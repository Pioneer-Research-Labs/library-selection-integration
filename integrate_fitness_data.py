### Integrate long read library data with short read fitness data

import pandas as pd
import numpy as np
from os.path import join
from tqdm import tqdm
from Levenshtein import distance
from scipy.sparse import vstack, lil_array, csr_array
from multiprocessing import Pool, cpu_count

LIBRARY_PATH = 's3://pioneer-sequencing/libraries'
MIN_BC_LENGTH = 42
MAX_BC_LENGTH = 46
N_CORES = 32
DISTANCE_CUTOFF = 2 # This corresponds to ~0.0001 probability of distances >=3 for a 44bp barcode

def compute_wrapper(args):
    '''
    Wrapper function to unpack arguments for multiprocessing
    '''
    return compute_distances_row(*args)

def compute_distances_row(lr, sr_bcs, distance_cutoff):    
    '''
    Compute Levenshtein distances between a long read and all short reads in the short read barcode list
    '''
    distances = csr_array((1, len(sr_bcs)))

    for i in range(len(sr_bcs)):
        d = distance(lr, sr_bcs[i], score_cutoff=distance_cutoff)
        if d <= distance_cutoff:
            distances[0, i] = d

    return distances

def load_library_data(library_id):
    '''
    Load object for long read library data
    '''
    # Library data
    library = pd.read_parquet(join(LIBRARY_PATH, library_id + '.parquet'))

    return library

def get_filtered_barcodes_to_correct(library, df_fitness):
    '''
    Get barcodes from long read library data and short read fitness data for error correction.
    We won't bother correcting the empty barcodes because there are a lot of them.
    '''
    # Get the barcodes
    sr_bcs = list(df_fitness.barcode.unique()) # Short read
    df_sr_bcs = pd.DataFrame(index=sr_bcs, 
                            columns=['length'],
                            data=[len(x) for x in sr_bcs]) # dataframe for short read
    
    # For now let's only look at barcodes between a certain length (i.e. 2 Levenshtein distances away from 44)
    lr_filter = list(library[library.bc_length.between(MIN_BC_LENGTH, MAX_BC_LENGTH)].bc_sequence)
    sr_filter = list(df_sr_bcs[df_sr_bcs.length.between(MIN_BC_LENGTH, MAX_BC_LENGTH)].index)

    # We only need to correct barcodes that are not exact matches between the two datasets
    intersection = set(lr_filter).intersection(set(sr_filter))

    lr_filter = list(set(lr_filter) - intersection)
    sr_filter = list(set(sr_filter) - intersection)
    
    return lr_filter, sr_filter, intersection

def calculate_distance_matrix(lr_filter, sr_filter):
    '''
    Calculate Levenshtein distance matrix between correctable long read and short read barcodes
    using multiprocessing with a progress bar
    '''
    # Generate inputs
    inputs = lr_filter

    # Prepare arguments as tuples for each row
    task_inputs = [(row, sr_filter, DISTANCE_CUTOFF) for row in inputs]

    # Set number of processes
    n_cores = min(cpu_count(), N_CORES)

    # Use multiprocessing Pool
    with Pool(processes=n_cores) as pool:
        results = list(tqdm(pool.imap(
            compute_wrapper, task_inputs), total=len(task_inputs), desc="Calculating Levenshtein distances"))
    # results = list(map(compute_wrapper, task_inputs))

    # Stack sparse arrays to get the final distance matrix
    distances = vstack(results)

    return distances

def error_correct_barcodes(lr_filter, sr_filter, distances):
    '''
    Error correct long read barcodes based on Levenshtein distance matrix
    '''
    # Get indices of nonzero distances less than cutoff
    nz_indices = distances.nonzero()
    np_lr_sub = np.array(lr_filter) # Numpy array of long read barcodes
    np_sr_sub = np.array(sr_filter) # Numpy array of short read barcodes

    # Make dataframe of correctable barcodes
    correctable_bcs = pd.DataFrame(index=np.unique(np_lr_sub[nz_indices[0]]),
                                columns=['distance','possible_mapped_bc'],
                                data=[])

    # Pull out distances, mapped barcodes, and downstream stats
    for i, j in zip(nz_indices[0], nz_indices[1]):
        try:
            correctable_bcs.loc[np_lr_sub[i],'distance'].append(distances[i,j])
            correctable_bcs.loc[np_lr_sub[i],'possible_mapped_bc'].append(np_sr_sub[j])
        except AttributeError:
            correctable_bcs.loc[np_lr_sub[i],'distance'] = [distances[i,j]]
            correctable_bcs.loc[np_lr_sub[i],'possible_mapped_bc'] = [np_sr_sub[j]]
    correctable_bcs['n'] = [len(x) for x in correctable_bcs['distance']]
    correctable_bcs['min_dist'] = [min(x) for x in correctable_bcs['distance']]
    correctable_bcs['n_min_dist'] = [
        sum(np.array(x)==y) for x,y in zip(correctable_bcs['distance'],correctable_bcs['min_dist'])]

    # Filter out barcodes that have more than one mappable target based on minimum distance
    correctable_bcs = correctable_bcs[correctable_bcs.n_min_dist == 1]
    correctable_bcs['mapped_bc'] = [np.array(x)[np.array(y) == z][0] for x,y,z in zip(correctable_bcs['possible_mapped_bc'],
                                                                        correctable_bcs['distance'],
                                                                        correctable_bcs['min_dist'])]

    # Save the mapping dictionary
    correction_map = correctable_bcs['mapped_bc'].to_dict()

    # Verify that the mapping are the minimum distance
    assert np.array([distance(x,y)==z for x,y,z in zip(
        correctable_bcs.index, 
        correctable_bcs.mapped_bc,
        correctable_bcs.min_dist)]).all()

    return correction_map, correctable_bcs

def create_integrated_dataframe(library, df_fitness, intersection, correction_map, correctable_bcs):
    '''
    Create an integrated dataframe with error corrected barcodes
    '''
    # Create a new corrected library dataframe
    library_corr = library[(library.bc_sequence.isin(list(intersection))) |
                           (library.bc_sequence.isin(correctable_bcs.index))]

    # Create a new column for correction status
    library_corr.loc[:,'corrected'] = False
    library_corr['corrected_bc_sequence'] = library_corr['bc_sequence'].copy()
    library_corr.loc[library_corr.bc_sequence.isin(correctable_bcs.index),'corrected'] = True
    library_corr.loc[library_corr.bc_sequence.isin(correctable_bcs.index),'bc_sequence'] = library_corr.loc[
        library_corr.bc_sequence.isin(correctable_bcs.index),'corrected_bc_sequence'].map(correction_map)
    library_corr.set_index('corrected_bc_sequence', inplace=True)

    # Merge with fitness data
    df_fitness.set_index('barcode', inplace=True)
    merge = pd.merge(library_corr, df_fitness, left_index=True, right_index=True, how='inner').reset_index()

    return merge