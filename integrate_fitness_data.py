### Integrate long read library data with short read fitness data

import pandas as pd
import numpy as np
from os.path import join
from tqdm.auto import tqdm
from rapidfuzz import process, fuzz
from Levenshtein import distance
from scipy.sparse import vstack, lil_array, csr_array
from multiprocessing import Pool, cpu_count

LIBRARY_PATH = 's3://pioneer-sequencing/libraries'
MIN_BC_LENGTH = 42
MAX_BC_LENGTH = 46
N_CORES = 90 # Number of cores to use for multiprocessing
DISTANCE_CUTOFF = 3 # Use edit distance of 3 as cutoff for correction

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

def calculate_correction_map(lr_filter, sr_filter):
    '''
    Calculate correction map between short read to long read barcodes based on RapidFuzz
    '''
    def split_given_size(a, size):
        return np.split(a, np.arange(size,len(a),size))

    split_size = 1000 # Use chunks of 1000 barcodes
    queries = sr_filter
    choices = lr_filter
    splits = split_given_size(queries, split_size)

    min_dist = [] # Minimum pairwise edit distance between queries and choices
    arg_min = [] # Index of choice barcode with minimum pairwise edit distance

    # Process in batches to reduce memory overhead
    for split in tqdm(splits):
        distance_split = process.cdist(split, 
                                       choices, 
                                       scorer=distance, 
                                       score_cutoff=DISTANCE_CUTOFF, 
                                       dtype=np.uint8, workers=N_CORES)
        min_dist.append(distance_split.min(axis=1))
        arg_min.append(np.argmin(distance_split, axis=1))
    min_dist = np.concatenate(min_dist)
    arg_min = np.concatenate(arg_min)

    # Calculate correctable queries and choices
    queries_to_correct = np.where(min_dist <= DISTANCE_CUTOFF)[0]
    target_choices = arg_min[queries_to_correct]

    # Calculate barcode correction map
    correction_map = {
        queries[queries_to_correct[i]]: choices[target_choices[i]] for i in range(len(queries_to_correct))
    }

    # Verify that the mapping are the minimum distance
    assert np.array([distance(x,y)==z for x,y,z in zip(
        np.array(sr_filter)[queries_to_correct], 
        np.array(lr_filter)[target_choices],
        min_dist[queries_to_correct])]).all()
    
    # Verify the mapping
    assert np.array([correction_map[x] == y for x,y in zip(
        np.array(sr_filter)[queries_to_correct], 
        np.array(lr_filter)[target_choices])]).all()

    return correction_map

def create_integrated_dataframe(library, df_fitness, intersection, correction_map):
    '''
    Create an integrated dataframe with error corrected barcodes
    '''
    # Create corrected short-read barcode column
    df_fitness['bc_sequence'] = df_fitness['barcode'].map(correction_map)
    df_fitness.loc[df_fitness['barcode'].isin(intersection), 'bc_sequence'] = df_fitness['barcode']

    # Create correction status label
    df_fitness['correction_status'] = 'corrected'
    df_fitness.loc[df_fitness['barcode'].isin(intersection), 'correction_status'] = 'exact_match'
    df_fitness.loc[pd.isnull(df_fitness['bc_sequence']), 'correction_status'] = 'uncorrected'

    # Merge library data with fitness data
    df_fitness = df_fitness.rename(columns={'barcode':'uncorrected_bc_sequence'})
    merge = pd.merge(library, df_fitness, on='bc_sequence', how='right')

    return merge