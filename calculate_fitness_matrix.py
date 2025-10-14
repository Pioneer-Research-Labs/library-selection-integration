### Calculate fold change trajectories of barcodes
import pandas as pd
import numpy as np
from tqdm import tqdm
from os.path import join

# Pseudocode
# Pivot table on nonzero days to get lr_corr wide format
# Get unique indices to subset the day 0 data (otherwise we are inefficient and use lots of useless barcodes)
# Calculate lr_corr - lr0_corr

def load_metadata(results_path, metadata_file):
    '''
    Load metadata from file. Returns a tuple of length two comprised of:
    1. A dictionary of sample groupings based on library, environment and base library
    2. A pandas DataFrame containing the metadata loaded from 'metadata_file'
    '''
    # Load metadata
    metadata = pd.read_csv(join(results_path, metadata_file))

    # Get unique library/environment groups and corresponding sample IDs
    groups = metadata.groupby(['library','environment','base_library']).groups
    metadata.set_index(['library','environment','base_library'], inplace=True)
    metadata.sort_index(inplace=True)

    # Convert sample IDs to string
    metadata['sample'] = metadata['sample'].astype(str)

    sample_dict = {}

    for g in groups:
        sample_dict[g] = metadata.loc[g, 'sample'].values

    metadata.reset_index(inplace=True)
    
    return sample_dict, metadata
def load_and_merge_data(results_path: str, metadata: pd.DataFrame, samples: list, min_counts: int) -> pd.DataFrame:
    '''
    Load counts and metadata and merge
    '''
    # Load counts from each sample and merge
    counts = []
    for sample in tqdm(samples, desc='Loading data'):
        try:
            counts.append(pd.read_csv(join(
                results_path, sample, sample + '_' + str(min_counts) + '_barcodes_freq.csv')))
        except FileNotFoundError:
            print('No counts file found for sample: ', sample)
            continue
    counts = pd.concat(counts)

    # Convert samples to string if not already
    counts['sample'] = counts['sample'].astype(str)
    
    # Merge counts with metadata
    counts_merge = pd.merge(counts, metadata, on=['sample'])
    
    return counts_merge

def calculate_psi_freq(counts_merge: pd.DataFrame, base_timepoint=0) -> pd.DataFrame:
    '''
    Given a barcode frequency pandas DataFrame with columns library, environment, replicate, and N,
    return a pandas DataFrame with row indexes of library,environment,replicate and columns corresponding
    to non-baseline timepoints and value of psi in frequency as calculated in the Boba-Seq paper.

    The original counts_merge table is also modified to add N0, psi_freq, and lr_corr (psi-corrected log2 abundance)
    '''
    # Calculate psi-freq and psi-corrected lr
    # psi-freq is the derived pseudofrequency from the boba-seq paper, which works out to 1/sqrt(N0*N)
    # N0 for each sample

    ### Grab N from the baseline timepoint and make a single dciontary in the form of {N: {lib,env,rep: N0}}
    N0_dict = counts_merge.loc[
        counts_merge['timepoint'] == base_timepoint, ['library','environment','replicate','N']].set_index(
        ['library','environment','replicate']).to_dict()

    # Calculate N0, psi, and corrected lr
    counts_merge['N0'] = [N0_dict['N'][x] for x in list(zip(
        counts_merge['library'],counts_merge['environment'],counts_merge['replicate']))]
    counts_merge['psi_freq'] = 1 / np.sqrt(counts_merge['N0'] * counts_merge['N'])
    counts_merge['lr_corr'] = np.log2(counts_merge['freq'] + counts_merge['psi_freq'])

    # Calculate psi_freq pivot table
    df_psi_freq = pd.pivot_table(counts_merge,
                            values='psi_freq',
                            index=['library','environment','replicate'],
                            columns='timepoint')
    df_psi_freq = df_psi_freq.loc[:,df_psi_freq.columns != base_timepoint]

    return df_psi_freq

def calculate_lr_df(df_counts, df_psi_freq, fill_value):
    '''
    Calculate pivoted dataframe with columns for days and lr values for a given environment (library x media x replicate)
    '''
    # Pivot data without TMPT 0. This is where we lose barcodes that are not present in at least 1 post-baseline timepoint!
    df_pivot = pd.pivot_table(df_counts,
                          values='freq',
                          index=['library','environment','replicate','barcode'],
                          columns='timepoint',
                          fill_value=fill_value)

    # Merge psi_freq to generate a table where each timepoint has a <tmpt_freq> and <tmpt_psi_freq> column
    # using the suffix keys
    df_merge = pd.merge(df_pivot.reset_index(),
                        df_psi_freq.reset_index(),
                        on=['library','environment','replicate'],
                        suffixes=('_freq','_psi_freq')
                       ).set_index(['library','environment','replicate','barcode'])
    
    # Generate lr dataframe
    df_lr = pd.DataFrame(index=df_merge.index,
                         columns=df_pivot.columns)
    for d in df_lr.columns:
        df_lr[d] = np.log2(df_merge[str(d) + '_freq'] + df_merge[str(d) + '_psi_freq'])

    return df_lr

def calculate_fitness(counts_merge, df_psi_freq, keep_dropouts=True, base_timepoint=0):
    '''
    Calculate barcodes x fitness df
    Note that we only keep barcodes that appear at least once in the non-base timepoint
    '''
    if keep_dropouts:
        fill_value=0
    else:
        fill_value=np.nan
        
    df_fitness = []

    ### Pull out the base_timepoint data and set the index as the unique comination of library,env,replicate,barcode
    counts_day0 = counts_merge[
        counts_merge['timepoint'] == base_timepoint].set_index(['library','environment','replicate','barcode'])
    ### Grab all sets of library, environment, replicate that aren't baseline
    g = counts_merge[
        counts_merge['timepoint'] != base_timepoint].groupby(['library','environment','replicate']) # Skip base timepoint
    keys = list(g.groups.keys())
    ### Loop through each
    for key in tqdm(keys):
        ### g.groups[key] returns the set of indexes in counts_merge that are in that group which allows the subset by loc
        ### This could presumably be done by operating on the groups in the groupby?
        counts_sub = counts_merge.loc[g.groups[key]]
        if len(counts_sub) > 0:
            ### Return a dataframe with index of library, env, rep, barcode with columns of timepoints and values of
            ### psi corrected log abundance. Will be missing all barcodes without at least one post-baseline detection.
            df_fitness_sub = calculate_lr_df(counts_sub, df_psi_freq, fill_value)
            ### In this merge we lose the baseline barcodes not detected at later timepoint
            ### The existence of lr_corr on counts_day0 is a consequence of running calculate_psi_freq on the same table and modifying in place
            counts_sub_merge = pd.merge(df_fitness_sub, counts_day0['lr_corr'],
                                        on=['library','environment','replicate','barcode'],
                                        how='left').rename(columns={'lr_corr':base_timepoint})
            counts_sub_merge = counts_sub_merge.reset_index().set_index(
                ['library','environment','replicate','barcode'])
            df_fitness.append(counts_sub_merge)
    df_fitness = pd.concat(df_fitness)
    
    # Subtract lr by lr0 to get fitness
    for c in df_fitness.columns:
        df_fitness[c] = df_fitness[c] - df_fitness[base_timepoint]
    
    # Melt to reshape
    df_fitness = pd.melt(df_fitness.reset_index(),
            id_vars=['library','environment','replicate','barcode'],
            var_name='timepoint',
            value_name='fitness')
    
    # Drop nan fitnesses which occur when a barcode is detected at a later timepoint but not baseline
    df_fitness = df_fitness.dropna(subset=['fitness'])

    # Drop duplicate rows ### unclear why this would occur
    df_fitness.drop_duplicates(inplace=True)

    # Merge frequency information with fitness data -- change 
    merge_cols = ['library', 'environment', 'replicate','timepoint','barcode','n','freq']
    df_fitness = pd.merge(
        df_fitness, counts_merge[merge_cols], on=['library', 'environment', 'replicate', 'timepoint', 'barcode'], how='left')

    ### Set un
    # Keep only barcodes that exist in t=1 timepoint
    print('Filtering for barcodes that exist in first nonzero timepoint.')
    print('Number of barcodes before filtering: ', len(df_fitness.barcode.unique()))

    timepoint_counts = df_fitness[df_fitness.timepoint == base_timepoint + 1].groupby(
        ['library','environment','replicate','barcode'])['n'].count()
    bcs_to_keep = timepoint_counts[timepoint_counts > 0].index
    df_fitness.set_index(['library','environment','replicate','barcode'], inplace=True)
    df_fitness = df_fitness.loc[bcs_to_keep]
    df_fitness.reset_index(inplace=True)

    print('Number of barcodes after filtering: ', len(df_fitness.barcode.unique()))

    return df_fitness



def calculate_psi_freq_v2(counts_merge: pd.DataFrame, base_timepoint=0) -> pd.DataFrame:
    '''
    Given a barcode frequency pandas DataFrame with columns library, environment, replicate, and N,
    return a pandas DataFrame with row indexes of library,environment,replicate and columns corresponding
    to non-baseline timepoints and value of psi in frequency as calculated in the Boba-Seq paper.

    The original counts_merge table is also modified to add N0, psi_freq, and lr_corr (psi-corrected log2 abundance)
    '''
    # Calculate psi-freq and psi-corrected lr
    # psi-freq is the derived pseudofrequency from the boba-seq paper, which works out to 1/sqrt(N0*N)
    # N0 for each sample

    ### Grab N from the baseline timepoint and make a single dciontary in the form of {N: {lib,env,rep: N0}}
    N0_dict = counts_merge.loc[
        counts_merge['timepoint'] == base_timepoint, ['library','environment','replicate','N']].set_index(
        ['library','environment','replicate']).to_dict()

    # Calculate N0, psi, and corrected lr
    counts_merge['N0'] = [N0_dict['N'][x] for x in list(zip(
        counts_merge['library'],counts_merge['environment'],counts_merge['replicate']))]
    counts_merge['psi_freq'] = 1 / np.sqrt(counts_merge['N0'] * counts_merge['N'])
    #counts_merge['lr_corr'] = np.log2(counts_merge['freq'] + counts_merge['psi_freq'])

    # Calculate psi_freq pivot table
    df_psi_freq = pd.pivot_table(counts_merge,
                            values='psi_freq',
                            index=['library','environment','replicate'],
                            columns='timepoint').reset_index()
    
    df_psi_freq_long = df_psi_freq.melt(id_vars = ['library','environment','replicate'], value_name="psi_freq")

    return df_psi_freq_long


def prep_and_filter_freq_table(counts_merge: pd.DataFrame, keep_dropouts=True, base_timepoint=0):
    '''
    Prep the counts table for downstream calculations by:
    1. Merge on baseline frequency values by the combination library, environment, replicate, barcode
    2. Fill in barcode data such that barcodes detected at a single timepoint are measured as frequency 0 at other timepoints rather than missing (if keep_dropouts = True)
    3. Filter out barcodes that are not detected in the baseline and baseline + 1 timepoint
    
    Returns a pandas DataFrame akin to counts_merged but with additional columns, rows, and filtering as described above.
    '''
    if keep_dropouts:
        fill_value=0
    else:
        fill_value=np.nan
        
    df_out = []

    ### Pull out the base_timepoint data and set the index as the unique comination of library,env,replicate,barcode
    baseline_freqs = counts_merge.loc[
        counts_merge['timepoint'] == base_timepoint,
        ['library','environment','replicate','barcode','freq']
        ].rename(columns={'freq':'baseline_freq'})
    
    ### Grab all sets of library, environment
    g = counts_merge.groupby(['library','environment','replicate'])
    keys = list(g.groups.keys())
    ### Loop through each
    for key in tqdm(keys):
        ### g.groups[key] returns the set of indexes in counts_merge that are in that group which allows the subset by loc
        ### This could presumably be done by operating on the groups in the groupby?
        counts_sub = counts_merge.loc[g.groups[key]]
        if len(counts_sub) > 0:

            # Pivot data so that we add 0 counts for everything that isn't found in all timepoints and then make long again
            df_pivot_wide = pd.pivot_table(counts_sub,
                          values='freq',
                          index=['library','environment','replicate','barcode'],
                          columns='timepoint',
                          fill_value=fill_value).reset_index()

                
            df_melt_long = df_pivot_wide.melt(id_vars=['library','environment','replicate','barcode'], value_name="freq")

            #df_merge = pd.merge(df_melt_long,
            #        df_psi_freq,
            #            on=['library','environment','replicate','timepoint']
            #           )
            
            ### Here we drop all things at later time points that aren't found in the baseline timepoint.
            df_merge = pd.merge(df_melt_long,
                    baseline_freqs,
                        on=['library','environment','replicate', 'barcode']
                       )
                        
            df_out.append(df_merge)

    df_out = pd.concat(df_out)

    # Drop duplicate rows ### unclear why this would occur
    df_out.dropna(inplace=True)
    df_out.drop_duplicates(inplace=True)

    # Merge n information with fitness data
    merge_cols = ['library', 'environment', 'replicate','timepoint','barcode','n']
    df_out = pd.merge(
        df_out, counts_merge[merge_cols], on=['library', 'environment', 'replicate', 'timepoint', 'barcode'], how='left')
    
    df_out['n'] = df_out['n'].fillna(0)
    bigN_table = counts_merge[['library', 'environment', 'replicate', 'timepoint', 'N']].drop_duplicates()
    df_out = pd.merge(df_out, bigN_table, on=['library', 'environment', 'replicate', 'timepoint'])

    # Keep only barcodes that exist in t=1 timepoint
    print('Filtering for barcodes that exist in first nonzero timepoint.')
    print('Number of barcodes found in at least one baseline sample: ', len(df_out.barcode.unique()))

    timepoint_counts = df_out[df_out.timepoint == base_timepoint + 1].groupby(
        ['library','environment','replicate','barcode'])['n'].sum()
    bcs_to_keep = timepoint_counts[timepoint_counts > 0].index
    df_out.set_index(['library','environment','replicate','barcode'], inplace=True)
    df_out = df_out.loc[bcs_to_keep]
    df_out.reset_index(inplace=True)

    print('Number of barcodes after filtering: ', len(df_out.barcode.unique()))

    return df_out


def calculate_fitness_final(frequency_table, psi_freq_table, freq_column='total_freq', bl_freq_column='total_bl_freq'):
    fitness_dataframe = pd.merge(frequency_table, psi_freq_table, on=['library','environment','replicate','timepoint'])
    
    fitness_dataframe["fitness"] =  np.log2((fitness_dataframe[freq_column] + fitness_dataframe["psi_freq"])) - np.log2((fitness_dataframe[bl_freq_column] + fitness_dataframe["psi_freq"]))

    return fitness_dataframe