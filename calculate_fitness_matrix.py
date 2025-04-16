### Calculate fold change trajectories of barcodes
import pandas as pd
import numpy as np
from tqdm import tqdm
from os.path import join

# Pseudocode
# Pivot table on nonzero days to get lr_corr wide format
# Get unique indices to subset the day 0 data (otherwise we are inefficient and use lots of useless barcodes)
# Calculate lr_corr - lr0_corr

def load_and_merge_data(results_path, metadata_file, min_counts=5):
    '''
    Load counts and metadata and merge
    '''
    # Load metadata
    metadata = pd.read_csv(join(results_path, metadata_file))

    # Load corrected counts from metadata samples
    samples = metadata['sample']

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
    
    # Merge counts with metadata
    counts_merge = pd.merge(counts, metadata, on=['sample'])
    
    return counts_merge

def calculate_psi_freq(counts_merge, base_timepoint=0):
    # Calculate psi-freq and psi-corrected lr
    # psi-freq is the derived pseudofrequency from the boba-seq paper, which works out to 1/sqrt(N0*N)
    # N0 for each sample
    N0_dict = counts_merge.loc[
        counts_merge['timepoint'] == base_timepoint, ['library','condition','replicate','N']].set_index(
        ['library','condition','replicate']).to_dict()

    # Calculate N0, psi, and corrected lr
    counts_merge['N0'] = [N0_dict['N'][x] for x in list(zip(
        counts_merge['library'],counts_merge['condition'],counts_merge['replicate']))]
    counts_merge['psi_freq'] = 1 / np.sqrt(counts_merge['N0'] * counts_merge['N'])
    counts_merge['lr_corr'] = np.log2(counts_merge['freq'] + counts_merge['psi_freq'])

    # Calculate psi_freq pivot table
    df_psi_freq = pd.pivot_table(counts_merge,
                            values='psi_freq',
                            index=['library','condition','replicate'],
                            columns='timepoint')
    df_psi_freq = df_psi_freq.loc[:,df_psi_freq.columns != base_timepoint]

    return df_psi_freq

def calculate_lr_df(df_counts, df_psi_freq, fill_value):
    '''
    Calculate pivoted dataframe with columns for days and lr values for a given condition (library x media x replicate)
    '''
    # Pivot
    df_pivot = pd.pivot_table(df_counts,
                          values='freq',
                          index=['library','condition','replicate','barcode'],
                          columns='timepoint',
                          fill_value=fill_value)

    # Merge psi_freq
    df_merge = pd.merge(df_pivot.reset_index(),
                        df_psi_freq.reset_index(),
                        on=['library','condition','replicate'],
                        suffixes=('_freq','_psi_freq')
                       ).set_index(['library','condition','replicate','barcode'])

    # Generate lr dataframe
    df_lr = pd.DataFrame(index=df_merge.index,
                         columns=df_pivot.columns)
    for d in df_lr.columns:
        df_lr[d] = np.log2(df_merge[str(d) + '_freq'] + df_merge[str(d) + '_psi_freq'])

    return df_lr

def calculate_fitness(counts_merge, df_psi_freq, keep_dropouts=True, base_timepoint=0):
    '''
    Calculate barcodes x fitness df
    '''
    if keep_dropouts:
        fill_value=0
    else:
        fill_value=np.nan
        
    df_fitness = []
    counts_day0 = counts_merge[
        counts_merge['timepoint'] == base_timepoint].set_index(['library','condition','replicate','barcode'])
    g = counts_merge[
        counts_merge['timepoint'] != base_timepoint].groupby(['library','condition','replicate']) # Skip base timepoint
    keys = list(g.groups.keys())
    
    for key in tqdm(keys):
        counts_sub = counts_merge.loc[g.groups[key]]
        if len(counts_sub) > 0:
            df_fitness_sub = calculate_lr_df(counts_sub, df_psi_freq, fill_value)
            counts_sub_merge = pd.merge(df_fitness_sub, counts_day0['lr_corr'],
                                        on=['library','condition','replicate','barcode'],
                                        how='left').rename(columns={'lr_corr':base_timepoint})
            counts_sub_merge = counts_sub_merge.reset_index().set_index(['library','condition','replicate','barcode'])
            df_fitness.append(counts_sub_merge)
    df_fitness = pd.concat(df_fitness)
    
    # Subtract lr by lr0 to get fitness
    for c in df_fitness.columns:
        df_fitness[c] = df_fitness[c] - df_fitness[base_timepoint]

    # Remove base timepoint column
    df_fitness = df_fitness.drop(columns=[base_timepoint])
    
    # Melt to reshape
    df_fitness = pd.melt(df_fitness.reset_index(),
            id_vars=['library','condition','replicate','barcode'],
            var_name='timepoint',
            value_name='fitness')

    return df_fitness

