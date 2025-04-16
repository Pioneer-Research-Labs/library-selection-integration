### Master pipeline
### TODO: integration with long read data
import argparse
from os.path import join
from calculate_fitness_matrix import (
    load_and_merge_data,
    calculate_psi_freq,
    calculate_fitness
)

if __name__ == '__main__':
    # Get short read results path and metadata file path
    parser = argparse.ArgumentParser(
        description="Pipeline for calculating fitness and integrating with library data.")
    parser.add_argument("--short_read_path", type=str, required=True, help="Path to short read results.")
    parser.add_argument("--metadata", type=str, required=False, help="Name of metadata file.")
    parser.add_argument("--base_timepoint", type=int, required=False, help="Base timepoint to use as reference.")

    args = parser.parse_args()

    short_read_path = args.short_read_path
    if args.metadata:
        metadata_path = args.metadata
    else:
        metadata_path = 'metadata.csv'
    if args.base_timepoint:
        base_timepoint = args.base_timepoint
    else:
        base_timepoint = 0

    # Load data
    print('Loading data...')
    counts_merge = load_and_merge_data(short_read_path, metadata_file=metadata_path)
    print('Data loaded.')

    # Calculate psi-freq
    print('Calculating psi-freq...')
    df_psi_freq = calculate_psi_freq(counts_merge, base_timepoint=base_timepoint)

    # Calculate fitness
    print('Calculating fitness...')
    df_fitness = calculate_fitness(counts_merge, df_psi_freq, base_timepoint=base_timepoint)
    
    # Save fitness data
    print('Saving fitness data...')
    df_fitness.to_parquet(join(short_read_path, 'fitness.parquet'), index=False)