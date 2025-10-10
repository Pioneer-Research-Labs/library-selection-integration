### Master pipeline
### Calculates fitness for each unique library/environment combination
### TODO: integration with long read data
import argparse
import pandas as pd
from os.path import join
from calculate_fitness_matrix import (
    load_metadata,
    load_and_merge_data,
    calculate_psi_freq,
    calculate_fitness
)
from integrate_fitness_data import (
    load_library_data,
    get_filtered_barcodes_to_correct,
    calculate_correction_map,
    create_integrated_dataframe  
)

if __name__ == '__main__':
    # Get short read results path and metadata filename
    parser = argparse.ArgumentParser(
        description="Pipeline for calculating fitness and integrating with library data.")
    parser.add_argument("--short_read_path", type=str, required=True, help="Path to short read results.")
    parser.add_argument("--metadata", type=str, required=False, help="Name of metadata file (default='metadata.csv').")
    parser.add_argument("--out_prefix", type=str, required=False, help="Name of output fitness file prefix (default='fitness').")
    parser.add_argument("--base_timepoint", type=int, required=False, help="Base timepoint to use as reference (default=0).")
    parser.add_argument("--min_counts", type=int, required=False, default=0, help="Min counts set by short read pipeline (default = 0)")

    args = parser.parse_args()

    short_read_path = args.short_read_path
    if args.metadata:
        metadata_file = args.metadata
    else:
        metadata_file = 'metadata.csv'
    if args.out_prefix:
        out_prefix = args.out_prefix
    else:
        out_prefix = 'fitness'
    if args.base_timepoint:
        base_timepoint = args.base_timepoint
    else:
        base_timepoint = 0

    # Load metadata
    print('Loading metadata...')
    sample_dict, metadata = load_metadata(short_read_path, metadata_file)
    print('Metadata loaded.')

    # For each library/environment/ combination, run the pipeline
    groups = sample_dict.keys()
    for g in groups:
        samples = sample_dict[g]
        print('Processing group: ', g)

        # Load data
        print('Loading data...')
        counts_merge = load_and_merge_data(short_read_path, metadata, samples, args.min_counts)
        print('Data loaded.')

        # Calculate psi-freq
        print('Calculating psi-freq...')
        df_psi_freq = calculate_psi_freq(counts_merge, base_timepoint=base_timepoint)

        # Calculate fitness
        print('Calculating fitness...')
        df_fitness = calculate_fitness(counts_merge, df_psi_freq, base_timepoint=base_timepoint)
        
        # Save fitness data
        out_name = out_prefix + '_' + str(g[0]) + '_' + str(g[1])
        print('Saving raw fitness data...')
        df_fitness.to_parquet(join(short_read_path, out_name + '.parquet'), index=False)
        print('Data saved.')

        # Load base library data
        library_id = g[2].split('-')[0] # Only use the prefix library ID 
        print('Loading base library data...')
        base_library = load_library_data(library_id)
        print('Base library data loaded.')

        # Merge library data with fitness data
        print('Correcting barcodes...')
        
        # Get filtered barcodes to correct
        lr_filter, sr_filter, intersection = get_filtered_barcodes_to_correct(base_library, df_fitness)

        
        # Calculate correction map
        correction_map = calculate_correction_map(lr_filter, sr_filter)

        print('Exact match barcodes: ', len(intersection))
        print('Correctable barcodes: ', len(correction_map))

        # Create integrated dataframe
        merge = create_integrated_dataframe(
            base_library, df_fitness, intersection, correction_map)
        
        print('Barcodes corrected.')

        print('Uncorrected barcodes: ', 
              len(merge[merge['ngs_correction_status'] == 'uncorrected'].uncorrected_bc_sequence.unique()))

        # Save merged data
        out_name = out_prefix + '_integrated_' + str(g[0]) + '_' + str(g[1])
        print('Saving integrated fitness data...')
        merge.to_parquet(join(short_read_path, out_name + '.parquet'), index=False)
        print('Data saved.')