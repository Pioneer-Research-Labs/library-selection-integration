### Master pipeline
### Calculates fitness for each unique library/environment combination
### TODO: integration with long read data
import argparse
import pandas as pd
from os.path import join
from calculate_fitness_matrix import (
    load_metadata,
    load_and_merge_data,
    calculate_psi_freq_v2,
    prep_and_filter_freq_table,
    calculate_fitness_final
)
from integrate_fitness_data import (
    load_library_data,
    filter_corrected_library_data,
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
    parser.add_argument("--out_path", type=str, required=False, default=None, help="Path to write out files, if not specified is short_read_path")
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
    if args.out_path is None:
        out_path = short_read_path
    else:
        out_path = args.out_path

    # Load metadata
    print('Loading metadata...')
    sample_dict, metadata = load_metadata(short_read_path, metadata_file)
    print('Metadata loaded.')

    # For each library/environment/base_library combination, run the pipeline
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
        df_psi_freq = calculate_psi_freq_v2(counts_merge, base_timepoint=base_timepoint)

        # Calculate fitness
        print('Calculating fitness...')
        prepped_freq_table = prep_and_filter_freq_table(counts_merge, df_psi_freq, base_timepoint=base_timepoint)


        # Load base library data
        library_id = g[2].split('-')[0] # Only use the prefix library ID 
        print('Loading base library data...')
        base_library = load_library_data(library_id)
        print('Base library data loaded.')
        
        if 'library_correction_status' in base_library.columns:
            base_library = filter_corrected_library_data(base_library)

        # Merge library data with fitness data
        print('Correcting barcodes...')
        
        # Get filtered barcodes to correct
        lr_filter, sr_filter, intersection = get_filtered_barcodes_to_correct(base_library, prepped_freq_table)
        
        # Calculate correction map
        correction_map = calculate_correction_map(lr_filter, sr_filter)

        print('Exact match barcodes: ', len(intersection))
        print('Correctable barcodes: ', len(correction_map))

        # Create integrated dataframe
        merge = create_integrated_dataframe(
            base_library, prepped_freq_table, intersection, correction_map)
        
        print('Barcodes corrected.')

        print('Uncorrected barcodes: ', 
              merge["ngs_correction_status"].apply(lambda x:x == ["uncorrected"]).sum())

        merge_fitness = calculate_fitness_final(frequency_table = merge, psi_freq_table = df_psi_freq)

        # Save merged data
        out_name = out_prefix + '_integrated_' + str(g[0]) + '_' + str(g[1])
        print('Saving integrated fitness data...')
        merge_fitness.to_parquet(join(out_path, out_name + '.parquet'), index=False)
        print('Data saved.')