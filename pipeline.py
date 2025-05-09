### Master pipeline
### Calculates fitness for each unique library/strain/condition combination
### TODO: integration with long read data
import argparse
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
    calculate_distance_matrix,
    error_correct_barcodes,
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

    # For each library/condition/strain combination, run the pipeline
    groups = sample_dict.keys()
    for g in groups:
        samples = sample_dict[g]
        print('Processing group: ', g)

        # Load data
        print('Loading data...')
        counts_merge = load_and_merge_data(short_read_path, metadata, samples)
        print('Data loaded.')

        # Calculate psi-freq
        print('Calculating psi-freq...')
        df_psi_freq = calculate_psi_freq(counts_merge, base_timepoint=base_timepoint)

        # Calculate fitness
        print('Calculating fitness...')
        df_fitness = calculate_fitness(counts_merge, df_psi_freq, base_timepoint=base_timepoint)
        
        # Save fitness data
        out_name = out_prefix + '_' + str(g[0]) + '_' + str(g[1]) + '_' + str(g[2])
        print('Saving raw fitness data...')
        df_fitness.to_parquet(join(short_read_path, out_name + '.parquet'), index=False)
        print('Data saved.')

        # Load library data
        library_id = g[0].split('-')[0] # Only use the prefix library ID 
        print('Loading library data...')
        library = load_library_data(library_id)
        print('Library data loaded.')

        # Merge library data with fitness data
        print('Correcting barcodes...')
        
        # Get filtered barcodes to correct
        lr_filter, sr_filter, intersection = get_filtered_barcodes_to_correct(library, df_fitness)
        
        # Calculate distance matrix
        distances = calculate_distance_matrix(lr_filter, sr_filter)

        # Do the actual error correction
        correction_map, correctable_bcs = error_correct_barcodes(lr_filter, sr_filter, distances)

        # Create integrated dataframe
        merge = create_integrated_dataframe(
            library, df_fitness, intersection, correction_map, correctable_bcs)
        
        print('Barcodes corrected.')

        # Save merged data
        out_name = out_prefix + '_integrated_' + str(g[0]) + '_' + str(g[1]) + '_' + str(g[2])
        print('Saving integrated fitness data...')
        merge.to_parquet(join(short_read_path, out_name + '.parquet'), index=False)
        print('Data saved.')