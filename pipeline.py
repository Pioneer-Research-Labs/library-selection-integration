### Master pipeline
### Calculates fitness for each unique library/environment combination
### TODO: integration with long read data
import argparse
import pandas as pd
from os.path import join
from os import makedirs
from datetime import datetime
from calculate_fitness_matrix import (
    load_metadata,
    load_and_merge_data,
    calculate_psi_freq_v2,
    prep_and_filter_freq_table,
    calculate_fitness_final,
    generate_per_sample_QC_metrics
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
    parser.add_argument("--experiment_ID", type=str, required = True, help="Experiment ID. Used to sort results in selection_parent_dir and as a prefix for output files.")
    parser.add_argument("--short_read_path", type=str, required=True, help="Path to short read results.")
    parser.add_argument("--metadata", type=str, required=False, default='metadata.csv', help="Name of metadata file (default='metadata.csv').")
    parser.add_argument("--selection_parent_dir", type=str, required=False, default='s3://pioneer-analysis/library-selection-output', 
                        help="Parent directory for selection outputs. default='s3://pioneer-analysis/library-selection-output')")
    parser.add_argument("--base_timepoint", type=int, required=False, default=0, help="Base timepoint to use as reference (default=0).")
    parser.add_argument("--min_counts", type=int, required=False, default=0, help="Min counts set by short read pipeline (default = 0)")

    args = parser.parse_args()

    short_read_path = args.short_read_path
    metadata_file = args.metadata
    base_timepoint = args.base_timepoint
    experiment_id = args.experiment_ID

    # Get the current datetime object
    now = datetime.now()
    # Format the datetime object into the desired string format
    formatted_datetime = now.strftime("%Y_%m_%d_%H_%M_%S")

    output_path = join(args.selection_parent_dir, experiment_id, formatted_datetime)
    ### if you aren't on S3, make you parent directory
    if not output_path.startswith("s3:"):
        makedirs(output_path)

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

        # Prep table by removing barcodes not found at baseline and baseline + 1
        # adding 0 frequency for other timepoints where not detected
        print('Preparing frequency table for barcode correction...')
        prepped_freq_table = prep_and_filter_freq_table(counts_merge, base_timepoint=base_timepoint)

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
        out_name = experiment_id + '_selection_' + str(g[0]) + '_' + str(g[1])
        print('Saving integrated fitness data...')
        merge_fitness.to_parquet(join(output_path, out_name + '.parquet'), index=False)
        print('Data saved.')

        # Generate some exploratory QC metrics that can be used to assess sample quality
        print('Generating per-sample QC metrics...')
        qc_table = generate_per_sample_QC_metrics(merge_fitness)
        qc_table.to_csv(join(output_path, out_name + '_qc_metrics.csv'), index=False)
        print('QC table saved.')

