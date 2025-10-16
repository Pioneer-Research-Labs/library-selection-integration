# 📈 Fitness Calculation & Integration Pipeline

This Python-based pipeline calculates **fitness metrics** for sequencing data across multiple **library**, **strain**, **replicate**, and **condition** combinations for a selection experiment. It integrates fitness results with PacBio library data and performs error correction on barcodes to the libraries.

---

## 🧬 Overview

This pipeline performs the following core functions:

1. **Loads and merges short-read data and metadata**
2. **Filters short barcodes within replicate-condition combination to those found in baseline and baseline + 1 timepoint**
3. **Corrects short read barcodes to library barcodes**
4. **Combines count and frequency data for short read barcodes that map to the same library barcode**
5. **Merges on library insert and mapping information**
6. **Calculates fitness for each combination barcode/replicate/condition for each timepoint comparing to baseline using the BOBA-seq method**

---

## 🛠️ Requirements

Python 3.7+

Install dependencies (if not already handled):

```bash
pip install pandas pyarrow 
```

Ensure the following scripts are in the same directory or importable:

* `calculate_fitness_matrix.py`
* `integrate_fitness_data.py`

---

## ▶️ Usage

```bash
python pipeline.py --short_read_path <SHORT_READ_DIR> [--metadata <METADATA_CSV>] /
    [--out_prefix <OUTPUT_PREFIX>] [--out_path <OUTPUT_PATH>] /
    [--base_timepoint <TIMEPOINT>] [--min_counts <MIN_COUNTS>]
```

### Arguments:

| Argument            | Type  | Required | Description                                             |
| ------------------- | ----- | -------- | ------------------------------------------------------- |
| `--short_read_path` | `str` | ✅        | Path to directory containing short-read results (this should usually be `results`).        |
| `--metadata`        | `str` | ❌        | Metadata CSV filename (default: `metadata.csv`).        |
| `--out_prefix`      | `str` | ❌        | Prefix for output files (default: `fitness`).           |
| `--out_path`        | `str` | ❌        | Output path for parquet files (default: `short_read_path`).  |
| `--base_timepoint`  | `int` | ❌        | Timepoint used as the fitness reference (default: `0`). |
| `--min_counts`      | `int` | ❌        | Min count for barcode reporting as set in short-read pipeline. (default: `0`) |
---

## 📂 Input Files

1. A metadata file containing information linking samples to experimental conditions. Must contain at minimum the following columns in any order:
   * sample (str) -- a sample identifier that was used for the short-read-pipeline
   * library (str) -- a library identifier corresponding the the library that was in the selection experiment (i.e. the base library inserted into a host)
   * environment (str) -- the condition in which the selection occurred
   * base_library(str) -- the donor library ID that has been characterized by long read sequencing and is deposited in 's3://pioneer-sequencing/libraries'
   * timepoint (int) -- must contain the value of base_timepoint and base_timepoint + 1 based on current filtering strategies
   * replicate (int/str) -- which replicate

2. Barcode count files generated from running `short-read-pipeline`. The output files are located in: short_read_path/sample/sample_min_counts_barcodes_freq.csv

---

## 📤 Output

For each unique `(library, environment)` group, a single output files is generated:

`fitness_integrated_<library>_<environment>.parquet`
   → Fitness data with library integration and barcode correction

Output columns:
* bc_length -> library_correction_status -- Data derived from the underlying base_library file that is merged on to short read data.
* bc_sequence: The error-corrected sequence used to merge short read barcodes to the long read library data. Occurs exactly once per combination of library, environment, replicate, and timepoint.
* library -- Library name from metadata file
* environment -- Environment from metadata file
* replicate -- Replicate from metadata file
* timepoint -- Timepoint from metadata file
* N -- total number of reads found in this library/environment/replicate/timepoint
* total_n -- Sum of read counts for all short read barcodes collapsed into this bc_sequence
* total_freq -- Sum of frequencies for all short read barcodes collapsed into this bc_sequence
* total_bl_freq -- Sum of baseline frequences for all short read barcodes collapsed into this bc_sequence
* n -- List of individual read counts for each short read barcode collapsed into the bc_sequence
* freq -- List of individual frequencies for each short read barcode collapsed into the bc_sequence
* bl_freq -- List of individual baseline frequences for each short read barcode collapsed into the bc_sequence
* uncorrected_bc -- List of original uncorrected BCs for each short read barcode collapsed into the bc_sequence
* ngs_correction_status -- List of NGS correction status for each short read barcode collapsed into the bc_sequence.
   * exact_match -- barcode found directly in library bc_sequence matches uncorrected_bc for this entry
   * corrected == was corrected to the bc_sequence. bc_sequence does NOT match uncorrected_bc for this entry
   * uncorrected == was not observed in the library and could not be corrected. The bc_sequence matches the uncorrected_bc for this entry, however there will be no corresponding library level information.
* psi_freq: The pseudocount frequency value derived from the psi calculation comparing this sample and its correponding baseline. Added to total_freq and total_bl_freq before fitness is calculated.
* Fitness: The log2 fold change of a bc_sequence frequncy compared to its baseline frequency after addition of psi_freq.

`fitness_integrated_<library>_<environment>_qc_metrics.csv` -> Per-sample QC metrics including information about total barcodes, 
number of barcodes detected, and basic frequency and fitness information.

Output columns:
* Sequencing info
   * n_short_reads -- total N of reads that came out of the short read pipeline, used for calculating frequencies
* All error-corrected barcodes (bc_sequences)
   * bc_n -- number barcodes
   * bc_n_unique -- number of unique barcodes. Expected to match bc_n but is left in for pipeline debugging
   * bc_n_detected -- how many barcodes are detected above 0
   * bc_median_fitness -- median fitness value for all barcodes
   * bc_freq_sum -- sum of frequencies for all barcodes surviving analysis. Should be near 1 unless lots of data lost via filtering.
* Barcodes that match the underlying long read library
   * bc_lib_matched_n -- number of barcodes that are matched to the underlying long read library
   * bc_lib_matched_n_detected -- number of barcodes that match the underlying long read library and are detected
   * bc_lib_matched_freq_sum -- summed frequencye of barcodes that match the underlying long read library
* Empty insert barcodes
   * empty_bc_n -- number of barcodes
   * empty_bc_n_detected -- number of barcodes detected
   * empty_bc_median_fitness -- median fitness of empty inserts
   * empty_bc_freq_sum -- summed frequency of empty inserts

## Using the output

The resulting parquet file can be read directly into python from s3 using:

```python
import pandas as pd
myData = pd.read_parquet("path/to/parquet/output")
```

To unnest the short-read barcode level data that has been collapsed into corresponding bc_sequences and
calculate individual uncorrected_bc level fitness values:

```python
import pandas as pd
import numpy as np

### Read in Data
myData = pd.read_parquet("path/to/parquet/output")

### Explode data to make it long by nested columns
myData_uncorrected_bc_level = myData.explode(column = ["n", "freq", "bl_freq", "uncorrected_bc", "ngs_correction_status"])

### Fix data types
myData_uncorrected_bc_level['freq'] = myData_uncorrected_bc_level['freq'].astype('float64')
myData_uncorrected_bc_level['bl_freq'] = myData_uncorrected_bc_level['bl_freq'].astype('float64')

### Calculate per-uncorrected barcode fitness
myData_uncorrected_bc_level["uncorrected_bc_fitness"] = np.log2(myData_uncorrected_bc_level["freq"] +\
   myData_uncorrected_bc_level["psi_freq"]) - np.log2(myData_uncorrected_bc_level["bl_freq"] +\
   myData_uncorrected_bc_level["psi_freq"])
```

---

## 🧩 Code Structure

* `pipeline.py`: Main pipeline script (this repo)
* `calculate_fitness_matrix.py`: Functions to load data, compute psi-freq and fitness
* `integrate_fitness_data.py`: Functions to integrate fitness with barcode data