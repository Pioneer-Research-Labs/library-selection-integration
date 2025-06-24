# 📈 Fitness Calculation & Integration Pipeline

This Python-based pipeline calculates **fitness metrics** for sequencing data across multiple **library**, **strain**, and **condition** combinations for a selection experiment. It integrates fitness results with PacBio library data and performs error correction on barcdoes to the libraries.

---

## 🧬 Overview

This pipeline performs the following core functions:

1. **Loads and merges short-read data and metadata**
2. **Calculates fitness scores**
3. **Integrates calculated fitness with library long read sequencing data**
4. **Corrects barcodes to library barcodes**
5. **Saves both raw and integrated fitness results**

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
python pipeline.py --short_read_path <SHORT_READ_DIR> [--metadata <METADATA_CSV>] [--out_prefix <OUTPUT_PREFIX>] [--base_timepoint <TIMEPOINT>]
```

### Arguments:

| Argument            | Type  | Required | Description                                             |
| ------------------- | ----- | -------- | ------------------------------------------------------- |
| `--short_read_path` | `str` | ✅        | Path to directory containing short-read results (this should usually be `results`).        |
| `--metadata`        | `str` | ❌        | Metadata CSV filename (default: `metadata.csv`).        |
| `--out_prefix`      | `str` | ❌        | Prefix for output files (default: `fitness`).           |
| `--base_timepoint`  | `int` | ❌        | Timepoint used as the fitness reference (default: `0`). |

---

## 📂 Input Files
These are generated from running `short-read-pipeline`.

---

## 📤 Output

For each unique `(library, strain, condition)` group, two output files are generated:

1. `fitness_<library>_<strain>_<condition>.parquet`
   → Raw fitness data

2. `fitness_integrated_<library>_<strain>_<condition>.parquet`
   → Fitness data with library integration and barcode correction

---

## 🧩 Code Structure

* `pipeline.py`: Main pipeline script (this repo)
* `calculate_fitness_matrix.py`: Functions to load data, compute psi-freq and fitness
* `integrate_fitness_data.py`: Functions to integrate fitness with barcode data

