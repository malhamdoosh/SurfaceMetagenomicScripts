# Bioinformatic Scripts for Metagenomics

This is repo is focused on Shotgun Metagenomics for Pathogen Detection on High-Touch Surfaces


## Clinical Validation Assay Analysis & Heatmap Generator

`generate_assay_analysis.py` is a Python tool designed to perform clinical validation by comparing metagenomic sequencing classification results (Kraken2 report files) against standard clinical diagnostic assay results (e.g., QIAstat-Dx assays). 

It generates a numeric presence/absence matrix CSV and visualizes the concordance between metagenomic reads and clinical assay findings using a categorized heatmap (True Positive, True Negative, False Positive, False Negative, and No Data).

---

## Table of Contents
- [Dependencies & Installation](#dependencies--installation)
- [Command Line Arguments](#command-line-arguments)
- [Input Files & Formatting](#input-files--formatting)
- [Execution Usage & Examples](#execution-usage--examples)
- [Outputs & Heatmap Categorization Key](#outputs--heatmap-categorization-key)

---

## Dependencies & Installation

This script requires **Python 3.x** and the following data analysis and visualization libraries:
- `pandas=3.0.0`
- `numpy=2.4.1`
- `seaborn=0.13.2`
- `matplotlib-base=3.10.8`

### Installation
Install all dependencies using `pip` or `conda`:

```bash
pip install pandas numpy seaborn matplotlib
```

*Or via Conda:*
```bash
conda install pandas numpy seaborn matplotlib
```

---

## Command Line Arguments

The script uses Python's standard `argparse` module. All arguments have built-in default values:

| Parameter | Type | Default Value | Description |
| :--- | :--- | :--- | :--- |
| `--targets` | `str` | `target_list.txt` | Path to the target list TSV file |
| `--assay` | `str` | `assay_results.txt` | Path to the clinical assay results TSV file |
| `--reports_dir` | `str` | `.` | Directory containing Kraken report files (`*_report.txt`) |
| `--out_matrix` | `str` | `pathogen_presence_absence_matrix.csv` | Output filename for presence/absence matrix |
| `--out_png` | `str` | `heatmap.png` | Output filename for heatmap PNG image |

---

## Input Files & Formatting

### 1. Target List TSV (`--targets`)
A tab-separated file mapping NCBI Taxonomy IDs from Kraken reports to clinical assay target names.
- **Header**: Skipped by the script (`skiprows=1`), but expects 4 tab-separated columns:
  1. `rank`: Taxonomic rank (e.g. `S`, `G`, `S1`)
  2. `TaxID`: NCBI Taxonomy ID (e.g. `2697049`, `11320`)
  3. `name_kraken_report`: Kraken reference organism name
  4. `target_name_assay`: Clinical assay target name (e.g. `SARS-CoV-2`, `Influenza A`)

**Example `target_list.txt` snippet:**
```tsv
rank	tax_ID	name_kraken_report	target_name_assay
S1	2697049	Severe acute respiratory syndrome coronavirus 2	SARS-CoV-2
S	11320	Influenza A virus	Influenza A
G	12059	Enterovirus	Rhinovirus/Enterovirus
```

### 2. Assay Results TSV (`--assay`)
A tab-separated file detailing positive assay targets for each clinical sample.
- **Format**: 2 tab-separated columns:
  1. `sample_name`: Sample identifier matching the Kraken report prefix (e.g. `CP17`, `HP6`)
  2. `assay_targets`: Comma-separated list of positive target names for that sample (or empty/NA if none)

**Example `assay_results.txt` snippet:**
```tsv
CP18	Adenovirus
AP9	Rhinovirus/Enterovirus,SARS-CoV-2
CP17	SARS-CoV-2,Vibrio cholerae
```

### 3. Kraken Report Directory (`--reports_dir`)
- Directory containing individual Kraken2 report files ending with `_report.txt` (e.g., `CP17_report.txt` or `Sample1_kraken2_trimmed_report.txt`).
- **Sample ID Parsing**: Sample names are automatically parsed from the file prefix before the first underscore `_`.
- **Columns Used**: Column 2 (Read count) and Column 5 (Taxonomy ID).

---

## Execution Usage & Examples

Run the script from your terminal specifying your input paths:

### Basic Execution (Using Defaults):
```bash
python generate_assay_analysis.py
```

### Full Execution with Custom Arguments:
```bash
python generate_assay_analysis.py \
  --targets target_list.txt \
  --assay assay_results_with_only_mgx_samples.txt \
  --reports_dir ./kraken_reports \
  --out_matrix pathogen_presence_absence_matrix.csv \
  --out_png clinical_validation_heatmap.png
```

---

## Outputs & Heatmap Categorization Key

The script produces two main output artifacts:

### 1. Presence/Absence Matrix CSV (`--out_matrix`)
Contains raw read counts mapped from Kraken2 reports for each target organism across all processed samples alongside target metadata (`rank`, `TaxID`, `name_kraken_report`, `target_name_assay`).

### 2. Clinical Validation Heatmap (`--out_png`)
Generates a high-resolution (300 DPI) Seaborn/Matplotlib heatmap plot illustrating concordance between Metagenomic sequencing (Kraken2) and Clinical Assay results:

| Category | Label | Color | Criteria |
| :--- | :--- | :--- | :--- |
| **TN** | True Negative | 🟢 Light Green (`#a8e6cf`) | Negative in Assay AND 0 Metagenomic reads |
| **TP** | True Positive | 🌲 Dark Green (`#2e7d32`) | Positive in Assay AND $\ge 1$ Metagenomic read |
| **FP** | Unexpected Hit | 🟡 Orange (`#ffc107`) | Negative in Assay BUT $\ge 1$ Metagenomic read detected |
| **FN** | Missing Hit | 🔴 Red (`#f44336`) | Positive in Assay BUT 0 Metagenomic reads detected |
| **NA** | No Data | ⚪ Gray (`#d1d1d1`) | Missing sample or missing assay data |
