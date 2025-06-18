# LDLR MPRA Full Dataset Analysis

This directory contains a shell script and all necessary files to run the fiber-seq analyis for the LDLR (Low-Density Lipoprotein Receptor) MPRA (Massively Parallel Reporter Assay) data using the `fiberseq_MPRA` package.

## Overview

The analysis pipeline processes fiber-seq data to identify significant differences between wild-type (WT) and variant sequences in the LDLR gene. 

## Directory Structure

```
test/
├── ldlr_mpra_analysis.sh          # Main analysis script
├── GRCh38_LDLR.tsv                # Reference TSV file
├── fp_WT_bed/
│   └── LDLR_WT_fp.bed             # Wild-type BED file
└── LDLR_pure_SNVs/
    └── fp_beds/                    # Directory containing variant BED files
        ├── variant1_fp.bed
        ├── variant2_fp.bed
        └── ...                     # Additional variant files
```

## Prerequisites

### Software Requirements
- Bash shell environment
- Python 3.7+
- `fiberseq_MPRA` Python package installed
- Access to high-performance computing resources (recommended)

### System Requirements
- **Memory**: 32+ GB RAM recommended
- **CPU**: Multi-core system (script uses up to 128 processes)
- **Storage**: Several GB of free space for output files
- **Time**: 2-8 hours depending on data size and system specs

## Input Files

### Required Files
1. **Wild-type BED file** (`fp_WT_bed/LDLR_WT_fp.bed`)
   - Contains fiber-seq footprint data for the wild-type LDLR sequence
   - Format: FiberHMM footprint BED format

2. **Variant BED files** (`LDLR_pure_SNVs/fp_beds/*_fp.bed`)
   - Individual BED files for each LDLR variant
   - Files must end with `_fp.bed` suffix
   - Same format as wild-type file

3. **Reference TSV file** (`GRCh38_LDLR.tsv`)
   - Contains reference information for LDLR variants (from Kircher et al. https://doi.org/10.1038/s41467-019-11526-w)
   - Used for annotation and filtering

## Usage

### Basic Execution

Navigate to the `test` directory and run:

```bash
cd test/
chmod +x ldlr_mpra_analysis.sh
./ldlr_mpra_analysis.sh
```

### Script Parameters

The script uses the following key parameters (hardcoded, can change as desired, but leave Reference sequence length):

- **Reference sequence length**: 4,718 bp
- **Subsample size**: 5,000 reads
- **Iterations**: 10,000 (for robust statistics)
- **Footprint size range**: 1-400 bp
- **Genomic position range**: 3,000-3,600 bp
- **Bin size**: 10 bp
- **Minimum coverage**: 2,000 reads
- **Processes**: 128 (adjust based on your system)

## Analysis Pipeline

The script runs four main steps:

### Step 1: Control Data Generation
- Generates control distributions from wild-type data
- Creates `control/control.pkl`
- **Runtime**: ~30-60 minutes

### Step 2: Null Distribution Generation
- Calculates null distributions for statistical testing
- Creates `null/null_means.pkl` and `null/null_stds.pkl`
- **Runtime**: ~30-60 minutes

### Step 3: Standard Analysis
- Performs main statistical analysis with standard p-values
- Generates plots and results in `results_standard/`
- **Runtime**: ~1-4 hours

### Step 4: FDR-Corrected Analysis
- Repeats analysis with False Discovery Rate (FDR) correction
- Generates plots and results in `results_fdr/`
- Uses FDR threshold of 0.05
- **Runtime**: ~1-4 hours

## Output Structure

After completion, the following directories will be created:

```
test/
├── fiberseq_MPRA_full_analysis/    # Main analysis directory
│   ├── control/
│   │   └── control.pkl             # Control data
│   ├── null/
│   │   ├── null_means.pkl          # Null distribution means
│   │   └── null_stds.pkl           # Null distribution standard deviations
│   ├── results_standard/           # Standard p-value results
│   │   ├── *.png                   # Heatmap plots
│   │   └── *.csv                   # Statistical results
│   ├── results_fdr/                # FDR-corrected results
│   │   ├── *.png                   # Heatmap plots
│   │   └── *.csv                   # Statistical results
│   └── analysis_log.txt            # Summary log file
```

## Interpreting Results

### Heatmap Plots
- **Color scale**: magma colormap (2-10 range)
- **X-axis**: Genomic positions (50 bp tick marks)
- **Y-axis**: Footprint sizes
- **Intensity**: Statistical significance

### Statistical Files
- CSV files contain p-values and effect sizes for each variant
- Compare standard vs. FDR-corrected results
- Lower p-values indicate more significant differences from WT

## Resuming Analysis

The script automatically detects existing preprocessing files:
- If `control/control.pkl` exists, control generation is skipped
- If null distribution files exist, null generation is skipped
- To force regeneration, delete the respective `.pkl` files

## Citation

If using this analysis pipeline, please cite the relevant publications for:
- https://www.biorxiv.org/content/10.1101/2025.05.27.656122v1
