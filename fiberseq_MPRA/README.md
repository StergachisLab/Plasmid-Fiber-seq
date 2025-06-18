# Footprint Analysis

A modular Python package for SNP footprint analysis with multiprocessing support.

## Overview

This package processes BED files containing SNP data, groups them by position, and creates comparison heatmaps for each position group. It uses multiprocessing to accelerate data processing and plotting.

## Installation

Clone the repository and install the package:

```bash
git clone https://github.com/yourusername/footprint_analysis.git
cd footprint_analysis
pip install -e .
```

## Package Structure

```
footprint_analysis/
├── __init__.py
├── core/
│   ├── __init__.py
│   ├── data_processing.py  # Core data processing functions
│   └── file_io.py          # File reading/writing functions
├── analysis/
│   ├── __init__.py
│   └── statistics.py       # Statistical analysis functions
├── visualization/
│   ├── __init__.py
│   └── plots.py            # Plotting functions
├── utils/
│   ├── __init__.py
│   └── helpers.py          # Helper functions and utilities
└── main.py                 # Command-line interface and main execution
```

## Usage

### Command-line Interface

The package provides a command-line interface for processing BED files and generating heatmaps:

```bash
python -m footprint_analysis.main \
    --input-dir /path/to/bed/files \
    --control-pkl /path/to/control.pkl \
    --null-means-pkl /path/to/null_means.pkl \
    --null-stds-pkl /path/to/null_stds.pkl \
    --output-dir /path/to/output
```

### Required Arguments

- `--input-dir`: Directory containing BED files to process
- `--control-pkl`: Path to the pickled control dataframe
- `--null-means-pkl`: Path to the pickled dictionary with means of null distributions
- `--null-stds-pkl`: Path to the pickled dictionary with standard deviations of null distributions
- `--output-dir`: Directory where the output plots will be saved

### Optional Arguments

- `--metadata-pkl`: Path to the pickled metadata with dimension information
- `--ref-seq-length`: Length of the reference sequence (default: 4718)
- `--row-min`: Minimum row value for the range (default: 1)
- `--row-max`: Maximum row value for the range (default: 400)
- `--col-min`: Minimum column value for the range (default: 3000)
- `--col-max`: Maximum column value for the range (default: 3600)
- `--bin-size`: Bin size for footprint sizes (default: 10)
- `--vmin`: Minimum value for the heatmap color scale (default: 2)
- `--vmax`: Maximum value for the heatmap color scale (default: 10)
- `--xticks`: Interval for x-axis ticks (default: 50)
- `--cmap`: Colormap for the heatmap (default: "magma")
- `--plot-wt`: Plot WT enrichment (default is variant)
- `--plot-both`: Plot both WT and variant enrichment
- `--processes`: Number of processes to use for parallel processing
- `--no-mp`: Disable multiprocessing and use sequential processing
- `--tsv-file`: Path to a TSV file with log2-FC values to display in the plots
- `--min-coverage`: Minimum number of reads required to process a BED file
- `--apply-fdr`: Apply Benjamini-Hochberg FDR correction to p-values
- `--fdr-threshold`: FDR threshold for significance (default: 0.05)

## Dependencies

- numpy
- pandas
- matplotlib
- seaborn
- scipy

## Implementation Notes

### External Dependencies

This package relies on functions from the FiberHMM_functions module for some core functionality. The following functions need to be implemented or imported:

- `grab_circular_reads`
- `prep_dfs_for_subtraction`
- `filter_fp_df`
- `read_ft_data`

These functions are referenced but not fully implemented in this package as they depend on external code.

## License

[MIT License](LICENSE)
