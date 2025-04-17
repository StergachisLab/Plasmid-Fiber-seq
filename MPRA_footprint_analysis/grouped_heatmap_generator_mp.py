#!/usr/bin/env python3
"""
Grouped Heatmap Generator for SNP Footprint Analysis with Multiprocessing

This script processes BED files containing SNP data, groups them by position,
and creates comparison heatmaps for each position group in a single PDF.
Multiprocessing is used to accelerate data processing and plotting.
"""

import os
import sys
import argparse
import warnings
import traceback
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for multiprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle
import numpy as np
sys.path.append('/gscratch/stergachislab/bmallo/large_home/python_scripts/')
from FiberHMM_functions import *
from scipy import stats
from scipy.stats import rankdata
from matplotlib.colors import Normalize
import multiprocessing as mp
from functools import partial
import time

def parse_bed_filename(filename):
    """
    Parse a BED filename to extract SNP position and base change information.
    Handles both single SNP format (3188_G_T_single_fp.bed) and 
    merged format (merged_3239_3242_fp.bed).
    
    Parameters:
    -----------
    filename : str
        Filename in format like "3244_G_A_fp.bed" or "merged_3239_3242_fp.bed"
        
    Returns:
    --------
    tuple
        For single SNP: (original_position, adjusted_position, complemented_ref_base, complemented_var_base, is_merged, second_position, second_adjusted)
        For merged: (first_position, first_adjusted, None, None, True, second_position, second_adjusted)
    """
    basename = os.path.basename(filename)
    
    # Check if this is a merged file
    if basename.startswith("merged_"):
        try:
            # Extract positions from merged filename
            parts = basename.replace("merged_", "").split("_")
            first_position = int(parts[0])
            second_position = int(parts[1])
            
            # Calculate adjusted positions
            first_adjusted = 11092732 - first_position
            second_adjusted = 11092732 - second_position
            
            # Return tuple indicating merged format
            return first_position, first_adjusted, None, None, True, second_position, second_adjusted
        except (IndexError, ValueError):
            return None, None, None, None, False, None, None
    
    # Original format processing
    try:
        parts = basename.split('_')
        original_position = int(parts[0])
        original_ref_base = parts[1]
        original_var_base = parts[2]
        
        # Calculate the adjusted position
        adjusted_position = 11092732 - original_position
        
        # Complement the bases
        base_complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        complemented_ref_base = base_complements.get(original_ref_base)
        complemented_var_base = base_complements.get(original_var_base)
        
        # Return tuple indicating regular format
        return original_position, adjusted_position, complemented_ref_base, complemented_var_base, False, None, None
    except (IndexError, ValueError):
        return None, None, None, None, False, None, None
    
def count_bed_file_rows(bed_file_path):
    """
    Count the number of rows (reads) in a BED file.
    
    Parameters:
    -----------
    bed_file_path : str
        Path to the BED file
        
    Returns:
    --------
    int
        Number of reads (rows) in the BED file
    """
    try:
        with open(bed_file_path, 'r') as file:
            # Count the lines in the file
            read_count = sum(1 for _ in file)
            return read_count
    except Exception as e:
        print(f"Error counting reads in {bed_file_path}: {str(e)}")
        return 0

def format_sample_name(sample_name):
    """
    Format a sample name like "3191_G_T" to show the complemented bases as "C→T".
    
    Parameters:
    -----------
    sample_name : str
        The sample name in format like "3191_G_T"
        
    Returns:
    --------
    str
        The formatted name showing complemented bases like "C→T"
    """
    try:
        # Check if the sample name follows the expected pattern
        parts = sample_name.split('_')
        if len(parts) >= 3:
            # Extract the reference and variant bases
            ref_base = parts[1].upper()
            var_base = parts[2].upper()
            
            # Define complement mapping
            base_complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            
            # Get complements
            complemented_ref = base_complements.get(ref_base, ref_base)
            complemented_var = base_complements.get(var_base, var_base)
            
            # Format as "C→T"
            return f"{complemented_ref}→{complemented_var}"
        else:
            # If the sample name doesn't follow the expected pattern, return as is
            return sample_name
    except Exception as e:
        print(f"Error formatting sample name {sample_name}: {str(e)}")
        return sample_name

# =============================================================================
# Load analysis dictionaries and process data
# =============================================================================

def load_analysis_dictionaries(control_pickle_path, null_means_dict_path, null_stds_dict_path, metadata_path=None):
    """
    Load control dataframe and precomputed null statistics dictionaries.
    
    Parameters:
    -----------
    control_pickle_path : str
        Path to the pickled control dataframe
    null_means_dict_path : str
        Path to the pickled dictionary with means of null distributions
    null_stds_dict_path : str
        Path to the pickled dictionary with standard deviations of null distributions
    metadata_path : str, optional
        Path to the pickled metadata with dimension information
        
    Returns:
    --------
    tuple of (control_df, null_means_dict, null_stds_dict) objects
    """
    # Load the control dataframe
    try:
        print(f"Loading control dataframe from {control_pickle_path}")
        control_df = pd.read_pickle(control_pickle_path)
    except Exception as e:
        raise Exception(f"Failed to load control dataframe from {control_pickle_path}: {e}")
    
    # Load precomputed null distribution statistics dictionaries
    try:
        print(f"Loading precomputed null distribution means dictionary from {null_means_dict_path}")
        with open(null_means_dict_path, 'rb') as f:
            null_means_dict = pickle.load(f)
        
        print(f"Loading precomputed null distribution standard deviations dictionary from {null_stds_dict_path}")
        with open(null_stds_dict_path, 'rb') as f:
            null_stds_dict = pickle.load(f)
        
        # Load metadata if provided
        if metadata_path:
            print(f"Loading metadata from {metadata_path}")
            with open(metadata_path, 'rb') as f:
                metadata = pickle.load(f)
            
            # Check if the dictionaries have the same dimensions as the control dataframe
            if (metadata['shape'][0] != control_df.shape[0] or 
                metadata['shape'][1] != control_df.shape[1]):
                warnings.warn("Precomputed null statistics have different dimensions than control dataframe. This may cause issues.")
    except Exception as e:
        raise Exception(f"Failed to load precomputed statistics: {e}")
    
    print("All data loaded successfully")
    return control_df, null_means_dict, null_stds_dict

def process_var_df(input_df, ref_seq_length, row_range, col_range, bin_size):
    """
    Process a single sample iteration.

    Parameters:
    -----------
    input_df : pandas.DataFrame
        The bed DataFrame
    ref_seq_length : int
        Length of the reference sequence
    row_range : tuple or None
        Range of rows to include
    col_range : tuple or None
        Range of columns to include
    bin_size : int or None
        Bin size for footprint sizes

    Returns:
    --------
    The filtered footprint count DataFrame for this iteration
    """
    # Get number of reads in the sample
    read_count = len(input_df)

    # Process bed file into expanded footprint dataframe
    circular_footprint_df = grab_circular_reads(input_df, ref_seq_length)

    # Transform expanded footprint dataframe into count of footprint sizes dataframe
    footprint_count_dfs = prep_dfs_for_subtraction(circular_footprint_df)

    # Filter footprint_count_dfs for desired footprint sizes and positions
    # Keep using the original filter_fp_df function
    filtered_footprint_count_df = filter_fp_df(footprint_count_dfs[0], bin_size, row_range, col_range)
    
    # Just apply a numeric sort to fix the bin ordering issues
    if bin_size and filtered_footprint_count_df.shape[0] > 0:
        try:
            # Extract the numeric starting value from each bin label
            bin_starts = [int(idx.split('-')[0]) for idx in filtered_footprint_count_df.index]
            # Sort the dataframe by these numeric starting values
            filtered_footprint_count_df = filtered_footprint_count_df.iloc[np.argsort(bin_starts)]
        except Exception as e:
            print(f"Warning: Could not sort bin indices: {e}")
    
    # Normalize by read count and return
    return filtered_footprint_count_df / read_count

def apply_multiple_testing_correction(p_value_df, fdr_threshold=0.05):
    """
    Apply Benjamini-Hochberg multiple testing correction to p-values.
    Only non-zero p-values will be considered in the correction.
    Values that don't meet the FDR threshold will be set to 0.
    
    Parameters:
    -----------
    p_value_df : pandas.DataFrame
        DataFrame with raw p-values (not -log10 transformed)
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with adjusted -log10(p-values), values not meeting threshold set to 0
    """
    # Extract non-zero p-values
    non_zero_pvalues = []
    for row in p_value_df.index:
        for col in p_value_df.columns:
            p_val = p_value_df.loc[row, col]
            if p_val < 1.0:  # Only include non-1.0 p-values (which correspond to positive enrichment)
                non_zero_pvalues.append((row, col, p_val))
    
    if len(non_zero_pvalues) == 0:
        return p_value_df.copy() * 0  # No p-values to correct, return zeros
    
    # Sort by p-value
    non_zero_pvalues.sort(key=lambda x: x[2])
    
    # Get the total number of tests
    n = len(non_zero_pvalues)
    
    # Calculate adjusted p-values
    adjusted_pvalues = {}
    for i, (row, col, p_val) in enumerate(non_zero_pvalues):
        rank = i + 1
        adjusted_p = min(p_val * n / rank, 1.0)
        adjusted_pvalues[(row, col)] = adjusted_p
    
    # Apply cumulative minimum from highest to lowest rank
    for i in range(n-2, -1, -1):
        row, col, _ = non_zero_pvalues[i]
        next_row, next_col, _ = non_zero_pvalues[i+1]
        adjusted_pvalues[(row, col)] = min(
            adjusted_pvalues[(row, col)],
            adjusted_pvalues[(next_row, next_col)]
        )
    
    # Create output DataFrame with -log10 transformation
    result_df = p_value_df.copy() * 0  # Initialize with zeros
    for (row, col), adj_p in adjusted_pvalues.items():
        # Only include values that meet the threshold
        if adj_p <= fdr_threshold:
            # Use a small value to avoid log10(0)
            min_p_value = 1e-16
            adj_p = max(adj_p, min_p_value)
            result_df.loc[row, col] = -np.log10(adj_p)
        # Values not meeting threshold remain 0
    
    return result_df

def calculate_pvalues_one_tailed_dict(enrichment_df, null_means_dict, null_stds_dict, 
                                     apply_fdr=False, fdr_threshold=0.05):
    """
    Calculate -log10(p-value) for each cell in the enrichment dataframe based on precomputed null statistics.
    Only calculates p-values for positive enrichment (values > 0).
    Uses dictionary lookups with safety checks for faster performance.
    Optionally applies FDR correction.
    
    Parameters:
    -----------
    enrichment_df : pandas.DataFrame
        Dataframe with enrichment values
    null_means_dict : dict
        Dictionary with (row, col) keys and mean values
    null_stds_dict : dict
        Dictionary with (row, col) keys and standard deviation values
    apply_fdr : bool
        Whether to apply FDR correction to the p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    pandas.DataFrame
        Dataframe with -log10(p-value) for each cell
    """
    # Initialize the raw p-value dataframe with ones (representing p-value = 1)
    raw_pvalues_df = pd.DataFrame(1.0, index=enrichment_df.index, columns=enrichment_df.columns)
    
    # Iterate only over cells with positive enrichment
    for row in enrichment_df.index:
        for col in enrichment_df.columns:
            # Get the observed value
            observed_value = enrichment_df.loc[row, col]
            
            # Skip cells with non-positive enrichment
            if observed_value <= 0:
                continue
            
            # Check if the key exists in the dictionaries
            key = (row, col)
            if key not in null_means_dict or key not in null_stds_dict:
                print(f"Warning: Key {key} not found in null statistics dictionaries")
                continue
                
            # Get the null statistics from dictionaries with faster lookup
            null_mean = null_means_dict[key]
            null_std = null_stds_dict[key]
            
            # Check if std is 0 to avoid division by zero
            if null_std == 0:
                # If std is 0, the distribution is a constant
                if observed_value > null_mean:
                    # If observed is greater than the constant value, p-value approaches 0
                    raw_pvalues_df.loc[row, col] = 1e-16
                continue
            
            # Calculate z-score
            z_score = (observed_value - null_mean) / null_std
            
            # Calculate one-tailed p-value (upper tail)
            p_value = 1 - stats.norm.cdf(z_score)
            
            # Apply floor to avoid computational issues
            min_p_value = 1e-16
            p_value = max(p_value, min_p_value)
            
            # Store the raw p-value
            raw_pvalues_df.loc[row, col] = p_value
    
    # Apply FDR correction if requested
    if apply_fdr:
        result_df = apply_multiple_testing_correction(raw_pvalues_df, fdr_threshold)
    else:
        # Convert raw p-values to -log10 scale
        result_df = pd.DataFrame(0, index=raw_pvalues_df.index, columns=raw_pvalues_df.columns)
        for row in raw_pvalues_df.index:
            for col in raw_pvalues_df.columns:
                p_value = raw_pvalues_df.loc[row, col]
                if p_value < 1.0:  # Skip p-value of 1.0 (non-significant)
                    min_p_value = 1e-16
                    p_value = max(p_value, min_p_value)
                    result_df.loc[row, col] = -np.log10(p_value)
    
    return result_df

# Function to safely compare samples without multiprocessing
def compare_samples_sequential(sample_directory_path, control_df, null_means_dict, null_stds_dict,
                               ref_seq_length, row_range=None, col_range=None, bin_size=None,
                               apply_fdr=False, fdr_threshold=0.05):
    """
    Compare variant and control samples using sequential processing instead of multiprocessing.
    This is a fallback function to use if multiprocessing encounters errors.
    
    Parameters:
    -----------
    Same as compare_samples_with_dictionaries
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue) DataFrames
    """
    # Setup dictionary to store processed dataframes
    pvalue_df_dict = {}
    
    # Read and process the sample dataframes
    bed_dict = read_ft_data(sample_directory_path)
    total_samples = len(bed_dict)
    
    print(f"Processing {total_samples} samples sequentially...")
    for i, (sample_name, input_df) in enumerate(bed_dict.items(), 1):
        print(f"Processing sample {i}/{total_samples}: {sample_name}")
        try:
            # Process the sample
            processed_count_df = process_var_df(input_df, ref_seq_length, row_range, col_range, bin_size)

            # Create enrichment dataframes
            WT_enrichment = control_df - processed_count_df  # Positive where WT is enriched
            Var_enrichment = processed_count_df - control_df  # Positive where Var is enriched

            # Convert enrichment values to -log10(p-values)
            print(f"  Calculating p-values for {sample_name}...")
            if apply_fdr:
                print(f"  Applying FDR correction with threshold {fdr_threshold}")
                
            WT_pvalue = calculate_pvalues_one_tailed_dict(WT_enrichment, null_means_dict, null_stds_dict,
                                                        apply_fdr=apply_fdr, fdr_threshold=fdr_threshold)
            Var_pvalue = calculate_pvalues_one_tailed_dict(Var_enrichment, null_means_dict, null_stds_dict,
                                                         apply_fdr=apply_fdr, fdr_threshold=fdr_threshold)

            # Store the p-value dataframes
            pvalue_df_dict[sample_name] = (WT_pvalue, Var_pvalue)
            print(f"  Finished processing {sample_name}")
        except Exception as e:
            print(f"Error processing sample {sample_name}: {str(e)}")
            print(traceback.format_exc())
    
    print(f"Successfully processed {len(pvalue_df_dict)} out of {total_samples} samples")
    return pvalue_df_dict

# Function to process a single sample for multiprocessing
def process_single_sample(sample_data, control_df, null_means_dict, null_stds_dict, 
                          ref_seq_length, row_range, col_range, bin_size,
                          apply_fdr=False, fdr_threshold=0.05):
    """
    Process a single sample for multiprocessing.
    
    Parameters:
    -----------
    sample_data : tuple
        Tuple of (sample_name, dataframe)
    control_df : pandas.DataFrame
        Control dataframe
    null_means_dict : dict
        Dictionary with means of null distributions
    null_stds_dict : dict
        Dictionary with standard deviations of null distributions
    ref_seq_length : int
        Length of the reference sequence
    row_range : tuple or None
        Range of rows to include
    col_range : tuple or None
        Range of columns to include
    bin_size : int or None
        Bin size for footprint sizes
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    tuple
        (sample_name, (WT_pvalue, Var_pvalue))
    """
    sample_name, input_df = sample_data
    print(f"Processing sample: {sample_name}")
    
    try:
        # Process the sample
        processed_count_df = process_var_df(input_df, ref_seq_length, row_range, col_range, bin_size)

        # Create enrichment dataframes
        WT_enrichment = control_df - processed_count_df  # Positive where WT is enriched
        Var_enrichment = processed_count_df - control_df  # Positive where Var is enriched

        # Convert enrichment values to -log10(p-values)
        print(f"  Calculating p-values for {sample_name}...")
        if apply_fdr:
            print(f"  Applying FDR correction with threshold {fdr_threshold}")
            
        WT_pvalue = calculate_pvalues_one_tailed_dict(WT_enrichment, null_means_dict, null_stds_dict, 
                                                     apply_fdr=apply_fdr, fdr_threshold=fdr_threshold)
        Var_pvalue = calculate_pvalues_one_tailed_dict(Var_enrichment, null_means_dict, null_stds_dict,
                                                      apply_fdr=apply_fdr, fdr_threshold=fdr_threshold)

        print(f"  Finished processing {sample_name}")
        return (sample_name, (WT_pvalue, Var_pvalue))
    
    except Exception as e:
        print(f"Error processing sample {sample_name}: {str(e)}")
        print(traceback.format_exc())
        return (sample_name, None)

def compare_samples_with_dictionaries(sample_directory_path, control_df, null_means_dict, null_stds_dict,
                                     ref_seq_length, row_range=None, col_range=None, bin_size=None,
                                     n_processes=None, use_multiprocessing=True,
                                     apply_fdr=False, fdr_threshold=0.05):
    """
    Compare variant and control samples using preloaded control dataframe and null statistics dictionaries.
    Uses parallel processing for faster analysis if multiprocessing is enabled and available.
    
    Parameters:
    -----------
    sample_directory_path : str
        Path to the directory containing BED files
    control_df : pandas.DataFrame
        Preloaded control dataframe
    null_means_dict : dict
        Dictionary with (row, col) keys and mean values
    null_stds_dict : dict
        Dictionary with (row, col) keys and standard deviation values
    ref_seq_length : int
        Length of the reference sequence
    row_range : tuple or None
        Range of rows to include
    col_range : tuple or None
        Range of columns to include
    bin_size : int or None
        Bin size for footprint sizes
    n_processes : int or None
        Number of processes to use for parallel processing
    use_multiprocessing : bool
        Whether to use multiprocessing or fall back to sequential processing
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue) DataFrames
    """
    # Read the sample dataframes first
    bed_dict = read_ft_data(sample_directory_path)
    total_samples = len(bed_dict)
    
    if total_samples == 0:
        print("No samples found. Check the input directory.")
        return {}
    
    # If multiprocessing is not requested or only 1 sample, use sequential processing
    if not use_multiprocessing or total_samples <= 1:
        return compare_samples_sequential(
            sample_directory_path, control_df, null_means_dict, null_stds_dict,
            ref_seq_length, row_range, col_range, bin_size,
            apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
        )
    
    # Determine number of processes to use
    if n_processes is None:
        n_processes = max(1, min(mp.cpu_count() - 1, total_samples))
    else:
        n_processes = max(1, min(n_processes, mp.cpu_count(), total_samples))
    
    print(f"Using {n_processes} processes for parallel sample processing")
    
    # Setup dictionary to store processed dataframes
    pvalue_df_dict = {}
    
    # Create a partial function with fixed arguments
    process_func = partial(
        process_single_sample,
        control_df=control_df,
        null_means_dict=null_means_dict,
        null_stds_dict=null_stds_dict,
        ref_seq_length=ref_seq_length,
        row_range=row_range,
        col_range=col_range,
        bin_size=bin_size,
        apply_fdr=apply_fdr,
        fdr_threshold=fdr_threshold
    )
    
    try:
        # Process samples in parallel
        print(f"Processing {total_samples} samples in parallel...")
        start_time = time.time()
        
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(process_func, bed_dict.items())
        
        # Store results in the dictionary
        for sample_name, pvalue_tuple in results:
            if pvalue_tuple is not None:
                pvalue_df_dict[sample_name] = pvalue_tuple
        
        end_time = time.time()
        print(f"All samples processed in {end_time - start_time:.2f} seconds")
        print(f"Successfully processed {len(pvalue_df_dict)} out of {total_samples} samples")
        
    except Exception as e:
        print(f"Error in parallel processing: {str(e)}")
        print(traceback.format_exc())
        print("Falling back to sequential processing...")
        
        return compare_samples_sequential(
            sample_directory_path, control_df, null_means_dict, null_stds_dict,
            ref_seq_length, row_range, col_range, bin_size,
            apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
        )
    
    # If no samples were processed successfully with multiprocessing, fall back to sequential
    if len(pvalue_df_dict) == 0:
        print("No samples were processed successfully with multiprocessing. Falling back to sequential processing...")
        return compare_samples_sequential(
            sample_directory_path, control_df, null_means_dict, null_stds_dict,
            ref_seq_length, row_range, col_range, bin_size,
            apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
        )
    
    return pvalue_df_dict

def compare_samples(sample_directory_path, control_pickle_path, 
                   null_means_dict_path, null_stds_dict_path, metadata_path=None,
                   ref_seq_length=4718, row_range=None, col_range=None, bin_size=None,
                   n_processes=None, use_multiprocessing=True,
                   apply_fdr=False, fdr_threshold=0.05):
    """
    Compare variant and control samples by processing them and creating p-value dataframes.
    Uses multiprocessing for faster performance if enabled.
    
    Parameters:
    -----------
    sample_directory_path : str
        Path to the directory containing BED files
    control_pickle_path : str
        Path to a pickled control dataframe
    null_means_dict_path : str
        Path to a pickled dictionary with means of null distributions
    null_stds_dict_path : str
        Path to a pickled dictionary with standard deviations of null distributions
    metadata_path : str, optional
        Path to the pickled metadata with dimension information
    ref_seq_length : int
        Length of the reference sequence
    row_range : tuple or None
        Range of rows to include
    col_range : tuple or None
        Range of columns to include
    bin_size : int or None
        Bin size for footprint sizes
    n_processes : int or None
        Number of processes to use for parallel processing
    use_multiprocessing : bool
        Whether to use multiprocessing or fall back to sequential processing
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue) DataFrames
    """
    # Load all required data
    control_df, null_means_dict, null_stds_dict = load_analysis_dictionaries(
        control_pickle_path, null_means_dict_path, null_stds_dict_path, metadata_path
    )
    
    # Process using the loaded data
    return compare_samples_with_dictionaries(
        sample_directory_path, control_df, null_means_dict, null_stds_dict,
        ref_seq_length, row_range, col_range, bin_size, n_processes, use_multiprocessing,
        apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
    )

# =============================================================================
# Plot functions
# =============================================================================

def plot_fp_heatmap(df, vmin=None, vmax=None, cmap="viridis", title="", xticks=0, 
                   output_directory=None, legend_label=None, snp_position=None, 
                   snp_ref_base=None, snp_var_base=None, read_count=None):
    """
    Creates a heatmap for footprint enrichment data with an optional vertical line indicating SNP location.

    Parameters:
    df (pd.DataFrame): Data for the heatmap, where rows represent footprint sizes and columns represent base positions.
    vmin (float): Minimum value for the heatmap color scale.
    vmax (float): Maximum value for the heatmap color scale.
    cmap (str): Colormap for the heatmap. Default is "viridis".
    title (str): Title for the plot.
    xticks (int): Interval for x-axis ticks.
    output_directory (str): Directory where the plot will be saved.
    legend_label (str): Label for the SNP location in the legend.
    snp_position (int): Position of the SNP on the x-axis. If None, no SNP line will be drawn.
    snp_ref_base (str): Reference base at the SNP position.
    snp_var_base (str): Variant base at the SNP position.
    read_count (int): Number of reads for this variant (added parameter).
    """

    # Set up figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(df, ax=ax, cbar=True, annot=False, fmt='.2f', vmin=vmin, vmax=vmax, cmap=cmap, xticklabels=xticks)

    # Customize colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_label("-log10(p-value)", fontsize=14)

    # Create legend items for both SNP location and read count
    legend_items = []
    legend_labels = []
    
    # Add vertical line for SNP location if provided
    if snp_position is not None and isinstance(snp_position, (int, float)):
        # Create a descriptive SNP label if base information is available
        snp_description = legend_label
        if not snp_description and snp_ref_base and snp_var_base:
            snp_description = f"{snp_ref_base}→{snp_var_base}"
        
        # Add the SNP indicator line
        snp_line = ax.axvline(x=snp_position, color='white', linestyle='--', linewidth=1, alpha=0.8)
        
        # Add to legend items
        if snp_description:
            legend_items.append(snp_line)
            legend_labels.append(snp_description)
    
    # Add read count to legend if available
    if read_count is not None:
        # Create a blank line with no visibility for the read count entry
        # Using a small line that's effectively invisible
        read_count_line = ax.plot([], [], ' ', alpha=0)[0]
        legend_items.append(read_count_line)
        legend_labels.append(f"Reads: {read_count:,}")
    
    # Add the legend if we have items to show
    if legend_items:
        legend = ax.legend(legend_items, legend_labels, 
                         loc='upper right', bbox_to_anchor=(1.2, 1.12),
                         facecolor='dimgray', edgecolor='black')
        for text in legend.get_texts():
            text.set_color("white")

    # Set title and labels (without modifying title to include read count)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_ylabel("Footprint Size", fontsize=14)

    # Save figure if output directory is provided
    if output_directory:
        # Create the output directory if it doesn't exist
        try:
            os.makedirs(output_directory, exist_ok=True)
            print(f"Output directory confirmed: {output_directory}")
            
            # Clean up the title for use in a filename
            safe_title = "".join([c if c.isalnum() or c in [' ', '_', '-'] else '_' for c in title])
            safe_title = safe_title.replace(' ', '_')
            
            # Save the figure with explicit path
            save_path = os.path.join(output_directory, f"{safe_title}.png")
            print(f"Attempting to save to: {save_path}")
            
            plt.savefig(save_path, format='png', bbox_inches='tight', dpi=300)
            print(f"Figure successfully saved to {save_path}")
            
            # Verify file exists
            if os.path.exists(save_path):
                print(f"Verified file exists: {save_path}")
            else:
                print(f"WARNING: File was not created at {save_path}")
                
        except Exception as e:
            print(f"Error saving figure: {str(e)}")
            # Try alternative approach with different path
            try:
                alternative_path = os.path.join(os.getcwd(), f"{safe_title}.png")
                print(f"Trying alternative save path: {alternative_path}")
                plt.savefig(alternative_path, format='png', dpi=300)
                print(f"Figure saved to alternative path: {alternative_path}")
            except Exception as e2:
                print(f"Alternative save also failed: {str(e2)}")

    # Return the figure and axis for further customization if needed
    return fig, ax

# Function to plot a single position group for multiprocessing
def plot_position_group(pos_data, output_directory, vmin, vmax, cmap, xticks, plot_both, plot_index, apply_fdr=False, fdr_threshold=0.05):
    """
    Plot a single position group for multiprocessing.
    
    Parameters:
    -----------
    pos_data : tuple
        Tuple of (position, sample_names, sample_info)
    output_directory : str
        Directory where the output plots will be saved
    vmin : float
        Minimum value for the heatmap color scale
    vmax : float
        Maximum value for the heatmap color scale
    cmap : str
        Colormap for the heatmap
    xticks : int
        Interval for x-axis ticks
    plot_both : bool
        If True, plot both WT and variant enrichment
    plot_index : int
        Index of the dataframe to plot (0 for WT, 1 for Variant)
    apply_fdr : bool
        Whether FDR correction was applied
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    tuple
        (position, success)
    """
    """
    Plot a single position group for multiprocessing.
    
    Parameters:
    -----------
    pos_data : tuple
        Tuple of (position, sample_names, sample_info)
    output_directory : str
        Directory where the output plots will be saved
    vmin : float
        Minimum value for the heatmap color scale
    vmax : float
        Maximum value for the heatmap color scale
    cmap : str
        Colormap for the heatmap
    xticks : int
        Interval for x-axis ticks
    plot_both : bool
        If True, plot both WT and variant enrichment
    plot_index : int
        Index of the dataframe to plot (0 for WT, 1 for Variant)
        
    Returns:
    --------
    tuple
        (position, success)
    """
    pos, sample_names, sample_info = pos_data
    
    try:
        print(f"Creating grouped plot for position {pos} with {len(sample_names)} samples")
        
        # Determine plot type and count
        n_plots = len(sample_names)
        if plot_both:
            n_plots *= 2  # Double if plotting both WT and variant
        
        # Calculate grid dimensions (try to make it somewhat square)
        n_cols = min(3, n_plots)  # Maximum 3 columns to leave more room for legends
        n_rows = (n_plots + n_cols - 1) // n_cols  # Ceiling division
        
        # Create figure with subplots - adjusted width to provide space for legends
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        
        # Flatten axes array for easier indexing
        if n_rows == 1 and n_cols == 1:
            axes = [axes]
        elif n_rows == 1 or n_cols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        # Keep track of plot index
        plot_idx = 0
        
        # Extract the original (plasmid) position for display
        sample_data = sample_info[sample_names[0]] if sample_names else None
        original_pos = sample_data.get("original_pos") if sample_data else None
        
        # Add a suptitle for the entire figure - including both positions and FDR info if applicable
        if apply_fdr:
            if original_pos:
                fig.suptitle(f"Footprint Enrichment at Position {pos} (Plasmid: {original_pos})\nFDR < {fdr_threshold}", 
                            fontsize=16, y=1.06)
            else:
                fig.suptitle(f"Footprint Enrichment at Position {pos}\nFDR < {fdr_threshold}", 
                            fontsize=16, y=1.06)
        else:
            if original_pos:
                fig.suptitle(f"Footprint Enrichment at Position {pos} (Plasmid: {original_pos})", 
                            fontsize=16, y=1.06)
            else:
                fig.suptitle(f"Footprint Enrichment at Position {pos}", 
                            fontsize=16, y=1.06)
        
        # Create a single colorbar for the entire figure
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)
        sm.set_array([])
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("-log10(p-value)", fontsize=14)
        
        # Plot each sample
        for sample_name in sample_names:
            sample_data = sample_info[sample_name]
            
            if plot_both:
                # Plot WT enrichment
                if plot_idx < len(axes):
                    ax = axes[plot_idx]
                    sns.heatmap(sample_data["wt_df"], ax=ax, cbar=False, 
                               annot=False, fmt='.2f', vmin=vmin, vmax=vmax, 
                               cmap=cmap, xticklabels=xticks)
                    
                    # Format the sample name to show complemented bases
                    formatted_name = format_sample_name(sample_name)
                    
                    # Create legend items for both SNP location, read count, and log2-FC value
                    legend_items = []
                    legend_labels = []
                    
                    # Add vertical line for SNP position if available
                    if sample_data["original_pos_plot"] is not None:
                        snp_line = ax.axvline(x=sample_data["original_pos_plot"], color='white', 
                                             linestyle='--', linewidth=1, alpha=0.8)
                        legend_items.append(snp_line)
                        legend_labels.append("SNP location")
                    
                    # Add read count to legend
                    if sample_data["read_count"] is not None:
                        # Create a blank line with no visibility for the read count entry
                        read_count_line = ax.plot([], [], ' ', alpha=0)[0]
                        legend_items.append(read_count_line)
                        legend_labels.append(f"Reads: {sample_data['read_count']:,}")
                    
                    # Add log2-FC value and P-value to legend if available
                    if sample_data.get("log2_fc_value") is not None:
                        # Create a blank line with no visibility for the log2-FC entry
                        log2_fc_line = ax.plot([], [], ' ', alpha=0)[0]
                        legend_items.append(log2_fc_line)
                        
                        # Format log2-FC value
                        log2_fc = sample_data["log2_fc_value"]
                        if isinstance(log2_fc, (int, float)):
                            legend_labels.append(f"log2-FC: {log2_fc:.2f}")
                            print(f"  Adding log2-FC to legend: {log2_fc:.2f}")
                        else:
                            legend_labels.append(f"log2-FC: {log2_fc}")
                            print(f"  Adding log2-FC to legend: {log2_fc} (not a number)")
                        
                        # Add P-value if available
                        if sample_data.get("pvalue") is not None:
                            pvalue_line = ax.plot([], [], ' ', alpha=0)[0]
                            legend_items.append(pvalue_line)
                            
                            # Format P-value
                            pvalue = sample_data["pvalue"]
                            if isinstance(pvalue, (int, float)):
                                # Display the original p-value without rounding
                                legend_labels.append(f"P-value: {pvalue}")
                                print(f"  Adding P-value to legend: {pvalue}")
                            else:
                                legend_labels.append(f"P-value: {pvalue}")
                                print(f"  Adding P-value to legend: {pvalue} (not a number)")
                    else:
                        print(f"  No log2-FC value found for sample (WT): {sample_name}")
                    
                    # Add the legend if we have items to show
                    if legend_items:
                        legend = ax.legend(legend_items, legend_labels, 
                                         loc='lower right', bbox_to_anchor=(1.05, 1.25),
                                         fontsize=8, framealpha=0.7,
                                         title_fontsize=9,
                                         facecolor='dimgray', edgecolor='black')
                        for text in legend.get_texts():
                            text.set_color("white")
                    
                    # Set title for subplot - now only showing WT status and base
                    ax.set_title(f"{formatted_name}\nWT", fontsize=10)
                    
                    # Only show y-label on leftmost plots
                    if plot_idx % n_cols == 0:
                        ax.set_ylabel("Footprint Size", fontsize=12)
                    else:
                        ax.set_ylabel("")
                    
                    plot_idx += 1
                
                # Plot Variant enrichment
                if plot_idx < len(axes):
                    ax = axes[plot_idx]
                    sns.heatmap(sample_data["var_df"], ax=ax, cbar=False, 
                               annot=False, fmt='.2f', vmin=vmin, vmax=vmax, 
                               cmap=cmap, xticklabels=xticks)
                    
                    # Format the sample name to show complemented bases
                    formatted_name = format_sample_name(sample_name)
                    
                    # Create legend items for both SNP location, read count, and log2-FC value
                    legend_items = []
                    legend_labels = []
                    
                    # Add vertical line for SNP position if available
                    if sample_data["original_pos_plot"] is not None:
                        snp_line = ax.axvline(x=sample_data["original_pos_plot"], color='white', 
                                             linestyle='--', linewidth=1, alpha=0.8)
                        legend_items.append(snp_line)
                        legend_labels.append("SNP location")
                    
                    # Add read count to legend
                    if sample_data["read_count"] is not None:
                        # Create a blank line with no visibility for the read count entry
                        read_count_line = ax.plot([], [], ' ', alpha=0)[0]
                        legend_items.append(read_count_line)
                        legend_labels.append(f"Reads: {sample_data['read_count']:,}")
                    
                    # Add log2-FC value and P-value to legend if available
                    if sample_data.get("log2_fc_value") is not None:
                        # Create a blank line with no visibility for the log2-FC entry
                        log2_fc_line = ax.plot([], [], ' ', alpha=0)[0]
                        legend_items.append(log2_fc_line)
                        
                        # Format log2-FC value
                        log2_fc = sample_data["log2_fc_value"]
                        if isinstance(log2_fc, (int, float)):
                            legend_labels.append(f"log2-FC: {log2_fc:.2f}")
                            print(f"  Adding log2-FC to legend: {log2_fc:.2f}")
                        else:
                            legend_labels.append(f"log2-FC: {log2_fc}")
                            print(f"  Adding log2-FC to legend: {log2_fc} (not a number)")
                        
                        # Add P-value if available
                        if sample_data.get("pvalue") is not None:
                            pvalue_line = ax.plot([], [], ' ', alpha=0)[0]
                            legend_items.append(pvalue_line)
                            
                            # Format P-value
                            pvalue = sample_data["pvalue"]
                            if isinstance(pvalue, (int, float)):
                                # Display the original p-value without rounding
                                legend_labels.append(f"P-value: {pvalue}")
                                print(f"  Adding P-value to legend: {pvalue}")
                            else:
                                legend_labels.append(f"P-value: {pvalue}")
                                print(f"  Adding P-value to legend: {pvalue} (not a number)")
                    else:
                        print(f"  No log2-FC value found for sample (Var): {sample_name}")
                    
                    # Add the legend if we have items to show
                    if legend_items:
                        legend = ax.legend(legend_items, legend_labels, 
                                         loc='lower right', bbox_to_anchor=(1.05, 1.25),
                                         fontsize=8, framealpha=0.7,
                                         title_fontsize=9,
                                         facecolor='dimgray', edgecolor='black')
                        for text in legend.get_texts():
                            text.set_color("white")
                    
                    # Set title for subplot - now only showing Variant status and base
                    ax.set_title(f"{formatted_name}\nVariant", fontsize=10)
                    
                    # Only show y-label on leftmost plots
                    if plot_idx % n_cols == 0:
                        ax.set_ylabel("Footprint Size", fontsize=12)
                    else:
                        ax.set_ylabel("")
                    
                    plot_idx += 1
            else:
                # Select which dataframe to plot based on plot_index
                df_to_plot = sample_data["wt_df"] if plot_index == 0 else sample_data["var_df"]
                base_type = sample_data["ref_base"] if plot_index == 0 else sample_data["var_base"]
                enrichment_type = "WT" if plot_index == 0 else "Variant"
                
                if plot_idx < len(axes):
                    ax = axes[plot_idx]
                    sns.heatmap(df_to_plot, ax=ax, cbar=False, 
                               annot=False, fmt='.2f', vmin=vmin, vmax=vmax, 
                               cmap=cmap, xticklabels=xticks)
                    
                    # Format the sample name to show complemented bases
                    formatted_name = format_sample_name(sample_name)
                    
                    # Select enrichment type label
                    enrichment_label = "WT" if plot_index == 0 else "Variant"
                    
                    # Create legend items for both SNP location, read count, and log2-FC value
                    legend_items = []
                    legend_labels = []
                    
                    # Add vertical line for SNP position if available
                    if sample_data["original_pos_plot"] is not None:
                        snp_line = ax.axvline(x=sample_data["original_pos_plot"], color='white', 
                                             linestyle='--', linewidth=1, alpha=0.8)
                        legend_items.append(snp_line)
                        legend_labels.append("SNP location")
                    
                    # Add read count to legend
                    if sample_data["read_count"] is not None:
                        # Create a blank line with no visibility for the read count entry
                        read_count_line = ax.plot([], [], ' ', alpha=0)[0]
                        legend_items.append(read_count_line)
                        legend_labels.append(f"Reads: {sample_data['read_count']:,}")
                    
                    # Add log2-FC value and P-value to legend if available
                    if sample_data.get("log2_fc_value") is not None:
                        # Create a blank line with no visibility for the log2-FC entry
                        log2_fc_line = ax.plot([], [], ' ', alpha=0)[0]
                        legend_items.append(log2_fc_line)
                        
                        # Format log2-FC value
                        log2_fc = sample_data["log2_fc_value"]
                        if isinstance(log2_fc, (int, float)):
                            legend_labels.append(f"log2-FC: {log2_fc:.2f}")
                            print(f"  Adding log2-FC to legend: {log2_fc:.2f}")
                        else:
                            legend_labels.append(f"log2-FC: {log2_fc}")
                            print(f"  Adding log2-FC to legend: {log2_fc} (not a number)")
                        
                        # Add P-value if available
                        if sample_data.get("pvalue") is not None:
                            pvalue_line = ax.plot([], [], ' ', alpha=0)[0]
                            legend_items.append(pvalue_line)
                            
                            # Format P-value
                            pvalue = sample_data["pvalue"]
                            if isinstance(pvalue, (int, float)):
                                # Display the original p-value without rounding
                                legend_labels.append(f"P-value: {pvalue}")
                                print(f"  Adding P-value to legend: {pvalue}")
                            else:
                                legend_labels.append(f"P-value: {pvalue}")
                                print(f"  Adding P-value to legend: {pvalue} (not a number)")
                    else:
                        print(f"  No log2-FC value found for sample (single plot): {sample_name}")
                    
                    # Add the legend if we have items to show
                    if legend_items:
                        legend = ax.legend(legend_items, legend_labels, 
                                         loc='lower right', bbox_to_anchor=(1.05, 1.25),
                                         fontsize=8, framealpha=0.7,
                                         title_fontsize=9,
                                         facecolor='dimgray', edgecolor='black')
                        for text in legend.get_texts():
                            text.set_color("white")
                    
                    # Set title for subplot - now only showing enrichment type and formatted name
                    ax.set_title(f"{formatted_name}\n{enrichment_label}", fontsize=10)
                    
                    # Only show y-label on leftmost plots
                    if plot_idx % n_cols == 0:
                        ax.set_ylabel("Footprint Size", fontsize=12)
                    else:
                        ax.set_ylabel("")
                    
                    plot_idx += 1
        
        # Hide unused subplots
        for i in range(plot_idx, len(axes)):
            axes[i].axis('off')
        
        # Make layout tight but leave room for legends and title
        plt.tight_layout()
        # Set axes for titles - leave more space at the top for legends and right for colorbar
        plt.subplots_adjust(top=0.80, right=0.87)
        
        # Save the figure
        if output_directory:
            # Create a clean filename based on the position
            is_merged = sample_info[sample_names[0]].get("is_merged", False) if sample_names else False
            
            if is_merged and sample_names:
                # For merged files, use both positions in the filename
                first_pos = sample_info[sample_names[0]].get("original_pos", "")
                second_pos = sample_info[sample_names[0]].get("second_original_pos", "")
                filename = f"position_{pos}_merged_{first_pos}_{second_pos}_comparison.png"
            else:
                # For regular files, use ref and var bases
                ref_base = sample_info[sample_names[0]].get("ref_base", "") if sample_names else ""
                var_base = sample_info[sample_names[0]].get("var_base", "") if sample_names else ""
                filename = f"position_{pos}_{ref_base}_{var_base}_comparison.png"
                
            save_path = os.path.join(output_directory, filename)
            
            try:
                plt.savefig(save_path, format='png', bbox_inches='tight', dpi=300)
                print(f"Figure successfully saved to {save_path}")
            except Exception as e:
                print(f"Error saving figure: {str(e)}")
                # Try alternative approach with different path
                try:
                    alternative_path = os.path.join(os.getcwd(), filename)
                    print(f"Trying alternative save path: {alternative_path}")
                    plt.savefig(alternative_path, format='png', dpi=300)
                    print(f"Figure saved to alternative path: {alternative_path}")
                except Exception as e2:
                    print(f"Alternative save also failed: {str(e2)}")
        
        # Close the figure to free memory
        plt.close(fig)
        
        return (pos, True)
    
    except Exception as e:
        print(f"Error plotting position group {pos}: {str(e)}")
        print(traceback.format_exc())
        return (pos, False)

# Function to plot position groups sequentially (fallback)
def plot_position_groups_sequential(position_groups, sample_info, output_directory, 
                                   vmin, vmax, cmap, xticks, plot_both, plot_index,
                                   apply_fdr=False, fdr_threshold=0.05):
    """
    Plot position groups sequentially (fallback method if multiprocessing fails).
    
    Parameters:
    -----------
    position_groups : dict
        Dictionary mapping positions to lists of sample names
    sample_info : dict
        Dictionary with sample information
    output_directory : str
        Directory where the output plots will be saved
    vmin : float
        Minimum value for the heatmap color scale
    vmax : float
        Maximum value for the heatmap color scale
    cmap : str
        Colormap for the heatmap
    xticks : int
        Interval for x-axis ticks
    plot_both : bool
        If True, plot both WT and variant enrichment
    plot_index : int
        Index of the dataframe to plot (0 for WT, 1 for Variant)
        
    Returns:
    --------
    int
        Number of successfully plotted position groups
    """
    successful_plots = 0
    total_positions = len(position_groups)
    
    print(f"Plotting {total_positions} position groups sequentially...")
    
    for i, (pos, sample_names) in enumerate(position_groups.items(), 1):
        print(f"Plotting position group {i}/{total_positions}: {pos}")
        try:
            pos_data = (pos, sample_names, sample_info)
            # Add debug output
            for sample_name in sample_names:
                if "log2_fc_value" in sample_info[sample_name] and sample_info[sample_name]["log2_fc_value"] is not None:
                    print(f"  Found log2-FC value for {sample_name}: {sample_info[sample_name]['log2_fc_value']}")
                
            _, success = plot_position_group(pos_data, output_directory, vmin, vmax, 
                                         cmap, xticks, plot_both, plot_index,
                                         apply_fdr=apply_fdr, fdr_threshold=fdr_threshold)
            if success:
                successful_plots += 1
        except Exception as e:
            print(f"Error plotting position group {pos}: {str(e)}")
            print(traceback.format_exc())
    
    return successful_plots

def read_values_from_tsv(tsv_file_path):
    """
    Read values from a TSV file containing genomic information.
    
    Parameters:
    -----------
    tsv_file_path : str
        Path to the TSV file with columns: chromosome, position, Ref, Alt, Tags, DNA, RNA, Value, P-Value
        
    Returns:
    --------
    dict
        Dictionary with keys as tuples of (position, ref_base, alt_base) and values as dict with 'log2_fc' and 'pvalue'
    """
    try:
        print(f"Reading values from TSV file: {tsv_file_path}")
        value_dict = {}
        
        # Determine if file is CSV or TSV by checking extension and content
        delimiter = '\t'  # Default to tab
        if tsv_file_path.lower().endswith('.csv'):
            delimiter = ','
        else:
            # Check the first line to determine delimiter
            with open(tsv_file_path, 'r') as f:
                first_line = f.readline()
                if ',' in first_line and '\t' not in first_line:
                    delimiter = ','
        
        print(f"Using delimiter: '{delimiter}'")
        
        with open(tsv_file_path, 'r') as f:
            # Read header line
            header = f.readline().strip().split(delimiter)
            print(f"Found header: {header}")
            
            # Find the indices of the required columns (case-insensitive)
            header_lower = [h.lower() for h in header]
            try:
                pos_idx = next((i for i, h in enumerate(header_lower) if 'position' in h), None)
                ref_idx = next((i for i, h in enumerate(header_lower) if 'ref' in h), None)
                alt_idx = next((i for i, h in enumerate(header_lower) if 'alt' in h), None)
                value_idx = next((i for i, h in enumerate(header_lower) if 'value' in h and 'p' not in h.lower()), None)
                pvalue_idx = next((i for i, h in enumerate(header_lower) if 'p-value' in h or 'pvalue' in h), None)
                
                if any(idx is None for idx in [pos_idx, ref_idx, alt_idx, value_idx]):
                    missing_cols = []
                    if pos_idx is None: missing_cols.append("position")
                    if ref_idx is None: missing_cols.append("Ref")
                    if alt_idx is None: missing_cols.append("Alt")
                    if value_idx is None: missing_cols.append("Value")
                    print(f"Error: Required columns not found in TSV file: {', '.join(missing_cols)}")
                    print(f"Available columns: {', '.join(header)}")
                    return {}
                else:
                    print(f"Found columns - Position: {header[pos_idx]}, Ref: {header[ref_idx]}, Alt: {header[alt_idx]}, Value: {header[value_idx]}")
                    if pvalue_idx is not None:
                        print(f"Found P-value column: {header[pvalue_idx]}")
                    else:
                        print("P-value column not found in TSV file")
                
            except ValueError as e:
                print(f"Error: Required column not found in TSV file: {e}")
                return {}
            
            # Read data and create lookup dictionary
            line_count = 0
            for line in f:
                line_count += 1
                fields = line.strip().split(delimiter)
                if len(fields) > max(pos_idx, ref_idx, alt_idx, value_idx):
                    try:
                        position = int(fields[pos_idx])
                        ref_base = fields[ref_idx].upper()  # Normalize case
                        alt_base = fields[alt_idx].upper()  # Normalize case
                        log2_fc = float(fields[value_idx])
                        
                        # Get p-value if available
                        pvalue = None
                        if pvalue_idx is not None and len(fields) > pvalue_idx:
                            try:
                                pvalue = float(fields[pvalue_idx])
                            except (ValueError, TypeError):
                                # Handle non-numeric p-values
                                pvalue = None
                        
                        # Store the entry with both values
                        key = (position, ref_base, alt_base)
                        value_dict[key] = {
                            'log2_fc': log2_fc,
                            'pvalue': pvalue
                        }
                        
                        # Print some debug info for the first few rows
                        if line_count <= 5:
                            pval_str = f", P-value={pvalue}" if pvalue is not None else ""
                            print(f"  Row {line_count}: Pos={position}, Ref={ref_base}, Alt={alt_base}, Value={log2_fc}{pval_str}")
                        
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse line {line_count}: {line.strip()}, Error: {e}")
                        continue
                        
        print(f"Successfully read {len(value_dict)} entries from TSV file")
        return value_dict
        
    except Exception as e:
        print(f"Error reading TSV file {tsv_file_path}: {str(e)}")
        return {}

def process_and_plot_grouped_samples(
    input_directory,
    control_pickle_path,
    null_means_dict_path,
    null_stds_dict_path,
    output_directory,
    tsv_file_path=None,  # Add parameter for TSV file
    metadata_path=None,
    ref_seq_length=4718,
    row_range=(1, 400),
    col_range=(3000, 3600),
    bin_size=10,
    vmin=2,
    vmax=10,
    xticks=50,
    cmap="magma",
    plot_index=1,  # By default, plot the second dataframe (variant enrichment) from each sample
    plot_both=False,  # Option to plot both WT and variant enrichment
    n_processes=None,  # Number of processes to use for multiprocessing
    use_multiprocessing=True,  # Whether to use multiprocessing or sequential processing
    min_coverage=None,  # Minimum number of reads required to process a BED file
    apply_fdr=False,  # Whether to apply FDR correction to p-values
    fdr_threshold=0.05  # FDR threshold for significance
):
    """
    Process bed files in a directory using compare_samples, group them by SNP position,
    and plot multiple heatmaps for each position group in a single PDF.
    Uses multiprocessing for faster performance if enabled and available.
    Files with fewer reads than the specified minimum coverage will be skipped.
    Can apply FDR correction to p-values to control for multiple testing.
    
    Parameters:
    -----------
    input_directory : str
        Path to the directory containing BED files to process
    control_pickle_path : str
        Path to the pickled control dataframe
    null_means_dict_path : str
        Path to the pickled dictionary with means of null distributions
    null_stds_dict_path : str
        Path to the pickled dictionary with standard deviations of null distributions
    output_directory : str
        Directory where the output plots will be saved
    tsv_file_path : str, optional
        Path to a TSV file with log2-FC values to display in the plots
    metadata_path : str, optional
        Path to the pickled metadata with dimension information
    ref_seq_length : int
        Length of the reference sequence
    row_range : tuple
        Range of rows to include (min, max)
    col_range : tuple
        Range of columns to include (min, max)
    bin_size : int
        Bin size for footprint sizes
    vmin : float
        Minimum value for the heatmap color scale
    vmax : float
        Maximum value for the heatmap color scale
    xticks : int
        Interval for x-axis ticks
    cmap : str
        Colormap for the heatmap
    plot_index : int
        Index of the dataframe to plot (0 for WT enrichment, 1 for Variant enrichment)
    plot_both : bool
        If True, plot both WT and variant enrichment for each sample
    n_processes : int or None
        Number of processes to use for multiprocessing
    use_multiprocessing : bool
        Whether to use multiprocessing or sequential processing
    min_coverage : int or None
        Minimum number of reads required to process a BED file. Files with fewer reads will be skipped.
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance. Only values passing this threshold will be displayed.
    
    Returns:
    --------
    None (saves plots to the specified output directory)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    print(f"Processing bed files from directory: {input_directory}")
    
    # Adjust vmin if FDR correction is enabled
    if apply_fdr:
        fdr_vmin = -np.log10(fdr_threshold)
        if vmin < fdr_vmin:
            original_vmin = vmin
            vmin = fdr_vmin
            print(f"FDR correction enabled. Adjusted vmin from {original_vmin} to {vmin:.2f} based on FDR threshold {fdr_threshold}")
        else:
            print(f"FDR correction enabled with threshold {fdr_threshold} (vmin = {vmin})")
    
    # Load values from TSV file if provided
    value_dict = {}
    if tsv_file_path:
        value_dict = read_values_from_tsv(tsv_file_path)
    
    # Step 1: Process all samples using compare_samples
    start_time = time.time()
    delta_df_dict = compare_samples(
        input_directory, 
        control_pickle_path,
        null_means_dict_path, 
        null_stds_dict_path, 
        metadata_path=metadata_path,
        ref_seq_length=ref_seq_length, 
        row_range=row_range, 
        col_range=col_range, 
        bin_size=bin_size,
        n_processes=n_processes,
        use_multiprocessing=use_multiprocessing,
        apply_fdr=apply_fdr,
        fdr_threshold=fdr_threshold
    )
    
    processing_time = time.time() - start_time
    print(f"Processing complete. Found {len(delta_df_dict)} samples in {processing_time:.2f} seconds.")
    
    # Step 2: Extract SNP information for all samples and group by position
    sample_info = {}
    position_groups = {}
    
    for sample_name, (wt_df, var_df) in delta_df_dict.items():
        print(f"Processing sample metadata: {sample_name}")
        
        # Extract SNP position and base change information from the available data
        original_pos, genomic_pos, ref_base, var_base = None, None, None, None
        is_merged = False
        second_original_pos, second_genomic_pos = None, None
        bed_file_path = None
        
        # Try to extract from sample name directly or from any related filenames
        sample_base_name = None
        for filename in os.listdir(input_directory):
            if filename.endswith("_fp.bed"):
                result = parse_bed_filename(filename)
                if result[1] is not None and filename.split('_fp.bed')[0] in sample_name:
                    if len(result) >= 5 and result[4]:  # It's a merged file
                        original_pos, genomic_pos = result[0], result[1]
                        is_merged = True
                        second_original_pos, second_genomic_pos = result[5], result[6]
                        ref_base, var_base = None, None  # No base info for merged files
                    else:  # Regular single SNP file
                        original_pos, genomic_pos, ref_base, var_base = result[:4]
                    sample_base_name = filename
                    bed_file_path = os.path.join(input_directory, filename)
                    break
        
        if sample_base_name:
            if is_merged:
                print(f"  Extracted merged SNP info from {sample_base_name}: original_pos={original_pos}/{second_original_pos}, genomic_pos={genomic_pos}/{second_genomic_pos}")
            else:
                print(f"  Extracted SNP info from {sample_base_name}: original_pos={original_pos}, genomic_pos={genomic_pos}, ref={ref_base}, var={var_base}")
        else:
            print("  Could not determine SNP position from filenames")
            # Skip samples with unknown position
            continue
        
        # Count reads in the BED file if found
        read_count = None
        if bed_file_path:
            read_count = count_bed_file_rows(bed_file_path)
            print(f"  Read count for {sample_base_name}: {read_count}")
            
        # Skip samples with insufficient coverage if min_coverage is specified
        if min_coverage is not None and (read_count is None or read_count < min_coverage):
            print(f"  Skipping {sample_name}: Read count {read_count} is below the minimum coverage threshold of {min_coverage}")
            continue
        
        # Adjust original position to be relative to the col_range for plotting
        original_pos_plot = None
        if original_pos is not None and col_range is not None:
            if col_range[0] <= original_pos <= col_range[1]:
                # Convert absolute position to position relative to the filtered range
                original_pos_plot = original_pos - col_range[0]
                print(f"  Adjusted original position for plotting: {original_pos_plot} (original: {original_pos}, range: {col_range})")
            else:
                print(f"  Original position {original_pos} is outside the plotted range {col_range}")
        
        # Look up the value from the TSV file if available
        log2_fc_value = None
        if value_dict and genomic_pos is not None and ref_base is not None and var_base is not None:
            # Try different key combinations to find a match
            possible_keys = [
                (genomic_pos, ref_base, var_base),  # Exact match
                (genomic_pos, ref_base.upper(), var_base.upper()),  # Case normalized
                (genomic_pos, ref_base.lower(), var_base.lower())   # Lowercase
            ]
            
            print(f"  Looking up log2-FC for position: {genomic_pos}, ref: {ref_base}, var: {var_base}")
            for key in possible_keys:
                if key in value_dict:
                    log2_fc_value = value_dict[key]
                    print(f"  Found log2-FC value for {key}: {log2_fc_value}")
                    break
            
            if log2_fc_value is None:
                print(f"  No log2-FC value found for position: {genomic_pos}, ref: {ref_base}, var: {var_base}")
                print(f"  Available keys in value_dict: {list(value_dict.keys())[:5]}...")
        else:
            if not value_dict:
                print(f"  No log2-FC dictionary available")
            else:
                print(f"  Missing lookup data - genomic_pos: {genomic_pos}, ref_base: {ref_base}, var_base: {var_base}")
        # Store the sample information
        sample_info[sample_name] = {
            "wt_df": wt_df,
            "var_df": var_df,
            "original_pos": original_pos,  # Original position from filename
            "genomic_pos": genomic_pos,    # Genomic position (11092732 - original_pos)
            "ref_base": ref_base,          # Complemented reference base
            "var_base": var_base,          # Complemented variant base
            "read_count": read_count,      # Number of reads
            "original_pos_plot": original_pos_plot,  # Original position adjusted for plotting
            "log2_fc_value": log2_fc_value  # log2 fold change value from TSV
        }
        
        # Group by position
        if genomic_pos is not None:
            if genomic_pos not in position_groups:
                position_groups[genomic_pos] = []
            position_groups[genomic_pos].append(sample_name)
    
    print(f"Grouped samples into {len(position_groups)} position groups")
    
    if len(position_groups) == 0:
        print("No valid position groups found. Check the input files and parsing logic.")
        return
    
    # Step 3: Create a multi-panel figure for each position group
    if use_multiprocessing and n_processes != 1 and len(position_groups) > 1:
        # Determine number of processes to use for plotting
        if n_processes is None:
            plot_processes = max(1, min(mp.cpu_count() - 1, len(position_groups)))
        else:
            plot_processes = max(1, min(n_processes, mp.cpu_count(), len(position_groups)))
        
        print(f"Creating plots for {len(position_groups)} position groups using {plot_processes} processes...")
        
        # Prepare data for parallel processing
        position_data = [
            (pos, sample_names, sample_info) 
            for pos, sample_names in position_groups.items()
        ]
        
        # Create a partial function with fixed parameters
        plot_func = partial(
            plot_position_group,
            output_directory=output_directory,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            xticks=xticks,
            plot_both=plot_both,
            plot_index=plot_index,
            apply_fdr=apply_fdr,
            fdr_threshold=fdr_threshold
        )
        
        try:
            # Plot position groups in parallel
            start_time = time.time()
            
            with mp.Pool(processes=plot_processes) as pool:
                results = pool.map(plot_func, position_data)
            
            # Count successful plots
            successful_plots = sum(1 for _, success in results if success)
            
            plotting_time = time.time() - start_time
            print(f"Plotting complete. Successfully created {successful_plots} out of {len(position_groups)} plots in {plotting_time:.2f} seconds.")
        
        except Exception as e:
            print(f"Error in parallel plotting: {str(e)}")
            print(traceback.format_exc())
            print("Falling back to sequential plotting...")
            
            # Fall back to sequential plotting
            start_time = time.time()
            successful_plots = plot_position_groups_sequential(
                position_groups, sample_info, output_directory, 
                vmin, vmax, cmap, xticks, plot_both, plot_index,
                apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
            )
            plotting_time = time.time() - start_time
            print(f"Sequential plotting complete. Successfully created {successful_plots} out of {len(position_groups)} plots in {plotting_time:.2f} seconds.")
    else:
        # Use sequential plotting
        print("Using sequential plotting...")
        start_time = time.time()
        successful_plots = plot_position_groups_sequential(
            position_groups, sample_info, output_directory, 
            vmin, vmax, cmap, xticks, plot_both, plot_index,
            apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
        )
        plotting_time = time.time() - start_time
        print(f"Sequential plotting complete. Successfully created {successful_plots} out of {len(position_groups)} plots in {plotting_time:.2f} seconds.")
    
    print(f"All position groups processed and plotted. Results saved to {output_directory}")

def parse_arguments():
    """Parse command-line arguments for the script."""
    parser = argparse.ArgumentParser(description="Generate grouped heatmaps for SNP footprint analysis.")
    
    # Required arguments
    parser.add_argument("--input-dir", required=True, 
                        help="Directory containing BED files to process")
    parser.add_argument("--control-pkl", required=True, 
                        help="Path to the pickled control dataframe")
    parser.add_argument("--null-means-pkl", required=True, 
                        help="Path to the pickled dictionary with means of null distributions")
    parser.add_argument("--null-stds-pkl", required=True, 
                        help="Path to the pickled dictionary with standard deviations of null distributions")
    parser.add_argument("--output-dir", required=True, 
                        help="Directory where the output plots will be saved")
    
    # Optional arguments
    parser.add_argument("--metadata-pkl", 
                        help="Path to the pickled metadata with dimension information")
    parser.add_argument("--ref-seq-length", type=int, default=4718, 
                        help="Length of the reference sequence")
    parser.add_argument("--row-min", type=int, default=1, 
                        help="Minimum row value for the range")
    parser.add_argument("--row-max", type=int, default=400, 
                        help="Maximum row value for the range")
    parser.add_argument("--col-min", type=int, default=3000, 
                        help="Minimum column value for the range")
    parser.add_argument("--col-max", type=int, default=3600, 
                        help="Maximum column value for the range")
    parser.add_argument("--bin-size", type=int, default=10, 
                        help="Bin size for footprint sizes")
    parser.add_argument("--vmin", type=float, default=2, 
                        help="Minimum value for the heatmap color scale")
    parser.add_argument("--vmax", type=float, default=10, 
                        help="Maximum value for the heatmap color scale")
    parser.add_argument("--xticks", type=int, default=50, 
                        help="Interval for x-axis ticks")
    parser.add_argument("--cmap", default="magma", 
                        help="Colormap for the heatmap")
    parser.add_argument("--plot-wt", action="store_true", 
                        help="Plot WT enrichment (default is variant)")
    parser.add_argument("--plot-both", action="store_true", 
                        help="Plot both WT and variant enrichment")
    parser.add_argument("--processes", type=int, default=None,
                        help="Number of processes to use for parallel processing")
    parser.add_argument("--no-mp", action="store_true",
                        help="Disable multiprocessing and use sequential processing")
    parser.add_argument("--tsv-file", 
                        help="Path to a TSV file with log2-FC values to display in the plots")
    parser.add_argument("--min-coverage", type=int, default=None,
                        help="Minimum number of reads required to process a BED file")
    parser.add_argument("--apply-fdr", action="store_true",
                        help="Apply Benjamini-Hochberg FDR correction to p-values")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                        help="FDR threshold for significance (default: 0.05)")
    
    return parser.parse_args()

def main():
    """Main function to execute the script."""
    start_time = time.time()
    print("Grouped Heatmap Generator for SNP Footprint Analysis")
    print("---------------------------------------------------")
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set plot_index based on --plot-wt flag
    plot_index = 0 if args.plot_wt else 1
    
    # Process data and create plots
    process_and_plot_grouped_samples(
        input_directory=args.input_dir,
        control_pickle_path=args.control_pkl,
        null_means_dict_path=args.null_means_pkl,
        null_stds_dict_path=args.null_stds_pkl,
        output_directory=args.output_dir,
        tsv_file_path=args.tsv_file,
        metadata_path=args.metadata_pkl,
        ref_seq_length=args.ref_seq_length,
        row_range=(args.row_min, args.row_max),
        col_range=(args.col_min, args.col_max),
        bin_size=args.bin_size,
        vmin=args.vmin,
        vmax=args.vmax,
        xticks=args.xticks,
        cmap=args.cmap,
        plot_index=plot_index,
        plot_both=args.plot_both,
        n_processes=args.processes,
        use_multiprocessing=not args.no_mp,
        min_coverage=args.min_coverage,
        apply_fdr=args.apply_fdr,
        fdr_threshold=args.fdr_threshold
    )
    
    total_time = time.time() - start_time
    print(f"Script execution completed in {total_time:.2f} seconds.")

if __name__ == "__main__":
    main()