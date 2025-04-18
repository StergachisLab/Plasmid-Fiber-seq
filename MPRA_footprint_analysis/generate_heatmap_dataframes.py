#!/usr/bin/env python3
"""
Generate Heatmap Dataframes for SNP Footprint Analysis

This script processes BED files containing SNP data, groups them by position,
and generates dataframes for WT and Variant enrichment analysis without plotting.
It uses core functions from the original grouped heatmap generator but focuses
on data processing rather than visualization.
"""

import os
import sys
import argparse
import warnings
import traceback
import pandas as pd
import pickle
import numpy as np
from scipy import stats
from scipy.stats import rankdata
import multiprocessing as mp
from functools import partial
import time
sys.path.append('/gscratch/stergachislab/bmallo/large_home/python_scripts/')
from FiberHMM_functions import *

# Import necessary functions from the original script
# Ensure these functions are in the same directory or in the Python path
from FiberHMM_functions import read_ft_data, grab_circular_reads, prep_dfs_for_subtraction, filter_fp_df

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
        (sample_name, (WT_pvalue, Var_pvalue, processed_count_df, WT_enrichment, Var_enrichment))
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
        return (sample_name, (WT_pvalue, Var_pvalue, processed_count_df, WT_enrichment, Var_enrichment))
    
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
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue, processed_count_df, WT_enrichment, Var_enrichment) DataFrames
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
    result_dict = {}
    
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
        for sample_name, result_tuple in results:
            if result_tuple is not None:
                result_dict[sample_name] = result_tuple
        
        end_time = time.time()
        print(f"All samples processed in {end_time - start_time:.2f} seconds")
        print(f"Successfully processed {len(result_dict)} out of {total_samples} samples")
        
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
    if len(result_dict) == 0:
        print("No samples were processed successfully with multiprocessing. Falling back to sequential processing...")
        return compare_samples_sequential(
            sample_directory_path, control_df, null_means_dict, null_stds_dict,
            ref_seq_length, row_range, col_range, bin_size,
            apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
        )
    
    return result_dict

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
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue, processed_count_df, WT_enrichment, Var_enrichment) DataFrames
    """
    # Setup dictionary to store processed dataframes
    result_dict = {}
    
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

            # Store all the dataframes
            result_dict[sample_name] = (WT_pvalue, Var_pvalue, processed_count_df, WT_enrichment, Var_enrichment)
            print(f"  Finished processing {sample_name}")
        except Exception as e:
            print(f"Error processing sample {sample_name}: {str(e)}")
            print(traceback.format_exc())
    
    print(f"Successfully processed {len(result_dict)} out of {total_samples} samples")
    return result_dict

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
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue, processed_count_df, WT_enrichment, Var_enrichment) DataFrames
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

def extract_sample_metadata(sample_directory_path, result_dict, min_coverage=None):
    """
    Extract metadata for each sample from BED filenames.
    Optionally filter samples based on minimum read coverage.
    
    Parameters:
    -----------
    sample_directory_path : str
        Path to the directory containing BED files
    result_dict : dict
        Dictionary of processed sample dataframes
    min_coverage : int, optional
        Minimum number of reads required for a sample to be included
        
    Returns:
    --------
    tuple
        (sample_metadata, filtered_result_dict) where:
        - sample_metadata: Dictionary mapping sample names to metadata dictionaries
        - filtered_result_dict: Dictionary of processed sample dataframes after coverage filtering
    """
    sample_metadata = {}
    filtered_result_dict = {}
    filtered_count = 0
    
    for sample_name in result_dict.keys():
        print(f"Extracting metadata for sample: {sample_name}")
        
        # Extract SNP position and base change information from the available data
        original_pos, genomic_pos, ref_base, var_base = None, None, None, None
        is_merged = False
        second_original_pos, second_genomic_pos = None, None
        bed_file_path = None
        
        # Try to extract from sample name directly or from any related filenames
        sample_base_name = None
        for filename in os.listdir(sample_directory_path):
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
                    bed_file_path = os.path.join(sample_directory_path, filename)
                    break
        
        # Count reads in the BED file if found
        read_count = None
        if bed_file_path:
            read_count = count_bed_file_rows(bed_file_path)
            print(f"  Read count for {sample_name}: {read_count}")
        
        # Check if sample meets minimum coverage requirement
        if min_coverage is not None and (read_count is None or read_count < min_coverage):
            print(f"  Filtering out {sample_name}: Read count {read_count} is below the minimum coverage threshold of {min_coverage}")
            filtered_count += 1
            continue
        
        # Store the sample metadata
        sample_metadata[sample_name] = {
            "original_pos": original_pos,  # Original position from filename
            "genomic_pos": genomic_pos,    # Genomic position (11092732 - original_pos)
            "ref_base": ref_base,          # Complemented reference base
            "var_base": var_base,          # Complemented variant base
            "read_count": read_count,      # Number of reads
            "is_merged": is_merged,
            "second_original_pos": second_original_pos,
            "second_genomic_pos": second_genomic_pos
        }
        
        # Include this sample in the filtered result dictionary
        filtered_result_dict[sample_name] = result_dict[sample_name]
    
    if min_coverage is not None:
        print(f"Coverage filtering: {filtered_count} samples filtered out, {len(filtered_result_dict)} samples retained")
    
    return sample_metadata, filtered_result_dict

def create_merged_pvalue_dataframes(result_dict):
    """
    Create merged dataframes containing the maximum p-value at each position
    for all WT and variant samples.
    
    Parameters:
    -----------
    result_dict : dict
        Dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue, processed_df, WT_enrichment, Var_enrichment)
    
    Returns:
    --------
    tuple
        (merged_wt_pvalue_df, merged_var_pvalue_df)
    """
    if not result_dict:
        print("No samples available to create merged dataframes")
        return None, None
    
    # Get the first sample to initialize the merged dataframes
    first_sample = next(iter(result_dict.values()))
    wt_pvalue_df = first_sample[0]
    var_pvalue_df = first_sample[1]
    
    # Initialize merged dataframes with zeros
    merged_wt_pvalue_df = pd.DataFrame(0, index=wt_pvalue_df.index, columns=wt_pvalue_df.columns)
    merged_var_pvalue_df = pd.DataFrame(0, index=var_pvalue_df.index, columns=var_pvalue_df.columns)
    
    print(f"Creating merged p-value dataframes from {len(result_dict)} samples...")
    
    # Iterate through all samples and update the merged dataframes with the maximum values
    for sample_name, (wt_pvalue, var_pvalue, _, _, _) in result_dict.items():
        # For WT p-values
        for row in wt_pvalue.index:
            for col in wt_pvalue.columns:
                merged_wt_pvalue_df.loc[row, col] = max(
                    merged_wt_pvalue_df.loc[row, col],
                    wt_pvalue.loc[row, col]
                )
        
        # For variant p-values
        for row in var_pvalue.index:
            for col in var_pvalue.columns:
                merged_var_pvalue_df.loc[row, col] = max(
                    merged_var_pvalue_df.loc[row, col],
                    var_pvalue.loc[row, col]
                )
    
    print("Successfully created merged p-value dataframes")
    return merged_wt_pvalue_df, merged_var_pvalue_df

def save_dataframes_to_pickle(result_dict, sample_metadata, output_directory):
    """
    Save the processed dataframes to pickle files in the specified output directory.
    Also creates and saves merged dataframes containing maximum p-values.
    
    Parameters:
    -----------
    result_dict : dict
        Dictionary mapping sample names to tuples of dataframes
    sample_metadata : dict
        Dictionary mapping sample names to metadata
    output_directory : str
        Directory where the pickle files will be saved
    
    Returns:
    --------
    tuple
        (merged_wt_pvalue_df, merged_var_pvalue_df)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    # Create merged dataframes with maximum p-values
    merged_wt_pvalue_df, merged_var_pvalue_df = create_merged_pvalue_dataframes(result_dict)
    
    # Create a dictionary to store all data
    output_data = {
        'result_dict': result_dict,
        'sample_metadata': sample_metadata,
        'merged_wt_pvalue': merged_wt_pvalue_df,
        'merged_var_pvalue': merged_var_pvalue_df
    }
    
    # Save the entire data structure to a single pickle file
    output_path = os.path.join(output_directory, 'heatmap_dataframes.pkl')
    print(f"Saving all dataframes to {output_path}")
    
    try:
        with open(output_path, 'wb') as f:
            pickle.dump(output_data, f)
        print(f"Successfully saved all dataframes to {output_path}")
    except Exception as e:
        print(f"Error saving dataframes: {str(e)}")
        print(traceback.format_exc())
    
    # Save merged dataframes separately
    try:
        if merged_wt_pvalue_df is not None and merged_var_pvalue_df is not None:
            merged_dir = os.path.join(output_directory, 'merged')
            os.makedirs(merged_dir, exist_ok=True)
            
            merged_wt_path = os.path.join(merged_dir, 'merged_wt_pvalue.pkl')
            merged_var_path = os.path.join(merged_dir, 'merged_var_pvalue.pkl')
            
            merged_wt_pvalue_df.to_pickle(merged_wt_path)
            merged_var_pvalue_df.to_pickle(merged_var_path)
            
            print(f"Saved merged WT p-value dataframe to {merged_wt_path}")
            print(f"Saved merged variant p-value dataframe to {merged_var_path}")
    except Exception as e:
        print(f"Error saving merged dataframes: {str(e)}")
        print(traceback.format_exc())
    
    # Also save individual sample data for easier access
    for sample_name, (wt_pvalue, var_pvalue, processed_df, wt_enrichment, var_enrichment) in result_dict.items():
        try:
            # Create a clean filename based on the sample name
            safe_name = "".join([c if c.isalnum() or c in [' ', '_', '-'] else '_' for c in sample_name])
            sample_dir = os.path.join(output_directory, safe_name)
            os.makedirs(sample_dir, exist_ok=True)
            
            # Save each dataframe separately
            dataframes = {
                'wt_pvalue': wt_pvalue,
                'var_pvalue': var_pvalue, 
                'processed_counts': processed_df,
                'wt_enrichment': wt_enrichment,
                'var_enrichment': var_enrichment
            }
            
            for df_name, df in dataframes.items():
                df_path = os.path.join(sample_dir, f"{df_name}.pkl")
                df.to_pickle(df_path)
            
            # Save the metadata
            if sample_name in sample_metadata:
                metadata_path = os.path.join(sample_dir, 'metadata.pkl')
                with open(metadata_path, 'wb') as f:
                    pickle.dump(sample_metadata[sample_name], f)
            
            print(f"Successfully saved dataframes for {sample_name}")
        except Exception as e:
            print(f"Error saving dataframes for {sample_name}: {str(e)}")
    
    print(f"All dataframes saved to {output_directory}")
    return merged_wt_pvalue_df, merged_var_pvalue_df
def generate_heatmap_dataframes(
    input_directory,
    control_pickle_path,
    null_means_dict_path,
    null_stds_dict_path,
    output_directory,
    metadata_path=None,
    ref_seq_length=4718,
    row_range=(1, 400),
    col_range=(3000, 3600),
    bin_size=10,
    min_coverage=None,
    n_processes=None,
    use_multiprocessing=True,
    apply_fdr=False,
    fdr_threshold=0.05
):
    """
    Process BED files to generate dataframes for SNP footprint analysis without plotting.
    
    Parameters:
    -----------
    input_directory : str
        Path to the directory containing BED files
    control_pickle_path : str
        Path to the pickled control dataframe
    null_means_dict_path : str
        Path to the pickled dictionary with means of null distributions
    null_stds_dict_path : str
        Path to the pickled dictionary with standard deviations of null distributions
    output_directory : str
        Directory where the output pickle files will be saved
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
    min_coverage : int, optional
        Minimum number of reads required to process a BED file
    n_processes : int or None
        Number of processes to use for parallel processing
    use_multiprocessing : bool
        Whether to use multiprocessing or sequential processing
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    tuple
        (result_dict, sample_metadata, merged_wt_pvalue_df, merged_var_pvalue_df)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    print(f"Processing bed files from directory: {input_directory}")
    
    if min_coverage is not None:
        print(f"Minimum read coverage threshold: {min_coverage}")
    
    # Step 1: Process all samples using compare_samples
    start_time = time.time()
    result_dict = compare_samples(
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
    print(f"Processing complete. Found {len(result_dict)} samples in {processing_time:.2f} seconds.")
    
    # Step 2: Extract metadata and filter by coverage if needed
    sample_metadata, filtered_result_dict = extract_sample_metadata(input_directory, result_dict, min_coverage)
    
    # Step 3: Save dataframes to pickle files and get merged dataframes
    merged_wt_pvalue_df, merged_var_pvalue_df = save_dataframes_to_pickle(
        filtered_result_dict, sample_metadata, output_directory
    )
    
    return filtered_result_dict, sample_metadata, merged_wt_pvalue_df, merged_var_pvalue_df

def parse_arguments():
    """Parse command-line arguments for the script."""
    parser = argparse.ArgumentParser(
        description="Generate dataframes for SNP footprint analysis without plotting."
    )
    
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
                        help="Directory where the output pickle files will be saved")
    
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
    parser.add_argument("--min-coverage", type=int, default=None,
                        help="Minimum number of reads required to process a BED file")
    parser.add_argument("--processes", type=int, default=None,
                        help="Number of processes to use for parallel processing")
    parser.add_argument("--no-mp", action="store_true",
                        help="Disable multiprocessing and use sequential processing")
    parser.add_argument("--apply-fdr", action="store_true",
                        help="Apply Benjamini-Hochberg FDR correction to p-values")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                        help="FDR threshold for significance (default: 0.05)")
    
    return parser.parse_args()

def main():
    """Main function to execute the script."""
    start_time = time.time()
    print("Heatmap Dataframe Generator for SNP Footprint Analysis")
    print("---------------------------------------------------")
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Generate dataframes
    result_dict, sample_metadata, merged_wt_pvalue_df, merged_var_pvalue_df = generate_heatmap_dataframes(
        input_directory=args.input_dir,
        control_pickle_path=args.control_pkl,
        null_means_dict_path=args.null_means_pkl,
        null_stds_dict_path=args.null_stds_pkl,
        output_directory=args.output_dir,
        metadata_path=args.metadata_pkl,
        ref_seq_length=args.ref_seq_length,
        row_range=(args.row_min, args.row_max),
        col_range=(args.col_min, args.col_max),
        bin_size=args.bin_size,
        min_coverage=args.min_coverage,
        n_processes=args.processes,
        use_multiprocessing=not args.no_mp,
        apply_fdr=args.apply_fdr,
        fdr_threshold=args.fdr_threshold
    )
    
    total_time = time.time() - start_time
    print(f"Script execution completed in {total_time:.2f} seconds.")
    print(f"Processed {len(result_dict)} samples.")
    
    # Print information about merged dataframes
    if merged_wt_pvalue_df is not None and merged_var_pvalue_df is not None:
        print(f"Created merged WT p-value dataframe with shape {merged_wt_pvalue_df.shape}")
        print(f"Created merged variant p-value dataframe with shape {merged_var_pvalue_df.shape}")
    
    print(f"Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()