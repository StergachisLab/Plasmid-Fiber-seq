#!/usr/bin/env python3
"""
Data Processing Functions for SNP Footprint Analysis

This module contains core functions for processing BED files and footprint data.
"""

import os
import sys
import warnings
import traceback
import numpy as np
import pandas as pd
from functools import partial
import multiprocessing as mp
import time

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

def grab_circular_reads(input_df, ref_seq_length):
    """
    Process a bed DataFrame to create an expanded footprint DataFrame,
    handling circular reads that wrap around the reference sequence.
    
    Note: This is a placeholder function that should be implemented 
    based on the FiberHMM_functions module referenced in the original script.
    
    Parameters:
    -----------
    input_df : pandas.DataFrame
        The input bed DataFrame
    ref_seq_length : int
        Length of the reference sequence
        
    Returns:
    --------
    pandas.DataFrame
        The expanded footprint DataFrame with circular reads
    """
    # Implementation required based on FiberHMM_functions.py
    raise NotImplementedError("This function needs to be implemented based on FiberHMM_functions module")

def prep_dfs_for_subtraction(circular_footprint_df):
    """
    Transform expanded footprint dataframe into count of footprint sizes dataframe.
    
    Note: This is a placeholder function that should be implemented 
    based on the FiberHMM_functions module referenced in the original script.
    
    Parameters:
    -----------
    circular_footprint_df : pandas.DataFrame
        The circular footprint DataFrame
        
    Returns:
    --------
    list of pandas.DataFrame
        A list containing the count of footprint sizes DataFrame
    """
    # Implementation required based on FiberHMM_functions.py
    raise NotImplementedError("This function needs to be implemented based on FiberHMM_functions module")

def filter_fp_df(footprint_df, bin_size, row_range, col_range):
    """
    Filter footprint_count_dfs for desired footprint sizes and positions.
    
    Note: This is a placeholder function that should be implemented 
    based on the FiberHMM_functions module referenced in the original script.
    
    Parameters:
    -----------
    footprint_df : pandas.DataFrame
        The footprint count DataFrame
    bin_size : int or None
        Bin size for footprint sizes
    row_range : tuple or None
        Range of rows to include
    col_range : tuple or None
        Range of columns to include
        
    Returns:
    --------
    pandas.DataFrame
        The filtered footprint count DataFrame
    """
    # Implementation required based on FiberHMM_functions.py
    raise NotImplementedError("This function needs to be implemented based on FiberHMM_functions module")

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

        # Import the statistics module for calculating p-values
        from footprint_analysis.analysis.statistics import calculate_pvalues_one_tailed_dict

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

def compare_samples_sequential(sample_directory_path, control_df, null_means_dict, null_stds_dict,
                               ref_seq_length, row_range=None, col_range=None, bin_size=None,
                               apply_fdr=False, fdr_threshold=0.05):
    """
    Compare variant and control samples using sequential processing instead of multiprocessing.
    This is a fallback function to use if multiprocessing encounters errors.
    
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
    apply_fdr : bool
        Whether to apply FDR correction to p-values
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    dictionary mapping sample names to tuples of (WT_pvalue, Var_pvalue) DataFrames
    """
    # Import the file I/O module
    from footprint_analysis.core.file_io import read_ft_data
    
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

            # Import the statistics module
            from footprint_analysis.analysis.statistics import calculate_pvalues_one_tailed_dict

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
    # Import the file I/O module
    from footprint_analysis.core.file_io import read_ft_data
    
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
