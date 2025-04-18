#!/usr/bin/env python3
"""
Statistical Analysis Functions for SNP Footprint Analysis

This module contains functions for statistical analysis, including calculation of p-values
and multiple testing correction.
"""

import numpy as np
import pandas as pd
from scipy import stats

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
