#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to precompute mean and standard deviation for each cell in a null distribution dataframe
and save them as pickled dictionaries for faster lookups.
"""

import os
import time
import pickle
import argparse
import numpy as np
import pandas as pd


def precompute_null_stats_as_dict(null_distribution_pickle_path, output_dir=None):
    """
    Precompute mean and standard deviation for each cell in the null distribution dataframe
    and save them as pickled dictionaries for faster lookups.
    
    Parameters:
    -----------
    null_distribution_pickle_path : str
        Path to the pickled dataframe containing null distributions for each cell
    output_dir : str, optional
        Directory to save the pickled dictionaries. If None, saves in the same directory as the input file.
        
    Returns:
    --------
    tuple of (means_dict_path, stds_dict_path, metadata_path) as strings
    """
    try:
        # Load the null distribution dataframe
        print(f"Loading null distribution from {null_distribution_pickle_path}...")
        null_distribution_df = pd.read_pickle(null_distribution_pickle_path)
        print(f"Loaded null distribution with shape {null_distribution_df.shape}")
    except Exception as e:
        raise Exception(f"Failed to load null distribution from {null_distribution_pickle_path}: {e}")
    
    # Determine output directory
    if output_dir is None:
        output_dir = os.path.dirname(null_distribution_pickle_path)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate base filename from the input file
    base_filename = os.path.splitext(os.path.basename(null_distribution_pickle_path))[0]
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    
    # Initialize dictionaries to store means and standard deviations
    means_dict = {}
    stds_dict = {}
    
    # Also store index and column information for dimension validation later
    index_list = list(null_distribution_df.index)
    columns_list = list(null_distribution_df.columns)
    
    # Calculate mean and std for each cell and store in dictionaries
    print("Calculating means and standard deviations...")
    for row in null_distribution_df.index:
        for col in null_distribution_df.columns:
            null_distribution = null_distribution_df.loc[row, col]
            means_dict[(row, col)] = np.mean(null_distribution)
            stds_dict[(row, col)] = np.std(null_distribution)
    
    # Create metadata dictionary with dimension information
    metadata = {
        'index': index_list,
        'columns': columns_list,
        'shape': null_distribution_df.shape
    }
    
    # Save the dictionaries and metadata as pickled files
    means_dict_path = os.path.join(output_dir, f"{base_filename}_means_dict_{timestamp}.pkl")
    stds_dict_path = os.path.join(output_dir, f"{base_filename}_stds_dict_{timestamp}.pkl")
    metadata_path = os.path.join(output_dir, f"{base_filename}_metadata_{timestamp}.pkl")
    
    print(f"Saving means dictionary to {means_dict_path}")
    with open(means_dict_path, 'wb') as f:
        pickle.dump(means_dict, f)
    
    print(f"Saving standard deviations dictionary to {stds_dict_path}")
    with open(stds_dict_path, 'wb') as f:
        pickle.dump(stds_dict, f)
    
    print(f"Saving metadata to {metadata_path}")
    with open(metadata_path, 'wb') as f:
        pickle.dump(metadata, f)
    
    print("Precomputation completed successfully")
    return means_dict_path, stds_dict_path, metadata_path


def main():
    """
    Main function to parse command line arguments and run the precomputation.
    """
    parser = argparse.ArgumentParser(
        description='Precompute means and standard deviations from null distribution.'
    )
    parser.add_argument(
        'null_distribution_path',
        type=str,
        help='Path to the pickled dataframe containing null distributions'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Directory to save the pickled dictionaries (default: same as input file)'
    )
    
    args = parser.parse_args()
    
    try:
        means_path, stds_path, metadata_path = precompute_null_stats_as_dict(
            args.null_distribution_path,
            args.output_dir
        )
        print("\nSummary:")
        print(f"- Means dictionary saved to: {means_path}")
        print(f"- Standard deviations dictionary saved to: {stds_path}")
        print(f"- Metadata saved to: {metadata_path}")
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
