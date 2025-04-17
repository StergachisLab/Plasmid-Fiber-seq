#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single Sample Average Distribution Script
-----------------------------------
Author: Ben Mallory
Last Updated: 2025-03-15
----------------------------------- 

This script takes in a footprint dataframe from FiberHMM and generates
an average distribution of footprint occupancy across multiple subsamples
of a dataset. Unlike the null distribution script which calculates differences
between pairs of subsamples, this script calculates the average occupancy
from single subsamples across iterations.
"""

import os
import pandas as pd
import tqdm
import time
import multiprocessing
from functools import partial
from collections import defaultdict
import sys
sys.path.append('/gscratch/stergachislab/bmallo/large_home/python_scripts/')
from FiberHMM_functions import *
import argparse

def read_ft_data_single(bed_file):
        df = pd.read_csv(bed_file,
            sep='\t',
            header=None,
            usecols=[0, 1, 2, 3, 8, 9]
        )  # Read the BED file
        df = df.reset_index(drop=True)
        df.columns = ['chrom', 'start', 'end', 'name', 'blockStarts', 'blockSizes']
        df = df.loc[df['blockStarts'] != '.']  # Filter rows where blockStarts is not '.'

        # Remove duplicate reads with worse alignment
        df['length'] = df['end'] - df['start']
        df = df.sort_values('length', ascending=False).drop_duplicates('name')

        return df  # Return the bed file as a dataframe

def process_single_sample(input_df, subsample_size, ref_seq_length, row_range, col_range, bin_size, iteration_idx=None):
    """
    Process a single sample iteration.
    This function is designed to be called by multiprocessing.
    
    Parameters:
    -----------
    input_df : pandas.DataFrame
        The bed DataFrame
    subsample_size : int
        Number of reads to include in the subsample
    ref_seq_length : int
        Length of the reference sequence
    row_range : tuple or None
        Range of rows to include
    col_range : tuple or None
        Range of columns to include
    bin_size : int or None
        Bin size for footprint sizes
    iteration_idx : int, optional
        Index of the iteration (for tracking purposes)
        
    Returns:
    --------
    The filtered footprint count DataFrame for this iteration
    """

    # Create a unique seed using iteration index, process ID, and current time
    if iteration_idx is not None:
        seed = int(iteration_idx) + os.getpid() + int(time.time() * 1000) % 10000
        np.random.seed(seed)
        random_state = np.random.RandomState(seed)
    else:
        random_state = np.random.RandomState()

    # Take a single subsample
    subsample = input_df.sample(n=subsample_size, replace=False, random_state=random_state)
    
    # Process bed file into expanded footprint dataframe
    circular_footprint_df = grab_circular_reads(subsample, ref_seq_length)

    # Transform expanded footprint dataframe into count of footprint sizes dataframe
    footprint_count_dfs = prep_dfs_for_subtraction(circular_footprint_df)

    # Filter footprint_count_dfs for desired footprint sizes and positions
    filtered_footprint_count_df = filter_fp_df(footprint_count_dfs[0], bin_size, row_range, col_range) / subsample_size
    
    return filtered_footprint_count_df

def collect_single_sample_average(input_df, num_iterations=1000, subsample_size=5000, ref_seq_length=4718, row_range=None, col_range=None, bin_size=None):
    """
    Processes multiple single samples and collects the average footprint occupancy.
    
    Parameters:
    -----------
    input_df : pandas.DataFrame
        The output df from read_ft_data
    num_iterations : int, default=1000
        Number of iterations to perform
    subsample_size : int, default=5000
        Number of reads to include in each subsample
    ref_seq_length : int, default=4718
        Length of the DNA read to pass to grab_circular_reads
    row_range : tuple, default=None
        Range of rows to include in the subsample
    col_range : tuple, default=None
        Range of columns to include in the subsample
    bin_size : int, default=None
        Bin size to pass to filter_fp_df    
    
    Returns:
    --------
    A pandas DataFrame with the same structure as a filtered footprint count DataFrame,
    but where each cell contains the average value across all iterations.
    """
    # Use defaultdict to store lists of values for each (row, col) pair
    results = defaultdict(list)
    
    # Store the structure of the first df
    first_df = None
    
    # Generate samples for all iterations
    progress_bar = tqdm.tqdm(
        range(num_iterations),
        desc="Processing iterations",
        unit="iteration",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
    )
    
    for i in progress_bar:
        # Update progress bar description with current iteration
        if i > 0 and i % 10 == 0:
            progress_bar.set_description(f"Completed {i}/{num_iterations} iterations")
            
        # Process a single sample
        filtered_footprint_count_df = process_single_sample(
            input_df,
            subsample_size,
            ref_seq_length,
            row_range,
            col_range,
            bin_size
        )
        
        # Store the structure of the first df
        if first_df is None:
            first_df = filtered_footprint_count_df
        
        # Add the values to the results dictionary
        for row in filtered_footprint_count_df.index:
            for col in filtered_footprint_count_df.columns:
                results[(row, col)].append(filtered_footprint_count_df.at[row, col])
    
    # Create the final DataFrame
    result_df = pd.DataFrame(index=first_df.index, columns=first_df.columns)
    
    # Fill the final DataFrame with the average values
    for (row, col), values in results.items():
        result_df.at[row, col] = sum(values) / len(values)
    
    return result_df

def collect_single_sample_average_parallel(input_df, num_iterations=1000, subsample_size=5000, ref_seq_length=4718, 
                                          row_range=None, col_range=None, bin_size=None, num_processes=None):
    """
    Parallel version of collect_single_sample_average that distributes iterations across multiple processes.
    
    Parameters:
    -----------
    Same as collect_single_sample_average, plus:
    num_processes : int, optional
        Number of processes to use. If None, uses the number of CPU cores.
    
    Returns:
    --------
    Same as collect_single_sample_average
    """
    # Determine number of processes to use
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    
    # Limit processes to a reasonable number (to avoid system overload)
    num_processes = min(num_processes, multiprocessing.cpu_count(), num_iterations)
    
    print(f"Running with {num_processes} parallel processes")
    
    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Create a partial function with all parameters except the iteration index
        worker_func = partial(
            process_single_sample,
            input_df,
            subsample_size,
            ref_seq_length,
            row_range,
            col_range,
            bin_size
        )
        
        # Process all iterations in parallel with progress tracking
        iterations = list(range(num_iterations))
        sample_dfs = list(tqdm.tqdm(
            pool.imap(worker_func, iterations),
            total=num_iterations,
            desc=f"Processing iterations in parallel",
            unit="iteration",
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
            mininterval=30  # Update log every 30 seconds
        ))
    
    # Use the first df to determine structure
    first_df = sample_dfs[0]
    
    # Initialize dictionary to store lists of values
    results = defaultdict(list)
    
    # Collect values from all sample_dfs
    for df in sample_dfs:
        for row in df.index:
            for col in df.columns:
                results[(row, col)].append(df.at[row, col])
    
    # Create the final DataFrame
    result_df = pd.DataFrame(index=first_df.index, columns=first_df.columns)
    
    # Fill the final DataFrame with the average values
    for (row, col), values in results.items():
        result_df.at[row, col] = sum(values) / len(values)
    
    return result_df

def main():
    """
    Main function to execute the footprint analysis workflow.
    """
    start_time = time.time()
    
    # Create example command for help text
    example_text = '''
Example Usage:
-------------
python single_sample_average.py \\
  --sample_path /path/to/control_sample.bed \\
  --output results/footprint_average.pkl \\
  --ref_seq_length 4718 \\
  --subsample_size 5000 \\
  --num_iterations 1000 \\
  --row_min 10 \\
  --row_max 50 \\
  --col_min 100 \\
  --col_max 200 \\
  --bin_size 10 \\
  --num_processes 8 \\
  --parallel
'''
    
    parser = argparse.ArgumentParser(
        description='Generate average footprint distribution from BED file.',
        epilog=example_text,
        formatter_class=argparse.RawDescriptionHelpFormatter  # This preserves the formatting of the epilog
    )
    parser.add_argument('--sample_path', type=str, required=True, help='path to FiberHMM bed file')
    parser.add_argument('--output', type=str, required=True, help='Output file path')
    parser.add_argument('--ref_seq_length', type=int, default=4718, help='Reference sequence length')
    parser.add_argument('--subsample_size', type=int, default=5000, help='Size of each subsample (should be similar to experimental read count)')
    parser.add_argument('--num_iterations', type=int, default=100, help='Number of iterations')
    parser.add_argument('--row_min', type=int, default=None, help='Minimum row value (default: None)')
    parser.add_argument('--row_max', type=int, default=None, help='Maximum row value (default: None)')
    parser.add_argument('--col_min', type=int, default=None, help='Minimum column value (default: None)')
    parser.add_argument('--col_max', type=int, default=None, help='Maximum column value (default: None)')
    parser.add_argument('--bin_size', type=int, default=None, help='Bin size for footprint sizes (default: None)')
    parser.add_argument('--num_processes', type=int, default=None, 
                      help='Number of processes to use for parallel processing (default: number of CPU cores)')
    parser.add_argument('--parallel', action='store_true', 
                      help='Use parallel processing for faster computation')
    
    args = parser.parse_args()
    
    # Check if the bed file exists
    if not os.path.exists(args.sample_path):
        print(f"Error: The specified BED file '{args.sample_path}' does not exist.")
        sys.exit(1)
        
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Read the data
    print(f"Reading data from {args.sample_path}")
    bed_df = read_ft_data_single(args.sample_path)
    
    # Set row and column ranges
    row_range = None if args.row_min is None or args.row_max is None else (args.row_min, args.row_max)
    col_range = None if args.col_min is None or args.col_max is None else (args.col_min, args.col_max)
    
    print(f"Starting processing of {args.num_iterations} iterations with subsample size {args.subsample_size}")
    print(f"This may take a while. Progress will be displayed below:")
    
    # Determine whether to use parallel processing
    if args.parallel:
        results_df = collect_single_sample_average_parallel(
            bed_df,
            num_iterations=args.num_iterations,
            subsample_size=args.subsample_size,
            ref_seq_length=args.ref_seq_length,
            row_range=row_range,
            col_range=col_range,
            bin_size=args.bin_size,
            num_processes=args.num_processes
        )
    else:
        results_df = collect_single_sample_average(
            bed_df,
            num_iterations=args.num_iterations,
            subsample_size=args.subsample_size,
            ref_seq_length=args.ref_seq_length,
            row_range=row_range,
            col_range=col_range,
            bin_size=args.bin_size
        )
    
    # Save the results
    print(f"Saving results to {args.output}")
    results_df.to_pickle(args.output)
    
    # Calculate and display total runtime
    total_time = time.time() - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Processing complete! Total runtime: {int(hours)}h {int(minutes)}m {int(seconds)}s")

if __name__ == "__main__":
    # This guard is important for multiprocessing to work correctly
    main()