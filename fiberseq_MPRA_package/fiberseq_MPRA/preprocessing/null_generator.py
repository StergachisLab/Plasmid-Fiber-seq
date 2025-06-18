#!/usr/bin/env python3
"""
Null Distribution Generator for Fiberseq MPRA Analysis Package
"""

import os
import sys
import argparse
import time
import pandas as pd
import numpy as np
import tqdm
import multiprocessing
import pickle
from functools import partial
from collections import defaultdict

# Import from fiberseq_MPRA package
from fiberseq_MPRA.core.data_processing import grab_circular_reads, prep_dfs_for_subtraction, filter_fp_df

def read_ft_data_single(bed_file):
    """Read a single BED file with the expected FiberHMM format."""
    try:
        df = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2, 3, 8, 9])
        df = df.reset_index(drop=True)
        df.columns = ['chrom', 'start', 'end', 'name', 'blockStarts', 'blockSizes']
        df = df.loc[df['blockStarts'] != '.']
        df['length'] = df['end'] - df['start']
        df = df.sort_values('length', ascending=False).drop_duplicates('name')
        return df
    except Exception as e:
        raise ValueError(f"Error reading BED file {bed_file}: {e}")

def get_subsampled_delta_df(input_df, subsample_size=5000, ref_seq_length=4718, 
                           row_range=None, col_range=None, bin_size=None):
    """Generate differential footprint data between two random subsamples."""
    random_state = np.random.RandomState()
    
    subsample_a = input_df.sample(n=subsample_size, replace=False, random_state=random_state)
    remaining_reads = input_df.drop(subsample_a.index)
    subsample_b = remaining_reads.sample(n=subsample_size, replace=False, random_state=random_state)
    
    circular_footprint_df_a = grab_circular_reads(subsample_a, ref_seq_length)
    circular_footprint_df_b = grab_circular_reads(subsample_b, ref_seq_length)

    footprint_count_dfs_a = prep_dfs_for_subtraction(circular_footprint_df_a)
    footprint_count_dfs_b = prep_dfs_for_subtraction(circular_footprint_df_b)

    filtered_a = filter_fp_df(footprint_count_dfs_a[0], bin_size, row_range, col_range) / subsample_size
    filtered_b = filter_fp_df(footprint_count_dfs_b[0], bin_size, row_range, col_range) / subsample_size

    return filtered_a - filtered_b

def process_single_iteration(bed_df, subsample_size, ref_seq_length, row_range, col_range, bin_size, iteration_idx=None):
    """Process a single iteration for parallel execution."""
    if iteration_idx is not None:
        seed = int(iteration_idx) + os.getpid() + int(time.time() * 1000) % 10000
        np.random.seed(seed)

    return get_subsampled_delta_df(bed_df, subsample_size, ref_seq_length, row_range, col_range, bin_size)

def generate_null_statistics_parallel(input_df, num_iterations=1000, subsample_size=5000, 
                                    ref_seq_length=4718, row_range=None, col_range=None, 
                                    bin_size=None, num_processes=None):
    """Generate null distribution statistics using parallel processing."""
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    
    num_processes = min(num_processes, multiprocessing.cpu_count(), num_iterations)
    
    print(f"Generating null statistics with {num_processes} parallel processes")
    print(f"Processing {num_iterations} iterations with subsample size {subsample_size}")
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        worker_func = partial(
            process_single_iteration, input_df, subsample_size, ref_seq_length, 
            row_range, col_range, bin_size
        )
        
        iterations = list(range(num_iterations))
        delta_dfs = list(tqdm.tqdm(
            pool.imap(worker_func, iterations),
            total=num_iterations,
            desc="Generating null distributions",
            unit="iteration"
        ))
    
    print("Computing statistics from null distributions...")
    
    first_df = delta_dfs[0]
    position_values = defaultdict(list)
    
    for delta_df in delta_dfs:
        for row in delta_df.index:
            for col in delta_df.columns:
                position_values[(row, col)].append(delta_df.at[row, col])
    
    means_dict = {}
    stds_dict = {}
    
    for (row, col), values in position_values.items():
        means_dict[(row, col)] = np.mean(values)
        stds_dict[(row, col)] = np.std(values)
    
    return means_dict, stds_dict, first_df.shape

def main():
    """Main function for null distribution statistics generation."""
    parser = argparse.ArgumentParser(description='Generate null distribution statistics for footprint analysis.')
    
    parser.add_argument('--input-bed', required=True, help='Path to control BED file')
    parser.add_argument('--output-means', required=True, help='Output path for null means pickle file')
    parser.add_argument('--output-stds', required=True, help='Output path for null stds pickle file')
    parser.add_argument('--ref-seq-length', type=int, default=4718, help='Reference sequence length')
    parser.add_argument('--subsample-size', type=int, default=5000, help='Number of reads per subsample')
    parser.add_argument('--num-iterations', type=int, default=1000, help='Number of iterations')
    parser.add_argument('--row-min', type=int, help='Minimum footprint size')
    parser.add_argument('--row-max', type=int, help='Maximum footprint size')
    parser.add_argument('--col-min', type=int, help='Minimum genomic position')
    parser.add_argument('--col-max', type=int, help='Maximum genomic position')
    parser.add_argument('--bin-size', type=int, help='Bin size for footprint sizes')
    parser.add_argument('--processes', type=int, help='Number of processes')
    parser.add_argument('--parallel', action='store_true', help='Use parallel processing')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_bed):
        print(f"Error: Input BED file '{args.input_bed}' does not exist.")
        sys.exit(1)
    
    print(f"Reading control BED file: {args.input_bed}")
    bed_df = read_ft_data_single(args.input_bed)
    print(f"Loaded {len(bed_df)} reads")
    
    min_required_reads = args.subsample_size * 2
    if len(bed_df) < min_required_reads:
        print(f"Error: Need at least {min_required_reads} reads for subsample size {args.subsample_size}")
        sys.exit(1)
    
    row_range = None if args.row_min is None or args.row_max is None else (args.row_min, args.row_max)
    col_range = None if args.col_min is None or args.col_max is None else (args.col_min, args.col_max)
    
    start_time = time.time()
    
    if args.parallel:
        means_dict, stds_dict, shape = generate_null_statistics_parallel(
            bed_df, args.num_iterations, args.subsample_size, args.ref_seq_length,
            row_range, col_range, args.bin_size, args.processes
        )
    else:
        print("Sequential processing not implemented yet. Use --parallel")
        sys.exit(1)
    
    print(f"Saving means dictionary to: {args.output_means}")
    with open(args.output_means, 'wb') as f:
        pickle.dump(means_dict, f)
    
    print(f"Saving stds dictionary to: {args.output_stds}")
    with open(args.output_stds, 'wb') as f:
        pickle.dump(stds_dict, f)
    
    total_time = time.time() - start_time
    print(f"Null distribution generation complete! Runtime: {total_time:.1f}s")
    print(f"Generated statistics for {len(means_dict):,} positions")

if __name__ == "__main__":
    main()
