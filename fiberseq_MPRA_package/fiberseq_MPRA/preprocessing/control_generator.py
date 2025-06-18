#!/usr/bin/env python3
"""
Control Data Generator for Fiberseq MPRA Analysis Package
"""

import os
import sys
import argparse
import time
import pandas as pd
import numpy as np
import tqdm
import multiprocessing
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

def process_single_sample(input_df, subsample_size, ref_seq_length, row_range, col_range, bin_size, iteration_idx=None):
    """Process a single sample iteration for control data generation."""
    if iteration_idx is not None:
        seed = int(iteration_idx) + os.getpid() + int(time.time() * 1000) % 10000
        np.random.seed(seed)
        random_state = np.random.RandomState(seed)
    else:
        random_state = np.random.RandomState()

    subsample = input_df.sample(n=subsample_size, replace=False, random_state=random_state)
    circular_footprint_df = grab_circular_reads(subsample, ref_seq_length)
    footprint_count_dfs = prep_dfs_for_subtraction(circular_footprint_df)
    filtered_footprint_count_df = filter_fp_df(footprint_count_dfs[0], bin_size, row_range, col_range)
    
    return filtered_footprint_count_df / subsample_size

def generate_control_data_parallel(input_df, num_iterations=1000, subsample_size=5000, 
                                 ref_seq_length=4718, row_range=None, col_range=None, 
                                 bin_size=None, num_processes=None):
    """Generate control data using parallel processing."""
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    
    num_processes = min(num_processes, multiprocessing.cpu_count(), num_iterations)
    
    print(f"Generating control data with {num_processes} parallel processes")
    print(f"Processing {num_iterations} iterations with subsample size {subsample_size}")
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        worker_func = partial(
            process_single_sample, input_df, subsample_size, ref_seq_length, 
            row_range, col_range, bin_size
        )
        
        iterations = list(range(num_iterations))
        sample_dfs = list(tqdm.tqdm(
            pool.imap(worker_func, iterations),
            total=num_iterations,
            desc="Generating control data",
            unit="iteration"
        ))
    
    first_df = sample_dfs[0]
    results = defaultdict(list)
    
    for df in sample_dfs:
        for row in df.index:
            for col in df.columns:
                results[(row, col)].append(df.at[row, col])
    
    result_df = pd.DataFrame(index=first_df.index, columns=first_df.columns)
    
    for (row, col), values in results.items():
        result_df.at[row, col] = sum(values) / len(values)
    
    return result_df

def main():
    """Main function for control data generation."""
    parser = argparse.ArgumentParser(description='Generate control data for footprint analysis from WT BED file.')
    
    parser.add_argument('--input-bed', required=True, help='Path to WT BED file')
    parser.add_argument('--output', required=True, help='Output path for control pickle file')
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
    
    print(f"Reading WT BED file: {args.input_bed}")
    bed_df = read_ft_data_single(args.input_bed)
    print(f"Loaded {len(bed_df)} reads")
    
    row_range = None if args.row_min is None or args.row_max is None else (args.row_min, args.row_max)
    col_range = None if args.col_min is None or args.col_max is None else (args.col_min, args.col_max)
    
    start_time = time.time()
    
    if args.parallel:
        control_df = generate_control_data_parallel(
            bed_df, args.num_iterations, args.subsample_size, args.ref_seq_length,
            row_range, col_range, args.bin_size, args.processes
        )
    else:
        print("Sequential processing not implemented yet. Use --parallel")
        sys.exit(1)
    
    print(f"Saving control data to: {args.output}")
    control_df.to_pickle(args.output)
    
    total_time = time.time() - start_time
    print(f"Control data generation complete! Runtime: {total_time:.1f}s")

if __name__ == "__main__":
    main()
