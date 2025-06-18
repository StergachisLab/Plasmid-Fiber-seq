#!/usr/bin/env python3
"""
Helper Functions for SNP Footprint Analysis

This module contains utility functions and helpers used across the project.
"""

import os
import time
import traceback
import multiprocessing as mp
from functools import partial

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
    # Import required modules
    from fiberseq_MPRA.core.file_io import load_analysis_dictionaries
    from fiberseq_MPRA.core.data_processing import compare_samples_with_dictionaries
    
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
    # Import required modules
    import os
    import numpy as np
    from functools import partial
    import multiprocessing as mp
    import time
    
    from fiberseq_MPRA.core.file_io import read_values_from_tsv, parse_bed_filename, count_bed_file_rows
    from fiberseq_MPRA.visualization.plots import plot_position_group, plot_position_groups_sequential
    
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
        pvalue = None
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
                    log2_fc_value = value_dict[key]['log2_fc']
                    pvalue = value_dict[key].get('pvalue')
                    print(f"  Found log2-FC value for {key}: {log2_fc_value}")
                    if pvalue is not None:
                        print(f"  Found P-value for {key}: {pvalue}")
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
            "log2_fc_value": log2_fc_value,  # log2 fold change value from TSV
            "pvalue": pvalue,               # P-value from TSV
            "is_merged": is_merged,         # Whether this is a merged file
            "second_original_pos": second_original_pos  # Second position for merged files
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
