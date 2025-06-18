#!/usr/bin/env python3
"""
Visualization Functions for SNP Footprint Analysis

This module contains functions for creating and saving heatmap plots and figures
for footprint enrichment data.
"""

import os
import traceback
import time
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial

from fiberseq_MPRA.core.file_io import format_sample_name

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
    
    Returns:
    tuple: (fig, ax) - the figure and axis objects for further customization
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
    apply_fdr : bool
        Whether FDR correction was applied
    fdr_threshold : float
        FDR threshold for significance
        
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

def parallel_plot_position_groups(position_groups, sample_info, output_directory, 
                                 vmin, vmax, cmap, xticks, plot_both, plot_index,
                                 n_processes=None, apply_fdr=False, fdr_threshold=0.05):
    """
    Plot position groups in parallel using multiprocessing.
    
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
    n_processes : int or None
        Number of processes to use for parallel processing
    apply_fdr : bool
        Whether FDR correction was applied
    fdr_threshold : float
        FDR threshold for significance
        
    Returns:
    --------
    int
        Number of successfully plotted position groups
    """
    # Determine number of processes to use
    if n_processes is None:
        n_processes = max(1, min(mp.cpu_count() - 1, len(position_groups)))
    else:
        n_processes = max(1, min(n_processes, mp.cpu_count(), len(position_groups)))
    
    print(f"Plotting {len(position_groups)} position groups using {n_processes} processes...")
    
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
        
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(plot_func, position_data)
        
        # Count successful plots
        successful_plots = sum(1 for _, success in results if success)
        
        plotting_time = time.time() - start_time
        print(f"Parallel plotting complete. Successfully created {successful_plots} out of {len(position_groups)} plots in {plotting_time:.2f} seconds.")
        
        return successful_plots
        
    except Exception as e:
        print(f"Error in parallel plotting: {str(e)}")
        print(traceback.format_exc())
        print("Falling back to sequential plotting...")
        
        # Fall back to sequential plotting
        return plot_position_groups_sequential(
            position_groups, sample_info, output_directory, 
            vmin, vmax, cmap, xticks, plot_both, plot_index,
            apply_fdr=apply_fdr, fdr_threshold=fdr_threshold
        )
