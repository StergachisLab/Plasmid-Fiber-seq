#!/usr/bin/env python3
"""
File I/O Functions for SNP Footprint Analysis

This module contains functions for reading and writing files, including BED files,
pickle files, and TSV files with genomic information.
"""

import os
import sys
import pickle
import pandas as pd
import traceback

def read_ft_data(directory_path):
    """
    Read all BED files in a directory into a dictionary.
    
    This uses the FiberHMM_functions.read_ft_data function directly.
    
    Parameters:
    -----------
    directory_path : str
        Path to the directory containing BED files
        
    Returns:
    --------
    dict
        Dictionary mapping sample names to BED dataframes
    """
    # Ensure FiberHMM_functions is available
    from .fiberhmm_integration import ensure_fiberhmm_in_path
    ensure_fiberhmm_in_path()
    
    # Import the actual function from FiberHMM_functions
    import FiberHMM_functions as fhmm
    
    # The function signatures match perfectly
    return fhmm.read_ft_data(directory_path)

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
                import warnings
                warnings.warn("Precomputed null statistics have different dimensions than control dataframe. This may cause issues.")
    except Exception as e:
        raise Exception(f"Failed to load precomputed statistics: {e}")
    
    print("All data loaded successfully")
    return control_df, null_means_dict, null_stds_dict

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
                        
                    except (ValueError, TypeError):
                        # Skip invalid rows
                        continue
            print(f"Successfully read {len(value_dict)} entries from TSV file")
            return value_dict
        
    except Exception as e:
        print(f"Error reading TSV file {tsv_file_path}: {str(e)}")
        return {}
