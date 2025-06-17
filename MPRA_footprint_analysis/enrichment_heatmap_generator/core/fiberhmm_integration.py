#!/usr/bin/env python3
"""
FiberHMM Integration Helper

This module handles the integration between the footprint_analysis package
and the FiberHMM_functions module.
"""

import sys
import os

def ensure_fiberhmm_in_path():
    """
    Ensure FiberHMM_functions is available for import.
    This function should be called before importing FiberHMM_functions.
    """
    # Get the directory containing this file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Look for FiberHMM_functions.py in various locations
    possible_locations = [
        current_dir,  # Same directory as this file
        os.path.dirname(current_dir),  # Parent directory (package root)
        os.path.dirname(os.path.dirname(current_dir)),  # Two levels up (enrichment_heatmap_generator)
        os.path.dirname(os.path.dirname(os.path.dirname(current_dir))),  # Three levels up (MPRA_footprint_analysis)
    ]
    
    for location in possible_locations:
        fiberhmm_path = os.path.join(location, 'FiberHMM_functions.py')
        if os.path.exists(fiberhmm_path):
            if location not in sys.path:
                sys.path.insert(0, location)
            print(f"Found FiberHMM_functions.py at: {fiberhmm_path}")
            return True
    
    # If not found, raise an informative error
    raise ImportError(
        "FiberHMM_functions.py not found. Please ensure it's in one of these locations:\n" +
        "\n".join(f"  - {loc}" for loc in possible_locations)
    )

def test_fiberhmm_import():
    """
    Test that FiberHMM_functions can be imported and has the required functions.
    """
    try:
        ensure_fiberhmm_in_path()
        import FiberHMM_functions as fhmm
        
        required_functions = [
            'grab_circular_reads',
            'prep_dfs_for_subtraction', 
            'filter_fp_df',
            'read_ft_data'
        ]
        
        missing_functions = []
        for func_name in required_functions:
            if not hasattr(fhmm, func_name):
                missing_functions.append(func_name)
        
        if missing_functions:
            raise ImportError(f"FiberHMM_functions is missing required functions: {missing_functions}")
        
        print("✓ FiberHMM_functions imported successfully with all required functions")
        return True
        
    except ImportError as e:
        print(f"✗ Failed to import FiberHMM_functions: {e}")
        return False
