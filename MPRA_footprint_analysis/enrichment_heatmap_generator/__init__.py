"""
Grouped Heatmap Generator for SNP Footprint Analysis with Multiprocessing

This package contains modules for processing BED files containing SNP data,
and creating comparison heatmaps for SNP footprint analysis.
"""

__version__ = '1.0.0'

# Ensure FiberHMM_functions is available
try:
    from .core.fiberhmm_integration import ensure_fiberhmm_in_path
    ensure_fiberhmm_in_path()
    print("FiberHMM integration initialized successfully")
except ImportError as e:
    import warnings
    warnings.warn(f"Could not setup FiberHMM_functions integration: {e}")
except Exception as e:
    import warnings
    warnings.warn(f"FiberHMM integration warning: {e}")
