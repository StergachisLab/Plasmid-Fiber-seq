"""
Fiberseq MPRA Analysis Package

A package for SNP footprint analysis with multiprocessing support.
"""

__version__ = '1.0.0'

# Try to ensure FiberHMM_functions is available, but don't fail if it's not
try:
    from .core.fiberhmm_integration import ensure_fiberhmm_in_path
    ensure_fiberhmm_in_path(verbose=False)
    print("FiberHMM integration initialized successfully")
except ImportError:
    # This is expected during testing - FiberHMM_functions might not be available
    import warnings
    warnings.warn("FiberHMM_functions not available - using mock functions for testing")
except Exception as e:
    # Don't fail the entire import if FiberHMM setup fails
    import warnings
    warnings.warn(f"FiberHMM integration warning: {e}")

# Make key functions easily accessible
try:
    from .utils.helpers import process_and_plot_grouped_samples
    from .core.file_io import parse_bed_filename, format_sample_name
    from .analysis.statistics import calculate_pvalues_one_tailed_dict
except ImportError:
    # Don't fail the entire package if individual modules have issues
    pass
