#!/usr/bin/env python3
"""
Quick integration test for FiberHMM functions
"""

def test_integration():
    print("Testing FiberHMM integration...")
    
    # Test 1: Package import
    try:
        import footprint_analysis
        print("✓ Package imported successfully")
    except ImportError as e:
        print(f"✗ Package import failed: {e}")
        return False
    
    # Test 2: FiberHMM integration
    try:
        from footprint_analysis.core.fiberhmm_integration import test_fiberhmm_import
        if test_fiberhmm_import():
            print("✓ FiberHMM integration working")
        else:
            print("✗ FiberHMM integration failed")
            return False
    except Exception as e:
        print(f"✗ FiberHMM integration error: {e}")
        return False
    
    # Test 3: Function imports
    try:
        from footprint_analysis.core.data_processing import grab_circular_reads, prep_dfs_for_subtraction, filter_fp_df
        from footprint_analysis.core.file_io import read_ft_data
        print("✓ All functions imported successfully")
    except ImportError as e:
        print(f"✗ Function import failed: {e}")
        return False
    
    # Test 4: CLI interface
    try:
        import subprocess
        import sys
        result = subprocess.run([
            sys.executable, '-m', 'footprint_analysis.main', '--help'
        ], capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0:
            print("✓ CLI interface working")
        else:
            print(f"✗ CLI interface failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"✗ CLI test failed: {e}")
        return False
    
    print("\n🎉 All integration tests passed!")
    print("Your package is ready to use with FiberHMM functions.")
    return True

if __name__ == "__main__":
    test_integration()
