"""
Command-line interface for fiberseq_MPRA package subcommands.
"""

import sys

def main():
    """Main entry point for fiberseq_MPRA package."""
    
    if len(sys.argv) < 2 or sys.argv[1] in ['--help', '-h', 'help']:
        print("Usage: python -m fiberseq_MPRA <command> [args]")
        print("\nAvailable commands:")
        print("  main             Run the main footprint analysis")
        print("  generate-control Generate control data from WT BED file")
        print("  generate-null    Generate null distribution statistics")
        print("\nFor help on a specific command:")
        print("  python -m fiberseq_MPRA <command> --help")
        print("\nExample workflow:")
        print("  1. python -m fiberseq_MPRA generate-control --input-bed WT.bed --output control.pkl")
        print("  2. python -m fiberseq_MPRA generate-null --input-bed WT.bed --output-means null_means.pkl --output-stds null_stds.pkl")
        print("  3. python -m fiberseq_MPRA main --input-dir variants/ --control-pkl control.pkl --null-means-pkl null_means.pkl --null-stds-pkl null_stds.pkl --output-dir results/")
        sys.exit(0)
    
    command = sys.argv[1]
    
    # Remove the command from sys.argv so submodules see normal arguments
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    
    if command == "main":
        from fiberseq_MPRA.main import main as main_analysis
        main_analysis()
    elif command == "generate-control":
        from fiberseq_MPRA.preprocessing.control_generator import main as generate_control
        generate_control()
    elif command == "generate-null":
        from fiberseq_MPRA.preprocessing.null_generator import main as generate_null
        generate_null()
    else:
        print(f"Unknown command: {command}")
        print("Use 'python -m fiberseq_MPRA --help' for available commands")
        sys.exit(1)

if __name__ == "__main__":
    main()
