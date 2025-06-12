#!/usr/bin/env python3

import argparse
import pyft
import pysam
import sys
import os
import subprocess
import logging
from pathlib import Path
from tqdm import tqdm


def setup_logging(verbose=False):
    """
    Set up logging configuration.
    
    Parameters:
    verbose (bool): If True, set logging level to DEBUG, otherwise INFO.
    """
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.CRITICAL + 1  # disables all logging output
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_dependencies():
    """
    Check if required external dependencies are available.
    
    Raises:
    SystemExit: If required dependencies are not found.
    """
    try:
        # Check if samtools is available
        result = subprocess.run(
            ["samtools", "--version"], 
            capture_output=True, 
            text=True, 
            check=True
        )
        logging.debug(f"samtools version: {result.stdout.split()[1]}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error("samtools not found. Please install samtools and ensure it's in your PATH.")
        sys.exit(1)


def validate_inputs(input_bam_path, output_bam_path, methylation_range, chromosomes_list):
    """
    Validate input parameters and file paths.
    
    Parameters:
    input_bam_path (str): Path to input BAM file
    output_bam_path (str): Path to output BAM file
    methylation_range (list): [min, max] methylation percentage range
    chromosomes_list (list): List of chromosomes to filter (can be None)
    
    Raises:
    SystemExit: If validation fails
    """
    # Check if input file exists
    input_path = Path(input_bam_path)
    if not input_path.exists():
        logging.error(f"Input BAM file '{input_bam_path}' not found")
        sys.exit(1)
    
    # Check if input file is readable
    if not input_path.is_file():
        logging.error(f"Input path '{input_bam_path}' is not a file")
        sys.exit(1)
    
    # Validate output directory exists
    output_path = Path(output_bam_path)
    if not output_path.parent.exists():
        logging.error(f"Output directory '{output_path.parent}' does not exist")
        sys.exit(1)
    
    # Check if output file already exists and warn user
    if output_path.exists():
        logging.warning(f"Output file '{output_bam_path}' already exists and will be overwritten")
    
    # Validate methylation range
    if len(methylation_range) != 2:
        logging.error("Methylation range must have exactly 2 values")
        sys.exit(1)
    
    if methylation_range[0] > methylation_range[1]:
        logging.error("Minimum methylation percentage must be <= maximum methylation percentage")
        sys.exit(1)
    
    if methylation_range[0] < 0 or methylation_range[1] > 1:
        logging.error("Methylation percentages must be between 0 and 1")
        sys.exit(1)
    
    # Log validation success
    logging.info("Input validation passed")
    if chromosomes_list:
        logging.info(f"Will filter for chromosomes: {', '.join(chromosomes_list)}")
    logging.info(f"Methylation range: {methylation_range[0]:.2f} - {methylation_range[1]:.2f}")


def count_reads(input_bam_path):
    """
    Counts the total number of reads in a BAM file using samtools.
    
    Parameters:
    input_bam_path (str): Path to the input BAM file
    
    Returns:
    int: Total number of reads in the BAM file
    """
    try:
        # Use samtools view -c to count the total number of reads in the BAM file
        result = subprocess.run(
            ["samtools", "view", "-c", input_bam_path],
            capture_output=True,
            text=True,
            check=True
        )
        read_count = int(result.stdout.strip())
        logging.info(f"Total reads in input BAM: {read_count:,}")
        return read_count
    except subprocess.CalledProcessError as e:
        logging.error(f"Error counting reads with samtools: {e.stderr}")
        return 0
    except ValueError as e:
        logging.error(f"Error parsing read count: {e}")
        return 0


def calculate_methylation_percentage(fiber):
    """
    Calculate the percentage of methylated AT basepairs for a given fiber.
    
    Parameters:
    fiber: A fiber object from pyft
    
    Returns:
    float: Percentage of methylated AT basepairs (0.0 to 1.0)
    """
    # Count total A and T bases in the sequence
    at_count = fiber.seq.count('A') + fiber.seq.count('T')
    
    # Handle edge case where there are no A or T bases
    if at_count == 0:
        logging.debug(f"Read {fiber.query_name} has no A or T bases, setting methylation to 0%")
        return 0.0
    
    # Calculate methylation percentage
    m6a_count = len(fiber.m6a.starts)
    percent_methylation = m6a_count / at_count
    
    return percent_methylation


def passes_filters(fiber, chromosomes, methylation_range):
    """
    Check if a fiber passes the specified filters.
    
    Parameters:
    fiber: A fiber object from pyft
    chromosomes (list or None): List of chromosomes to include, None for all
    methylation_range (list): [min, max] methylation percentage range
    
    Returns:
    bool: True if fiber passes all filters, False otherwise
    """
    # Calculate methylation percentage safely
    percent_methylation = calculate_methylation_percentage(fiber)
    
    # Check methylation range filter
    if not (methylation_range[0] <= percent_methylation <= methylation_range[1]):
        return False
    
    # Check chromosome filter (if specified)
    if chromosomes is not None and fiber.chrom not in chromosomes:
        return False
    
    return True


def filter_by_m6a(input_bam_path, output_bam_path, chromosomes=None, methylation_range=[0.1, 1]):
    """
    Filters reads based on the percent of methylated AT basepairs and writes the filtered reads to a new BAM file.

    Parameters:
    input_bam_path (str): Path to the input BAM file.
    output_bam_path (str): Path (including new file name) to the output BAM file.
    chromosomes (list): A list of chromosomes you want included in the output BAM file. 
                       The default is None, which means all chromosomes in the input BAM will be included.
    methylation_range (list): A list with two elements specifying the lower and upper bounds 
                             for the methylation percentage. Default is [0.1, 1], meaning it will 
                             include reads with at least 10% and up to 100% of methylated AT basepairs.

    Returns:
    tuple: (total_reads, reads_written) - counts of total and filtered reads
    """
    logging.info("Starting BAM filtering process")
    
    # Count total reads for progress tracking
    total_reads = count_reads(input_bam_path)
    if total_reads == 0:
        logging.error("No reads found in input BAM file or error occurred")
        return 0, 0

    # Initialize objects for reading and writing
    fiberbam = None
    writer = None
    reads_written = 0
    reads_processed = 0
    
    try:
        # Read in input bam as a Fiberbam object
        logging.info("Opening input BAM file")
        fiberbam = pyft.Fiberbam(input_bam_path)

        # Open a Fiberwriter object to write the filtered reads to the output BAM file
        logging.info("Creating output BAM file")
        writer = pyft.Fiberwriter(output_bam_path, input_bam_path)

        # Process each read with progress tracking
        logging.info("Processing reads...")
        for fiber in tqdm(fiberbam, total=total_reads, desc="Filtering reads"):
            reads_processed += 1
            
            # Check if fiber passes all filters
            if passes_filters(fiber, chromosomes, methylation_range):
                writer.write(fiber)
                reads_written += 1
            
            # Log progress for very large files (every 100k reads)
            if reads_processed % 100000 == 0:
                logging.debug(f"Processed {reads_processed:,} reads, written {reads_written:,}")

    except Exception as e:
        logging.error(f"Error during BAM filtering: {e}")
        raise

    # Calculate and log final statistics
    if total_reads > 0:
        pass_percentage = (reads_written / total_reads) * 100
        logging.info(f"Filtering complete:")
        logging.info(f"  Total reads processed: {total_reads:,}")
        logging.info(f"  Reads passing filter: {reads_written:,}")
        logging.info(f"  Filter pass rate: {pass_percentage:.2f}%")
        
        # Print summary for user (in addition to logging)
        print(f'Total reads = {total_reads:,}')
        print(f'Reads passing filter = {reads_written:,}')
        print(f'Percentage of reads passing filter = {pass_percentage:.2f}%')
    
    return total_reads, reads_written


def add_command_to_header(bam_path, command, make_optional=True):
    """
    Opens the BAM file using pysam and adds a tag to the header with the command used to filter the reads.
    
    Parameters:
    bam_path (str): Path to the BAM file to modify
    command (str): Command line string to add to header
    make_optional (bool): If True, continue even if header modification fails
    """
    try:
        logging.info("Adding command to BAM header")
        
        # Create temporary file path
        temp_output_path = str(Path(bam_path).with_suffix('.temp.bam'))
        
        # Read the existing header
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            header = bam_file.header.to_dict()

        # Ensure 'PG' (program) section exists in the header
        if 'PG' not in header:
            header['PG'] = []

        # Create a unique program ID to avoid conflicts
        pg_id = 'm6a_filter'
        
        # Add the new command to the 'PG' section
        header['PG'].append({
            'ID': pg_id,
            'PN': 'filter_by_m6a.py',  # Program name
            'VN': '1.0',              # Version
            'CL': command             # Command line
        })

        # Write new BAM file with updated header
        with pysam.AlignmentFile(temp_output_path, "wb", header=header) as out_bam:
            with pysam.AlignmentFile(bam_path, "rb") as bam_file:
                for read in bam_file:
                    out_bam.write(read)

        # Replace original file with the updated BAM file
        os.replace(temp_output_path, bam_path)
        
        # Generate index for the new BAM file
        logging.info("Indexing output BAM file")
        pysam.index(bam_path)
        
        logging.info("Successfully updated BAM header and created index")
        
    except Exception as e:
        # Clean up temporary file if it exists
        if os.path.exists(temp_output_path):
            try:
                os.remove(temp_output_path)
            except:
                pass
        
        if make_optional:
            logging.warning(f"Failed to update BAM header (continuing anyway): {e}")
        else:
            logging.error(f"Failed to update BAM header: {e}")
            raise


def main():
    """Main function to handle command-line interface and coordinate the filtering process."""
    
    # Set up argument parser with detailed help
    parser = argparse.ArgumentParser(
        description="Filter reads based on m6A methylation percentage and output BAM with filtered reads.",
        epilog="Example command: python filter_by_m6a.py input.bam output.bam --chromosomes chr1,chr2 --methylation_range 0.1 1 --verbose",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Required arguments for input and output BAM paths
    parser.add_argument("input_bam_path", type=str,
                        help="Path to the input BAM file.")
    parser.add_argument("output_bam_path", type=str,
                        help="Path to the output BAM file.")

    # Optional argument for chromosome filtering
    parser.add_argument("--chromosomes", type=str,
                        help="Comma-separated list of chromosomes to include in the output BAM, "
                             "or leave blank to include all reads in the input BAM.")

    # Optional argument for methylation range filtering
    parser.add_argument("--methylation_range", nargs=2, type=float,
                        default=[0.1, 1.0], metavar=('MIN', 'MAX'), 
                        help="Methylation percentage range (default: 0.1 to 1.0). "
                             "Values should be between 0 and 1.")

    # Optional argument to skip header modification
    parser.add_argument("--skip_header", action='store_true',
                        help="Skip adding command information to BAM header (faster for large files).")

    # Optional argument for verbose logging
    parser.add_argument("--verbose", "-v", action='store_true',
                        help="Enable verbose logging output.")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Set up logging based on verbosity level
    setup_logging(args.verbose)
    
    logging.info("Starting m6A BAM filtering script")
    logging.info(f"Input: {args.input_bam_path}")
    logging.info(f"Output: {args.output_bam_path}")

    try:
        # Check dependencies first
        check_dependencies()

        # Convert comma-separated chromosomes string to a list if provided, else set to None
        chromosomes_list = None
        if args.chromosomes:
            chromosomes_list = [chrom.strip() for chrom in args.chromosomes.split(',')]
            # Remove empty strings that might result from extra commas
            chromosomes_list = [chrom for chrom in chromosomes_list if chrom]

        # Validate all inputs
        validate_inputs(args.input_bam_path, args.output_bam_path, 
                       args.methylation_range, chromosomes_list)

        # Capture the full command-line input as a string for header
        command_line_input = ' '.join(sys.argv)

        # Run the filter_by_m6a function
        total_reads, reads_written = filter_by_m6a(
            args.input_bam_path, 
            args.output_bam_path,
            chromosomes=chromosomes_list, 
            methylation_range=args.methylation_range
        )

        # Add the command-line input to the BAM header (unless skipped)
        if not args.skip_header and reads_written > 0:
            add_command_to_header(args.output_bam_path, command_line_input, make_optional=True)
        elif args.skip_header:
            logging.info("Skipped BAM header modification as requested")
        
        # Final success message
        if reads_written > 0:
            logging.info("BAM filtering completed successfully!")
        else:
            logging.warning("No reads passed the filtering criteria")
            
    except KeyboardInterrupt:
        logging.info("Process interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        if args.verbose:
            logging.exception("Full traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()