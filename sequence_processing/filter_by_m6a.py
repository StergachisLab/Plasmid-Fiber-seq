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
    verbose (bool): If True, set logging level to DEBUG, otherwise WARNING.
    """
    level = logging.DEBUG if verbose else logging.WARNING
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
        # Safely get read name for logging (handle case where query_name might not exist)
        read_name = getattr(fiber, 'query_name', 'unknown')
        logging.debug(f"Read {read_name} has no A or T bases, setting methylation to 0%")
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
    
    finally:
        # Clean up resources - pyft objects are automatically cleaned up when they go out of scope
        # No explicit close() methods are available for Fiberbam and Fiberwriter objects
        logging.debug("Cleaning up resources (automatic cleanup for pyft objects)")

    # Calculate and log final statistics
    if total_reads > 0:
        pass_percentage = (reads_written / total_reads) * 100
        logging.info(f"Filtering complete:")
        logging.info(f"  Total reads processed: {total_reads:,}")
        logging.info(f"  Reads passing filter: {reads_written:,}")
        logging.info(f"  Filter pass rate: {pass_percentage:.2f}%")
        
        # Print summary for user (in addition to logging)
        print(f'Total reads = {total_reads:,}, Reads passing filter = {reads_written:,}, '
              f'percentage of reads passing filter = {pass_percentage:.2f}%')
    
    return total_reads, reads_written


def create_command_info(args, chromosomes_list):
    """
    Create comprehensive command information including defaults for header tracking.
    
    Parameters:
    args: Parsed command line arguments
    chromosomes_list: Processed list of chromosomes (can be None)
    
    Returns:
    str: Formatted command information showing all parameters used
    """
    import datetime
    import getpass
    import socket
    
    # Start with the original command as entered
    original_command = ' '.join(sys.argv)
    
    # Create detailed parameter information
    param_info = []
    
    # Input/Output files (use basename for privacy/brevity)
    param_info.append(f"input={Path(args.input_bam_path).name}")
    param_info.append(f"output={Path(args.output_bam_path).name}")
    
    # Methylation range (show whether default or custom)
    if args.methylation_range == [0.1, 1.0]:
        param_info.append(f"methylation_range=0.1-1.0 (default)")
    else:
        param_info.append(f"methylation_range={args.methylation_range[0]}-{args.methylation_range[1]} (custom)")
    
    # Chromosomes (show whether all or filtered)
    if chromosomes_list is None:
        param_info.append("chromosomes=all (default)")
    else:
        chrom_str = ','.join(chromosomes_list[:3])  # Show first 3 chromosomes
        if len(chromosomes_list) > 3:
            chrom_str += f"... ({len(chromosomes_list)} total)"
        param_info.append(f"chromosomes={chrom_str} (filtered)")
    
    # Flags
    flags = []
    if args.verbose:
        flags.append("verbose")
    if args.skip_header:
        flags.append("skip_header")
    else:
        flags.append("enhanced_header")
    
    if flags:
        param_info.append(f"flags={','.join(flags)}")
    else:
        param_info.append("flags=none")
    
    # Add execution metadata
    try:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        param_info.append(f"timestamp={timestamp}")
    except:
        pass
    
    try:
        username = getpass.getuser()
        param_info.append(f"user={username}")
    except:
        pass
    
    try:
        hostname = socket.gethostname()
        param_info.append(f"host={hostname}")
    except:
        pass
    
    # Combine original command with parameter details
    detailed_info = f"{original_command} | PARAMS: {'; '.join(param_info)}"
    
    return detailed_info


def add_command_to_header_samtools(bam_path, command_info, make_optional=True):
    """
    Use samtools reheader to modify BAM header - much faster than manual rewriting.
    
    Parameters:
    bam_path (str): Path to the BAM file to modify
    command_info (str): Detailed command information including parameters to add to header
    make_optional (bool): If True, continue even if header modification fails
    """
    try:
        logging.info("Adding command to BAM header using samtools reheader")
        
        # Create temporary paths
        temp_header_path = str(Path(bam_path).with_suffix('.temp_header.sam'))
        temp_output_path = str(Path(bam_path).with_suffix('.temp.bam'))
        
        # Extract current header
        result = subprocess.run([
            "samtools", "view", "-H", bam_path
        ], capture_output=True, text=True, check=True)
        
        current_header = result.stdout
        
        # Parse and modify header
        header_lines = current_header.strip().split('\n')
        
        # Generate unique PG ID
        existing_pg_ids = []
        for line in header_lines:
            if line.startswith('@PG\tID:'):
                # Extract ID from @PG line
                for field in line.split('\t'):
                    if field.startswith('ID:'):
                        existing_pg_ids.append(field[3:])
                        break
        
        # Create unique PG ID
        base_id = 'm6a_filter'
        pg_id = base_id
        counter = 1
        while pg_id in existing_pg_ids:
            pg_id = f"{base_id}_{counter}"
            counter += 1
        
        # Create new PG line with comprehensive information
        # Escape any problematic characters in command_info for SAM format
        escaped_command = command_info.replace('\t', ' ').replace('\n', ' ')
        new_pg_line = f"@PG\tID:{pg_id}\tPN:filter_by_m6a.py\tVN:1.0\tCL:{escaped_command}"
        
        # Insert new PG line in the correct location (after existing @PG lines, or before @CO lines)
        modified_lines = []
        pg_added = False
        
        # Track where we are in the header structure
        past_sq_lines = False
        past_rg_lines = False
        
        for line in header_lines:
            # Check if we've passed the @SQ and @RG sections
            if line.startswith('@SQ'):
                past_sq_lines = True
            elif line.startswith('@RG'):
                past_rg_lines = True
            elif line.startswith('@CO') and not pg_added:
                # Insert before @CO lines if we haven't added PG yet
                modified_lines.append(new_pg_line)
                pg_added = True
            elif line.startswith('@PG'):
                # Add our PG line after the last existing @PG line
                modified_lines.append(line)
                if not pg_added:
                    modified_lines.append(new_pg_line)
                    pg_added = True
                continue
            elif past_sq_lines and past_rg_lines and not line.startswith(('@PG', '@CO')) and not pg_added:
                # If we're past @SQ and @RG but haven't seen @PG or @CO, add our PG line here
                modified_lines.append(new_pg_line)
                pg_added = True
            
            modified_lines.append(line)
        
        # If we still haven't added the PG line, add it at the end
        if not pg_added:
            modified_lines.append(new_pg_line)
        
        # Write modified header to temp file
        with open(temp_header_path, 'w') as f:
            f.write('\n'.join(modified_lines) + '\n')
        
        # Use samtools reheader (faster than manual BAM rewriting)
        with open(temp_output_path, 'wb') as outfile:
            subprocess.run([
                "samtools", "reheader", temp_header_path, bam_path
            ], stdout=outfile, check=True)
        
        # Replace original file atomically
        os.replace(temp_output_path, bam_path)
        
        # Regenerate index
        logging.info("Regenerating BAM index")
        subprocess.run(["samtools", "index", bam_path], check=True)
        
        # Clean up temp files
        if os.path.exists(temp_header_path):
            os.remove(temp_header_path)
        
        logging.info("Successfully updated BAM header using samtools reheader")
        
    except Exception as e:
        # Clean up temp files on error
        for temp_file in [temp_header_path, temp_output_path]:
            if 'temp_file' in locals() and os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except:
                    pass
        
        if make_optional:
            logging.warning(f"Failed to update BAM header using samtools reheader: {e}")
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

    # Optional argument to skip enhanced header modification
    parser.add_argument("--skip_header", action='store_true',
                        help="Skip adding enhanced command information to BAM header (basic info may still be added by pyft).")

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

        # Run the filter_by_m6a function
        total_reads, reads_written = filter_by_m6a(
            args.input_bam_path, 
            args.output_bam_path,
            chromosomes=chromosomes_list, 
            methylation_range=args.methylation_range
        )

        # Add enhanced command information to the BAM header (unless skipped)
        if not args.skip_header and reads_written > 0:
            # Create comprehensive command information for header tracking
            command_info = create_command_info(args, chromosomes_list)
            # Add filtering results to command info
            pass_rate = (reads_written / total_reads) * 100 if total_reads > 0 else 0
            command_info_with_stats = f"{command_info} | RESULTS: total_reads={total_reads:,}; passed={reads_written:,}; pass_rate={pass_rate:.2f}%"
            add_command_to_header_samtools(args.output_bam_path, command_info_with_stats, make_optional=True)
            logging.info("Added comprehensive command information to BAM header")
        else:
            logging.info("Skipped enhanced BAM header modification (basic info may still be added by pyft)")
        
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