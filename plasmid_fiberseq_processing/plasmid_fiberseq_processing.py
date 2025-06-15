#!/usr/bin/env python3

import os
import argparse
import pyft
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import ast
import json
import subprocess
import logging
import sys
import shlex
import platform
import datetime
from pathlib import Path
from tqdm import tqdm
from typing import Dict, List, Tuple, Optional


def extract_sample_name_from_bam(bam_path: str) -> str:
    """
    Extract sample name from BAM file path.
    
    Args:
        bam_path: Path to the BAM file
        
    Returns:
        Sample name extracted from filename
    """
    # Get the filename without directory path
    filename = os.path.basename(bam_path)
    
    # Remove .bam extension (and .sam if present)
    sample_name = filename
    if sample_name.endswith('.bam'):
        sample_name = sample_name[:-4]
    elif sample_name.endswith('.sam'):
        sample_name = sample_name[:-4]
    
    # Log the extracted sample name
    logging.info(f"Extracted sample name from BAM filename: '{sample_name}'")
    
    return sample_name


def setup_logging(verbose=False):
    """Set up logging configuration."""
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_dependencies():
    """Check if required external dependencies are available."""
    try:
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


def extract_reference_from_bam_header(bam_path: str) -> Optional[str]:
    """
    Extract reference FASTA path from BAM header @PG tags.
    
    Returns:
        Path to reference FASTA file if found, None otherwise
    """
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            header = bam_file.header.to_dict()
            
            # Look for @PG tags
            if 'PG' not in header:
                logging.warning("No @PG tags found in BAM header")
                return None
            
            # Search through all PG entries for alignment commands
            for pg_entry in header['PG']:
                if 'CL' not in pg_entry:
                    continue
                    
                command_line = pg_entry['CL']
                program_name = pg_entry.get('PN', '')
                
                # Look for common alignment tools
                if any(aligner in program_name.lower() for aligner in ['pbmm2', 'minimap2', 'bwa', 'bowtie']):
                    logging.info(f"Found alignment command from {program_name}: {command_line}")
                    
                    # Parse command line to extract reference path
                    reference_path = parse_reference_from_command(command_line, program_name)
                    if reference_path:
                        return reference_path
            
            logging.warning("No alignment commands found in BAM header @PG tags")
            return None
            
    except Exception as e:
        logging.error(f"Error reading BAM header: {e}")
        return None


def parse_reference_from_command(command_line: str, program_name: str) -> Optional[str]:
    """
    Parse reference FASTA path from alignment command line.
    
    Args:
        command_line: The command line string from @PG CL tag
        program_name: The program name from @PG PN tag
    
    Returns:
        Path to reference FASTA if found, None otherwise
    """
    try:
        # Split command line respecting quotes
        parts = shlex.split(command_line)
        
        if program_name.lower() == 'pbmm2':
            # pbmm2 format: pbmm2 align [options] reference.fa input.bam output.bam
            # Find position of 'align' subcommand
            if 'align' in parts:
                align_idx = parts.index('align')
                # Reference is typically the first positional argument after 'align' and any options
                for i in range(align_idx + 1, len(parts)):
                    arg = parts[i]
                    # Skip options (start with -)
                    if arg.startswith('-'):
                        # Skip option and its value if it takes one
                        continue
                    # First non-option argument should be reference
                    if arg.endswith(('.fa', '.fasta', '.fna')):
                        logging.info(f"Extracted reference from pbmm2 command: {arg}")
                        return arg
                        
        elif program_name.lower() in ['minimap2']:
            # minimap2 format: minimap2 [options] reference.fa input.fastq > output.sam
            for i, arg in enumerate(parts[1:], 1):  # Skip program name
                if not arg.startswith('-') and arg.endswith(('.fa', '.fasta', '.fna')):
                    logging.info(f"Extracted reference from minimap2 command: {arg}")
                    return arg
                    
        elif program_name.lower() in ['bwa']:
            # bwa format: bwa mem [options] reference.fa input.fastq
            if 'mem' in parts:
                mem_idx = parts.index('mem')
                for i in range(mem_idx + 1, len(parts)):
                    arg = parts[i]
                    if not arg.startswith('-') and arg.endswith(('.fa', '.fasta', '.fna')):
                        logging.info(f"Extracted reference from bwa command: {arg}")
                        return arg
        
        # Generic fallback - look for any .fa/.fasta/.fna file in command
        for arg in parts:
            if not arg.startswith('-') and arg.endswith(('.fa', '.fasta', '.fna')):
                logging.info(f"Extracted reference from generic parsing: {arg}")
                return arg
                
        logging.warning(f"Could not parse reference from command: {command_line}")
        return None
        
    except Exception as e:
        logging.error(f"Error parsing command line '{command_line}': {e}")
        return None


def should_skip_chromosome(chrom_name: str) -> bool:
    """
    Check if a chromosome should be skipped based on naming patterns.
    
    Args:
        chrom_name: Name of the chromosome
        
    Returns:
        True if chromosome should be skipped, False otherwise
    """
    return chrom_name.startswith('chrUn') or chrom_name.endswith('_random')


def classify_chromosome(chrom_name: str) -> str:
    """
    Classify a chromosome as either 'genomic' or 'plasmid' based on naming patterns.
    
    Args:
        chrom_name: Name of the chromosome
        
    Returns:
        'genomic' for standard genomic chromosomes, 'plasmid' for others
    """
    # Convert to lowercase for case-insensitive matching
    chrom_lower = chrom_name.lower()
    
    # Standard human chromosomes (including sex chromosomes and mitochondrial)
    genomic_patterns = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrx', 'chry', 'chrm', 'chrmt', 'chrebv'
    ]
    
    # Check for exact matches first (case-insensitive)
    if chrom_lower in genomic_patterns:
        logging.debug(f"Chromosome {chrom_name} classified as genomic (exact match: {chrom_lower})")
        return 'genomic'
    
    # Check for chromosome names that start with standard patterns but have suffixes
    # (e.g., chr1_alt, chr2_patch, etc.)
    for pattern in genomic_patterns:
        if chrom_lower.startswith(pattern + '_'):
            logging.debug(f"Chromosome {chrom_name} classified as genomic (prefix match: {pattern})")
            return 'genomic'
    
    # Check for non-chr prefixed standard chromosomes (1, 2, 3, etc.)
    if chrom_name.isdigit() and 1 <= int(chrom_name) <= 22:
        logging.debug(f"Chromosome {chrom_name} classified as genomic (numeric autosome)")
        return 'genomic'
    
    # Check for sex chromosomes without chr prefix (case-insensitive)
    if chrom_lower in ['x', 'y', 'm', 'mt', 'ebv']:
        logging.debug(f"Chromosome {chrom_name} classified as genomic (sex/mito/viral chromosome)")
        return 'genomic'
    
    # Check for common alternative naming patterns
    if chrom_lower.startswith('chromosome'):
        # Handle "chromosome1", "chromosome_1", etc.
        remaining = chrom_lower.replace('chromosome', '').lstrip('_')
        if remaining.isdigit() and 1 <= int(remaining) <= 22:
            logging.debug(f"Chromosome {chrom_name} classified as genomic (chromosome prefix)")
            return 'genomic'
        if remaining.lower() in ['x', 'y', 'm', 'mt', 'ebv']:
            logging.debug(f"Chromosome {chrom_name} classified as genomic (chromosome prefix sex/mito/viral)")
            return 'genomic'
    
    # Everything else is considered plasmid
    logging.debug(f"Chromosome {chrom_name} classified as plasmid (no genomic pattern match)")
    return 'plasmid'


def parse_fasta_file(fasta_path: str) -> Dict[str, int]:
    """Parse FASTA file and extract chromosome names and their lengths."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    reference_length_dict = {}
    skipped_chromosomes = []
    genomic_chromosomes = []
    plasmid_chromosomes = []
    
    try:
        with pysam.FastaFile(fasta_path) as fasta:
            for reference in fasta.references:
                # Skip chromosomes that start with "chrUn" or end with "_random"
                if should_skip_chromosome(reference):
                    skipped_chromosomes.append(reference)
                    continue
                    
                length = fasta.get_reference_length(reference)
                reference_length_dict[reference] = length
                
                # Classify chromosome for logging
                chrom_type = classify_chromosome(reference)
                if chrom_type == 'genomic':
                    genomic_chromosomes.append(reference)
                else:
                    plasmid_chromosomes.append(reference)
                
                logging.info(f"Found {chrom_type} chromosome: {reference} (length: {length:,} bp)")
        
        # Log summary of chromosome types
        if genomic_chromosomes:
            logging.info(f"Found {len(genomic_chromosomes)} genomic chromosomes: {', '.join(genomic_chromosomes)}")
        if plasmid_chromosomes:
            logging.info(f"Found {len(plasmid_chromosomes)} plasmid chromosomes: {', '.join(plasmid_chromosomes)}")
        
        if skipped_chromosomes:
            chrUn_count = sum(1 for chrom in skipped_chromosomes if chrom.startswith('chrUn'))
            random_count = sum(1 for chrom in skipped_chromosomes if chrom.endswith('_random'))
            
            skip_msg = []
            if chrUn_count > 0:
                skip_msg.append(f"{chrUn_count} chrUn chromosomes")
            if random_count > 0:
                skip_msg.append(f"{random_count} _random chromosomes")
            
            logging.info(f"Skipped {len(skipped_chromosomes)} chromosomes ({', '.join(skip_msg)}): " +
                        f"{', '.join(skipped_chromosomes[:5])}" + 
                        (f" and {len(skipped_chromosomes)-5} more..." if len(skipped_chromosomes) > 5 else ""))
        
        if not reference_length_dict:
            raise ValueError("No valid sequences found in FASTA file (all were chrUn/_random or empty)")
            
        logging.info(f"Successfully parsed {len(reference_length_dict)} chromosomes from FASTA file")
        return reference_length_dict
        
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Error parsing FASTA file: {str(e)}")


def parse_reference_length_dict(dict_str: str) -> Dict[str, int]:
    """Parse reference length dictionary from string format (legacy support)."""
    try:
        reference_length_dict = ast.literal_eval(dict_str)
        if not isinstance(reference_length_dict, dict):
            raise ValueError(
                "The reference length dictionary must be a valid Python dictionary.")
        return reference_length_dict
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"Invalid dictionary format: {str(e)}")


def check_bam_chromosomes(input_bam_path: str, reference_chromosomes: set) -> Tuple[set, set]:
    """
    Check which chromosomes are present in the BAM file and compare with reference.
    
    Returns:
        Tuple of (bam_chromosomes, missing_from_reference)
    """
    try:
        with pysam.AlignmentFile(input_bam_path, "rb") as bam_file:
            bam_chromosomes = set(bam_file.references)
            
        # Filter out chrUn and _random chromosomes from BAM chromosomes for comparison
        bam_chromosomes_filtered = {chrom for chrom in bam_chromosomes if not should_skip_chromosome(chrom)}
        skipped_chromosomes = {chrom for chrom in bam_chromosomes if should_skip_chromosome(chrom)}
        
        if skipped_chromosomes:
            chrUn_count = sum(1 for chrom in skipped_chromosomes if chrom.startswith('chrUn'))
            random_count = sum(1 for chrom in skipped_chromosomes if chrom.endswith('_random'))
            
            skip_msg = []
            if chrUn_count > 0:
                skip_msg.append(f"{chrUn_count} chrUn")
            if random_count > 0:
                skip_msg.append(f"{random_count} _random")
            
            logging.info(f"Found {len(skipped_chromosomes)} chromosomes in BAM file that will be ignored ({', '.join(skip_msg)})")
            
        missing_from_reference = bam_chromosomes_filtered - reference_chromosomes
        
        logging.info(f"Found {len(bam_chromosomes_filtered)} processable chromosomes in BAM file")
        
        if missing_from_reference:
            logging.warning(f"Found {len(missing_from_reference)} chromosome(s) in BAM file that are NOT in the reference FASTA:")
            for chrom in sorted(missing_from_reference):
                logging.warning(f"  - {chrom}")
            logging.warning("These chromosomes will be IGNORED during processing.")
            logging.warning("Consider updating your FASTA file if these should be included.")
            
            # Also print warning for user visibility
            print(f"\n⚠️  WARNING: Found chromosomes in BAM file not present in reference FASTA:")
            for chrom in sorted(missing_from_reference):
                print(f"    - {chrom}")
            print(f"   These {len(missing_from_reference)} chromosome(s) will be ignored during processing.")
            print(f"   Consider updating your FASTA file if these should be included.\n")
        
        # Check for chromosomes in reference but not in BAM
        missing_from_bam = reference_chromosomes - bam_chromosomes_filtered
        if missing_from_bam:
            logging.info(f"Reference chromosomes not found in BAM file ({len(missing_from_bam)}):")
            for chrom in sorted(missing_from_bam):
                logging.info(f"  - {chrom}")
            logging.info("This is normal if your BAM file doesn't contain reads for all references.")
        
        return bam_chromosomes_filtered, missing_from_reference
        
    except Exception as e:
        logging.error(f"Error checking BAM chromosomes: {e}")
        return set(), set()


def validate_inputs(input_bam_path: str, output_folder: str, 
                   genomic_methylation_range: List[float], plasmid_methylation_range: List[float], 
                   reference_length_dict: Dict[str, int], fasta_file: Optional[str] = None):
    """Validate input parameters and file paths."""
    # Check if input file exists
    input_path = Path(input_bam_path)
    if not input_path.exists():
        logging.error(f"Input BAM file '{input_bam_path}' not found")
        sys.exit(1)
    
    if not input_path.is_file():
        logging.error(f"Input path '{input_bam_path}' is not a file")
        sys.exit(1)
    
    # Check FASTA file if provided
    if fasta_file:
        fasta_path = Path(fasta_file)
        if not fasta_path.exists():
            logging.error(f"FASTA file '{fasta_file}' not found")
            sys.exit(1)
        if not fasta_path.is_file():
            logging.error(f"FASTA path '{fasta_file}' is not a file")
            sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Validate genomic methylation range
    if len(genomic_methylation_range) != 2:
        logging.error("Genomic methylation range must have exactly 2 values")
        sys.exit(1)
    
    if genomic_methylation_range[0] > genomic_methylation_range[1]:
        logging.error("Minimum genomic methylation percentage must be <= maximum genomic methylation percentage")
        sys.exit(1)
    
    if genomic_methylation_range[0] < 0 or genomic_methylation_range[1] > 1:
        logging.error("Genomic methylation percentages must be between 0 and 1")
        sys.exit(1)
    
    # Validate plasmid methylation range
    if len(plasmid_methylation_range) != 2:
        logging.error("Plasmid methylation range must have exactly 2 values")
        sys.exit(1)
    
    if plasmid_methylation_range[0] > plasmid_methylation_range[1]:
        logging.error("Minimum plasmid methylation percentage must be <= maximum plasmid methylation percentage")
        sys.exit(1)
    
    if plasmid_methylation_range[0] < 0 or plasmid_methylation_range[1] > 1:
        logging.error("Plasmid methylation percentages must be between 0 and 1")
        sys.exit(1)
    
    logging.info("Input validation passed")
    logging.info(f"Genomic methylation range: {genomic_methylation_range[0]:.2f} - {genomic_methylation_range[1]:.2f}")
    logging.info(f"Plasmid methylation range: {plasmid_methylation_range[0]:.2f} - {plasmid_methylation_range[1]:.2f}")
    
    # Log chromosome information
    logging.info(f"Processing {len(reference_length_dict)} chromosomes:")
    for chrom, length in reference_length_dict.items():
        logging.info(f"  {chrom}: {length:,} bp")
    
    # Check for chromosome mismatches between BAM and reference
    reference_chromosomes = set(reference_length_dict.keys())
    bam_chromosomes, missing_from_reference = check_bam_chromosomes(input_bam_path, reference_chromosomes)
    
    return bam_chromosomes, missing_from_reference


def count_reads(input_bam_path: str) -> int:
    """Count the total number of reads in a BAM file using samtools."""
    try:
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


def calculate_methylation_percentage(fiber) -> float:
    """Calculate the percentage of methylated AT basepairs for a given fiber."""
    at_count = fiber.seq.count('A') + fiber.seq.count('T')
    
    if at_count == 0:
        logging.debug(f"Read {fiber.qname} has no A or T bases, setting methylation to 0%")
        return 0.0
    
    m6a_count = len(fiber.m6a.starts)
    percent_methylation = m6a_count / at_count
    
    return percent_methylation


def passes_methylation_filter(fiber, genomic_methylation_range: List[float], 
                             plasmid_methylation_range: List[float]) -> bool:
    """Check if a fiber passes the methylation filter based on chromosome type."""
    percent_methylation = calculate_methylation_percentage(fiber)
    
    # Determine chromosome type and apply appropriate filter
    chrom_type = classify_chromosome(fiber.chrom)
    
    if chrom_type == 'genomic':
        return genomic_methylation_range[0] <= percent_methylation <= genomic_methylation_range[1]
    else:  # plasmid
        return plasmid_methylation_range[0] <= percent_methylation <= plasmid_methylation_range[1]


def filter_by_methylation(input_bam_path: str, output_bam_path: str, 
                         genomic_methylation_range: List[float] = [0.1, 1.0],
                         plasmid_methylation_range: List[float] = [0.1, 1.0]) -> Tuple[int, int, Dict[str, int]]:
    """Filter reads based on m6A methylation percentage with separate ranges for genomic and plasmid chromosomes."""
    logging.info("Step 1: Starting methylation filtering")
    
    total_reads = count_reads(input_bam_path)
    if total_reads == 0:
        logging.error("No reads found in input BAM file")
        return 0, 0, {}

    reads_written = 0
    reads_processed = 0
    reads_failed_methylation = 0
    reads_failed_genomic = 0
    reads_failed_plasmid = 0
    reads_passed_genomic = 0
    reads_passed_plasmid = 0
    writer = None
    
    try:
        logging.info("Opening input BAM file for methylation filtering")
        fiberbam = pyft.Fiberbam(input_bam_path)
        writer = pyft.Fiberwriter(output_bam_path, input_bam_path)

        logging.info(f"Processing reads for methylation filtering:")
        logging.info(f"  Genomic range: {genomic_methylation_range[0]:.2f} - {genomic_methylation_range[1]:.2f}")
        logging.info(f"  Plasmid range: {plasmid_methylation_range[0]:.2f} - {plasmid_methylation_range[1]:.2f}")
        
        for fiber in tqdm(fiberbam, total=total_reads, desc="Methylation filtering"):
            reads_processed += 1
            
            if passes_methylation_filter(fiber, genomic_methylation_range, plasmid_methylation_range):
                writer.write(fiber)
                reads_written += 1
                
                # Track statistics by chromosome type
                chrom_type = classify_chromosome(fiber.chrom)
                if chrom_type == 'genomic':
                    reads_passed_genomic += 1
                else:
                    reads_passed_plasmid += 1
            else:
                reads_failed_methylation += 1
                
                # Track failed reads by chromosome type
                chrom_type = classify_chromosome(fiber.chrom)
                if chrom_type == 'genomic':
                    reads_failed_genomic += 1
                else:
                    reads_failed_plasmid += 1
                
                # Log details for first few failed reads
                if reads_failed_methylation <= 5:
                    methylation_pct = calculate_methylation_percentage(fiber)
                    expected_range = genomic_methylation_range if chrom_type == 'genomic' else plasmid_methylation_range
                    logging.debug(f"Read {fiber.qname} from {chrom_type} {fiber.chrom} failed methylation filter: {methylation_pct:.3f} (expected: {expected_range[0]:.2f}-{expected_range[1]:.2f})")
            
            if reads_processed % 100000 == 0:
                logging.debug(f"Processed {reads_processed:,} reads, written {reads_written:,}, failed methylation: {reads_failed_methylation:,}")

    except Exception as e:
        logging.error(f"Error during methylation filtering: {e}")
        raise
    finally:
        # Clean up writer - pyft.Fiberwriter closes automatically when deleted
        if writer is not None:
            del writer

    # Log detailed methylation statistics
    logging.info(f"Methylation filtering statistics:")
    logging.info(f"  Reads processed: {reads_processed:,}")
    logging.info(f"  Reads written: {reads_written:,}")
    logging.info(f"  Reads failed methylation: {reads_failed_methylation:,}")
    logging.info(f"    - Genomic reads passed: {reads_passed_genomic:,}")
    logging.info(f"    - Genomic reads failed: {reads_failed_genomic:,}")
    logging.info(f"    - Plasmid reads passed: {reads_passed_plasmid:,}")
    logging.info(f"    - Plasmid reads failed: {reads_failed_plasmid:,}")
    
    if reads_processed != total_reads:
        logging.warning(f"  Read count mismatch! Expected {total_reads:,}, processed {reads_processed:,}")

    # Index the methylation-filtered BAM file for subsequent processing
    if reads_written > 0:
        logging.info("Indexing methylation-filtered BAM file...")
        try:
            subprocess.run(["samtools", "index", output_bam_path], check=True)
            logging.info("Successfully indexed methylation-filtered BAM file")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to index methylation-filtered BAM file: {e}")
            raise

    if total_reads > 0:
        pass_percentage = (reads_written / total_reads) * 100
        genomic_pass_rate = (reads_passed_genomic / (reads_passed_genomic + reads_failed_genomic) * 100) if (reads_passed_genomic + reads_failed_genomic) > 0 else 0
        plasmid_pass_rate = (reads_passed_plasmid / (reads_passed_plasmid + reads_failed_plasmid) * 100) if (reads_passed_plasmid + reads_failed_plasmid) > 0 else 0
        
        logging.info(f"Methylation filtering complete:")
        logging.info(f"  Total reads processed: {total_reads:,}")
        logging.info(f"  Reads passing filter: {reads_written:,}")
        logging.info(f"  Overall pass rate: {pass_percentage:.2f}%")
        logging.info(f"  Genomic pass rate: {genomic_pass_rate:.2f}%")
        logging.info(f"  Plasmid pass rate: {plasmid_pass_rate:.2f}%")
        
        print(f'Methylation filtering results:')
        print(f'  Total reads = {total_reads:,}')
        print(f'  Reads passing filter = {reads_written:,}')
        print(f'  Overall percentage passing filter = {pass_percentage:.2f}%')
        print(f'  Genomic reads: {reads_passed_genomic:,} passed, {reads_failed_genomic:,} failed ({genomic_pass_rate:.2f}% pass rate)')
        print(f'  Plasmid reads: {reads_passed_plasmid:,} passed, {reads_failed_plasmid:,} failed ({plasmid_pass_rate:.2f}% pass rate)')
    
    # Return statistics for potential use by calling functions
    methylation_stats = {
        'total_processed': reads_processed,
        'total_passed': reads_written,
        'genomic_passed': reads_passed_genomic,
        'genomic_failed': reads_failed_genomic,
        'plasmid_passed': reads_passed_plasmid,
        'plasmid_failed': reads_failed_plasmid
    }
    
    return total_reads, reads_written, methylation_stats


def plot_read_counts(read_counts: Dict[str, int], sample_name: str, 
                    output_folder: str, plot_type: str = "chromosome") -> None:
    """Create and save read count plots."""
    plt.figure(figsize=(12, 8))
    
    if plot_type == "chromosome":
        sns.set_style("whitegrid")
        sns.barplot(x=list(read_counts.keys()), y=list(read_counts.values()), 
                   hue=list(read_counts.keys()), palette="viridis", 
                   edgecolor='black', legend=False)
        plt.xlabel('Chromosome', fontsize=16)
        plt.ylabel('Read Count', fontsize=16)
        plt.title('Number of Reads Written to Each BAM File', fontsize=20)
        plt.xticks(rotation=45, fontsize=12)
        filename = f'{sample_name}_plasmid_read_counts.png'
    else:  # nucleosome count
        plt.bar(list(read_counts.keys()), list(read_counts.values()), 
               color='blue', edgecolor='black')
        plt.xlabel('Nucleosome Count', fontsize=16)
        plt.ylabel('Number of Reads', fontsize=16)
        plt.title(f'Reads vs Nucleosome Count - {sample_name}', fontsize=20)
        filename = f'reads_vs_nucleosome_count-{sample_name}.png'
    
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, filename), dpi=300)
    plt.close()  # This is plt.close(), not writer.close() - this is correct for matplotlib


def sort_reads_by_plasmid(input_bam: str, output_folder: str, 
                         reference_length_dict: Dict[str, int], sample_name: str,
                         length_range: Tuple[int, int], exact_alignment: bool) -> Dict[str, str]:
    """Step 2: Sort reads by plasmid chromosome."""
    logging.info("Step 2: Starting plasmid sorting")
    
    if not os.path.exists(input_bam):
        raise FileNotFoundError(f"Input BAM file not found: {input_bam}")
    
    os.makedirs(output_folder, exist_ok=True)
    
    output_bams = {}
    output_writers = {}
    read_counts = {chrom: 0 for chrom in reference_length_dict.keys()}
    skipped_reads = 0
    
    # Debug: Log all chromosomes and their classifications
    logging.info("Chromosome classifications:")
    for chrom in reference_length_dict.keys():
        chrom_type = classify_chromosome(chrom)
        logging.info(f"  {chrom} -> {chrom_type}")

    try:
        # Create BAM writers for each chromosome with organized directory structure
        for chrom in reference_length_dict.keys():
            # Classify chromosome and create appropriate directory structure
            chrom_type = classify_chromosome(chrom)
            chrom_output_folder = os.path.join(output_folder, chrom_type, chrom)
            os.makedirs(chrom_output_folder, exist_ok=True)
            logging.debug(f"Created directory: {chrom_output_folder}")

            output_bam_path = os.path.join(chrom_output_folder, f'{chrom}_{sample_name}.bam')
            output_writer = pyft.Fiberwriter(output_bam_path, input_bam)

            output_writers[chrom] = output_writer
            output_bams[chrom] = output_bam_path

        # Process reads
        fiberbam = pyft.Fiberbam(input_bam)
        logging.info("Processing reads for plasmid sorting...")
        
        # Debug: Track which chromosomes we see in the BAM
        chromosomes_seen = set()
        
        # Debug: Track read processing for detailed analysis
        total_reads_processed = 0
        total_reads_written = 0
        skipped_by_chromosome_filter = 0
        skipped_by_writer_missing = 0
        skipped_by_plasmid_length = 0
        skipped_by_plasmid_alignment = 0
        
        for fiber in tqdm(fiberbam, desc="Plasmid sorting"):
            total_reads_processed += 1
            
            # Track chromosomes we encounter
            chromosomes_seen.add(fiber.chrom)
            
            # Skip chrUn and _random chromosomes
            if should_skip_chromosome(fiber.chrom):
                skipped_reads += 1
                skipped_by_chromosome_filter += 1
                continue
                
            if fiber.chrom in output_writers:
                chrom_type = classify_chromosome(fiber.chrom)
                
                # Apply different filtering based on chromosome type
                if chrom_type == 'genomic':
                    # For genomic chromosomes: NO FILTERS AT ALL - just write the read
                    output_writers[fiber.chrom].write(fiber)
                    read_counts[fiber.chrom] += 1
                    total_reads_written += 1
                    logging.debug(f"Wrote genomic read from {fiber.chrom} (length: {fiber.get_seq_length()}, start: {fiber.start}, end: {fiber.end})")
                else:
                    # For plasmid chromosomes: apply both length and exact alignment filters
                    reference_length = reference_length_dict[fiber.chrom]
                    
                    # Apply length filter for plasmids
                    if (reference_length - length_range[0] <= fiber.get_seq_length() <= 
                        reference_length + length_range[1]):
                        
                        # Apply exact alignment filter if enabled
                        if not exact_alignment or (fiber.start == 0 and fiber.end == reference_length):
                            output_writers[fiber.chrom].write(fiber)
                            read_counts[fiber.chrom] += 1
                            total_reads_written += 1
                            logging.debug(f"Wrote plasmid read from {fiber.chrom} (length: {fiber.get_seq_length()}, expected: {reference_length})")
                        else:
                            skipped_by_plasmid_alignment += 1
                            logging.debug(f"Filtered out plasmid read from {fiber.chrom} due to alignment (start: {fiber.start}, end: {fiber.end}, expected end: {reference_length})")
                    else:
                        skipped_by_plasmid_length += 1
                        logging.debug(f"Filtered out plasmid read from {fiber.chrom} due to length (length: {fiber.get_seq_length()}, expected: {reference_length} ± {length_range})")
            else:
                skipped_by_writer_missing += 1
                logging.debug(f"Read from {fiber.chrom} skipped - not in output writers (not in reference)")

        # Log detailed read processing statistics
        logging.info(f"Detailed read processing statistics:")
        logging.info(f"  Total reads processed: {total_reads_processed:,}")
        logging.info(f"  Total reads written: {total_reads_written:,}")
        logging.info(f"  Reads skipped by chromosome filter (chrUn/_random): {skipped_by_chromosome_filter:,}")
        logging.info(f"  Reads skipped - chromosome not in reference: {skipped_by_writer_missing:,}")
        logging.info(f"  Reads skipped - plasmid length filter: {skipped_by_plasmid_length:,}")
        logging.info(f"  Reads skipped - plasmid alignment filter: {skipped_by_plasmid_alignment:,}")
        
        total_accounted = (total_reads_written + skipped_by_chromosome_filter + 
                          skipped_by_writer_missing + skipped_by_plasmid_length + skipped_by_plasmid_alignment)
        missing_reads = total_reads_processed - total_accounted
        if missing_reads != 0:
            logging.warning(f"  UNACCOUNTED READS: {missing_reads:,} - this suggests internal filtering!")


        # Debug: Log chromosomes we saw vs. what we expected
        logging.info(f"Chromosomes encountered in BAM: {sorted(chromosomes_seen)}")
        expected_chroms = set(reference_length_dict.keys())
        missing_chroms = expected_chroms - chromosomes_seen
        extra_chroms = chromosomes_seen - expected_chroms
        
        if missing_chroms:
            logging.warning(f"Expected chromosomes not found in BAM: {sorted(missing_chroms)}")
        if extra_chroms:
            logging.info(f"Additional chromosomes in BAM (not in reference): {sorted(extra_chroms)}")

    except Exception as e:
        logging.error(f"Error during plasmid sorting: {e}")
        raise e
    
    finally:
        # Always clean up writers - pyft.Fiberwriter closes automatically when deleted
        for chrom, writer in list(output_writers.items()):
            try:
                del writer
            except:
                pass
        output_writers.clear()

    # Log skipped reads
    if skipped_reads > 0:
        logging.info(f"Skipped {skipped_reads:,} reads from chrUn and _random chromosomes")

    # Print read counts and remove chromosomes with 0 reads from output_bams
    chromosomes_with_reads = {}
    genomic_reads = 0
    plasmid_reads = 0
    
    for chrom in read_counts:
        chrom_type = classify_chromosome(chrom)
        logging.info(f'Reads written to {chrom_type}/{chrom}/{os.path.basename(output_bams[chrom])}: {read_counts[chrom]}')
        
        if read_counts[chrom] > 0:
            chromosomes_with_reads[chrom] = output_bams[chrom]
            if chrom_type == 'genomic':
                genomic_reads += read_counts[chrom]
            else:
                plasmid_reads += read_counts[chrom]
        else:
            # Remove empty BAM files and their directories if empty
            try:
                os.remove(output_bams[chrom])
                # Try to remove empty chromosome directory
                chrom_dir = os.path.dirname(output_bams[chrom])
                try:
                    os.rmdir(chrom_dir)
                except OSError:
                    pass  # Directory not empty, that's fine
                logging.debug(f"Removed empty BAM file: {output_bams[chrom]}")
            except:
                pass

    # Log summary of reads by type
    total_reads = genomic_reads + plasmid_reads
    if total_reads > 0:
        logging.info(f"Read distribution summary:")
        logging.info(f"  Genomic reads: {genomic_reads:,} ({genomic_reads/total_reads*100:.1f}%)")
        logging.info(f"  Plasmid reads: {plasmid_reads:,} ({plasmid_reads/total_reads*100:.1f}%)")
        logging.info(f"  Total reads: {total_reads:,}")
    else:
        logging.warning("No reads found for any chromosomes after filtering")

    # Index all plasmid-sorted BAM files that have reads
    if chromosomes_with_reads:
        logging.info("Indexing plasmid-sorted BAM files...")
        for chrom, bam_path in chromosomes_with_reads.items():
            try:
                subprocess.run(["samtools", "index", bam_path], check=True)
                logging.debug(f"Successfully indexed {bam_path}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to index {bam_path}: {e}")
                raise
    else:
        logging.warning("No reads found for any chromosomes after filtering")

    # Create figures folder and plot (only for chromosomes with reads)
    figures_folder = os.path.join(output_folder, 'figures')
    os.makedirs(figures_folder, exist_ok=True)
    
    # Filter read_counts to only include chromosomes with reads for plotting
    filtered_read_counts = {chrom: count for chrom, count in read_counts.items() if count > 0}
    if filtered_read_counts:
        plot_read_counts(filtered_read_counts, sample_name, figures_folder, "chromosome")

    logging.info("Plasmid sorting complete")
    return chromosomes_with_reads  # Return only chromosomes with reads


def sort_by_nucleosome_count(input_bam_paths: Dict[str, str], 
                           chromosomes: List[str], 
                           output_folder: str,
                           sample_name: str) -> List[str]:
    """Step 3: Sort reads by nucleosome count for each chromosome."""
    logging.info("Step 3: Starting nucleosome count sorting")
    indexed_bam_files = []
    
    # Filter to only process chromosomes that have BAM files (i.e., had reads)
    chromosomes_to_process = [chrom for chrom in chromosomes if chrom in input_bam_paths]
    
    if not chromosomes_to_process:
        logging.warning("No chromosomes with reads found for nucleosome sorting")
        return indexed_bam_files
    
    logging.info(f"Processing nucleosome sorting for {len(chromosomes_to_process)} chromosomes with reads")
    
    for chrom in chromosomes_to_process:
        logging.info(f"Processing nucleosome sorting for chromosome: {chrom}")
        
        # Maintain the same directory structure (genomic/plasmid)
        chrom_type = classify_chromosome(chrom)
        destination_path = os.path.join(output_folder, chrom_type, chrom, 'sorted_by_nuc_count')
        os.makedirs(destination_path, exist_ok=True)
        
        input_path = input_bam_paths.get(chrom)
        if input_path is None:
            logging.warning(f'Input BAM file not found for chromosome {chrom}. Skipping...')
            continue
            
        if not os.path.exists(input_path):
            logging.warning(f'Input BAM file {input_path} does not exist. Skipping...')
            continue
            
        # Check if BAM file has an index
        index_path = input_path + '.bai'
        if not os.path.exists(index_path):
            logging.warning(f'Index file not found for {input_path}. Skipping...')
            continue

        output_bams = {}
        read_counts = {}
        
        try:
            fiberbam = pyft.Fiberbam(input_path)

            for fiber in tqdm(fiberbam, desc=f"Nucleosome sorting {chrom}"):
                if fiber.chrom != chrom:
                    continue
                    
                nuc_count = len(fiber.nuc.starts)
                
                if nuc_count not in output_bams:
                    output_bam_path = os.path.join(
                        destination_path, f'{nuc_count}_nucleosomes.bam')
                    output_bams[nuc_count] = pyft.Fiberwriter(output_bam_path, input_path)
                    
                output_bams[nuc_count].write(fiber)
                read_counts[nuc_count] = read_counts.get(nuc_count, 0) + 1

        except Exception as e:
            logging.error(f'Error processing chromosome {chrom}: {str(e)}')
            continue
        
        finally:
            # Clean up writers - pyft.Fiberwriter closes automatically when deleted
            for nuc_count, writer in list(output_bams.items()):
                try:
                    del writer
                    bam_path = os.path.join(destination_path, f'{nuc_count}_nucleosomes.bam')
                    indexed_bam_files.append(bam_path)
                except:
                    pass

        # Create plot for nucleosome counts (only if we have data)
        if read_counts:
            figures_folder = os.path.join(output_folder, 'figures')
            os.makedirs(figures_folder, exist_ok=True)
            plot_read_counts(read_counts, chrom, figures_folder, "nucleosome")
            
            logging.info(f'Nucleosome sorting for {chrom_type} chromosome {chrom} complete. Found {len(read_counts)} different nucleosome counts.')
        else:
            logging.warning(f'No reads found for chromosome {chrom} during nucleosome sorting')

    return indexed_bam_files


def collect_chromosome_statistics(chromosomes_with_reads: Dict[str, str], 
                                reference_length_dict: Dict[str, int]) -> Dict[str, Dict[str, int]]:
    """Collect statistics for each chromosome."""
    chromosome_stats = {}
    
    for chrom in reference_length_dict.keys():
        stats = {
            'final_written': 0,
            'length_bp': reference_length_dict[chrom]
        }
        
        # Count reads if BAM file exists
        if chrom in chromosomes_with_reads:
            bam_path = chromosomes_with_reads[chrom]
            if os.path.exists(bam_path):
                try:
                    # Count reads using samtools
                    result = subprocess.run(
                        ["samtools", "view", "-c", bam_path],
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    stats['final_written'] = int(result.stdout.strip())
                except:
                    stats['final_written'] = 0
        
        chromosome_stats[chrom] = stats
    
    return chromosome_stats


def generate_filtering_report(
    input_bam: str,
    output_folder: str,
    sample_name: str,
    reference_length_dict: Dict[str, int],
    fasta_source: str,
    genomic_methylation_range: List[float],
    plasmid_methylation_range: List[float],
    length_range: Tuple[int, int],
    exact_alignment: bool,
    sort_by_nucleosomes: bool,
    skip_methylation: bool,
    methylation_stats: Optional[Dict[str, int]],
    chromosome_stats: Dict[str, Dict[str, int]],
    output_bams: Dict[str, str],
    total_input_reads: int,
    command_line: str
) -> str:
    """Generate a comprehensive filtering report."""
    
    # Create report directory
    report_dir = os.path.join(output_folder, 'report')
    os.makedirs(report_dir, exist_ok=True)
    
    report_path = os.path.join(report_dir, f"{sample_name}_filtering_report.txt")
    
    # Get timestamp
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Classify chromosomes
    genomic_chroms = []
    plasmid_chroms = []
    for chrom in reference_length_dict.keys():
        if classify_chromosome(chrom) == 'genomic':
            genomic_chroms.append(chrom)
        else:
            plasmid_chroms.append(chrom)
    
    # Calculate totals
    total_output_reads = sum(chromosome_stats.get(chrom, {}).get('final_written', 0) 
                           for chrom in reference_length_dict.keys())
    genomic_reads = sum(chromosome_stats.get(chrom, {}).get('final_written', 0) 
                       for chrom in genomic_chroms)
    plasmid_reads = sum(chromosome_stats.get(chrom, {}).get('final_written', 0) 
                       for chrom in plasmid_chroms)
    
    # Get script version
    script_version = "comprehensive_plasmid_processor_v2.0"
    
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("COMPREHENSIVE PLASMID BAM PROCESSOR - FILTERING REPORT\n")
        f.write("=" * 80 + "\n")
        
        # RUN INFORMATION
        f.write("RUN INFORMATION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Timestamp:           {timestamp}\n")
        f.write(f"Script Version:      {script_version}\n")
        f.write(f"Input BAM:           {os.path.basename(input_bam)}\n")
        f.write(f"Input Path:          {os.path.abspath(input_bam)}\n")
        f.write(f"Output Directory:    {os.path.abspath(output_folder)}\n")
        f.write(f"Reference Source:    {fasta_source}\n")
        f.write(f"Sample Name:         {sample_name}\n")
        f.write(f"Host System:         {platform.node()}\n")
        f.write(f"Python Version:      {platform.python_version()}\n")
        f.write("\n")
        
        # FILTERING PARAMETERS
        f.write("FILTERING PARAMETERS:\n")
        f.write("-" * 40 + "\n")
        
        if skip_methylation:
            f.write("Methylation Filter:  DISABLED\n")
        else:
            if genomic_methylation_range == plasmid_methylation_range:
                f.write(f"Methylation Filter:  ENABLED (range: {genomic_methylation_range[0]:.2f} - {genomic_methylation_range[1]:.2f})\n")
                f.write("                     Both genomic and plasmid chromosomes use same range\n")
            else:
                f.write("Methylation Filter:  ENABLED (separate ranges)\n")
                f.write(f"                     Genomic: {genomic_methylation_range[0]:.2f} - {genomic_methylation_range[1]:.2f}\n")
                f.write(f"                     Plasmid: {plasmid_methylation_range[0]:.2f} - {plasmid_methylation_range[1]:.2f}\n")
            f.write("                     Reads must have m6A methylation within specified range\n")
        
        f.write(f"Length Filter:       PLASMID ONLY (±{length_range[0]}-{length_range[1]} bp from reference)\n")
        f.write("                     Genomic chromosomes: NO length restrictions\n")
        
        alignment_status = "exact alignment required" if exact_alignment else "any alignment accepted"
        f.write(f"Alignment Filter:    PLASMID ONLY ({alignment_status})\n")
        f.write("                     Genomic chromosomes: NO alignment restrictions\n")
        
        f.write("Chromosome Types:    GENOMIC, PLASMID\n")
        f.write("Excluded:            chrUn_* and *_random chromosomes\n")
        
        nucleosome_status = "ENABLED" if sort_by_nucleosomes else "DISABLED"
        f.write(f"Nucleosome Sorting:  {nucleosome_status}\n")
        f.write("\n")
        
        # CHROMOSOME CLASSIFICATIONS
        f.write("CHROMOSOME CLASSIFICATIONS:\n")
        f.write("-" * 40 + "\n")
        
        if genomic_chroms:
            genomic_display = genomic_chroms[:10]
            if len(genomic_chroms) > 10:
                remaining = len(genomic_chroms) - 10
                genomic_str = ", ".join(genomic_display) + f" ... and {remaining} more"
            else:
                genomic_str = ", ".join(genomic_chroms)
            f.write(f"Genomic Chromosomes ({len(genomic_chroms)}):\n")
            f.write(f"  {genomic_str}\n")
        
        if plasmid_chroms:
            f.write(f"Plasmid Chromosomes ({len(plasmid_chroms)}):\n")
            f.write(f"  {', '.join(plasmid_chroms)}\n")
        f.write("\n")
        
        # PROCESSING RESULTS
        f.write("PROCESSING RESULTS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total Input Reads:      {total_input_reads:,}\n")
        f.write(f"Total Output BAM Files: {len(output_bams)}\n")
        f.write(f"Total Reads Written:    {total_output_reads:,}\n")
        
        if total_output_reads > 0:
            genomic_pct = (genomic_reads / total_output_reads) * 100
            plasmid_pct = (plasmid_reads / total_output_reads) * 100
            f.write(f"  Genomic Reads:        {genomic_reads:,} ({genomic_pct:.1f}%)\n")
            f.write(f"  Plasmid Reads:        {plasmid_reads:,} ({plasmid_pct:.1f}%)\n")
        
        if total_input_reads > 0:
            overall_pass_rate = (total_output_reads / total_input_reads) * 100
            f.write(f"Overall Pass Rate:      {overall_pass_rate:.2f}%\n")
        f.write("\n")
        
        # DETAILED FILTERING STATISTICS
        if not skip_methylation and methylation_stats:
            f.write("METHYLATION FILTERING STATISTICS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Reads Processed:        {methylation_stats.get('total_processed', 0):,}\n")
            f.write(f"Reads Passed Filter:    {methylation_stats.get('total_passed', 0):,}\n")
            f.write(f"Genomic Reads Passed:   {methylation_stats.get('genomic_passed', 0):,}\n")
            f.write(f"Genomic Reads Failed:   {methylation_stats.get('genomic_failed', 0):,}\n")
            f.write(f"Plasmid Reads Passed:   {methylation_stats.get('plasmid_passed', 0):,}\n")
            f.write(f"Plasmid Reads Failed:   {methylation_stats.get('plasmid_failed', 0):,}\n")
            
            # Calculate pass rates
            genomic_total = methylation_stats.get('genomic_passed', 0) + methylation_stats.get('genomic_failed', 0)
            plasmid_total = methylation_stats.get('plasmid_passed', 0) + methylation_stats.get('plasmid_failed', 0)
            
            if genomic_total > 0:
                genomic_pass_rate = (methylation_stats.get('genomic_passed', 0) / genomic_total) * 100
                f.write(f"Genomic Pass Rate:      {genomic_pass_rate:.2f}%\n")
            
            if plasmid_total > 0:
                plasmid_pass_rate = (methylation_stats.get('plasmid_passed', 0) / plasmid_total) * 100
                f.write(f"Plasmid Pass Rate:      {plasmid_pass_rate:.2f}%\n")
            f.write("\n")
        
        # PER-CHROMOSOME BREAKDOWN
        f.write("PER-CHROMOSOME BREAKDOWN:\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Chromosome':<20} {'Type':<8} {'Reads':<12} {'Length (bp)':<12} {'Status'}\n")
        f.write("-" * 75 + "\n")
        
        for chrom in sorted(reference_length_dict.keys()):
            chrom_type = classify_chromosome(chrom)
            chrom_length = reference_length_dict[chrom]
            reads_written = chromosome_stats.get(chrom, {}).get('final_written', 0)
            
            if reads_written > 0:
                status = "✓ Written"
            else:
                status = "- No reads"
            
            f.write(f"{chrom:<20} {chrom_type:<8} {reads_written:<12,} {chrom_length:<12,} {status}\n")
        f.write("\n")
        
        # OUTPUT BAM FILES
        if output_bams:
            f.write("OUTPUT BAM FILES:\n")
            f.write("-" * 40 + "\n")
            f.write(f"{'Chromosome':<15} {'Type':<8} {'Reads':<10} {'File Path'}\n")
            f.write("-" * 85 + "\n")
            
            for chrom, bam_path in sorted(output_bams.items()):
                chrom_type = classify_chromosome(chrom)
                reads = chromosome_stats.get(chrom, {}).get('final_written', 0)
                rel_path = os.path.relpath(bam_path, output_folder)
                f.write(f"{chrom:<15} {chrom_type:<8} {reads:<10,} {rel_path}\n")
            f.write("\n")
        
        # FILTERING RULES APPLIED
        f.write("FILTERING RULES APPLIED:\n")
        f.write("-" * 40 + "\n")
        f.write("Genomic Chromosomes (chr1-22, chrX, chrY, chrM):\n")
        
        if not skip_methylation:
            f.write(f"  ✓ Methylation: {genomic_methylation_range[0]:.2f}-{genomic_methylation_range[1]:.2f}\n")
        else:
            f.write("  - Methylation: SKIPPED\n")
        
        f.write("  ✓ Length: NO restrictions (all read lengths accepted)\n")
        f.write("  ✓ Alignment: NO restrictions (any alignment position)\n")
        f.write("\n")
        
        f.write("Plasmid Chromosomes (all others):\n")
        
        if not skip_methylation:
            f.write(f"  ✓ Methylation: {plasmid_methylation_range[0]:.2f}-{plasmid_methylation_range[1]:.2f}\n")
        else:
            f.write("  - Methylation: SKIPPED\n")
        
        f.write(f"  ✓ Length: Within ±{length_range[0]}-{length_range[1]} bp of reference length\n")
        
        if exact_alignment:
            f.write("  ✓ Alignment: Must span entire plasmid (start=0, end=reference_length)\n")
        else:
            f.write("  ✓ Alignment: Any alignment position accepted\n")
        f.write("\n")
        
        # QUICK REFERENCE
        f.write("QUICK REFERENCE:\n")
        f.write("-" * 40 + "\n")
        f.write("To reproduce this analysis:\n")
        f.write(f"  {command_line}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("End of Report\n")
        f.write("=" * 80 + "\n")
    
    return report_path


def process_bam_pipeline(input_bam: str, output_folder: str, 
                        reference_length_dict: Dict[str, int], sample_name: str,
                        genomic_methylation_range: List[float], plasmid_methylation_range: List[float], 
                        length_range: Tuple[int, int], exact_alignment: bool, sort_by_nucleosomes: bool,
                        skip_methylation: bool) -> Tuple[Dict[str, str], int, Optional[Dict[str, int]]]:
    """Main pipeline function that coordinates all processing steps."""
    
    methylation_stats = None
    total_input_reads = count_reads(input_bam)
    
    # Step 1: Methylation filtering (optional)
    if not skip_methylation:
        methylation_filtered_bam = os.path.join(output_folder, f"{sample_name}_methylation_filtered.bam")
        total_reads, filtered_reads, methylation_stats = filter_by_methylation(
            input_bam, methylation_filtered_bam, genomic_methylation_range, plasmid_methylation_range
        )
        if filtered_reads == 0:
            logging.error("No reads passed methylation filtering. Pipeline stopped.")
            return {}, total_input_reads, methylation_stats
        current_bam = methylation_filtered_bam
    else:
        logging.info("Skipping methylation filtering as requested")
        current_bam = input_bam
    
    # Step 2: Plasmid sorting
    output_bams = sort_reads_by_plasmid(
        current_bam, output_folder, reference_length_dict, 
        sample_name, length_range, exact_alignment
    )
    
    # Step 3: Nucleosome count sorting (optional)
    if sort_by_nucleosomes:
        chromosomes = list(reference_length_dict.keys())
        indexed_bam_files = sort_by_nucleosome_count(output_bams, chromosomes, output_folder, sample_name)
        
        # Index nucleosome-sorted BAM files
        logging.info("Indexing nucleosome-sorted BAM files...")
        for bam_path in indexed_bam_files:
            try:
                subprocess.run(["samtools", "index", bam_path], check=True)
            except subprocess.CalledProcessError as e:
                logging.warning(f"Failed to index {bam_path}: {e}")
    
    return output_bams, total_input_reads, methylation_stats


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive BAM processing pipeline: methylation filtering → plasmid sorting → nucleosome count sorting",
        epilog="""Example commands: 
        # Using BAM header reference with separate methylation ranges (default - most convenient):
        python comprehensive_plasmid_processor.py input.bam --genomic_methylation_range 0.1 0.8 --plasmid_methylation_range 0.3 1.0 --sort_by_nucleosomes
        
        # Using legacy single methylation range for both:
        python comprehensive_plasmid_processor.py input.bam --methylation_range 0.1 1.0 --sort_by_nucleosomes
        
        # Using specific FASTA file (override):
        python comprehensive_plasmid_processor.py input.bam --fasta references.fasta --genomic_methylation_range 0.2 0.9 --sort_by_nucleosomes
        
        # Using manual dictionary (legacy):
        python comprehensive_plasmid_processor.py input.bam --reference_dict '{"pUC19": 2682, "SV40_Full": 4723}' --plasmid_methylation_range 0.4 1.0
        
        # Full pipeline with all options:
        python comprehensive_plasmid_processor.py input.bam --genomic_methylation_range 0.1 0.8 --plasmid_methylation_range 0.3 1.0 --length_range 100 150 --sort_by_nucleosomes --verbose""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('input_bam', help='Path to the input BAM file')
    
    # Reference input (mutually exclusive, with default behavior)
    ref_group = parser.add_mutually_exclusive_group(required=False)
    ref_group.add_argument('--fasta', '--reference_fasta', dest='fasta_file',
                          help='Path to FASTA file containing reference chromosomes (overrides BAM header detection)')
    ref_group.add_argument('--reference_dict', type=parse_reference_length_dict,
                          help='Dictionary of reference lengths for chromosomes in Python dict format (legacy)')
    
    # Optional arguments
    parser.add_argument('--output_folder', 
                        default=os.path.join(os.getcwd(), 'plasmid_bams_processed'),
                        help='Folder to save output files (default: plasmid_bams_processed)')
    parser.add_argument('--genomic_methylation_range', nargs=2, type=float, default=[0.1, 1.0],
                        metavar=('MIN', 'MAX'),
                        help='Methylation percentage range for GENOMIC chromosomes (default: 0.1 to 1.0)')
    parser.add_argument('--plasmid_methylation_range', nargs=2, type=float, default=[0.1, 1.0],
                        metavar=('MIN', 'MAX'),
                        help='Methylation percentage range for PLASMID chromosomes (default: 0.1 to 1.0)')
    parser.add_argument('--methylation_range', nargs=2, type=float, default=None,
                        metavar=('MIN', 'MAX'),
                        help='Methylation percentage range for BOTH genomic and plasmid (legacy - overrides separate ranges if specified)')
    parser.add_argument('--length_range', type=int, nargs=2, default=[50, 50], 
                        metavar=('LOWER_RANGE', 'UPPER_RANGE'),
                        help='Range of acceptable fiber lengths for PLASMID chromosomes only (default: [50, 50]). Genomic chromosomes accept all read lengths.')
    parser.add_argument('--exact_alignment', type=int, choices=[0, 1], default=1,
                        help='Apply exact alignment filter to PLASMID chromosomes only - read must span entire plasmid (0 for False, 1 for True, default: 1). Genomic chromosomes have no alignment restrictions.')
    parser.add_argument('--sort_by_nucleosomes', action='store_true',
                        help='Sort reads by nucleosome count after plasmid sorting')
    parser.add_argument('--skip_methylation', action='store_true',
                        help='Skip methylation filtering step')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging output')

    args = parser.parse_args()
    
    # Extract sample name from input BAM filename
    sample_name = extract_sample_name_from_bam(args.input_bam)
    
    # Set up logging
    setup_logging(args.verbose)
    
    # Handle methylation range arguments - support both legacy and new separate ranges
    if args.methylation_range is not None:
        # Legacy mode - use the same range for both genomic and plasmid
        genomic_methylation_range = args.methylation_range
        plasmid_methylation_range = args.methylation_range
        logging.info(f"Using legacy methylation range for both genomic and plasmid: {args.methylation_range}")
    else:
        # New mode - use separate ranges
        genomic_methylation_range = args.genomic_methylation_range
        plasmid_methylation_range = args.plasmid_methylation_range
        logging.info(f"Using separate methylation ranges - Genomic: {genomic_methylation_range}, Plasmid: {plasmid_methylation_range}")
    
    logging.info("Starting comprehensive BAM processing pipeline")
    logging.info(f"Input: {args.input_bam}")
    logging.info(f"Sample name: {sample_name}")
    logging.info(f"Output folder: {args.output_folder}")

    try:
        # Check dependencies
        check_dependencies()
        
        # Get reference length dictionary - priority order:
        # 1. Manual dictionary (--reference_dict)
        # 2. Specified FASTA file (--fasta)  
        # 3. Extract from BAM header (default)
        reference_length_dict = None
        fasta_source = None
        
        if args.reference_dict:
            logging.info("Using manually provided reference dictionary")
            reference_length_dict = args.reference_dict
            fasta_source = "manual_dict"
            
        elif args.fasta_file:
            logging.info(f"Using specified FASTA file: {args.fasta_file}")
            reference_length_dict = parse_fasta_file(args.fasta_file)
            fasta_source = args.fasta_file
            
        else:
            # Default behavior: extract from BAM header
            logging.info("Attempting to extract reference FASTA from BAM header...")
            extracted_fasta_path = extract_reference_from_bam_header(args.input_bam)
            
            if extracted_fasta_path:
                logging.info(f"Successfully extracted reference path from BAM header: {extracted_fasta_path}")
                
                # Check if the extracted path exists
                if os.path.exists(extracted_fasta_path):
                    logging.info("Reference file exists, parsing...")
                    reference_length_dict = parse_fasta_file(extracted_fasta_path)
                    fasta_source = extracted_fasta_path
                else:
                    logging.error(f"Reference file extracted from BAM header does not exist: {extracted_fasta_path}")
                    logging.error("Please provide a valid reference using --fasta or --reference_dict")
                    sys.exit(1)
            else:
                logging.error("Could not extract reference FASTA path from BAM header")
                logging.error("Please provide a reference using --fasta or --reference_dict")
                logging.error("Example: python script.py input.bam --fasta references.fasta")
                sys.exit(1)
        
        if not reference_length_dict:
            logging.error("No reference sequences found")
            sys.exit(1)
            
        logging.info(f"Reference source: {fasta_source}")
        
        # Validate inputs
        bam_chromosomes, missing_from_reference = validate_inputs(
            args.input_bam, args.output_folder, 
            genomic_methylation_range, plasmid_methylation_range, reference_length_dict, fasta_source if fasta_source != "manual_dict" else None)
        
        # Get the full command for the report
        command_line = ' '.join(sys.argv)
        
        # Convert exact_alignment from int to bool
        exact_alignment = bool(args.exact_alignment)

        # Run the complete pipeline
        output_bams, total_input_reads, methylation_stats = process_bam_pipeline(
            args.input_bam,
            args.output_folder,
            reference_length_dict,
            sample_name,
            genomic_methylation_range,
            plasmid_methylation_range,
            args.length_range,
            exact_alignment,
            args.sort_by_nucleosomes,
            args.skip_methylation
        )

        # Generate comprehensive statistics
        chromosome_stats = collect_chromosome_statistics(output_bams, reference_length_dict)

        # Generate the comprehensive report
        report_path = generate_filtering_report(
            input_bam=args.input_bam,
            output_folder=args.output_folder,
            sample_name=sample_name,
            reference_length_dict=reference_length_dict,
            fasta_source=fasta_source,
            genomic_methylation_range=genomic_methylation_range,
            plasmid_methylation_range=plasmid_methylation_range,
            length_range=args.length_range,
            exact_alignment=exact_alignment,
            sort_by_nucleosomes=args.sort_by_nucleosomes,
            skip_methylation=args.skip_methylation,
            methylation_stats=methylation_stats,
            chromosome_stats=chromosome_stats,
            output_bams=output_bams,
            total_input_reads=total_input_reads,
            command_line=command_line
        )

        # Index the final BAM files
        if output_bams:
            logging.info("Indexing final BAM files...")
            for bam_path in output_bams.values():
                try:
                    subprocess.run(["samtools", "index", bam_path], check=True)
                except subprocess.CalledProcessError as e:
                    logging.warning(f"Failed to index {bam_path}: {e}")

        print(f"\nPipeline complete! Output saved to: {args.output_folder}")
        print(f"Sample name used: {sample_name}")
        print(f"Comprehensive report: {report_path}")
        if fasta_source != "manual_dict":
            print(f"Reference used: {fasta_source}")
        logging.info("All processing steps completed successfully!")
        
    except KeyboardInterrupt:
        logging.info("Process interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        if args.verbose:
            logging.exception("Full traceback:")
        sys.exit(1)
    
    return 0


if __name__ == '__main__':
    exit(main())