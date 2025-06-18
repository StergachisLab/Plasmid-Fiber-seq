#!/usr/bin/env python3

import argparse
import pysam
import json
import sys
import logging
import os
import re
from collections import defaultdict, Counter
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import tempfile

# Set up logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def is_transition(ref, alt):
    """
    Determine if a SNV is a transition (A<->G or C<->T).

    Parameters:
    - ref: reference base
    - alt: alternate base

    Returns:
    - True if transition, False if transversion
    """
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    return (ref.upper(), alt.upper()) in transitions


def process_snv_bam(snv_data, read_to_variants, bam_template, reference_chrom, start_0based, end_0based, 
                   snv_bam_dir, target_variant_position=None, bam_log=None, quiet=False):
    """
    Process a single SNV and create its BAM file - designed for parallel processing.
    
    Parameters:
    - snv_data: tuple of (snv_id, read_names)
    - read_to_variants: dictionary mapping read names to their variants
    - bam_template: template BAM file to use for creating new BAM files
    - reference_chrom, start_0based, end_0based: region coordinates
    - snv_bam_dir: directory to save the BAM file
    - target_variant_position: position of the target variant (optional)
    - bam_log: file object for logging BAM details
    - quiet: whether to suppress informational messages
    
    Returns:
    - dict with analysis results and SNV information
    """
    snv_id, read_names = snv_data
    result = {
        'snv_id': snv_id,
        'success': False,
        'error': None,
        'analysis': {},
    }
    
    try:
        # Create a safe filename from the SNV ID
        safe_snv_id = re.sub(r'[^\w.-]', '_', snv_id)
        output_snv_bam = os.path.join(snv_bam_dir, f"{safe_snv_id}.bam")
        
        # Count reads with only this specific variant
        single_variant_count = 0
        single_variant_read_names = set()
        
        # Track variant combinations
        variant_combinations = Counter()
        
        # Track positions of all variants in these reads
        variant_positions = Counter()
        
        # Create a set of read names for faster lookup
        read_name_set = set(read_names)
        total_reads_for_snv = len(read_name_set)
        
        # Parse position from SNV_ID if not provided
        if target_variant_position is None:
            target_variant_position = int(snv_id.split(':')[0])
            
        variant_positions[target_variant_position] = total_reads_for_snv  # 100% by definition
        
        # Analyze variants in reads
        for read_name in read_name_set:
            variants = read_to_variants.get(read_name, [])
            
            # If only one variant and it's our target SNV
            if len(variants) == 1 and variants[0] == snv_id:
                single_variant_count += 1
                single_variant_read_names.add(read_name)
            
            # Add this combination to our counter
            variant_combinations[frozenset(variants)] += 1
            
            # Extract positions from all variants
            for var in variants:
                if var != snv_id:  # Skip the target SNV since we already counted it
                    # Extract position from variant ID
                    if ">" in var:  # SNV format
                        pos = int(var.split(':')[0])
                        variant_positions[pos] += 1
                    elif "+" in var:  # Insertion format
                        pos = int(var.split(':')[0])
                        variant_positions[pos] += 1
                    else:  # Deletion format
                        parts = var.split(':')
                        pos = int(parts[0])
                        variant_positions[pos] += 1
        
        # Find most common variant combination
        most_common_combo, most_common_count = variant_combinations.most_common(1)[0]
        
        # Calculate number of unique combinations
        unique_combo_count = len(variant_combinations)
        
        # Calculate position frequencies as percentages
        position_percentages = {pos: (count / total_reads_for_snv * 100) 
                               for pos, count in variant_positions.items()}
        
        # Convert most common combination to string for tag
        most_common_combo_str = json.dumps(sorted(list(most_common_combo)))
        
        # Convert position percentages to string for tag
        position_percentages_str = json.dumps(position_percentages)
        
        # Calculate the single variant proportion
        single_variant_proportion = single_variant_count / total_reads_for_snv if total_reads_for_snv > 0 else 0
        
        # Store analysis results
        result['analysis'] = {
            'total_reads': total_reads_for_snv,
            'single_variant_count': single_variant_count,
            'single_variant_proportion': single_variant_proportion,
            'unique_combo_count': unique_combo_count,
            'most_common_combo': sorted(list(most_common_combo)),
            'most_common_count': most_common_count,
            'position_percentages': position_percentages,
            'single_variant_read_names': single_variant_read_names,
            'bam_path': output_snv_bam
        }
        
        # Create a temporary BAM file for faster writing
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as temp_file:
            temp_bam_path = temp_file.name
        
        # Open BAM files
        with pysam.AlignmentFile(temp_bam_path, "wb", template=bam_template) as snv_bam:
            # Write all reads with this SNV to the temp BAM
            for read in read_to_variants.keys():
                if read in read_name_set:
                    # Get the variants for this read
                    variants = read_to_variants.get(read, [])
                    
                    # Convert to the same format as used for the VR tag
                    var = variants if variants else "WT"
                    var_serialized = json.dumps(var)
                    
                    # Create a new AlignmentSegment
                    new_read = pysam.AlignedSegment()
                    new_read.query_name = read
                    
                    # Add the tags
                    new_read.set_tag("VR", var_serialized)
                    new_read.set_tag("SV", single_variant_proportion)
                    new_read.set_tag("UC", unique_combo_count)
                    new_read.set_tag("MC", most_common_combo_str)
                    new_read.set_tag("VP", position_percentages_str)
                    
                    # Write the read with all tags to the SNV-specific BAM
                    snv_bam.write(new_read)
        
        # Move the temp BAM to the final location
        os.rename(temp_bam_path, output_snv_bam)
        
        # Log information about this BAM file
        if bam_log:
            variant_positions_str = ",".join([f"{pos}:{percent:.1f}%" for pos, percent in position_percentages.items()])
            most_common_str = ",".join(sorted(list(most_common_combo)))
            log_entry = f"SNV\t{snv_id}\t{output_snv_bam}\t{total_reads_for_snv}\t{single_variant_count}\t{single_variant_proportion:.4f}\t{unique_combo_count}\t{most_common_str}\t{most_common_count}\t{most_common_count/total_reads_for_snv:.4f}\t{variant_positions_str}\n"
            with open(bam_log, 'a') as f:
                f.write(log_entry)
        
        result['success'] = True
        
    except Exception as e:
        result['error'] = str(e)
        
    return result


def process_single_variant_bam(snv_data, bam_template, single_variant_bam_dir, bam_log=None, quiet=False):
    """
    Process a single SNV and create a BAM file with only reads that have this single variant.
    
    Parameters:
    - snv_data: tuple of (snv_id, read_names)
    - bam_template: template BAM file to use for creating new BAM files
    - single_variant_bam_dir: directory to save the BAM file
    - bam_log: file object for logging BAM details
    - quiet: whether to suppress informational messages
    
    Returns:
    - dict with results
    """
    snv_id, single_variant_read_names = snv_data
    result = {
        'snv_id': snv_id,
        'success': False,
        'error': None,
        'read_count': 0,
    }
    
    try:
        # Create a safe filename
        safe_snv_id = re.sub(r'[^\w.-]', '_', snv_id)
        output_single_variant_bam = os.path.join(single_variant_bam_dir, f"{safe_snv_id}_single.bam")
        
        # Create a temporary BAM file for faster writing
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as temp_file:
            temp_bam_path = temp_file.name
        
        # Write only reads that have this single variant
        with pysam.AlignmentFile(temp_bam_path, "wb", template=bam_template) as single_variant_bam:
            read_count = 0
            
            # Get target variant position
            target_variant_position = int(snv_id.split(':')[0])
            
            for read_name in single_variant_read_names:
                # Create a new AlignmentSegment
                new_read = pysam.AlignedSegment()
                new_read.query_name = read_name
                
                # Add the VR tag
                var_serialized = json.dumps([snv_id])
                new_read.set_tag("VR", var_serialized)
                
                # Write to the single-variant BAM
                single_variant_bam.write(new_read)
                read_count += 1
        
        # Move the temp BAM to the final location
        os.rename(temp_bam_path, output_single_variant_bam)
        
        # Log information about this single-variant BAM file
        if bam_log:
            variant_positions_str = f"{target_variant_position}:100.0%"  # By definition, only this variant
            log_entry = f"SingleVariant\t{snv_id}\t{output_single_variant_bam}\t{read_count}\t{read_count}\t1.0000\t1\t{snv_id}\t{read_count}\t1.0000\t{variant_positions_str}\n"
            with open(bam_log, 'a') as f:
                f.write(log_entry)
        
        result['success'] = True
        result['read_count'] = read_count
        
    except Exception as e:
        result['error'] = str(e)
        
    return result


def add_var_tag(bam_file, output_bam_file, reference_fasta, reference_chrom, start, end,
             snv=False, insertion=False, deletion=False,
             snv_output=None, vr_output=None, report_output=None,
             write_snv_bams=False, snv_bam_dir="sorted_by_SNV",
             min_coverage_threshold=500, min_snv_bam_coverage=1,
             single_variant_cutoff=None, single_variant_bam_dir=None,
             quiet=False, bam_log_file=None, threads=1, chunk_size=1000):
    """
    Extract SNVs, insertions, and deletions from reads in a BAM file within a reference range by parsing the CIGAR string.
    Add the variants as a new BAM tag for each read. Track unique SNVs and VR tags with their counts.

    Parameters:
    - bam_file: str, path to the input BAM file.
    - output_bam_file: str, path to the output BAM file.
    - reference_fasta: str, path to the reference FASTA file.
    - reference_chrom: str, reference chromosome to fetch.
    - start: int, start coordinate of the reference range (1-based, inclusive).
    - end: int, end coordinate of the reference range (1-based, inclusive).
    - snv: bool, whether to include SNVs in the output tag.
    - insertion: bool, whether to include insertions in the output tag.
    - deletion: bool, whether to include deletions in the output tag.
    - snv_output: str, optional path to output file for unique SNVs and their counts.
    - vr_output: str, optional path to output file for unique VR tags and their counts.
    - min_coverage_threshold: int, threshold for reporting high-coverage variants (default: 500).
    - min_snv_bam_coverage: int, minimum read coverage required to write SNV specific bam.
    - threads: int, number of parallel processes to use (default: 1).
    - chunk_size: int, number of reads to process at once (default: 1000).
    """
    try:
        # Set logging level based on quiet flag
        if quiet:
            logger.setLevel(logging.WARNING)
        else:
            logger.setLevel(logging.INFO)
            
        # Set up BAM log file if specified
        if bam_log_file:
            with open(bam_log_file, 'w') as bam_log:
                bam_log.write("BAM_Type\tSNV_ID\tBAM_Path\tTotal_Reads\tSingle_Variant_Reads\tSingle_Variant_Proportion\tUnique_Combinations\tMost_Common_Combo\tMost_Common_Count\tMost_Common_Proportion\tVariant_Positions\n")
            
        # Open the BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
        output_bam = pysam.AlignmentFile(output_bam_file, "wb", template=bam)

        # Open the reference FASTA file
        fasta = pysam.FastaFile(reference_fasta)

        # Get chromosome length to validate coordinates
        chrom_lengths = dict(zip(fasta.references, fasta.lengths))
        if reference_chrom not in chrom_lengths:
            raise ValueError(f"Chromosome {reference_chrom} not found in reference FASTA")

        chrom_length = chrom_lengths[reference_chrom]

        # Validate coordinates
        if start <= 0 or end <= 0:
            raise ValueError("Start and end positions must be positive (1-based)")
        if start > chrom_length or end > chrom_length:
            raise ValueError(f"Coordinates exceed chromosome length ({chrom_length})")
        if start > end:
            raise ValueError("Start position must be less than or equal to end position")

        # Convert to 0-based coordinates for internal processing
        start_0based = start - 1
        end_0based = end - 1

        # Fetch only the required region of the reference
        ref_region = fasta.fetch(reference_chrom, start_0based, end_0based + 1)

        # Count processed reads and variants
        processed_reads = 0
        total_snvs = 0
        total_insertions = 0
        total_deletions = 0
        wt_count = 0  # Count of wild-type reads (no variants)

        # Dictionaries to track unique SNVs and VR tags with their counts
        unique_snvs = {}  # Format: {snv_id: count}
        unique_vr_tags = {}  # Format: {vr_tag_json: count}

        # For tracking transitions and transversions
        transitions = 0
        transversions = 0
        transition_snvs = set()
        transversion_snvs = set()

        # For tracking SNV-specific reads if write_snv_bams is True
        snv_to_reads = defaultdict(list) if write_snv_bams else None
        
        # New: For tracking read variants for each SNV-specific BAM
        if write_snv_bams:
            # Store full variant information for each read
            read_to_variants = {}  # Format: {read_name: list_of_variants}

        # Create SNV BAM directory if needed
        if write_snv_bams:
            os.makedirs(snv_bam_dir, exist_ok=True)
            if not quiet:
                logger.info(f"Created directory for SNV-specific BAM files: {snv_bam_dir}")
            
        # Create single-variant BAM directory if needed
        if single_variant_cutoff is not None and single_variant_bam_dir is not None:
            os.makedirs(single_variant_bam_dir, exist_ok=True)
            if not quiet:
                logger.info(f"Created directory for single-variant BAM files: {single_variant_bam_dir}")
                
        # Get total number of reads for progress bar
        total_reads = 0
        if not quiet:
            try:
                # Count reads in region for progress bar
                for _ in bam.fetch(reference_chrom, start_0based, end_0based + 1):
                    total_reads += 1
                bam.reset()
            except:
                # If we can't count, just use a placeholder
                total_reads = None

        # First BAM pass: Process each read and collect SNV information
        reads_iterator = bam.fetch(reference_chrom, start_0based, end_0based + 1)
        if not quiet:
            reads_iterator = tqdm(reads_iterator, total=total_reads, desc="Processing reads", unit="reads")
            
        for read in reads_iterator:
            # Skip unmapped, secondary, or supplementary alignments
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Extract relevant fields
            cigar = read.cigar
            query_seq = read.query_sequence
            ref_pos = read.reference_start  # 0-based

            # Skip if read does not overlap with the region of interest
            read_end = read.reference_end  # 0-based, exclusive
            if read_end is None or ref_pos > end_0based or read_end <= start_0based:
                continue

            # Containers for variants
            snvs = []
            insertions = []
            deletions = []

            # Parse CIGAR string for variants
            query_index = 0  # 0-based index for the query sequence
            current_ref_pos = ref_pos  # 0-based

            for op, length in cigar:
                if current_ref_pos > end_0based:  # Past our region of interest
                    break

                # PacBio CIGAR operations
                if op == 1:  # Insertion (I)
                    if start_0based <= current_ref_pos <= end_0based:
                        inserted_seq = query_seq[query_index:query_index + length]
                        # Record as 1-based for output
                        insertions.append(f"{current_ref_pos + 1}:+{inserted_seq}")
                    query_index += length

                elif op == 2:  # Deletion (D)
                    if current_ref_pos + length > start_0based and current_ref_pos <= end_0based:
                        # Calculate overlap with our region of interest
                        overlap_start = max(current_ref_pos, start_0based)
                        overlap_end = min(current_ref_pos + length, end_0based + 1)

                        if overlap_start < overlap_end:
                            # Get the deleted reference sequence
                            rel_start = overlap_start - start_0based
                            rel_end = overlap_end - start_0based
                            deleted_bases = ref_region[rel_start:rel_end]

                            # Record as 1-based for output
                            deletions.append(f"{overlap_start + 1}:{overlap_end - overlap_start}{deleted_bases}")

                    current_ref_pos += length

                elif op == 4:  # Soft clipping (S)
                    query_index += length

                elif op == 5:  # Hard clipping (H)
                    # No sequence in read, just skip
                    pass

                elif op == 7:  # Sequence match (=)
                    # Move both pointers forward
                    if current_ref_pos + length > start_0based and current_ref_pos <= end_0based:
                        # Calculate overlap with our region of interest
                        overlap_start = max(current_ref_pos, start_0based)
                        overlap_end = min(current_ref_pos + length, end_0based + 1)

                        # Adjust query position for the overlap start
                        adjusted_query_index = query_index + (overlap_start - current_ref_pos)

                        # No SNVs to record for matched sequences
                        pass

                    query_index += length
                    current_ref_pos += length

                elif op == 8:  # Sequence mismatch (X)
                    if current_ref_pos + length > start_0based and current_ref_pos <= end_0based:
                        # Calculate overlap with our region of interest
                        overlap_start = max(current_ref_pos, start_0based)
                        overlap_end = min(current_ref_pos + length, end_0based + 1)

                        # Adjust query position for the overlap start
                        adjusted_query_index = query_index + (overlap_start - current_ref_pos)

                        for i in range(overlap_end - overlap_start):
                            rel_pos = overlap_start + i - start_0based
                            if rel_pos < len(ref_region):
                                ref_base = ref_region[rel_pos]
                                query_base = query_seq[adjusted_query_index + i]

                                # Record SNV as 1-based for output
                                snv_id = f"{overlap_start + i + 1}:{ref_base}>{query_base}"
                                snvs.append(snv_id)

                                # Determine if this is a transition or transversion
                                if is_transition(ref_base, query_base):
                                    transition_snvs.add(snv_id)
                                else:
                                    transversion_snvs.add(snv_id)

                    query_index += length
                    current_ref_pos += length

                elif op == 0:  # ALIGNMENT_MATCH (M) - forbidden in PacBio BAM spec
                    logger.warning(f"Read {read.query_name} contains forbidden ALIGNMENT_MATCH (M) CIGAR operation. "
                                  f"This violates the PacBio BAM spec. Skipping this read.")
                    break

                else:  # Other operations (N, P) - skip for now
                    if op == 3:  # REFERENCE_SKIP (N)
                        current_ref_pos += length
                    # For PADDING (P), we don't move either pointer

            # If we encountered a forbidden operation, skip this read
            if op == 0:
                continue

            # Collect all variants in a single flat list
            all_variants = []
            if snv and snvs:
                all_variants.extend(snvs)
                total_snvs += len(snvs)
            if insertion and insertions:
                all_variants.extend(insertions)
                total_insertions += len(insertions)
            if deletion and deletions:
                all_variants.extend(deletions)
                total_deletions += len(deletions)

            # Add variants as a new BAM tag
            if not all_variants:  # If no variants are included
                var = "WT"
                wt_count += 1
            else:
                var = all_variants

            # Convert list to JSON string
            var_serialized = json.dumps(var)

            # Track unique VR tags
            if var_serialized in unique_vr_tags:
                unique_vr_tags[var_serialized] += 1
            else:
                unique_vr_tags[var_serialized] = 1

            # Track individual SNVs
            if snv and snvs:
                for snv_id in snvs:
                    if snv_id in unique_snvs:
                        unique_snvs[snv_id] += 1
                    else:
                        unique_snvs[snv_id] = 1

                    # If write_snv_bams is True, store read object for each SNV
                    if write_snv_bams:
                        # Store the read name instead of the read object to avoid modifying original read
                        snv_to_reads[snv_id].append(read.query_name)
                        
                        # New: Store all variants for this read
                        read_to_variants[read.query_name] = all_variants

            read.set_tag("VR", var_serialized)  # Add the 'VR' tag with variant information
            output_bam.write(read)
            processed_reads += 1

        # Log summary statistics only if not quiet
        if not quiet:
            logger.info(f"Processed {processed_reads} reads in region {reference_chrom}:{start}-{end}")
            logger.info(f"Found {total_snvs} SNVs, {total_insertions} insertions, and {total_deletions} deletions")
            logger.info(f"Wild-type reads (no variants): {wt_count} ({(wt_count/processed_reads*100):.2f}% of total)")

        # Filter SNVs with coverage below the threshold
        filtered_snvs = {}
        for snv_id, read_names in snv_to_reads.items():
            if len(read_names) >= min_snv_bam_coverage:
                filtered_snvs[snv_id] = read_names
        
        snv_bams_written = 0
        snv_bams_skipped = len(snv_to_reads) - len(filtered_snvs)
        
        # Process SNV-specific BAMs in parallel if requested
        if write_snv_bams and filtered_snvs:
            snv_items = list(filtered_snvs.items())
            if not quiet:
                logger.info(f"Creating {len(snv_items)} SNV-specific BAM files using {threads} threads...")
                
            # Process SNVs in parallel if multiple threads are available
            if threads > 1 and len(snv_items) > 1:
                # Create a pool of workers
                with mp.Pool(processes=threads) as pool:
                    # Create a partial function with fixed arguments
                    process_func = partial(
                        process_snv_bam,
                        read_to_variants=read_to_variants,
                        bam_template=bam,
                        reference_chrom=reference_chrom,
                        start_0based=start_0based,
                        end_0based=end_0based,
                        snv_bam_dir=snv_bam_dir,
                        bam_log=bam_log_file,
                        quiet=quiet
                    )
                    
                    # Process SNVs in chunks for better progress reporting
                    results = []
                    for i in range(0, len(snv_items), chunk_size):
                        chunk = snv_items[i:i+chunk_size]
                        
                        # Show progress if not quiet
                        if not quiet:
                            with tqdm(total=len(chunk), desc=f"Processing SNV BAMs {i+1}-{i+len(chunk)}/{len(snv_items)}", unit="BAMs") as pbar:
                                chunk_results = []
                                for result in pool.imap_unordered(process_func, chunk):
                                    pbar.update()
                                    chunk_results.append(result)
                                results.extend(chunk_results)
                        else:
                            # Process without progress bar
                            chunk_results = list(pool.imap_unordered(process_func, chunk))
                            results.extend(chunk_results)
                            
                    # Count successful BAMs
                    snv_bams_written = sum(1 for r in results if r['success'])
                    
                    # Get analysis data for single variant filtering
                    snvs_meeting_cutoff = []
                    reads_with_single_variant = {}
                    
                    for result in results:
                        if result['success']:
                            snv_id = result['snv_id']
                            analysis = result['analysis']
                            
                            # Check if this SNV meets the single-variant cutoff
                            meets_cutoff = False
                            if single_variant_cutoff is not None:
                                # If the cutoff is a float between 0-1, treat it as a proportion
                                if 0 <= single_variant_cutoff <= 1:
                                    meets_cutoff = analysis['single_variant_proportion'] >= single_variant_cutoff
                                # Otherwise treat it as an absolute count
                                else:
                                    meets_cutoff = analysis['single_variant_count'] >= single_variant_cutoff
                                    
                            if meets_cutoff and single_variant_bam_dir is not None:
                                snvs_meeting_cutoff.append(snv_id)
                                reads_with_single_variant[snv_id] = analysis['single_variant_read_names']
            else:
                # Process SNVs sequentially if only one thread or just one SNV
                results = []
                if not quiet:
                    snv_progress = tqdm(snv_items, desc="Creating SNV BAMs", unit="BAMs")
                else:
                    snv_progress = snv_items
                    
                # Store SNV data for the single variant BAM creation
                snvs_meeting_cutoff = []
                reads_with_single_variant = {}
                
                for snv_data in snv_progress:
                    snv_id = snv_data[0]
                    result = process_snv_bam(
                        snv_data=snv_data,
                        read_to_variants=read_to_variants,
                        bam_template=bam,
                        reference_chrom=reference_chrom,
                        start_0based=start_0based,
                        end_0based=end_0based,
                        snv_bam_dir=snv_bam_dir,
                        bam_log=bam_log_file,
                        quiet=quiet
                    )
                    
                    results.append(result)
                    
                    if result['success']:
                        snv_bams_written += 1
                        
                        # Check if this SNV meets the single-variant cutoff
                        if single_variant_cutoff is not None and single_variant_bam_dir is not None:
                            analysis = result['analysis']
                            meets_cutoff = False
                            
                            # If the cutoff is a float between 0-1, treat it as a proportion
                            if 0 <= single_variant_cutoff <= 1:
                                meets_cutoff = analysis['single_variant_proportion'] >= single_variant_cutoff
                            # Otherwise treat it as an absolute count
                            else:
                                meets_cutoff = analysis['single_variant_count'] >= single_variant_cutoff
                                
                            if meets_cutoff:
                                snvs_meeting_cutoff.append(snv_id)
                                reads_with_single_variant[snv_id] = analysis['single_variant_read_names']
                        
                        # Update progress bar description if not quiet
                        if not quiet:
                            snv_progress.set_description(f"Created {snv_bams_written} SNV BAMs")
            
            if not quiet:
                logger.info(f"Wrote {snv_bams_written} SNV-specific BAM files to {snv_bam_dir}/")
                if snv_bams_skipped > 0:
                    logger.info(f"Skipped {snv_bams_skipped} SNVs with coverage below {min_snv_bam_coverage}x")
                
            # Third pass: Create BAM files for SNVs that meet the single-variant cutoff
            if single_variant_cutoff is not None and single_variant_bam_dir is not None and snvs_meeting_cutoff:
                single_variant_bams_written = 0
                
                if not quiet:
                    logger.info(f"Creating BAM files for {len(snvs_meeting_cutoff)} SNVs that meet the single-variant cutoff: {single_variant_cutoff}")
                
                # Process single-variant BAMs in parallel if multiple threads are available
                if threads > 1 and len(snvs_meeting_cutoff) > 1:
                    # Create data for parallel processing
                    single_variant_items = [(snv_id, reads_with_single_variant[snv_id]) for snv_id in snvs_meeting_cutoff]
                    
                    # Create a pool of workers
                    with mp.Pool(processes=threads) as pool:
                        # Create a partial function with fixed arguments
                        process_func = partial(
                            process_single_variant_bam,
                            bam_template=bam,
                            single_variant_bam_dir=single_variant_bam_dir,
                            bam_log=bam_log_file,
                            quiet=quiet
                        )
                        
                        # Process SNVs in chunks for better progress reporting
                        results = []
                        for i in range(0, len(single_variant_items), chunk_size):
                            chunk = single_variant_items[i:i+chunk_size]
                            
                            # Show progress if not quiet
                            if not quiet:
                                with tqdm(total=len(chunk), desc=f"Processing single-variant BAMs {i+1}-{i+len(chunk)}/{len(single_variant_items)}", unit="BAMs") as pbar:
                                    chunk_results = []
                                    for result in pool.imap_unordered(process_func, chunk):
                                        pbar.update()
                                        chunk_results.append(result)
                                    results.extend(chunk_results)
                            else:
                                # Process without progress bar
                                chunk_results = list(pool.imap_unordered(process_func, chunk))
                                results.extend(chunk_results)
                                
                        # Count successful BAMs
                        single_variant_bams_written = sum(1 for r in results if r['success'])
                else:
                    # Process single-variant BAMs sequentially
                    if not quiet:
                        single_variant_progress = tqdm(snvs_meeting_cutoff, desc="Creating single-variant BAMs", unit="BAMs")
                    else:
                        single_variant_progress = snvs_meeting_cutoff
                    
                    for snv_id in single_variant_progress:
                        result = process_single_variant_bam(
                            (snv_id, reads_with_single_variant[snv_id]),
                            bam_template=bam,
                            single_variant_bam_dir=single_variant_bam_dir,
                            bam_log=bam_log_file,
                            quiet=quiet
                        )
                        
                        if result['success']:
                            single_variant_bams_written += 1
                            
                            # Update progress bar description if not quiet
                            if not quiet:
                                single_variant_progress.set_description(f"Created {single_variant_bams_written} single-variant BAMs")
                
                if not quiet:
                    logger.info(f"Wrote {single_variant_bams_written} single-variant BAM files to {single_variant_bam_dir}/")

        # Calculate transitions and transversions from the unique SNVs
        transitions = len(transition_snvs)
        transversions = len(transversion_snvs)
        total_unique_snvs = len(unique_snvs)

        # Verify the counts
        if transitions + transversions != total_unique_snvs:
            logger.warning(f"Transition ({transitions}) + Transversion ({transversions}) count doesn't match total unique SNVs ({total_unique_snvs})")

        # Calculate transition/transversion ratio
        ti_tv_ratio = transitions / transversions if transversions > 0 else float('inf')

        # Write unique SNVs and their counts to file if specified
        if snv_output and unique_snvs:
            with open(snv_output, 'w') as f:
                f.write("SNV_ID\tCount\tType\n")
                for snv_id, count in sorted(unique_snvs.items(), key=lambda x: x[1], reverse=True):
                    snv_type = "Transition" if snv_id in transition_snvs else "Transversion"
                    f.write(f"{snv_id}\t{count}\t{snv_type}\n")
            if not quiet:
                logger.info(f"Wrote {len(unique_snvs)} unique SNVs to {snv_output}")

            # Calculate statistics for SNVs with coverage > threshold
            high_coverage_snvs = sum(1 for count in unique_snvs.values() if count > min_coverage_threshold)
            percent_high_coverage = (high_coverage_snvs / len(unique_snvs)) * 100 if unique_snvs else 0
            if not quiet:
                logger.info(f"SNV Statistics: {len(unique_snvs)} unique SNVs, {high_coverage_snvs} ({percent_high_coverage:.2f}%) with coverage > {min_coverage_threshold}x")

        # Write unique VR tags and their counts to file if specified
        if vr_output and unique_vr_tags:
            with open(vr_output, 'w') as f:
                f.write("VR_Tag\tCount\n")
                for vr_tag, count in sorted(unique_vr_tags.items(), key=lambda x: x[1], reverse=True):
                    f.write(f"{vr_tag}\t{count}\n")
            if not quiet:
                logger.info(f"Wrote {len(unique_vr_tags)} unique VR tags to {vr_output}")

            # Calculate statistics for VR tags with coverage > threshold
            high_coverage_vr = sum(1 for count in unique_vr_tags.values() if count > min_coverage_threshold)
            percent_high_coverage_vr = (high_coverage_vr / len(unique_vr_tags)) * 100 if unique_vr_tags else 0
            if not quiet:
                logger.info(f"VR Tag Statistics: {len(unique_vr_tags)} unique VR tags, {high_coverage_vr} ({percent_high_coverage_vr:.2f}%) with coverage > {min_coverage_threshold}x")

        # Print summary report to stdout and a file if specified
        report_lines = []
        report_lines.append("===== SUMMARY REPORT =====")
        report_lines.append(f"Processed {processed_reads} reads in region {reference_chrom}:{start}-{end}")
        report_lines.append(f"Found {total_snvs} total SNV occurrences, {total_insertions} insertions, and {total_deletions} deletions")
        report_lines.append(f"Wild-type reads (no variants): {wt_count} ({(wt_count/processed_reads*100):.2f}% of total)")

        if write_snv_bams and snv_to_reads:
            if snv_bams_skipped > 0:
                report_lines.append(f"Wrote BAM files for {snv_bams_written} unique SNVs to {snv_bam_dir}/ (skipped {snv_bams_skipped} SNVs below {min_snv_bam_coverage}x coverage)")
            else:
                report_lines.append(f"Wrote BAM files for {snv_bams_written} unique SNVs to {snv_bam_dir}/")
                
        # Add single-variant BAM information to the report
        if single_variant_cutoff is not None and single_variant_bam_dir is not None:
            if 0 <= single_variant_cutoff <= 1:
                cutoff_desc = f"{single_variant_cutoff*100:.1f}% proportion"
            else:
                cutoff_desc = f"{single_variant_cutoff} count"
                
            if single_variant_bams_written > 0:
                report_lines.append(f"Wrote {single_variant_bams_written} single-variant BAM files to {single_variant_bam_dir}/ (variants meeting {cutoff_desc} cutoff)")
                report_lines.append(f"These BAMs contain only reads with a single specific variant and no other variants")
            else:
                report_lines.append(f"No SNVs met the single-variant cutoff of {cutoff_desc}")

        if unique_snvs:
            high_coverage_snvs = sum(1 for count in unique_snvs.values() if count > min_coverage_threshold)
            percent_high_coverage = (high_coverage_snvs / len(unique_snvs)) * 100
            report_lines.append(f"\nSNV Statistics:")
            report_lines.append(f"  - {len(unique_snvs)} unique SNVs")
            report_lines.append(f"  - {transitions} transitions ({(transitions/len(unique_snvs)*100):.2f}%)")
            report_lines.append(f"  - {transversions} transversions ({(transversions/len(unique_snvs)*100):.2f}%)")
            report_lines.append(f"  - Transition/Transversion ratio: {ti_tv_ratio:.2f}")
            report_lines.append(f"  - {high_coverage_snvs} SNVs ({percent_high_coverage:.2f}%) with coverage > {min_coverage_threshold}x")

        if unique_vr_tags:
            high_coverage_vr = sum(1 for count in unique_vr_tags.values() if count > min_coverage_threshold)
            percent_high_coverage_vr = (high_coverage_vr / len(unique_vr_tags)) * 100
            report_lines.append(f"\nVR Tag Statistics:")
            report_lines.append(f"  - {len(unique_vr_tags)} unique VR tags")
            report_lines.append(f"  - {high_coverage_vr} VR tags ({percent_high_coverage_vr:.2f}%) with coverage > {min_coverage_threshold}x")

        report_lines.append("=========================\n")

        # Save to report file if specified
        if report_output:
            with open(report_output, 'w') as f:
                for line in report_lines:
                    f.write(f"{line}\n")
            if not quiet:
                logger.info(f"Wrote summary report to {report_output}")
            
        # Print to stdout only if not quiet
        if not quiet:
            for line in report_lines:
                print(line)

        bam.close()
        output_bam.close()
        fasta.close()
        
        # Final summary message
        if not quiet:
            logger.info("Analysis complete.")

    except Exception as e:
        logger.error(f"Error in add_var_tag: {str(e)}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add SNVs, insertions, and deletions as a tag to a BAM file.")
    parser.add_argument("bam_file", type=str, help="Input BAM file")
    parser.add_argument("output_bam_file", type=str, help="Output BAM file")
    parser.add_argument("reference_fasta", type=str, help="Reference FASTA file")
    parser.add_argument("reference_chrom", type=str, help="Reference chromosome")
    parser.add_argument("start", type=int, help="Start position (1-based, inclusive)")
    parser.add_argument("end", type=int, help="End position (1-based, inclusive)")
    parser.add_argument("--snv", action="store_true", help="Include SNVs in the output tag")
    parser.add_argument("--insertion", action="store_true", help="Include insertions in the output tag")
    parser.add_argument("--deletion", action="store_true", help="Include deletions in the output tag")
    parser.add_argument("--snv-output", type=str, help="Output file for unique SNVs and their counts (default: <output_bam_prefix>_snvs.txt)")
    parser.add_argument("--vr-output", type=str, help="Output file for unique VR tags and their counts (default: <output_bam_prefix>_vr_tags.txt)")
    parser.add_argument("--report-output", type=str, help="Output file for summary report (default: <output_bam_prefix>_report.txt)")
    parser.add_argument("--write-snv-bams", action="store_true", help="Write separate BAM files for each unique SNV")
    parser.add_argument("--snv-bam-dir", type=str, default="sorted_by_SNV", help="Directory to write SNV-specific BAM files (default: sorted_by_SNV)")
    parser.add_argument("--min-coverage", type=int, default=500, help="Minimum coverage threshold for reporting (default: 500)")
    parser.add_argument("--min-snv-bam-coverage", type=int, default=1, help="Minimum coverage required to write a BAM file for an SNV (default: 1)")
    parser.add_argument("--single-variant-cutoff", type=float, help="Cutoff for filtering single-variant reads (absolute count if >= 1, proportion if 0-1)")
    parser.add_argument("--single-variant-bam-dir", type=str, help="Directory to write BAM files containing only reads with a single specific variant")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--quiet", action="store_true", help="Minimize output, show only progress bars and errors")
    parser.add_argument("--bam-log", type=str, help="File to log detailed information about each BAM file created")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel processes to use (default: 1)")
    parser.add_argument("--chunk-size", type=int, default=100, help="Number of SNVs to process per chunk in parallel mode (default: 100)")

    args = parser.parse_args()

    # Set logging level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    else:
        logging.getLogger().setLevel(logging.INFO)

    # Ensure at least one variant type is included
    if not (args.snv or args.insertion or args.deletion):
        sys.exit("Error: At least one of --snv, --insertion, or --deletion must be specified.")

    # Create default output filenames based on the output BAM file if not specified
    output_prefix = args.output_bam_file.rsplit('.', 1)[0]  # Remove extension

    if args.snv and args.snv_output is None:
        args.snv_output = f"{output_prefix}_snvs.txt"
        if not args.quiet:
            print(f"Using default SNV output file: {args.snv_output}")

    if args.vr_output is None:
        args.vr_output = f"{output_prefix}_vr_tags.txt"
        if not args.quiet:
            print(f"Using default VR tag output file: {args.vr_output}")

    if args.report_output is None:
        args.report_output = f"{output_prefix}_report.txt"
        if not args.quiet:
            print(f"Using default report output file: {args.report_output}")

    # Create default BAM log file if write-snv-bams is enabled but no log file specified
    if args.write_snv_bams and args.bam_log is None:
        args.bam_log = f"{output_prefix}_bam_log.tsv"
        if not args.quiet:
            print(f"Using default BAM log file: {args.bam_log}")

    try:
        add_var_tag(
            bam_file=args.bam_file,
            output_bam_file=args.output_bam_file,
            reference_fasta=args.reference_fasta,
            reference_chrom=args.reference_chrom,
            start=args.start,
            end=args.end,
            snv=args.snv,
            insertion=args.insertion,
            deletion=args.deletion,
            snv_output=args.snv_output,
            vr_output=args.vr_output,
            report_output=args.report_output,
            write_snv_bams=args.write_snv_bams,
            snv_bam_dir=args.snv_bam_dir,
            min_coverage_threshold=args.min_coverage,
            min_snv_bam_coverage=args.min_snv_bam_coverage,
            single_variant_cutoff=args.single_variant_cutoff,
            single_variant_bam_dir=args.single_variant_bam_dir,
            quiet=args.quiet,
            bam_log_file=args.bam_log,
            threads=args.threads,
            chunk_size=args.chunk_size
        )
    except Exception as e:
        sys.exit(f"Error: {str(e)}")
