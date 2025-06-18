from statistics import mean, median, mode
import matplotlib.pyplot as plt
import pysam
import pyft
import pandas as pd

def get_nuc_lengths(fiberbam):
    nuc_length_list = []
    read_count = 0
    for idx, fiber in enumerate(fiberbam):
        read_count += 1
        nuc_length_list.extend(fiber.nuc.lengths)

    print(f'read_count = {read_count}',f'nuc_count = {len(nuc_length_list)}, '
        f'mean = {round(mean(nuc_length_list))}, '
        f'median = {median(nuc_length_list)}, '
        f'mode = {mode(nuc_length_list)}')
    
    return nuc_length_list
    
# Need to use pysam.AlignmentFile to to get nuc lengths by nuc count
def get_nuc_lengths_by_nuc_count(samfile):
    # Ensure the 'nc' tag is present in the BAM file (e.g., added using add_nc_tag.py)
    nuc_length_dict = {}

    for idx, fiber in enumerate(samfile):
        # check if both 'nc' and 'nl' tags are present
        if fiber.has_tag('nc') and fiber.has_tag('nl'):
            # get the 'nc' and 'nl' tags
            nc = fiber.get_tag('nc')
            nl = [ item for item in fiber.get_tag('nl')]
            # Add 'nc' and 'nl' to the dictionary
            if nc not in nuc_length_dict:
                nuc_length_dict[nc] = []
            nuc_length_dict[nc].extend(nl)

    return nuc_length_dict

def get_npr(fiberbam):
    npr_list = []
    read_count = 0
    for idx, fiber in enumerate(fiberbam):
        read_count += 1
        npr_list.append(len(fiber.nuc.lengths))

    print(f'read_count = {read_count}',f'mean = {round(mean(npr_list))}',f'mode = {mode(npr_list)}')
    
    return npr_list

def get_read_lengths(fiberbam):
    read_length_list = []
    read_count = 0
    for idx, fiber in enumerate(fiberbam):
        read_count += 1
        read_length_list.append(fiber.get_seq_length())
    
    print(f'read_count = {read_count}',f'mean = {round(mean(read_length_list))}',f'median = {median(read_length_list)}',f'mode = {mode(read_length_list)}')
    
    return read_length_list

def get_percent_methylated(fiberbam):
    percent_methylated = []
    for idx, fiber in enumerate(fiberbam):
        AT_count = fiber.seq.count('A') + fiber.seq.count('T')
        m6a_count = len(fiber.m6a.starts)
        percent_methylated.append(m6a_count / AT_count)

    return percent_methylated

def get_bpn(fiberbam):
    bpn_list = []
    read_count = 0
    for idx, fiber in enumerate(fiberbam):
        read_length = fiber.get_seq_length()
        nuc_count = len(fiber.nuc.lengths)
        if nuc_count == 0:
            continue
        bpn_list.append(read_length/nuc_count)

    print(f'read_count = {read_count}',f'mean = {round(mean(bpn_list))}',f'median = {median(bpn_list)}',f'mode = {mode(bpn_list)}')
    
    return bpn_list

def get_msp_lengths(fiberbam):
    msp_length_list = []
    read_count = 0
    for idx, fiber in enumerate(fiberbam):
        msp_length_list.extend(fiber.msp.lengths)

    print(f'read_count = {read_count}',f'mean = {round(mean(msp_length_list))}',f'median = {median(msp_length_list)}',f'mode = {mode(msp_length_list)}')
    
    return msp_length_list

def get_msp_counts(fiberbam):
    msp_count_list = []
    read_count = 0
    for idx, fiber in enumerate(fiberbam):
        msp_count_list.append(len(fiber.msp.lengths))

    print(f'read_count = {read_count}',f'mean = {round(mean(msp_count_list))}',f'median = {median(msp_count_list)}',f'mode = {mode(msp_count_list)}')

    return msp_count_list

def analyze_bam(input_bam, analysis_module):
    """ 
    Runs the provided analysis module on a BAM file.

    Args:
        input_bam (str): Full path to the input BAM file.
        analysis_module (function): A function that will perform analysis on the BAM file. 
                                    Example: `get_nuc_lengths` from `pyft_analysis`.
    """   

    # Read in input_bam as pyft fiberbam object
    fiberbam = pyft.Fiberbam(input_bam)

    # Run analysis module on bam file
    result = analysis_module(fiberbam)

    return result

# Get a dataframe containing nucleosome per read data from an untagged bam file
def get_npr_df(bamfile):
    npr_df = pd.DataFrame(columns=["read_length", "nucleosome_count", "nucleosome_lengths", "mean_nucleosome_length"])
    bam_in = pysam.AlignmentFile(bamfile, "rb")
    rows = []

    for idx, read in enumerate(bam_in.fetch(until_eof=True)):
        # Skip reads without an "nl" tag
        if not read.has_tag("nl"):
            continue
        
        read_length = read.query_length
        nucleosome_lengths = list(read.get_tag("nl"))
        nucleosome_count = len(nucleosome_lengths)

        # Handle potential errors in mean calculation
        try:
            mean_nucleosome_length = round(sum(nucleosome_lengths) / len(nucleosome_lengths))
        except (ZeroDivisionError, TypeError):
            mean_nucleosome_length = "N/A"

        rows.append({
            "read_length": read_length, 
            "nucleosome_count": nucleosome_count, 
            "nucleosome_lengths": nucleosome_lengths, 
            "mean_nucleosome_length": mean_nucleosome_length
        })

    npr_df = pd.concat([npr_df, pd.DataFrame(rows)], ignore_index=True)
    return npr_df
