import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import sys
import ast
import tqdm.notebook
from scipy.ndimage import gaussian_filter1d

def blocks_to_array(row):
    #turn sizes and starts into alternating footprint length and gap length
    #encoded as -10, 40, -30, 90, -32, etc.
    fp_sizes=ast.literal_eval(row['blockSizes'])
    fp_starts=ast.literal_eval(row['blockStarts'])
    
    # check if it's a single footprint, which won't be read by ast correctly
    if isinstance(fp_sizes, int):
        fp_sizes=[fp_sizes]
        fp_starts=[fp_starts]
        
    fp_starts=np.array(fp_starts)
    fp_sizes=np.array(fp_sizes)

    length=row['end']-row['start']

    #calculate gap sizes and starts
    gap_starts=np.append([0],fp_starts+fp_sizes)
    gap_sizes=np.append(fp_starts,[length])
    gap_sizes=-(gap_sizes-gap_starts)

    #add to an array
    read=np.zeros(len(gap_sizes)+len(fp_sizes)).astype(int)    
    read[1::2] = fp_sizes
    read[0::2] = gap_sizes
    read=read[read!=0]
    
    return read

def unpack(read):
    #unpack compressed read into expanded read
    #from -1, 3, -2, 1, etc. to -1, 3, 3, 3, -2, -2, 1, etc.
    return np.repeat(read, np.abs(read))


def pack(read):
    #convert expanded to compressed version of read
    #from -1, 3, 3, 3, -2, -2, 1, etc. to -1, 3, -2, 1, etc.
    return read[np.insert(np.diff(read) != 0, 0, True)]

def slice_read(read, start, end, strand):
    # for a given read, cut to the dimensions in query
    # the input is the expanded version of the read (i.e. unpack())
    
   # # add 2s to encode portions of read not within query
    if start < 0:
        pre_read = np.full(-start, 2)
    else:
        pre_read = np.array([])
        
    if end > 0:
        post_read = np.full(end, 2)
    else:
        post_read = np.array([]) 
        
    # check if either is outside of query
    start = max(start, 0)
    end = min(end, 0)
    if end == 0:
        end = len(read)
   # print(start,end,len(read))

    # trim read
    read = read[start:end]
    
    # add 2s to outer part
    read = np.concatenate((pre_read, read, post_read))
    if strand=='-':
        read=read[::-1]

    return read.astype(int)
    
def grab_reads(reads, query):
    #given a dataframe of reads (the bed file dataframe) and a set of query coordinates
    #return a dataframe of reads with overlap (and the relative coordinates)
    
    #convert columns into arrays
    #reads
    chrom_reads = reads['chrom'].values
    start_reads = reads['start'].values
    end_reads = reads['end'].values
    rid_reads = reads['name'].values
    blockStarts = dict(zip(reads['name'],reads['blockStarts']))
    blockSizes = dict(zip(reads['name'],reads['blockSizes']))

    #query
    chrom_query = query['chrom'].values
    start_query = query['start'].values
    end_query = query['end'].values
    rid_query = query['name'].values
    strand = dict(zip(query['name'],query['strand']))

    overlaps = []
    
    for chrom in tqdm.notebook.tqdm(np.unique(chrom_reads), leave=False):
        # find indices for current chromosome in both bed files
        idxr = np.where(chrom_reads == chrom)[0]
        idxq = np.where(chrom_query == chrom)[0]
        
        if len(idxr) == 0 or len(idxq) == 0:
            continue  
        
        # filter start and end positions for current chromosome
        start_reads_chrom = start_reads[idxr]
        end_reads_chrom = end_reads[idxr]
        start_query_chrom = start_query[idxq]
        end_query_chrom = end_query[idxq]
        
        # find entries where query is contained within reads
        overlap_matrix = (
            (start_reads_chrom[:, None] <= end_query_chrom) & 
            (end_reads_chrom[:, None] >= start_query_chrom)
        )
        
        # grab overlapping indices
        overlap_idxr, overlap_idxq = np.nonzero(overlap_matrix)

        # grab overlapping rows
        overlaps.extend(
            zip(
                chrom_reads[idxr[overlap_idxr]],
                start_reads_chrom[overlap_idxr],
                end_reads_chrom[overlap_idxr],
                rid_reads[idxr[overlap_idxr]],
                
                start_query_chrom[overlap_idxq],
                end_query_chrom[overlap_idxq],
            
                rid_query[idxq[overlap_idxq]]
            )
        )

    # combine columns again
    overlap_df = pd.DataFrame(overlaps, columns=['chrom', 'start', 'end', 'name', 'q_start', 'q_end', 'q_id'])
    overlap_df['blockStarts'] = overlap_df['name'].map(blockStarts)
    overlap_df['blockSizes'] = overlap_df['name'].map(blockSizes)
    overlap_df['rel_start'] = overlap_df['q_start']-overlap_df['start']
    overlap_df['rel_end'] = overlap_df['q_end']-overlap_df['end']
    overlap_df['strand'] = overlap_df['q_id'].map(strand)
    return overlap_df

def read_ft_data_single(bed_file):
        df = pd.read_csv(bed_file,
            sep='\t',
            header=None,
            usecols=[0, 1, 2, 3, 8, 9]
        )  # Read the BED file
        df = df.reset_index(drop=True)
        df.columns = ['chrom', 'start', 'end', 'name', 'blockStarts', 'blockSizes']
        df = df.loc[df['blockStarts'] != '.']  # Filter rows where blockStarts is not '.'

        # Remove duplicate reads with worse alignment
        df['length'] = df['end'] - df['start']
        df = df.sort_values('length', ascending=False).drop_duplicates('name')

        return df  # Return the bed file as a dataframe

def read_ft_data(sample_directory_path):
    bed_dict = {}  # Initialize a dictionary to store DataFrames
    for f in os.listdir(sample_directory_path):
        if 'fp.bed' in f:
            bed_name = '_'.join(f.split('_')[:-1])  # Extract the key for the dictionary
            bed = pd.read_csv(
                os.path.join(sample_directory_path, f),
                sep='\t',
                header=None,
                usecols=[0, 1, 2, 3, 8, 9]
            )  # Read the BED file
            bed = bed.reset_index(drop=True)
            bed.columns = ['chrom', 'start', 'end', 'name', 'blockStarts', 'blockSizes']
            bed = bed.loc[bed['blockStarts'] != '.']  # Filter rows where blockStarts is not '.'

            # Remove duplicate reads with worse alignment
            bed['length'] = bed['end'] - bed['start']
            bed = bed.sort_values('length', ascending=False).drop_duplicates('name')

            # Save the DataFrame to the dictionary
            bed_dict[bed_name] = bed

    return bed_dict  # Return the dictionary containing all DataFrames

def grab_circular_reads(sample_name, ref_seq_length, number_of_reads=None):
    """
    Extract circular reads from a sample DataFrame.
    
    Parameters:
    sample_name (pd.DataFrame): DataFrame containing read data.
    ref_seq_length (int): Length of the reference sequence.
    number_of_reads (int, optional): Number of reads to sample. If None, use all reads.
    
    Returns:
    pd.DataFrame: DataFrame with circular reads.
    """
    if number_of_reads is None:
        number_of_reads = sample_name.shape[0]
    tmp = sample_name.sample(n=number_of_reads)  # subset the bed as desired

    # initialize reads_df
    reads_df = []
    # iterate, unpack, trim
    for index, row in tmp.iterrows():  # Removed tqdm wrapper
        read = unpack(blocks_to_array(row))
        read = read[ref_seq_length:ref_seq_length*2]
        reads_df.append(read)
    # turn into dataframe
    reads_df = pd.DataFrame(reads_df)
    return reads_df

def prep_dfs_for_subtraction(*dataframes):
    """
    Processes multiple DataFrames to ensure rows cover a continuous range of values, 
    filling gaps with 0, and removes any rows with negative indices.

    Parameters:
    - *dataframes (pd.DataFrame): One or more DataFrames to process.

    Returns:
    - list of pd.DataFrame: Processed DataFrames, reindexed and aligned.
    """
    # Helper function to process a single DataFrame
    def process_dataframe(df):
        # Initialize an empty dictionary to store value counts
        value_counts_dict = {}
        
        # Iterate over each column in the original DataFrame
        for col in df.columns:
            # Get the counts of unique values in the column
            value_counts = df[col].value_counts()
            # Add the counts as a Series to the dictionary
            value_counts_dict[col] = value_counts
        
        # Convert the dictionary to a DataFrame, where rows are footprint sizes, columns are base positions
        new_df = pd.DataFrame(value_counts_dict)
        
        # Fill NaN values with 0
        new_df.fillna(0, inplace=True)
        
        return new_df

    # Process all DataFrames
    processed_dfs = [process_dataframe(df) for df in dataframes]

    # Determine the range of the new index across all DataFrames
    min_value = 0  # Start at 0 to ensure no negative indices
    max_value = max(df.index.max() for df in processed_dfs)
    max_value = int(max_value)  # Ensure the maximum value is an integer

    # Create a new index covering the entire range
    new_index = range(min_value, max_value + 1)

    # Reindex all DataFrames to fill gaps with 0 and explicitly remove negative indices
    aligned_dfs = [
        df[df.index >= 0].reindex(new_index, fill_value=0).loc[new_index] 
        for df in processed_dfs
    ]

    return aligned_dfs

def filter_fp_df(df, bin_size, row_range=None, col_range=None):
    """
    Processes a footprint DataFrame that has been preped for subtraction
    and optionally selects specific rows and columns to plot on the hatmap.

    Parameters:
    df (pd.DataFrame): The input DataFrame where rows are footprint sizes, columns are base positions,
                       and values are the number of footprints of that size over each base position.
    bin_size (int): The size of bins to group footprint sizes.
    row_range (tuple, optional): A tuple (start, end) specifying the range of rows (index values) to include.
                                 If None, all rows are included.
    col_range (tuple, optional): A tuple (start, end) specifying the range of columns to include.
                                 If None, all columns are included.

    Returns:
    pd.DataFrame: A DataFrame with rows binned by the specified bin size, and
                  selected rows and columns if specified.
    """
    # Subset the DataFrame based on row_range and col_range if provided
    if row_range is not None:
        df = df[(df.index >= row_range[0]) & (df.index <= row_range[1])]
    if col_range is not None:
        df = df.loc[:, col_range[0]:col_range[1]]
    
    # If no bin_size is provided, return the filtered DataFrame without binning
    if bin_size is None or bin_size <= 0:
        return df
    
    # Determine min and max values for binning
    min_val = row_range[0] if row_range is not None else df.index.min() 
    max_val = row_range[1] if row_range is not None else df.index.max()
    
    # Ensure min_val is rounded down to nearest bin_size multiple
    min_bin = (min_val // bin_size) * bin_size
    
    # Ensure max_val is rounded up to nearest bin_size multiple + bin_size
    max_bin = ((max_val // bin_size) + 1) * bin_size
    
    # Create bin boundaries and labels
    bins = range(min_bin, max_bin + 1, bin_size)
    labels = [f'{i+1}-{i+bin_size}' for i in range(min_bin, max_bin, bin_size)]
    
    # Create bin labels for each row
    bin_labels = pd.cut(df.index, bins=bins, labels=labels)
    
    # Group by the bin labels and sum
    binned_df = df.groupby(bin_labels).sum()
    
    return binned_df

def plot_fp_heatmap(df, snp_location_1=100, snp_location_2=200, vmin=None, vmax=None, cmap="viridis", title="", xticks=0, output_directory=None, legend_label=None):
    """
    Creates a heatmap for footprint enrichment data with a vertical line indicating SNP location.

    Parameters:
    df (pd.DataFrame): Data for the heatmap, where rows represent footprint sizes and columns represent base positions.
    snp_location (int): X-coordinate for the SNP location vertical line. Default is vmin.
    vmin (float): Minimum value for the heatmap color scale.
    vmax (float): Maximum value for the heatmap color scale.
    cmap (str): Colormap for the heatmap. Default is "viridis".
    """
    # Set up figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(df, ax=ax, cbar=True, annot=False, fmt='.2f', vmin=vmin, vmax=vmax, cmap=cmap, xticklabels=xticks)

    # Customize colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_label("Footprint Enrichment on WT", fontsize=14)

    # Add vertical line for SNP location
    ax.axvline(x=snp_location_1, color='white', linestyle='--', linewidth=1, 
               label='SNP Location', alpha=0.5)
    
    # Add vertical line for second SNP location
    ax.axvline(x=snp_location_2, color='white', linestyle='--', linewidth=1, 
               label='SNP Location', alpha=0.5)

    # Set title and labels
    ax.set_title(f'{title}', 
                 fontsize=16, fontweight='bold')
    ax.set_ylabel("Footprint Size", fontsize=14)

    # Add legend
    legend = ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.12), 
                       labels=[legend_label], facecolor='dimgray', edgecolor='black')
    for text in legend.get_texts():
        text.set_color("white")

    if output_directory:
        save_directory = os.path.join(
            '/gscratch/stergachislab/bmallo/large_home/ft_data/plasmid_fiberseq/time_course/figures/FiberHMM_heatmaps',
            output_directory
        )
        if not os.path.exists(save_directory):
            os.makedirs(save_directory)
        
        save_path = os.path.join(save_directory, f"{title}.pdf")
        plt.savefig(save_path, dpi=300)
        print(f"Figure saved to {save_path}")

    plt.show()

def meta_line_plot(reads_df, fp_size_min, fp_size_max, smoothing=None):
    """
    Counts the frequency of a given size range of footprints at positions in a DataFrame of unpacked reads.
    Reads must be trimmed and centered in the DataFrame. Output is an array with fractional enrichment.

    Parameters:
        reads_df (pd.DataFrame): DataFrame of unpacked reads.
        fp_size_min (int): Minimum size of footprints to include.
        fp_size_max (int): Maximum size of footprints to include.
        smoothing (float, optional): Standard deviation for Gaussian smoothing. Default is None.

    Returns:
        np.ndarray: Array with fractional enrichment.
    """
    # This accounts for the "2" values masking the ends of reads
    count_df = reads_df.copy()
    count_df[count_df != 2] = 1
    count_df[count_df == 2] = 0
    counts = count_df.sum(axis=0)
    
    df = reads_df.copy()
    df[df <= 2] = 0

    if fp_size_min:
        df[df < fp_size_min] = 0
    if fp_size_max:
        df[df > fp_size_max] = 0
    df[df > 0] = 1

    df_a = df.sum(axis=0)
    df_a = np.array(df_a) / counts.to_numpy()

    # Apply smoothing if specified
    if smoothing is not None:
        df_a = gaussian_filter1d(df_a, sigma=smoothing)

    return df_a