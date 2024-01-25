import pandas as pd
import numpy as np
import os
import argparse
import sys
from matplotlib import pyplot as plt

# this converts bed files to split parquet files which minimized memory impact of the downstream steps
# only takes an input of the directory containing the bed files, and the chromosome info file

def options():

    chromlist = []

    parser = argparse.ArgumentParser(description='process options')
    parser.add_argument('-i', '--infile', dest='indir', required=True,
                        help='Input directory')
    parser.add_argument('-c', '--chrom_info', dest='chrom_info', required=True,
                        help='TSV file with chromosome names in first column')
    args = parser.parse_args()

    indir = args.indir
    chrom_info = args.chrom_info

    chrom_info_path = os.path.join(indir, "Reference", chrom_info)

    with open(chrom_info_path, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            chromlist.append(line[0])

    return indir, chromlist

def bed_to_parquet(infile, outdir, prefix, chromlist):
    # makes a directory
    os.system('mkdir '+outdir)

    # reads in the reads as a dataframe, columns could be wrong down the line if the m6a caller changes or something
    reads = pd.read_csv(infile, usecols=[0, 1, 2, 3, 11], names=[
                        'chrom', 'start', 'end', 'rid', 'me'], sep='\t')
    reads = reads.loc[reads['me'] != '.']
    reads = reads.sort_values(by=['chrom', 'start'])
    reads = reads.reset_index(drop=True)
    reads.index = reads.index.astype(str)

    # only grab the chromosomes we actually want, I'll add an option to input a chrominfo file
    for chrom in chromlist:
        tmp = reads.loc[reads['chrom'] == chrom]
        tmp = tmp['me'].str.split(pat=',', expand=True).T
        tmp.to_parquet(outdir+'/'+prefix+'_'+chrom+'.pq')

    # write out to parquet
    reads = reads.drop(columns=['me'])
    reads = reads.loc[reads['chrom'].isin(chromlist)]
    reads.to_parquet(outdir+'/'+prefix+'_read-info.pq')


indir, chromlist = options()

indir = indir+'/Infiles/bed'

# iterates through all bed files

parquet_dir = indir.replace('bed', 'parquet')
context_dir = indir.replace('bed', 'context')

if not os.path.exists(parquet_dir):
    os.mkdir(parquet_dir)
else:
    print(f"The directory {parquet_dir} already exists, skipping creation.")

if not os.path.exists(context_dir):
    os.mkdir(context_dir)
else:
    print(f"The directory {context_dir} already exists, skipping creation.")

for f in os.listdir(indir):
    print(f)
    infile = indir+'/'+f
    prefix = f.replace('.bed', '')
    outdir1 = indir.replace('bed', 'parquet')+'/'+prefix
    print('converting to split methylation parquet file')
    bed_to_parquet(infile, outdir1, prefix, chromlist)
