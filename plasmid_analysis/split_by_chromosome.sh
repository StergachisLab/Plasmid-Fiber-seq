#!/bin/bash

# Check if the input BAM file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 input_bam_file"
    exit 1
fi

# Get the input BAM file from the command-line argument
input_bam="$1"

# Check if the input BAM file exists
if [ ! -f "$input_bam" ]; then
    echo "Error: Input BAM file not found."
    exit 1
fi

# Iterate over each chromosome in the BAM file
for chr in $(samtools view -H $input_bam | awk '/^@SQ/{print $2}' | sed 's/SN://'); do
    # Extract reads for each chromosome into a new BAM file
    samtools view -b $input_bam $chr > ${chr}.bam
done

