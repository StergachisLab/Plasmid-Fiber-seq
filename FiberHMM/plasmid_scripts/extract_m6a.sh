#!/bin/bash

# Check if the correct number of command-line arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Set the input and output directory paths from the command-line arguments
input_directory="$1"
output_directory="$2"

# Check if the specified input directory exists
if [ ! -d "$input_directory" ]; then
    echo "Error: Input directory not found - $input_directory"
    exit 1
fi

# Check if the specified output directory exists; if not, create it
if [ ! -d "$output_directory" ]; then
    mkdir -p "$output_directory"
fi

# Iterate through each file in the input directory
for input_file in "$input_directory"/*_nucleosomes.bam; do
    # Extract the number from the input file name
    file_number=$(basename "$input_file" | cut -d'_' -f1)

    # Define the output file name with the specified output directory
    output_file="${output_directory}/${file_number}_nucleosomes.m6a.bed.gz"

    # Execute the ft extract command
    ft extract -r -v --m6a "$output_file" "$input_file"

    # Optionally, you can print a message indicating the completion of each file
    echo "Processed $input_file and saved result as $output_file"
done

