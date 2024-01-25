#!/bin/bash

# Check if the required command line arguments are provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 experiment_name1 [experiment_name2 ...]"
    exit 1
fi

# Iterate over each experiment name provided in the command line arguments
for experiment_name in "$@"; do
    # Initialize the command line for each experiment
    command_line=("python" "/gscratch/stergachislab/bmallo/large_home/hmm-footprint-caller/multi_cell_line/apply_model_exact.py" "-i" "$experiment_name" "-e" "multi_cell_hg38.h5" "-f")

    # Set the sample directory based on the experiment name
    sample_directory="$experiment_name/Infiles/parquet"

    # Iterate over the files in the directory and add file as sample
    for sample_file in "$sample_directory"/*; do
        sample_names+=("$(basename "$sample_file")")
    done

    # Use IFS to join sample names with commas
    IFS=',' file_name_str="${sample_names[*]}"

    # Add formatted sample names to command line code
    command_line+=("$file_name_str")

    # Print the command line for the current experiment
    echo "Running command for experiment $experiment_name: ${command_line[@]}"

    # Run the command
    "${command_line[@]}"

    # Clear sample_names array for the next iteration
    sample_names=()
done

