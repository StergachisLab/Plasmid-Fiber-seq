#!/bin/bash
#SBATCH --job-name=precompute_null_stats
#SBATCH --output=/gscratch/stergachislab/bmallo/large_home/slurm_jobs/log_files/precompute_null_stats_%j.out
#SBATCH --error=/gscratch/stergachislab/bmallo/large_home/slurm_jobs/log_files/precompute_null_stats_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=cpu-g2

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_JOB_NODELIST"
echo "Start Time: $(date)"

# Initialize conda (may be required on some clusters)
eval "$(conda shell.bash hook)"

# Activate the conda environment
conda activate /mmfs1/home/bmallo/miniconda3/envs/general

# Define input and output paths
NULL_DIST_PATH="/gscratch/stergachislab/bmallo/large_home/slurm_jobs/results/null_WT_distribution_2000_reads.pkl"
OUTPUT_DIR="/gscratch/stergachislab/bmallo/large_home/slurm_jobs/results"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the precomputation script
echo "Starting precomputation..."
python /gscratch/stergachislab/bmallo/large_home/git_repos/Plasmid-Fiber-seq/MPRA_footprint_analysis/precompute_null_stats.py "$NULL_DIST_PATH" --output-dir "$OUTPUT_DIR"

# Check if the script executed successfully
if [ $? -eq 0 ]; then
    echo "Precomputation completed successfully."
else
    echo "Precomputation failed. Check error logs."
    exit 1
fi

# Print completion information
echo "End Time: $(date)"
echo "Job completed."

# Optional: Save environment and module information for reproducibility
module list
python --version
pip list > "${OUTPUT_DIR}/pip_packages_${SLURM_JOB_ID}.txt"

exit 0
