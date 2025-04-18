#!/bin/bash

#SBATCH --job-name=mean_WT_df      # Job name
#SBATCH --partition=cpu-g2         # Partition (queue) name
#SBATCH --nodes=1                  # Run on a single node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --cpus-per-task=64          # Request 16 CPU cores
#SBATCH --mem=200G                 # Memory limit
#SBATCH --time=48:00:00            # Time limit hrs:min:sec
#SBATCH --output=/gscratch/stergachislab/bmallo/large_home/slurm_jobs/log_files/fp_avg_%j.log     # Standard output and error log (%j is replaced by job ID)
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=bmallo@uw.edu  # Email address for notifications

# Print job information
echo "Job started on $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on host $(hostname)"
echo "Using $SLURM_CPUS_PER_TASK CPU cores"

# Initialize conda
eval "$(conda shell.bash hook)"

# Activate the conda environment
conda activate /mmfs1/home/bmallo/miniconda3/envs/new_general

# Verify the environment is active
echo "Using Python: $(which python)"
echo "Python version: $(python --version)"

# Set up paths and parameters
SAMPLE_PATH="/gscratch/stergachislab/bmallo/large_home/FiberHMM/FiberHMM_samples/LDLR_MPRA/SNV_sorted_reads/fp_WT_bed/LDLR_WT_fp.bed"
OUTPUT_PATH="/gscratch/stergachislab/bmallo/large_home/slurm_jobs/results/WT_footprint_avg_df_2000_reads.pkl"
REF_SEQ_LENGTH=4718
SUBSAMPLE_SIZE=2000
NUM_ITERATIONS=10000
ROW_MIN=1
ROW_MAX=400
COL_MIN=3000
COL_MAX=3600
BIN_SIZE=10

# Run the single sample average script with parallel processing on 8 cores
# Set PYTHONUNBUFFERED to ensure output is written immediately to the log file
PYTHONUNBUFFERED=1 python /gscratch/stergachislab/bmallo/large_home/python_scripts/generate_average_WT_footprint_df.py \
  --sample_path $SAMPLE_PATH \
  --output $OUTPUT_PATH \
  --ref_seq_length $REF_SEQ_LENGTH \
  --subsample_size $SUBSAMPLE_SIZE \
  --num_iterations $NUM_ITERATIONS \
  --row_min $ROW_MIN \
  --row_max $ROW_MAX \
  --col_min $COL_MIN \
  --col_max $COL_MAX \
  --bin_size $BIN_SIZE \
  --num_processes $SLURM_CPUS_PER_TASK \
  --parallel

# Send detailed completion email
if [ $? -eq 0 ]; then
    mail -s "Footprint Avg Job $SLURM_JOB_ID Completed Successfully" your.email@domain.com << EOF
Job Details:
------------
Job ID: $SLURM_JOB_ID
Completed: $(date)
Output file: $OUTPUT_PATH
CPU cores used: $SLURM_CPUS_PER_TASK
Iterations: $NUM_ITERATIONS

The analysis has completed successfully and results are available.
EOF
else
    mail -s "Footprint Avg Job $SLURM_JOB_ID FAILED" your.email@domain.com << EOF
Job Details:
------------
Job ID: $SLURM_JOB_ID
Failed: $(date)
Check log file: fp_avg_${SLURM_JOB_ID}.log for error details.
EOF
fi

# Print completion info
echo "Job finished on $(date)"
