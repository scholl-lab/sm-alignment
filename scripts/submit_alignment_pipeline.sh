#!/bin/bash
#
#SBATCH --job-name=sm_alignment_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log

#
# Usage:
#   sbatch submit_alignment_pipeline.sh [SNAKEMAKE_FILE] [CONFIG_FILE] [MAX_JOBS]
#
#   - SNAKEMAKE_FILE:   Path to the Snakemake workflow file (default: "alignment_pipeline.smk")
#   - CONFIG_FILE:      Path to a Snakemake config file (default: "config.yaml")
#   - MAX_JOBS:         Number of Snakemake jobs (in parallel) to use (default: 20)
#

# ------------------------------------------------------------------------------------
# Parse command-line arguments with defaults
# ------------------------------------------------------------------------------------
SNAKEMAKE_FILE=${1:-"alignment_pipeline.smk"}
CONFIG_FILE=${2:-"config.yaml"}
MAX_JOBS=${3:-20}

# ------------------------------------------------------------------------------------
# HPC environment setup
# ------------------------------------------------------------------------------------
export TMPDIR="$HOME/scratch/tmp"
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Create the slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

echo "DEBUG: Running Snakemake pipeline with:"
echo "       Snakefile:    $SNAKEMAKE_FILE"
echo "       Config file:  $CONFIG_FILE"
echo "       Max jobs:     $MAX_JOBS"
echo "       TMPDIR:       $TMPDIR"

date

# ------------------------------------------------------------------------------------
# Run Snakemake via srun (or directly if your cluster allows)
# ------------------------------------------------------------------------------------
srun snakemake \
    -s "$SNAKEMAKE_FILE" \
    --use-conda \
    --profile=cubi-v1 \
    -j "$MAX_JOBS" \
    --configfile "$CONFIG_FILE"

date
