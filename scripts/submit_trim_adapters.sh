#!/bin/bash
#
#SBATCH --job-name=sm_trim_adapters_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log

#
# Usage:
#   sbatch submit_trim_adapters.sh [CONFIG_FILE] [MAX_JOBS] [SNAKEMAKE_FILE]
#
#   - CONFIG_FILE:      Path to a Snakemake config file (default: "config.yaml")
#   - MAX_JOBS:         Number of Snakemake jobs (in parallel) to use (default: 20)
#   - SNAKEMAKE_FILE:   Path to the Snakemake workflow file (default: "trim_adapters.smk")
#
# ------------------------------------------------------------------------------------
# Parse command-line arguments with defaults
# ------------------------------------------------------------------------------------
CONFIG_FILE=${1:-"config.yaml"}
MAX_JOBS=${2:-20}
SNAKEMAKE_FILE=${3:-"trim_adapters.smk"}

# ------------------------------------------------------------------------------------
# HPC environment setup
# ------------------------------------------------------------------------------------
export TMPDIR=$HOME/scratch/tmp
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Create the slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

# Export default SBATCH outputs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date

# ------------------------------------------------------------------------------------
# Run Snakemake via srun
# ------------------------------------------------------------------------------------
srun snakemake \
    -s "$SNAKEMAKE_FILE" \
    --use-conda \
    --profile=cubi-v1 \
    -j "$MAX_JOBS" \
    --configfile "$CONFIG_FILE"

date
