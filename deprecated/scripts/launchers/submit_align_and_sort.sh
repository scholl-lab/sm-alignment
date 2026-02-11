#!/bin/bash
#
#SBATCH --job-name=sm_align_and_sort_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log

# Based on:
# https://hpc-docs.cubi.bihealth.org/best-practice/temp-files/#tmpdir-and-the-scheduler
# https://bihealth.github.io/bih-cluster/slurm/snakemake/#custom-logging-directory

# First, point TMPDIR to the scratch in your home as mktemp will use this
export TMPDIR=$HOME/scratch/tmp
# Second, create another unique temporary directory within this directory
export TMPDIR=$(mktemp -d)
# Finally, setup the cleanup trap
trap "rm -rf $TMPDIR" EXIT

# Create the slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

# Export SBATCH_DEFAULTS for any sbatch commands within the script
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

# Record the start time
date

# Execute the Snakemake workflow
srun snakemake -s align_and_sort.smk --use-conda --profile=cubi-v1 -j100

# Record the end time
date
