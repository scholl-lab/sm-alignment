#!/bin/bash
#SBATCH --job-name=sm_pipeline
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log
#
# Generic Snakemake launcher for SLURM.
#
# Usage:
#   sbatch scripts/run_snakemake.sh <SNAKEFILE> [CONFIG_FILE] [PROFILE] [EXTRA_ARGS...]
#
# Examples:
#   sbatch scripts/run_snakemake.sh workflow/Snakefile
#   sbatch --job-name=sm_trim scripts/run_snakemake.sh workflow/Snakefile config/config_trim.yaml
#   sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml cubi-v1 --forceall

set -euo pipefail

# ---- Arguments with defaults ----
SNAKEFILE="${1:?Error: SNAKEFILE path required as first argument}"
CONFIG_FILE="${2:-config/config.yaml}"
PROFILE="${3:-cubi-v1}"
shift 3 2>/dev/null || true

# ---- TMPDIR setup ----
BASE_TMPDIR="${HOME}/scratch/tmp"
mkdir -p "${BASE_TMPDIR}"
export TMPDIR=$(mktemp -d "${BASE_TMPDIR}/sm.XXXXXX")
trap 'rm -rf "${TMPDIR}"' EXIT

# ---- SLURM logging ----
mkdir -p slurm_logs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

# ---- Launch ----
echo "=== Snakemake Launch ==="
echo "  Snakefile:  ${SNAKEFILE}"
echo "  Config:     ${CONFIG_FILE}"
echo "  Profile:    ${PROFILE}"
echo "  TMPDIR:     ${TMPDIR}"
echo "  Extra args: $*"
echo "  Start:      $(date)"
echo "========================"

srun snakemake \
    -s "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile profiles/default \
    --profile="${PROFILE}" \
    "$@"

echo "=== Finished: $(date) ==="
