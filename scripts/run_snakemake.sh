#!/bin/bash
#SBATCH --job-name=sm_pipeline
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_logs/%x-%j.log
#
# Snakemake SLURM launcher — auto-detects BIH and Charité HPC clusters.
#
# Usage:
#   sbatch scripts/run_snakemake.sh <SNAKEFILE> [CONFIG_FILE] [EXTRA_ARGS...]
#
# Examples:
#   sbatch scripts/run_snakemake.sh workflow/Snakefile
#   sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml
#   sbatch --job-name=sm_exomes scripts/run_snakemake.sh workflow/Snakefile
#   sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml --forceall

set -euo pipefail

# ---- Cluster auto-detection ----
detect_cluster() {
    local fqdn
    fqdn=$(hostname -f 2>/dev/null || hostname)
    if [[ -d "/etc/xdg/snakemake/cubi-v1" ]] || [[ "$fqdn" =~ cubi|bihealth ]]; then
        echo "bih"
    elif [[ -f "/etc/profile.d/conda.sh" ]] || [[ "$fqdn" =~ charite|\.sc- ]]; then
        echo "charite"
    else
        echo "local"
    fi
}

CLUSTER=$(detect_cluster)

# ---- Arguments with defaults ----
SNAKEFILE="${1:?Error: SNAKEFILE path required as first argument}"
CONFIG_FILE="${2:-config/config.yaml}"
shift 2 2>/dev/null || shift 1 2>/dev/null || true

# ---- Conda activation (Charité requires explicit sourcing) ----
if [[ "$CLUSTER" == "charite" ]] && [[ -f /etc/profile.d/conda.sh ]]; then
    source /etc/profile.d/conda.sh
fi
conda activate snakemake 2>/dev/null || true

# ---- TMPDIR setup ----
if [[ "$CLUSTER" == "bih" ]]; then
    BASE_TMPDIR="${HOME}/scratch/tmp"
else
    BASE_TMPDIR="${TMPDIR:-/tmp}/snakemake"
fi
mkdir -p "${BASE_TMPDIR}"
TMPDIR=$(mktemp -d "${BASE_TMPDIR}/sm.XXXXXX")
export TMPDIR
trap 'rm -rf "${TMPDIR}"' EXIT

# ---- SLURM logging ----
mkdir -p slurm_logs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

# ---- Cluster-specific Snakemake arguments ----
CLUSTER_ARGS=()
if [[ "$CLUSTER" == "bih" ]]; then
    CLUSTER_ARGS+=(--profile=cubi-v1)
elif [[ "$CLUSTER" == "charite" ]]; then
    CLUSTER_ARGS+=(--profile profiles/charite)
fi

# ---- Launch ----
echo "=== Snakemake Launch ==="
echo "  Cluster:    ${CLUSTER}"
echo "  Snakefile:  ${SNAKEFILE}"
echo "  Config:     ${CONFIG_FILE}"
echo "  TMPDIR:     ${TMPDIR}"
echo "  Extra args: $*"
echo "  Start:      $(date)"
echo "========================"

snakemake \
    -s "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile profiles/default \
    "${CLUSTER_ARGS[@]}" \
    "$@"

echo "=== Finished: $(date) ==="
