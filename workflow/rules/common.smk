import os
import pandas as pd
from rules.helpers import (
    get_java_opts as _get_java_opts_impl,
    get_samples as _get_samples_impl,
    get_all_basenames as _get_all_basenames_impl,
    get_basenames_for_sample as _get_basenames_impl,
    resolve_fastq_path as _resolve_fastq_path_impl,
)


# =============================================================================
# Config shortcuts
# =============================================================================
REF = config["ref"]["genome"]
REF_GZ = config["ref"].get("genome_gz", "") or REF
REF_BUILD = config["ref"]["build"]
KNOWN_SITES = config["ref"]["known_sites"]

FASTQ_DIR = config["paths"]["fastq_folder"]
OUTPUT_DIR = config["paths"]["output_folder"]
LOG_SUBDIR = config["paths"].get("log_subdir", "logs")
LOG_DIR = os.path.join(OUTPUT_DIR, LOG_SUBDIR)

PLATFORM = config.get("read_group", {}).get("platform", "ILLUMINA")

# --- FASTQ suffixes (raw and trimmed) ---
R1_SUFFIX = config.get("fastq", {}).get("r1_suffix", "_R1_001.fastq.gz")
R2_SUFFIX = config.get("fastq", {}).get("r2_suffix", "_R2_001.fastq.gz")
TRIMMED_R1_SUFFIX = config.get("fastq", {}).get("trimmed_r1_suffix", ".bbduk_R1_001.fastq.gz")
TRIMMED_R2_SUFFIX = config.get("fastq", {}).get("trimmed_r2_suffix", ".bbduk_R2_001.fastq.gz")

# --- Trimming ---
TRIMMING_ENABLED = config.get("trimming", {}).get("enabled", False)
TRIMMED_DIR = os.path.join(OUTPUT_DIR, "bbduk_trimmed")

MERGED_SUFFIX = config.get("bam", {}).get("merged_suffix", ".merged.bam")
DEDUP_SUFFIX = config.get("bam", {}).get("dedup_suffix", ".merged.dedup.bam")
DEDUP_METRICS_SUFFIX = config.get("bam", {}).get(
    "dedup_metrics_suffix", ".merged.dedup_metrics.txt"
)
RECAL_TABLE_SUFFIX = config.get("bam", {}).get(
    "recal_table_suffix", ".merged.dedup.recal_data.table"
)
FINAL_BAM_SUFFIX = config.get("bam", {}).get("final_suffix", ".merged.dedup.bqsr.bam")

COMPRESSION_LEVEL = config.get("processing", {}).get("compression_level", 6)

ALIGNED_DIR = os.path.join(OUTPUT_DIR, "aligned")
MERGED_DIR = os.path.join(OUTPUT_DIR, "merged")
DEDUP_DIR = os.path.join(OUTPUT_DIR, "dedup")
BQSR_DIR = os.path.join(OUTPUT_DIR, "bqsr")

# --- Quality control ---
QC_CFG = config.get("qc", {})
QC_ENABLED = QC_CFG.get("enabled", False)
QC_DIR = os.path.join(OUTPUT_DIR, "qc")


# =============================================================================
# Samples metadata
# =============================================================================
samples_df = pd.read_table(config["paths"]["samples"]).set_index(
    "fastq_files_basename", drop=False
)


def get_samples():
    """Return sorted list of unique sample names from metadata."""
    return _get_samples_impl(samples_df)


def get_all_basenames():
    """Return sorted list of all FASTQ basenames from metadata."""
    return _get_all_basenames_impl(samples_df)


def get_basenames_for_sample(sample):
    """Return all FASTQ basenames belonging to a given sample."""
    return _get_basenames_impl(samples_df, sample)


def _resolve_fastq_path(wildcards, suffix):
    """Build FASTQ path respecting trimming mode and optional subfolder."""
    return _resolve_fastq_path_impl(
        wildcards.basename, suffix, samples_df, TRIMMING_ENABLED, TRIMMED_DIR, FASTQ_DIR
    )


def get_fastq_r1(wildcards):
    """Return R1 FASTQ path - trimmed or raw depending on config."""
    suffix = TRIMMED_R1_SUFFIX if TRIMMING_ENABLED else R1_SUFFIX
    return _resolve_fastq_path(wildcards, suffix)


def get_fastq_r2(wildcards):
    """Return R2 FASTQ path - trimmed or raw depending on config."""
    suffix = TRIMMED_R2_SUFFIX if TRIMMING_ENABLED else R2_SUFFIX
    return _resolve_fastq_path(wildcards, suffix)


# =============================================================================
# Resource helpers
# =============================================================================
def get_java_opts(wildcards, resources):
    """
    Derive GATK --java-options from allocated resources.
    Reserves 20% of mem_mb for JVM non-heap overhead.
    """
    return _get_java_opts_impl(resources.mem_mb, resources.tmpdir)
