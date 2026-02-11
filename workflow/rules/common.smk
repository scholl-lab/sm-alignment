import os
import pandas as pd


# =============================================================================
# Config shortcuts
# =============================================================================
REF = config["ref"]["genome"]
REF_GZ = config["ref"].get("genome_gz", "")
REF_BUILD = config["ref"]["build"]
KNOWN_SITES = config["ref"]["known_sites"]

FASTQ_DIR = config["paths"]["fastq_folder"]
OUTPUT_DIR = config["paths"]["output_folder"]
LOG_SUBDIR = config["paths"].get("log_subdir", "logs")
LOG_DIR = os.path.join(OUTPUT_DIR, LOG_SUBDIR)

PLATFORM = config.get("read_group", {}).get("platform", "ILLUMINA")

R1_SUFFIX = config.get("fastq", {}).get("trimmed_r1_suffix", ".bbduk_R1_001.fastq.gz")
R2_SUFFIX = config.get("fastq", {}).get("trimmed_r2_suffix", ".bbduk_R2_001.fastq.gz")

MERGED_SUFFIX = config.get("bam", {}).get("merged_suffix", ".merged.bam")
DEDUP_SUFFIX = config.get("bam", {}).get("dedup_suffix", ".merged.dedup.bam")
DEDUP_METRICS_SUFFIX = config.get("bam", {}).get("dedup_metrics_suffix", ".merged.dedup_metrics.txt")
RECAL_TABLE_SUFFIX = config.get("bam", {}).get("recal_table_suffix", ".merged.dedup.recal_data.table")
FINAL_BAM_SUFFIX = config.get("bam", {}).get("final_suffix", ".merged.dedup.bqsr.bam")

COMPRESSION_LEVEL = config.get("processing", {}).get("compression_level", 6)

ALIGNED_DIR = os.path.join(OUTPUT_DIR, "aligned")
MERGED_DIR = os.path.join(OUTPUT_DIR, "merged")
DEDUP_DIR = os.path.join(OUTPUT_DIR, "dedup")
BQSR_DIR = os.path.join(OUTPUT_DIR, "bqsr")


# =============================================================================
# Samples metadata
# =============================================================================
samples_df = pd.read_table(config["paths"]["samples"]).set_index(
    "fastq_files_basename", drop=False
)


def get_samples():
    """Return sorted list of unique sample names from metadata."""
    return sorted(samples_df["project_sample"].unique().tolist())


def get_basenames_for_sample(sample):
    """Return all FASTQ basenames belonging to a given sample."""
    return samples_df.loc[
        samples_df["project_sample"] == sample, "fastq_files_basename"
    ].tolist()


# =============================================================================
# Resource helpers
# =============================================================================
def get_java_opts(wildcards, resources):
    """
    Derive GATK --java-options from allocated resources.
    Reserves 20% of mem_mb for JVM non-heap overhead.
    """
    xmx = int(resources.mem_mb * 0.8)
    xms = int(resources.mem_mb * 0.2)
    tmpdir = resources.tmpdir
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={tmpdir}"


# =============================================================================
# Ensure output directories exist
# =============================================================================
for _d in [ALIGNED_DIR, MERGED_DIR, DEDUP_DIR, BQSR_DIR, LOG_DIR]:
    os.makedirs(_d, exist_ok=True)
