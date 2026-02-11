"""Pure-Python helper functions for the sm-alignment pipeline.

Extracted from common.smk so they can be unit-tested independently.
"""

from __future__ import annotations

import os


def get_java_opts(resources_mem_mb: int, tmpdir: str) -> str:
    """Derive GATK --java-options from allocated resources.

    Reserves 20% of mem_mb for JVM non-heap overhead.

    Args:
        resources_mem_mb: Allocated memory in megabytes.
        tmpdir: Temporary directory path.

    Returns:
        Java options string for -Xms, -Xmx, and tmpdir.
    """
    xmx = int(resources_mem_mb * 0.8)
    xms = int(resources_mem_mb * 0.2)
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={tmpdir}"


def get_samples(samples_df) -> list[str]:
    """Return sorted list of unique sample names from metadata.

    Args:
        samples_df: DataFrame with a 'project_sample' column.

    Returns:
        Sorted list of unique sample names.
    """
    return sorted(samples_df["project_sample"].unique().tolist())


def get_basenames_for_sample(samples_df, sample: str) -> list[str]:
    """Return all FASTQ basenames belonging to a given sample.

    Args:
        samples_df: DataFrame indexed by fastq_files_basename.
        sample: Sample name to look up.

    Returns:
        List of FASTQ basenames for the sample.
    """
    result: list[str] = samples_df.loc[
        samples_df["project_sample"] == sample, "fastq_files_basename"
    ].tolist()
    return result


def resolve_fastq_path(
    basename: str,
    suffix: str,
    samples_df,
    trimming_enabled: bool,
    trimmed_dir: str,
    fastq_dir: str,
) -> str:
    """Build FASTQ path respecting trimming mode and optional subfolder.

    Args:
        basename: FASTQ file basename (used as index into samples_df).
        suffix: File suffix (e.g. '_R1_001.fastq.gz').
        samples_df: DataFrame indexed by fastq_files_basename.
        trimming_enabled: Whether trimming mode is active.
        trimmed_dir: Directory for trimmed FASTQ files.
        fastq_dir: Directory for raw FASTQ files.

    Returns:
        Full path to the FASTQ file.
    """
    row = samples_df.loc[basename]
    subfolder = row.get("subfolder", "") if "subfolder" in samples_df.columns else ""
    if trimming_enabled:
        base_dir = trimmed_dir
    else:
        base_dir = fastq_dir
    if subfolder:
        return os.path.join(base_dir, str(subfolder), f"{basename}{suffix}")
    return os.path.join(base_dir, f"{basename}{suffix}")
