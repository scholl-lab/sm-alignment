"""Pure-Python helper functions for the sm-alignment pipeline.

Extracted from common.smk so they can be unit-tested independently.
"""

from __future__ import annotations

import os

import pandas as pd


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


def get_all_basenames(samples_df) -> list[str]:
    """Return sorted list of all FASTQ basenames from metadata.

    Args:
        samples_df: DataFrame with a 'fastq_files_basename' column.

    Returns:
        Sorted list of all FASTQ basenames.
    """
    return sorted(samples_df["fastq_files_basename"].unique().tolist())


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
    if trimming_enabled:
        return os.path.join(trimmed_dir, f"{basename}{suffix}")
    row = samples_df.loc[basename]
    subfolder = row.get("subfolder", "") if "subfolder" in samples_df.columns else ""
    if subfolder and str(subfolder).strip() and str(subfolder).strip().lower() != "nan":
        return os.path.join(fastq_dir, str(subfolder).strip(), f"{basename}{suffix}")
    return os.path.join(fastq_dir, f"{basename}{suffix}")


def load_samples(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Load and validate samples.tsv with string columns and no NaN coercion."""
    df = pd.read_table(path, dtype=str, keep_default_na=False)
    if df.empty:
        raise ValueError(f"Samples table is empty: {path}")
    if "fastq_files_basename" not in df.columns:
        raise ValueError(f"Missing required column 'fastq_files_basename' in {path}")
    basenames = df["fastq_files_basename"]
    if (basenames.str.strip() == "").any():
        raise ValueError(f"Empty fastq_files_basename entry in {path}")
    if basenames.str.contains(r"[/\\]", regex=True).any():
        raise ValueError(f"fastq_files_basename must not contain path separators in {path}")
    if basenames.duplicated().any():
        dups = sorted(basenames[basenames.duplicated()].unique().tolist())
        raise ValueError(f"Duplicate fastq_files_basename entries in {path}: {dups}")
    return df.set_index("fastq_files_basename", drop=False)


def format_read_group(
    samples_df,
    basename: str,
    platform: str = "ILLUMINA",
) -> str:
    """Format @RG header string for a FASTQ basename row in samples_df.

    Produces ID:{flowcell}.{lane}.{library}, SM:{sample}, LB:{library},
    PL:{platform}, PU:{flowcell}.{lane}.{library}, falling back to
    subfolder/mdc_project and project_sample when flowcell/library are omitted.
    """
    row = samples_df.loc[basename]
    sample = str(row["project_sample"]).strip()
    lane = str(row["lane"]).strip()

    flowcell = ""
    for col in ("flowcell", "subfolder", "mdc_project"):
        if col in samples_df.columns:
            val = row.get(col, "")
            if val is not None and str(val).strip() and str(val).strip().lower() != "nan":
                flowcell = str(val).strip()
                break

    library = sample
    if "library" in samples_df.columns:
        val = row.get("library", "")
        if val is not None and str(val).strip() and str(val).strip().lower() != "nan":
            library = str(val).strip()

    rg_id = f"{flowcell}.{lane}.{library}"
    return f"@RG\\tID:{rg_id}\\tSM:{sample}\\tLB:{library}\\tPL:{platform}\\tPU:{rg_id}"
