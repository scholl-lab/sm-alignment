"""Pure-Python helper functions for the sm-alignment pipeline.

Extracted from common.smk so they can be unit-tested independently.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

FGBIO_READ_STRUCTURES_RE = re.compile(
    r"^(?=.*M)(?=[^\s]*T)(?:[1-9][0-9]*[TBMCS])*(?:[1-9][0-9]*|\+)[TBMCS]"
    r"\s+(?=[^\s]*T)(?:[1-9][0-9]*[TBMCS])*(?:[1-9][0-9]*|\+)[TBMCS]$"
)

DEFAULT_UMI_PRESETS_PATH = (
    Path(__file__).resolve().parent.parent.parent / "config" / "umi_presets.yaml"
)


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
    umi_enabled: bool = False,
    umi_dir: str = "",
) -> str:
    """Build FASTQ path respecting trimming/UMI mode and optional subfolder.

    Args:
        basename: FASTQ file basename (used as index into samples_df).
        suffix: File suffix (e.g. '_R1_001.fastq.gz').
        samples_df: DataFrame indexed by fastq_files_basename.
        trimming_enabled: Whether trimming mode is active.
        trimmed_dir: Directory for trimmed FASTQ files.
        fastq_dir: Directory for raw FASTQ files.
        umi_enabled: Whether UMI extraction is active.
        umi_dir: Directory for UMI-extracted FASTQ files.

    Returns:
        Full path to the FASTQ file.
    """
    if trimming_enabled:
        return os.path.join(trimmed_dir, f"{basename}{suffix}")
    if umi_enabled:
        return os.path.join(umi_dir, f"{basename}{suffix}")
    row = samples_df.loc[basename]
    subfolder = row.get("subfolder", "") if "subfolder" in samples_df.columns else ""
    if subfolder and str(subfolder).strip():
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
    if "kit" in df.columns and "project_sample" in df.columns:
        for sample, grp in df.groupby("project_sample"):
            kits = {str(k).strip() for k in grp["kit"]}
            if len(kits) > 1:
                raise ValueError(
                    f"Sample '{sample}' has conflicting kits in {path}: {sorted(kits)}"
                )
    return df.set_index("fastq_files_basename", drop=False)


def get_sample_kit(samples_df, sample: str) -> str:
    """Return the capture kit identifier for a sample (or '' if unset)."""
    if "kit" not in samples_df.columns:
        return ""
    vals = {
        str(k).strip()
        for k in samples_df.loc[samples_df["project_sample"] == sample, "kit"]
        if str(k).strip()
    }
    return next(iter(vals)) if vals else ""


def get_read_group_fields(
    samples_df,
    basename: str,
    platform: str = "ILLUMINA",
) -> dict[str, str]:
    """Extract read group fields (id, sample, library, platform, pu) for a basename."""
    row = samples_df.loc[basename]
    sample = str(row["project_sample"]).strip()
    lane = str(row["lane"]).strip()

    flowcell = ""
    for col in ("flowcell", "subfolder", "mdc_project"):
        if col in samples_df.columns:
            val = row.get(col, "")
            if val is not None and str(val).strip():
                flowcell = str(val).strip()
                break

    library = sample
    if "library" in samples_df.columns:
        val = row.get("library", "")
        if val is not None and str(val).strip():
            library = str(val).strip()

    rg_id = f"{flowcell}.{lane}.{library}"
    return {
        "id": rg_id,
        "sample": sample,
        "library": library,
        "platform": platform,
        "pu": rg_id,
    }


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
    rg = get_read_group_fields(samples_df, basename, platform)
    return (
        f"@RG\\tID:{rg['id']}\\tSM:{rg['sample']}\\tLB:{rg['library']}"
        f"\\tPL:{rg['platform']}\\tPU:{rg['pu']}"
    )


def load_umi_presets(path: str | os.PathLike[str] | None = None) -> dict[str, dict[str, Any]]:
    """Load verified UMI presets from config/umi_presets.yaml."""
    preset_path = Path(path) if path is not None else DEFAULT_UMI_PRESETS_PATH
    with open(preset_path) as f:
        data = yaml.safe_load(f) or {}
    return data


def _validate_read_structures_syntax(read_structures: str) -> str:
    """Validate paired-end fgbio read-structure syntax and return normalized string."""
    rs = " ".join(str(read_structures).split())
    if not FGBIO_READ_STRUCTURES_RE.match(rs):
        raise ValueError(
            f"Invalid fgbio read_structures '{read_structures}': expected two space-separated "
            "read structures with <int|finally '+'>[TBMCS] segments including M and T "
            "(e.g. '3M2S+T 3M2S+T')."
        )
    return rs


def validate_umi_config(
    umi_cfg: dict[str, Any] | None,
    presets: dict[str, dict[str, Any]] | None = None,
) -> None:
    """Validate the umi configuration block with clear error messages."""
    if not umi_cfg:
        return
    known_presets = presets if presets is not None else load_umi_presets()
    enabled = bool(umi_cfg.get("enabled", False))
    raw_preset = umi_cfg.get("preset")
    preset = str(raw_preset).strip() if raw_preset is not None else ""
    raw_rs = umi_cfg.get("read_structures")
    rs = str(raw_rs).strip() if raw_rs is not None else ""

    if preset and preset != "custom" and preset not in known_presets:
        raise ValueError(
            f"Unknown umi.preset '{preset}'. Valid presets: {sorted(known_presets.keys())} or 'custom'."
        )
    if rs:
        _validate_read_structures_syntax(rs)

    umi_tag = str(umi_cfg.get("umi_tag") or "RX").strip()
    dedup_cfg = umi_cfg.get("dedup") or {}
    raw_bt = dedup_cfg.get("barcode_tag")
    if raw_bt is not None and str(raw_bt).strip() and str(raw_bt).strip() != umi_tag:
        raise ValueError(
            f"umi.dedup.barcode_tag '{str(raw_bt).strip()}' must match umi.umi_tag '{umi_tag}' "
            "(or be null/empty for position-only deduplication)."
        )

    if enabled:
        has_kit_preset = bool(preset and preset != "custom")
        has_rs = bool(rs)
        if has_kit_preset and has_rs:
            raise ValueError(
                "umi config conflict: specify either 'umi.preset' or 'umi.read_structures', not both."
            )
        if not has_kit_preset and not has_rs:
            raise ValueError(
                "umi.enabled is true: must specify either 'umi.preset' or 'umi.read_structures'."
            )


def resolve_umi_read_structures(
    umi_cfg: dict[str, Any] | None,
    presets: dict[str, dict[str, Any]] | None = None,
) -> str:
    """Resolve effective fgbio read_structures from umi config and presets.

    A custom read_structures value overrides preset when provided.
    Returns '' when umi.enabled is false.
    """
    if not umi_cfg or not umi_cfg.get("enabled", False):
        return ""
    known_presets = presets if presets is not None else load_umi_presets()
    raw_rs = umi_cfg.get("read_structures")
    rs = str(raw_rs).strip() if raw_rs is not None else ""
    if rs:
        return _validate_read_structures_syntax(rs)

    raw_preset = umi_cfg.get("preset")
    preset = str(raw_preset).strip() if raw_preset is not None else ""
    if not preset or preset == "custom":
        raise ValueError(
            "umi.enabled is true: must specify either 'umi.preset' or 'umi.read_structures'."
        )
    if preset not in known_presets:
        raise ValueError(
            f"Unknown umi.preset '{preset}'. Valid presets: {sorted(known_presets.keys())} or 'custom'."
        )
    return _validate_read_structures_syntax(str(known_presets[preset]["read_structures"]))


def write_umi_provenance(
    log_dir: str | os.PathLike[str],
    umi_cfg: dict[str, Any] | None,
    presets: dict[str, dict[str, Any]] | None = None,
) -> Path | None:
    """Log the resolved UMI read structure and write umi_provenance.yaml."""
    if not umi_cfg or not umi_cfg.get("enabled", False):
        return None
    known_presets = presets if presets is not None else load_umi_presets()
    read_structures = resolve_umi_read_structures(umi_cfg, known_presets)
    raw_preset = umi_cfg.get("preset")
    preset = str(raw_preset).strip() if raw_preset is not None and str(raw_preset).strip() else ""
    raw_rs = umi_cfg.get("read_structures")
    has_custom_rs = bool(raw_rs is not None and str(raw_rs).strip())

    effective_preset = "custom" if has_custom_rs else (preset or "custom")
    source = (
        str(known_presets[preset].get("source", ""))
        if (not has_custom_rs and preset in known_presets)
        else "custom (config.umi.read_structures)"
    )
    umi_tag = str(umi_cfg.get("umi_tag") or "RX")
    duplex = bool(umi_cfg.get("duplex", True))
    dedup_cfg = umi_cfg.get("dedup") or {}
    barcode_tag = dedup_cfg.get("barcode_tag", umi_tag)

    provenance = {
        "umi": {
            "enabled": True,
            "preset": effective_preset,
            "read_structures": read_structures,
            "umi_tag": umi_tag,
            "duplex": duplex,
            "barcode_tag": barcode_tag,
            "source": source,
        }
    }
    logger.info(
        "UMI enabled: preset=%s read_structures='%s' umi_tag=%s duplex=%s source=%s",
        effective_preset,
        read_structures,
        umi_tag,
        duplex,
        source,
    )
    out_dir = Path(log_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "umi_provenance.yaml"
    out_path.write_text(yaml.safe_dump(provenance, sort_keys=False))
    return out_path


def format_bwa_extra(bwa_extra: str, umi_enabled: bool = False) -> str:
    """Return bwa mem extra arguments, appending -C when UMI mode is enabled."""
    base = bwa_extra.strip()
    if not umi_enabled:
        return base
    tokens = base.split()
    if "-C" not in tokens:
        return f"{base} -C".strip()
    return base


def format_mark_duplicates_extra(
    base_extra: str,
    umi_cfg: dict[str, Any] | None = None,
) -> str:
    """Return GATK MarkDuplicates extra arguments with optional UMI flags."""
    cleaned = re.sub(r"--CREATE_INDEX(?:\s+|=)\S+", "", base_extra).strip()
    if not umi_cfg or not umi_cfg.get("enabled", False):
        return cleaned

    umi_tag = str(umi_cfg.get("umi_tag") or "RX")
    dedup_cfg = umi_cfg.get("dedup") or {}
    raw_barcode_tag = dedup_cfg.get("barcode_tag", umi_tag)
    barcode_tag = str(raw_barcode_tag).strip() if raw_barcode_tag is not None else ""

    parts = [cleaned] if cleaned else []
    if barcode_tag:
        parts.append(f"--BARCODE_TAG {barcode_tag}")
        if umi_cfg.get("duplex", True):
            parts.append("--DUPLEX_UMI true")

    dedup_extra = str(dedup_cfg.get("extra") or "").strip()
    if dedup_extra:
        parts.append(dedup_extra)

    return " ".join(parts)
