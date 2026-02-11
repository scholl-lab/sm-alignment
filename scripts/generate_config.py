#!/usr/bin/env python3
"""Generate config/samples.tsv from Illumina SampleSheet.csv and FASTQ directories.

Parses sequencing facility deliverables (SampleSheet.csv + FASTQ files) to produce
the samples.tsv metadata file required by the sm-alignment Snakemake pipeline.

Supports two SampleSheet formats:
  1. BIH/Charite minimal (bare CSV rows, no section headers)
  2. Standard Illumina with [Header]/[Data] sections

Usage:
    python scripts/generate_config.py --fastq-dir /path/to/fastqs
    python scripts/generate_config.py --fastq-dir /path/to/fastqs --samplesheet SampleSheet.csv
    python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template --dry-run
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pandas as pd
except ImportError:
    sys.exit(
        "Error: pandas is required but not installed.\n"
        "Install it with: pip install pandas\n"
        "(pandas is already available in Snakemake conda environments.)"
    )


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ILLUMINA_RE = re.compile(
    r"^(?P<sample>.+?)_S(?P<snum>\d+)_L(?P<lane>\d{3})_(?P<read>R[12])_001\.fastq\.gz$"
)

# Pattern to extract project ID from folder/sample names (e.g., A5297)
PROJECT_ID_RE = re.compile(r"(A\d{4,})")

# DNA bases only (for detecting index columns in minimal SampleSheet format)
INDEX_RE = re.compile(r"^[ACGTNacgtn]+$")

# Known section headers in Illumina SampleSheet v2
SECTION_HEADERS = {"[Data]", "[BCLConvert_Data]", "[BCLConvert_Settings]", "[Header]",
                   "[Reads]", "[Settings]"}


# ---------------------------------------------------------------------------
# SampleSheet parsing
# ---------------------------------------------------------------------------

def detect_samplesheet_format(lines: List[str]) -> str:
    """Detect whether a SampleSheet uses Illumina sections or minimal format.

    Args:
        lines: Non-empty, stripped lines from the SampleSheet file.

    Returns:
        'illumina' if section headers like [Data] are found, 'minimal' otherwise.
    """
    for line in lines:
        stripped = line.strip()
        if stripped and stripped.split(",")[0].strip() in SECTION_HEADERS:
            return "illumina"
    return "minimal"


def parse_illumina_samplesheet(lines: List[str]) -> List[Dict[str, str]]:
    """Parse a standard Illumina SampleSheet with [Header]/[Data] sections.

    Skips everything until [Data] or [BCLConvert_Data] section is found,
    then reads the CSV header and data rows.

    Args:
        lines: All lines from the SampleSheet file.

    Returns:
        List of dicts with normalized keys: lane, sample_name, index_i7,
        index_i5, sample_project.
    """
    data_start = None
    for i, line in enumerate(lines):
        stripped = line.strip().split(",")[0].strip()
        if stripped in ("[Data]", "[BCLConvert_Data]"):
            data_start = i + 1
            break

    if data_start is None:
        print("Warning: No [Data] or [BCLConvert_Data] section found in SampleSheet.",
              file=sys.stderr)
        return []

    # Read the header row
    if data_start >= len(lines):
        return []

    header_line = lines[data_start].strip()
    reader = csv.DictReader(lines[data_start:])

    # Normalize column names: Illumina uses varying capitalization
    column_map = {
        "lane": "lane",
        "sample_id": "sample_id",
        "sample_name": "sample_name",
        "index": "index_i7",
        "index2": "index_i5",
        "sample_project": "sample_project",
    }

    samples = []
    for row in reader:
        if not any(v.strip() for v in row.values() if v):
            continue  # skip empty rows

        normalized = {}
        for key, value in row.items():
            if key is None:
                continue
            norm_key = key.strip().lower().replace(" ", "_")
            if norm_key in column_map:
                normalized[column_map[norm_key]] = value.strip() if value else ""

        # Use Sample_Name if available, fall back to Sample_ID
        sample_name = normalized.get("sample_name", "") or normalized.get("sample_id", "")

        samples.append({
            "lane": normalized.get("lane", ""),
            "sample_name": sample_name,
            "index_i7": normalized.get("index_i7", ""),
            "index_i5": normalized.get("index_i5", ""),
            "sample_project": normalized.get("sample_project", ""),
        })

    return samples


def parse_minimal_samplesheet(lines: List[str]) -> List[Dict[str, str]]:
    """Parse a BIH/Charite minimal SampleSheet (bare CSV, no section headers).

    Expected columns (positional): Lane, Sample_Name, index_i7, index_i5, Sample_Project

    Heuristic detection: Lane is numeric, indices are ACGT sequences.

    Args:
        lines: All lines from the SampleSheet file.

    Returns:
        List of dicts with keys: lane, sample_name, index_i7, index_i5, sample_project.
    """
    samples = []

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = [f.strip() for f in stripped.split(",")]
        if len(fields) < 4:
            continue

        # Check if first line is a header row
        if fields[0].lower() in ("lane", "sample_lane"):
            continue

        # Heuristic: detect column layout
        # Look for: numeric lane, sample name, two index columns, optional project
        lane_val = None
        sample_name = None
        index_i7 = None
        index_i5 = None
        project = ""

        if len(fields) >= 5:
            # Standard 5-column: Lane, Sample_Name, index_i7, index_i5, Sample_Project
            if fields[0].isdigit() and INDEX_RE.match(fields[2]) and INDEX_RE.match(fields[3]):
                lane_val = fields[0]
                sample_name = fields[1]
                index_i7 = fields[2]
                index_i5 = fields[3]
                project = fields[4]
            # Alternative: Sample_Name, Lane, index_i7, index_i5, Sample_Project
            elif fields[1].isdigit() and INDEX_RE.match(fields[2]) and INDEX_RE.match(fields[3]):
                sample_name = fields[0]
                lane_val = fields[1]
                index_i7 = fields[2]
                index_i5 = fields[3]
                project = fields[4]
        elif len(fields) == 4:
            # 4-column: Lane, Sample_Name, index_i7, index_i5 (no project)
            if fields[0].isdigit() and INDEX_RE.match(fields[2]) and INDEX_RE.match(fields[3]):
                lane_val = fields[0]
                sample_name = fields[1]
                index_i7 = fields[2]
                index_i5 = fields[3]

        if lane_val is not None and sample_name is not None:
            samples.append({
                "lane": lane_val,
                "sample_name": sample_name,
                "index_i7": index_i7 or "",
                "index_i5": index_i5 or "",
                "sample_project": project,
            })
        else:
            print(f"Warning: Could not parse SampleSheet line: {stripped}", file=sys.stderr)

    return samples


def parse_samplesheet(path: Path) -> List[Dict[str, str]]:
    """Parse a SampleSheet.csv file, auto-detecting the format.

    Args:
        path: Path to the SampleSheet.csv file.

    Returns:
        List of dicts with keys: lane, sample_name, index_i7, index_i5, sample_project.
    """
    with open(path, "r", encoding="utf-8-sig") as fh:
        lines = fh.readlines()

    if not lines:
        print(f"Warning: SampleSheet is empty: {path}", file=sys.stderr)
        return []

    fmt = detect_samplesheet_format(lines)

    if fmt == "illumina":
        samples = parse_illumina_samplesheet(lines)
    else:
        samples = parse_minimal_samplesheet(lines)

    return samples


def find_samplesheet(fastq_dir: Path) -> Optional[Path]:
    """Auto-detect a SampleSheet.csv in the FASTQ directory or its parent.

    Searches for common SampleSheet file names in the FASTQ directory
    and one level up.

    Args:
        fastq_dir: Directory containing FASTQ files.

    Returns:
        Path to the SampleSheet if found, None otherwise.
    """
    candidates = [
        "SampleSheet.csv",
        "samplesheet.csv",
        "SampleSheet_v2.csv",
    ]
    search_dirs = [fastq_dir]
    if fastq_dir.parent != fastq_dir:
        search_dirs.append(fastq_dir.parent)

    for search_dir in search_dirs:
        for name in candidates:
            path = search_dir / name
            if path.is_file():
                return path
    return None


# ---------------------------------------------------------------------------
# FASTQ discovery
# ---------------------------------------------------------------------------

def discover_fastq_files(fastq_dir: Path) -> List[Dict[str, str]]:
    """Scan a directory for Illumina-named FASTQ file pairs.

    Finds all *_R1_001.fastq.gz files matching the Illumina naming convention,
    extracts metadata, and verifies that a matching R2 file exists.

    Args:
        fastq_dir: Directory to scan (non-recursive).

    Returns:
        List of dicts with keys: basename, lane, lane_str, sample_name, r1_path, r2_path.
    """
    if not fastq_dir.is_dir():
        sys.exit(f"Error: FASTQ directory does not exist: {fastq_dir}")

    all_files = sorted(fastq_dir.iterdir())
    fastq_gz_files = [f for f in all_files if f.name.endswith(".fastq.gz")]

    if not fastq_gz_files:
        sys.exit(
            f"Error: No .fastq.gz files found in {fastq_dir}\n"
            "Please check the directory path and ensure FASTQ files are present."
        )

    pairs = []
    r1_count = 0
    skipped_r2 = 0

    for fpath in fastq_gz_files:
        match = ILLUMINA_RE.match(fpath.name)
        if not match:
            continue
        if match.group("read") != "R1":
            continue

        r1_count += 1
        sample = match.group("sample")
        snum = match.group("snum")
        lane = match.group("lane")

        basename = f"{sample}_S{snum}_L{lane}"
        lane_str = f"L{lane}"
        r2_name = fpath.name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        r2_path = fpath.parent / r2_name

        if not r2_path.is_file():
            print(f"Warning: R2 missing for {fpath.name}, skipping pair.", file=sys.stderr)
            skipped_r2 += 1
            continue

        pairs.append({
            "basename": basename,
            "lane": lane,
            "lane_str": lane_str,
            "sample_name": sample,
            "r1_path": str(fpath),
            "r2_path": str(r2_path),
        })

    if not pairs:
        sys.exit(
            f"Error: Found {r1_count} R1 files but no valid pairs in {fastq_dir}\n"
            f"(Skipped {skipped_r2} files with missing R2 mates.)\n"
            "Expected Illumina naming: <sample>_S<N>_L<NNN>_R[12]_001.fastq.gz"
        )

    return pairs


# ---------------------------------------------------------------------------
# Merge SampleSheet + FASTQ data into samples.tsv rows
# ---------------------------------------------------------------------------

def infer_project(fastq_dir: Path, project_arg: Optional[str]) -> str:
    """Determine the project identifier.

    Priority:
        1. Explicit --project CLI argument
        2. Extract A{NNNN} pattern from directory name
        3. Use directory name as-is

    Args:
        fastq_dir: The FASTQ directory path.
        project_arg: Value of --project CLI argument, or None.

    Returns:
        Project identifier string.
    """
    if project_arg:
        return project_arg

    # Try to extract A{NNNN+} pattern from the directory path
    dir_name = fastq_dir.resolve().name
    match = PROJECT_ID_RE.search(str(fastq_dir.resolve()))
    if match:
        return match.group(1)

    return dir_name


def build_samples_table(
    fastq_pairs: List[Dict[str, str]],
    samplesheet_entries: Optional[List[Dict[str, str]]],
    project: str,
) -> pd.DataFrame:
    """Combine FASTQ discovery results with SampleSheet metadata.

    For each FASTQ pair, attempts to match it with a SampleSheet entry
    by sample name. Falls back to FASTQ-derived metadata if no SampleSheet
    or no match found.

    Args:
        fastq_pairs: Results from discover_fastq_files().
        samplesheet_entries: Results from parse_samplesheet(), or None.
        project: Project identifier for the mdc_project column.

    Returns:
        DataFrame with columns: fastq_files_basename, lane, project_sample, mdc_project.
    """
    # Build a lookup from SampleSheet sample_name to entry
    ss_lookup: Dict[str, Dict[str, str]] = {}
    if samplesheet_entries:
        for entry in samplesheet_entries:
            name = entry["sample_name"]
            if name:
                ss_lookup[name] = entry

    rows = []
    matched_ss = set()

    for pair in fastq_pairs:
        fastq_sample = pair["sample_name"]
        lane_str = pair["lane_str"]
        basename = pair["basename"]

        # Try to match SampleSheet entry by sample name
        # The FASTQ sample name might be a superset of the SampleSheet name
        # (e.g., SampleSheet has "A5297_DNA_01_STREAM_P1_L1" and FASTQ has the same)
        project_sample = fastq_sample  # default: derive from FASTQ
        entry_project = project

        if fastq_sample in ss_lookup:
            ss_entry = ss_lookup[fastq_sample]
            project_sample = ss_entry["sample_name"]
            if ss_entry.get("sample_project"):
                entry_project = ss_entry["sample_project"] or project
            matched_ss.add(fastq_sample)
        else:
            # Try partial matching: SampleSheet name contained in FASTQ name
            for ss_name, ss_entry in ss_lookup.items():
                if ss_name and (ss_name in fastq_sample or fastq_sample.startswith(ss_name)):
                    project_sample = ss_entry["sample_name"]
                    if ss_entry.get("sample_project"):
                        entry_project = ss_entry["sample_project"] or project
                    matched_ss.add(ss_name)
                    break

        rows.append({
            "fastq_files_basename": basename,
            "lane": lane_str,
            "project_sample": project_sample,
            "mdc_project": entry_project,
        })

    # Warn about SampleSheet entries that had no matching FASTQ files
    if samplesheet_entries:
        for entry in samplesheet_entries:
            name = entry["sample_name"]
            if name and name not in matched_ss:
                print(
                    f"Warning: SampleSheet sample '{name}' has no matching FASTQ files, skipping.",
                    file=sys.stderr,
                )

    df = pd.DataFrame(rows, columns=["fastq_files_basename", "lane", "project_sample", "mdc_project"])
    df = df.sort_values(["project_sample", "lane", "fastq_files_basename"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Config template generation
# ---------------------------------------------------------------------------

CONFIG_TEMPLATE = """\
# =============================================================================
# config/config.yaml -- sm-alignment pipeline configuration
# =============================================================================
# Generated by: scripts/generate_config.py
# Edit paths and parameters below for your project.

# --- Reference genome & known variant sites ---
ref:
  genome: "analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  genome_gz: "analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
  build: "GRCh38"
  known_sites:
    - "analysis/GATK_resource_bundle/af-only-gnomad.hg38.vcf.gz"
    - "analysis/GATK_resource_bundle/af-only-gnomad.hg38.common_biallelic.vcf.gz"
    - "analysis/GATK_resource_bundle/1000g_pon.hg38.vcf.gz"

# --- Paths ---
paths:
  samples: "{samples_path}"
  fastq_folder: "{fastq_folder}"
  output_folder: "{output_folder}"
  log_subdir: "logs"

# --- Read group defaults ---
read_group:
  platform: "ILLUMINA"

# --- FASTQ file naming patterns ---
fastq:
  r1_suffix: "_R1_001.fastq.gz"
  r2_suffix: "_R2_001.fastq.gz"
  trimmed_r1_suffix: ".bbduk_R1_001.fastq.gz"
  trimmed_r2_suffix: ".bbduk_R2_001.fastq.gz"

# --- BAM file naming suffixes ---
bam:
  merged_suffix: ".merged.bam"
  dedup_suffix: ".merged.dedup.bam"
  dedup_metrics_suffix: ".merged.dedup_metrics.txt"
  recal_table_suffix: ".merged.dedup.recal_data.table"
  final_suffix: ".merged.dedup.bqsr.bam"

# --- Processing options ---
processing:
  compression_level: 6

# --- Tool-specific extra CLI arguments (passthrough strings) ---
params:
  bwa_mem:
    extra: ""
  samtools:
    sort_extra: ""
    merge_extra: ""
  gatk:
    MarkDuplicates: "--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"
    BaseRecalibrator: ""
    ApplyBQSR: ""

# --- BBDuk trimming parameters ---
trimming:
  enabled: false
  bbduk_ref: "adapters,artifacts"
  ktrim: "r"
  k: 23
  mink: 11
  hdist: 1
  tpe: "t"
  tbo: "t"
  ftl: 5
  trimpolyg: 3
  trimpolya: 3
  qtrim: "t"
  trimq: 10
  ziplevel: 5
  quantize: "0,10,20,30,40,50,60"

# --- Subset BAM (optional) ---
subset:
  bed_file: ""
  output_suffix: ".subset.bam"
"""


def generate_config_template(
    config_output: Path,
    fastq_dir: Path,
    samples_path: str,
    project: str,
    dry_run: bool = False,
    force: bool = False,
) -> None:
    """Generate a starter config.yaml from template.

    Args:
        config_output: Where to write the config file.
        fastq_dir: FASTQ directory (fills paths.fastq_folder).
        samples_path: Value for paths.samples.
        project: Project identifier (used in output folder path).
        dry_run: If True, print content but do not write.
        force: If True, overwrite existing file without asking.
    """
    fastq_folder = str(fastq_dir).replace("\\", "/")
    output_folder = f"results/{project}"

    content = CONFIG_TEMPLATE.format(
        samples_path=samples_path,
        fastq_folder=fastq_folder,
        output_folder=output_folder,
    )

    if dry_run:
        print(f"\n--- Config template ({config_output}) ---")
        print(content)
        return

    if config_output.is_file() and not force:
        response = input(f"Config file already exists: {config_output}. Overwrite? [y/N] ")
        if response.lower() not in ("y", "yes"):
            print(f"Skipped writing {config_output}")
            return

    config_output.parent.mkdir(parents=True, exist_ok=True)
    with open(config_output, "w", encoding="utf-8") as fh:
        fh.write(content)

    print(f"Written: {config_output}")


# ---------------------------------------------------------------------------
# Output formatting and writing
# ---------------------------------------------------------------------------

def print_summary(
    df: pd.DataFrame,
    samplesheet_entries: Optional[List[Dict[str, str]]],
    fastq_pairs: List[Dict[str, str]],
) -> None:
    """Print a human-readable summary of the generated samples table.

    Args:
        df: The samples DataFrame.
        samplesheet_entries: Parsed SampleSheet entries (or None).
        fastq_pairs: Discovered FASTQ pairs.
    """
    print()
    if samplesheet_entries:
        print(f"Found SampleSheet.csv with {len(samplesheet_entries)} sample(s)")
    else:
        print("No SampleSheet found; using FASTQ filenames only")
    print(f"Matched {len(fastq_pairs) * 2} FASTQ files ({len(fastq_pairs)} pair(s))")
    print()

    # Print table with aligned columns
    col_widths = {}
    for col in df.columns:
        max_val = df[col].astype(str).str.len().max()
        col_widths[col] = max(len(col), max_val)

    header = "  ".join(col.ljust(col_widths[col]) for col in df.columns)
    print(header)
    for _, row in df.iterrows():
        line = "  ".join(str(row[col]).ljust(col_widths[col]) for col in df.columns)
        print(line)
    print()


def write_samples_tsv(
    df: pd.DataFrame,
    output_path: Path,
    dry_run: bool = False,
    force: bool = False,
) -> None:
    """Write the samples DataFrame to a TSV file.

    Args:
        df: The samples DataFrame.
        output_path: Path for the output TSV.
        dry_run: If True, do not write the file.
        force: If True, overwrite without asking.
    """
    if dry_run:
        print(f"[dry-run] Would write: {output_path} ({len(df)} sample(s))")
        return

    if output_path.is_file() and not force:
        response = input(f"Output file already exists: {output_path}. Overwrite? [y/N] ")
        if response.lower() not in ("y", "yes"):
            print(f"Skipped writing {output_path}")
            return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"Written: {output_path} ({len(df)} sample(s))")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point: parse arguments, discover files, generate output."""
    parser = argparse.ArgumentParser(
        description="Generate config/samples.tsv from sequencing facility deliverables.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s --fastq-dir /data/project/fastqs
  %(prog)s --fastq-dir /data/project/fastqs --samplesheet SampleSheet.csv
  %(prog)s --fastq-dir /data/project/fastqs --config-template --dry-run
  %(prog)s --fastq-dir /data/project/fastqs --project A5297 --force
        """,
    )
    parser.add_argument(
        "--fastq-dir",
        required=True,
        help="Directory containing FASTQ files (and optionally SampleSheet.csv)",
    )
    parser.add_argument(
        "--samplesheet",
        help="Path to SampleSheet.csv (auto-detected in fastq-dir if not specified)",
    )
    parser.add_argument(
        "--project",
        help="Project identifier for mdc_project column (default: inferred from folder name)",
    )
    parser.add_argument(
        "--output",
        default="config/samples.tsv",
        help="Output samples.tsv path (default: config/samples.tsv)",
    )
    parser.add_argument(
        "--config-template",
        action="store_true",
        help="Also generate a config.yaml template",
    )
    parser.add_argument(
        "--config-output",
        default="config/config.yaml",
        help="Output config.yaml path (default: config/config.yaml)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be generated without writing files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files without asking",
    )

    args = parser.parse_args()

    fastq_dir = Path(args.fastq_dir).resolve()

    # ---- Discover FASTQ files ----
    fastq_pairs = discover_fastq_files(fastq_dir)

    # ---- Parse SampleSheet (if available) ----
    samplesheet_entries: Optional[List[Dict[str, str]]] = None

    if args.samplesheet:
        ss_path = Path(args.samplesheet)
        if not ss_path.is_file():
            sys.exit(f"Error: SampleSheet not found: {ss_path}")
        samplesheet_entries = parse_samplesheet(ss_path)
    else:
        ss_path = find_samplesheet(fastq_dir)
        if ss_path:
            samplesheet_entries = parse_samplesheet(ss_path)

    # ---- Determine project ----
    project = infer_project(fastq_dir, args.project)

    # ---- Build samples table ----
    df = build_samples_table(fastq_pairs, samplesheet_entries, project)

    # ---- Print summary ----
    print_summary(df, samplesheet_entries, fastq_pairs)

    # ---- Write samples.tsv ----
    output_path = Path(args.output)
    write_samples_tsv(df, output_path, dry_run=args.dry_run, force=args.force)

    # ---- Generate config template (optional) ----
    if args.config_template:
        config_output = Path(args.config_output)
        generate_config_template(
            config_output=config_output,
            fastq_dir=fastq_dir,
            samples_path=args.output,
            project=project,
            dry_run=args.dry_run,
            force=args.force,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
