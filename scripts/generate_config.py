#!/usr/bin/env python3
"""Generate config/samples.tsv from Illumina SampleSheet.csv and FASTQ directories.

Parses sequencing facility deliverables (SampleSheet.csv + FASTQ files) to produce
the samples.tsv metadata file required by the sm-alignment Snakemake pipeline.

Supports two SampleSheet formats:
  1. BIH/Charite minimal (bare CSV rows, no section headers)
  2. Standard Illumina with [Header]/[Data] sections

Two modes:
  Interactive:  python scripts/generate_config.py          (guided wizard)
  Flags:        python scripts/generate_config.py --fastq-dir /path/to/fastqs

Usage:
    python scripts/generate_config.py
    python scripts/generate_config.py --fastq-dir /path/to/fastqs
    python scripts/generate_config.py --fastq-dir /path/to/fastqs --samplesheet SampleSheet.csv
    python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template --dry-run
    python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template --output-dir /data/results
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any

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
SECTION_HEADERS = {
    "[Data]",
    "[BCLConvert_Data]",
    "[BCLConvert_Settings]",
    "[Header]",
    "[Reads]",
    "[Settings]",
}

# Reference genome file extensions
GENOME_EXTENSIONS = ("*.fna", "*.fa", "*.fasta")

# BWA index companion extensions (classic and bwtsw/64-bit)
BWA_INDEX_EXTS = (
    ".amb",
    ".ann",
    ".bwt",
    ".pac",
    ".sa",
    ".64.amb",
    ".64.ann",
    ".64.bwt",
    ".64.pac",
    ".64.sa",
)

# BQSR-appropriate known-sites VCF patterns (GATK BaseRecalibrator)
BQSR_KNOWN_SITES_PATTERNS = (
    "*dbsnp*.vcf*",
    "*known_indels*.vcf*",
    "*Mills_and_1000G*.vcf*",
    "*1000G_phase1.snps.high_confidence*.vcf*",
)

# Additional VCFs found but NOT used for BQSR (retained for reference):
# gnomAD af-only → Mutect2 --germline-resource
# 1000g_pon      → Mutect2 --panel-of-normals

# Legacy alias kept for backward compatibility
KNOWN_SITES_PATTERNS = BQSR_KNOWN_SITES_PATTERNS

# Standard search directories for reference data (relative to project root)
REF_SEARCH_DIRS = [
    "resources/ref",
    "resources/ref/GRCh38",
    "analysis/ref/GRCh38",
    "../resources/ref/GRCh38",
    "../resources/ref",
]

# Known-sites search directories
KNOWN_SITES_SEARCH_DIRS = [
    "resources/gatk_bundle/hg38",
    "resources/gatk_bundle",
    "analysis/GATK_resource_bundle",
    "../resources/gatk_bundle/hg38",
    "../resources/gatk_bundle",
]

# Well-known shared locations on BIH and Charité HPC
SHARED_REF_DIRS = [
    "/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38",
    "/data/cephfs-1/work/groups/scholl/shared/ref",
]

SHARED_KNOWN_SITES_DIRS = [
    "/data/cephfs-1/work/projects/apa-sequencing/analysis/GATK_resource_bundle",
]


# ---------------------------------------------------------------------------
# SampleSheet parsing
# ---------------------------------------------------------------------------


def detect_samplesheet_format(lines: list[str]) -> str:
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


def parse_illumina_samplesheet(lines: list[str]) -> list[dict[str, str]]:
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
        print(
            "Warning: No [Data] or [BCLConvert_Data] section found in SampleSheet.",
            file=sys.stderr,
        )
        return []

    # Read the header row
    if data_start >= len(lines):
        return []

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

        samples.append(
            {
                "lane": normalized.get("lane", ""),
                "sample_name": sample_name,
                "index_i7": normalized.get("index_i7", ""),
                "index_i5": normalized.get("index_i5", ""),
                "sample_project": normalized.get("sample_project", ""),
            }
        )

    return samples


def parse_minimal_samplesheet(lines: list[str]) -> list[dict[str, str]]:
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
            samples.append(
                {
                    "lane": lane_val,
                    "sample_name": sample_name,
                    "index_i7": index_i7 or "",
                    "index_i5": index_i5 or "",
                    "sample_project": project,
                }
            )
        else:
            print(
                f"Warning: Could not parse SampleSheet line: {stripped}",
                file=sys.stderr,
            )

    return samples


def parse_samplesheet(path: Path) -> list[dict[str, str]]:
    """Parse a SampleSheet.csv file, auto-detecting the format.

    Args:
        path: Path to the SampleSheet.csv file.

    Returns:
        List of dicts with keys: lane, sample_name, index_i7, index_i5, sample_project.
    """
    with open(path, encoding="utf-8-sig") as fh:
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


def find_samplesheet(fastq_dir: Path) -> Path | None:
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


def discover_fastq_files(fastq_dir: Path) -> list[dict[str, str]]:
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

        pairs.append(
            {
                "basename": basename,
                "lane": lane,
                "lane_str": lane_str,
                "sample_name": sample,
                "r1_path": str(fpath),
                "r2_path": str(r2_path),
            }
        )

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


def infer_project(fastq_dir: Path, project_arg: str | None) -> str:
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
    fastq_pairs: list[dict[str, str]],
    samplesheet_entries: list[dict[str, str]] | None,
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
    ss_lookup: dict[str, dict[str, str]] = {}
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

        rows.append(
            {
                "fastq_files_basename": basename,
                "lane": lane_str,
                "project_sample": project_sample,
                "mdc_project": entry_project,
            }
        )

    # Warn about SampleSheet entries that had no matching FASTQ files
    if samplesheet_entries:
        for entry in samplesheet_entries:
            name = entry["sample_name"]
            if name and name not in matched_ss:
                print(
                    f"Warning: SampleSheet sample '{name}' has no matching FASTQ files, skipping.",
                    file=sys.stderr,
                )

    df = pd.DataFrame(
        rows, columns=["fastq_files_basename", "lane", "project_sample", "mdc_project"]
    )
    df = df.sort_values(["project_sample", "lane", "fastq_files_basename"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Reference data discovery
# ---------------------------------------------------------------------------


def _find_genome_fastas(search_dir: Path) -> list[dict[str, Any]]:
    """Find reference genome FASTA files in a directory."""
    results: list[dict[str, Any]] = []
    if not search_dir.is_dir():
        return results
    for ext in GENOME_EXTENSIONS:
        for fasta in search_dir.glob(ext):
            if fasta.name.startswith("."):
                continue
            # Check companion files
            has_fai = (fasta.parent / (fasta.name + ".fai")).is_file()
            has_dict = any(
                (fasta.parent / fasta.name.rsplit(".", 1)[0]).with_suffix(".dict").is_file()
                for _ in [None]
            )
            has_bwa = (
                any(
                    (fasta.parent / (fasta.name + ext)).is_file()
                    for ext in BWA_INDEX_EXTS[:5]  # check classic first
                )
                or any(
                    (fasta.parent / (fasta.name + ext)).is_file()
                    for ext in BWA_INDEX_EXTS[5:]  # check 64-bit
                )
            )
            # Check for .gz companion
            gz_path = fasta.parent / (fasta.name + ".gz")
            has_gz = gz_path.is_file()
            results.append(
                {
                    "path": str(fasta),
                    "gz_path": str(gz_path) if has_gz else "",
                    "has_fai": has_fai,
                    "has_dict": has_dict,
                    "has_bwa": has_bwa,
                    "has_gz": has_gz,
                    "name": fasta.name,
                }
            )
    # Also check for .gz-only references (BWA index built from .gz)
    for ext in GENOME_EXTENSIONS:
        for fasta_gz in search_dir.glob(ext + ".gz"):
            uncompressed = fasta_gz.parent / fasta_gz.name[:-3]
            if uncompressed.is_file():
                continue  # already found above
            has_bwa = any(
                (fasta_gz.parent / (fasta_gz.name + ext)).is_file() for ext in BWA_INDEX_EXTS
            )
            if has_bwa:
                results.append(
                    {
                        "path": "",
                        "gz_path": str(fasta_gz),
                        "has_fai": False,
                        "has_dict": False,
                        "has_bwa": has_bwa,
                        "has_gz": True,
                        "name": fasta_gz.name,
                    }
                )
    return results


def _find_known_sites_vcfs(
    search_dir: Path,
    patterns: tuple[str, ...] = BQSR_KNOWN_SITES_PATTERNS,
) -> list[str]:
    """Find known-sites VCF files matching GATK resource bundle patterns."""
    results: list[str] = []
    if not search_dir.is_dir():
        return results
    for pattern in patterns:
        for vcf in search_dir.glob(pattern):
            # Only include .vcf.gz or .vcf (not .vcf.gz.tbi)
            if vcf.name.endswith(".tbi") or vcf.name.endswith(".idx"):
                continue
            # Check for tabix index
            has_index = (vcf.parent / (vcf.name + ".tbi")).is_file() or (
                vcf.parent / (vcf.name + ".idx")
            ).is_file()
            if not has_index:
                print(f"  Warning: No index for {vcf.name}", file=sys.stderr)
            results.append(str(vcf))
    # Deduplicate (patterns may overlap)
    return sorted(set(results))


def discover_reference_data(
    ref_dir: Path | None,
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Scan for reference genome and known-sites VCFs.

    Searches in order: explicit --ref-dir, relative project paths, shared HPC locations.

    Returns:
        Dict with keys: genome, genome_gz, build, known_sites, search_log.
    """
    search_log: list[str] = []
    genomes: list[dict[str, Any]] = []
    known_sites: list[str] = []

    # Build search order for reference genomes
    ref_search = []
    if ref_dir:
        ref_search.append(ref_dir)
    if project_root:
        for rel in REF_SEARCH_DIRS:
            ref_search.append(project_root / rel)
    for shared in SHARED_REF_DIRS:
        ref_search.append(Path(shared))

    # Search for genomes
    for search_path in ref_search:
        if not search_path.is_dir():
            continue
        found = _find_genome_fastas(search_path)
        if found:
            search_log.append(f"  Found {len(found)} genome(s) in {search_path}")
            for g in found:
                status = []
                if g["has_bwa"]:
                    status.append("BWA")
                if g["has_fai"]:
                    status.append("FAI")
                if g["has_dict"]:
                    status.append("Dict")
                search_log.append(f"    {g['name']}  [{', '.join(status) or 'no indexes'}]")
            genomes.extend(found)

    # Build search order for known-sites
    ks_search = []
    if ref_dir:
        ks_search.append(ref_dir)
        # Also check parent and sibling dirs
        if ref_dir.parent.is_dir():
            for sibling in ("gatk_bundle", "GATK_resource_bundle", "known_sites"):
                candidate = ref_dir.parent / sibling
                if candidate.is_dir():
                    ks_search.append(candidate)
    if project_root:
        for rel in KNOWN_SITES_SEARCH_DIRS:
            ks_search.append(project_root / rel)
    for shared in SHARED_KNOWN_SITES_DIRS:
        ks_search.append(Path(shared))

    # Search for known-sites
    for search_path in ks_search:
        if not search_path.is_dir():
            continue
        found_vcfs = _find_known_sites_vcfs(search_path)
        if found_vcfs:
            search_log.append(f"  Found {len(found_vcfs)} known-sites VCF(s) in {search_path}")
            for vcf_path in found_vcfs:
                search_log.append(f"    {Path(vcf_path).name}")
            known_sites.extend(found_vcfs)

    # Deduplicate known-sites (same file found via different search paths)
    seen_names = set()
    unique_ks = []
    for ks in known_sites:
        name = Path(ks).name
        if name not in seen_names:
            seen_names.add(name)
            unique_ks.append(ks)
    known_sites = unique_ks

    # Pick best genome (prefer one with BWA index + uncompressed)
    genome_path = ""
    genome_gz_path = ""
    build = "GRCh38"
    if genomes:
        # Sort: prefer BWA-indexed, then with FAI, then with uncompressed path
        best = sorted(
            genomes,
            key=lambda g: (g["has_bwa"], g["has_fai"], g["has_dict"], bool(g["path"])),
            reverse=True,
        )[0]
        genome_path = str(Path(best["path"]).resolve()) if best["path"] else ""
        # Verify BWA index exists for the .gz before setting genome_gz (#25)
        if best.get("has_gz") and best.get("gz_path"):
            gz = Path(best["gz_path"])
            gz_has_bwa = any((gz.parent / (gz.name + ext)).is_file() for ext in BWA_INDEX_EXTS[:5])
            if gz_has_bwa:
                genome_gz_path = str(gz.resolve())
            else:
                search_log.append(f"  Warning: No BWA index for {gz.name}, genome_gz left empty")
                genome_gz_path = ""
        # Infer build from filename
        name_lower = best["name"].lower()
        if "grch37" in name_lower or "hg19" in name_lower or "hs37" in name_lower:
            build = "GRCh37"

    # Resolve known-sites paths (#23)
    known_sites = [str(Path(ks).resolve()) for ks in known_sites]

    if not genomes:
        search_log.append("  No reference genome found")
    if not known_sites:
        search_log.append("  No known-sites VCFs found")

    return {
        "genome": genome_path,
        "genome_gz": genome_gz_path,
        "build": build,
        "known_sites": known_sites,
        "search_log": search_log,
    }


# ---------------------------------------------------------------------------
# Config template generation
# ---------------------------------------------------------------------------


def _resolve_path(p: str) -> str:
    """Resolve a path string to canonical form, skipping placeholders."""
    if not p or "EDIT_ME" in p:
        return p
    return str(Path(p).resolve())


def _build_config_yaml(
    ref_data: dict[str, Any],
    fastq_folder: str,
    output_folder: str,
    samples_path: str,
) -> str:
    """Build config.yaml content from discovered reference data and paths."""
    genome = ref_data.get("genome", "") or "EDIT_ME: /path/to/reference.fna"
    genome_gz = ref_data.get("genome_gz", "") or ""
    build = ref_data.get("build", "GRCh38")
    known_sites = ref_data.get("known_sites", [])

    # Canonicalize paths (#23) — resolve ../ segments
    genome = _resolve_path(genome)
    genome_gz = _resolve_path(genome_gz)
    fastq_folder = _resolve_path(fastq_folder)
    known_sites = [_resolve_path(ks) for ks in known_sites]

    # Format known-sites list
    if known_sites:
        ks_lines = "\n".join(f'    - "{ks}"' for ks in known_sites)
    else:
        ks_lines = '    - "EDIT_ME: /path/to/known_sites.vcf.gz"'

    return f"""\
# =============================================================================
# config/config.yaml -- sm-alignment pipeline configuration
# =============================================================================
# Generated by: scripts/generate_config.py
# Review all paths below before running the pipeline.

# --- Reference genome & known variant sites ---
ref:
  genome: "{genome}"
  genome_gz: "{genome_gz}"
  build: "{build}"
  known_sites:
{ks_lines}

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
    ref_data: dict[str, Any],
    dry_run: bool = False,
    force: bool = False,
    output_dir: str | None = None,
) -> None:
    """Generate config.yaml with discovered reference paths.

    Args:
        config_output: Where to write the config file.
        fastq_dir: FASTQ directory (fills paths.fastq_folder).
        samples_path: Value for paths.samples.
        project: Project identifier (used in output folder path).
        ref_data: Discovered reference data from discover_reference_data().
        dry_run: If True, print content but do not write.
        force: If True, overwrite existing file without asking.
        output_dir: Pipeline output directory (default: results/<project>).
    """
    fastq_folder = str(fastq_dir).replace("\\", "/")
    output_folder = output_dir if output_dir else f"results/{project}"

    content = _build_config_yaml(ref_data, fastq_folder, output_folder, samples_path)

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
    samplesheet_entries: list[dict[str, str]] | None,
    fastq_pairs: list[dict[str, str]],
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


def _prompt(prompt: str, default: str = "") -> str:
    """Prompt the user for input with an optional default value."""
    if default:
        result = input(f"{prompt} [{default}]: ").strip()
        return result if result else default
    return input(f"{prompt}: ").strip()


def _prompt_yn(prompt: str, default: bool = True) -> bool:
    """Prompt for yes/no with a default."""
    suffix = "[Y/n]" if default else "[y/N]"
    result = input(f"{prompt} {suffix}: ").strip().lower()
    if not result:
        return default
    return result in ("y", "yes")


def _prompt_path(prompt: str, default: str = "", must_exist: bool = True) -> str:
    """Prompt for a filesystem path, re-prompting on invalid input."""
    while True:
        raw = _prompt(prompt, default)
        if not raw:
            if not must_exist:
                return ""
            print("  Path cannot be empty. Please try again.")
            continue
        p = Path(raw).expanduser()
        if must_exist and not p.exists():
            print(f"  Path does not exist: {p}")
            retry = _prompt_yn("  Try again?", default=True)
            if not retry:
                return str(p)
            continue
        return str(p)


def interactive_mode() -> None:
    """Guided wizard for generating pipeline config files."""
    print()
    print("=" * 60)
    print("  sm-alignment — Config Generator (interactive)")
    print("=" * 60)
    print()

    # 1. FASTQ directory
    fastq_dir_str = _prompt_path("FASTQ directory")
    fastq_dir = Path(fastq_dir_str).resolve()

    # 2. SampleSheet
    ss_path: Path | None = find_samplesheet(fastq_dir)
    if ss_path:
        print(f"\n  Auto-detected SampleSheet: {ss_path}")
        use_detected = _prompt_yn("  Use this SampleSheet?", default=True)
        if not use_detected:
            custom = _prompt_path(
                "  Path to SampleSheet.csv (leave empty to skip)", must_exist=False
            )
            ss_path = Path(custom) if custom else None
    else:
        print("\n  No SampleSheet.csv found in FASTQ directory.")
        custom = _prompt_path("  Path to SampleSheet.csv (leave empty to skip)", must_exist=False)
        ss_path = Path(custom) if custom else None

    # 3. Project name
    auto_project = infer_project(fastq_dir, None)
    project = _prompt("Project identifier", default=auto_project)

    # 4. Output paths
    samples_output = _prompt("Output samples.tsv path", default="config/samples.tsv")

    # 5. Generate config.yaml?
    gen_config = _prompt_yn("\nAlso generate config/config.yaml?", default=True)

    ref_dir: Path | None = None
    config_output = "config/config.yaml"
    output_dir: str | None = None
    if gen_config:
        config_output = _prompt("Output config.yaml path", default="config/config.yaml")
        default_output_dir = f"results/{project}"
        output_dir_str = _prompt(
            "Pipeline output directory (paths.output_folder)", default=default_output_dir
        )
        if output_dir_str != default_output_dir:
            output_dir = output_dir_str

        # 6. Reference data
        print("\n  Reference data discovery")
        print("  The script can scan for reference genome and known-sites VCFs.")
        custom_ref = _prompt_path(
            "  Reference data directory (leave empty to auto-scan)",
            must_exist=False,
        )
        if custom_ref:
            ref_dir = Path(custom_ref).resolve()

    # 7. Dry-run or write?
    dry_run = not _prompt_yn("\nWrite files now?", default=True)
    force = False
    if not dry_run:
        force = _prompt_yn("Overwrite existing files without asking?", default=False)

    # ---- Execute with collected parameters ----
    print()

    # Discover FASTQs
    fastq_pairs = discover_fastq_files(fastq_dir)

    # Parse SampleSheet
    samplesheet_entries: list[dict[str, str]] | None = None
    if ss_path and ss_path.is_file():
        samplesheet_entries = parse_samplesheet(ss_path)

    # Build table
    df = build_samples_table(fastq_pairs, samplesheet_entries, project)
    print_summary(df, samplesheet_entries, fastq_pairs)

    # Write samples.tsv
    write_samples_tsv(df, Path(samples_output), dry_run=dry_run, force=force)

    # Generate config.yaml
    if gen_config:
        project_root = fastq_dir.parent
        print("Scanning for reference data...")
        ref_data = discover_reference_data(ref_dir, project_root)
        for line in ref_data.get("search_log", []):
            print(line)
        print()

        generate_config_template(
            config_output=Path(config_output),
            fastq_dir=fastq_dir,
            samples_path=samples_output,
            project=project,
            ref_data=ref_data,
            dry_run=dry_run,
            force=force,
            output_dir=output_dir,
        )

    print("\nDone.")


def main() -> None:
    """Main entry point: parse arguments, discover files, generate output."""
    parser = argparse.ArgumentParser(
        description="Generate config/samples.tsv from sequencing facility deliverables.\n"
        "Run without arguments for an interactive guided wizard.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s                                                   # interactive wizard
  %(prog)s --fastq-dir /data/project/fastqs                  # flags mode
  %(prog)s --fastq-dir /data/project/fastqs --config-template --dry-run
  %(prog)s --fastq-dir /data/project/fastqs --project A5297 --force
        """,
    )
    parser.add_argument(
        "--fastq-dir",
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
        "--output-dir",
        help="Pipeline output directory for BAM files etc. "
        "(default: results/<project>). Written into paths.output_folder in config.yaml",
    )
    parser.add_argument(
        "--ref-dir",
        help="Directory containing reference genome and/or known-sites VCFs "
        "(auto-scans common locations if not specified)",
    )
    parser.add_argument(
        "--config-template",
        action="store_true",
        help="Also generate a config.yaml with discovered reference paths",
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

    # If no --fastq-dir provided, launch interactive wizard
    if not args.fastq_dir:
        interactive_mode()
        return

    fastq_dir = Path(args.fastq_dir).resolve()

    # ---- Discover FASTQ files ----
    fastq_pairs = discover_fastq_files(fastq_dir)

    # ---- Parse SampleSheet (if available) ----
    samplesheet_entries: list[dict[str, str]] | None = None

    if args.samplesheet:
        ss_explicit = Path(args.samplesheet)
        if not ss_explicit.is_file():
            sys.exit(f"Error: SampleSheet not found: {ss_explicit}")
        samplesheet_entries = parse_samplesheet(ss_explicit)
    else:
        ss_auto = find_samplesheet(fastq_dir)
        if ss_auto:
            samplesheet_entries = parse_samplesheet(ss_auto)

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
        ref_dir = Path(args.ref_dir).resolve() if args.ref_dir else None
        # Use parent of fastq_dir as project root for relative searches
        project_root = fastq_dir.parent

        print("Scanning for reference data...")
        ref_data = discover_reference_data(ref_dir, project_root)
        for line in ref_data.get("search_log", []):
            print(line)
        print()

        config_output = Path(args.config_output)
        generate_config_template(
            config_output=config_output,
            fastq_dir=fastq_dir,
            samples_path=args.output,
            project=project,
            ref_data=ref_data,
            dry_run=args.dry_run,
            force=args.force,
            output_dir=args.output_dir,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
