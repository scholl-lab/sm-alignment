"""Shared fixtures for sm-alignment test suite."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# SampleSheet content fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def illumina_samplesheet_lines():
    """Standard Illumina SampleSheet with [Header]/[Data] sections."""
    return [
        "[Header]\n",
        "IEMFileVersion,5\n",
        "Investigator Name,Test\n",
        "[Reads]\n",
        "151\n",
        "[Data]\n",
        "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project\n",
        "1,S1,SampleA,ATCACG,TTAGGC,TestProject\n",
        "1,S2,SampleB,CGATGT,TTAGGC,TestProject\n",
        "2,S1,SampleA,ATCACG,TTAGGC,TestProject\n",
    ]


@pytest.fixture()
def minimal_5col_lines():
    """BIH/Charite minimal SampleSheet (5 columns, no section headers)."""
    return [
        "1,SampleA,ATCACG,TTAGGC,TestProject\n",
        "1,SampleB,CGATGT,TTAGGC,TestProject\n",
        "2,SampleA,ATCACG,TTAGGC,TestProject\n",
    ]


@pytest.fixture()
def minimal_4col_lines():
    """Minimal SampleSheet with 4 columns (no project)."""
    return [
        "1,SampleA,ATCACG,TTAGGC\n",
        "1,SampleB,CGATGT,TTAGGC\n",
    ]


# ---------------------------------------------------------------------------
# FASTQ directory fixtures
# ---------------------------------------------------------------------------

FASTQ_CONTENT = b"@read1\nACGTACGT\n+\nIIIIIIII\n"


def _create_fastq(path: Path) -> None:
    """Create a tiny valid gzip FASTQ file."""
    with gzip.open(path, "wb") as f:
        f.write(FASTQ_CONTENT)


@pytest.fixture()
def fastq_dir(tmp_path):
    """Temporary directory with 6 valid FASTQ files (3 pairs)."""
    fq_dir = tmp_path / "fastqs"
    fq_dir.mkdir()
    names = [
        "SampleA_S1_L001_R1_001.fastq.gz",
        "SampleA_S1_L001_R2_001.fastq.gz",
        "SampleA_S1_L002_R1_001.fastq.gz",
        "SampleA_S1_L002_R2_001.fastq.gz",
        "SampleB_S2_L001_R1_001.fastq.gz",
        "SampleB_S2_L001_R2_001.fastq.gz",
    ]
    for name in names:
        _create_fastq(fq_dir / name)
    return fq_dir


@pytest.fixture()
def fastq_dir_with_orphan(fastq_dir):
    """FASTQ directory with an extra R1 that has no R2 mate."""
    _create_fastq(fastq_dir / "OrphanSample_S9_L001_R1_001.fastq.gz")
    return fastq_dir


# ---------------------------------------------------------------------------
# DataFrame fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def samples_df():
    """Pre-built samples DataFrame (no subfolder column)."""
    data = {
        "fastq_files_basename": [
            "SampleA_S1_L001",
            "SampleA_S1_L002",
            "SampleB_S2_L001",
        ],
        "lane": ["L001", "L002", "L001"],
        "project_sample": ["SampleA", "SampleA", "SampleB"],
        "mdc_project": ["TestProject", "TestProject", "TestProject"],
    }
    return pd.DataFrame(data).set_index("fastq_files_basename", drop=False)


@pytest.fixture()
def samples_df_with_subfolder():
    """Pre-built samples DataFrame with subfolder column."""
    data = {
        "fastq_files_basename": [
            "SampleA_S1_L001",
            "SampleA_S1_L002",
            "SampleB_S2_L001",
        ],
        "lane": ["L001", "L002", "L001"],
        "project_sample": ["SampleA", "SampleA", "SampleB"],
        "mdc_project": ["TestProject", "TestProject", "TestProject"],
        "subfolder": ["run1", "run1", "run2"],
    }
    return pd.DataFrame(data).set_index("fastq_files_basename", drop=False)


# ---------------------------------------------------------------------------
# Reference / known-sites directory fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def ref_dir(tmp_path):
    """Temporary directory with skeleton reference genome files."""
    d = tmp_path / "ref"
    d.mkdir()
    base = d / "dummy.fna"
    base.touch()
    for ext in (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"):
        (d / f"dummy.fna{ext}").touch()
    (d / "dummy.dict").touch()
    return d


@pytest.fixture()
def known_sites_dir(tmp_path):
    """Temporary directory with skeleton known-sites VCF files."""
    d = tmp_path / "known_sites"
    d.mkdir()
    for name in ("dbsnp_mini.vcf.gz", "Mills_and_1000G_mini.vcf.gz"):
        (d / name).touch()
        (d / f"{name}.tbi").touch()
    return d


# ---------------------------------------------------------------------------
# FASTQ pairs fixture (as returned by discover_fastq_files)
# ---------------------------------------------------------------------------


@pytest.fixture()
def fastq_pairs(fastq_dir):
    """Pre-built list matching discover_fastq_files() return format."""
    return [
        {
            "basename": "SampleA_S1_L001",
            "lane": "001",
            "lane_str": "L001",
            "sample_name": "SampleA",
            "r1_path": str(fastq_dir / "SampleA_S1_L001_R1_001.fastq.gz"),
            "r2_path": str(fastq_dir / "SampleA_S1_L001_R2_001.fastq.gz"),
        },
        {
            "basename": "SampleA_S1_L002",
            "lane": "002",
            "lane_str": "L002",
            "sample_name": "SampleA",
            "r1_path": str(fastq_dir / "SampleA_S1_L002_R1_001.fastq.gz"),
            "r2_path": str(fastq_dir / "SampleA_S1_L002_R2_001.fastq.gz"),
        },
        {
            "basename": "SampleB_S2_L001",
            "lane": "001",
            "lane_str": "L001",
            "sample_name": "SampleB",
            "r1_path": str(fastq_dir / "SampleB_S2_L001_R1_001.fastq.gz"),
            "r2_path": str(fastq_dir / "SampleB_S2_L001_R2_001.fastq.gz"),
        },
    ]
