"""Unit tests for workflow/rules/helpers.py — extracted pipeline helpers."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import pytest

# Add workflow/rules to path so we can import helpers directly
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "workflow" / "rules"))
from helpers import (
    get_basenames_for_sample,
    get_java_opts,
    get_samples,
    resolve_fastq_path,
)

# ============================================================================
# TestGetJavaOpts
# ============================================================================


class TestGetJavaOpts:
    """Tests for get_java_opts()."""

    def test_standard_allocation(self):
        result = get_java_opts(8000, "/tmp")
        assert "-Xmx6400m" in result
        assert "-Xms1600m" in result
        assert "-Djava.io.tmpdir=/tmp" in result

    def test_large_allocation(self):
        result = get_java_opts(32000, "/scratch/tmp")
        assert "-Xmx25600m" in result
        assert "-Xms6400m" in result
        assert "-Djava.io.tmpdir=/scratch/tmp" in result

    def test_small_allocation(self):
        result = get_java_opts(1000, "/tmp")
        assert "-Xmx800m" in result
        assert "-Xms200m" in result

    def test_return_type(self):
        result = get_java_opts(4000, "/tmp")
        assert isinstance(result, str)

    @pytest.mark.parametrize(
        "mem_mb,expected_xmx,expected_xms",
        [
            (1000, 800, 200),
            (4000, 3200, 800),
            (8000, 6400, 1600),
            (16000, 12800, 3200),
        ],
    )
    def test_80_20_split(self, mem_mb, expected_xmx, expected_xms):
        result = get_java_opts(mem_mb, "/tmp")
        assert f"-Xmx{expected_xmx}m" in result
        assert f"-Xms{expected_xms}m" in result


# ============================================================================
# TestGetSamples
# ============================================================================


class TestGetSamples:
    """Tests for get_samples()."""

    def test_unique_sorted(self, samples_df):
        result = get_samples(samples_df)
        assert result == ["SampleA", "SampleB"]

    def test_single_sample(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["X_S1_L001"],
                "project_sample": ["OnlySample"],
            }
        )
        result = get_samples(df)
        assert result == ["OnlySample"]

    def test_sort_order(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["z_S1_L001", "a_S1_L001"],
                "project_sample": ["Zebra", "Alpha"],
            }
        )
        result = get_samples(df)
        assert result == ["Alpha", "Zebra"]


# ============================================================================
# TestGetBasenamesForSample
# ============================================================================


class TestGetBasenamesForSample:
    """Tests for get_basenames_for_sample()."""

    def test_multi_lane(self, samples_df):
        result = get_basenames_for_sample(samples_df, "SampleA")
        assert sorted(result) == ["SampleA_S1_L001", "SampleA_S1_L002"]

    def test_single_lane(self, samples_df):
        result = get_basenames_for_sample(samples_df, "SampleB")
        assert result == ["SampleB_S2_L001"]

    def test_nonexistent_sample(self, samples_df):
        result = get_basenames_for_sample(samples_df, "NoSuchSample")
        assert result == []


# ============================================================================
# TestResolveFastqPath
# ============================================================================


class TestResolveFastqPath:
    """Tests for resolve_fastq_path()."""

    def test_raw_mode(self, samples_df):
        result = resolve_fastq_path(
            "SampleA_S1_L001",
            "_R1_001.fastq.gz",
            samples_df,
            trimming_enabled=False,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
        )
        assert result == os.path.join("/data/fastqs", "SampleA_S1_L001_R1_001.fastq.gz")

    def test_trimmed_mode(self, samples_df):
        result = resolve_fastq_path(
            "SampleA_S1_L001",
            ".bbduk_R1_001.fastq.gz",
            samples_df,
            trimming_enabled=True,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
        )
        assert result == os.path.join("/out/trimmed", "SampleA_S1_L001.bbduk_R1_001.fastq.gz")

    def test_subfolder(self, samples_df_with_subfolder):
        result = resolve_fastq_path(
            "SampleA_S1_L001",
            "_R1_001.fastq.gz",
            samples_df_with_subfolder,
            trimming_enabled=False,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
        )
        assert result == os.path.join("/data/fastqs", "run1", "SampleA_S1_L001_R1_001.fastq.gz")

    def test_subfolder_with_trimming(self, samples_df_with_subfolder):
        result = resolve_fastq_path(
            "SampleA_S1_L001",
            ".bbduk_R1_001.fastq.gz",
            samples_df_with_subfolder,
            trimming_enabled=True,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
        )
        assert result == os.path.join(
            "/out/trimmed", "run1", "SampleA_S1_L001.bbduk_R1_001.fastq.gz"
        )

    def test_r2_suffix(self, samples_df):
        result = resolve_fastq_path(
            "SampleB_S2_L001",
            "_R2_001.fastq.gz",
            samples_df,
            trimming_enabled=False,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
        )
        assert result == os.path.join("/data/fastqs", "SampleB_S2_L001_R2_001.fastq.gz")
