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
    format_bwa_extra,
    format_mark_duplicates_extra,
    format_read_group,
    get_all_basenames,
    get_basenames_for_sample,
    get_java_opts,
    get_sample_kit,
    get_samples,
    load_samples,
    load_umi_presets,
    resolve_fastq_path,
    resolve_umi_read_structures,
    validate_umi_config,
    write_umi_provenance,
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
# TestGetAllBasenames
# ============================================================================


class TestGetAllBasenames:
    """Tests for get_all_basenames()."""

    def test_returns_all_sorted(self, samples_df):
        result = get_all_basenames(samples_df)
        assert result == ["SampleA_S1_L001", "SampleA_S1_L002", "SampleB_S2_L001"]

    def test_single_entry(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["X_S1_L001"],
                "project_sample": ["OnlySample"],
            }
        )
        result = get_all_basenames(df)
        assert result == ["X_S1_L001"]

    def test_no_duplicates(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["A_S1_L001", "A_S1_L001", "B_S2_L001"],
                "project_sample": ["A", "A", "B"],
            }
        )
        result = get_all_basenames(df)
        assert result == ["A_S1_L001", "B_S2_L001"]


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
        assert result == os.path.join("/out/trimmed", "SampleA_S1_L001.bbduk_R1_001.fastq.gz")

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


# ============================================================================
# TestFormatReadGroup
# ============================================================================


class TestFormatReadGroup:
    """Tests for format_read_group()."""

    def test_with_flowcell_and_library(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["S1_L001", "S1_L001_run2"],
                "lane": ["L001", "L001"],
                "project_sample": ["APA1-T", "APA1-T"],
                "mdc_project": ["A5297", "A5297"],
                "flowcell": ["FC1", "FC2"],
                "library": ["LIB1", "LIB2"],
            }
        ).set_index("fastq_files_basename", drop=False)
        rg1 = format_read_group(df, "S1_L001", "ILLUMINA")
        rg2 = format_read_group(df, "S1_L001_run2", "ILLUMINA")
        assert rg1 == r"@RG\tID:FC1.L001.LIB1\tSM:APA1-T\tLB:LIB1\tPL:ILLUMINA\tPU:FC1.L001.LIB1"
        assert rg2 == r"@RG\tID:FC2.L001.LIB2\tSM:APA1-T\tLB:LIB2\tPL:ILLUMINA\tPU:FC2.L001.LIB2"
        assert rg1 != rg2

    def test_fallback_to_subfolder_and_sample(self, samples_df_with_subfolder):
        rg = format_read_group(samples_df_with_subfolder, "SampleA_S1_L001", "ILLUMINA")
        assert (
            rg
            == r"@RG\tID:run1.L001.SampleA\tSM:SampleA\tLB:SampleA\tPL:ILLUMINA\tPU:run1.L001.SampleA"
        )

    def test_fallback_to_mdc_project_and_sample(self, samples_df):
        rg = format_read_group(samples_df, "SampleA_S1_L001", "ILLUMINA")
        assert (
            rg
            == r"@RG\tID:TestProject.L001.SampleA\tSM:SampleA\tLB:SampleA\tPL:ILLUMINA\tPU:TestProject.L001.SampleA"
        )


# ============================================================================
# TestLoadSamples
# ============================================================================


class TestLoadSamples:
    """Tests for load_samples()."""

    def test_numeric_ids_and_empty_subfolder(self, tmp_path):
        tsv = tmp_path / "samples.tsv"
        tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\tsubfolder\n"
            "1001_S1_L001\t1\t1001\t5297\t\n"
            "1002_S2_L001\t2\t1002\t5297\trun2\n"
        )
        df = load_samples(tsv)
        assert df.loc["1001_S1_L001", "project_sample"] == "1001"
        assert df.loc["1001_S1_L001", "lane"] == "1"
        assert df.loc["1001_S1_L001", "subfolder"] == ""
        path = resolve_fastq_path(
            "1001_S1_L001",
            "_R1_001.fastq.gz",
            df,
            trimming_enabled=False,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
        )
        assert path == "/data/fastqs/1001_S1_L001_R1_001.fastq.gz"
        assert "nan" not in path

    def test_duplicate_basenames_rejected(self, tmp_path):
        tsv = tmp_path / "samples.tsv"
        tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\n"
            "dup_S1_L001\tL001\tSampleA\tProj\n"
            "dup_S1_L001\tL002\tSampleB\tProj\n"
        )
        with pytest.raises(ValueError, match="Duplicate fastq_files_basename"):
            load_samples(tsv)

    def test_path_separator_in_basename_rejected(self, tmp_path):
        tsv = tmp_path / "samples.tsv"
        tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\n"
            "sub/S1_L001\tL001\tSampleA\tProj\n"
        )
        with pytest.raises(ValueError, match="path separators"):
            load_samples(tsv)

    def test_kit_per_sample_and_conflict_rejected(self, tmp_path):
        tsv = tmp_path / "samples.tsv"
        tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\tkit\n"
            "S1_L001\tL001\tSampleA\tProj\tS33266340\n"
            "S1_L002\tL002\tSampleA\tProj\tS33266340\n"
            "S2_L001\tL001\tSampleB\tProj\tS33266436\n"
        )
        df = load_samples(tsv)
        assert get_sample_kit(df, "SampleA") == "S33266340"
        assert get_sample_kit(df, "SampleB") == "S33266436"

        bad_tsv = tmp_path / "bad_samples.tsv"
        bad_tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\tkit\n"
            "S1_L001\tL001\tSampleA\tProj\tS33266340\n"
            "S1_L002\tL002\tSampleA\tProj\tS33266436\n"
        )
        with pytest.raises(ValueError, match="conflicting kits"):
            load_samples(bad_tsv)


# ============================================================================
# TestUmiHelpers
# ============================================================================


class TestUmiHelpers:
    """Tests for UMI preset resolution, validation, provenance, and CLI flag helpers."""

    def test_load_umi_presets_contains_agilent_xths2(self):
        presets = load_umi_presets()
        assert "agilent_xths2" in presets
        assert presets["agilent_xths2"]["read_structures"] == "3M2S+T 3M2S+T"
        assert "G9985-90010" in presets["agilent_xths2"]["source"]

    def test_preset_resolves_to_right_structure(self):
        cfg = {"enabled": True, "preset": "agilent_xths2", "read_structures": None}
        assert resolve_umi_read_structures(cfg) == "3M2S+T 3M2S+T"

    def test_custom_structure_overrides_preset(self):
        for preset in (None, "custom", "agilent_xths2"):
            cfg = {
                "enabled": True,
                "preset": preset,
                "read_structures": "5M2S+T 5M2S+T",
            }
            assert resolve_umi_read_structures(cfg) == "5M2S+T 5M2S+T"
        assert (
            resolve_umi_read_structures({"enabled": True, "read_structures": "3M2S+T\n3M2S+T"})
            == "3M2S+T 3M2S+T"
        )

    def test_disabled_returns_empty(self):
        assert resolve_umi_read_structures({"enabled": False}) == ""
        assert resolve_umi_read_structures(None) == ""

    def test_validate_umi_config_rejects_conflict(self):
        with pytest.raises(ValueError, match="umi config conflict"):
            validate_umi_config(
                {
                    "enabled": True,
                    "preset": "agilent_xths2",
                    "read_structures": "5M+T 5M+T",
                }
            )

    def test_validate_umi_config_rejects_missing_both(self):
        with pytest.raises(ValueError, match="must specify either"):
            validate_umi_config({"enabled": True, "preset": None, "read_structures": None})
        with pytest.raises(ValueError, match="must specify either"):
            validate_umi_config({"enabled": True, "preset": "custom", "read_structures": ""})

    def test_validate_umi_config_rejects_invalid_structure(self):
        with pytest.raises(ValueError, match="Invalid fgbio read_structures"):
            validate_umi_config({"enabled": True, "read_structures": "3X+T 3X+T"})

    def test_validate_umi_config_rejects_unknown_preset(self):
        with pytest.raises(ValueError, match=r"Unknown umi\.preset"):
            validate_umi_config({"enabled": True, "preset": "unknown_kit"})

    def test_validate_umi_config_rejects_barcode_tag_mismatch(self):
        with pytest.raises(
            ValueError, match=r"umi\.dedup\.barcode_tag 'RX' must match umi\.umi_tag 'OX'"
        ):
            validate_umi_config(
                {
                    "enabled": True,
                    "preset": "agilent_xths2",
                    "umi_tag": "OX",
                    "dedup": {"barcode_tag": "RX"},
                }
            )

    def test_write_umi_provenance_logs_and_writes_yaml(self, tmp_path, caplog):
        import logging

        import yaml

        cfg = {"enabled": True, "preset": "agilent_xths2", "umi_tag": "RX", "duplex": True}
        with caplog.at_level(logging.INFO):
            out = write_umi_provenance(tmp_path / "logs", cfg)
        assert out is not None and out.exists()
        assert "read_structures='3M2S+T 3M2S+T'" in caplog.text
        data = yaml.safe_load(out.read_text())
        assert data["umi"]["preset"] == "agilent_xths2"
        assert data["umi"]["read_structures"] == "3M2S+T 3M2S+T"
        assert data["umi"]["umi_tag"] == "RX"
        assert data["umi"]["duplex"] is True
        assert "G9985-90010" in data["umi"]["source"]

        assert write_umi_provenance(tmp_path / "logs", {"enabled": False}) is None

    def test_format_bwa_extra(self):
        assert format_bwa_extra("-K 100000000 -Y", umi_enabled=False) == "-K 100000000 -Y"
        assert format_bwa_extra("-K 100000000 -Y", umi_enabled=True) == "-K 100000000 -Y -C"
        assert format_bwa_extra("-K 100000000 -Y -C", umi_enabled=True) == "-K 100000000 -Y -C"

    def test_format_mark_duplicates_extra(self):
        import yaml
        from snakemake.utils import validate as snakemake_validate

        base = "--VALIDATION_STRINGENCY LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
        assert format_mark_duplicates_extra(base, {"enabled": False}) == base
        assert (
            format_mark_duplicates_extra(base, {"enabled": True, "duplex": True, "umi_tag": "RX"})
            == f"{base} --BARCODE_TAG RX --DUPLEX_UMI true"
        )
        assert (
            format_mark_duplicates_extra(base, {"enabled": True, "duplex": False, "umi_tag": "RX"})
            == f"{base} --BARCODE_TAG RX"
        )
        assert (
            format_mark_duplicates_extra(
                base,
                {
                    "enabled": True,
                    "duplex": True,
                    "dedup": {"barcode_tag": None, "extra": "--TAGGING_POLICY All"},
                },
            )
            == f"{base} --TAGGING_POLICY All"
        )

        repo_root = Path(__file__).resolve().parent.parent
        cfg = yaml.safe_load((repo_root / "tests" / "data" / "config" / "config.yaml").read_text())
        cfg["umi"] = {"enabled": True, "umi_tag": "OX", "preset": "agilent_xths2"}
        snakemake_validate(cfg, str(repo_root / "workflow" / "schemas" / "config.schema.yaml"))
        validate_umi_config(cfg["umi"])
        assert format_mark_duplicates_extra("", cfg["umi"]) == "--BARCODE_TAG OX --DUPLEX_UMI true"

    def test_resolve_fastq_path_umi_mode(self, samples_df):
        result = resolve_fastq_path(
            "SampleA_S1_L001",
            ".umi_R1_001.fastq.gz",
            samples_df,
            trimming_enabled=False,
            trimmed_dir="/out/trimmed",
            fastq_dir="/data/fastqs",
            umi_enabled=True,
            umi_dir="/out/umi",
        )
        assert result == os.path.join("/out/umi", "SampleA_S1_L001.umi_R1_001.fastq.gz")
