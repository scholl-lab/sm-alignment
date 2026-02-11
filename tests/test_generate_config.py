"""Unit tests for scripts/generate_config.py — config generator."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

# Ensure the scripts directory is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import generate_config as gc

# ============================================================================
# TestDetectSamplesheetFormat
# ============================================================================


class TestDetectSamplesheetFormat:
    """Tests for detect_samplesheet_format()."""

    def test_illumina_data_section(self, illumina_samplesheet_lines):
        assert gc.detect_samplesheet_format(illumina_samplesheet_lines) == "illumina"

    def test_bclconvert_data_section(self):
        lines = ["[BCLConvert_Data]\n", "Lane,Sample_ID\n", "1,SampleA\n"]
        assert gc.detect_samplesheet_format(lines) == "illumina"

    def test_minimal_format(self, minimal_5col_lines):
        assert gc.detect_samplesheet_format(minimal_5col_lines) == "minimal"

    def test_empty_lines(self):
        assert gc.detect_samplesheet_format([]) == "minimal"

    def test_header_section_only(self):
        lines = ["[Header]\n", "Key,Value\n"]
        assert gc.detect_samplesheet_format(lines) == "illumina"

    def test_data_with_surrounding_commas(self):
        lines = ["[Data],\n", "Lane,Sample_ID\n"]
        assert gc.detect_samplesheet_format(lines) == "illumina"


# ============================================================================
# TestParseIlluminaSamplesheet
# ============================================================================


class TestParseIlluminaSamplesheet:
    """Tests for parse_illumina_samplesheet()."""

    def test_standard_parsing(self, illumina_samplesheet_lines):
        result = gc.parse_illumina_samplesheet(illumina_samplesheet_lines)
        assert len(result) == 3
        assert result[0]["sample_name"] == "SampleA"
        assert result[0]["lane"] == "1"
        assert result[0]["sample_project"] == "TestProject"

    def test_bclconvert_section(self):
        lines = [
            "[BCLConvert_Data]\n",
            "Lane,Sample_ID,index,index2,Sample_Project\n",
            "1,MySample,ATCACG,TTAGGC,Proj\n",
        ]
        result = gc.parse_illumina_samplesheet(lines)
        assert len(result) == 1
        # Falls back to Sample_ID when Sample_Name not present
        assert result[0]["sample_name"] == "MySample"

    def test_empty_rows_skipped(self):
        lines = [
            "[Data]\n",
            "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project\n",
            "1,S1,SampleA,ATCACG,TTAGGC,Proj\n",
            ",,,,\n",
            "2,S2,SampleB,CGATGT,AAGCTA,Proj\n",
        ]
        result = gc.parse_illumina_samplesheet(lines)
        assert len(result) == 2

    def test_sample_id_fallback(self):
        lines = [
            "[Data]\n",
            "Lane,Sample_ID,index,index2,Sample_Project\n",
            "1,FallbackName,ATCACG,TTAGGC,Proj\n",
        ]
        result = gc.parse_illumina_samplesheet(lines)
        assert result[0]["sample_name"] == "FallbackName"

    def test_no_data_section(self):
        lines = ["[Header]\n", "Key,Value\n"]
        result = gc.parse_illumina_samplesheet(lines)
        assert result == []

    def test_data_section_at_end(self):
        lines = [
            "[Header]\n",
            "Key,Value\n",
            "[Data]\n",
        ]
        result = gc.parse_illumina_samplesheet(lines)
        assert result == []


# ============================================================================
# TestParseMinimalSamplesheet
# ============================================================================


class TestParseMinimalSamplesheet:
    """Tests for parse_minimal_samplesheet()."""

    def test_5col_standard(self, minimal_5col_lines):
        result = gc.parse_minimal_samplesheet(minimal_5col_lines)
        assert len(result) == 3
        assert result[0]["lane"] == "1"
        assert result[0]["sample_name"] == "SampleA"
        assert result[0]["sample_project"] == "TestProject"

    def test_4col_no_project(self, minimal_4col_lines):
        result = gc.parse_minimal_samplesheet(minimal_4col_lines)
        assert len(result) == 2
        assert result[0]["sample_project"] == ""

    def test_header_skip(self):
        lines = [
            "Lane,Sample_Name,index_i7,index_i5,Project\n",
            "1,SampleA,ATCACG,TTAGGC,Proj\n",
        ]
        result = gc.parse_minimal_samplesheet(lines)
        assert len(result) == 1

    def test_comment_skip(self):
        lines = [
            "# This is a comment\n",
            "1,SampleA,ATCACG,TTAGGC,Proj\n",
        ]
        result = gc.parse_minimal_samplesheet(lines)
        assert len(result) == 1

    def test_alt_column_order(self):
        """Sample_Name first, Lane second."""
        lines = ["SampleA,1,ATCACG,TTAGGC,Proj\n"]
        result = gc.parse_minimal_samplesheet(lines)
        assert len(result) == 1
        assert result[0]["sample_name"] == "SampleA"
        assert result[0]["lane"] == "1"

    def test_too_few_columns(self):
        lines = ["1,SampleA,ATCACG\n"]
        result = gc.parse_minimal_samplesheet(lines)
        assert len(result) == 0

    def test_empty_lines_skipped(self):
        lines = ["\n", "  \n", "1,SampleA,ATCACG,TTAGGC,Proj\n"]
        result = gc.parse_minimal_samplesheet(lines)
        assert len(result) == 1


# ============================================================================
# TestParseSamplesheet
# ============================================================================


class TestParseSamplesheet:
    """Tests for parse_samplesheet() file-level wrapper."""

    def test_illumina_file(self, tmp_path):
        ss = tmp_path / "SampleSheet.csv"
        ss.write_text(
            "[Data]\nLane,Sample_ID,Sample_Name,index,index2,Sample_Project\n"
            "1,S1,SampleA,ATCACG,TTAGGC,Proj\n",
            encoding="utf-8",
        )
        result = gc.parse_samplesheet(ss)
        assert len(result) == 1
        assert result[0]["sample_name"] == "SampleA"

    def test_empty_file(self, tmp_path):
        ss = tmp_path / "SampleSheet.csv"
        ss.write_text("", encoding="utf-8")
        result = gc.parse_samplesheet(ss)
        assert result == []

    def test_bom_handling(self, tmp_path):
        ss = tmp_path / "SampleSheet.csv"
        # Write BOM-prefixed content; utf-8-sig encoding adds BOM automatically
        content = (
            "[Data]\nLane,Sample_ID,Sample_Name,index,index2,Sample_Project\n"
            "1,S1,SampleA,ATCACG,TTAGGC,Proj\n"
        )
        ss.write_text(content, encoding="utf-8-sig")
        result = gc.parse_samplesheet(ss)
        assert len(result) == 1


# ============================================================================
# TestFindSamplesheet
# ============================================================================


class TestFindSamplesheet:
    """Tests for find_samplesheet()."""

    def test_found_in_dir(self, tmp_path):
        (tmp_path / "SampleSheet.csv").write_text("data", encoding="utf-8")
        result = gc.find_samplesheet(tmp_path)
        assert result is not None
        assert result.name == "SampleSheet.csv"

    def test_found_in_parent(self, tmp_path):
        subdir = tmp_path / "fastqs"
        subdir.mkdir()
        (tmp_path / "SampleSheet.csv").write_text("data", encoding="utf-8")
        result = gc.find_samplesheet(subdir)
        assert result is not None

    def test_lowercase(self, tmp_path):
        (tmp_path / "samplesheet.csv").write_text("data", encoding="utf-8")
        result = gc.find_samplesheet(tmp_path)
        assert result is not None
        assert result.name.lower() == "samplesheet.csv"

    def test_not_found(self, tmp_path):
        result = gc.find_samplesheet(tmp_path)
        assert result is None

    def test_preference_order(self, tmp_path):
        """SampleSheet.csv preferred over samplesheet.csv."""
        (tmp_path / "SampleSheet.csv").write_text("preferred", encoding="utf-8")
        (tmp_path / "samplesheet.csv").write_text("fallback", encoding="utf-8")
        result = gc.find_samplesheet(tmp_path)
        assert result is not None
        assert result.name == "SampleSheet.csv"


# ============================================================================
# TestDiscoverFastqFiles
# ============================================================================


class TestDiscoverFastqFiles:
    """Tests for discover_fastq_files()."""

    def test_pair_discovery(self, fastq_dir):
        pairs = gc.discover_fastq_files(fastq_dir)
        assert len(pairs) == 3

    def test_orphan_skip(self, fastq_dir_with_orphan):
        pairs = gc.discover_fastq_files(fastq_dir_with_orphan)
        assert len(pairs) == 3
        basenames = [p["basename"] for p in pairs]
        assert "OrphanSample_S9_L001" not in basenames

    def test_lane_extraction(self, fastq_dir):
        pairs = gc.discover_fastq_files(fastq_dir)
        lanes = sorted(set(p["lane_str"] for p in pairs))
        assert "L001" in lanes
        assert "L002" in lanes

    def test_sample_extraction(self, fastq_dir):
        pairs = gc.discover_fastq_files(fastq_dir)
        samples = sorted(set(p["sample_name"] for p in pairs))
        assert samples == ["SampleA", "SampleB"]

    def test_nonexistent_dir(self, tmp_path):
        with pytest.raises(SystemExit):
            gc.discover_fastq_files(tmp_path / "nonexistent")

    def test_empty_dir(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(SystemExit):
            gc.discover_fastq_files(empty)

    def test_basename_format(self, fastq_dir):
        pairs = gc.discover_fastq_files(fastq_dir)
        for p in pairs:
            assert "_S" in p["basename"]
            assert "_L" in p["basename"]


# ============================================================================
# TestInferProject
# ============================================================================


class TestInferProject:
    """Tests for infer_project()."""

    def test_explicit_arg(self, tmp_path):
        result = gc.infer_project(tmp_path, "MyProject")
        assert result == "MyProject"

    def test_a_pattern(self, tmp_path):
        d = tmp_path / "delivery_A5297_exomes"
        d.mkdir()
        result = gc.infer_project(d, None)
        assert result == "A5297"

    def test_fallback_dirname(self, tmp_path):
        d = tmp_path / "no_pattern"
        d.mkdir()
        result = gc.infer_project(d, None)
        assert result == "no_pattern"

    def test_override(self, tmp_path):
        d = tmp_path / "delivery_A5297"
        d.mkdir()
        result = gc.infer_project(d, "CustomProject")
        assert result == "CustomProject"


# ============================================================================
# TestBuildSamplesTable
# ============================================================================


class TestBuildSamplesTable:
    """Tests for build_samples_table()."""

    def test_merge(self, fastq_pairs):
        ss_entries = [
            {
                "sample_name": "SampleA",
                "lane": "1",
                "index_i7": "",
                "index_i5": "",
                "sample_project": "SSProject",
            },
        ]
        df = gc.build_samples_table(fastq_pairs, ss_entries, "DefaultProj")
        assert len(df) == 3
        assert "fastq_files_basename" in df.columns
        assert "project_sample" in df.columns

    def test_no_samplesheet(self, fastq_pairs):
        df = gc.build_samples_table(fastq_pairs, None, "TestProject")
        assert len(df) == 3
        assert all(df["mdc_project"] == "TestProject")

    def test_project_override(self, fastq_pairs):
        ss_entries = [
            {
                "sample_name": "SampleA",
                "lane": "1",
                "index_i7": "",
                "index_i5": "",
                "sample_project": "SSProject",
            },
        ]
        df = gc.build_samples_table(fastq_pairs, ss_entries, "FallbackProj")
        # SampleA rows should use SSProject, SampleB should use FallbackProj
        sample_a = df[df["project_sample"] == "SampleA"]
        sample_b = df[df["project_sample"] == "SampleB"]
        assert all(sample_a["mdc_project"] == "SSProject")
        assert all(sample_b["mdc_project"] == "FallbackProj")

    def test_sorting(self, fastq_pairs):
        df = gc.build_samples_table(fastq_pairs, None, "Proj")
        # Should be sorted by project_sample, lane, basename
        samples = df["project_sample"].tolist()
        assert samples == sorted(samples)

    def test_warnings_unmatched(self, fastq_pairs, capsys):
        ss_entries = [
            {
                "sample_name": "NoSuchSample",
                "lane": "1",
                "index_i7": "",
                "index_i5": "",
                "sample_project": "Proj",
            },
        ]
        gc.build_samples_table(fastq_pairs, ss_entries, "Proj")
        captured = capsys.readouterr()
        assert "NoSuchSample" in captured.err

    def test_partial_match(self, fastq_pairs):
        """SampleSheet name is a prefix of the FASTQ sample name."""
        ss_entries = [
            {
                "sample_name": "Sample",
                "lane": "1",
                "index_i7": "",
                "index_i5": "",
                "sample_project": "PartialProj",
            },
        ]
        df = gc.build_samples_table(fastq_pairs, ss_entries, "FallbackProj")
        # Should match via partial matching (fastq_sample.startswith(ss_name))
        assert any(df["mdc_project"] == "PartialProj")


# ============================================================================
# TestFindGenomeFastas
# ============================================================================


class TestFindGenomeFastas:
    """Tests for _find_genome_fastas()."""

    def test_with_companions(self, ref_dir):
        results = gc._find_genome_fastas(ref_dir)
        assert len(results) == 1
        assert results[0]["has_bwa"]
        assert results[0]["has_fai"]

    def test_nonexistent_dir(self, tmp_path):
        results = gc._find_genome_fastas(tmp_path / "nope")
        assert results == []

    def test_hidden_files_excluded(self, tmp_path):
        d = tmp_path / "ref"
        d.mkdir()
        (d / ".hidden.fna").touch()
        results = gc._find_genome_fastas(d)
        assert len(results) == 0

    def test_gz_only_reference(self, tmp_path):
        d = tmp_path / "ref"
        d.mkdir()
        gz = d / "genome.fna.gz"
        gz.touch()
        # Add BWA index for the .gz
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            (d / f"genome.fna.gz{ext}").touch()
        results = gc._find_genome_fastas(d)
        assert len(results) == 1
        assert results[0]["has_gz"]


# ============================================================================
# TestFindKnownSitesVcfs
# ============================================================================


class TestFindKnownSitesVcfs:
    """Tests for _find_known_sites_vcfs()."""

    def test_matching(self, known_sites_dir):
        results = gc._find_known_sites_vcfs(known_sites_dir)
        assert len(results) >= 2

    def test_tbi_excluded(self, known_sites_dir):
        results = gc._find_known_sites_vcfs(known_sites_dir)
        for r in results:
            assert not r.endswith(".tbi")

    def test_nonexistent(self, tmp_path):
        results = gc._find_known_sites_vcfs(tmp_path / "nope")
        assert results == []

    def test_missing_index_warning(self, tmp_path, capsys):
        d = tmp_path / "ks"
        d.mkdir()
        (d / "dbsnp_test.vcf.gz").touch()
        # No .tbi file
        gc._find_known_sites_vcfs(d)
        captured = capsys.readouterr()
        assert "Warning" in captured.err

    def test_dedup(self, tmp_path):
        d = tmp_path / "ks"
        d.mkdir()
        # Same file matched by multiple patterns
        (d / "dbsnp_known_indels.vcf.gz").touch()
        (d / "dbsnp_known_indels.vcf.gz.tbi").touch()
        results = gc._find_known_sites_vcfs(d)
        # Should be deduplicated
        assert len(results) == len(set(results))

    def test_gnomad_excluded(self, tmp_path):
        """gnomAD and panel-of-normals VCFs should not be found by default patterns."""
        d = tmp_path / "ks"
        d.mkdir()
        # BQSR-relevant
        (d / "dbsnp_138.vcf.gz").touch()
        (d / "dbsnp_138.vcf.gz.tbi").touch()
        # NOT BQSR-relevant
        (d / "af-only-gnomad.hg38.vcf.gz").touch()
        (d / "af-only-gnomad.hg38.vcf.gz.tbi").touch()
        (d / "1000g_pon.hg38.vcf.gz").touch()
        (d / "1000g_pon.hg38.vcf.gz.tbi").touch()
        results = gc._find_known_sites_vcfs(d)
        names = [Path(r).name for r in results]
        assert "dbsnp_138.vcf.gz" in names
        assert "af-only-gnomad.hg38.vcf.gz" not in names
        assert "1000g_pon.hg38.vcf.gz" not in names

    def test_custom_patterns(self, tmp_path):
        """Custom patterns parameter should override defaults."""
        d = tmp_path / "ks"
        d.mkdir()
        (d / "af-only-gnomad.hg38.vcf.gz").touch()
        (d / "af-only-gnomad.hg38.vcf.gz.tbi").touch()
        results = gc._find_known_sites_vcfs(d, patterns=("*gnomad*.vcf.gz",))
        assert len(results) == 1


# ============================================================================
# TestDiscoverReferenceData
# ============================================================================


class TestDiscoverReferenceData:
    """Tests for discover_reference_data()."""

    def test_explicit_dir(self, ref_dir, known_sites_dir):
        # Put known-sites in sibling dir to ref_dir
        result = gc.discover_reference_data(ref_dir)
        assert result["genome"]
        assert result["build"] == "GRCh38"

    def test_nothing_found(self, tmp_path):
        result = gc.discover_reference_data(tmp_path / "empty")
        assert result["genome"] == ""

    def test_grch37_detection(self, tmp_path):
        d = tmp_path / "ref"
        d.mkdir()
        (d / "hg19.fna").touch()
        (d / "hg19.fna.fai").touch()
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            (d / f"hg19.fna{ext}").touch()
        result = gc.discover_reference_data(d)
        assert result["build"] == "GRCh37"

    def test_bwa_preference(self, tmp_path):
        d = tmp_path / "ref"
        d.mkdir()
        # Genome without BWA index
        (d / "no_bwa.fna").touch()
        (d / "no_bwa.fna.fai").touch()
        # Genome with BWA index
        (d / "with_bwa.fna").touch()
        (d / "with_bwa.fna.fai").touch()
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            (d / f"with_bwa.fna{ext}").touch()
        result = gc.discover_reference_data(d)
        assert "with_bwa" in result["genome"]

    def test_genome_gz_without_bwa_index(self, tmp_path):
        """genome_gz should be empty when .gz exists but has no BWA index."""
        d = tmp_path / "ref"
        d.mkdir()
        (d / "ref.fna").touch()
        (d / "ref.fna.fai").touch()
        (d / "ref.fna.gz").touch()  # .gz exists but no BWA index for it
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            (d / f"ref.fna{ext}").touch()  # BWA index for uncompressed only
        result = gc.discover_reference_data(d)
        assert result["genome"]  # uncompressed found
        assert result["genome_gz"] == ""  # .gz has no BWA index
        assert any("Warning" in line for line in result["search_log"])

    def test_genome_gz_with_bwa_index(self, tmp_path):
        """genome_gz should be set when .gz has a BWA index."""
        d = tmp_path / "ref"
        d.mkdir()
        (d / "ref.fna").touch()
        (d / "ref.fna.fai").touch()
        (d / "ref.fna.gz").touch()
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            (d / f"ref.fna{ext}").touch()
            (d / f"ref.fna.gz{ext}").touch()  # BWA index for .gz too
        result = gc.discover_reference_data(d)
        assert result["genome_gz"]  # .gz has BWA index

    def test_paths_are_resolved(self, tmp_path):
        """Discovered paths should be canonical (no ../ segments)."""
        # Create ref dir one level deep, then search from parent with ..
        d = tmp_path / "sub" / "ref"
        d.mkdir(parents=True)
        (d / "ref.fna").touch()
        (d / "ref.fna.fai").touch()
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            (d / f"ref.fna{ext}").touch()
        result = gc.discover_reference_data(d)
        assert ".." not in result["genome"]


# ============================================================================
# TestBuildConfigYaml
# ============================================================================


class TestBuildConfigYaml:
    """Tests for _build_config_yaml()."""

    def test_required_sections(self):
        ref_data = {
            "genome": "/path/to/ref.fna",
            "genome_gz": "",
            "build": "GRCh38",
            "known_sites": ["/path/to/dbsnp.vcf.gz"],
        }
        result = gc._build_config_yaml(ref_data, "/fastqs", "/output", "samples.tsv")
        assert "ref:" in result
        assert "paths:" in result
        assert "trimming:" in result
        assert "params:" in result

    def test_edit_me_placeholders(self):
        ref_data = {"genome": "", "genome_gz": "", "build": "GRCh38", "known_sites": []}
        result = gc._build_config_yaml(ref_data, "/fastqs", "/output", "samples.tsv")
        assert "EDIT_ME" in result

    def test_trimming_section(self):
        ref_data = {
            "genome": "/ref.fna",
            "genome_gz": "",
            "build": "GRCh38",
            "known_sites": ["/ks.vcf.gz"],
        }
        result = gc._build_config_yaml(ref_data, "/fastqs", "/output", "samples.tsv")
        assert "trimming:" in result
        assert "enabled: false" in result

    def test_edit_me_not_resolved(self):
        """EDIT_ME placeholders should not be resolved by _resolve_path."""
        ref_data = {"genome": "", "genome_gz": "", "build": "GRCh38", "known_sites": []}
        result = gc._build_config_yaml(ref_data, "/fastqs", "/output", "samples.tsv")
        assert "EDIT_ME" in result


# ============================================================================
# TestWriteSamplesTsv
# ============================================================================


class TestWriteSamplesTsv:
    """Tests for write_samples_tsv()."""

    def test_dry_run(self, tmp_path, capsys):
        df = pd.DataFrame({"col": ["val"]})
        output = tmp_path / "samples.tsv"
        gc.write_samples_tsv(df, output, dry_run=True)
        assert not output.exists()
        captured = capsys.readouterr()
        assert "dry-run" in captured.out

    def test_write(self, tmp_path):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["S1"],
                "lane": ["L001"],
                "project_sample": ["Sample"],
                "mdc_project": ["Proj"],
            }
        )
        output = tmp_path / "samples.tsv"
        gc.write_samples_tsv(df, output, dry_run=False, force=True)
        assert output.exists()
        content = output.read_text()
        assert "S1" in content

    def test_force_overwrite(self, tmp_path):
        df = pd.DataFrame({"col": ["new"]})
        output = tmp_path / "samples.tsv"
        output.write_text("old content")
        gc.write_samples_tsv(df, output, dry_run=False, force=True)
        content = output.read_text()
        assert "new" in content

    def test_no_force_prompt(self, tmp_path):
        df = pd.DataFrame({"col": ["val"]})
        output = tmp_path / "samples.tsv"
        output.write_text("old")
        with patch("builtins.input", return_value="n"):
            gc.write_samples_tsv(df, output, dry_run=False, force=False)
        assert output.read_text() == "old"
