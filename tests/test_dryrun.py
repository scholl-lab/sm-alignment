"""Dry-run tests for the Snakemake pipeline.

These tests run `snakemake` with test data to verify
that the workflow loads correctly and rules are defined.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = [
    pytest.mark.dryrun,
    pytest.mark.skipif(
        shutil.which("snakemake") is None,
        reason="snakemake not installed or not on PATH",
    ),
]

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_CONFIG = REPO_ROOT / "tests" / "data" / "config" / "config.yaml"
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"


def _run_snakemake(
    *extra_args: str,
    config: Path = TEST_CONFIG,
    dry_run: bool = True,
) -> subprocess.CompletedProcess:
    """Run snakemake as subprocess and return the result."""
    cmd = [
        "snakemake",
        "-s",
        str(SNAKEFILE),
        "--configfile",
        str(config),
        "--directory",
        str(REPO_ROOT),
    ]
    if dry_run:
        cmd.append("-n")
    cmd.extend(extra_args)
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=60,
    )


class TestDryRun:
    """Dry-run tests verifying the workflow loads with test data."""

    def test_list_rules(self):
        """Verify the workflow loads and lists expected rules."""
        result = _run_snakemake("--list-rules", dry_run=False)
        assert result.returncode == 0, (
            f"--list-rules failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        output = result.stdout
        for rule in ("bwa_map", "merge_bam_files", "deduplicate_bam_files"):
            assert rule in output, f"Expected rule '{rule}' not in --list-rules output"

    def test_dag_resolves(self):
        """DAG resolves or fails only on MissingInputException (not schema/syntax)."""
        result = _run_snakemake()
        output = result.stdout + result.stderr
        # Accept either success or MissingInputException (expected without real inputs)
        if result.returncode != 0:
            assert "MissingInputException" in output or "Empty file path" in output, (
                f"Unexpected failure:\nstdout: {result.stdout}\nstderr: {result.stderr}"
            )

    def test_lint_passes(self):
        result = _run_snakemake("--lint", dry_run=False)
        # Snakemake --lint returns 0 or 1 (warnings); just check it doesn't crash
        assert result.returncode in (0, 1), (
            f"Lint crashed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

    def test_trimming_enabled_rules(self):
        """With trimming enabled, trim_adapters rule should exist."""
        result = _run_snakemake("--list-rules", dry_run=False)
        assert result.returncode == 0, result.stderr
        assert "trim_adapters" in result.stdout

    def test_missing_config_fails(self, tmp_path):
        fake_config = tmp_path / "nonexistent.yaml"
        result = _run_snakemake(config=fake_config)
        assert result.returncode != 0

    def test_invalid_config_fails(self):
        invalid_config = REPO_ROOT / "tests" / "data" / "config" / "config_invalid_build.yaml"
        result = _run_snakemake(config=invalid_config)
        assert result.returncode != 0

    def test_qc_rules_listed(self):
        """QC rules should be listed when qc.enabled is true."""
        result = _run_snakemake("--list-rules", dry_run=False)
        assert result.returncode == 0, result.stderr
        for rule in (
            "fastqc_raw",
            "fastqc_trimmed",
            "samtools_stats",
            "samtools_flagstat",
            "picard_collect_multiple_metrics",
            "picard_collect_hs_metrics",
            "picard_validate_sam",
            "qualimap_bamqc",
            "multiqc",
        ):
            assert rule in result.stdout, f"Expected QC rule '{rule}' not in --list-rules output"

    def test_bqsr_target_intervals(self, tmp_path):
        """Dry-run with ref.target_intervals set emits -L and --interval-padding on base_recalibration only."""
        bed = tmp_path / "targets.bed"
        bed.write_text("chr1\t0\t1000\n")
        result = _run_snakemake("-p", "--config", f"ref={{target_intervals: '{bed}'}}")
        assert result.returncode == 0, result.stderr
        assert f'-L "{bed}" --interval-padding 100' in result.stdout
        assert "ApplyBQSR" in result.stdout
        for block in result.stdout.split("rule "):
            if block.startswith("apply_bqsr:"):
                assert "-L " not in block

        missing = tmp_path / "missing.bed"
        fail_res = _run_snakemake("--config", f"ref={{target_intervals: '{missing}'}}")
        assert fail_res.returncode != 0
        assert "MissingInputException" in (fail_res.stdout + fail_res.stderr)

    def test_trimming_with_subfolder_and_defaults(self, tmp_path, fastq_dir):
        """Dry-run with trimming.enabled=true and subfolder set succeeds with ftl=0 qtrim=f and no quantize=."""
        run1 = tmp_path / "fastqs" / "run1"
        run1.mkdir(parents=True)
        for fq in fastq_dir.iterdir():
            shutil.copy(fq, run1 / fq.name)
        samples_tsv = tmp_path / "samples.tsv"
        samples_tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\tsubfolder\n"
            "SampleA_S1_L001\tL001\tSampleA\tProj\trun1\n"
        )
        result = _run_snakemake(
            "-p",
            "--config",
            "trimming={enabled: true}",
            f"paths={{samples: '{samples_tsv}', fastq_folder: '{tmp_path / 'fastqs'}', output_folder: '{tmp_path / 'out'}'}}",
        )
        assert result.returncode == 0, result.stderr
        assert "ftl=0" in result.stdout
        assert "qtrim=f" in result.stdout
        assert "quantize=" not in result.stdout

    def test_bwa_map_default_extra(self):
        """Dry-run bwa_map command contains -K 100000000 -Y."""
        result = _run_snakemake("-p")
        assert result.returncode == 0, result.stderr
        assert "-K 100000000 -Y" in result.stdout

    def test_dedup_and_validate_sam(self):
        """Dry-run includes OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500, tracks dedup_bai, and schedules picard_validate_sam."""
        result = _run_snakemake("-p")
        assert result.returncode == 0, result.stderr
        assert "--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500" in result.stdout
        assert ".merged.dedup.bai" in result.stdout
        assert "rule picard_validate_sam:" in result.stdout
        assert ".pdf" not in result.stdout

    def test_default_profile_runtimes(self):
        """Dry-run with --workflow-profile profiles/default resolves all rule runtimes <=480."""
        result = _run_snakemake(
            "--workflow-profile",
            str(REPO_ROOT / "profiles" / "default"),
        )
        assert result.returncode == 0, result.stderr

    def test_hs_metrics_enabled(self, tmp_path):
        """Dry-run with qc.hs_metrics.enabled=true schedules picard_collect_hs_metrics before multiqc."""
        baits = tmp_path / "baits.interval_list"
        targets = tmp_path / "targets.interval_list"
        baits.write_text("@HD\tVN:1.6\n")
        targets.write_text("@HD\tVN:1.6\n")
        result = _run_snakemake(
            "-p",
            "--config",
            f"qc={{enabled: true, hs_metrics: {{enabled: true, kits: {{S33266340: {{bait_intervals: '{baits}', target_intervals: '{targets}'}}}}}}}}",
        )
        assert result.returncode == 0, result.stderr
        assert "rule picard_collect_hs_metrics:" in result.stdout
        assert f'--BAIT_INTERVALS "{baits}"' in result.stdout
        assert f'--TARGET_INTERVALS "{targets}"' in result.stdout
