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
