"""Dry-run tests for the Snakemake pipeline.

These tests run `snakemake` with test data to verify
that the workflow loads correctly and rules are defined.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

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


def _write_config(tmp_path: Path, name: str = "config.yaml", **overrides) -> Path:
    """Write a temporary config file derived from TEST_CONFIG with top-level overrides."""
    cfg = yaml.safe_load(TEST_CONFIG.read_text())
    for key, val in overrides.items():
        if isinstance(val, dict) and isinstance(cfg.get(key), dict):
            cfg[key] = {**cfg[key], **val}
        else:
            cfg[key] = val
    out = tmp_path / name
    out.write_text(yaml.safe_dump(cfg))
    return out


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
        cfg = _write_config(tmp_path, "cfg_ok.yaml", ref={"target_intervals": str(bed)})
        result = _run_snakemake("-p", config=cfg)
        assert result.returncode == 0, result.stderr
        assert f'-L "{bed}" --interval-padding 100' in result.stdout
        assert "ApplyBQSR" in result.stdout
        for block in result.stdout.split("rule "):
            if block.startswith("apply_bqsr:"):
                assert "-L " not in block

        missing = tmp_path / "missing.bed"
        bad_cfg = _write_config(tmp_path, "cfg_bad.yaml", ref={"target_intervals": str(missing)})
        fail_res = _run_snakemake(config=bad_cfg)
        assert fail_res.returncode != 0
        assert "MissingInputException" in (fail_res.stdout + fail_res.stderr)

    def test_trimming_with_subfolder_and_defaults(self, tmp_path, fastq_dir):
        """Dry-run with trimming.enabled=true and subfolder set succeeds with ftl=0 qtrim=f and no quantize=."""
        run1 = fastq_dir / "run1"
        run1.mkdir()
        for fq in fastq_dir.glob("*.fastq.gz"):
            shutil.copy(fq, run1 / fq.name)
        samples_tsv = tmp_path / "samples.tsv"
        samples_tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\tsubfolder\n"
            "SampleA_S1_L001\tL001\tSampleA\tProj\trun1\n"
        )
        cfg = _write_config(
            tmp_path,
            trimming={"enabled": True},
            paths={
                "samples": str(samples_tsv),
                "fastq_folder": str(tmp_path / "fastqs"),
                "output_folder": str(tmp_path / "out"),
            },
        )
        result = _run_snakemake("-p", config=cfg)
        assert result.returncode == 0, result.stderr
        assert "-Xmx3400m" in result.stdout
        assert "ftl=0" in result.stdout
        assert "qtrim=f" in result.stdout
        assert "quantize=" not in result.stdout

    def test_bwa_map_default_extra(self):
        """Dry-run bwa_map command contains -K 100000000 -Y."""
        result = _run_snakemake("-p")
        assert result.returncode == 0, result.stderr
        assert "-K 100000000 -Y" in result.stdout

    def test_dedup_and_validate_sam(self, tmp_path):
        """Dry-run includes OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500, tracks dedup_bai, and schedules picard_validate_sam."""
        result = _run_snakemake("-p")
        assert result.returncode == 0, result.stderr
        assert "--CREATE_INDEX true" in result.stdout
        assert "--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500" in result.stdout
        assert ".merged.dedup.bai" in result.stdout
        assert "rule picard_validate_sam:" in result.stdout
        assert "--IGNORE_WARNINGS true" in result.stdout
        assert ".pdf" not in result.stdout

        legacy_cfg = _write_config(
            tmp_path,
            params={
                "gatk": {"MarkDuplicates": "--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"}
            },
        )
        legacy_res = _run_snakemake("-p", config=legacy_cfg)
        assert legacy_res.returncode == 0, legacy_res.stderr
        assert "--VALIDATION_STRINGENCY SILENT" in legacy_res.stdout
        for block in legacy_res.stdout.split("rule "):
            if block.startswith("deduplicate_bam_files:"):
                assert block.count("--CREATE_INDEX") == 1

    def test_default_profile_runtimes(self, tmp_path):
        """Dry-run with --workflow-profile profiles/default and profiles/bih resolves runtimes <=480 and slurm_account."""
        bih_dir = tmp_path / "bih"
        bih_dir.mkdir()
        bih_lines = [
            line
            for line in (REPO_ROOT / "profiles" / "bih" / "config.yaml").read_text().splitlines()
            if not line.strip().startswith("executor:")
        ]
        (bih_dir / "config.yaml").write_text("\n".join(bih_lines) + "\n")
        result = _run_snakemake(
            "--workflow-profile",
            str(REPO_ROOT / "profiles" / "default"),
            "--profile",
            str(bih_dir),
        )
        assert result.returncode == 0, result.stderr
        res_lines = [
            line.strip()
            for line in result.stdout.splitlines()
            if line.strip().startswith("resources:")
        ]
        assert res_lines, "Expected resources: lines in dry-run output"
        for line in res_lines:
            assert "slurm_account=hpc-ag-scholl" in line, f"Missing slurm_account in {line}"
            assert "slurm_partition=standard" in line, f"Missing slurm_partition in {line}"
            m = re.search(r"\bruntime=(\d+)\b", line)
            assert m is not None, f"Missing runtime in {line}"
            assert int(m.group(1)) <= 480, f"runtime {m.group(1)} > 480 in {line}"

    def test_hs_metrics_enabled(self, tmp_path):
        """Dry-run with qc.hs_metrics.enabled=true schedules picard_collect_hs_metrics before multiqc."""
        baits = tmp_path / "baits.interval_list"
        targets = tmp_path / "targets.interval_list"
        baits.write_text("@HD\tVN:1.6\n")
        targets.write_text("@HD\tVN:1.6\n")
        samples_tsv = tmp_path / "samples.tsv"
        samples_tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\tkit\n"
            "SampleA_S1_L001\tL001\tSampleA\tTestProject\tS33266340\n"
            "SampleA_S1_L002\tL002\tSampleA\tTestProject\tS33266340\n"
            "SampleB_S2_L001\tL001\tSampleB\tTestProject\tS33266340\n"
        )
        cfg = _write_config(
            tmp_path,
            paths={"samples": str(samples_tsv)},
            qc={
                "enabled": True,
                "hs_metrics": {
                    "enabled": True,
                    "kits": {
                        "S33266340": {
                            "bait_intervals": str(baits),
                            "target_intervals": str(targets),
                        }
                    },
                },
            },
        )
        result = _run_snakemake("-p", config=cfg)
        assert result.returncode == 0, result.stderr
        assert "rule picard_collect_hs_metrics:" in result.stdout
        assert "SampleA.hs_metrics.txt" in result.stdout
        assert "SampleB.hs_metrics.txt" in result.stdout
        assert f'--BAIT_INTERVALS "{baits}"' in result.stdout
        assert f'--TARGET_INTERVALS "{targets}"' in result.stdout

    def test_external_staging_directory(self, tmp_path):
        """Invoking snakemake -s workflow/Snakefile --configfile ... -n from an empty temp dir succeeds."""
        stage_dir = tmp_path / "stage"
        stage_dir.mkdir()
        cfg = _write_config(
            tmp_path,
            "stage_config.yaml",
            ref={
                "genome": str(REPO_ROOT / "tests" / "data" / "ref" / "dummy.fna"),
                "known_sites": [
                    str(REPO_ROOT / "tests" / "data" / "known_sites" / "dbsnp_mini.vcf.gz"),
                    str(
                        REPO_ROOT
                        / "tests"
                        / "data"
                        / "known_sites"
                        / "Mills_and_1000G_mini.vcf.gz"
                    ),
                ],
            },
            paths={
                "samples": str(REPO_ROOT / "tests" / "data" / "config" / "samples.tsv"),
                "fastq_folder": str(REPO_ROOT / "tests" / "data" / "fastqs"),
                "output_folder": str(stage_dir / "out"),
            },
        )
        result = subprocess.run(
            ["snakemake", "-s", str(SNAKEFILE), "--configfile", str(cfg), "-n"],
            cwd=stage_dir,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, (
            f"External stage dry-run failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        assert "SampleA.merged.dedup.bqsr.bam" in result.stdout

        spooled_launcher = tmp_path / "slurm_script"
        shutil.copy(REPO_ROOT / "scripts" / "run_snakemake.sh", spooled_launcher)
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        args_file = tmp_path / "snakemake_args.txt"
        conda_file = tmp_path / "conda_called.txt"
        (bin_dir / "snakemake").write_text(f'#!/bin/bash\nprintf "%s\\n" "$@" > "{args_file}"\n')
        (bin_dir / "snakemake").chmod(0o755)
        (bin_dir / "conda").write_text(f'#!/bin/bash\necho "$@" > "{conda_file}"\n')
        (bin_dir / "conda").chmod(0o755)
        launcher_res = subprocess.run(
            ["bash", str(spooled_launcher), str(SNAKEFILE), str(cfg), "-n"],
            cwd=stage_dir,
            env={**os.environ, "PATH": f"{bin_dir}:/usr/bin:/bin"},
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert launcher_res.returncode == 0, launcher_res.stderr
        assert not conda_file.exists()
        logged_args = args_file.read_text().splitlines()
        assert str(REPO_ROOT / "profiles" / "default") in logged_args
