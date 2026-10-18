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

pytestmark = [pytest.mark.dryrun]

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

        cli_only_res = subprocess.run(
            [
                "snakemake",
                "-s",
                str(SNAKEFILE),
                "--list-rules",
                "--config",
                f'paths={{"samples": "{REPO_ROOT / "tests/data/config/samples.tsv"}"}}',
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert cli_only_res.returncode == 0, (
            f"CLI-only --config failed:\nstdout: {cli_only_res.stdout}\nstderr: {cli_only_res.stderr}"
        )
        assert "bwa_map" in cli_only_res.stdout

    def test_umi_disabled_dag_and_commands_identical(self, tmp_path):
        """With umi.enabled=false, DAG and rendered shell commands are identical to config without umi."""
        cfg_without_umi = _write_config(tmp_path, "cfg_no_umi.yaml", trimming={"enabled": True})
        raw = yaml.safe_load(cfg_without_umi.read_text())
        raw.pop("umi", None)
        cfg_without_umi.write_text(yaml.safe_dump(raw))

        cfg_umi_false = _write_config(
            tmp_path,
            "cfg_umi_false.yaml",
            trimming={"enabled": True},
            umi={"enabled": False, "preset": None, "read_structures": None},
        )

        dag_no_umi = _run_snakemake("--dag", config=cfg_without_umi, dry_run=False)
        dag_umi_false = _run_snakemake("--dag", config=cfg_umi_false, dry_run=False)
        assert dag_no_umi.returncode == 0, dag_no_umi.stderr
        assert dag_umi_false.returncode == 0, dag_umi_false.stderr
        assert "umi_extract" not in dag_umi_false.stdout
        assert dag_no_umi.stdout == dag_umi_false.stdout

        dry_no_umi = _run_snakemake("-p", "--quiet", "progress", config=cfg_without_umi)
        dry_umi_false = _run_snakemake("-p", "--quiet", "progress", config=cfg_umi_false)
        assert dry_no_umi.returncode == 0, dry_no_umi.stderr
        assert dry_umi_false.returncode == 0, dry_umi_false.stderr
        assert "rule umi_extract:" not in dry_umi_false.stdout
        assert "--BARCODE_TAG" not in dry_umi_false.stdout
        assert "--DUPLEX_UMI" not in dry_umi_false.stdout
        for block in dry_umi_false.stdout.split("rule "):
            if block.startswith("bwa_map:"):
                assert " -C" not in block

    def test_umi_enabled_dryrun_commands(self, tmp_path):
        """Dry-run wiring test: umi.enabled=true schedules umi_extract -> trim_adapters -> bwa_map (-C) -> MarkDuplicates (BARCODE_TAG)."""
        cfg = _write_config(
            tmp_path,
            "cfg_umi_on.yaml",
            trimming={"enabled": True},
            umi={"enabled": True, "preset": "agilent_xths2", "duplex": True},
        )
        result = _run_snakemake("-p", config=cfg)
        assert result.returncode == 0, result.stderr
        assert "rule umi_extract:" in result.stdout
        assert "--read-structures 3M2S+T 3M2S+T" in result.stdout
        assert "--umi-tag RX" in result.stdout
        assert "samtools fastq" in result.stdout
        assert "-T RX" in result.stdout
        assert "-K 100000000 -Y -C" in result.stdout
        assert "--BARCODE_TAG RX --DUPLEX_UMI true" in result.stdout

        cfg_custom = _write_config(
            tmp_path,
            "cfg_umi_custom.yaml",
            trimming={"enabled": True},
            umi={"enabled": True, "umi_tag": "OX", "read_structures": "3M2S+T\n3M2S+T"},
        )
        res_custom = _run_snakemake("-p", config=cfg_custom)
        assert res_custom.returncode == 0, res_custom.stderr
        assert "--read-structures 3M2S+T 3M2S+T" in res_custom.stdout
        assert "--umi-tag OX" in res_custom.stdout
        assert "-T OX" in res_custom.stdout
        assert "--BARCODE_TAG OX --DUPLEX_UMI true" in res_custom.stdout

    def test_umi_enabled_fixture_end_to_end(self, tmp_path):
        """Execute umi_extract -> trim_adapters -> bwa_map -> merge -> dedup on a paired FASTQ fixture with MBC + dark bases + read-through adapter."""
        import gzip

        def revcomp(seq: str) -> str:
            return seq.translate(str.maketrans("ACGTB", "TGCAV"))[::-1]

        ref_seq = (
            "GCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCG"
            "ATCGTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCG"
            "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
            "CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA"
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        )
        ref_dir = tmp_path / "ref"
        ref_dir.mkdir()
        ref_fa = ref_dir / "tiny.fna"
        ref_fa.write_text(f">chr1\n{ref_seq}\n")
        for ext in (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"):
            (ref_dir / f"tiny.fna{ext}").touch()
        (ref_dir / "tiny.dict").touch()

        insert = ref_seq[20:90]  # 70 bp genomic template
        r1_mbc, r1_dark = "AAA", "TG"
        r2_mbc, r2_dark = "BBB", "CA"
        r1_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        r1_3prime_remnant = "TGGGG" + r1_adapter[:20]
        r2_3prime_remnant = "CATTT" + r2_adapter[:20]

        r1_full = r1_mbc + r1_dark + insert + r1_3prime_remnant
        r2_full = r2_mbc + r2_dark + revcomp(insert) + r2_3prime_remnant

        fq_dir = tmp_path / "fastqs"
        fq_dir.mkdir()
        with gzip.open(fq_dir / "SampleA_S1_L001_R1_001.fastq.gz", "wt") as f:
            f.write(f"@read1 1:N:0:ATCACG\n{r1_full}\n+\n" + "I" * len(r1_full) + "\n")
        with gzip.open(fq_dir / "SampleA_S1_L001_R2_001.fastq.gz", "wt") as f:
            f.write(f"@read1 2:N:0:ATCACG\n{r2_full}\n+\n" + "I" * len(r2_full) + "\n")

        samples_tsv = tmp_path / "samples.tsv"
        samples_tsv.write_text(
            "fastq_files_basename\tlane\tproject_sample\tmdc_project\n"
            "SampleA_S1_L001\tL001\tSampleA\tTestProject\n"
        )

        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        pkgs_root = Path.home() / "miniforge3" / "pkgs"

        def _find_tool(name: str, pkg_glob: str) -> str | None:
            found = shutil.which(name)
            if found:
                return found
            if pkgs_root.is_dir():
                matches = sorted(pkgs_root.glob(f"{pkg_glob}/bin/{name}"))
                if matches:
                    return str(matches[-1])
            return None

        tools = {
            "fgbio": _find_tool("fgbio", "fgbio-*"),
            "samtools": _find_tool("samtools", "samtools-1.2*"),
            "bbduk.sh": _find_tool("bbduk.sh", "bbmap-*"),
            "bwa": _find_tool("bwa", "bwa-*"),
        }
        missing = [name for name, path in tools.items() if not path]
        assert not missing, f"Missing required bioinformatics tools: {missing}"
        for name, target in tools.items():
            assert target is not None
            (bin_dir / name).symlink_to(target)
        assert tools["bwa"] is not None
        subprocess.run([tools["bwa"], "index", str(ref_fa)], check=True, capture_output=True)

        markdup_args_file = tmp_path / "markdup_args.txt"
        gatk_shim = bin_dir / "gatk"
        gatk_shim.write_text(
            "#!/bin/bash\n"
            f'printf "%s " "$@" > "{markdup_args_file}"\n'
            "while [[ $# -gt 0 ]]; do\n"
            '  case "$1" in\n'
            '    -I) in_bam="$2"; shift 2 ;;\n'
            '    -O) out_bam="$2"; shift 2 ;;\n'
            '    -M) metrics="$2"; shift 2 ;;\n'
            "    *) shift ;;\n"
            "  esac\n"
            "done\n"
            'cp "$in_bam" "$out_bam"\n'
            'touch "${out_bam%.bam}.bai" "$metrics"\n'
        )
        gatk_shim.chmod(0o755)

        out_dir = tmp_path / "out"
        cfg = _write_config(
            tmp_path,
            "cfg_e2e.yaml",
            ref={"genome": str(ref_fa)},
            paths={
                "samples": str(samples_tsv),
                "fastq_folder": str(fq_dir),
                "output_folder": str(out_dir),
            },
            trimming={"enabled": True},
            qc={"enabled": False},
            umi={"enabled": True, "preset": "agilent_xths2", "duplex": True},
        )

        target_dedup = out_dir / "dedup" / "SampleA.merged.dedup.bam"
        snk_bin = shutil.which("snakemake")
        assert snk_bin is not None
        snk_dir = str(Path(snk_bin).parent)
        env = {
            **os.environ,
            "PATH": f"{bin_dir}:{snk_dir}:{os.environ.get('PATH', '/usr/bin:/bin')}",
        }

        run_res = subprocess.run(
            [
                snk_bin,
                "-s",
                str(SNAKEFILE),
                "--configfile",
                str(cfg),
                "--workflow-profile",
                "none",
                "--directory",
                str(tmp_path),
                "--cores",
                "2",
                "--notemp",
                str(target_dedup),
            ],
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert run_res.returncode == 0, (
            f"End-to-end UMI run failed:\nstdout: {run_res.stdout}\nstderr: {run_res.stderr}"
        )

        # 1. Verify umi_extract stripped first 5 bases (3M+2S) and wrote RX:Z:AAA-BBB comment
        umi_r1 = out_dir / "umi_extracted" / "SampleA_S1_L001.umi_R1_001.fastq.gz"
        umi_r2 = out_dir / "umi_extracted" / "SampleA_S1_L001.umi_R2_001.fastq.gz"
        with gzip.open(umi_r1, "rt") as f:
            u1_lines = [line.strip() for line in f]
        with gzip.open(umi_r2, "rt") as f:
            u2_lines = [line.strip() for line in f]
        assert u1_lines[0] == "@read1\tRX:Z:AAA-BBB"
        assert u2_lines[0] == "@read1\tRX:Z:AAA-BBB"
        assert u1_lines[1] == r1_full[5:]
        assert u2_lines[1] == r2_full[5:]

        # 2. Verify bbduk preserved RX comment and removed 3' mate-MBC remnant + adapter
        trim_r1 = out_dir / "bbduk_trimmed" / "SampleA_S1_L001.bbduk_R1_001.fastq.gz"
        trim_r2 = out_dir / "bbduk_trimmed" / "SampleA_S1_L001.bbduk_R2_001.fastq.gz"
        with gzip.open(trim_r1, "rt") as f:
            t1_lines = [line.strip() for line in f]
        with gzip.open(trim_r2, "rt") as f:
            t2_lines = [line.strip() for line in f]
        assert t1_lines[0] == "@read1\tRX:Z:AAA-BBB"
        assert t2_lines[0] == "@read1\tRX:Z:AAA-BBB"
        assert t1_lines[1] == insert
        assert t2_lines[1] == revcomp(insert)

        # 3. Verify aligned BAM has RX:Z:AAA-BBB on all mapped reads and 70M CIGAR
        aligned_bam = out_dir / "aligned" / "SampleA_S1_L001.bam"
        sam_view = subprocess.run(
            [str(bin_dir / "samtools"), "view", str(aligned_bam)],
            env=env,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
        assert len(sam_view) == 2
        for record in sam_view:
            assert "RX:Z:AAA-BBB" in record
            assert "\t70M\t" in record

        # 4. Verify MarkDuplicates received BARCODE_TAG RX and DUPLEX_UMI true
        markdup_args = markdup_args_file.read_text()
        assert "--BARCODE_TAG RX" in markdup_args
        assert "--DUPLEX_UMI true" in markdup_args

        # 5. Verify UMI provenance file written with resolved read structure
        prov_file = out_dir / "logs" / "umi_provenance.yaml"
        assert prov_file.exists()
        prov = yaml.safe_load(prov_file.read_text())
        assert prov["umi"]["read_structures"] == "3M2S+T 3M2S+T"
        assert prov["umi"]["preset"] == "agilent_xths2"
