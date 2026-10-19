"""Schema validation tests for config and samples schemas.

These tests validate config/samples data against the JSON schemas
in workflow/schemas/ using the same validation mechanism as Snakemake.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pandas as pd
import pytest
from snakemake.exceptions import WorkflowError as _SchemaError
from snakemake.utils import validate as snakemake_validate

SCHEMA_DIR = Path(__file__).resolve().parent.parent / "workflow" / "schemas"
DATA_DIR = Path(__file__).resolve().parent / "data"


def _load_yaml(path: Path) -> dict:
    """Load a YAML file as a dict."""
    import yaml

    with open(path) as f:
        return yaml.safe_load(f)


# ============================================================================
# Config schema tests
# ============================================================================


class TestConfigSchema:
    """Tests for config.schema.yaml validation."""

    @pytest.fixture()
    def valid_config(self):
        return _load_yaml(DATA_DIR / "config" / "config.yaml")

    def test_valid_minimal(self, valid_config):
        snakemake_validate(valid_config, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_missing_ref(self, valid_config):
        cfg = copy.deepcopy(valid_config)
        del cfg["ref"]
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_missing_paths(self, valid_config):
        cfg = copy.deepcopy(valid_config)
        del cfg["paths"]
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_invalid_build_enum(self):
        cfg = _load_yaml(DATA_DIR / "config" / "config_invalid_build.yaml")
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_empty_known_sites(self, valid_config):
        cfg = copy.deepcopy(valid_config)
        cfg["ref"]["known_sites"] = []
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_compression_bounds(self, valid_config):
        cfg = copy.deepcopy(valid_config)
        cfg["processing"] = {"compression_level": 0}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_actual_repo_config(self):
        """Validate the actual config/config.yaml from the repo."""
        repo_config = Path(__file__).resolve().parent.parent / "config" / "config.yaml"
        assert repo_config.exists()
        cfg = _load_yaml(repo_config)
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_qc_defaults_valid(self, valid_config):
        """QC section is optional — config without it should pass."""
        cfg = copy.deepcopy(valid_config)
        cfg.pop("qc", None)
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_qc_all_flags(self, valid_config):
        """QC section with all flags set should pass."""
        cfg = copy.deepcopy(valid_config)
        cfg["qc"] = {
            "enabled": True,
            "fastqc": True,
            "samtools_stats": True,
            "samtools_flagstat": True,
            "picard_collect_metrics": True,
            "qualimap": False,
            "qualimap_feature_file": "",
        }
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_qc_disabled(self, valid_config):
        """QC disabled should pass."""
        cfg = copy.deepcopy(valid_config)
        cfg["qc"] = {"enabled": False}
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_qc_invalid_type(self, valid_config):
        """String where boolean expected should fail."""
        cfg = copy.deepcopy(valid_config)
        cfg["qc"] = {"enabled": "yes"}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_umi_disabled_defaults_valid(self, valid_config):
        """UMI disabled with null preset and read_structures should pass."""
        cfg = copy.deepcopy(valid_config)
        cfg["umi"] = {
            "enabled": False,
            "preset": None,
            "read_structures": None,
            "umi_tag": "RX",
            "duplex": True,
            "dedup": {"barcode_tag": "RX", "extra": ""},
            "fgbio": {"fastq_to_bam_extra": "", "mem_mb": 8000},
        }
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_umi_preset_agilent_xths2_valid(self, valid_config):
        """UMI enabled with preset=agilent_xths2 and null read_structures should pass."""
        cfg = copy.deepcopy(valid_config)
        cfg["umi"] = {"enabled": True, "preset": "agilent_xths2", "read_structures": None}
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_umi_custom_read_structures_valid(self, valid_config):
        """UMI enabled with valid custom read_structures (preset=null or custom) should pass."""
        for preset in (None, "", "custom"):
            cfg = copy.deepcopy(valid_config)
            cfg["umi"] = {
                "enabled": True,
                "preset": preset,
                "read_structures": "3M2S+T 3M2S+T",
            }
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    @pytest.mark.parametrize(
        "bad_structure",
        [
            "invalid",
            "3X2S+T 3X2S+T",
            "+M+T +M+T",
            "3M2S+T",
            "+T +T",
            "3M2S 3M2S",
            "0M+T 0M+T",
        ],
    )
    def test_umi_invalid_read_structures_rejected(self, valid_config, bad_structure):
        """Invalid fgbio read_structures syntax must be rejected by schema."""
        cfg = copy.deepcopy(valid_config)
        cfg["umi"] = {"enabled": True, "read_structures": bad_structure}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_umi_preset_and_read_structures_conflict_rejected(self, valid_config):
        """Setting both a kit preset and read_structures when enabled must be rejected."""
        cfg = copy.deepcopy(valid_config)
        cfg["umi"] = {
            "enabled": True,
            "preset": "agilent_xths2",
            "read_structures": "5M+T 5M+T",
        }
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    @pytest.mark.parametrize("preset", [None, "", "custom"])
    def test_umi_enabled_missing_both_rejected(self, valid_config, preset):
        """Enabling UMI without preset or read_structures must be rejected."""
        cfg = copy.deepcopy(valid_config)
        cfg["umi"] = {"enabled": True, "preset": preset, "read_structures": None}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_umi_unknown_preset_rejected(self, valid_config):
        """Unknown UMI preset name must be rejected by schema enum."""
        cfg = copy.deepcopy(valid_config)
        cfg["umi"] = {"enabled": True, "preset": "nonexistent_kit"}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_params_multiqc_extra_valid(self, valid_config):
        """params.multiqc.extra string should pass schema validation."""
        cfg = copy.deepcopy(valid_config)
        cfg.setdefault("params", {})["multiqc"] = {"extra": "--title MyReport"}
        snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))

    def test_params_multiqc_extra_invalid_type(self, valid_config):
        """Non-string params.multiqc.extra must be rejected by schema."""
        cfg = copy.deepcopy(valid_config)
        cfg.setdefault("params", {})["multiqc"] = {"extra": 123}
        with pytest.raises(_SchemaError):
            snakemake_validate(cfg, str(SCHEMA_DIR / "config.schema.yaml"))


# ============================================================================
# Samples schema tests
# ============================================================================


class TestSamplesSchema:
    """Tests for samples.schema.yaml validation."""

    def test_valid_samples(self):
        df = pd.read_table(DATA_DIR / "config" / "samples.tsv")
        snakemake_validate(df, str(SCHEMA_DIR / "samples.schema.yaml"))

    def test_missing_required_column(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["S1"],
                "lane": ["L001"],
                # missing project_sample and mdc_project
            }
        )
        with pytest.raises(_SchemaError):
            snakemake_validate(df, str(SCHEMA_DIR / "samples.schema.yaml"))

    def test_optional_subfolder(self):
        df = pd.DataFrame(
            {
                "fastq_files_basename": ["S1"],
                "lane": ["L001"],
                "project_sample": ["Sample"],
                "mdc_project": ["Proj"],
                "subfolder": ["run1"],
            }
        )
        snakemake_validate(df, str(SCHEMA_DIR / "samples.schema.yaml"))
