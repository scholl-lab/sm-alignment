"""Schema validation tests for config and samples schemas.

These tests validate config/samples data against the JSON schemas
in workflow/schemas/ using the same validation mechanism as Snakemake.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pandas as pd
import pytest

try:
    from snakemake.utils import validate as snakemake_validate

    HAS_SNAKEMAKE = True
except ImportError:
    HAS_SNAKEMAKE = False

# snakemake.validate() may raise WorkflowError or jsonschema.ValidationError
# depending on the version; we catch the broadest applicable base.
try:
    from snakemake.exceptions import WorkflowError as _SchemaError
except Exception:
    _SchemaError = Exception  # type: ignore[assignment,misc]

pytestmark = pytest.mark.skipif(not HAS_SNAKEMAKE, reason="snakemake not installed")

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
        if not repo_config.exists():
            pytest.skip("No repo config/config.yaml")
        cfg = _load_yaml(repo_config)
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
