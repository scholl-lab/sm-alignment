# Deprecated Files

These files are superseded by the modular `workflow/` pipeline structure introduced in February 2025. They are retained here temporarily for reference during migration.

**Planned removal: August 2025**

## What was replaced

| Old location | Replaced by |
|---|---|
| `scripts/snakemake/*.smk` (12 standalone workflows) | `workflow/rules/*.smk` (7 modular rules included by `workflow/Snakefile`) |
| `scripts/run_*.sh`, `scripts/submit_*.sh` (11 SLURM launchers) | `scripts/run_snakemake.sh` (single generic launcher) |
| `configs/*.yaml`, root `config.yaml` (3 config files) | `config/config.yaml` (unified hierarchical config with schema validation) |

## Migration

```bash
# Old way (one script per step):
sbatch scripts/run_alignment.sh
sbatch scripts/run_merge_bams.sh
sbatch scripts/run_dedup_bams.sh

# New way (single launcher, full pipeline):
sbatch scripts/run_snakemake.sh workflow/Snakefile
```

See `README.md` in the repository root for full usage instructions.

## Structure

```
deprecated/
├── configs/                    # Old configuration files
│   ├── config_alignment.yaml
│   ├── config_legacy_root.yaml # Was config.yaml at repo root
│   └── config_trim_adapters.yaml
└── scripts/
    ├── launchers/              # Old per-step SLURM submission scripts
    └── snakemake/              # Old standalone Snakemake workflows
```
