# sm-alignment

Snakemake pipeline for DNA sequence alignment and BAM processing on SLURM HPC clusters.

**FASTQ → BWA → merge lanes → MarkDuplicates → BQSR → analysis-ready BAM**

## Quick Start

```bash
# 1. Configure
cp config/samples.tsv config/samples.tsv        # edit with your sample info
vim config/config.yaml                           # set reference paths, known-sites, FASTQ location

# 2. Dry-run (verify DAG without executing)
snakemake -s workflow/Snakefile --configfile config/config.yaml -n

# 3. Submit to SLURM
sbatch scripts/run_snakemake.sh workflow/Snakefile
```

## Configuration

### `config/config.yaml`

| Section | Key settings |
|---|---|
| `ref` | `genome`, `build`, `known_sites` (list of VCFs for BQSR) |
| `paths` | `fastq_folder`, `output_folder`, `samples` (path to TSV) |
| `trimming` | Set `enabled: true` to run BBDuk adapter trimming before alignment |
| `params` | Extra CLI flags passed through to bwa, samtools, GATK tools |

### `config/samples.tsv`

One row per FASTQ pair (per lane). Columns:

| Column | Example | Description |
|---|---|---|
| `fastq_files_basename` | `SampleA_S1_L001` | FASTQ filename prefix (before `_R1_001.fastq.gz`) |
| `lane` | `L001` | Sequencing lane (used for read groups) |
| `project_sample` | `SampleA` | Sample name (lane BAMs are merged by this) |
| `mdc_project` | `ProjectX` | Project identifier (read group PU tag) |

### `profiles/default/config.yaml`

Per-rule resource allocation (threads, memory, walltime). Adjust for your cluster without touching workflow code.

## Project Structure

```
workflow/
├── Snakefile              # Entry point (schema validation, rule all)
├── rules/                 # Modular rule files
├── envs/                  # Conda environments (pinned versions)
└── schemas/               # JSON Schema for config + samples validation
config/                    # Pipeline config + sample metadata
profiles/default/          # Resource allocation profile
scripts/run_snakemake.sh   # Single SLURM launcher
```

## Usage Examples

```bash
# Default: config/config.yaml, cubi-v1 cluster profile
sbatch scripts/run_snakemake.sh workflow/Snakefile

# Custom config
sbatch scripts/run_snakemake.sh workflow/Snakefile config/my_config.yaml

# Custom job name
sbatch --job-name=sm_exomes scripts/run_snakemake.sh workflow/Snakefile

# Override resources for a single run
snakemake -s workflow/Snakefile --set-threads bwa_map=8 --set-resources bwa_map:mem_mb=9600

# Pass extra Snakemake flags
sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml cubi-v1 --forceall
```

## Tools & Environments

| Conda env | Tools | Used by |
|---|---|---|
| `workflow/envs/bwa_samtools.yaml` | bwa 0.7.18, samtools 1.21, samblaster 0.1.26 | alignment, merge, utilities |
| `workflow/envs/gatk.yaml` | gatk4 4.6.1.0, samtools 1.21 | dedup, BQSR |
| `workflow/envs/bbtools.yaml` | bbmap 39.06 | trimming |

## Deprecated Files

Old standalone workflows and SLURM launchers are in `deprecated/` for reference during migration. See `deprecated/README.md` for details. Scheduled for removal August 2025.

## License

MIT
