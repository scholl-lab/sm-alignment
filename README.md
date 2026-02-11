# sm-alignment

Snakemake pipeline for DNA sequence alignment and BAM processing (BWA + GATK).

## Quick Start

1. Edit `config/config.yaml` with your reference paths, known-sites, and sample metadata
2. Edit `config/samples.tsv` with your sample information
3. Adjust resource allocation in `profiles/default/config.yaml` if needed
4. Submit to SLURM:

```bash
sbatch scripts/run_snakemake.sh workflow/Snakefile
```

## Pipeline Stages

FASTQ → BWA alignment → merge lanes → GATK MarkDuplicates → GATK BQSR → analysis-ready BAM

Optional: BBDuk adapter trimming (set `trimming.enabled: true` in config)

## Dry Run

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml -n
```

## Configuration

| File | Purpose |
|---|---|
| `config/config.yaml` | Pipeline settings (references, paths, tool parameters) |
| `config/samples.tsv` | Sample metadata (basename, lane, sample name, project) |
| `profiles/default/config.yaml` | Resource allocation (threads, memory, time per rule) |

## Project Structure

```
workflow/           Snakemake workflow (Snakefile, rules/, envs/, schemas/)
config/             Configuration and sample metadata
profiles/default/   Resource allocation profile
scripts/            SLURM submission script
```

## Legacy Workflows

The standalone workflows in `scripts/snakemake/` and individual submission scripts in `scripts/run_*.sh` / `scripts/submit_*.sh` are retained for reference but superseded by the modular `workflow/` structure.
