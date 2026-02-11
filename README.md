# sm-alignment

Snakemake pipeline for DNA sequence alignment and BAM processing on SLURM HPC clusters.

**FASTQ → (BBDuk trim) → BWA → merge lanes → MarkDuplicates → BQSR → analysis-ready BAM**

Supports both **BIH HPC** (`cubi-v1` profile) and **Charité HPC** (auto-detected).

---

## Setup

### 1. Clone into your project directory

The pipeline should live inside each project's directory on the cluster, not in your home directory (quota is too small for results).

```bash
# Charité HPC
cd /sc-projects/<your-project>
git clone https://github.com/scholl-lab/sm-alignment.git
cd sm-alignment

# BIH HPC
cd /data/cephfs-1/work/projects/<your-project>
git clone https://github.com/scholl-lab/sm-alignment.git
cd sm-alignment
```

Recommended project layout:

```
/sc-projects/<your-project>/         # or /data/.../projects/<your-project>/
├── sm-alignment/                    # this pipeline (git clone)
│   ├── workflow/
│   ├── config/
│   │   ├── config.yaml              # ← edit for your project
│   │   └── samples.tsv              # ← generate from SampleSheet
│   ├── profiles/
│   └── scripts/
├── resources/                       # reference data (or symlinks to shared)
│   ├── ref/GRCh38/
│   └── gatk_bundle/hg38/
└── data/                            # raw FASTQ deliveries (or symlinks)
```

### 2. Set up reference data

Follow the [lab handbook: Reference Data Setup](https://github.com/scholl-lab/lab-handbook/blob/main/docs/reference-data-setup.md) to download the reference genome and GATK known-sites.

**BIH HPC** — shared data already exists:

```
/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38/
/data/cephfs-1/work/projects/apa-sequencing/analysis/GATK_resource_bundle/
```

**Charité HPC** — set up per project:

```
/sc-projects/<your-project>/resources/ref/GRCh38/
/sc-projects/<your-project>/resources/gatk_bundle/hg38/
```

---

## Generate Config Files

### Interactive wizard (recommended for first-time setup)

Run without arguments to get a step-by-step guided setup:

```bash
python scripts/generate_config.py
```

The wizard prompts for FASTQ directory, SampleSheet, project name, reference data location, and output paths. It auto-detects defaults where possible.

### From a sequencing facility delivery (flags mode)

When you receive FASTQ files with a `SampleSheet.csv` from the core facility:

```bash
# Generate samples.tsv from the delivery folder
python scripts/generate_config.py \
    --fastq-dir /path/to/250903_LH00253_0332_B232J72LT4_A5297_FASTQ/

# Output:
#   Found SampleSheet.csv with 2 samples
#   Matched 4 FASTQ files (2 pairs)
#   Written: config/samples.tsv (2 samples)
```

The script auto-detects `SampleSheet.csv` in the FASTQ directory, parses Illumina filenames, and generates `config/samples.tsv`. It handles both standard Illumina `[Data]` format and the BIH/Charité minimal CSV format.

**Options:**

```bash
# Preview without writing (dry-run)
python scripts/generate_config.py --fastq-dir /path/to/fastqs --dry-run

# Generate samples.tsv AND config.yaml (scans for reference genome + known-sites)
python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template

# Point to a specific reference data directory
python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template \
    --ref-dir /sc-projects/<project>/resources/ref/GRCh38

# Override project name
python scripts/generate_config.py --fastq-dir /path/to/fastqs --project A5297

# Overwrite existing files
python scripts/generate_config.py --fastq-dir /path/to/fastqs --force
```

When `--config-template` is used, the script scans for reference data:

1. **Explicit** `--ref-dir` (if provided)
2. **Relative** paths near the FASTQ directory (`resources/ref/GRCh38/`, `resources/gatk_bundle/hg38/`, etc.)
3. **Shared BIH HPC** locations (`/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38/`)

It checks for companion files (BWA index, FAI, dict, tabix index) and reports what it finds. Discovered paths are written directly into `config/config.yaml`. Any paths not found are marked with `EDIT_ME:` placeholders.

### Review config/config.yaml

After generating, review the config and adjust paths if needed:

```yaml
ref:
  genome: "/sc-projects/<project>/resources/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  known_sites:
    - "/sc-projects/<project>/resources/gatk_bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
    - "/sc-projects/<project>/resources/gatk_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
    - "/sc-projects/<project>/resources/gatk_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

paths:
  fastq_folder: "/path/to/your/fastq/directory"   # auto-filled from --fastq-dir
  output_folder: "results/A5297"                   # auto-filled from project name

trimming:
  enabled: false    # set to true to auto-trim before alignment
```

---

## Run the Pipeline

### Dry-run first

Always verify the execution plan before submitting:

```bash
# On a compute node (Charité) or login node (BIH)
conda activate snakemake

snakemake -s workflow/Snakefile --configfile config/config.yaml -n
```

### Submit to SLURM

```bash
# Create log directory
mkdir -p slurm_logs

# Submit — cluster is auto-detected
sbatch scripts/run_snakemake.sh workflow/Snakefile
```

The launcher auto-detects whether you're on BIH HPC or Charité HPC:

| Cluster | Detection | Behavior |
|---------|-----------|----------|
| **BIH HPC** | `cubi-v1` profile exists | Uses `--profile=cubi-v1` for job submission |
| **Charité HPC** | `/etc/profile.d/conda.sh` exists | Uses SLURM executor plugin (`profiles/charite`), sources conda automatically |
| **Other/local** | Fallback | Runs without cluster submission |

### Usage examples

```bash
# Custom config file
sbatch scripts/run_snakemake.sh workflow/Snakefile config/my_config.yaml

# Custom job name (visible in squeue)
sbatch --job-name=sm_exomes scripts/run_snakemake.sh workflow/Snakefile

# Pass extra Snakemake flags
sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml --forceall

# Override resources for a single run
sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml \
    --set-threads bwa_map=8 --set-resources bwa_map:mem_mb=9600
```

### Monitor progress

```bash
# Check running jobs
squeue -u $USER

# Follow coordinator log
tail -f slurm_logs/sm_pipeline-*.log

# Follow individual rule logs
tail -f slurm_logs/slurm-*.log
```

---

## Configuration Reference

### `config/config.yaml`

| Section | Key settings |
|---|---|
| `ref` | `genome`, `build`, `known_sites` (list of VCFs for BQSR) |
| `paths` | `fastq_folder`, `output_folder`, `samples` (path to TSV) |
| `trimming` | Set `enabled: true` to run BBDuk adapter trimming before alignment |
| `params` | Extra CLI flags passed through to bwa, samtools, GATK tools |
| `subset` | Optional BAM subsetting by BED regions |

### `config/samples.tsv`

One row per FASTQ pair (per lane):

| Column | Example | Description |
|---|---|---|
| `fastq_files_basename` | `A5297_DNA_01_STREAM_P1_L1_S1_L008` | FASTQ filename prefix (before `_R1_001.fastq.gz`) |
| `lane` | `L008` | Sequencing lane (used for read groups) |
| `project_sample` | `A5297_DNA_01_STREAM_P1_L1` | Sample name (lane BAMs are merged by this) |
| `mdc_project` | `A5297` | Project identifier (read group PU tag) |
| `subfolder` | `run1` | *(optional)* Subdirectory within FASTQ folder |

### `profiles/default/config.yaml`

Per-rule resource allocation (threads, memory, walltime). Adjust for your cluster without touching workflow code. Runtimes are capped at 2880 min (48h) for Charité compatibility.

### `profiles/charite/config.yaml`

Charité-specific SLURM submission settings. Used automatically when the launcher detects the Charité cluster.

---

## Tools & Conda Environments

| Conda env | Tools | Used by |
|---|---|---|
| `workflow/envs/bwa_samtools.yaml` | bwa 0.7.18, samtools 1.21, samblaster 0.1.26 | alignment, merge, utilities |
| `workflow/envs/gatk.yaml` | gatk4 4.6.1.0, samtools 1.21 | dedup, BQSR |
| `workflow/envs/bbtools.yaml` | bbmap 39.06 | trimming |

Conda environments are created automatically by Snakemake on first run (`software-deployment-method: conda` in the workflow profile).

### Software Deployment

The pipeline supports two software deployment strategies:

**1. Per-rule conda environments (default)**

Snakemake creates isolated conda environments from `workflow/envs/*.yaml` on the first run. This is fully reproducible but slow initially. Use `--conda-prefix` to share envs across projects:

```bash
snakemake --conda-prefix /shared/conda-envs ...
```

**2. Pre-installed tools (skip per-rule conda)**

Install all tools directly into your `snakemake` environment and disable per-rule conda. This avoids the first-run env creation overhead and works around mamba 2.x incompatibilities:

```bash
# Install tools into the snakemake environment
mamba install -n snakemake -c bioconda -c conda-forge \
    bwa=0.7.18 samtools=1.21 samblaster=0.1.26 gatk4=4.6.1.0 bbmap=39.06
```

Then comment out `software-deployment-method` in `profiles/default/config.yaml` or pass `--sdm none` on the CLI:

```bash
sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml --sdm none
```

---

## Deprecated Files

Old standalone workflows and SLURM launchers are in `deprecated/` for reference during migration. See `deprecated/README.md` for details. Scheduled for removal August 2025.

## License

MIT
