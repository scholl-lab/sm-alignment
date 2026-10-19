# sm-alignment

Snakemake pipeline for DNA sequence alignment and BAM processing on SLURM HPC clusters.

Designed for SLURM HPC clusters with automatic cluster detection.

```mermaid
flowchart LR
    FASTQ["Paired-end\nFASTQ files"]
    UMI["fgbio + samtools\nUMI extract"]
    TRIM["BBDuk\ntrim"]
    ALIGN["BWA MEM\nalign"]
    MERGE["samtools merge\nper sample"]
    DEDUP["GATK\nMarkDuplicates"]
    BQSR["GATK\nBQSR"]
    BAM["Analysis-ready\nBAM"]
    FQC_RAW["FastQC\nraw"]
    FQC_TRIM["FastQC\ntrimmed"]
    QC["samtools stats\nPicard metrics\nQualimap"]
    MULTIQC["MultiQC\nreport"]

    FASTQ -.->|umi\nenabled| UMI
    UMI -.-> TRIM
    FASTQ --> TRIM
    TRIM -.->|trimming\nenabled| ALIGN
    FASTQ -->|per lane| ALIGN
    ALIGN --> MERGE
    MERGE --> DEDUP
    DEDUP --> BQSR
    BQSR --> BAM

    FASTQ -.-> FQC_RAW
    TRIM -.-> FQC_TRIM
    BAM -.-> QC
    FQC_RAW -.-> MULTIQC
    FQC_TRIM -.-> MULTIQC
    QC -.-> MULTIQC
    DEDUP -.->|dedup\nmetrics| MULTIQC

    style UMI stroke-dasharray: 5 5
    style TRIM stroke-dasharray: 5 5
    style FQC_RAW stroke-dasharray: 5 5
    style FQC_TRIM stroke-dasharray: 5 5
    style QC stroke-dasharray: 5 5
    style MULTIQC stroke-dasharray: 5 5
```

---

## Setup

### Prerequisites

You need **Snakemake 8+** and **conda/mamba**. Install via [Miniforge](https://github.com/conda-forge/miniforge) (do not use Anaconda):

```bash
mamba create -n snakemake -y -c conda-forge -c bioconda python=3.11 snakemake=8
conda activate snakemake
```

Pipeline tools (BWA, GATK, samtools, BBDuk) are installed automatically as per-rule conda environments on first run — no manual tool installation needed. See [Software Deployment](#software-deployment) for alternatives.

### 1. Clone the repository

```bash
git clone https://github.com/scholl-lab/sm-alignment.git
cd sm-alignment
```

Recommended project layout:

```
<your-project>/
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

The pipeline requires a reference genome and known-sites VCFs for base quality score recalibration (BQSR):

| Data | Source | Maps to `config.yaml` |
|------|--------|-----------------------|
| **GRCh38 no-alt analysis set** (FASTA + BWA index + .fai + .dict) | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) | `ref.genome` |
| **dbSNP 138** | [Google Cloud — Broad references](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0) | `ref.known_sites[]` |
| **Known indels** | same bucket | `ref.known_sites[]` |
| **Mills & 1000G gold-standard indels** | same bucket | `ref.known_sites[]` |

Download the files and place them under `resources/` as shown in the project layout above. Then point `config/config.yaml` at the paths, or use `generate_config.py --config-template` to auto-detect them (see [Generate Config Files](#generate-config-files)).

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
    --fastq-dir /path/to/delivery_folder/

# Output:
#   Found SampleSheet.csv with 2 samples
#   Matched 4 FASTQ files (2 pairs)
#   Written: config/samples.tsv (2 samples)
```

The script auto-detects `SampleSheet.csv` in the FASTQ directory, parses Illumina filenames, and generates `config/samples.tsv`. It handles both standard Illumina `[Data]` format and minimal CSV formats.

**Options:**

```bash
# Preview without writing (dry-run)
python scripts/generate_config.py --fastq-dir /path/to/fastqs --dry-run

# Generate samples.tsv AND config.yaml (scans for reference genome + known-sites)
python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template

# Point to a specific reference data directory
python scripts/generate_config.py --fastq-dir /path/to/fastqs --config-template \
    --ref-dir /path/to/resources/ref/GRCh38

# Override project name
python scripts/generate_config.py --fastq-dir /path/to/fastqs --project MyProject

# Overwrite existing files
python scripts/generate_config.py --fastq-dir /path/to/fastqs --force
```

When `--config-template` is used, the script scans for reference data:

1. **Explicit** `--ref-dir` (if provided)
2. **Relative** paths near the FASTQ directory (`resources/ref/GRCh38/`, `resources/gatk_bundle/hg38/`, etc.)

It checks for companion files (BWA index, FAI, dict, tabix index) and reports what it finds. Discovered paths are written directly into `config/config.yaml`. Any paths not found are marked with `EDIT_ME:` placeholders.

### Review config/config.yaml

After generating, review the config and adjust paths if needed:

```yaml
ref:
  genome: "/path/to/resources/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  known_sites:
    - "/path/to/resources/gatk_bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
    - "/path/to/resources/gatk_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
    - "/path/to/resources/gatk_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

paths:
  fastq_folder: "/path/to/your/fastq/directory"   # auto-filled from --fastq-dir
  output_folder: "results/MyProject"               # auto-filled from project name

trimming:
  enabled: false    # set to true to auto-trim before alignment
```

---

## Run the Pipeline

### Dry-run first

Always verify the execution plan before submitting:

```bash
# On a compute or login node
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

The launcher auto-detects your cluster environment and selects the appropriate SLURM submission method. If no known cluster is detected, it falls back to local execution.

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
| `umi` | UMI / molecular barcode extraction (`fgbio FastqToBam`) and UMI-aware duplicate marking (`--BARCODE_TAG`, `--DUPLEX_UMI`) |
| `qc` | Quality control: `enabled`, per-tool toggles (`fastqc`, `samtools_stats`, `samtools_flagstat`, `picard_collect_metrics`, `qualimap`) |
| `params` | Extra CLI flags passed through to bwa, samtools, GATK, and MultiQC |
| `subset` | Optional BAM subsetting by BED regions |

### `config/samples.tsv`

One row per FASTQ pair (per lane). Samples sequenced across multiple lanes have multiple rows — the pipeline aligns each lane independently and merges by `project_sample`.

| Column | Example | Required | Description |
|---|---|---|---|
| `fastq_files_basename` | `Sample1_S1_L001` | yes | FASTQ filename prefix (before `_R1_001.fastq.gz`). Must be unique — used as index. |
| `lane` | `L001` | yes | Sequencing lane. Used in the read group ID (`@RG ID:{lane}-{sample}`) and PU tag. |
| `project_sample` | `Sample1` | yes | Logical sample name. All rows sharing this value are merged into one BAM after alignment. |
| `mdc_project` | `ProjectX` | yes | Project identifier. Used in the read group PU tag (`PU:{lane}-{project}`). |
| `subfolder` | `run1` | no | Subdirectory within the FASTQ folder, for deliveries split across sub-directories. |

**Multi-lane example** — one sample sequenced on two lanes produces two rows:

```tsv
fastq_files_basename	lane	project_sample	mdc_project
Sample1_S1_L001	L001	Sample1	ProjectX
Sample1_S1_L002	L002	Sample1	ProjectX
```

This generates two lane-level BAMs that are merged into a single `Sample1.bam` before deduplication.

**Generating `samples.tsv`:**

- **Automatic** (recommended): `python scripts/generate_config.py --fastq-dir /path/to/fastqs` parses the Illumina `SampleSheet.csv` and FASTQ filenames. Run with `--help` for all options, or without arguments for the interactive wizard.
- **Manual**: create a tab-separated file with the columns above. Ensure `fastq_files_basename` matches the actual FASTQ filenames in your `paths.fastq_folder` directory.

### `profiles/default/config.yaml`

Per-rule resource allocation (threads, memory, walltime). Adjust for your cluster without touching workflow code.

### `profiles/charite/config.yaml`

Cluster-specific SLURM executor plugin settings. Used automatically when the launcher detects a matching cluster.

---

## UMI / Molecular Barcodes

Libraries with inline molecular barcodes (UMIs) require extracting the barcode and skipping stem/dark bases **before** 3' adapter/overlap trimming and alignment. Otherwise, 5' barcode bases produce soft-clips or read-end mismatches, and short inserts that read through into the opposite mate retain a 3' mate-barcode remnant ahead of the standard adapter sequence.

### Data flow

When `umi.enabled: true`, the pipeline runs:

```
fgbio FastqToBam -> samtools fastq -T RX -> bbduk -> bwa -C -> MarkDuplicates BARCODE_TAG
```

1. **`umi_extract`** (`fgbio FastqToBam` piped to `samtools fastq -T RX`): extracts inline UMI bases into the `RX` tag, strips 5' UMI and skip/dark bases from both R1 and R2 according to the resolved fgbio read structure, and emits paired FASTQs with `RX:Z:<R1_UMI>-<R2_UMI>` in the FASTQ header comment (dropping the original Illumina read comment so `bwa mem -C` emits a valid SAM tag).
2. **`trim_adapters`** (`bbduk.sh` with `tbo=t tpe=t`): because both 5' barcodes and dark bases have already been stripped, mate-overlap detection (`tbo=t`) identifies the exact genomic insert on short read-through fragments and trims both the 3' reverse-complement mate-MBC remnant and the 3' adapter while preserving the `RX:Z:` read comment.
3. **`bwa_map`** (`bwa mem -C`): appends the FASTQ header comment (`RX:Z:...`) directly to every aligned BAM record.
4. **`deduplicate_bam_files`** (`gatk MarkDuplicates`): passes `--BARCODE_TAG RX` (and `--DUPLEX_UMI true` when `umi.duplex: true`) for UMI-aware duplicate marking.
5. **Run provenance**: the resolved preset, fgbio read structure, UMI tag, and protocol source are logged and written to `{output_folder}/logs/umi_provenance.yaml`.

### Configuration reference

| Key | Type | Default | Description |
|---|---|---|---|
| `umi.enabled` | `bool` | `false` | Master switch. When `false`, DAG and commands are identical to standard non-UMI mode. |
| `umi.preset` | `str \| null` | `null` | Verified kit preset from `config/umi_presets.yaml` (`agilent_xths2` or `custom`). Mutually exclusive with `umi.read_structures` when enabled (unless `custom`). |
| `umi.read_structures` | `str \| null` | `null` | Paired-end [fgbio read structures](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) (e.g. `"3M2S+T 3M2S+T"`). Overrides preset. |
| `umi.umi_tag` | `str` | `"RX"` | Two-character SAM tag used by `fgbio FastqToBam` and `samtools fastq -T`. |
| `umi.duplex` | `bool` | `true` | Pass `--DUPLEX_UMI true` to GATK `MarkDuplicates`. |
| `umi.dedup.barcode_tag` | `str \| null` | `"RX"` | Tag passed to `MarkDuplicates --BARCODE_TAG` (`null` = position-only duplicate marking). |
| `umi.dedup.extra` | `str` | `""` | Extra CLI arguments appended to `MarkDuplicates` when UMI mode is enabled. |
| `umi.fgbio.fastq_to_bam_extra` | `str` | `""` | Extra CLI arguments passed through to `fgbio FastqToBam`. |
| `umi.fgbio.mem_mb` | `int` | `8000` | Default memory allocation (MB) for `umi_extract` (overridable via workflow profile). |

### Worked example: Agilent SureSelect XT HS2 (G9984A)

Agilent SureSelect XT HS2 DNA libraries carry a **3-nt molecular barcode (MBC) followed by 2 dark bases at the 5' end of both Read 1 and Read 2** (*Agilent SureSelect XT HS2 DNA System Protocol*, publication G9985-90010, pp. 128, 134). All 5 leading bases (`3M2S+T 3M2S+T`) must be removed prior to alignment:

```yaml
trimming:
  enabled: true

umi:
  enabled: true
  preset: agilent_xths2          # resolves to "3M2S+T 3M2S+T" (Agilent G9985-90010, pp. 128, 134)
  umi_tag: RX
  duplex: true
  dedup:
    barcode_tag: RX
```

### How to add a kit

1. **Ad-hoc custom kit (no code changes):** set `umi.preset: custom` (or `null`) and specify the paired-end fgbio read structure directly in `config/config.yaml`:
   ```yaml
   umi:
     enabled: true
     preset: custom
     read_structures: "5M2S+T 5M2S+T"
   ```
2. **Permanent verified preset:**
   - Verify the exact R1/R2 barcode and skip/stem structure in the vendor protocol documentation.
   - Add an entry to `config/umi_presets.yaml` with `description`, `read_structures`, `duplex`, and `source` (vendor document number and page).
   - Add the preset key to the `umi.preset` enum (and `then.oneOf` preset list) in `workflow/schemas/config.schema.yaml`.

---

## Tools & Conda Environments

| Conda env | Tools | Used by |
|---|---|---|
| `workflow/envs/bwa_samtools.yaml` | bwa 0.7.18, samtools 1.21, samblaster 0.1.26, fgbio 3.1.2 | UMI extraction, alignment, merge, utilities, samtools QC |
| `workflow/envs/gatk.yaml` | gatk4 4.6.1.0, samtools 1.21 | dedup, BQSR, Picard metrics |
| `workflow/envs/bbtools.yaml` | bbmap 39.06 | trimming |
| `workflow/envs/fastqc.yaml` | fastqc 0.12.1 | FASTQ quality control |
| `workflow/envs/multiqc.yaml` | multiqc 1.33 | QC report aggregation |
| `workflow/envs/qualimap.yaml` | qualimap 2.3 | alignment coverage analysis (opt-in) |

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
    bwa=0.7.18 samtools=1.21 samblaster=0.1.26 fgbio=3.1.2 gatk4=4.6.1.0 bbmap=39.06
```

Then comment out `software-deployment-method` in `profiles/default/config.yaml` or pass `--sdm none` on the CLI:

```bash
sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml --sdm none
```

---

## Quality Control

QC runs automatically alongside the main pipeline (controlled by `qc.enabled` in config). The final output is a single MultiQC HTML report aggregating all metrics.

### What gets collected

| Tool | Runs on | Metrics |
|---|---|---|
| **FastQC** | Raw FASTQs (+ trimmed if enabled) | Per-base quality, adapter content, GC distribution, overrepresented sequences |
| **samtools stats** | Final BAM | Error rate per cycle, base composition, insert size, indel distribution |
| **samtools flagstat** | Final BAM | Total/mapped/paired/duplicate read counts |
| **Picard CollectMultipleMetrics** | Final BAM | Alignment summary, insert size, GC bias, quality by cycle, sequencing artifacts |
| **Picard MarkDuplicates** | (already produced) | Duplication rate, optical duplicates |
| **Qualimap** (opt-in) | Final BAM | Coverage distribution, mean depth, mapping quality distribution |

### Configuration

```yaml
qc:
  enabled: true                 # master switch
  fastqc: true                  # FastQC on FASTQs
  samtools_stats: true          # samtools stats on final BAM
  samtools_flagstat: true       # samtools flagstat on final BAM
  picard_collect_metrics: true  # Picard CollectMultipleMetrics
  qualimap: false               # opt-in (requires ~16 GB RAM on GRCh38)
  qualimap_feature_file: ""     # BED for exome on-target stats
```

Set `qc.enabled: false` to skip all QC and produce only BAM files (pipeline behaves exactly as before).

### Output

```
{output_folder}/qc/
├── fastqc/raw/          # FastQC reports per FASTQ file
├── fastqc/trimmed/      # (when trimming enabled)
├── samtools/            # *.stats.txt, *.flagstat.txt
├── picard/              # CollectMultipleMetrics output files
├── qualimap/            # (when qc.qualimap: true)
├── multiqc_report.html  # ← aggregated report
└── multiqc_report_data/
```

---

## Deprecated Files

Old standalone workflows and SLURM launchers are in `deprecated/` for reference during migration. See `deprecated/README.md` for details. Scheduled for removal August 2025.

## License

MIT
