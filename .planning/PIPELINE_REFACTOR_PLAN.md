# Plan: Unified Pipeline Refactor + Config Generator

## Problems to Fix

### 1. Trimming is disconnected from the main DAG

`rule all` requests BQSR BAMs but never triggers trimming. `common.smk` always uses `trimmed_r1_suffix` regardless of `trimming.enabled`. The trim rule uses a `{sample}` wildcard while alignment uses `{basename}` — fragile implicit coupling.

### 2. No way to generate config from sequencing facility deliverables

Users receive a folder like:
```
250903_LH00253_0332_B232J72LT4_A5297_FASTQ/
├── SampleSheet.csv
├── A5297_DNA_01_STREAM_P1_L1_S1_L008_R1_001.fastq.gz
├── A5297_DNA_01_STREAM_P1_L1_S1_L008_R2_001.fastq.gz
├── A5297_DNA_02_STREAM_P1_N1_S2_L008_R1_001.fastq.gz
├── A5297_DNA_02_STREAM_P1_N1_S2_L008_R2_001.fastq.gz
├── *.md5
├── SampleSheet.csv
└── *_multiqc_report.html
```

SampleSheet.csv format (minimal, from BIH/Charité core facility):
```
8,A5297_DNA_01_STREAM_P1_L1,GTTATCGA,ACTACTTC,Project
8,A5297_DNA_02_STREAM_P1_N1,ATAGTGAC,TGTATCGA,Project
```

Columns: `Lane, Sample_Name, index_i7, index_i5, Sample_Project`

Users must currently hand-build `config/samples.tsv` from this. Error-prone and tedious.

### 3. FASTQ input path is rigid

The pipeline assumes all FASTQs are in one flat directory (`paths.fastq_folder`). Real deliverables may span multiple run folders or subfolders.

---

## Refactoring Plan

### Phase 1: Fix conditional trimming in the DAG

**Goal**: `rule all` drives the entire pipeline. When `trimming.enabled: true`, trimming runs automatically before alignment. When `false`, alignment reads raw FASTQs directly.

**Approach** (input function pattern from snakemake-workflows/dna-seq-gatk-variant-calling):

#### 1a. Add input function in `common.smk`

```python
TRIMMING_ENABLED = config.get("trimming", {}).get("enabled", False)

def get_fastq_r1(wildcards):
    """Return R1 FASTQ path — trimmed or raw depending on config."""
    if TRIMMING_ENABLED:
        return os.path.join(TRIMMED_DIR, f"{wildcards.basename}{TRIMMED_R1_SUFFIX}")
    return os.path.join(FASTQ_DIR, f"{wildcards.basename}{R1_SUFFIX}")

def get_fastq_r2(wildcards):
    """Return R2 FASTQ path — trimmed or raw depending on config."""
    if TRIMMING_ENABLED:
        return os.path.join(TRIMMED_DIR, f"{wildcards.basename}{TRIMMED_R2_SUFFIX}")
    return os.path.join(FASTQ_DIR, f"{wildcards.basename}{R2_SUFFIX}")
```

#### 1b. Define both raw and trimmed suffixes in `common.smk`

```python
# Raw FASTQ suffixes
R1_SUFFIX = config.get("fastq", {}).get("r1_suffix", "_R1_001.fastq.gz")
R2_SUFFIX = config.get("fastq", {}).get("r2_suffix", "_R2_001.fastq.gz")

# Trimmed FASTQ suffixes
TRIMMED_R1_SUFFIX = config.get("fastq", {}).get("trimmed_r1_suffix", ".bbduk_R1_001.fastq.gz")
TRIMMED_R2_SUFFIX = config.get("fastq", {}).get("trimmed_r2_suffix", ".bbduk_R2_001.fastq.gz")

TRIMMED_DIR = os.path.join(OUTPUT_DIR, "bbduk_trimmed")
```

#### 1c. Update `alignment.smk` to use input functions

```python
rule bwa_map:
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2,
    ...
```

#### 1d. Align trim rule wildcards with alignment wildcards

Rename `{sample}` → `{basename}` in `trim.smk` so the wildcard matches. Update `_get_trim_samples()` to use `samples_df["fastq_files_basename"]` instead of filesystem globbing. This way the samples.tsv is the single source of truth for both trimming and alignment.

```python
rule trim_adapters:
    input:
        r1=lambda wc: os.path.join(FASTQ_DIR, f"{wc.basename}{R1_SUFFIX}"),
        r2=lambda wc: os.path.join(FASTQ_DIR, f"{wc.basename}{R2_SUFFIX}"),
    output:
        trimmed_r1=os.path.join(TRIMMED_DIR, "{basename}" + TRIMMED_R1_SUFFIX),
        trimmed_r2=os.path.join(TRIMMED_DIR, "{basename}" + TRIMMED_R2_SUFFIX),
    ...
```

#### 1e. Remove `trim_all` as separate target

With the DAG now chaining through input functions, `trim_all` is unnecessary. Snakemake will auto-trigger `trim_adapters` when `bwa_map` needs trimmed FASTQs.

**Files changed**: `common.smk`, `alignment.smk`, `trim.smk`

---

### Phase 2: Config generator script (`scripts/generate_config.py`)

**Goal**: Parse a SampleSheet.csv (+ optionally scan for FASTQs) and generate `config/samples.tsv` with correct column values. Minimal dependencies (stdlib + pandas, already available).

#### 2a. SampleSheet.csv parsing

Handle the BIH/Charité minimal format (no section headers, just CSV data rows):
```
Lane,Sample_Name,index_i7,index_i5,Sample_Project
```

Also handle standard Illumina format with `[Header]`, `[Data]` sections (skip lines until data rows).

```python
def parse_samplesheet(path):
    """Parse Illumina SampleSheet.csv, return list of dicts.

    Handles both:
    - Full format with [Header]/[Data] sections
    - Minimal format (bare CSV rows)
    """
```

#### 2b. FASTQ discovery and basename construction

Given SampleSheet data + a FASTQ directory, match each sample to its actual FASTQ files:

```python
ILLUMINA_RE = re.compile(
    r'^(?P<sample>.+?)_S(?P<snum>\d+)_L(?P<lane>\d{3})_(?P<read>R[12])_001\.fastq\.gz$'
)
```

For the example data:
- SampleSheet row: `8, A5297_DNA_01_STREAM_P1_L1, GTTATCGA, ACTACTTC, Project`
- FASTQ file: `A5297_DNA_01_STREAM_P1_L1_S1_L008_R1_001.fastq.gz`
- Extracted basename: `A5297_DNA_01_STREAM_P1_L1_S1_L008`
- Extracted lane: `L008`
- project_sample: `A5297_DNA_01_STREAM_P1_L1` (from SampleSheet Sample_Name)

#### 2c. Generate samples.tsv

Map SampleSheet fields to pipeline columns:

| SampleSheet field | samples.tsv column | Derivation |
|---|---|---|
| — | `fastq_files_basename` | FASTQ filename before `_R1_001.fastq.gz` |
| Lane | `lane` | `L{lane:03d}` from SampleSheet or parsed from filename |
| Sample_Name | `project_sample` | Direct mapping |
| Sample_Project or folder name | `mdc_project` | From SampleSheet or CLI arg `--project` |

#### 2d. CLI interface

```
usage: generate_config.py [-h] --fastq-dir DIR [--samplesheet CSV]
                          [--project NAME] [--output PATH]

Generate config/samples.tsv from sequencing facility deliverables.

Required:
  --fastq-dir DIR        Directory containing FASTQ files

Optional:
  --samplesheet CSV      Path to SampleSheet.csv (auto-detected if in fastq-dir)
  --project NAME         Project identifier for mdc_project column
                         (default: inferred from folder name, e.g. A5297)
  --output PATH          Output samples.tsv path (default: config/samples.tsv)
  --config-template      Also generate a config.yaml template
  --dry-run              Show what would be generated without writing files
```

#### 2e. Example usage

```bash
# Simplest: point at the delivery folder
python scripts/generate_config.py \
    --fastq-dir /charite-store-f/f-bih-gp/processed/A5297/250903_LH00253_0332_B232J72LT4_A5297_FASTQ/

# Output:
# ✓ Found SampleSheet.csv with 2 samples
# ✓ Matched 4 FASTQ files (2 pairs)
#
# fastq_files_basename                    lane   project_sample              mdc_project
# A5297_DNA_01_STREAM_P1_L1_S1_L008      L008   A5297_DNA_01_STREAM_P1_L1   A5297
# A5297_DNA_02_STREAM_P1_N1_S2_L008      L008   A5297_DNA_02_STREAM_P1_N1   A5297
#
# ✓ Written: config/samples.tsv (2 samples)
```

**Files created**: `scripts/generate_config.py`

---

### Phase 3: Multi-folder FASTQ support

**Goal**: Support FASTQs spread across multiple delivery folders (e.g., re-sequencing runs) via a `subfolder` column in samples.tsv.

#### 3a. Add `subfolder` support to samples.tsv

The `subfolder` column already exists in the current samples.tsv template. Wire it into the input functions:

```python
def get_fastq_r1(wildcards):
    row = samples_df.loc[wildcards.basename]
    subfolder = row.get("subfolder", "")
    base_dir = TRIMMED_DIR if TRIMMING_ENABLED else FASTQ_DIR
    suffix = TRIMMED_R1_SUFFIX if TRIMMING_ENABLED else R1_SUFFIX
    if subfolder:
        return os.path.join(base_dir, subfolder, f"{wildcards.basename}{suffix}")
    return os.path.join(base_dir, f"{wildcards.basename}{suffix}")
```

#### 3b. Update schema

Add `subfolder` as optional field in `samples.schema.yaml`.

#### 3c. Update trim rule

If trimming is enabled, the trim rule input must also respect `subfolder` to find raw FASTQs in the correct delivery folder.

**Files changed**: `common.smk`, `trim.smk`, `samples.schema.yaml`

---

### Phase 4: Config template generator

**Goal**: `generate_config.py --config-template` also produces a starter `config/config.yaml` with paths pre-filled.

Scan for:
- Reference genome (look for `*.fna`, `*.fa`, `*.fasta` with BWA index files)
- Known-sites VCFs (look for `*dbsnp*`, `*gnomad*`, `*1000g*` with `.tbi` index)
- Pre-fill `paths.fastq_folder` from `--fastq-dir`
- Pre-fill `paths.output_folder` based on project convention

This builds on the discovery logic from `.planning/CONFIG_HELPER_PLAN.md` but as a non-interactive script (suitable for HPC where interactive terminals are limited).

**Files changed**: `scripts/generate_config.py` (extend)

---

## Implementation Order

| Phase | Scope | Files changed |
|---|---|---|
| 1 | Fix conditional trimming DAG | `common.smk`, `alignment.smk`, `trim.smk` |
| 2 | Config generator script | `scripts/generate_config.py` (new) |
| 3 | Multi-folder FASTQ support | `common.smk`, `trim.smk`, `samples.schema.yaml` |
| 4 | Config template generation | `scripts/generate_config.py` (extend) |

Phases 1 and 2 are independent and can be done in parallel. Phase 3 depends on Phase 1. Phase 4 depends on Phase 2.

---

## Validation Checklist

After implementation, verify:

- [ ] `trimming.enabled: false` → pipeline reads raw FASTQs, no trim rule fires
- [ ] `trimming.enabled: true` → trim rule fires, alignment reads trimmed FASTQs
- [ ] `snakemake -n` shows correct DAG for both modes
- [ ] `generate_config.py` correctly parses BIH/Charité SampleSheet.csv format
- [ ] `generate_config.py` correctly parses standard Illumina SampleSheet.csv with `[Data]` header
- [ ] Generated `samples.tsv` validates against `samples.schema.yaml`
- [ ] `subfolder` column works for multi-folder inputs
- [ ] Dry-run succeeds with generated config against real data paths
