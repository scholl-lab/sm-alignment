# Refactoring Plan: Configurable sm-alignment Pipeline

This plan addresses all hardcoded values identified in `HARDCODED_VALUES_AUDIT.md` and aligns the pipeline with current Snakemake best practices (v8+/v9).

---

## 1. Target Directory Structure

Migrate from the current flat layout to the [official Snakemake standardized structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html):

```
sm-alignment/
├── workflow/
│   ├── Snakefile                     # Main entry point (replaces alignment_pipeline.smk)
│   ├── rules/
│   │   ├── common.smk               # Shared helpers: get_samples(), metadata loading
│   │   ├── alignment.smk            # bwa_map rule
│   │   ├── merge.smk                # merge_bam_files rule
│   │   ├── dedup.smk                # deduplicate_bam_files rule
│   │   ├── bqsr.smk                 # base_recalibration + apply_bqsr rules
│   │   ├── trim.smk                 # BBDuk trimming rule
│   │   └── utilities.smk            # subset_bam, bam_to_fastq, md5sum
│   ├── envs/
│   │   ├── bwa_samtools.yaml         # bwa + samtools with pinned versions
│   │   ├── gatk.yaml                 # gatk4 with pinned version
│   │   └── bbtools.yaml              # bbmap/bbtools with pinned version
│   └── schemas/
│       ├── config.schema.yaml        # JSON Schema for config validation
│       └── samples.schema.yaml       # JSON Schema for metadata TSV validation
├── config/
│   ├── config.yaml                   # Main experiment/pipeline config
│   └── samples.tsv                   # Sample metadata (was metadata/metadata_exomes.tsv)
├── profiles/
│   └── default/
│       └── config.yaml               # Workflow profile: resource overrides, SLURM settings
├── scripts/
│   └── run_snakemake.sh              # Single generic SLURM launcher (replaces 11 scripts)
└── resources/                        # Static resources (BED files, adapter refs, etc.)
```

The old `scripts/snakemake/` and individual `scripts/run_*.sh`/`scripts/submit_*.sh` files are retired.

---

## 2. Unified Configuration File

Replace the current scattered config files (`config.yaml`, `configs/config_alignment.yaml`, `configs/config_trim_adapters.yaml`) with a single hierarchical `config/config.yaml`.

Design follows the pattern used by [snakemake-workflows/dna-seq-gatk-variant-calling](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling):

```yaml
# =============================================================================
# config/config.yaml — Pipeline configuration
# =============================================================================

# --- Reference genome & known sites ---
ref:
  genome: "analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  genome_gz: "analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
  build: "GRCh38"
  known_sites:
    - "analysis/GATK_resource_bundle/af-only-gnomad.hg38.vcf.gz"
    - "analysis/GATK_resource_bundle/af-only-gnomad.hg38.common_biallelic.vcf.gz"
    - "analysis/GATK_resource_bundle/1000g_pon.hg38.vcf.gz"

# --- Paths ---
paths:
  samples: "config/samples.tsv"
  fastq_folder: "results/exomes/bbduk_trimmed"
  output_folder: "results/exomes"
  log_subdir: "logs"

# --- Read group defaults ---
read_group:
  platform: "ILLUMINA"

# --- FASTQ file patterns ---
fastq:
  r1_suffix: "_R1_001.fastq.gz"
  r2_suffix: "_R2_001.fastq.gz"
  trimmed_r1_suffix: ".bbduk_R1_001.fastq.gz"
  trimmed_r2_suffix: ".bbduk_R2_001.fastq.gz"

# --- BAM file naming ---
bam:
  merged_suffix: ".merged.bam"
  dedup_suffix: ".merged.dedup.bam"
  dedup_metrics_suffix: ".merged.dedup_metrics.txt"
  recal_table_suffix: ".merged.dedup.recal_data.table"
  final_suffix: ".merged.dedup.bqsr.bam"

# --- Processing options ---
processing:
  remove_duplicates: true
  compression_level: 6

# --- Tool-specific extra CLI arguments (passthrough strings) ---
params:
  bwa_mem:
    extra: ""
  samtools:
    sort_extra: ""
    merge_extra: ""
  gatk:
    MarkDuplicates: "--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"
    BaseRecalibrator: ""
    ApplyBQSR: ""

# --- BBDuk trimming parameters ---
trimming:
  enabled: false
  bbduk_ref: "adapters,artifacts"
  ktrim: "r"
  k: 23
  mink: 11
  hdist: 1
  tpe: "t"
  tbo: "t"
  ftl: 5
  trimpolyg: 3
  trimpolya: 3
  qtrim: "t"
  trimq: 10
  ziplevel: 5

# --- Subset BAM (optional) ---
subset:
  bed_file: ""   # path to BED file; leave empty to skip
  output_suffix: ".subset.bam"
```

### What changed vs. current config

| Current | New | Rationale |
|---------|-----|-----------|
| Flat keys (`fastq_folder`, `aligned_folder`, ...) | Hierarchical (`ref`, `paths`, `params`, ...) | Organized by concern, follows catalog conventions |
| Reference paths only in config, but known-sites hardcoded in `bqsr_bams.smk` shell | `ref.known_sites` list in config, iterated in rule | Eliminates critical hardcoding bug |
| `ILLUMINA` platform hardcoded in 4 files | `read_group.platform` in config | Configurable |
| FASTQ suffixes hardcoded in 6 files | `fastq.r1_suffix` etc. in config | One place to change |
| BAM suffixes hardcoded in 8 files | `bam.*_suffix` in config | One place to change |
| Java heap `-Xms4000m -Xmx7g` in 6 shell blocks | Derived from `resources.mem_mb` (see Section 4) | Consistent with allocated memory |
| GATK flags (`--CREATE_INDEX`, `--VALIDATION_STRINGENCY`) hardcoded | `params.gatk.MarkDuplicates` passthrough string | User-configurable without touching Snakefile |
| Absolute BED path in `subset_bam.smk` | `subset.bed_file` in config | Portable |
| BBDuk params split across `config_trim_adapters.yaml` | `trimming.*` section in unified config | Single config file |

---

## 3. Conda Environment YAML Files

Replace named environment references (`conda: "base"`, `conda: "gatk"`) with portable YAML env files with pinned versions. Named envs break portability — anyone cloning the repo would need to manually recreate them.

Reference: [Snakemake Deployment docs](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)

### `workflow/envs/bwa_samtools.yaml`

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bwa=0.7.18
  - samtools=1.21
```

### `workflow/envs/gatk.yaml`

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - gatk4=4.6.1.0
  - samtools=1.21
```

### `workflow/envs/bbtools.yaml`

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bbmap=39.06
```

### Usage in rules

```python
rule bwa_map:
    conda: "../envs/bwa_samtools.yaml"   # relative to Snakefile

rule deduplicate_bam_files:
    conda: "../envs/gatk.yaml"
```

### Version pinning

After initial setup, generate platform-specific pin files for exact reproducibility:

```bash
snakedeploy pin-conda-envs workflow/envs/bwa_samtools.yaml
# produces workflow/envs/bwa_samtools.linux-64.pin.txt
```

---

## 4. Workflow Profile for Resources

Move all thread counts, memory values, time limits, and SLURM settings out of rule definitions and into a workflow profile at `profiles/default/config.yaml`.

This is the Snakemake-recommended approach — rules keep sensible defaults, but the profile overrides them per deployment. See [Snakemake CLI docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

### `profiles/default/config.yaml`

```yaml
# =============================================================================
# Workflow Profile — Resource allocation & SLURM settings
# =============================================================================
# Override precedence: CLI > workflow profile > rule defaults > default-resources
# =============================================================================

# --- Execution settings ---
latency-wait: 60
jobs: 20
use-conda: true
rerun-incomplete: true
printshellcmds: true

# --- Default resources (baseline for ALL rules) ---
default-resources:
  mem_mb: 4000
  runtime: 120            # 2 hours in minutes
  tmpdir: "system_tmpdir" # inherit from SLURM node $TMPDIR

# --- Per-rule thread counts ---
set-threads:
  bwa_map: 16
  merge_bam_files: 8
  deduplicate_bam_files: 4
  base_recalibration: 4
  apply_bqsr: 4
  trim_adapters: 8
  convert_bam_to_fastq: 8
  samtools_view_subset: 1

# --- Per-rule resource overrides ---
set-resources:
  bwa_map:
    mem_mb: 19200          # 16 * 1200
    runtime: 1440          # 24h
  merge_bam_files:
    mem_mb: 9600           # 8 * 1200
    runtime: 1440          # 24h
  deduplicate_bam_files:
    mem_mb: 17600          # 4 * 4400
    runtime: 4320          # 72h
  base_recalibration:
    mem_mb: 17600          # 4 * 4400
    runtime: 4320          # 72h
  apply_bqsr:
    mem_mb: 17600          # 4 * 4400
    runtime: 4320          # 72h
  trim_adapters:
    mem_mb: 8000           # 8 * 1000
    runtime: 720           # 12h
```

### What this replaces

- All `threads: N` hardcoded values in rules — rules can still set defaults, but the profile overrides
- All `get_mem_from_threads()` helper functions — memory is set directly per-rule
- All `time = "72:00:00"` resource values — use `runtime` in minutes (Snakemake standard)
- All `SCRATCH_DIR` / `TMPDIR` handling in rules — `default-resources.tmpdir` handles it

### CLI overrides for one-off runs

```bash
# Give BWA more threads on a large-memory node
snakemake --set-threads bwa_map=32 --set-resources bwa_map:mem_mb=38400
```

---

## 5. Java Heap Size: Derive from resources.mem_mb

Instead of hardcoding `-Xms4000m -Xmx7g` in every GATK rule's shell block, derive the JVM heap from the allocated memory. This follows the pattern used by [snakemake-wrapper-utils](https://github.com/snakemake/snakemake-wrapper-utils/blob/master/snakemake_wrapper_utils/java.py).

### Helper function in `workflow/rules/common.smk`

```python
def get_java_opts(wildcards, resources):
    """
    Generate GATK --java-options string from allocated resources.
    Uses 80% of mem_mb for -Xmx (20% reserved for JVM overhead).
    """
    xmx = int(resources.mem_mb * 0.8)
    xms = int(resources.mem_mb * 0.2)
    tmpdir = resources.tmpdir
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={tmpdir}"
```

### Usage in GATK rules

```python
rule deduplicate_bam_files:
    params:
        java_opts = get_java_opts,
        extra = config["params"]["gatk"]["MarkDuplicates"],
    shell:
        r"""
        gatk --java-options '{params.java_opts}' MarkDuplicates \
            {params.extra} \
            -I "{input.merged_bam}" \
            -O "{output.dedup_bam}" \
            -M "{output.metrics}" \
            2> {log}
        """
```

This means changing `mem_mb` in the profile automatically adjusts JVM heap size — no disconnect between SLURM-allocated memory and what Java actually uses.

---

## 6. Single Generic SLURM Launcher

Replace all 11 submission scripts with one parameterized `scripts/run_snakemake.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=sm_pipeline
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log
#
# Usage:
#   sbatch scripts/run_snakemake.sh <SNAKEFILE> [CONFIG] [PROFILE] [EXTRA_ARGS...]
#
# Examples:
#   sbatch scripts/run_snakemake.sh workflow/Snakefile
#   sbatch --job-name=sm_trim scripts/run_snakemake.sh workflow/Snakefile config/config_trim.yaml
#   sbatch scripts/run_snakemake.sh workflow/Snakefile config/config.yaml cubi-v1 --forceall

set -euo pipefail

SNAKEFILE="${1:?Error: SNAKEFILE path required}"
CONFIG_FILE="${2:-config/config.yaml}"
PROFILE="${3:-cubi-v1}"
shift 3 2>/dev/null || true

# TMPDIR setup
export TMPDIR="${HOME}/scratch/tmp"
mkdir -p "${TMPDIR}"
export TMPDIR=$(mktemp -d "${TMPDIR}/sm.XXXXXX")
trap 'rm -rf "${TMPDIR}"' EXIT

# Logging
mkdir -p slurm_logs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

echo "=== Snakemake Launch ==="
echo "  Snakefile:  ${SNAKEFILE}"
echo "  Config:     ${CONFIG_FILE}"
echo "  Profile:    ${PROFILE}"
echo "  TMPDIR:     ${TMPDIR}"
echo "  Extra args: $*"
echo "  Start:      $(date)"

srun snakemake \
    -s "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile profiles/default \
    --profile="${PROFILE}" \
    "$@"

echo "=== Finished: $(date) ==="
```

### What this replaces

| Old (11 scripts) | New (1 script) |
|---|---|
| `sbatch scripts/submit_alignment_pipeline.sh` | `sbatch scripts/run_snakemake.sh workflow/Snakefile` |
| `sbatch scripts/submit_trim_adapters.sh` | `sbatch --job-name=sm_trim scripts/run_snakemake.sh workflow/Snakefile config/config_trim.yaml` |
| `sbatch scripts/run_dedup_bams.sh` | No longer needed — integrated pipeline handles all stages |
| etc. | |

Key: the `--workflow-profile profiles/default` flag loads the resource profile from the repo automatically. The `--profile=cubi-v1` continues to use the BIH cluster's global profile for SLURM executor settings.

---

## 7. Config Schema Validation

Add JSON Schema files (written in YAML) to catch config errors at pipeline startup rather than mid-run.

### `workflow/schemas/config.schema.yaml`

```yaml
$schema: "https://json-schema.org/draft-04/schema#"
description: Configuration schema for sm-alignment pipeline
type: object
properties:
  ref:
    type: object
    properties:
      genome:
        type: string
      genome_gz:
        type: string
      build:
        type: string
        enum: ["GRCh37", "GRCh38"]
      known_sites:
        type: array
        items:
          type: string
        minItems: 1
    required: [genome, build, known_sites]

  paths:
    type: object
    properties:
      samples:
        type: string
      fastq_folder:
        type: string
      output_folder:
        type: string
      log_subdir:
        type: string
        default: "logs"
    required: [samples, fastq_folder, output_folder]

  read_group:
    type: object
    properties:
      platform:
        type: string
        default: "ILLUMINA"

  fastq:
    type: object
    properties:
      r1_suffix:
        type: string
        default: "_R1_001.fastq.gz"
      r2_suffix:
        type: string
        default: "_R2_001.fastq.gz"
      trimmed_r1_suffix:
        type: string
        default: ".bbduk_R1_001.fastq.gz"
      trimmed_r2_suffix:
        type: string
        default: ".bbduk_R2_001.fastq.gz"

  bam:
    type: object
    properties:
      final_suffix:
        type: string
        default: ".merged.dedup.bqsr.bam"

  params:
    type: object
    properties:
      bwa_mem:
        type: object
        properties:
          extra:
            type: string
            default: ""
      gatk:
        type: object
        properties:
          MarkDuplicates:
            type: string
            default: "--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"
          BaseRecalibrator:
            type: string
            default: ""
          ApplyBQSR:
            type: string
            default: ""

required: [ref, paths]
```

### `workflow/schemas/samples.schema.yaml`

```yaml
$schema: "https://json-schema.org/draft-04/schema#"
description: Row schema for sample metadata TSV
properties:
  fastq_files_basename:
    type: string
  lane:
    type: string
  project_sample:
    type: string
  mdc_project:
    type: string
required: [fastq_files_basename, lane, project_sample]
```

### Usage in `workflow/Snakefile`

```python
import pandas as pd
from snakemake.utils import validate

configfile: "config/config.yaml"
validate(config, "schemas/config.schema.yaml")

samples = pd.read_table(config["paths"]["samples"])
validate(samples, "schemas/samples.schema.yaml")
samples = samples.set_index("fastq_files_basename", drop=False)
```

---

## 8. Refactored Rule Examples

### `workflow/rules/common.smk` — Shared helpers

```python
import os
import pandas as pd

# --- Config shortcuts ---
REF         = config["ref"]["genome"]
REF_GZ      = config["ref"].get("genome_gz", "")
KNOWN_SITES = config["ref"]["known_sites"]
FASTQ_DIR   = config["paths"]["fastq_folder"]
OUTPUT_DIR  = config["paths"]["output_folder"]
LOG_SUBDIR  = config["paths"].get("log_subdir", "logs")
LOG_DIR     = os.path.join(OUTPUT_DIR, LOG_SUBDIR)
PLATFORM    = config.get("read_group", {}).get("platform", "ILLUMINA")
R1_SUFFIX   = config.get("fastq", {}).get("trimmed_r1_suffix", ".bbduk_R1_001.fastq.gz")
R2_SUFFIX   = config.get("fastq", {}).get("trimmed_r2_suffix", ".bbduk_R2_001.fastq.gz")

# --- Load and validate samples ---
samples_df = pd.read_table(config["paths"]["samples"]).set_index("fastq_files_basename", drop=False)

def get_samples():
    return sorted(samples_df["project_sample"].unique().tolist())

def get_java_opts(wildcards, resources):
    xmx = int(resources.mem_mb * 0.8)
    xms = int(resources.mem_mb * 0.2)
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={resources.tmpdir}"
```

### `workflow/rules/alignment.smk` — BWA alignment

```python
rule bwa_map:
    input:
        r1 = lambda wc: os.path.join(FASTQ_DIR, f"{wc.basename}{R1_SUFFIX}"),
        r2 = lambda wc: os.path.join(FASTQ_DIR, f"{wc.basename}{R2_SUFFIX}"),
    output:
        bam = temp(os.path.join(OUTPUT_DIR, "aligned", "{basename}.bam")),
    params:
        reference = REF_GZ,
        read_group = lambda wc: (
            '"@RG\\tID:{lane}-{sample}\\tSM:{sample}\\tLB:{sample}'
            '\\tPL:{platform}\\tPU:{lane}-{project}"'
            .format(
                lane=samples_df.loc[wc.basename, "lane"],
                sample=samples_df.loc[wc.basename, "project_sample"],
                project=samples_df.loc[wc.basename, "mdc_project"],
                platform=PLATFORM,
            )
        ),
        sort_threads = 2,
        extra = config.get("params", {}).get("bwa_mem", {}).get("extra", ""),
    threads: 16
    resources:
        runtime = 1440,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        bwa      = os.path.join(LOG_DIR, "map.bwa.{basename}.log"),
        samtools = os.path.join(LOG_DIR, "map.samtools.{basename}.log"),
    shell:
        r"""
        BWA_THREADS=$(({{threads}} - {{params.sort_threads}}))
        SORT_MEM=$(({{params.sort_threads}} * 2000))

        TMP_SORT_DIR=$(mktemp -p {{resources.tmpdir}} -d samtools-sort.XXXXXX)

        bwa mem -t $BWA_THREADS \
            {{params.extra}} \
            -R {{params.read_group}} \
            {{params.reference}} \
            {{input.r1}} {{input.r2}} \
            2> {{log.bwa}} \
        | samtools sort -@ {{params.sort_threads}} \
            -m ${{SORT_MEM}}M \
            -O BAM \
            -T "$TMP_SORT_DIR"/tmp \
            -o {{output.bam}} \
            2> {{log.samtools}}

        rm -rf "$TMP_SORT_DIR"
        """
```

### `workflow/rules/bqsr.smk` — BQSR (known-sites from config list)

```python
rule base_recalibration:
    input:
        dedup_bam = os.path.join(OUTPUT_DIR, "dedup", "{sample}.merged.dedup.bam"),
    output:
        recal_table = os.path.join(OUTPUT_DIR, "bqsr", "{sample}.merged.dedup.recal_data.table"),
    params:
        java_opts   = get_java_opts,
        known_sites = lambda wc: " ".join(
            f'--known-sites "{ks}"' for ks in KNOWN_SITES
        ),
        extra       = config.get("params", {}).get("gatk", {}).get("BaseRecalibrator", ""),
    threads: 4
    resources:
        runtime = 4320,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "recal.gatk.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' BaseRecalibrator \
            {params.extra} \
            -I "{input.dedup_bam}" \
            -R "{REF}" \
            {params.known_sites} \
            -O "{output.recal_table}" \
            2> {log}
        """
```

This eliminates the hardcoded known-sites paths from `bqsr_bams.smk` — the list is iterated from `config["ref"]["known_sites"]`.

---

## 9. Implementation Phases

### Phase 1: Foundation (non-breaking)

1. Create `workflow/envs/` with YAML env files (pin current tool versions)
2. Create `workflow/schemas/` with config and samples schemas
3. Create `config/config.yaml` with the new hierarchical structure
4. Create `profiles/default/config.yaml` with all resource overrides
5. Create `scripts/run_snakemake.sh` generic launcher

### Phase 2: Migrate integrated pipeline

6. Create `workflow/Snakefile` as the main entry point
7. Extract `workflow/rules/common.smk` with shared helpers
8. Migrate `alignment_pipeline.smk` rules into modular `workflow/rules/*.smk` files
9. Replace all hardcoded values with config references (per sections 2-5 above)
10. Replace `conda: "base"` / `conda: "gatk"` with YAML env file paths
11. Replace `time = "HH:MM:SS"` with `runtime = <minutes>` in all rules
12. Replace hardcoded Java opts with `get_java_opts()` helper
13. Move tool-internal params (`sort_threads`, `bwa_threads`) to `params:` block

### Phase 3: Migrate standalone workflows

14. Migrate `trim_adapters.smk` into `workflow/rules/trim.smk`
15. Migrate utility workflows (`subset_bam.smk`, `bam_to_fastq.smk`, `md5sum_files.smk`) into `workflow/rules/utilities.smk`
16. Fix `bqsr_bams.smk` hardcoded reference/known-sites paths (or retire in favor of integrated pipeline)
17. Fix `merge_bams.smk`, `dedup_bams.smk` hardcoded `results/` paths (or retire)

### Phase 4: Cleanup

18. Add config validation (`validate()`) to `workflow/Snakefile`
19. Move `config/samples.tsv` from `metadata/metadata_exomes.tsv`, update config path
20. Move static resources (BED files) to `resources/` and update config paths
21. Deprecate old `scripts/run_*.sh` and `scripts/submit_*.sh` (keep in git history)
22. Retire old standalone `.smk` files that are fully superseded by the integrated pipeline
23. Update `README.md` with new usage instructions
24. Dry-run test: `snakemake -s workflow/Snakefile --configfile config/config.yaml -n`

---

## 10. Migration Cheat Sheet

Quick reference for converting each category of hardcoded value:

| Hardcoded value | Where it moves | Access pattern |
|---|---|---|
| `conda: "base"` | `workflow/envs/bwa_samtools.yaml` | `conda: "../envs/bwa_samtools.yaml"` |
| `conda: "gatk"` | `workflow/envs/gatk.yaml` | `conda: "../envs/gatk.yaml"` |
| `threads: 16` | Keep as default in rule | Override via `profiles/default/config.yaml` `set-threads:` |
| `mem_mb = threads * 1200` | `profiles/default/config.yaml` `set-resources:` | `resources.mem_mb` (auto from profile) |
| `time = "24:00:00"` | `profiles/default/config.yaml` `set-resources:` | `resources.runtime` (integer minutes) |
| `-Xms4000m -Xmx7g` | Derived from `resources.mem_mb` | `get_java_opts()` helper in `params:` |
| `-Dsamjdk.compression_level=6` | `config["processing"]["compression_level"]` | `params:` block in rule |
| `--VALIDATION_STRINGENCY SILENT` | `config["params"]["gatk"]["MarkDuplicates"]` | `{params.extra}` in shell |
| `ILLUMINA` (platform) | `config["read_group"]["platform"]` | Variable `PLATFORM` in `common.smk` |
| `_R1_001.fastq.gz` | `config["fastq"]["r1_suffix"]` | Variable `R1_SUFFIX` in `common.smk` |
| `.merged.dedup.bqsr.bam` | `config["bam"]["final_suffix"]` | Variable from config |
| Known-sites paths in shell | `config["ref"]["known_sites"]` list | Iterated in `params:` lambda |
| Absolute BED file path | `config["subset"]["bed_file"]` | Input from config |
| `$HOME/scratch/tmp` | `scripts/run_snakemake.sh` + `default-resources.tmpdir` | `{resources.tmpdir}` |
| `cubi-v1` profile | `scripts/run_snakemake.sh` arg `$3` | `--profile=${PROFILE}` |
| `slurm_logs` | `scripts/run_snakemake.sh` | Centralized in launcher |
| `'metadata.tsv'` | `config["paths"]["samples"]` | Loaded in `common.smk` |
| `configfile: "config.yaml"` | `workflow/Snakefile` + `--configfile` CLI | Single entry point |
| `results/merged`, `results/dedup` etc. | Derived from `config["paths"]["output_folder"]` | `os.path.join(OUTPUT_DIR, "merged")` |

---

## 11. Compatibility Notes

- **Snakemake version**: The `runtime` resource (integer minutes), `--workflow-profile`, and `set-threads`/`set-resources` in profiles require Snakemake >= 8. Verify the cluster has SM 8+ or coordinate with HPC admins.
- **cubi-v1 profile**: The BIH cluster's `cubi-v1` profile remains as the global executor profile. The new `profiles/default/` is a *workflow* profile layered on top — they are complementary, not conflicting.
- **Backward compatibility**: The old standalone workflows (`alignment.smk`, `merge_bams.smk`, etc.) can be kept during migration but should not be actively maintained. The integrated pipeline in `workflow/Snakefile` supersedes them.
- **Conda env creation**: First run with new YAML env files will take time to create environments. Subsequent runs use cached envs. Consider running `snakemake --conda-create-envs-only` before the first production run.

---

## Sources

- [Snakemake Best Practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)
- [Snakemake Distribution & Reproducibility (standard directory layout)](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
- [Snakemake Configuration & Schema Validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)
- [Snakemake CLI: set-threads, set-resources, default-resources](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
- [Snakemake SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [snakemake-workflows/dna-seq-gatk-variant-calling](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling) — Reference GATK pipeline
- [snakemake-wrapper-utils: get_java_opts](https://github.com/snakemake/snakemake-wrapper-utils/blob/master/snakemake_wrapper_utils/java.py) — Canonical Java heap derivation pattern
- [smk-simple-slurm](https://github.com/jdblischak/smk-simple-slurm) — Minimal SLURM profile example
