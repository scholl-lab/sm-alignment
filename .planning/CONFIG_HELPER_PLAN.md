# Plan: Interactive Config Setup Helper

A Python script that auto-discovers files from standardized directory layouts and interactively guides the user through generating `config/config.yaml` and `config/samples.tsv`.

---

## Inspiration

The approach is modeled after **nf-core launch** — a schema-driven interactive wizard — but enhanced with filesystem auto-discovery, which nf-core does not do. The key insight: most bioinformatics projects follow predictable directory structures, so the tool can propose sensible defaults by scanning the filesystem before asking questions.

---

## User Experience Flow

```
$ python scripts/setup_config.py

=== sm-alignment Config Setup ===

[Phase 1: Auto-Discovery]

Scanning for reference genomes...
  Found: analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
         BWA index: Yes | FAI: Yes | Dict: Yes

Scanning for known variant sites...
  Found 3 VCFs in analysis/GATK_resource_bundle/:
  ┌──────────────────────────────────────────────────┬───────┐
  │ VCF Path                                         │ Index │
  ├──────────────────────────────────────────────────┼───────┤
  │ af-only-gnomad.hg38.vcf.gz                       │  Yes  │
  │ af-only-gnomad.hg38.common_biallelic.vcf.gz      │  Yes  │
  │ 1000g_pon.hg38.vcf.gz                            │  Yes  │
  └──────────────────────────────────────────────────┴───────┘

Scanning for FASTQ files...
  Found 24 files → 6 samples across 2 lanes in download/exomes/
  ┌─────────────────────┬──────────┬───────┬──────┐
  │ Basename            │ Sample   │ Lane  │ R2?  │
  ├─────────────────────┼──────────┼───────┼──────┤
  │ SampleA_S1_L001     │ SampleA  │ L001  │ Yes  │
  │ SampleA_S1_L002     │ SampleA  │ L002  │ Yes  │
  │ SampleB_S2_L001     │ SampleB  │ L001  │ Yes  │
  │ ...                 │ ...      │ ...   │ ...  │
  └─────────────────────┴──────────┴───────┴──────┘

[Phase 2: Interactive Prompts]

? Reference genome build: (Use arrow keys)
  » GRCh38
    GRCh37

? Use discovered reference genome? (Y/n) Y

? Select known variant sites for BQSR: (checkbox, space to toggle)
  » [x] af-only-gnomad.hg38.vcf.gz
    [x] af-only-gnomad.hg38.common_biallelic.vcf.gz
    [x] 1000g_pon.hg38.vcf.gz

? Output directory: (results/exomes) █

? Enable BBDuk adapter trimming? (y/N) N

? Subset BAM by BED regions? (y/N) N

[Phase 3: Generate Files]

✓ Written: config/config.yaml
✓ Written: config/samples.tsv (6 samples)
✓ Validated against schema

Done! Run with: sbatch scripts/run_snakemake.sh workflow/Snakefile
```

---

## Architecture

### Single-file script (no package install required)

```
scripts/
  setup_config.py       # Self-contained helper (~400 lines)
```

Dependencies are all available in a standard bioinformatics conda env or pip-installable:

| Package | Purpose | Already available? |
|---|---|---|
| `questionary` | Interactive prompts with path autocomplete, checkboxes | New (pip install) |
| `rich` | Formatted tables, colored output | New (pip install) |
| `PyYAML` | YAML generation | Yes (Snakemake dep) |
| `pandas` | samples.tsv generation | Yes (Snakemake dep) |

### Module-level structure (within the single file)

```python
# --- Auto-discovery functions ---
discover_reference_genomes(search_dirs) -> list[dict]
discover_known_sites(search_dirs) -> list[dict]
discover_fastq_files(search_dirs) -> list[dict]
parse_illumina_filename(filename) -> dict | None

# --- Display functions (rich) ---
display_reference_genomes(candidates)
display_known_sites(vcfs)
display_samples(samples)

# --- Interactive prompt functions (questionary) ---
prompt_reference(discovered) -> dict
prompt_known_sites(discovered) -> list[str]
prompt_paths(discovered_fastq_dir) -> dict
prompt_optional_features() -> dict
prompt_trimming() -> dict

# --- Generation functions ---
generate_config_yaml(answers, output_path)
generate_samples_tsv(samples_df, output_path)
validate_outputs(config_path, samples_path)

# --- Entry point ---
main()
```

---

## Auto-Discovery Details

### Reference genome scanning

Search directories (in order): `analysis/ref/`, `references/`, `ref/`, `resources/`, `.`

```python
GENOME_EXTENSIONS = ["*.fna", "*.fa", "*.fasta"]  # also .gz variants
```

For each candidate, check for companion files:
- BWA index: `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` (alongside the FASTA or `.gz`)
- FAI index: `.fai`
- Sequence dict: `.dict`

Present with status indicators so user can see which genomes are ready to use.

### Known-sites VCF scanning

Search directories: `analysis/GATK_resource_bundle/`, `resources/`, `ref/`

```python
KNOWN_SITE_PATTERNS = [
    "*dbsnp*.vcf.gz",
    "*Mills_and_1000G*.vcf.gz",
    "*known_indels*.vcf.gz",
    "*af-only-gnomad*.vcf.gz",
    "*1000g_pon*.vcf.gz",
    "*gnomad*common_biallelic*.vcf.gz",
]
```

Check for `.tbi` index alongside each VCF. Warn if missing.

### FASTQ file scanning and Illumina filename parsing

Search directories: `download/`, `data/`, `fastq/`, `raw/`, or user-specified.

Parse using regex that handles both bcl2fastq2 and BCL Convert formats:

```python
# Handles: SampleA_S1_L001_R1_001.fastq.gz  (bcl2fastq2, with lane)
#          SampleA_S1_R1_001.fastq.gz        (BCL Convert, no lane)
ILLUMINA_RE = re.compile(
    r'^(?P<sample>.+?)_S(?P<snum>\d+)(?:_L(?P<lane>\d{3}))?_(?P<read>R[12])_(?P<set>\d{3})\.fastq\.gz$'
)
```

Fallback regex for non-Illumina naming:

```python
# Handles: anything_R1_001.fastq.gz, anything_R1.fastq.gz
GENERIC_RE = re.compile(
    r'^(?P<basename>.+?)[_.](?P<read>R[12])(?:[_.](?P<set>\d{3}))?\.fastq\.gz$'
)
```

R1/R2 pairing by string substitution (`_R1_` → `_R2_`), then existence check.

Build a pandas DataFrame with columns matching `samples.schema.yaml`:
- `fastq_files_basename`: everything before the R1/R2 suffix
- `lane`: from parsed filename (or `L001` if not present)
- `project_sample`: sample name extracted from filename
- `mdc_project`: prompted from user (project identifier)

---

## Interactive Prompts Details

### Prompt types mapped to config sections

| Config section | Prompt type | Discovery source |
|---|---|---|
| `ref.genome` | `path()` with autocomplete | Pre-filled from discovery |
| `ref.build` | `select()` GRCh37/GRCh38 | Inferred from filename if possible |
| `ref.known_sites` | `checkbox()` multi-select | Pre-filled from discovery |
| `paths.fastq_folder` | `path()` | Pre-filled from discovery |
| `paths.output_folder` | `text()` | Default: `results/` |
| `read_group.platform` | `select()` | Default: ILLUMINA |
| `trimming.enabled` | `confirm()` | — |
| `trimming.*` (if enabled) | `text()` with defaults | From schema defaults |
| `subset.bed_file` | `path()` (if enabled) | — |
| `mdc_project` | `text()` | Not auto-discoverable |

### Conditional flow

- If trimming is disabled → skip all `trimming.*` prompts
- If no BED file → skip `subset.*` prompts
- If only one reference genome found → auto-select, just confirm
- If known-sites found → pre-check all, let user deselect

---

## Config YAML Generation

Use `yaml.safe_dump()` with `sort_keys=False` to preserve the logical section ordering from the schema. Only write non-default values for optional sections (trimming, subset) to keep the output clean.

Add a header comment:

```yaml
# =============================================================================
# config/config.yaml — Generated by setup_config.py
# =============================================================================
```

### Validation

After writing both files, validate against the existing schemas:

```python
from snakemake.utils import validate
import yaml, pandas as pd

config = yaml.safe_load(open("config/config.yaml"))
validate(config, "workflow/schemas/config.schema.yaml")

samples = pd.read_table("config/samples.tsv")
validate(samples, "workflow/schemas/samples.schema.yaml")
```

---

## CLI Arguments

```
usage: setup_config.py [-h] [--config-out PATH] [--samples-out PATH]
                       [--scan-dir DIR] [--no-interactive]

Options:
  --config-out PATH     Output config file (default: config/config.yaml)
  --samples-out PATH    Output samples file (default: config/samples.tsv)
  --scan-dir DIR        Base directory to scan (default: current directory)
  --no-interactive      Auto-discovery only, no prompts (use all defaults)
```

The `--no-interactive` flag enables CI/scripted usage: discover files, apply defaults, write config, validate.

---

## Edge Cases to Handle

| Scenario | Behavior |
|---|---|
| No reference genome found | Prompt for manual path entry |
| No FASTQ files found | Prompt for directory, re-scan |
| No known-sites VCFs found | Warn, prompt for manual paths or skip BQSR |
| Mixed FASTQ naming conventions | Try Illumina regex first, fall back to generic, warn |
| Underscores in sample names (Illumina converts to dashes) | Document this; parse as-is |
| Existing config.yaml | Ask to overwrite or merge |
| Existing samples.tsv | Ask to overwrite or append |
| R2 file missing for an R1 | Warn per file, exclude from samples.tsv |
| `.gz` vs uncompressed reference | Discover both, let user choose |

---

## Implementation Phases

### Phase 1: Core discovery + YAML generation (MVP)

- `discover_reference_genomes()`, `discover_known_sites()`, `discover_fastq_files()`
- `parse_illumina_filename()` with fallback regex
- `generate_config_yaml()` and `generate_samples_tsv()`
- Basic `argparse` CLI
- Print results to stdout (no rich dependency yet)

### Phase 2: Interactive prompts

- Add `questionary` prompts for each config section
- Conditional flow (trimming, subset)
- Path autocompletion for file prompts

### Phase 3: Rich display + polish

- Add `rich` tables for discovery results
- Progress indicators during scanning
- Colored status (green = ready, red = missing index)
- `--no-interactive` mode

### Phase 4: Advanced features

- Detect genome build from FASTA header or filename
- Parse Illumina SampleSheet.csv if found alongside FASTQs
- Support re-running to update an existing config (merge mode)
- Detect if BWA index needs rebuilding (mtime comparison)

---

## Sources

- [nf-core launch](https://nf-co.re/docs/nf-core-tools/pipelines/launch) — schema-driven interactive pipeline config wizard
- [questionary](https://questionary.readthedocs.io/en/stable/) — Python interactive prompt library (used by nf-core)
- [rich](https://rich.readthedocs.io/) — terminal formatting (tables, colors)
- [Illumina FASTQ naming convention](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)
- [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
- [Snakemake config validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)
