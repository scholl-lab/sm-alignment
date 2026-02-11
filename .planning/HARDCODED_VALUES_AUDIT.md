# Hardcoded Values Audit -- Legacy Scripts

This audit covers all hardcoded (non-configurable) values found in the **legacy** scripts
under `scripts/*.sh` and `scripts/snakemake/*.smk`. It does NOT cover the new `workflow/`
code. The purpose is to document what the refactoring plan needs to address and to serve
as a traceability record.

**Scope:**
- 12 shell scripts in `scripts/`
- 12 Snakemake files in `scripts/snakemake/`
- 3 legacy config files (`config.yaml`, `configs/config_alignment.yaml`, `configs/config_trim_adapters.yaml`)

---

## 1. Critical Issues (paths bypassing config, absolute paths)

These are the highest-priority findings: values embedded directly in shell blocks or
Python code that should have been read from configuration but were not.

| # | File | Value | Context | Configurable? |
|---|------|-------|---------|---------------|
| 1.1 | `snakemake/bqsr_bams.smk` (line 63-67) | `analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` | Reference genome path hardcoded directly in `BaseRecalibrator` shell block instead of using config variable | NO -- bypasses `config["reference"]` entirely |
| 1.2 | `snakemake/bqsr_bams.smk` (line 65) | `analysis/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf` | Known-sites VCF hardcoded in shell block | NO -- not in any config file |
| 1.3 | `snakemake/bqsr_bams.smk` (line 66) | `analysis/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz` | Known-sites VCF hardcoded in shell block | NO -- not in any config file |
| 1.4 | `snakemake/bqsr_bams.smk` (line 67) | `analysis/GATK_resource_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | Known-sites VCF hardcoded in shell block | NO -- not in any config file |
| 1.5 | `snakemake/subset_bam.smk` (line 35) | `/data/cephfs-1/work/groups/scholl/shared/target_files/apa_genes/hg38/apa_genes.genes2bed.GRCh38.S33266436_Regions.padding1000bp.bed` | Absolute path to BED file hardcoded as rule input | NO -- not configurable at all |
| 1.6 | `snakemake/subset_bam.smk` (line 7) | `results/exomes` | Base path hardcoded in `prefix_results` instead of coming from config | NO -- no configfile directive |
| 1.7 | `snakemake/bam_to_fastq.smk` (line 7) | `results/exomes` | Base path hardcoded in `prefix_results` instead of coming from config | NO -- no configfile directive |
| 1.8 | `snakemake/align_and_sort.smk` (line 11) | `results/exomes` | Base path hardcoded in `prefix_results` instead of coming from config | NO -- config loaded but path is still hardcoded |
| 1.9 | `snakemake/merge_bams.smk` (line 61) | `results/merged` | Hardcoded path in `expand()` target and `mkdir -p results/merged` shell command | Partially -- `MERGE_DIR` exists as variable but `expand()` and `mkdir` use literal |
| 1.10 | `snakemake/dedup_bams.smk` (line 44) | `results/dedup` | Hardcoded in `expand()` target and `mkdir -p results/dedup` shell command | Partially -- `DEDUP_DIR` variable exists but `expand()` uses literal string |
| 1.11 | `snakemake/bqsr_bams.smk` (line 39) | `results/bqsr` | Hardcoded in `expand()` target | Partially -- `BQSR_DIR` variable exists but `expand()` uses literal string |
| 1.12 | `snakemake/alignment_v2.smk` (line 27) | `metadata.tsv` | Metadata file path hardcoded as literal `'metadata.tsv'` | NO -- should be in config |
| 1.13 | `snakemake/merge_bams_v2.smk` (line 26) | `metadata.tsv` | Metadata file path hardcoded as literal `'metadata.tsv'` | NO -- should be in config |
| 1.14 | `snakemake/alignment_v2.smk` (line 16) | `config.yaml` | Config file re-opened with `open("config.yaml")` instead of using Snakemake's `config` dict | NO -- fragile; breaks if `--configfile` points elsewhere |
| 1.15 | `snakemake/merge_bams_v2.smk` (line 13) | `config.yaml` | Same re-open issue as above | NO |
| 1.16 | `snakemake/merge_bams.smk` (line 18) | `results/aligned/` | Fallback default for `ALIGNED_DIR` hardcoded as literal | Partially -- `config.get("aligned_folder", "results/aligned/")` with hardcoded default |

**Note:** Items 1.1-1.4 are the most dangerous finding. The `bqsr_bams.smk` standalone
workflow hardcodes different known-sites files (dbSNP, known indels, Mills) than those
referenced in `config.yaml` and `configs/config_alignment.yaml` (gnomAD, 1000G PoN).
This means the standalone BQSR workflow uses an entirely different set of calibration
resources than the integrated pipeline, with no configuration control.

---

## 2. Conda Environment Names

Named conda environments (e.g., `conda: "gatk"`) are not portable. Anyone cloning the
repo must manually create identically-named local environments.

| # | File | Value | Rule(s) | Configurable? |
|---|------|-------|---------|---------------|
| 2.1 | `snakemake/alignment_pipeline.smk` | `conda: "base"` | `bwa_map`, `merge_bam_files` | NO |
| 2.2 | `snakemake/alignment_pipeline.smk` | `conda: "gatk"` | `deduplicate_bam_files`, `base_recalibration`, `apply_bqsr` | NO |
| 2.3 | `snakemake/dedup_bams.smk` | `conda: "gatk"` | `deduplicate_bam_files` | NO |
| 2.4 | `snakemake/bqsr_bams.smk` | `conda: "gatk"` | `base_recalibration`, `apply_bqsr` | NO |
| 2.5 | `snakemake/trim_adapters.smk` | `conda: "base"` | `trim_adapters` | NO |
| 2.6 | `snakemake/subset_bam.smk` | `conda: "base"` | `samtools_view_apa_genes` | NO |
| 2.7 | `snakemake/bam_to_fastq.smk` | `conda: "base"` | `convert_bam_to_fastq` | NO |
| 2.8 | `snakemake/align_and_sort.smk` | `conda: "base"` | `align_and_sort`, `index_bam` | NO |

**Total occurrences:** 12 (8 using `"base"`, 4 using `"gatk"`)

---

## 3. SLURM Parameters

All shell submission scripts share the same boilerplate SBATCH headers. These values
cannot be overridden without editing the scripts directly.

| # | File | Parameter | Value | Configurable? |
|---|------|-----------|-------|---------------|
| 3.1 | `run_alignment.sh` | `--job-name` | `sm_alignment_main_job` | NO |
| 3.2 | `run_alignment_v2.sh` | `--job-name` | `sm_alignment_v2_main_job` | NO |
| 3.3 | `run_bqsr_bams.sh` | `--job-name` | `sm_bqsr_main_job` | NO |
| 3.4 | `run_dedup_bams.sh` | `--job-name` | `sm_dedup_main_job` | NO |
| 3.5 | `run_md5sum_files.sh` | `--job-name` | `sm_md5sum_files_main_job` | NO |
| 3.6 | `run_merge_bams.sh` | `--job-name` | `sm_merge_bams_main_job` | NO |
| 3.7 | `submit_bam_to_fastq.sh` | `--job-name` | `sm_bam_to_fastq_main_job` | NO |
| 3.8 | `submit_align_and_sort.sh` | `--job-name` | `sm_align_and_sort_main_job` | NO |
| 3.9 | `submit_alignment_pipeline.sh` | `--job-name` | `sm_alignment_main_job` | NO |
| 3.10 | `submit_trim_adapters.sh` | `--job-name` | `sm_trim_adapters_main_job` | NO |
| 3.11 | `run_subset_bam.sh` | `--job-name` | `sm_subset_bam_main_job` | NO |
| 3.12 | All 11 legacy `.sh` | `--ntasks` | `1` | NO |
| 3.13 | All 11 legacy `.sh` | `--nodes` | `1` | NO |
| 3.14 | 9 of 11 `.sh` | `--time` | `168:00:00` (7 days) | NO |
| 3.15 | `submit_alignment_pipeline.sh` | `--time` | `72:00:00` (3 days) | NO |
| 3.16 | `submit_trim_adapters.sh` | `--time` | `72:00:00` (3 days) | NO |
| 3.17 | All 11 legacy `.sh` | `--mem` | `2000M` | NO |
| 3.18 | All 11 legacy `.sh` | `--output` | `slurm_logs/%x-%j.log` | NO |
| 3.19 | All 11 legacy `.sh` | `SBATCH_DEFAULTS` | `--output=slurm_logs/%x-%j.log` | NO |
| 3.20 | All 11 legacy `.sh` | `TMPDIR` base | `$HOME/scratch/tmp` | NO -- hardcoded path |
| 3.21 | All 11 legacy `.sh` | `--profile` | `cubi-v1` | NO (except `submit_alignment_pipeline.sh`, `submit_trim_adapters.sh`) |
| 3.22 | `run_alignment.sh` | `-j` (max jobs) | `10` | NO |
| 3.23 | `run_alignment_v2.sh` | `-j` (max jobs) | `10` | NO |
| 3.24 | `run_bqsr_bams.sh` | `-j` (max jobs) | `20` | NO |
| 3.25 | `run_dedup_bams.sh` | `-j` (max jobs) | `20` | NO |
| 3.26 | `run_md5sum_files.sh` | `-j` (max jobs) | `20` | NO |
| 3.27 | `run_merge_bams.sh` | `-j` (max jobs) | `20` | NO |
| 3.28 | `submit_bam_to_fastq.sh` | `-j` (max jobs) | `20` | NO |
| 3.29 | `submit_align_and_sort.sh` | `-j` (max jobs) | `100` | NO |
| 3.30 | `run_subset_bam.sh` | `-j` (max jobs) | `20` | NO |
| 3.31 | `submit_alignment_pipeline.sh` | `-j` default | `20` | YES -- `$3` argument with default |
| 3.32 | `submit_trim_adapters.sh` | `-j` default | `20` | YES -- `$3` argument with default |

---

## 4. Thread Counts

Thread values set in Snakemake rule `threads:` directives.

| # | File | Rule | Value | Configurable? |
|---|------|------|-------|---------------|
| 4.1 | `snakemake/alignment.smk` | `bwa_map` | `threads: 16` | NO -- hardcoded in rule |
| 4.2 | `snakemake/alignment_v2.smk` | `bwa_map` | `threads: 16` | NO |
| 4.3 | `snakemake/alignment_pipeline.smk` | `bwa_map` | `threads: 16` | NO |
| 4.4 | `snakemake/alignment_pipeline.smk` | `merge_bam_files` | `threads: 8` | NO |
| 4.5 | `snakemake/alignment_pipeline.smk` | `deduplicate_bam_files` | `threads: 4` | NO |
| 4.6 | `snakemake/alignment_pipeline.smk` | `base_recalibration` | `threads: 4` | NO |
| 4.7 | `snakemake/alignment_pipeline.smk` | `apply_bqsr` | `threads: 4` | NO |
| 4.8 | `snakemake/merge_bams.smk` | `merge_bam_files` | `threads: 8` | NO |
| 4.9 | `snakemake/merge_bams_v2.smk` | `merge_bam_files` | `threads: 8` | NO |
| 4.10 | `snakemake/dedup_bams.smk` | `deduplicate_bam_files` | `threads: 4` | NO |
| 4.11 | `snakemake/bqsr_bams.smk` | `base_recalibration` | `threads: 4` | NO |
| 4.12 | `snakemake/bqsr_bams.smk` | `apply_bqsr` | `threads: 4` | NO |
| 4.13 | `snakemake/trim_adapters.smk` | `trim_adapters` | `threads: 4` | NO (note: the rule also has `bbduk_threads` resource set from config `THREADS=8`) |
| 4.14 | `snakemake/subset_bam.smk` | `samtools_view_apa_genes` | `threads: 1` | NO |
| 4.15 | `snakemake/bam_to_fastq.smk` | `convert_bam_to_fastq` | `threads: 8` | NO |
| 4.16 | `snakemake/align_and_sort.smk` | `align_and_sort` | `threads: 8` | NO |
| 4.17 | `snakemake/align_and_sort.smk` | `index_bam` | `threads: 1` | NO |

Additionally, several helper functions derive sub-thread counts from the total:

| # | File | Function | Formula | Configurable? |
|---|------|----------|---------|---------------|
| 4.18 | `snakemake/alignment.smk` | `get_bwa_threads()` | `threads - 2` | NO -- hardcoded offset |
| 4.19 | `snakemake/alignment.smk` | `get_sort_threads()` | `return 2` | NO -- literal constant |
| 4.20 | `snakemake/alignment_v2.smk` | `get_bwa_threads()` | `threads - 2` | NO |
| 4.21 | `snakemake/alignment_v2.smk` | `get_sort_threads()` | `return 2` | NO |
| 4.22 | `snakemake/alignment_pipeline.smk` | `get_bwa_threads()` | `max(1, threads - 2)` | NO |
| 4.23 | `snakemake/alignment_pipeline.smk` | `get_sort_threads()` | `return 2` | NO |

---

## 5. Memory Values

Memory allocation values, both in helper functions and direct resource declarations.

| # | File | Rule/Function | Value | Configurable? |
|---|------|---------------|-------|---------------|
| 5.1 | `snakemake/alignment.smk` | `get_mem_from_threads()` | `threads * 1200` (MB) | NO -- multiplier hardcoded |
| 5.2 | `snakemake/alignment.smk` | `get_sort_mem()` | `2 * 1200` = 2400 MB | NO |
| 5.3 | `snakemake/alignment_v2.smk` | `get_mem_from_threads()` | `threads * 1200` (MB) | NO |
| 5.4 | `snakemake/alignment_v2.smk` | `get_sort_mem()` | `2 * 1200` = 2400 MB | NO |
| 5.5 | `snakemake/alignment_pipeline.smk` | `get_mem_from_threads()` | `threads * 1200` (MB) | NO |
| 5.6 | `snakemake/alignment_pipeline.smk` | `get_sort_mem()` | `2 * 2000` = 4000 MB | NO -- note: different multiplier (2000) than alignment.smk (1200) |
| 5.7 | `snakemake/alignment_pipeline.smk` | `deduplicate_bam_files` | `threads * 4400` (MB) | NO |
| 5.8 | `snakemake/alignment_pipeline.smk` | `base_recalibration` | `threads * 4400` (MB) | NO |
| 5.9 | `snakemake/alignment_pipeline.smk` | `apply_bqsr` | `threads * 4400` (MB) | NO |
| 5.10 | `snakemake/merge_bams.smk` | `get_mem_from_threads()` | `threads * 1200` (MB) | NO |
| 5.11 | `snakemake/merge_bams_v2.smk` | `get_mem_from_threads()` | `threads * 1200` (MB) | NO |
| 5.12 | `snakemake/dedup_bams.smk` | `get_mem_from_threads()` | `threads * 4400` (MB) | NO |
| 5.13 | `snakemake/bqsr_bams.smk` | `get_mem_from_threads()` | `threads * 4400` (MB) | NO |
| 5.14 | `snakemake/trim_adapters.smk` | `get_mem_from_threads()` | `threads * 1000` (MB) | NO |
| 5.15 | `snakemake/align_and_sort.smk` | `align_and_sort` | `mem_mb = 15000` (literal) | NO |
| 5.16 | `snakemake/align_and_sort.smk` | `index_bam` | `mem_mb = 15000` (literal) | NO |

**Inconsistency:** The memory-per-thread multiplier varies across files:
- Alignment rules: 1200 MB/thread
- Alignment pipeline sort: 2000 MB/thread (different from standalone)
- GATK rules (dedup/BQSR): 4400 MB/thread
- Trim adapters: 1000 MB/thread
- Align-and-sort: flat 15000 MB regardless of threads

---

## 6. Java/GATK Parameters

JVM heap settings and GATK-specific flags embedded in shell blocks.

| # | File | Rule | Value | Configurable? |
|---|------|------|-------|---------------|
| 6.1 | `snakemake/alignment_pipeline.smk` | `deduplicate_bam_files` | `-Xms4000m -Xmx7g` | NO |
| 6.2 | `snakemake/alignment_pipeline.smk` | `deduplicate_bam_files` | `--CREATE_INDEX true` | NO |
| 6.3 | `snakemake/alignment_pipeline.smk` | `deduplicate_bam_files` | `--VALIDATION_STRINGENCY SILENT` | NO |
| 6.4 | `snakemake/alignment_pipeline.smk` | `base_recalibration` | `-Xms4000m -Xmx7g` | NO |
| 6.5 | `snakemake/alignment_pipeline.smk` | `apply_bqsr` | `-Xms4000m -Xmx7g` | NO |
| 6.6 | `snakemake/alignment_pipeline.smk` | `apply_bqsr` | `-Dsamjdk.compression_level=6` | NO |
| 6.7 | `snakemake/dedup_bams.smk` | `deduplicate_bam_files` | `-Xms4000m -Xmx7g` | NO |
| 6.8 | `snakemake/dedup_bams.smk` | `deduplicate_bam_files` | `--CREATE_INDEX` | NO |
| 6.9 | `snakemake/dedup_bams.smk` | `deduplicate_bam_files` | `--VALIDATION_STRINGENCY SILENT` | NO |
| 6.10 | `snakemake/bqsr_bams.smk` | `base_recalibration` | `-Xms4000m -Xmx7g` | NO |
| 6.11 | `snakemake/bqsr_bams.smk` | `apply_bqsr` | `-Xms4000m -Xmx7g` | NO |
| 6.12 | `snakemake/bqsr_bams.smk` | `apply_bqsr` | `-Dsamjdk.compression_level=6` | NO |

**Note:** The JVM heap (`-Xmx7g` = 7168 MB) is disconnected from the SLURM memory
allocation (`threads * 4400` = 17600 MB for 4 threads). The JVM uses only ~40% of
allocated memory, wasting cluster resources.

---

## 7. Rule Time Limits

Time limits set in Snakemake rule `resources: time` directives. All use the legacy
`HH:MM:SS` string format rather than the Snakemake 8+ `runtime` (integer minutes).

| # | File | Rule | Value | Configurable? |
|---|------|------|-------|---------------|
| 7.1 | `snakemake/alignment.smk` | `bwa_map` | `time = '24:00:00'` | NO |
| 7.2 | `snakemake/alignment_v2.smk` | `bwa_map` | `time = '24:00:00'` | NO |
| 7.3 | `snakemake/alignment_pipeline.smk` | `bwa_map` | `time = "24:00:00"` | NO |
| 7.4 | `snakemake/alignment_pipeline.smk` | `merge_bam_files` | `time = "24:00:00"` | NO |
| 7.5 | `snakemake/alignment_pipeline.smk` | `deduplicate_bam_files` | `time = "72:00:00"` | NO |
| 7.6 | `snakemake/alignment_pipeline.smk` | `base_recalibration` | `time = "72:00:00"` | NO |
| 7.7 | `snakemake/alignment_pipeline.smk` | `apply_bqsr` | `time = "72:00:00"` | NO |
| 7.8 | `snakemake/merge_bams.smk` | `merge_bam_files` | `time = '24:00:00'` | NO |
| 7.9 | `snakemake/merge_bams_v2.smk` | `merge_bam_files` | `time = '24:00:00'` | NO |
| 7.10 | `snakemake/dedup_bams.smk` | `deduplicate_bam_files` | `time = '72:00:00'` | NO |
| 7.11 | `snakemake/bqsr_bams.smk` | `base_recalibration` | `time = '72:00:00'` | NO |
| 7.12 | `snakemake/bqsr_bams.smk` | `apply_bqsr` | `time = '72:00:00'` | NO |
| 7.13 | `snakemake/trim_adapters.smk` | `trim_adapters` | `time = "12:00:00"` | NO |

---

## 8. File Naming Patterns

FASTQ suffixes, BAM suffixes, and output naming conventions embedded in code.

| # | File | Value | Context | Configurable? |
|---|------|-------|---------|---------------|
| 8.1 | `snakemake/alignment.smk` | `_R1_001.fastq.gz` | `get_input_files()` filter and `get_rg()` parser | NO |
| 8.2 | `snakemake/alignment.smk` | `_R1_001.fastq.gz` / `_R2_001.fastq.gz` | Input patterns in `bwa_map` rule | NO |
| 8.3 | `snakemake/alignment_v2.smk` | `_R1_001.fastq.gz` / `_R2_001.fastq.gz` | Input lambdas in `bwa_map` rule | NO |
| 8.4 | `snakemake/alignment_pipeline.smk` | `.bbduk_R1_001.fastq.gz` / `.bbduk_R2_001.fastq.gz` | `find_r1()` and `find_r2()` functions | NO -- trimmed FASTQ suffix hardcoded |
| 8.5 | `snakemake/alignment_pipeline.smk` | `.merged.bam` | Merge output suffix in `get_merged_bam()` | NO |
| 8.6 | `snakemake/alignment_pipeline.smk` | `.merged.dedup.bam` | Dedup output suffix in `get_dedup_bam()` | NO |
| 8.7 | `snakemake/alignment_pipeline.smk` | `.merged.dedup_metrics.txt` | Dedup metrics filename suffix | NO |
| 8.8 | `snakemake/alignment_pipeline.smk` | `.merged.dedup.recal_data.table` | Recal table suffix in `get_recal_table()` | NO |
| 8.9 | `snakemake/alignment_pipeline.smk` | `.merged.dedup.bqsr.bam` | Final BQSR BAM suffix in `get_bqsr_bam()` and `apply_bqsr` output | NO |
| 8.10 | `snakemake/merge_bams.smk` | `.merged.bam` | Merge output suffix in `expand()` and output | NO |
| 8.11 | `snakemake/merge_bams_v2.smk` | `.merged.bam` | Merge output suffix | NO |
| 8.12 | `snakemake/dedup_bams.smk` | `.merged.bam` | Input filename pattern | NO |
| 8.13 | `snakemake/dedup_bams.smk` | `.merged.dedup.bam` | Dedup output suffix | NO |
| 8.14 | `snakemake/dedup_bams.smk` | `.merged.dedup_metrics.txt` | Metrics suffix | NO |
| 8.15 | `snakemake/bqsr_bams.smk` | `.merged.dedup.bam` | Input suffix | NO |
| 8.16 | `snakemake/bqsr_bams.smk` | `.merged.dedup.recal_data.table` | Recal table suffix | NO |
| 8.17 | `snakemake/bqsr_bams.smk` | `.merged.dedup.bqsr.bam` | Final BQSR BAM suffix | NO |
| 8.18 | `snakemake/trim_adapters.smk` | `_R1_001.fastq.gz` / `_R2_001.fastq.gz` | Input FASTQ suffixes in `get_samples()`, `find_r1()`, `find_r2()` | NO |
| 8.19 | `snakemake/trim_adapters.smk` | `.bbduk_R1_001.fastq.gz` / `.bbduk_R2_001.fastq.gz` | Trimmed output suffixes | NO |
| 8.20 | `snakemake/align_and_sort.smk` | `.bbduk_R1.fastq.gz` / `.bbduk_R2.fastq.gz` | Input FASTQ suffix (note: different pattern from trim_adapters output!) | NO |
| 8.21 | `snakemake/align_and_sort.smk` | `.sorted.bam` / `.sorted.bam.bai` | Output BAM suffixes | NO |
| 8.22 | `snakemake/bam_to_fastq.smk` | `_R1.fastq.gz` / `_R2.fastq.gz` | Output FASTQ suffixes | NO |
| 8.23 | `snakemake/subset_bam.smk` | `.apa-genes.bam` | Output BAM suffix | NO |
| 8.24 | `snakemake/md5sum_files.smk` | `.fastq.gz` | Input file extension filter | NO |
| 8.25 | `snakemake/md5sum_files.smk` | `.md5sum` | Output checksum file extension | NO |
| 8.26 | `snakemake/md5sum_files.smk` | `all_md5sums.txt` | Concatenated output filename | NO |

**Inconsistency found:** `align_and_sort.smk` expects input files with suffix
`.bbduk_R1.fastq.gz` (no underscore before `001`), while `trim_adapters.smk` produces
files with suffix `.bbduk_R1_001.fastq.gz`. These two workflows would not connect
without manual renaming.

---

## 9. Read Group and Platform Values

| # | File | Value | Context | Configurable? |
|---|------|-------|---------|---------------|
| 9.1 | `snakemake/alignment.smk` | `ILLUMINA` | Hardcoded `\\tPL:ILLUMINA` in `get_rg()` function | NO |
| 9.2 | `snakemake/alignment.smk` | `ILLUMINA` | Hardcoded `\\tPL:ILLUMINA` in `bwa_map` params | NO |
| 9.3 | `snakemake/alignment_v2.smk` | `ILLUMINA` | Hardcoded `\\tPL:ILLUMINA` in `bwa_map` params lambda | NO |
| 9.4 | `snakemake/alignment_pipeline.smk` | `ILLUMINA` | Hardcoded `\\tPL:ILLUMINA` in `bwa_map` params lambda | NO |
| 9.5 | `snakemake/align_and_sort.smk` | `ILLUMINA` | Hardcoded `\\tPL:ILLUMINA` in read group params lambda | NO |

---

## 10. Configfile Directives

Every standalone `.smk` file that loads config uses a hardcoded path.

| # | File | Value | Configurable? |
|---|------|-------|---------------|
| 10.1 | `snakemake/alignment.smk` | `configfile: "config.yaml"` | NO -- fixed filename |
| 10.2 | `snakemake/alignment_v2.smk` | `configfile: "config.yaml"` | NO |
| 10.3 | `snakemake/alignment_pipeline.smk` | `configfile: "config.yaml"` | Partially -- can override with `--configfile` CLI |
| 10.4 | `snakemake/merge_bams.smk` | `configfile: "config.yaml"` | NO |
| 10.5 | `snakemake/merge_bams_v2.smk` | `configfile: "config.yaml"` | NO |
| 10.6 | `snakemake/dedup_bams.smk` | `configfile: "config.yaml"` | NO |
| 10.7 | `snakemake/bqsr_bams.smk` | `configfile: "config.yaml"` | NO |
| 10.8 | `snakemake/trim_adapters.smk` | `configfile: "config.yaml"` | Partially |
| 10.9 | `snakemake/md5sum_files.smk` | `configfile: "config.yaml"` | NO |
| 10.10 | `snakemake/align_and_sort.smk` | `configfile: "config.yaml"` | NO |
| 10.11 | `snakemake/subset_bam.smk` | (no configfile at all) | N/A -- entirely hardcoded |
| 10.12 | `snakemake/bam_to_fastq.smk` | (no configfile at all) | N/A -- entirely hardcoded |

---

## 11. Other Hardcoded Values

| # | File | Value | Context | Configurable? |
|---|------|-------|---------|---------------|
| 11.1 | `snakemake/alignment_pipeline.smk` | `/tmp` | Fallback default in `os.environ.get("TMPDIR", "/tmp")` | Partially -- env var, but default is `/tmp` |
| 11.2 | `snakemake/trim_adapters.smk` | `/tmp` | Same fallback default | Partially |
| 11.3 | `snakemake/subset_bam.smk` | `samtools_view_apa_genes_padding` | Output subdirectory name hardcoded | NO |
| 11.4 | `snakemake/bam_to_fastq.smk` | `samtools_view_apa_genes_padding` | Input/output subdirectory name hardcoded | NO |
| 11.5 | `snakemake/align_and_sort.smk` | `samtools_view_apa_genes_padding` | Input/output subdirectory name hardcoded | NO |
| 11.6 | All legacy `.sh` | `slurm_logs` | Log output directory name | NO |
| 11.7 | All legacy `.sh` | `--use-conda` | Snakemake conda flag always passed | NO |
| 11.8 | `submit_bam_to_fastq.sh` | `bam_to_fastq.smk` | Snakefile path (no `snakemake/` prefix unlike others) | NO -- likely a bug |
| 11.9 | `submit_align_and_sort.sh` | `align_and_sort.smk` | Snakefile path (no `snakemake/` prefix) | NO -- likely a bug |
| 11.10 | `snakemake/alignment.smk` | `download/` | Input directory structure assumed as `download/{run}/` | Partially -- `config["fastq_folder"]` exists but rule input uses literal `download/` |
| 11.11 | `snakemake/alignment.smk` | `_` (underscore) | Filename parsing assumes underscore-delimited fields at specific positions (split("_")[3], etc.) | NO -- fragile naming convention |
| 11.12 | `snakemake/merge_bams.smk` | regex `(.+)_DNA_(\d+)_(.+)_S\d+_L\d+_lane\d+\.bam` | Hardcoded filename regex for sample detection | NO -- breaks on different naming conventions |
| 11.13 | `configs/config_alignment.yaml` | `metadata/metadata_exomes.tsv` | Metadata file path | YES -- in config, but path is project-specific |

---

## Summary Table

| Category | Count | Configurable | Not Configurable |
|----------|------:|:------------:|:----------------:|
| 1. Critical issues (paths bypassing config) | 16 | 1 | 15 |
| 2. Conda environment names | 12 | 0 | 12 |
| 3. SLURM parameters | 32 | 2 | 30 |
| 4. Thread counts | 23 | 0 | 23 |
| 5. Memory values | 16 | 0 | 16 |
| 6. Java/GATK parameters | 12 | 0 | 12 |
| 7. Rule time limits | 13 | 0 | 13 |
| 8. File naming patterns | 26 | 0 | 26 |
| 9. Read group / platform | 5 | 0 | 5 |
| 10. Configfile directives | 12 | 2 | 10 |
| 11. Other | 13 | 2 | 11 |
| **TOTAL** | **180** | **7** | **173** |

Only **3.9%** of hardcoded values are configurable. The remaining **96.1%** require
editing source files to change.

---

## Key Observations

1. **Dangerous divergence in BQSR known-sites.** The standalone `bqsr_bams.smk` hardcodes
   dbSNP + known indels + Mills, while the integrated `alignment_pipeline.smk` uses
   gnomAD + 1000G PoN from config. These are entirely different calibration resource
   sets with no overlap.

2. **Inconsistent memory formulas.** The same logical step (e.g., samtools sort) uses
   different memory multipliers across files (1200 vs 2000 MB per thread), making
   resource usage unpredictable.

3. **JVM heap is disconnected from SLURM allocation.** GATK rules allocate 17.6 GB via
   SLURM but cap JVM at 7 GB, wasting 60% of reserved memory.

4. **FASTQ suffix mismatch between workflows.** `trim_adapters.smk` produces
   `.bbduk_R1_001.fastq.gz` but `align_and_sort.smk` expects `.bbduk_R1.fastq.gz`.

5. **11 nearly-identical submission scripts** differ only in job name and Snakefile path.
   Two of them (`submit_bam_to_fastq.sh`, `submit_align_and_sort.sh`) are missing the
   `snakemake/` directory prefix, which is likely a bug.

6. **No config validation.** No schema checking means typos in config keys fail silently
   at runtime, often deep into multi-day pipeline runs.

7. **Three workflows have no configfile at all** (`subset_bam.smk`, `bam_to_fastq.smk`,
   `align_and_sort.smk` only partially). Every path in these files is hardcoded.
