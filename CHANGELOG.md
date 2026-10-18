# Changelog

All notable changes to `sm-alignment` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- Optional UMI / molecular barcode extraction and UMI-aware duplicate marking (`umi` config block, off by default; #47):
  - Verified kit preset table in `config/umi_presets.yaml` including Agilent SureSelect XT HS2 (`agilent_xths2`: `3M2S+T 3M2S+T`, Agilent protocol G9985-90010, pp. 128, 134).
  - `umi_extract` rule (`fgbio FastqToBam` piped to `samtools fastq -T RX`) upstream of `trim_adapters` / `bwa_map`.
  - Automatic `-C` flag on `bwa mem` when `umi.enabled: true` to transfer `RX:Z:` FASTQ comments into aligned BAM records.
  - UMI-aware GATK `MarkDuplicates` flags (`--BARCODE_TAG`, `--DUPLEX_UMI`, and passthrough `umi.dedup.extra`).
  - JSON Schema and helper validation for `umi.preset`, fgbio read-structure syntax, and `preset` XOR `read_structures` conflict checking.
  - Run provenance logging and `{output_folder}/logs/umi_provenance.yaml` emission.
  - Pinned `fgbio=3.1.2` in `workflow/envs/bwa_samtools.yaml`.

## [0.6.0] - 2026-10-16

### Fixed
- Run Snakemake dry-run and schema validation tests in CI and make `Snakefile` and SLURM launcher stage-safe (#40, #46).

## [0.5.0] - 2026-10-16

### Fixed
- Use explicit time units for `runtime` in `profiles/default`, `profiles/bih`, and `profiles/charite` (#45).
- Pass explicit `-Xmx` derived from job `mem_mb` to `bbduk.sh` in `trim_adapters` (#43, #44).

## [0.4.0] - 2026-10-16

### Added
- Per-kit GATK `CollectHsMetrics` QC rule and `Picard ValidateSamFile` (#38, #41).
- BIH SLURM profile and capped rule runtimes (#39).

### Fixed
- Default `bwa mem` flags (`-K 100000000 -Y`), flowcell/library `@RG` formatting, padded BQSR target intervals, and `MarkDuplicates` optical pixel distance (#34, #35, #36, #37, #38).
