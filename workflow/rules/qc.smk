# =============================================================================
# workflow/rules/qc.smk — Quality control rules
# =============================================================================
# FastQC (raw + trimmed FASTQs), samtools stats/flagstat, Picard
# CollectMultipleMetrics, Qualimap bamqc, and MultiQC aggregation.
#
# Controlled by the qc section in config.yaml:
#   qc.enabled             — master switch (default: true)
#   qc.fastqc              — FastQC on FASTQs (default: true)
#   qc.samtools_stats      — samtools stats on final BAM (default: true)
#   qc.samtools_flagstat   — samtools flagstat on final BAM (default: true)
#   qc.picard_collect_metrics — Picard CollectMultipleMetrics (default: true)
#   qc.qualimap            — Qualimap bamqc, opt-in (default: false)
# =============================================================================

FASTQC_RAW_DIR = os.path.join(QC_DIR, "fastqc", "raw")
FASTQC_TRIMMED_DIR = os.path.join(QC_DIR, "fastqc", "trimmed")
SAMTOOLS_QC_DIR = os.path.join(QC_DIR, "samtools")
PICARD_QC_DIR = os.path.join(QC_DIR, "picard")
QUALIMAP_QC_DIR = os.path.join(QC_DIR, "qualimap")

# Picard CollectMultipleMetrics output prefix (contains {sample} wildcard)
_PICARD_PREFIX = os.path.join(PICARD_QC_DIR, "{sample}.multiple_metrics")


# -----------------------------------------------------------------------------
# FastQC helpers — build a lookup from output stem to source FASTQ path.
# FastQC names outputs by stripping compression/FASTQ extensions from the
# input filename, so we precompute exact output stems from the configured
# suffixes.  Handle .fastq.gz, .fq.gz, .fastq, .fq (and bare .gz fallback).
# -----------------------------------------------------------------------------
def _fastqc_stem(suffix):
    """Return the filename stem FastQC will use for a given FASTQ suffix."""
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if suffix.endswith(ext):
            return suffix[: -len(ext)]
    if suffix.endswith(".gz"):
        return suffix[: -len(".gz")]
    return suffix


_RAW_R1_STEM = _fastqc_stem(R1_SUFFIX)  # e.g. "_R1_001"
_RAW_R2_STEM = _fastqc_stem(R2_SUFFIX)  # e.g. "_R2_001"
_TRIM_R1_STEM = _fastqc_stem(TRIMMED_R1_SUFFIX)  # e.g. ".bbduk_R1_001"
_TRIM_R2_STEM = _fastqc_stem(TRIMMED_R2_SUFFIX)  # e.g. ".bbduk_R2_001"

# Maps: fq_stem -> source FASTQ path  (built at parse time)
_FASTQC_RAW_LOOKUP: dict[str, str] = {}
_FASTQC_TRIMMED_LOOKUP: dict[str, str] = {}

for _bn in get_all_basenames():
    _row = samples_df.loc[_bn]
    _sub = _row.get("subfolder", "") if "subfolder" in samples_df.columns else ""
    _raw_base = os.path.join(FASTQ_DIR, str(_sub)) if _sub else FASTQ_DIR

    _FASTQC_RAW_LOOKUP[_bn + _RAW_R1_STEM] = os.path.join(_raw_base, _bn + R1_SUFFIX)
    _FASTQC_RAW_LOOKUP[_bn + _RAW_R2_STEM] = os.path.join(_raw_base, _bn + R2_SUFFIX)

    if TRIMMING_ENABLED:
        _FASTQC_TRIMMED_LOOKUP[_bn + _TRIM_R1_STEM] = os.path.join(
            TRIMMED_DIR, _bn + TRIMMED_R1_SUFFIX
        )
        _FASTQC_TRIMMED_LOOKUP[_bn + _TRIM_R2_STEM] = os.path.join(
            TRIMMED_DIR, _bn + TRIMMED_R2_SUFFIX
        )


# =============================================================================
# FastQC — raw FASTQs
# =============================================================================
rule fastqc_raw:
    """Run FastQC on a raw FASTQ file."""
    input:
        lambda wc: _FASTQC_RAW_LOOKUP[wc.fq_stem],
    output:
        html=os.path.join(FASTQC_RAW_DIR, "{fq_stem}_fastqc.html"),
        zip=os.path.join(FASTQC_RAW_DIR, "{fq_stem}_fastqc.zip"),
    params:
        outdir=FASTQC_RAW_DIR,
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    log:
        os.path.join(LOG_DIR, "fastqc.raw.{fq_stem}.log"),
    shell:
        r"""
        fastqc --quiet --threads {threads} \
            --outdir "{params.outdir}" \
            "{input}" \
            2> {log}
        """


# =============================================================================
# FastQC — trimmed FASTQs (only when trimming is enabled)
# =============================================================================
def _get_trimmed_fastq(wildcards):
    """Resolve trimmed FASTQ path; raise a clear error when trimming is off."""
    if not TRIMMING_ENABLED:
        raise ValueError(
            f"fastqc_trimmed requested for '{wildcards.fq_stem}' but "
            "trimming.enabled is false in config.yaml"
        )
    return _FASTQC_TRIMMED_LOOKUP[wildcards.fq_stem]


rule fastqc_trimmed:
    """Run FastQC on a trimmed FASTQ file."""
    input:
        _get_trimmed_fastq,
    output:
        html=os.path.join(FASTQC_TRIMMED_DIR, "{fq_stem}_fastqc.html"),
        zip=os.path.join(FASTQC_TRIMMED_DIR, "{fq_stem}_fastqc.zip"),
    params:
        outdir=FASTQC_TRIMMED_DIR,
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    log:
        os.path.join(LOG_DIR, "fastqc.trimmed.{fq_stem}.log"),
    shell:
        r"""
        fastqc --quiet --threads {threads} \
            --outdir "{params.outdir}" \
            "{input}" \
            2> {log}
        """


# =============================================================================
# samtools stats
# =============================================================================
rule samtools_stats:
    """Compute alignment statistics on final BAM."""
    input:
        bam=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
    output:
        os.path.join(SAMTOOLS_QC_DIR, "{sample}.stats.txt"),
    threads: 1
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        os.path.join(LOG_DIR, "samtools_stats.{sample}.log"),
    shell:
        r"""
        samtools stats "{input.bam}" > {output} 2> {log}
        """


# =============================================================================
# samtools flagstat
# =============================================================================
rule samtools_flagstat:
    """Compute flag-based read counts on final BAM."""
    input:
        bam=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
    output:
        os.path.join(SAMTOOLS_QC_DIR, "{sample}.flagstat.txt"),
    threads: 1
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        os.path.join(LOG_DIR, "samtools_flagstat.{sample}.log"),
    shell:
        r"""
        samtools flagstat "{input.bam}" > {output} 2> {log}
        """


# =============================================================================
# Picard CollectMultipleMetrics (via GATK4)
# =============================================================================
rule picard_collect_multiple_metrics:
    """Collect comprehensive alignment QC metrics in a single pass."""
    input:
        bam=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
    output:
        multiext(
            _PICARD_PREFIX,
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".quality_distribution_metrics",
            ".quality_distribution.pdf",
            ".quality_by_cycle_metrics",
            ".quality_by_cycle.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
            ".bait_bias_detail_metrics",
            ".bait_bias_summary_metrics",
            ".error_summary_metrics",
            ".pre_adapter_detail_metrics",
            ".pre_adapter_summary_metrics",
            ".quality_yield_metrics",
        ),
    params:
        java_opts=get_java_opts,
        reference=REF,
        output_prefix=_PICARD_PREFIX,
        extra=config.get("params", {}).get("gatk", {}).get("CollectMultipleMetrics", ""),
    threads: 1
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "picard_collect_metrics.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' CollectMultipleMetrics \
            {params.extra} \
            -I "{input.bam}" \
            -R "{params.reference}" \
            -O "{params.output_prefix}" \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM QualityScoreDistribution \
            --PROGRAM MeanQualityByCycle \
            --PROGRAM CollectBaseDistributionByCycle \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectSequencingArtifactMetrics \
            --PROGRAM CollectQualityYieldMetrics \
            2> {log}
        """


# =============================================================================
# Qualimap bamqc (opt-in)
# =============================================================================
rule qualimap_bamqc:
    """Run Qualimap bamqc for coverage and alignment quality analysis."""
    input:
        bam=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
    output:
        directory(os.path.join(QUALIMAP_QC_DIR, "{sample}")),
    params:
        feature_file=(
            f'--feature-file "{QC_CFG.get("qualimap_feature_file")}"'
            if QC_CFG.get("qualimap_feature_file")
            else ""
        ),
    threads: 4
    conda:
        "../envs/qualimap.yaml"
    log:
        os.path.join(LOG_DIR, "qualimap.{sample}.log"),
    shell:
        r"""
        qualimap bamqc \
            -bam "{input.bam}" \
            -outdir "{output}" \
            --java-mem-size={resources.mem_mb}M \
            -nt {threads} \
            {params.feature_file} \
            2> {log}
        """


# =============================================================================
# MultiQC — aggregate all QC outputs into a single report
# =============================================================================
def _collect_multiqc_inputs(wildcards):
    """Build the list of QC outputs to aggregate based on config flags."""
    inputs = []
    all_basenames = get_all_basenames()
    all_samples = get_samples()

    if QC_CFG.get("fastqc", True):
        for bn in all_basenames:
            for stem_suffix in [_RAW_R1_STEM, _RAW_R2_STEM]:
                inputs.append(os.path.join(FASTQC_RAW_DIR, f"{bn}{stem_suffix}_fastqc.zip"))
            if TRIMMING_ENABLED:
                for stem_suffix in [_TRIM_R1_STEM, _TRIM_R2_STEM]:
                    inputs.append(
                        os.path.join(
                            FASTQC_TRIMMED_DIR,
                            f"{bn}{stem_suffix}_fastqc.zip",
                        )
                    )

    if QC_CFG.get("samtools_stats", True):
        inputs.extend(
            expand(
                os.path.join(SAMTOOLS_QC_DIR, "{sample}.stats.txt"),
                sample=all_samples,
            )
        )

    if QC_CFG.get("samtools_flagstat", True):
        inputs.extend(
            expand(
                os.path.join(SAMTOOLS_QC_DIR, "{sample}.flagstat.txt"),
                sample=all_samples,
            )
        )

    if QC_CFG.get("picard_collect_metrics", True):
        inputs.extend(
            expand(
                os.path.join(
                    PICARD_QC_DIR,
                    "{sample}.multiple_metrics.alignment_summary_metrics",
                ),
                sample=all_samples,
            )
        )

    # Always include MarkDuplicates metrics (produced by dedup rule)
    inputs.extend(
        expand(
            os.path.join(DEDUP_DIR, "{sample}" + DEDUP_METRICS_SUFFIX),
            sample=all_samples,
        )
    )

    if QC_CFG.get("qualimap", False):
        inputs.extend(
            expand(
                os.path.join(QUALIMAP_QC_DIR, "{sample}"),
                sample=all_samples,
            )
        )

    return inputs


rule multiqc:
    """Aggregate all QC outputs into a single MultiQC report."""
    input:
        _collect_multiqc_inputs,
    output:
        report=os.path.join(QC_DIR, "multiqc_report.html"),
        data=directory(os.path.join(QC_DIR, "multiqc_data")),
    params:
        outdir=QC_DIR,
        search_dirs=lambda wc: f'"{QC_DIR}" "{DEDUP_DIR}"',
        extra=config.get("params", {}).get("multiqc", {}).get("extra", ""),
    threads: 1
    conda:
        "../envs/multiqc.yaml"
    log:
        os.path.join(LOG_DIR, "multiqc.log"),
    shell:
        r"""
        multiqc \
            --force \
            --outdir "{params.outdir}" \
            --filename multiqc_report.html \
            {params.extra} \
            {params.search_dirs} \
            2> {log}
        """
