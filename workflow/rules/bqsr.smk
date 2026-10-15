rule base_recalibration:
    """Generate BQSR recalibration table from known variant sites."""
    input:
        dedup_bam=os.path.join(DEDUP_DIR, "{sample}" + DEDUP_SUFFIX),
        dedup_bai=os.path.join(
            DEDUP_DIR, "{sample}" + DEDUP_SUFFIX.removesuffix(".bam") + ".bai"
        ),
        reference=REF,
        known_sites=KNOWN_SITES,
        intervals=[TARGET_INTERVALS] if TARGET_INTERVALS else [],
    output:
        recal_table=os.path.join(BQSR_DIR, "{sample}" + RECAL_TABLE_SUFFIX),
    params:
        java_opts=get_java_opts,
        known_sites=lambda wc, input: " ".join(
            f'--known-sites "{ks}"' for ks in input.known_sites
        ),
        padding=BQSR_INTERVAL_PADDING,
        intervals=lambda wc, input: (
            f'-L "{input.intervals}" --interval-padding {BQSR_INTERVAL_PADDING}'
            if input.intervals
            else ""
        ),
        extra=config.get("params", {}).get("gatk", {}).get("BaseRecalibrator", ""),
    threads: 4
    resources:
        runtime=480,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "recal.gatk.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' BaseRecalibrator \
            {params.extra} \
            -I "{input.dedup_bam}" \
            -R "{input.reference}" \
            {params.known_sites} \
            {params.intervals} \
            -O "{output.recal_table}" \
            2> {log}
        """


rule apply_bqsr:
    """Apply base quality score recalibration to produce final BAMs."""
    input:
        dedup_bam=os.path.join(DEDUP_DIR, "{sample}" + DEDUP_SUFFIX),
        dedup_bai=os.path.join(
            DEDUP_DIR, "{sample}" + DEDUP_SUFFIX.removesuffix(".bam") + ".bai"
        ),
        recal_table=os.path.join(BQSR_DIR, "{sample}" + RECAL_TABLE_SUFFIX),
        reference=REF,
    output:
        bqsr_bam=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
        bqsr_bai=os.path.join(
            BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX.removesuffix(".bam") + ".bai"
        ),
    params:
        java_opts=lambda wc, resources: (
            f"-Xms{int(resources.mem_mb * 0.2)}m"
            f" -Xmx{int(resources.mem_mb * 0.8)}m"
            f" -Djava.io.tmpdir={resources.tmpdir}"
            f" -Dsamjdk.compression_level={COMPRESSION_LEVEL}"
        ),
        extra=config.get("params", {}).get("gatk", {}).get("ApplyBQSR", ""),
    threads: 4
    resources:
        runtime=480,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "apply_bqsr.gatk.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' ApplyBQSR \
            {params.extra} \
            -R "{input.reference}" \
            -I "{input.dedup_bam}" \
            -bqsr "{input.recal_table}" \
            -O "{output.bqsr_bam}" \
            2> {log}
        """
