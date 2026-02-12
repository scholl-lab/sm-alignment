rule base_recalibration:
    """Generate BQSR recalibration table from known variant sites."""
    input:
        dedup_bam=os.path.join(DEDUP_DIR, "{sample}" + DEDUP_SUFFIX),
    output:
        recal_table=os.path.join(BQSR_DIR, "{sample}" + RECAL_TABLE_SUFFIX),
    params:
        java_opts=get_java_opts,
        reference=REF,
        known_sites=lambda wc: " ".join(f'--known-sites "{ks}"' for ks in KNOWN_SITES),
        extra=config.get("params", {}).get("gatk", {}).get("BaseRecalibrator", ""),
    threads: 4
    resources:
        runtime=4320,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "recal.gatk.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' BaseRecalibrator \
            {params.extra} \
            -I "{input.dedup_bam}" \
            -R "{params.reference}" \
            {params.known_sites} \
            -O "{output.recal_table}" \
            2> {log}
        """


rule apply_bqsr:
    """Apply base quality score recalibration to produce final BAMs."""
    input:
        dedup_bam=os.path.join(DEDUP_DIR, "{sample}" + DEDUP_SUFFIX),
        recal_table=os.path.join(BQSR_DIR, "{sample}" + RECAL_TABLE_SUFFIX),
    output:
        bqsr_bam=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
    params:
        java_opts=lambda wc, resources: (
            f"-Xms{int(resources.mem_mb*0.2)}m"
            f" -Xmx{int(resources.mem_mb*0.8)}m"
            f" -Djava.io.tmpdir={resources.tmpdir}"
            f" -Dsamjdk.compression_level={COMPRESSION_LEVEL}"
        ),
        reference=REF,
        extra=config.get("params", {}).get("gatk", {}).get("ApplyBQSR", ""),
    threads: 4
    resources:
        runtime=4320,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "apply_bqsr.gatk.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' ApplyBQSR \
            {params.extra} \
            -R "{params.reference}" \
            -I "{input.dedup_bam}" \
            -bqsr "{input.recal_table}" \
            -O "{output.bqsr_bam}" \
            2> {log}
        """
