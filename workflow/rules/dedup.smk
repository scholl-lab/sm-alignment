rule deduplicate_bam_files:
    """Mark PCR duplicates with GATK MarkDuplicates."""
    input:
        merged_bam=os.path.join(MERGED_DIR, "{sample}" + MERGED_SUFFIX),
    output:
        dedup_bam=temp(os.path.join(DEDUP_DIR, "{sample}" + DEDUP_SUFFIX)),
        dedup_bai=temp(
            os.path.join(DEDUP_DIR, "{sample}" + DEDUP_SUFFIX.removesuffix(".bam") + ".bai")
        ),
        metrics=os.path.join(DEDUP_DIR, "{sample}" + DEDUP_METRICS_SUFFIX),
    params:
        java_opts=get_java_opts,
        extra=config.get("params", {})
        .get("gatk", {})
        .get(
            "MarkDuplicates",
            "--CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
            " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500",
        ),
    threads: 4
    resources:
        runtime=480,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "dedup.gatk.{sample}.log"),
    shell:
        r"""
        gatk --java-options '{params.java_opts}' MarkDuplicates \
            {params.extra} \
            -I "{input.merged_bam}" \
            -O "{output.dedup_bam}" \
            -M "{output.metrics}" \
            2> {log}
        """
