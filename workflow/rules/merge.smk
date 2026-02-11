rule merge_bam_files:
    """Merge all lane-level BAMs for a sample into a single BAM."""
    input:
        lambda wc: [
            os.path.join(ALIGNED_DIR, f"{basename}.bam")
            for basename in get_basenames_for_sample(wc.sample)
        ],
    output:
        merged_bam=temp(os.path.join(MERGED_DIR, "{sample}" + MERGED_SUFFIX)),
    params:
        list_file=os.path.join(MERGED_DIR, "{sample}.bamlist"),
        extra=config.get("params", {}).get("samtools", {}).get("merge_extra", ""),
    threads: 8
    resources:
        runtime=1440,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        os.path.join(LOG_DIR, "merge.samtools.{sample}.log"),
    shell:
        r"""
        echo "{input}" | tr " " "\n" > "{params.list_file}"

        samtools merge -@ {threads} -O BAM \
            {params.extra} \
            -b "{params.list_file}" \
            "{output.merged_bam}" \
            2> {log}

        rm -f "{params.list_file}"
        """
