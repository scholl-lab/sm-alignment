rule bwa_map:
    """Align paired-end FASTQ per basename, pipe to samtools sort."""
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2,
    output:
        bam=temp(os.path.join(ALIGNED_DIR, "{basename}.bam")),
    params:
        reference=REF_GZ,
        read_group=lambda wc: (
            '"@RG\\tID:{lane}-{sample}\\tSM:{sample}\\tLB:{sample}'
            '\\tPL:{platform}\\tPU:{lane}-{project}"'.format(
                lane=samples_df.loc[wc.basename, "lane"],
                sample=samples_df.loc[wc.basename, "project_sample"],
                project=samples_df.loc[wc.basename, "mdc_project"],
                platform=PLATFORM,
            )
        ),
        sort_threads=2,
        sort_mem=4000,
        extra=config.get("params", {}).get("bwa_mem", {}).get("extra", ""),
    threads: 16
    resources:
        runtime=1440,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        bwa=os.path.join(LOG_DIR, "map.bwa.{basename}.log"),
        samtools=os.path.join(LOG_DIR, "map.samtools.{basename}.log"),
    shell:
        r"""
        BWA_THREADS=$(({threads} - {params.sort_threads}))

        TMP_SORT_DIR=$(mktemp -p {resources.tmpdir} -d samtools-sort.XXXXXX)

        bwa mem -t $BWA_THREADS \
            {params.extra} \
            -R {params.read_group} \
            {params.reference} \
            {input.r1} {input.r2} \
            2> {log.bwa} \
        | samtools sort -@ {params.sort_threads} \
            -m {params.sort_mem}M \
            -O BAM \
            -T "$TMP_SORT_DIR"/tmp \
            -o {output.bam} \
            2> {log.samtools}

        rm -rf "$TMP_SORT_DIR"
        """
