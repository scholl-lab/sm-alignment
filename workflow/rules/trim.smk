TRIM_CFG = config.get("trimming", {})


def _resolve_raw_fastq(wildcards, suffix):
    """Locate raw FASTQ for UMI extraction or trimming, respecting optional subfolder."""
    return _resolve_fastq_path_impl(
        wildcards.basename, suffix, samples_df, False, TRIMMED_DIR, FASTQ_DIR
    )


def _resolve_trim_input(wildcards, raw_suffix, umi_suffix):
    """Locate input FASTQ for trimming (UMI-extracted when UMI enabled, else raw)."""
    if UMI_ENABLED:
        return os.path.join(UMI_DIR, f"{wildcards.basename}{umi_suffix}")
    return _resolve_raw_fastq(wildcards, raw_suffix)


rule umi_extract:
    """Extract UMI / molecular barcode and skip dark bases with fgbio FastqToBam + samtools fastq."""
    input:
        r1=lambda wc: _resolve_raw_fastq(wc, R1_SUFFIX),
        r2=lambda wc: _resolve_raw_fastq(wc, R2_SUFFIX),
    output:
        r1=temp(os.path.join(UMI_DIR, "{basename}" + UMI_R1_SUFFIX)),
        r2=temp(os.path.join(UMI_DIR, "{basename}" + UMI_R2_SUFFIX)),
    params:
        xmx=lambda wc, resources: int(getattr(resources, "mem_mb", UMI_FGBIO_MEM_MB) * 0.85),
        read_structures=UMI_READ_STRUCTURES,
        preset=UMI_CFG.get("preset") or ("custom" if UMI_CFG.get("read_structures") else "none"),
        umi_tag=UMI_TAG,
        sample=lambda wc: _get_read_group_fields_impl(samples_df, wc.basename, PLATFORM)["sample"],
        library=lambda wc: _get_read_group_fields_impl(samples_df, wc.basename, PLATFORM)[
            "library"
        ],
        rg_id=lambda wc: _get_read_group_fields_impl(samples_df, wc.basename, PLATFORM)["id"],
        platform=PLATFORM,
        rg_pu=lambda wc: _get_read_group_fields_impl(samples_df, wc.basename, PLATFORM)["pu"],
        extra=UMI_FGBIO_EXTRA,
    threads: 4
    resources:
        mem_mb=UMI_FGBIO_MEM_MB,
        runtime=240,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        os.path.join(LOG_DIR, "umi_extract.{basename}.log"),
    shell:
        r"""
        echo "umi_extract: preset={params.preset} read_structures='{params.read_structures}' umi_tag={params.umi_tag}" > {log}

        fgbio -Xmx{params.xmx}m \
            --compression 0 \
            --tmp-dir "{resources.tmpdir}" \
            FastqToBam \
            --input "{input.r1}" "{input.r2}" \
            --output /dev/stdout \
            --read-structures {params.read_structures} \
            --umi-tag {params.umi_tag} \
            --sample "{params.sample}" \
            --library "{params.library}" \
            --read-group-id "{params.rg_id}" \
            --platform "{params.platform}" \
            --platform-unit "{params.rg_pu}" \
            {params.extra} \
            2>> {log} \
        | samtools fastq \
            -@ {threads} \
            -T {params.umi_tag} \
            -1 "{output.r1}" \
            -2 "{output.r2}" \
            -0 /dev/null \
            -s /dev/null \
            -n \
            - \
            2>> {log}
        """


rule trim_adapters:
    """Adapter and quality trimming with BBDuk."""
    input:
        r1=lambda wc: _resolve_trim_input(wc, R1_SUFFIX, UMI_R1_SUFFIX),
        r2=lambda wc: _resolve_trim_input(wc, R2_SUFFIX, UMI_R2_SUFFIX),
    output:
        trimmed_r1=os.path.join(TRIMMED_DIR, "{basename}" + TRIMMED_R1_SUFFIX),
        trimmed_r2=os.path.join(TRIMMED_DIR, "{basename}" + TRIMMED_R2_SUFFIX),
    params:
        xmx=lambda wc, resources: int(getattr(resources, "mem_mb", 4000) * 0.85),
        ref=TRIM_CFG.get("bbduk_ref", "adapters,artifacts"),
        ziplevel=TRIM_CFG.get("ziplevel", 5),
        ktrim=TRIM_CFG.get("ktrim", "r"),
        k=TRIM_CFG.get("k", 23),
        mink=TRIM_CFG.get("mink", 11),
        hdist=TRIM_CFG.get("hdist", 1),
        tpe=TRIM_CFG.get("tpe", "t"),
        tbo=TRIM_CFG.get("tbo", "t"),
        ftl=TRIM_CFG.get("ftl", 0),
        trimpolyg=TRIM_CFG.get("trimpolyg", 0),
        trimpolya=TRIM_CFG.get("trimpolya", 0),
        qtrim=TRIM_CFG.get("qtrim", "f"),
        trimq=TRIM_CFG.get("trimq", 10),
        quantize=(f"quantize={TRIM_CFG['quantize']}" if TRIM_CFG.get("quantize", "") else ""),
    threads: 4
    resources:
        runtime=240,
    conda:
        "../envs/bbtools.yaml"
    log:
        os.path.join(LOG_DIR, "bbduk_trim.{basename}.log"),
    shell:
        r"""
        bbduk.sh \
            -Xmx{params.xmx}m \
            threads={threads} \
            in={input.r1} in2={input.r2} \
            out={output.trimmed_r1} out2={output.trimmed_r2} \
            ziplevel={params.ziplevel} \
            ref={params.ref} \
            ktrim={params.ktrim} k={params.k} mink={params.mink} hdist={params.hdist} \
            tpe={params.tpe} tbo={params.tbo} \
            ftl={params.ftl} trimpolyg={params.trimpolyg} trimpolya={params.trimpolya} \
            qtrim={params.qtrim} trimq={params.trimq} {params.quantize} \
            > {log} 2>&1
        """
