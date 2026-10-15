TRIM_CFG = config.get("trimming", {})


def _resolve_raw_fastq(wildcards, suffix):
    """Locate raw FASTQ for trimming, respecting optional subfolder."""
    return _resolve_fastq_path_impl(
        wildcards.basename, suffix, samples_df, False, TRIMMED_DIR, FASTQ_DIR
    )


rule trim_adapters:
    """Adapter and quality trimming with BBDuk."""
    input:
        r1=lambda wc: _resolve_raw_fastq(wc, R1_SUFFIX),
        r2=lambda wc: _resolve_raw_fastq(wc, R2_SUFFIX),
    output:
        trimmed_r1=os.path.join(TRIMMED_DIR, "{basename}" + TRIMMED_R1_SUFFIX),
        trimmed_r2=os.path.join(TRIMMED_DIR, "{basename}" + TRIMMED_R2_SUFFIX),
    params:
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
        quantize=(
            f"quantize={TRIM_CFG['quantize']}" if TRIM_CFG.get("quantize", "") else ""
        ),
    threads: 4
    resources:
        runtime=720,
    conda:
        "../envs/bbtools.yaml"
    log:
        os.path.join(LOG_DIR, "bbduk_trim.{basename}.log"),
    shell:
        r"""
        bbduk.sh \
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
