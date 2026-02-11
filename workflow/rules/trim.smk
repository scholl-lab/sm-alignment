import glob as _glob


# =============================================================================
# Trimming configuration
# =============================================================================
TRIM_CFG = config.get("trimming", {})
TRIM_FASTQ_DIRS = TRIM_CFG.get("fastq_dirs", [])
TRIM_FASTQ_PATTERN = TRIM_CFG.get("fastq_pattern", "*_R1_001.fastq.gz")
TRIMMED_DIR = os.path.join(OUTPUT_DIR, "bbduk_trimmed")

R1_RAW_SUFFIX = config.get("fastq", {}).get("r1_suffix", "_R1_001.fastq.gz")
R2_RAW_SUFFIX = config.get("fastq", {}).get("r2_suffix", "_R2_001.fastq.gz")


def _get_trim_samples():
    """Discover sample names from R1 FASTQ files across all input directories."""
    r1_files = []
    for d in TRIM_FASTQ_DIRS:
        r1_files.extend(_glob.glob(os.path.join(d, TRIM_FASTQ_PATTERN)))
    return [os.path.basename(f).replace(R1_RAW_SUFFIX, "") for f in r1_files]


def _find_raw_r1(wildcards):
    """Locate R1 FASTQ for a sample across input directories."""
    for d in TRIM_FASTQ_DIRS:
        candidate = os.path.join(d, f"{wildcards.sample}{R1_RAW_SUFFIX}")
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(
        f"No R1 FASTQ for '{wildcards.sample}' in: {TRIM_FASTQ_DIRS}"
    )


def _find_raw_r2(wildcards):
    """Locate R2 FASTQ for a sample across input directories."""
    for d in TRIM_FASTQ_DIRS:
        candidate = os.path.join(d, f"{wildcards.sample}{R2_RAW_SUFFIX}")
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(
        f"No R2 FASTQ for '{wildcards.sample}' in: {TRIM_FASTQ_DIRS}"
    )


# =============================================================================
# Rules
# =============================================================================
rule trim_all:
    input:
        expand(
            os.path.join(TRIMMED_DIR, "{sample}.bbduk_R1_001.fastq.gz"),
            sample=_get_trim_samples(),
        ),
        expand(
            os.path.join(TRIMMED_DIR, "{sample}.bbduk_R2_001.fastq.gz"),
            sample=_get_trim_samples(),
        ),


rule trim_adapters:
    """Adapter and quality trimming with BBDuk."""
    input:
        r1=_find_raw_r1,
        r2=_find_raw_r2,
    output:
        trimmed_r1=os.path.join(TRIMMED_DIR, "{sample}.bbduk_R1_001.fastq.gz"),
        trimmed_r2=os.path.join(TRIMMED_DIR, "{sample}.bbduk_R2_001.fastq.gz"),
    params:
        ref=TRIM_CFG.get("bbduk_ref", "adapters,artifacts"),
        ziplevel=TRIM_CFG.get("ziplevel", 5),
        ktrim=TRIM_CFG.get("ktrim", "r"),
        k=TRIM_CFG.get("k", 23),
        mink=TRIM_CFG.get("mink", 11),
        hdist=TRIM_CFG.get("hdist", 1),
        tpe=TRIM_CFG.get("tpe", "t"),
        tbo=TRIM_CFG.get("tbo", "t"),
        ftl=TRIM_CFG.get("ftl", 5),
        trimpolyg=TRIM_CFG.get("trimpolyg", 3),
        trimpolya=TRIM_CFG.get("trimpolya", 3),
        qtrim=TRIM_CFG.get("qtrim", "t"),
        trimq=TRIM_CFG.get("trimq", 10),
        quantize=TRIM_CFG.get("quantize", "0,10,20,30,40,50,60"),
    threads: 4
    resources:
        runtime=720,
    conda:
        "../envs/bbtools.yaml"
    log:
        os.path.join(LOG_DIR, "bbduk_trim.{sample}.log"),
    shell:
        r"""
        mkdir -p "$(dirname {output.trimmed_r1})"

        bbduk.sh \
            threads={threads} \
            in={input.r1} in2={input.r2} \
            out={output.trimmed_r1} out2={output.trimmed_r2} \
            ziplevel={params.ziplevel} \
            ref={params.ref} \
            ktrim={params.ktrim} k={params.k} mink={params.mink} hdist={params.hdist} \
            tpe={params.tpe} tbo={params.tbo} \
            ftl={params.ftl} trimpolyg={params.trimpolyg} trimpolya={params.trimpolya} \
            qtrim={params.qtrim} trimq={params.trimq} quantize={params.quantize} \
            > {log} 2>&1
        """
