#####################################################################
# trim_adapters.smk
#####################################################################
import os
import glob

# Load the config file
configfile: "config.yaml"

# Read parameters from config
FASTQ_DIRS       = config["params"]["fastq_dirs"]   # multiple input directories
PREFIX_RESULTS   = config["params"]["prefix_results"]
TRIMMED_DIR_SUB  = config["params"]["trimmed_dir_sub"]
LOG_DIR_SUB      = config["params"]["log_dir_sub"]
FASTQ_PATTERN    = config["params"]["fastq_pattern"]
BBDUK_REF        = config["params"]["bbduk_ref"]
THREADS          = config["params"]["threads"]

# Additional BBDuk parameters
ZIPLEVEL   = config["params"]["ziplevel"]
KTRIM      = config["params"]["ktrim"]
K          = config["params"]["k"]
MINK       = config["params"]["mink"]
HDIST      = config["params"]["hdist"]
TPE        = config["params"]["tpe"]
TBO        = config["params"]["tbo"]
FTL        = config["params"]["ftl"]
TRIMPOLYG  = config["params"]["trimpolyg"]
TRIMPOLYA  = config["params"]["trimpolya"]

#####################################################################
# Define output directories
#####################################################################
TRIMMED_DIR = os.path.join(PREFIX_RESULTS, TRIMMED_DIR_SUB)
LOG_DIR     = os.path.join(PREFIX_RESULTS, LOG_DIR_SUB)

#####################################################################
# Helper functions
#####################################################################
def get_fastq_files():
    """
    Return all R1 FASTQ files matching FASTQ_PATTERN
    in all directories listed in FASTQ_DIRS.
    """
    r1_files = []
    print("DEBUG: Looking for R1 files with pattern:", FASTQ_PATTERN)
    for d in FASTQ_DIRS:
        print(f"DEBUG: Searching directory: {d}")
        found = glob.glob(os.path.join(d, FASTQ_PATTERN))
        if found:
            print(f"DEBUG: Found {len(found)} file(s) in {d}: {found}")
        else:
            print(f"DEBUG: Found no files in {d} with pattern {FASTQ_PATTERN}")
        r1_files.extend(found)

    print("DEBUG: Combined R1 file list:", r1_files)
    return r1_files

def get_samples():
    """
    Identify sample names by removing '_R1_001.fastq.gz'
    from each R1 file's base name.
    """
    r1_files = get_fastq_files()
    samples = [
        os.path.basename(f).replace('_R1_001.fastq.gz', '')
        for f in r1_files
    ]
    print("DEBUG: Extracted sample names:", samples)
    return samples

def find_r1(wildcards):
    """
    Given a sample name, locate its R1 FASTQ file by
    searching each directory in FASTQ_DIRS.
    We expect the file to end with '_R1_001.fastq.gz'.
    """
    sample = wildcards.sample
    for d in FASTQ_DIRS:
        candidate = os.path.join(d, f"{sample}_R1_001.fastq.gz")
        if os.path.exists(candidate):
            print(f"DEBUG: Found R1 for {sample} in {d}: {candidate}")
            return candidate
    raise FileNotFoundError(
        f"No R1 FASTQ found for sample '{sample}' in directories: {FASTQ_DIRS}"
    )

def find_r2(wildcards):
    """
    Given a sample name, locate its R2 FASTQ file by
    searching each directory in FASTQ_DIRS.
    We expect the file to end with '_R2_001.fastq.gz'.
    """
    sample = wildcards.sample
    for d in FASTQ_DIRS:
        candidate = os.path.join(d, f"{sample}_R2_001.fastq.gz")
        if os.path.exists(candidate):
            print(f"DEBUG: Found R2 for {sample} in {d}: {candidate}")
            return candidate
    raise FileNotFoundError(
        f"No R2 FASTQ found for sample '{sample}' in directories: {FASTQ_DIRS}"
    )

#####################################################################
# Workflow rules
#####################################################################
all_samples = get_samples()
print("DEBUG: Final sample list (all_samples):", all_samples)

rule all:
    input:
        expand(os.path.join(TRIMMED_DIR, "{sample}.bbduk_R1_001.fastq.gz"), sample=all_samples),
        expand(os.path.join(TRIMMED_DIR, "{sample}.bbduk_R2_001.fastq.gz"), sample=all_samples)

rule trim_adapters:
    input:
        r1 = find_r1,
        r2 = find_r2
    output:
        trimmed_r1 = os.path.join(TRIMMED_DIR, "{sample}.bbduk_R1_001.fastq.gz"),
        trimmed_r2 = os.path.join(TRIMMED_DIR, "{sample}.bbduk_R2_001.fastq.gz")
    log:
        log = os.path.join(LOG_DIR, "bbduk_trim.{sample}.log")
    threads:
        THREADS
    params:
        ref = BBDUK_REF
    conda:
        "base"
    shell:
        """
        # Create trimmed output + log directories
        mkdir -p {TRIMMED_DIR} {LOG_DIR}

        # Run BBDuk
        bbduk.sh \
            threads={threads} \
            in={input.r1} in2={input.r2} \
            out={output.trimmed_r1} out2={output.trimmed_r2} \
            ziplevel={ZIPLEVEL} \
            ref={params.ref} \
            ktrim={KTRIM} k={K} mink={MINK} hdist={HDIST} \
            tpe={TPE} tbo={TBO} \
            ftl={FTL} trimpolyg={TRIMPOLYG} trimpolya={TRIMPOLYA} \
            > {log.log} 2>&1
        """
