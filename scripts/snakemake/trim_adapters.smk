import os
import glob
import functools

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths
prefix_results = functools.partial(os.path.join, 'results/exomes')
FASTQ_DIR = prefix_results('samtools_view_apa_genes_padding')  # Input directory containing FASTQ files
TRIMMED_DIR = prefix_results('samtools_view_apa_genes_padding')               # Output directory for trimmed FASTQ files
LOG_DIR = prefix_results('logs')                            # Directory for log files
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions to find all FASTQ files and extract sample names

def get_fastq_files():
    r1_files = glob.glob(os.path.join(FASTQ_DIR, "*_R1.fastq.gz"))
    return r1_files

def get_samples():
    r1_files = get_fastq_files()
    samples = [os.path.basename(f).replace('_R1.fastq.gz', '') for f in r1_files]
    return samples
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

rule all:
    input:
        expand(os.path.join(TRIMMED_DIR, "{sample}.bbduk_R1.fastq.gz"), sample=get_samples()),
        expand(os.path.join(TRIMMED_DIR, "{sample}.bbduk_R2.fastq.gz"), sample=get_samples())

rule trim_adapters:
    input:
        r1 = os.path.join(FASTQ_DIR, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(FASTQ_DIR, "{sample}_R2.fastq.gz")
    output:
        trimmed_r1 = os.path.join(TRIMMED_DIR, "{sample}.bbduk_R1.fastq.gz"),
        trimmed_r2 = os.path.join(TRIMMED_DIR, "{sample}.bbduk_R2.fastq.gz")
    log:
        log = os.path.join(LOG_DIR, "bbduk_trim.{sample}.log")
    threads: 8  # Number of threads to use for bbduk
    params:
        ref = "adapters"  # Reference file for adapters (ensure this file is accessible)
    conda:
        "base"
    shell:
        """
        # Create the output and log directories if they do not exist
        mkdir -p {TRIMMED_DIR}
        mkdir -p {LOG_DIR}

        # Run bbduk command
        bbduk.sh threads={threads} \
            in={input.r1} in2={input.r2} \
            out={output.trimmed_r1} out2={output.trimmed_r2} \
            ref={params.ref} \
            ktrim=r k=23 mink=11 hdist=1 tpe tbo \
            > {log.log} 2>&1
        """
# ----------------------------------------------------------------------------------- #
