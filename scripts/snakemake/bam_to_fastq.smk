import os
import glob
import functools

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths
prefix_results = functools.partial(os.path.join, 'results/exomes')
BAM_DIR = prefix_results('samtools_view_apa_genes_padding')  # Input directory containing BAM files
FASTQ_DIR = prefix_results('samtools_view_apa_genes_padding')  # Output directory for FASTQ files
LOG_DIR = prefix_results('logs')  # Directory for log files
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions to find all BAM files and extract sample names

def get_bam_files():
    return glob.glob(os.path.join(BAM_DIR, "*.bam"))

def get_samples():
    bam_files = get_bam_files()
    samples = [os.path.basename(bam).replace('.bam', '') for bam in bam_files]
    return samples
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

rule all:
    input:
        expand(os.path.join(FASTQ_DIR, "{sample}_R1.fastq.gz"), sample=get_samples()),
        expand(os.path.join(FASTQ_DIR, "{sample}_R2.fastq.gz"), sample=get_samples())

rule convert_bam_to_fastq:
    input:
        bam_file = os.path.join(BAM_DIR, '{sample}.bam')
    output:
        r1 = os.path.join(FASTQ_DIR, '{sample}_R1.fastq.gz'),
        r2 = os.path.join(FASTQ_DIR, '{sample}_R2.fastq.gz')
    log:
        log = os.path.join(LOG_DIR, 'samtools_fastq.{sample}.log')
    threads: 8
    conda:
        "base"
    shell:
        """
        # Create the output and log directories if they do not exist
        mkdir -p {FASTQ_DIR}
        mkdir -p {LOG_DIR}

        # Run samtools fastq command to convert BAM to FASTQ
        # first sort the bam file by name
        samtools sort -n -@ {threads} {input.bam_file} | \
        samtools fastq -@ {threads} - \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            2> {log.log}
        """
# ----------------------------------------------------------------------------------- #
