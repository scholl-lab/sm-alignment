import os
import glob
import functools

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths
prefix_results = functools.partial(os.path.join, 'results/exomes')
BAM_DIR = prefix_results('bqsr')
OUTPUT_DIR = prefix_results('samtools_view_apa_genes_padding')
LOG_DIR = prefix_results('logs')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions to find all BAM files

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
        expand(os.path.join(OUTPUT_DIR, "{sample}.apa-genes.bam"), sample=get_samples())

rule samtools_view_apa_genes:
    input:
        bam_file = os.path.join(BAM_DIR, '{sample}.bam'),
        bed_file = '_other/apa-genes.sorted.padding1000.bed'
    output:
        bam = os.path.join(OUTPUT_DIR, '{sample}.apa-genes.bam')
    log:
        log = os.path.join(LOG_DIR, 'samtools_view.{sample}.log')
    threads: 1
    conda:
        "base"
    shell:
        """
        # Create the output directory if it does not exist
        mkdir -p {OUTPUT_DIR}

        # Run samtools view command
        samtools view -b -h -P -M -L {input.bed_file} {input.bam_file} > {output.bam} 2> {log.log}
        """
# ----------------------------------------------------------------------------------- #
