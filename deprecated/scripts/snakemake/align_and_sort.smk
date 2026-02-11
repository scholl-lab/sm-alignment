import os
import glob
import functools

# ----------------------------------------------------------------------------------- #
# Read the configuration file
configfile: "config.yaml"

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths
prefix_results = functools.partial(os.path.join, 'results/exomes')
INPUT_DIR = prefix_results('samtools_view_apa_genes_padding')  # Directory containing trimmed FASTQ files
OUTPUT_DIR = prefix_results('samtools_view_apa_genes_padding')  # Output directory (same as input)
LOG_DIR = prefix_results('logs')  # Directory for log files
REFERENCE = config["reference"]  # Reference genome path from config.yaml
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions to find all trimmed FASTQ R1 files and extract sample names

def get_trimmed_fastq_files():
    # Match filenames like {sample}.bbduk_R1.fastq.gz
    return glob.glob(os.path.join(INPUT_DIR, "*.bbduk_R1.fastq.gz"))

def get_samples():
    r1_files = get_trimmed_fastq_files()
    # Extract sample names by removing '.bbduk_R1.fastq.gz'
    samples = [os.path.basename(f).replace('.bbduk_R1.fastq.gz', '') for f in r1_files]
    return samples
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{sample}.sorted.bam"), sample=get_samples()),
        expand(os.path.join(OUTPUT_DIR, "{sample}.sorted.bam.bai"), sample=get_samples())

rule align_and_sort:
    input:
        r1 = os.path.join(INPUT_DIR, "{sample}.bbduk_R1.fastq.gz"),
        r2 = os.path.join(INPUT_DIR, "{sample}.bbduk_R2.fastq.gz"),
        reference = REFERENCE
    output:
        sorted_bam = os.path.join(OUTPUT_DIR, "{sample}.sorted.bam"),
        bwa_log = os.path.join(LOG_DIR, "{sample}.bwa.log"),
        samblaster_log = os.path.join(LOG_DIR, "{sample}.samblaster.log"),
        sambamba_view_log = os.path.join(LOG_DIR, "{sample}.sambamba_view.log"),
        samtools_sort_log = os.path.join(LOG_DIR, "{sample}.samtools_sort.log")
    params:
        # Constructing the read group (RG) information based on sample name
        read_group = lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA\\tPU:{wildcards.sample}"
    threads: 8
    resources:
        mem_mb = 15000  # 15 GB per job
    conda:
        "base"
    shell:
        """
        # Run the alignment and sorting pipeline using samtools sort
        bwa mem -C -t {threads} \
            -R "{params.read_group}" \
            {input.reference} \
            {input.r1} {input.r2} \
            2> {output.bwa_log} \
        | samblaster 2> {output.samblaster_log} \
        | samtools view -@ {threads} -bS -o /dev/stdout /dev/stdin \
        | samtools sort -@ {threads} -o {output.sorted_bam} \
            2> {output.samtools_sort_log}

        # Placeholder log for samtools view (no actual stderr captured)
        echo "Samtools view completed." > {output.sambamba_view_log}
        """

rule index_bam:
    input:
        sorted_bam = os.path.join(OUTPUT_DIR, "{sample}.sorted.bam")
    output:
        index = os.path.join(OUTPUT_DIR, "{sample}.sorted.bam.bai")
    log:
        log = os.path.join(LOG_DIR, "{sample}.samtools_index.log")
    threads: 1
    resources:
        mem_mb = 15000  # 15 GB per job
    conda:
        "base"
    shell:
        """
        # Index the sorted BAM file using samtools
        samtools index -@ {threads} {input.sorted_bam} \
            2> {log.log}
        """
# ----------------------------------------------------------------------------------- #
