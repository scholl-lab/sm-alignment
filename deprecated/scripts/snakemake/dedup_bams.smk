import os
import glob
import functools

# Read the configuration file
configfile: "config.yaml"

# ----------------------------------------------------------------------------------- #
# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, config['output_folder'])
MERGE_DIR = prefix_results('merged')
DEDUP_DIR = prefix_results('dedup')
LOG_DIR = prefix_results('logs')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions to find all merged BAM files

def get_merged_files():
    return glob.glob(MERGE_DIR + "/*.merged.bam")

# Function to return the number of threads for MarkDuplicatesSpark
def get_dedup_threads(wildcards, threads):
    return threads

# Function to return the memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

# Define the rule to collect all the dedup targets
# issues with disk space fix:
# https://gatk.broadinstitute.org/hc/en-us/community/posts/360067258451-MarkDuplicatesSpark-doesn-t-work-for-large-bam-files
rule all:
    input:
        expand("results/dedup/{sample}.merged.dedup.bam", sample=[os.path.basename(x).replace('.merged.bam','') for x in get_merged_files()])

# Define the rule to deduplicate BAM files
# maybe have to set java options temp https://hpc.nih.gov/apps/GATK.html
rule deduplicate_bam_files:
    input:
        bam_file = os.path.join(MERGE_DIR, '{sample}.merged.bam')
    output:
        dedup_bam = os.path.join(DEDUP_DIR, '{sample}.merged.dedup.bam'),
        metrics = os.path.join(DEDUP_DIR, '{sample}.merged.dedup_metrics.txt')
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,                # Temporary directory
    conda:
        "gatk"  # This sets the Conda environment for this rule
    log:
        dedup = os.path.join(LOG_DIR, 'dedup.gatk.{sample}.log'),
    shell:
        """
        # Create the dedup directory if it does not exist
        mkdir -p results/dedup
        
        # Use GATK to deduplicate the BAM file
        gatk --java-options '-Xms4000m -Xmx7g -Djava.io.tmpdir={resources.tmpdir}' MarkDuplicates \
            -I "{input.bam_file}" \
            -O "{output.dedup_bam}" \
            -M "{output.metrics}" \
            --CREATE_INDEX \
            --VALIDATION_STRINGENCY SILENT 2> {log.dedup}
        """
# ----------------------------------------------------------------------------------- #
