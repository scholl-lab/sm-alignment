import os
import csv
import yaml
import functools

# Read the configuration file
configfile: "config.yaml"

# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')

# Extract user-defined input and output directories from the configuration file
with open("config.yaml", "r") as stream:
    config = yaml.safe_load(stream)

OUTPUT_DIR = config["output_folder"]

# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, OUTPUT_DIR)
ALIGNED_DIR = prefix_results('aligned')
MERGE_DIR = prefix_results('merged')
LOG_DIR = prefix_results('logs')

# Read metadata table into a dictionary from the TSV file
metadata_dict = {}
with open('metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        basename = row['fastq_files_basename']
        metadata_dict[basename] = row

# Helper function to get unique samples
def get_samples():
    samples = set(row['project_sample'] for row in metadata_dict.values())
    return list(samples)

# Function to return the number of threads for samtools merge
def get_merge_threads(wildcards, threads):
    return threads

# Function to return the memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 1200

# List of unique samples
samples = get_samples()

# Define the pipeline rules
rule all:
    input:
        expand(f"{MERGE_DIR}/{{sample}}.merged.bam", sample=samples)

rule merge_bam_files:
    input:
        bam_files = lambda wildcards: [
            f"{ALIGNED_DIR}/{metadata_dict[basename]['subfolder']}/{basename}.bam"
            for basename in metadata_dict
            if metadata_dict[basename]['project_sample'] == wildcards.sample
        ]
    output:
        bam = os.path.join(MERGE_DIR, '{sample}.merged.bam')
    params:
        list_file = os.path.join(MERGE_DIR, '{sample}.bamlist')
    threads: 8
    resources:
        mem_mb = get_mem_from_threads,
        time = '24:00:00',
        tmpdir = SCRATCH_DIR,
        merge_threads = get_merge_threads
    log:
        merge = os.path.join(LOG_DIR, 'merge.samtools.{sample}.log')
    shell:
        """
        # Write the BAM file names to merge into a list file
        echo "{input.bam_files}" | tr " " "\\n" > "{params.list_file}"
        
        # Create the merged directory if it does not exist
        mkdir -p {MERGE_DIR}
        
        # Use samtools to merge the BAM files
        samtools merge -@{resources.merge_threads} -O BAM -b "{params.list_file}" "{output.bam}" 2> {log.merge}
        
        # Remove the list file
        rm "{params.list_file}"
        """
