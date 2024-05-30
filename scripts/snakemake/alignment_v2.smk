import os
import functools
import csv
import yaml

# Load configuration file containing user-defined settings
configfile: "config.yaml"

# ----------------------------------------------------------------------------------- #
# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract user-defined input and output directories and reference file from the configuration file
with open("config.yaml", "r") as stream:
    config = yaml.safe_load(stream)

INPUT_DIR = config["fastq_folder"]
OUTPUT_DIR = config["output_folder"]
REFERENCE_FILE = config["reference"]
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Read metadata table into a dictionary from the TSV file
metadata_dict = {}
with open('metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        basename = row['fastq_files_basename']
        metadata_dict[basename] = row
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, config['output_folder'])
ALIGNED_DIR = prefix_results('aligned')
MERGE_DIR = prefix_results('merged')
MD_DIR = prefix_results('mark_duplicates')
LOG_DIR = prefix_results('logs')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions

# Function to return the number of threads for bwa mem
def get_bwa_threads(wildcards, threads):
    return threads - 2

# Function to return the number of threads for samtools sort
def get_sort_threads(wildcards, threads):
    return 2

# Function to return the memory for samtools sort in MB
def get_sort_mem(wildcards, threads):
    return 2 * 1200

# Function to return the memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 1200
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

# Rule "all": Defines the final output of the pipeline
rule all:
    input:
        expand(f"{ALIGNED_DIR}/{{fastq_files_basename}}_{{subfolder}}.bam", 
               zip, 
               fastq_files_basename=metadata_dict.keys(), 
               subfolder=[metadata_dict[key]['subfolder'] for key in metadata_dict.keys()])

# Rule "bwa_map": Performs BWA alignment and SAMtools sorting
rule bwa_map:
    input:
        fastq_f = lambda wildcards: f"{INPUT_DIR}/{metadata_dict[wildcards.fastq_files_basename]['subfolder']}/{wildcards.fastq_files_basename}_R1_001.fastq.gz",
        fastq_r = lambda wildcards: f"{INPUT_DIR}/{metadata_dict[wildcards.fastq_files_basename]['subfolder']}/{wildcards.fastq_files_basename}_R2_001.fastq.gz",
    output:
        bam_lane = f"{ALIGNED_DIR}/{{fastq_files_basename}}_{{subfolder}}.bam",
    params:
        reference = REFERENCE_FILE,
        read_group = lambda wildcards: f'"@RG\\tID:{metadata_dict[wildcards.fastq_files_basename]["lane"]}-{metadata_dict[wildcards.fastq_files_basename]["project_sample"]}\\tSM:{metadata_dict[wildcards.fastq_files_basename]["project_sample"]}\\tLB:{metadata_dict[wildcards.fastq_files_basename]["project_sample"]}\\tPL:ILLUMINA\\tPU:{metadata_dict[wildcards.fastq_files_basename]["lane"]}-{metadata_dict[wildcards.fastq_files_basename]["mdc_project"]}"',
    threads: 16,
    resources:
        mem_mb = get_mem_from_threads, # Memory in MB based on the number of threads
        time = '24:00:00',              # Time limit for the job
        tmpdir = SCRATCH_DIR,           # Temporary directory
        bwa_threads = get_bwa_threads,  # Number of threads for BWA
        sort_threads = get_sort_threads,# Number of threads for SAMtools sort
        sort_mem = get_sort_mem,        # Memory for SAMtools sort in MB
    log:
        bwa = f"{LOG_DIR}/map.bwa.{{fastq_files_basename}}_{{subfolder}}.log",
        samtools = f"{LOG_DIR}/map.samtools.{{fastq_files_basename}}_{{subfolder}}.log",
    shell:
        """
        bwa mem -t {resources.bwa_threads} -R {params.read_group} {params.reference} {input.fastq_f} {input.fastq_r} 2> {log.bwa} | \
        samtools sort -@{resources.sort_threads} -m{resources.sort_mem}M -O BAM -T {resources.tmpdir} -o {output.bam_lane} - 2> {log.samtools}
        """
# ----------------------------------------------------------------------------------- #
