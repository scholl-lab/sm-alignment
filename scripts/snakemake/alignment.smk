import os
import functools

# Load configuration file containing user-defined settings
configfile: "config.yaml"

# examples:
# https://github.com/Boyle-Lab/short-read-seq_snakemake/blob/main/Snakefile

# ----------------------------------------------------------------------------------- #
# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract user-defined input and output directories and reference file from the configuration file
INPUT_DIR = config["fastq_folder"]
OUTPUT_DIR = config["output_folder"]
REFERENCE_FILE = config["reference"]
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
# TODO: remove if not needed

# Helper function to obtain a list of all forward .fastq.gz files in a folder and its subfolders
def get_input_files(folder):
    input_files = []
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith('_R1_001.fastq.gz'):
                input_files.append(os.path.relpath(os.path.join(root, file), folder))
    return input_files

# Define a function to get the readgroup (RG) string from the filename and subfolder
def get_rg(filename):
    print(filename)
    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    sample_id = basename.split("_")[3]
    lane = "-".join(basename.split("_")[4:6]) + '-' + dirname
    read_group = '@RG\tID:' + lane + '-' + sample_id + '\tSM:' + sample_id + '\tLB:' + sample_id + '\tPL:ILLUMINA\tPU:' + lane + '-' + sample_id
    return read_group

# Function to get the lane bam from the input fastq filename and subfolder
def get_bam_lane_name(filename):
    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    sample_id = basename.split("_")[3]
    bam_lane_basename = basename.replace("_R1_001.fastq.gz", "") + '_' + dirname
    return bam_lane_basename

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
# Generate file lists using the get_input_files function
input_fq_f = get_input_files(INPUT_DIR)
# ----------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

# Rule "all": Defines the final output of the pipeline, i.e., the aligned BAM files
rule all:
    input:
        # lane BAM output
        [os.path.join(ALIGNED_DIR, '{}.bam'.format(get_bam_lane_name(x))) for x in input_fq_f]

# Rule "bwa_map": Performs BWA alignment and SAMtools sorting for the given paired-end reads
rule bwa_map:
    input:
        fastq_f = "download/{run}/{mdc_project}_DNA_{dna_number}_{project_sample}_{sample_sheet_number}_{lane}_R1_001.fastq.gz",
        fastq_r = "download/{run}/{mdc_project}_DNA_{dna_number}_{project_sample}_{sample_sheet_number}_{lane}_R2_001.fastq.gz",
    output:
        bam_lane = os.path.join(ALIGNED_DIR, '{mdc_project}_DNA_{dna_number}_{project_sample}_{sample_sheet_number}_{lane}_{run}.bam')
    params:
        reference = REFERENCE_FILE,
        read_group = '"@RG\\tID:{run}-{sample_sheet_number}-{lane}-{project_sample}\\tSM:{project_sample}\\tLB:{project_sample}\\tPL:ILLUMINA\\tPU:{run}-{sample_sheet_number}-{lane}"',
    threads: 16
    resources:
        mem_mb = get_mem_from_threads, # Memory in MB based on the number of threads
        time = '24:00:00',              # Time limit for the job
        tmpdir = SCRATCH_DIR,           # Temporary directory
        bwa_threads = get_bwa_threads,  # Number of threads for BWA
        sort_threads = get_sort_threads,# Number of threads for SAMtools sort
        sort_mem = get_sort_mem,        # Memory for SAMtools sort in MB
    log:
        bwa = os.path.join(LOG_DIR, 'map.bwa.{mdc_project}_DNA_{dna_number}_{project_sample}_{sample_sheet_number}_{lane}_{run}.log'),
        samtools = os.path.join(LOG_DIR, 'map.samtools.{mdc_project}_DNA_{dna_number}_{project_sample}_{sample_sheet_number}_{lane}_{run}.log')
    shell:
        """
        bwa mem -t {resources.bwa_threads} -R {params.read_group} {params.reference} {input.fastq_f} {input.fastq_r} 2> {log.bwa} | samtools sort -@{resources.sort_threads} -m{resources.sort_mem}M -O BAM -T {resources.tmpdir} -o {output.bam_lane} - 2> {log.samtools}
        """
# ----------------------------------------------------------------------------------- #