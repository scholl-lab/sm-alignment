import glob
import re
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
# TODO: has to be adapted to the config file
ALIGNED_DIR = config.get("aligned_folder", "results/aligned/")
MERGE_DIR = prefix_results('merged')
LOG_DIR = prefix_results('logs')
# ----------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------- #
# Helper functions

# Function to find all unique sample prefixes
def get_samples():
    samples = set()
    for file in glob.glob(ALIGNED_DIR + "*.bam"):
        # Extract just the filename without the directory
        filename = file.replace(ALIGNED_DIR, "")
        
        m = re.match(r"(.+)_DNA_(\d+)_(.+)_S\d+_L\d+_lane\d+\.bam", filename)
        if m:
            samples.add(f"{m.group(1)}_DNA_{m.group(2)}_{m.group(3)}")
    return list(samples)

# Function to return the number of threads for samtools sort
def get_merge_threads(wildcards, threads):
    return threads

# Function to return the memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 1200
# ----------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------- #
# List of unique samples
samples = get_samples()
# ----------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

# Define the rule to collect all the targets
rule all:
    input:
        expand("results/merged/{sample}.merged.bam", sample=samples)

# Define the rule to merge BAM files
rule merge_bam_files:
    input:
        bam_files = lambda wildcards: glob.glob(ALIGNED_DIR + wildcards.sample + "_*.bam")
    output:
        bam = os.path.join(MERGE_DIR, '{sample}.merged.bam'),
    params:
        list_file = os.path.join(MERGE_DIR, '{sample}.bamlist'),
    threads: 8
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '24:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory
        merge_threads = get_merge_threads,  # Number of threads for SAMtools sort
    log:
        merge = os.path.join(LOG_DIR, 'merge.samtools.{sample}.log'),
    shell:
        """
        # Write the BAM file names to merge into a list file
        echo "{input.bam_files}" | tr " " "\\n" > "{params.list_file}"
        
        # Create the merged directory if it does not exist
        mkdir -p results/merged
        
        # Use samtools to merge the BAM files
        samtools merge -@{resources.merge_threads} -O BAM -b "{params.list_file}" "{output.bam}" 2> {log.merge}
        
        # Remove the list file
        rm "{params.list_file}"
        """
# ----------------------------------------------------------------------------------- #