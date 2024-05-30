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
BQSR_DIR = prefix_results('bqsr')
LOG_DIR = prefix_results('logs')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper functions to find all deduplicated BAM files

def get_dedup_files():
    return glob.glob(DEDUP_DIR + "/*.merged.dedup.bam")

# Function to return the memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the pipeline rules

# Define the rule to collect all the BQSR targets
rule all:
    input:
        expand("results/bqsr/{sample}.merged.dedup.bqsr.bam", sample=[os.path.basename(x).replace('.merged.dedup.bam','') for x in get_dedup_files()])

# Define the rule to perform base recalibration using BaseRecalibrator
rule base_recalibration:
    input:
        bam_file = os.path.join(DEDUP_DIR, '{sample}.merged.dedup.bam')
    output:
        recal_table = os.path.join(BQSR_DIR, '{sample}.merged.dedup.recal_data.table')
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory
    conda:
        "gatk"
    log:
        recal = os.path.join(LOG_DIR, 'recal.gatk.{sample}.log'),
    shell:
        """
        # Create the BQSR directory if it does not exist
        mkdir -p {BQSR_DIR}
        
        # Use GATK's BaseRecalibrator to perform base recalibration
        gatk --java-options '-Xms4000m -Xmx7g -Djava.io.tmpdir={resources.tmpdir}' BaseRecalibrator \
            -I "{input.bam_file}" \
            -R analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
            --known-sites analysis/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
            --known-sites analysis/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
            --known-sites analysis/GATK_resource_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            -O "{output.recal_table}" 2> {log.recal}
        """

# Define the rule to apply BQSR using ApplyBQSR
rule apply_bqsr:
    input:
        bam_file = os.path.join(DEDUP_DIR, '{sample}.merged.dedup.bam'),
        recal_table = os.path.join(BQSR_DIR, '{sample}.merged.dedup.recal_data.table')
    output:
        bqsr_bam = os.path.join(BQSR_DIR, '{sample}.merged.dedup.bqsr.bam')
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory
    conda:
        "gatk"
    log:
        apply_bqsr = os.path.join(LOG_DIR, 'apply_bqsr.gatk.{sample}.log'),
    shell:
        """
        # Use GATK's ApplyBQSRS to apply base recalibration
        gatk --java-options '-Xms4000m -Xmx7g -Djava.io.tmpdir={resources.tmpdir} -Dsamjdk.compression_level=6' ApplyBQSR \
            -I "{input.bam_file}" \
            -bqsr "{input.recal_table}" \
            -O "{output.bqsr_bam}" 2> {log.apply_bqsr}
        """
# ----------------------------------------------------------------------------------- #