#####################################################################
# alignment_pipeline.smk
#####################################################################
import os
import glob
import csv
import yaml

##############################################################################
# 1) Load the config file
##############################################################################
configfile: "config.yaml"

# We'll assume your config has top-level keys:
#   fastq_folder
#   aligned_folder
#   output_folder
#   final_bam_folder
#   reference
#   reference_unpacked
#   panel_of_normals
#   af_only_gnomad
#   common_biallelic_gnomad
#   reference_version
#   mutect_scatter_by_chromosome
#   final_bam_file_extension
#   log_dir_sub

FASTQ_FOLDER            = config["fastq_folder"]
ALIGNED_FOLDER          = config["aligned_folder"]
OUTPUT_FOLDER           = config["output_folder"]
FINAL_BAM_FOLDER        = config["final_bam_folder"]
REFERENCE               = config["reference"]
REFERENCE_UNPACKED      = config["reference_unpacked"]
PANEL_OF_NORMALS        = config["panel_of_normals"]
AF_ONLY_GNOMAD          = config["af_only_gnomad"]
COMMON_BIALLELIC_GNOMAD = config["common_biallelic_gnomad"]
REFERENCE_VERSION       = config["reference_version"]
FINAL_BAM_EXTENSION     = config["final_bam_file_extension"]
LOG_SUBFOLDER           = config["log_dir_sub"]

##############################################################################
# 2) Optional environment-based scratch directory (e.g. for cluster jobs)
##############################################################################
SCRATCH_DIR = os.environ.get("TMPDIR", "/tmp")

##############################################################################
# 3) Read the metadata (metadata.tsv)
##############################################################################
metadata_dict = {}
METADATA_FILE = config["metadata_file"]

def load_metadata():
    print(f"DEBUG: Loading metadata from {METADATA_FILE}")
    with open(METADATA_FILE, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            basename = row["fastq_files_basename"]
            metadata_dict[basename] = row

load_metadata()
print(f"DEBUG: Number of entries in metadata: {len(metadata_dict)}")

##############################################################################
# 4) Define helper functions for sample-lists and resource usage
##############################################################################
def get_samples():
    """
    Return a list of unique sample names (project_sample) from metadata.
    """
    samples = set(md["project_sample"] for md in metadata_dict.values())
    samples_list = sorted(samples)
    print(f"DEBUG: Found {len(samples_list)} unique sample(s): {samples_list}")
    return samples_list

def get_mem_from_threads(wildcards, threads):
    """Default: 1200 MB per thread."""
    return threads * 1200

def get_bwa_threads(wildcards, threads):
    """
    e.g. if threads=9, then use 8 for BWA, 1 for samtools sort,
    so it doesn't slow down sorting too much.
    """
    return max(1, threads - 1)

def get_sort_threads(wildcards, threads):
    # With threads=9, this returns 1 for samtools sort
    return 1

def get_sort_mem(wildcards, threads):
    return 2 * 1200

all_samples = get_samples()

##############################################################################
# 5) Make sure directories exist
##############################################################################
os.makedirs(ALIGNED_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER,  exist_ok=True)
os.makedirs(FINAL_BAM_FOLDER, exist_ok=True)

LOG_DIR = os.path.join(OUTPUT_FOLDER, LOG_SUBFOLDER)
os.makedirs(LOG_DIR, exist_ok=True)

##############################################################################
# 6) Final outputs (the BQSR BAMs) & "all" rule
##############################################################################
def final_bam_path(sample):
    return os.path.join(FINAL_BAM_FOLDER, f"{sample}{FINAL_BAM_EXTENSION}")

rule all:
    input:
        [final_bam_path(s) for s in all_samples]

##############################################################################
# 7) ALIGNMENT RULE (BWA + SAMtools sort)
##############################################################################
def find_r1(basename):
    """ e.g. FASTQ_FOLDER/basename.bbduk_R1_001.fastq.gz """
    return os.path.join(FASTQ_FOLDER, f"{basename}.bbduk_R1_001.fastq.gz")

def find_r2(basename):
    """ e.g. FASTQ_FOLDER/basename.bbduk_R2_001.fastq.gz """
    return os.path.join(FASTQ_FOLDER, f"{basename}.bbduk_R2_001.fastq.gz")

rule bwa_map:
    """
    Align per 'basename' -> produce ALIGNED_FOLDER/{basename}.bam
    """
    input:
        lambda wc: [find_r1(wc.basename), find_r2(wc.basename)]
    output:
        bam_lane = os.path.join(ALIGNED_FOLDER, "{basename}.bam")
    params:
        reference = REFERENCE,
        read_group = lambda wc: (
            "\"@RG\\tID:{lane}-{sample}\\tSM:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tPU:{lane}-{mdc_project}\""
            .format(
                lane=metadata_dict[wc.basename]["lane"],
                sample=metadata_dict[wc.basename]["project_sample"],
                mdc_project=metadata_dict[wc.basename]["mdc_project"]
            )
        )
    threads: 9
    resources:
        mem_mb       = get_mem_from_threads,
        time         = "24:00:00",
        tmpdir       = SCRATCH_DIR,
        bwa_threads  = get_bwa_threads,
        sort_threads = get_sort_threads,
        sort_mem     = get_sort_mem
    conda:
        "base"
    log:
        bwa      = os.path.join(LOG_DIR, "map.bwa.{basename}.log"),
        samtools = os.path.join(LOG_DIR, "map.samtools.{basename}.log")
    shell:
        r"""
        echo "DEBUG: Starting alignment for {wildcards.basename}" >&2

        # Create a unique subfolder for samtools sort temp files:
        TMP_SORT_DIR=$(mktemp -p {resources.tmpdir} -d samtools-sort.XXXXXX)
        echo "DEBUG: Created sort tempdir: $TMP_SORT_DIR" >&2

        # Run BWA + pipe to samtools sort
        bwa mem -t {resources.bwa_threads} \
            -R {params.read_group} \
            {params.reference} \
            {input[0]} {input[1]} \
            2> {log.bwa} \
        | samtools sort -@ {resources.sort_threads} \
            -m {resources.sort_mem}M \
            -O BAM \
            -T "$TMP_SORT_DIR"/tmp \
            -o {output.bam_lane} \
            2> {log.samtools}

        # Clean up temporary directory
        rm -rf "$TMP_SORT_DIR"
        echo "DEBUG: Removed tmpdir: $TMP_SORT_DIR" >&2

        echo "DEBUG: Finished alignment for {wildcards.basename}" >&2
        """

##############################################################################
# 8) MERGE per sample
##############################################################################
def get_merged_bam(sample):
    return os.path.join(OUTPUT_FOLDER, "merged", f"{sample}.merged.bam")

rule merge_bam_files:
    """
    Merge all lane-level BAMs that correspond to the same 'project_sample'.
    """
    input:
        lambda wc: [
            os.path.join(ALIGNED_FOLDER, f"{basename}.bam")
            for basename in metadata_dict
            if metadata_dict[basename]["project_sample"] == wc.sample
        ]
    output:
        merged_bam = get_merged_bam("{sample}")
    params:
        list_file = os.path.join(OUTPUT_FOLDER, "merged", "{sample}.bamlist")
    threads: 8
    resources:
        mem_mb = get_mem_from_threads,  
        time   = "24:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "base"
    log:
        merge = os.path.join(LOG_DIR, "merge.samtools.{sample}.log")
    shell:
        r"""
        echo "DEBUG: Merging BAMs for sample {wildcards.sample}" >&2
        mkdir -p {OUTPUT_FOLDER}/merged
        echo "{input}" | tr " " "\n" > "{params.list_file}"

        samtools merge -@ {threads} -O BAM -b "{params.list_file}" "{output.merged_bam}" 2> {log.merge}

        rm "{params.list_file}"
        echo "DEBUG: Finished merging for {wildcards.sample}" >&2
        """

##############################################################################
# 9) DEDUPLICATION
##############################################################################
def get_dedup_bam(sample):
    dd_dir = os.path.join(OUTPUT_FOLDER, "dedup")
    os.makedirs(dd_dir, exist_ok=True)
    return os.path.join(dd_dir, f"{sample}.merged.dedup.bam")

rule deduplicate_bam_files:
    input:
        merged_bam = get_merged_bam("{sample}")
    output:
        dedup_bam = get_dedup_bam("{sample}"),
        metrics   = os.path.join(OUTPUT_FOLDER, "dedup", "{sample}.merged.dedup_metrics.txt")
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 4400,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        dedup = os.path.join(LOG_DIR, "dedup.gatk.{sample}.log")
    shell:
        r"""
        echo "DEBUG: MarkDuplicates for sample {wildcards.sample}" >&2
        gatk --java-options "-Xms4000m -Xmx7g -Djava.io.tmpdir={resources.tmpdir}" MarkDuplicates \
            -I "{input.merged_bam}" \
            -O "{output.dedup_bam}" \
            -M "{output.metrics}" \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            2> {log.dedup}
        echo "DEBUG: Finished dedup for {wildcards.sample}" >&2
        """

##############################################################################
# 10) BQSR
##############################################################################
def get_recal_table(sample):
    bqsr_dir = os.path.join(OUTPUT_FOLDER, "bqsr")
    os.makedirs(bqsr_dir, exist_ok=True)
    return os.path.join(bqsr_dir, f"{sample}.merged.dedup.recal_data.table")

def get_bqsr_bam(sample):
    bqsr_dir = os.path.join(OUTPUT_FOLDER, "bqsr")
    return os.path.join(bqsr_dir, f"{sample}.merged.dedup.bqsr.bam")

rule base_recalibration:
    input:
        dedup_bam = get_dedup_bam("{sample}")
    output:
        recal_table = get_recal_table("{sample}")
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 4400,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        recal = os.path.join(LOG_DIR, "recal.gatk.{sample}.log")
    shell:
        r"""
        echo "DEBUG: BaseRecalibrator for sample {wildcards.sample}" >&2
        gatk --java-options "-Xms4000m -Xmx7g -Djava.io.tmpdir={resources.tmpdir}" BaseRecalibrator \
            -I "{input.dedup_bam}" \
            -R "{REFERENCE_UNPACKED}" \
            --known-sites "{AF_ONLY_GNOMAD}" \
            --known-sites "{COMMON_BIALLELIC_GNOMAD}" \
            --known-sites "{PANEL_OF_NORMALS}" \
            -O "{output.recal_table}" \
            2> {log.recal}
        echo "DEBUG: Finished BaseRecalibrator for sample {wildcards.sample}" >&2
        """

rule apply_bqsr:
    input:
        dedup_bam    = get_dedup_bam("{sample}"),
        recal_table  = get_recal_table("{sample}")
    output:
        bqsr_bam     = get_bqsr_bam("{sample}")
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 4400,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        apply_bqsr = os.path.join(LOG_DIR, "apply_bqsr.gatk.{sample}.log")
    shell:
        r"""
        echo "DEBUG: ApplyBQSR for sample {wildcards.sample}" >&2
        gatk --java-options "-Xms4000m -Xmx7g -Djava.io.tmpdir={resources.tmpdir} -Dsamjdk.compression_level=6" \
            ApplyBQSR \
            -I "{input.dedup_bam}" \
            -bqsr "{input.recal_table}" \
            -O "{output.bqsr_bam}" \
            2> {log.apply_bqsr}
        echo "DEBUG: Finished ApplyBQSR for sample {wildcards.sample}" >&2
        """

##############################################################################
# 11) Final BQSR’d BAM => rule all
##############################################################################
# "rule all" is already defined above, collecting final BQSR bams.
##############################################################################
