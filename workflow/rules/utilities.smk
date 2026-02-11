import glob as _glob


# =============================================================================
# Subset BAM by BED regions
# =============================================================================
SUBSET_CFG = config.get("subset", {})
SUBSET_BED = SUBSET_CFG.get("bed_file", "")
SUBSET_SUFFIX = SUBSET_CFG.get("output_suffix", ".subset.bam")
SUBSET_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "subset")


rule subset_bam:
    """Extract BAM reads overlapping regions in a BED file."""
    input:
        bam_file=os.path.join(BQSR_DIR, "{sample}" + FINAL_BAM_SUFFIX),
        bed_file=SUBSET_BED,
    output:
        bam=os.path.join(SUBSET_OUTPUT_DIR, "{sample}" + SUBSET_SUFFIX),
    threads: 1
    resources:
        runtime=120,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        os.path.join(LOG_DIR, "samtools_view_subset.{sample}.log"),
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})"

        samtools view -b -h -P -M \
            -L {input.bed_file} \
            {input.bam_file} \
            > {output.bam} \
            2> {log}
        """


# =============================================================================
# BAM to FASTQ conversion
# =============================================================================
rule convert_bam_to_fastq:
    """Convert BAM back to paired-end FASTQ files."""
    input:
        bam_file="{prefix}.bam",
    output:
        r1="{prefix}_R1.fastq.gz",
        r2="{prefix}_R2.fastq.gz",
    threads: 8
    resources:
        runtime=1440,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        "{prefix}.bam_to_fastq.log",
    shell:
        r"""
        samtools sort -n -@ {threads} {input.bam_file} \
        | samtools fastq -@ {threads} - \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            2> {log}
        """


# =============================================================================
# MD5 checksums
# =============================================================================
def _get_checksum_files(folder):
    """Walk a directory tree and return all .fastq.gz relative paths."""
    result = []
    for root, _, files in os.walk(folder):
        for f in files:
            if f.endswith(".fastq.gz"):
                result.append(os.path.relpath(os.path.join(root, f), folder))
    return result


rule calculate_md5sum:
    """Calculate MD5 checksum for a single file."""
    input:
        fastq=os.path.join(FASTQ_DIR, "{file}"),
    output:
        md5sum=os.path.join(OUTPUT_DIR, "md5sum", "{file}.md5sum"),
    shell:
        r"""
        mkdir -p "$(dirname {output.md5sum})"
        md5sum {input.fastq} > {output.md5sum}
        """


rule join_md5sums:
    """Concatenate all individual checksums into a single file."""
    input:
        lambda wc: expand(
            os.path.join(OUTPUT_DIR, "md5sum", "{file}.md5sum"),
            file=_get_checksum_files(FASTQ_DIR),
        ),
    output:
        os.path.join(OUTPUT_DIR, "all_md5sums.txt"),
    shell:
        "cat {input} > {output}"


# =============================================================================
# Simplified align + sort (for pre-subset FASTQ)
# =============================================================================
rule align_and_sort:
    """BWA alignment with samblaster dedup piped to samtools sort."""
    input:
        r1="{prefix}.bbduk_R1.fastq.gz",
        r2="{prefix}.bbduk_R2.fastq.gz",
    output:
        sorted_bam="{prefix}.sorted.bam",
    params:
        reference=REF_GZ,
        read_group=lambda wc: (
            f"@RG\\tID:{os.path.basename(wc.prefix)}"
            f"\\tSM:{os.path.basename(wc.prefix)}"
            f"\\tLB:{os.path.basename(wc.prefix)}"
            f"\\tPL:{PLATFORM}"
            f"\\tPU:{os.path.basename(wc.prefix)}"
        ),
    threads: 8
    resources:
        runtime=1440,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        bwa="{prefix}.bwa.log",
        samtools="{prefix}.samtools_sort.log",
    shell:
        r"""
        bwa mem -C -t {threads} \
            -R "{params.read_group}" \
            {params.reference} \
            {input.r1} {input.r2} \
            2> {log.bwa} \
        | samblaster 2> /dev/null \
        | samtools view -@ {threads} -bS -o /dev/stdout /dev/stdin \
        | samtools sort -@ {threads} -o {output.sorted_bam} \
            2> {log.samtools}
        """


rule index_bam:
    """Index a sorted BAM file."""
    input:
        sorted_bam="{prefix}.sorted.bam",
    output:
        index="{prefix}.sorted.bam.bai",
    threads: 1
    resources:
        runtime=120,
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        "{prefix}.samtools_index.log",
    shell:
        r"""
        samtools index -@ {threads} {input.sorted_bam} 2> {log}
        """
