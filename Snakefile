# Snakefile

import os

configfile: "config.yaml"

# Extract variables from the config file
fastq_files = config["fastq_files"]
index_file = config["index_file"]
output_dir = config["output_dir"]
human_genome = config["human_genome"]
mouse_genome = config["mouse_genome"]
threads = config["threads"]

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Rule to define the conda environment
rule all:
    input:
        expand(f"{output_dir}/qc/{{sample}}_fastqc.html", sample=["merged_R1", "merged_R2"]),
        expand(f"{output_dir}/demultiplexed/{{sample}}_R1.fastq.gz", sample=["sample1", "sample2", "sample3"]),
        expand(f"{output_dir}/aligned/{{sample}}.bam", sample=["sample1", "sample2", "sample3"])

# Rule for initial QC of raw FASTQ files
rule fastqc_raw:
    input:
        fastq=fastq_files
    output:
        html=expand(f"{output_dir}/qc/{{sample}}_fastqc.html", sample=["raw_" + os.path.basename(f).split(".")[0] for f in fastq_files])
    conda:
        "environment.yaml"
    shell:
        """
        fastqc -o {output_dir}/qc {input.fastq}
        """

# Rule to merge FASTQ files into R1 and R2
rule merge_fastq:
    input:
        R1=[f for f in fastq_files if "_R1" in f],
        R2=[f for f in fastq_files if "_R2" in f]
    output:
        R1=f"{output_dir}/merged_R1.fastq.gz",
        R2=f"{output_dir}/merged_R2.fastq.gz"
    conda:
        "environment.yaml"
    shell:
        """
        cat {input.R1} | gzip > {output.R1}
        cat {input.R2} | gzip > {output.R2}
        """

# Rule for QC of merged FASTQ files
rule fastqc_merged:
    input:
        R1=f"{output_dir}/merged_R1.fastq.gz",
        R2=f"{output_dir}/merged_R2.fastq.gz"
    output:
        R1_html=f"{output_dir}/qc/merged_R1_fastqc.html",
        R2_html=f"{output_dir}/qc/merged_R2_fastqc.html"
    conda:
        "environment.yaml"
    shell:
        """
        fastqc -o {output_dir}/qc {input.R1} {input.R2}
        """

# Rule for demultiplexing merged FASTQ files
rule demultiplex:
    input:
        R1=f"{output_dir}/merged_R1.fastq.gz",
        R2=f"{output_dir}/merged_R2.fastq.gz",
        index=index_file
    output:
        R1=expand(f"{output_dir}/demultiplexed/{{sample}}_R1.fastq.gz", sample=["sample1", "sample2", "sample3"]),
        R2=expand(f"{output_dir}/demultiplexed/{{sample}}_R2.fastq.gz", sample=["sample1", "sample2", "sample3"])
    conda:
        "environment.yaml"
    script:
        "scripts/demultiplex.py"

# Rule for QC of demultiplexed FASTQ files
rule fastqc_demultiplexed:
    input:
        R1=expand(f"{output_dir}/demultiplexed/{{sample}}_R1.fastq.gz", sample=["sample1", "sample2", "sample3"]),
        R2=expand(f"{output_dir}/demultiplexed/{{sample}}_R2.fastq.gz", sample=["sample1", "sample2", "sample3"])
    output:
        html=expand(f"{output_dir}/qc/{{sample}}_fastqc.html", sample=["demux_" + s for s in ["sample1", "sample2", "sample3"]])
    conda:
        "environment.yaml"
    shell:
        """
        fastqc -o {output_dir}/qc {input.R1} {input.R2}
        """

# Rule for UMI deduplication
rule umi_deduplication:
    input:
        R1=lambda wildcards: f"{output_dir}/demultiplexed/{wildcards.sample}_R1.fastq.gz",
        R2=lambda wildcards: f"{output_dir}/demultiplexed/{wildcards.sample}_R2.fastq.gz"
    output:
        R1=f"{output_dir}/deduplicated/{{sample}}_R1.fastq.gz",
        R2=f"{output_dir}/deduplicated/{{sample}}_R2.fastq.gz"
    conda:
        "environment.yaml"
    script:
        "scripts/umi_dedup.py"

# Rule for alignment to mixed human/mouse genome
rule align_reads:
    input:
        R1=f"{output_dir}/deduplicated/{{sample}}_R1.fastq.gz",
        R2=f"{output_dir}/deduplicated/{{sample}}_R2.fastq.gz"
    output:
        bam=f"{output_dir}/aligned/{{sample}}.bam"
    conda:
        "environment.yaml"
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {human_genome},{mouse_genome} \
            --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {output.bam}
        """

