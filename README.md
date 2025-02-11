# BIH-Massspec Pipeline

This Snakemake pipeline is designed to process sequencing data through various stages, including quality control, merging, demultiplexing, UMI deduplication, and alignment to a mixed human/mouse genome.

## Pipeline Overview

1. **Initial QC of Raw FASTQ Files:**
   - Runs FastQC on raw sequencing files to generate quality control reports.

2. **Merging FASTQ Files:**
   - Merges R1 and R2 FASTQ files into combined files for further processing.

3. **QC of Merged FASTQ Files:**
   - Runs FastQC on merged FASTQ files to generate quality control reports.

4. **Demultiplexing Merged FASTQ Files:**
   - Demultiplexes merged FASTQ files into individual sample files using an index file.

5. **QC of Demultiplexed FASTQ Files:**
   - Runs FastQC on demultiplexed FASTQ files to generate quality control reports.

6. **UMI Deduplication:**
   - Removes duplicate reads using Unique Molecular Identifiers (UMIs).

7. **Alignment to Mixed Human/Mouse Genome:**
   - Aligns reads to a mixed human and mouse genome using the STAR aligner.

## Configuration

The pipeline uses a configuration file (`config.yaml`) to specify input files, output directory, genomes, and other parameters.

## Usage

1. **Install Dependencies:**
   - Ensure you have Snakemake and the necessary conda environment set up.

2. **Run the Pipeline:**
   - Execute the Snakemake pipeline with the following command:
     ```
     snakemake --use-conda
     ```

## Key Files

- `Snakefile`: Defines the pipeline rules and workflow.
- `config.yaml`: Configuration file with input and output specifications.
- `environment.yaml`: Conda environment file for dependency management.
- `scripts/demultiplex.py`: Script for demultiplexing FASTQ files.
- `scripts/umi_dedup.py`: Script for UMI deduplication.

## Output

The output directory will contain the following subdirectories and files:
- `qc/`: Quality control reports for raw, merged, and demultiplexed FASTQ files.
- `merged_R1.fastq.gz`, `merged_R2.fastq.gz`: Merged FASTQ files.
- `demultiplexed/`: Demultiplexed FASTQ files for each sample.
- `deduplicated/`: UMI deduplicated FASTQ files for each sample.
- `aligned/`: BAM files of aligned reads for each sample.

## Notes

- Ensure the paths and filenames in the configuration file are correctly specified before running the pipeline.
- Adjust the number of threads and other parameters as needed in the `config.yaml` file.
