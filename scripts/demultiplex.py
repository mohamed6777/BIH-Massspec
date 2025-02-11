import os
import sys
import gzip
import pandas as pd

index_file = sys.argv[1]
r1_file = sys.argv[2]
r2_file = sys.argv[3]
output_dir = sys.argv[4]

# Load index file containing barcodes and sample assignments
index_df = pd.read_csv(index_file, sep="\t", header=None, names=["barcode", "sample"])

sample_dict = {row["barcode"]: row["sample"] for _, row in index_df.iterrows()}

# Function to write reads into respective files
def demultiplex(fq1, fq2, sample_dict, output_dir):
    handles = {}
    for sample in set(sample_dict.values()):
        handles[sample] = (
            gzip.open(f"{output_dir}/{sample}_R1.fastq.gz", "wt"),
            gzip.open(f"{output_dir}/{sample}_R2.fastq.gz", "wt")
        )

    with gzip.open(fq1, "rt") as r1, gzip.open(fq2, "rt") as r2:
        while True:
            try:
                r1_header, r1_seq, _, r1_qual = [next(r1).strip() for _ in range(4)]
                r2_header, r2_seq, _, r2_qual = [next(r2).strip() for _ in range(4)]

                barcode = r1_seq[:8]  # Assuming first 8 bases contain the barcode
                if barcode in sample_dict:
                    sample = sample_dict[barcode]
                    handles[sample][0].write(f"{r1_header}\n{r1_seq}\n+\n{r1_qual}\n")
                    handles[sample][1].write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")
            except StopIteration:
                break

    for handle in handles.values():
        handle[0].close()
        handle[1].close()

demultiplex(r1_file, r2_file, sample_dict, output_dir)
