#!/usr/bin/env bash

#SBATCH --job-name=star_trial
#SBATCH --output=log_%x_%j.out
#SBATCH --error=log_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --qos=serial
#SBATCH --partition=bigmem
#SBATCH --mem=100GB

# === Define variables ===
FASTQ_FILE=$1  # pass the file as argument
FASTQ_DIR="/scratch/sgrazian/tmp/marco_data/RNA_seq/fastq_data"
GTF_FILE="/scratch/sgrazian/tmp/marco_data/reference/gencode_mouse_annotated_genome.gtf"
BOWTIE2_INDEX="/scratch/sgrazian/tmp/marco_data/RNA_seq/cufflinks_analysis/bowtie2_index/gencode_mouse_index"
OUT_DIR="/scratch/sgrazian/tmp/marco_data/RNA_seq/cufflinks_analysis"

# === Run Tophat ===
BASE=$(basename "$FASTQ_FILE" .fastq)
SAMPLE_OUT="$OUT_DIR/${BASE}_tophat_out"
mkdir -p "$SAMPLE_OUT"

echo "Running TopHat2 on $FASTQ_FILE..."
tophat2 -p 8 -G "$GTF_FILE" -o "$SAMPLE_OUT" --no-coverage-search "$BOWTIE2_INDEX" "$FASTQ_FILE" > "$SAMPLE_OUT/tophat.log" 2>&1
echo "Finished $BASE"