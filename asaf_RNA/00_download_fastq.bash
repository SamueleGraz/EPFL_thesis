#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --job-name=download
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=100GB

singularity exec ~/singularity_images/rna-seq-image_version_0.3.sif bash -c "
  source /opt/conda/etc/profile.d/conda.sh && \
  conda activate sra && \
  fasterq-dump SRR12045661 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Early_rep1.fastq \
  --bufsize 5
  \
    fasterq-dump SRR12045662 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Early_rep2.fastq \
  --bufsize 5
  \
   fasterq-dump SRR12045663 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Early_rep3.fastq \
  --bufsize 5
  \
    fasterq-dump SRR12045664 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Late_rep1.fastq \
  --bufsize 5
  \
    fasterq-dump SRR12045665 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Late_rep2.fastq \
  --bufsize 5
  \
   fasterq-dump SRR12045666 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Late_rep3.fastq \
  --bufsize 5
  \
    fasterq-dump SRR12045667 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Reactivated_rep1.fastq \
  --bufsize 5
  \
   fasterq-dump SRR12045668 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Reactivated_rep2.fastq \
  --bufsize 5
  \
    fasterq-dump SRR12045669 \
  --outdir /home/sgrazian/tmp/marco_data/RNA_seq/fastq_data \
  --progress \
  --verbose \
  --outfile RNA_Reactivated_rep3.fastq \
  --bufsize 5
"