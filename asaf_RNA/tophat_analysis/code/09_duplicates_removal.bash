#!/usr/bin/env bash

#SBATCH --job-name=dup_rem
#SBATCH --output=log_%x_%j.out
#SBATCH --error=log_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --qos=serial
#SBATCH --partition=bigmem
#SBATCH --mem=100GB

singularity exec ~/singularity_images/rna-seq-image_version_0.5.sif bash -c "
  source /opt/conda/etc/profile.d/conda.sh && \
  conda activate align && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Basal_rep1_tophat_sorted.bam RNA_Basal_rep1_tophat_dupremoved.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Basal_rep2_tophat_sorted.bam RNA_Basal_rep2_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Basal_rep3_tophat_sorted.bam RNA_Basal_rep3_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Early_rep1_tophat_sorted.bam RNA_Early_rep1_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Early_rep2_tophat_sorted.bam RNA_Early_rep2_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Early_rep3_tophat_sorted.bam RNA_Early_rep3_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Late_rep1_tophat_sorted.bam RNA_Late_rep1_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Late_rep2_tophat_sorted.bam RNA_Late_rep2_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Late_rep3_tophat_sorted.bam RNA_Late_rep3_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Reactivated_rep1_tophat_sorted.bam RNA_Reactivated_rep1_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Reactivated_rep2_tophat_sorted.bam RNA_Reactivated_rep2_tophat_removed.bam && \
  samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r RNA_Reactivated_rep3_tophat_sorted.bam RNA_Reactivated_rep3_tophat_removed.bam"