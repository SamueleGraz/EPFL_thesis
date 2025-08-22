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

featureCounts -a /scratch/sgrazian/tmp/marco_data/reference/gencode_mouse_annotated_genome.gtf \
 -o counts.txt \
 RNA_Basal_rep1_tophat_out/accepted_hits.bam \
 RNA_Basal_rep2_tophat_out/accepted_hits.bam \
 RNA_Basal_rep3_tophat_out/accepted_hits.bam \
 RNA_Early_rep1_tophat_out/accepted_hits.bam \
 RNA_Early_rep2_tophat_out/accepted_hits.bam \
 RNA_Early_rep3_tophat_out/accepted_hits.bam \
 RNA_Late_rep1_tophat_out/accepted_hits.bam \
 RNA_Late_rep2_tophat_out/accepted_hits.bam \
 RNA_Late_rep3_tophat_out/accepted_hits.bam \
 RNA_Reactivated_rep1_tophat_out/accepted_hits.bam \
 RNA_Reactivated_rep2_tophat_out/accepted_hits.bam \
 RNA_Reactivated_rep3_tophat_out/accepted_hits.bam \
 -T 8 \
 --verbose