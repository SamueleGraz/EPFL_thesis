#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --job-name=qualimap_trial
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --qos=serial
#SBATCH --partition=bigmem
#SBATCH --mem=100GB

#for bamfile in *.bam; do
 # sample_name=$(basename "$bamfile" .bam)
  #echo "Running Qualimap for $sample_name"
  
singularity exec ~/singularity_images/rna-seq-image_version_0.5.sif bash -c "
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate align && \
    qualimap rnaseq -bam Aligned_RNA_Reactivated_rep3_Aligned.sortedByCoord.out.bam -gtf /home/sgrazian/work/reference/gencode.vM33.primary_assembly.annotation.gtf --java-mem-size=4G -outdir qualimap_results/Aligned_RNA_Reactivated_rep3_Aligned.sortedByCoord.out -outformat PDF:HTML"
done