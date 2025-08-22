#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --job-name=star_trial
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --qos=serial
#SBATCH --partition=bigmem
#SBATCH --mem=100GB

singularity exec ~/singularity_images/rna-seq-image_version_0.5.sif bash -c "
  source /opt/conda/etc/profile.d/conda.sh && \
  conda activate align && \
  STAR \
  --runMode genomeGenerate \
  --genomeDir mouse_star_index \
  --genomeFastaFiles /home/sgrazian/work/reference/gencode_mouse_genome.fa \
  --sjdbGTFfile /home/sgrazian/work/reference/gencode.vM33.primary_assembly.annotation.gtf \
  --sjdbOverhang 50 \
  --outFileNamePrefix mouse 
"