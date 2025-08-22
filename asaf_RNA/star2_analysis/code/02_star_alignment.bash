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

cd /scratch/sgrazian/tmp/marco_data/RNA_seq

echo "Running STAR for file: ${FILE}"

singularity exec ~/singularity_images/rna-seq-image_version_0.5.sif bash -c "
  source /opt/conda/etc/profile.d/conda.sh && \
  conda activate align && \
  STAR \
    --genomeDir star_analysis/star2/star_codes/mouse_star_index \
    --readFilesIn fastq_data/${FILE}.fastq \
    --runMode alignReads \
    --runThreadN 8 \
    --quantMode GeneCounts \
    --outFileNamePrefix Aligned_${FILE}_ \
    --outSAMtype BAM SortedByCoordinate \
    --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
    --outWigType bedGraph \
    --outWigNorm RPM \
    --outWigStrand Unstranded"