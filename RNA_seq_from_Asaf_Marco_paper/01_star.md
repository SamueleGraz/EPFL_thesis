#creating the index used by star tool through calling the image made with the dockerfile

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

singularity exec ~/singularity_images/rna-seq-image_version_0.3.sif bash -c "
  source /opt/conda/etc/profile.d/conda.sh && \
  conda activate align && \
  STAR \
  --runMode genomeGenerate \
  --genomeDir mouse_star_index \
  --genomeFastaFiles /scratch/sgrazian/tmp/marco_data/RNA_seq/reference/gencode_mouse_genome.fa \
  --sjdbGTFfile scratch/sgrazian/tmp/marco_data/RNA_seq/reference/gencode_mouse_annotated_genome.gtf \
  --sjdbOverhang 50 \
  --outFileNamePrefix mouse 
"


#Run star command in order to have the gene counts. A loop was used in order to create a different job for each file was needed to be analysed

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

singularity exec ~/singularity_images/rna-seq-image_version_0.3.sif bash -c "
  source /opt/conda/etc/profile.d/conda.sh && \
  conda activate align && \
  STAR \
    --genomeDir star/mouse_star_index \
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

#text file (sequence.txt) with the parameters used from the loop:

RNA_Basal_rep1
RNA_Basal_rep2
RNA_Basal_rep3
RNA_Early_rep1
RNA_Early_rep2
RNA_Early_rep3
RNA_Late_rep1
RNA_Late_rep2
RNA_Late_rep3
RNA_Reactivated_rep1
RNA_Reactivated_rep2
RNA_Reactivated_rep3


#Here's the while loop used in the terminal:

while IFS= read -r i; do echo "Submitting job for file: $i"; sbatch --export=FILE="$i" trial_loop.bash; done < sequence.txt
