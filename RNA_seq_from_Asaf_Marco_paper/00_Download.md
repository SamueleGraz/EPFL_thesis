#Download fastq files using a computer cluster through SLURM

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
```
The following samples were downloaded:
+ RNA_Basal_rep1;
+ RNA_Basal_rep2;
+ RNA_Basal_rep3;
+ RNA_Early_rep1;
+ RNA_Early_rep2 ;
+ RNA_Early_rep3 ;
+ RNA_Late_rep1 ;
+ RNA_Late_rep2 ;
+ RNA_Late_rep3 ;
+ RNA_Reactivated_rep1 ;
+ RNA_Reactivated_rep2 ;
+ RNA_Reactivated_rep3 .
