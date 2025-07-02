In the next lines, you will see my analysis done in order to compare the available RNA-seq data from [**Asaf Marco 's paper**]([URL]https://www.nature.com/articles/s41593-020-00717-0#data-availability) with the data generated from Johannes Graeff lab.

First of all, I created a Dockerfile in order to use specific tools and to be reproducible. I used [**samuelegraz/rna-seq-image:version_0.3**]([URL]https://hub.docker.com/repository/docker/samuelegraz/rna-seq-image/tags/version_0.3/sha256-df5af533563bc42b992b953af69a3041f3d643828bdc1484d5e360ddb0284ccc) image. 

Then I used mostly slurm sbatch files in order to execute multiple tasks in the same time.
Here the code I used for downloading the fastq files. Note that the files from Marco, were already trimmed, so it was possible to start from the alignment using STAR.
```
#Ex. download fastq files

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
I downloaded the following samples:
+ RNA_Early_rep1 ;
+ RNA_Early_rep2 ;
+ RNA_Early_rep3 ;
+ RNA_Late_rep1 ;
+ RNA_Late_rep2 ;
+ RNA_Late_rep3 ;
+ RNA_Reactivated_rep1 ;
+ RNA_Reactivated_rep2 ;
+ RNA_Reactivated_rep3 .







