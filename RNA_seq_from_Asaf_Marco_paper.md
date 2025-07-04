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

04.07.25
I started to work on scratch because of disk memory problems.
I solved some issues with star and the slurm file to send multiple jobs. Here how I created the star_index:

```

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
```
In that step, I had problems on make the slurm file read the sequences, so I had to put all the files in the same folders because it seems like star command works just if the files and the directories you are calling are inside the same directory in which you run the sbatch or in a sub-folder.
To avoid this problem and keep separated the sbatch files and the outputs in different folders, I added one code in the next slurm file:
```
cd /scratch/sgrazian/tmp/marco_data/RNA_seq
```
In this way slurm will start to run STAR just from that common folder. Moreover, in order to run in paralle different jobs I run a loop code in order to re-run everytime the same slurm file with a different file in the parameter that I specified. Here the loop I wrote and the slurm file.
slurm file (trial_loop.bash):
```
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
```
text file (sequence.txt) with the paramters used from the loop:
```
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
```
Here the while loop:
```
while IFS= read -r i; do echo "Submitting job for file: $i"; sbatch --export=FILE="$i" trial_loop.bash; done < sequence.txt
```








