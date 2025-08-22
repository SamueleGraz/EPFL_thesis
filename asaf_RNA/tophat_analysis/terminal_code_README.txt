Here the steps I followed in order to create ReadsPerGene tables after using tophat as aligner and 
featurecounts as counter

create conda env with tools available in 2014 since it was the data in which the version cufflinks 2.2 was released (for this reason we need python 2.7)
\\\
conda activate -n tophat_cufflinks python=2.7 -y
\\\

config additional channels from conda
\\\
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
\\\

install the packages
\\\
conda install bowtie2=2.2.5 tophat=2.0.13 cufflinks=2.2.1 -y
conda install -c bioconda subread 
\\\

create index with bowtie2 (like the one you create with star)
\\\
bowtie2-build /path/to/genome.fa /path/to/index/genome_index
\\\

run the loop for tophat alignment (see 01_tophat.bash file)
use this code in the terminal:
\\\
for fq in /scratch/sgrazian/tmp/marco_data/RNA_seq/fastq_data/*.fastq; do
    sbatch 01_tophat.bash "$fq"
done
\\\
this code use the variable fq and it put is as second argument after the command sbatch. You can set this 
second argument inside the bash script with $1 (0=1st argument after the command, the script itself, 
1=2nd argument). Basically the fq = /scratch/sgrazian/tmp/marco_data/RNA_seq/fastq_data/*.fastq will be
set with the variable 1 inside the bash script

create a table for the genecount using the package subread from conda (featurecounts tool)
see 02_genecount.bash file

clean the table counts.txt from the useless columns (chr, start, end, strand, length)
\\\
cut -f1,7- counts.txt > final_unstranded_counts.tsv
\\\

Then I deleted manually the first row and I changed the names of the columns so that
ncol and the nrow of marco_matrix and marco_samplesheet correspond

I obtained with deseq2 more genes significantly and differentially expressed with featurecounts. 
I run flagstag on the bam files originated after Tophat and I saw that the duplicates were not marked
next steps will be to mark those duplicates, sort and index the files and then to count them again and repeat the analysis

Here's the steps:

Sort by coordinate:
samtools sort -o input_sorted.bam input.bam

Now you can run markdup:
samtools markdup -s --duplicate-count --threads 8 --verbosity 2 input_sorted.bam output_dupmarked.bam

Then I want to check the number of duplicates
samtools flagstat RNA_Basal_rep1_tophat_dupmarked.bam > RNA_Basal_rep1_tophat_dupmarked_stats.txt
there are a lot of duplicates

And lastly I want to remove the duplicates in order to do again the R analysis
samtools markdup -s --duplicate-count --threads 8 --verbosity 2 -r \
RNA_Basal_rep1_tophat_sorted.bam RNA_Basal_rep1_tophat_dupremoved.bam

check again the duplicates with flagstat
samtools flagstat RNA_Basal_rep1_tophat_dupremoved.bam > RNA_Basal_rep1_tophat_dupremoved_stats.txt

SThen I checked with flagstat again and now I have 0 duplicates

After that, I run 02_genecount2.bash in order to count again the reads without dups
If I compare counts_nodups.txt with counts.txt I can already see differences