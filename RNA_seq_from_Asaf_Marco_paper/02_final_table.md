#Readspergene filtration to have just the Gene ID and the unstranded reads.

#!/usr/bin/env bash

#SBATCH --job-name=filter_unstranded
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --qos=serial
#SBATCH --partition=bigmem
#SBATCH --mem=100GB

# Variable passed externally
FILE="$TABLE"

# Check that the variable has been passed
if [ -z "$FILE" ]; then
    echo "Error: No file specified. Use --export=TABLE=nomefile"
    exit 1
fi

# Extract GeneID and Unstranded
echo "Filtering the file: $FILE"
cut -f1,2 "$FILE" > "${FILE%.out.tab}_unstranded.txt"
echo "Done: created ${FILE%.out.tab}_unstranded.txt"

#Here's the loop in order to create different slurm files for each ReadsPerGene.out.tab file:

for f in *ReadsPerGene.out.tab; do sbatch --export=TABLE="$f" 03_table_filtering.bash; done

#All the files were moved in a new folder where they were sorted by the ID so as to merge all the files in one unique table

mkdir gene_count_unstranded

#move unstranded files
mv *unstranded* gene_count_unstranded/

#sort files by ID
for f in *_unstranded.txt; do     sort -k1,1 "$f" > "${f%.txt}_sorted.txt"; done

#create a header.txt file where to keep the names of the files. It will be useful in order to create the header for the final table:
echo -e "GeneID\t$(ls *_unstranded_sorted.txt | sed 's/_ReadsPerGene_unstranded_sorted.txt//g' | paste -sd '\t')" > header.txt

#A merged.txt was made where to merge all the files
for f in *_unstranded_sorted.txt; do     if [[ "$f" != "merged.txt" ]]; then         join -a1 -a2 -e 0 -o auto -t $'\t' merged.txt "$f" > tmp && mv tmp merged.txt;     fi; done

#The header was made
cat header.txt merged.txt > final_unstranded_counts.tsv
