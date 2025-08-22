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