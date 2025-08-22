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

bowtie2-build /scratch/sgrazian/tmp/marco_data/reference/gencode_mouse_genome.fa /scratch/sgrazian/tmp/marco_data/reference/gencode_mouse_index