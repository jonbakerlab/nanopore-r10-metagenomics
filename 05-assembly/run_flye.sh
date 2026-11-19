#!/bin/bash
#SBATCH --job-name=flye_ens
#SBATCH --partition=interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --output=flye_ens.out
#SBATCH --error=flye_ens.err

module load flye

# Execute Flye
flye --nano-raw filtered_metagenomics_reads.fastq --out-dir metagenomic_assembly --meta --threads 14
