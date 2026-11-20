#!/bin/bash

#SBATCH --job-name=dorado
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=3
#SBATCH --mem-per-gpu=20G
#SBATCH --cpus-per-gpu=4
#SBATCH --time=14:00:00
#SBATCH --export=ALL

absolute_path="/home/exacloud/gscratch/BakerLab"
pod5_path="pod5_skip"
working_dir=" "

mkdir ${absolute_path}/${working_dir}/dorado-output/

srun ${absolute_path}/shared_resources/dorado-0.9.1-linux-x64/bin/dorado basecaller ${absolute_path}/shared_resources/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
   ${absolute_path}/${pod5_path} \
   --batchsize 450 \
   --min-qscore 9 \
   --no-trim \
   > ${absolute_path}/${working_dir}/dorado-output/supv500.bam

srun ${absolute_path}/shared_resources/dorado-0.9.1-linux-x64/bin/dorado demux \
   ${absolute_path}/${working_dir}/dorado-output/supv500.bam \
   --kit-name SQK-NBD114-24 \
   --emit-fastq \
   --output-dir ${absolute_path}/${working_dir}/dorado-output/demux

for f in "${absolute_path}/${working_dir}/dorado-output/demux/"*_barcode*.fastq "${absolute_path}/${working_dir}/dorado-output/demux/"*_unclassified.fastq; do
  newname=$(echo "$f" | sed -E 's/.*_(barcode[0-9]+|unclassified)\.fastq/\1.fastq/')
  mv "$f" "${absolute_path}/${working_dir}/dorado-output/demux/${newname}"
done
