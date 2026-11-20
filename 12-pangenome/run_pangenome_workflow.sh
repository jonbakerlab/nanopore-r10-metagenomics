#!/bin/bash
# Run the anvi'o pangenomics workflow on SLURM. max_threads in pan-configs.json should match --cores.

anvi-run-workflow -w pangenomics -c pan-configs.json \
  --additional-params --jobs 24 --keep-going --rerun-incomplete \
  --cores 24 \
  --cluster 'sbatch --job-name=Rothia_pangenome --output={log} --error={log} --nodes=1 --ntasks=1 --cpus-per-task=1'
