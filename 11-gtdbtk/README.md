# 11. Classify your genome of interest using GTDB-Tk

Note: make sure to have the latest version, which uses *skani* to calculate Average Nucleotide Identity (ANI) and incorporates it into `classify_wf`.

```sh
gtdbtk --version
conda update -c bioconda gtdbtk
```

## Download and alias the GTDB-Tk reference data (R226)

The most recent database is R226 at the moment, but this may vary; check the GTDB-Tk website for updates.

```sh
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```
```sh
export GTDBTK_DATA_PATH="/home/exacloud/gscratch/BakerLab/charlom/gtdbtk_db/release226"
echo 'export GTDBTK_DATA_PATH="/home/exacloud/gscratch/BakerLab/charlom/gtdbtk_db/release226"' >> ~/.bashrc
source ~/.bashrc
```
```sh
gtdbtk check_install
```
Note: deactivate and reactivate the conda environment.

## Make output directories

```sh
mkdir output_gtdb-tk
mkdir scratch_dir
```

Note: set up an interactive job if needed:
```sh
srun --job-name=gtdbtk --nodes=1 -c 12 --mem-per-cpu=4G -t 03:00:00 --partition=interactive --ntasks-per-node=1 --pty --export=ALL /usr/bin/bash
```

## Run GTDB-Tk classify_wf

Add the path to your closed genome(s) of interest after `--genome_dir`.

```sh
gtdbtk classify_wf \
  --genome_dir \
  --out_dir output_gtdb-tk \
  --mash_db $GTDBTK_DATA_PATH/gtdbtk.mashdb \
  -x fasta \
  --cpus 12 \
  --pplacer_cpus 1 \
  --scratch_dir scratch_dir
```

OR

```sh
gtdbtk classify_wf \
  --genome_dir \
  --out_dir output_gtdb-tk \
  --skip_ani_screen \
  -x fasta \
  --cpus 12 \
  --pplacer_cpus 1 \
  --scratch_dir scratch_dir
```

The first command runs the GTDB-Tk classify workflow (`classify_wf`) using a MASH database (`--mash_db`) to quickly pre-filter the input genomes, identifying bacterial and archaeal genomes for classification. The second command also runs `classify_wf` but skips the MASH/ANI filtering step (`--skip_ani_screen`), meaning it classifies all genomes provided regardless of their content, which is useful if you are confident that all input genomes are valid.

```sh
gtdbtk ani_rep --genome_dir --out_dir -x fasta --cpus 12
```

`ani_rep` only calculates the ANI between your genomes and reference genomes (R226 GTDB-Tk database) without performing taxonomic classification. This is useful if you only need similarity information.

The main output is a `summary.tsv` file.
