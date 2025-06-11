# Complete Metagenome Assembled Genomes (MAGs) pipeline

## Introduction

High-molecular-weight genomic DNA was extracted from saliva and metagenomic long-read library was prepared using a ligation sequencing gDNA - Native Barcoding Kit 24 V14 (Oxford Nanopore Technologies) and sequenced on a MinION using an R10 flow cell (Oxford Nanopore Technologies) using the adaptive sampling option.

This analysis pipeline can take the pod5 files outputted by a minION Nanopore sequencer and do the super high accuracy basecalling using Dorado v5.0.0. with the aim to obtain complete, high quality microbial assembled genomes from saliva.
This work was done at the BakerLab with support from Dr. Jonathon Baker and his team, specially Br. Matthew Barbisan.

### A note on adaptive sampling:

If you are interested in obtaining a higher ratio of microbial DNA reads / Human DNA reads adaptive sampling (https://nanoporetech.com/es/document/adaptive-sampling) can help a lot. First you would have to download the human reference genome.
```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```
To set up your run selecting the "Deplete Sequences" options in the Adaptive sampling section as well as  attaching the human reference genome file as the aligment reference. This file should be in the same minKNOW directory for it to be accesible.

## 1 - Super high accuracy (sup) basecalling using dorado

You can do this either locally using MinKNOW or using dorado on the ARC. The last option is going to be faster.

### 1.1 - Installing dorado

Dorado can be downloaded [here](https://github.com/nanoporetech/dorado) onto a local device and copied via `scp` to its desired directory, for more information on installing dorado see the Bacterial-Methylation-Workflow (https://github.com/jonbakerlab/Baker-Lab-Scripts/commit/ee94e5f2da2773b992c020b29905b9baf3f4c267).

Alongside `dorado`, the sup model also needs to be downloaded.

```sh
./dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```

### 1.2 - Running the sup model with dorado on the ARC

With all of the .pod5 files, dorado can be used to basecall using the sup model by inputting the appropriate directory. Adapter sequences are trimmed by default once the demultimplexing step is completed. The following script will output, quality filtered (>Q9), demultiplexed .fastq files.

```sh
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
```

### Notes:
If using minKNOW, make sure to select the apropiate basecalling model as well as to check the options to do the quality filtering (e.g: >Q9), demultiplexing and trimming barcodes. The model of basecalling can be found in the first line of the .fastq files (you can check this using the command head), in basecall_model_version_id where the acronym sup stands for 'super high accuracy'. Other accronyms include: * fast == fast basecalling model * hac == high accuracy. 

## 2 - Install and set up conda environments

Creating a conda enviroments to compartimentalize all needed packages for our analysis is recommended. For more information on installing conda and mamba see the "Introduction to Conda and Mamba" repository (https://github.com/jonbakerlab/conda-microbiome-tutorial). Here is an example for creating a conda environment and installing minimap2. This would need to be repeated for each pogram needed.

### 2.1 - Create a conda enviornment

```sh
conda create -n minimap2 -c bioconda minimap2
conda activate minimap2
conda install bioconda::minimap2
```
### Programs needed:

minimap2, samtools, parallel, seqkit, Flye, Sylph, dnaapler, bakta, Anvio (developers version).

## 3 - Filter out human reads using *minimap2*

Since we are interested in microbial genomic DNA we will filter out all human reads.

### 3.1 - Download the human reference genome

```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

### 3.1 - Index the human genome

```sh
minimap2 -d human.mmi GRCh38.primary_assembly.genome.fa
```

### 3.3 - Align the reads to the human genome

Make sure to previously unzip you fastq.gz files.

```sh
minimap2 -ax lr:hq human.mmi barcode01/*.fastq > aligned.sam
```

### 3.4 - Filter out human reads

This command selects unmapped reads (reads that did not align with the human genome). The output is a bam file that needs to be converted back to fastq.

```sh
samtools view -@ 14 -b -F 4 aligned.sam > human_reads.bam
samtools view -@ 14 -b -f 4 aligned.sam > non_human_reads.bam
```
```sh
samtools fastq -@ 14 human_reads.bam > huaman_reads.fastq
samtools fastq -@ 14 non_human_reads.bam > filtered_metagenomics_reads.fastq
```

### Notes:

It is worth to keep the human_reads.fastq file to compare file sizes with the non_human_reads.fastq file, as well as to check that the sum of the size of both files makes sense with the unfiltered file you started with.

You can view file sizes using:
```sh
ls -lh 
```
or 
```sh
du -sh *yourfile*
```

Concatenate two or more metagenomics_filtered.fastq files from different runs if needed.

```sh
cat filtered_metagenomics_reads_1.fastq filtered_metagenomics_reads_2.fastq > filtered_metagenomics_reads_merged.fastq
```
Verify the integrity of the merged file

```sh
parallel -j 14 grep -c "^@" ::: filtered_metagenomics_reads_1.fastq filtered_metagenomics_reads_2.fastq
```
```sh
parallel -j 14 grep -c "^@" ::: filtered_metagenomics_reads_merged.fastq
```

## 4 - Genomic assembly using Flye

To assemble the long reads I used Flye with the --meta flag, making sure to use enough threads to speed up the process. The input file for Flye is going to be the filtered file that was the output from minimap2 which contains only microbial DNA reads.

```sh
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

#Excute Flye
flye --nano-raw filtered_metagenomics_reads.fastq --out-dir metagenomic_assembly --meta --threads 14
```
```sh
sbatch --nodes=1 --ntasks=1 --cpus-per-task=14 --mem=100G --time=12:00:00 --partition=interactive flye_metagenome.sh
```

### Visualize the .gfa file from the assembly using Bandage and easily identify circularized contigs.

## 5 - Profiling with sylph

For taxonomic assignment I used Sylph. Sylph can profile can calculate the abundances of genomes in a metagenomic sample by using a reference database (GTDB-R220) using nanopore reads with high precision.

```sh
conda install -c bioconda sylph
```
### 5.1 - Download GTDB-R220 pre-built database (~13 GB)
```sh
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb
```
### Note:
 There is no pre-built version of the latest GTDB-R226 database yet. .
 
 *.syldb* files (generated by sketching genomes) are aggregated together into one large database.

### 5.2 - Sketch

Compact representations of metagenomic reads (in .sylsp format) are generated by extracting k-mers. This representation is more manageable for analysis and comparison with reference databases (GTDB-R220 pre-built).

```sh
sylph sketch -t 14 filtered_metagenomics_reads.fastq -o metagenomic_sketch -d output_directory
```
The output is a *.sylsp* file which is located in the output_directory.

### 5.3 - Standard profiling

```sh
sylph profile gtdb-r220-c200-dbv1.syldb output_directory/filtered_metagenomics_reads.fastq.sylsp -t 14 > profiling_results.tsv
```
*.sylsp* file (generated by sketching reads) is queried/profiled against the combined database *.syldb*. Can also redirect using -o instead of >.

### 5.4 - Incorporating taxonomy into sylph

Sylph's profiling_results .tsv outputs do not have taxonomic information. sylph-tax can turn sylph's .tsv output into a taxonomic profile like Kraken or MetaPhlAn. sylph-tax does this by using custom taxonomy files to annotate sylph's output.

```sh
conda install -c bioconda sylph-tax
```

### 5.4.1 - Download the taxonomy files

```sh
mkdir taxonomy_file_folder
sylph-tax download --download-to taxonomy_file_folder
#You only have to do this once
```
### 5.4.2 - Use Sylph-tax

```sh
sylph-tax taxprof profiling_results.tsv -t taxonomy_file_folder/gtdb_r220_metadata.tsv.gz  -o output_prefix-sample1
#use the same meta_database
```
The output file is going to be a file *.sylphmpa* with the name asigned after the prefix. 

### 6 - Select closed contigs that represent bacterial genomes

We are now going to work with assembly_info.txt and assembly.fasta files which are the output from the Flye metagenome assembly.

### 6.1 - Select closed contigs from the *assembly_info.txt*

Closed contigs are assigend a "Y" in the 4th column of the file, which corresponds to circularized (Yes/No).

```sh
awk -F'\t' '$4 == "Y" {print $1}' assembly_info.txt > contigs_closed.txt
```

### 6.2 - Map and extract the circularized contigs in the *assembly.fasta* file from the Flye metagenomic assembly output.

```sh
seqkit grep -f contigs_closed.txt assembly.fasta > contigs_closed_assembly.fasta
```

### 6.3 - Filter by size 

Filter those contigs that have more than, e.g: 500.000 bp

```sh
seqkit seq -m 500000 -g contigs_closed_assembly.fasta > closed_genomes_filtered_500k.fasta
```
### Note: 

You can also extract particular contigs that are of interest to you. After running meta-Flye contigs are named as *"edge_XXX"* where the XXX represent a unique number for each contig. You can easily see the name that was assigned to your contig of interest selecting the contig in Bandage.

```sh
seqkit grep -p "contig_XXX" assembly.fasta > edge_XXX.fasta
```
```sh
cat edge_XXX.fasta closed_genomes_filtered_500k.fasta > closed_genomes_filtered_with_edge_XXX.fasta
```
### 6.4 - Verify all contigs included in your list

```sh
seqkit seq -n closed_genomes_filtered.fasta
```

### 6.5 - Split each contig in a different fasta file

```sh
seqkit split -i closed_genomes_filtered.fasta -O output_dir 
```

## 7 - Find missassebliesss

### 7.1 - Use minimap2 to map each genome assembly (e.g: closed_genome) back to filtered_metagenomics_reads.fastq file.

```sh
minimap2 -ax lr:hq -p 1 -o mapping-uper.sam -t 12 --secondary-seq closed_genome.fasta filtered_metagenomics_reads.fastq
```

### 7.2 - Convert the sam file to a bam file, sort it and index.

```sh
samtools view -bS -F 4 -@ 12 mapping-uper.sam -o mapping-super.bam
```
```sh
anvi-init-bam -T 12 -o mapping-super.bam mapping-super.bam
```

### 7.3 Run anvi-script-find-misassemblies

To run the command anvi-script-find-misassemblies you will need a bam-file of long reads mapped onto an assembly made from these long reads. The script creates two output .txt files: one for clipping events and one for contig regions with no coverage.

```sh
anvi-script-find-misassemblies -b mapping-super.bam -o misassemblies-super --just-do-it --clipping-ratio 0.5
```

### Note: Here you can learn  more about how to know more about the nature of the errors: https://merenlab.org/data/benchmarking-long-read-assemblers/

## 8 - Reorient each circular contig to start with the dnaA gene

```sh
dnaapler chromosome -i closed_genome.fasta -o dnaapler/closed_genome_reoriented.fasta --autocomplete nearest
```
The dnaapler chromosome command reorients your sequence to begin with the dnaA chromosomal replication initiator gene. The --autocomplete nearest option works if the tool cannot find a dnaA gene directly (e.g., if it is not annotated). It then searches for the gene most similar to dnaA and reorients the genome to begin with that gene.

## 9 - Check genome completeness

CheckM2 is a tool that evaluates the quality of your genome assembly by predicting completeness (the percentage of expected genes present in the genome) and contamination (the presence of sequences from other organisms).

```sh
checkm2 predict --input closed_genome_reoriented.fasta  --output-directory checkm2_output -x fasta
```

## 10 - Classify yout genome of interest using GTDB-Tk

### Note: Meke sure to have the last updated version which uses *skani* to calculate Average Nucleotide Identity (ANI) and incorporate it into the classify-wf.

```sh
gtdbtk --version
conda update -c bioconda gtdbtk
```

## 10.1 - Download and alias the GTDB-Tk reference data (226)

The most recent database is R226 at the moment, but this may vary, check the GTDB-Tk website for uptdates.

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
### Note: deactivate and activate again the conda enviornment.

### 10.2 - Make an output directory to save all the results

```sh
mkdir output_gtdb-tk
mkdir scratch_dir
```

### Note: set up an interactive job if needed
```sh
srun --job-name=gtdbtk --nodes=1 -c 12  --mem-per-cpu=4G -t 03:00:00 --partition=interactive --ntasks-per-node=1 --pty --export=ALL /usr/bin/bash
```

### 10.3 Run GTDB-Tk classify_wf

### Note: Add the path to your closed genome of interest after the flag --genome_dir, also add the correspoding paths to --output_dir and --scratch_dir.

```sh
gtdbtk classify_wf \
  --genome_dir \
  --output_dir \
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

The first command runs the GTDB-Tk classify workflow (classify_wf) using a MASH database (--mash_db) to quickly pre-filter the input genomes, identifying bacterial and archaeal genomes for classification. The second command also runs classify_wf but skips the MASH/ANI filtering step (--skip_ani_screen), meaning it will classify all genomes provided, regardless of their content, which is useful if you're confident that all input genomes are valid.

```sh
gtdbtk ani_rep --genome_dir --out_dir -x fasta --cpus 12
```

The command ani_rep only calculates the ANI between your genomes and reference genomes (R226 GTDB-Tk database) without performing taxonomic classification. This is useful if you only need similarity information without the full classification. 

The main output is a summary.tsv file.

## 11 - Pangenome analysis using Anvi'o

Once we identify to which species our genome of interest corresponds to, we can do a pangenome analysis. 

### 11.1 - Download species-cluster genomes

The first thing to do would be to download all genomes from the species cluster on GTDB-Tk. That can be done manually or using ncbi-genome-download.

```sh
ncbi-genome-download bacteria -s genbank -A GCA_963515415.1,GCA_927910885.1,GCA_937910925.1,GCA_963531205.1 --formats fasta --output-folder Species_cluster
```
You can always make a .txt file with all the accessions number of the genomes you want to download.

```sh
ACCESSIONS=$(paste -sd, accessions.txt)

ncbi-genome-download bacteria -s genbank -A "$ACCESSIONS" --formats fasta --output-folder .
```
Extract the .fna.gz files from each sub-folder created and move them to your current directory.

```sh
mv GCA_* ../../
```
```sh
for d in GCA_*; do
    fna_file=$(find "$d" -name "*.fna.gz" | head -n 1)
    if [ -f "$fna_file" ]; then
        mv "$fna_file" .
    else
        echo "No se encontró archivo .fna.gz en $d"
    fi
done
```
```sh
for f in *.fna; do
mv "$f" "${f%.fna}.fasta"
done
```
```sh
rm -rf genome
```

### 11.2 - Creat the fasta.txt file

The fasta.txt file that you need to create for running anvi pangenome workflow includes two columns, "name" and  "path" to the .fasta files and it is going to have as many lines as genomes included in the analysis. It is the main input for the pangenome workflow.

Remember to copy or move your genome assembly (closed_genome_reoriented.fasta) to the directory where you have the species genomes you donwloaded from NCBI.

```sh
echo -e "name\tpath" > fasta.txt
for f in *.fasta; do
    name=$(basename "$f" .fasta)
    path=$(realpath "$f")
    echo -e "${name}\t${path}" >> fasta.txt
done
```
Anvi'o does not accept all type of characters in the .txt file, if you get an error the first time running it, try to clean the .txt file with the following script.

```sh
awk 'BEGIN{OFS="\t"} NR==1 {print $1, $2; next} {gsub(/[^A-Za-z0-9_]/, "_", $1); print $1, $2}' fasta.txt > fasta.txt.cleaned
mv fasta.txt.cleaned fasta.txt
```
### 11.3 - Creating the pan-configs.json file

Create a .json file and edit it as needed.

```sh
nano pan-configs.json
```

```sh
{
    "fasta_txt": "fasta.txt",
    "anvi_gen_contigs_database": {
        "--project-name": "{group}",
        "--description": "",
        "--skip-gene-calling": "",
        "--ignore-internal-stop-codons": "",
        "--skip-mindful-splitting": "",
        "--contigs-fasta": "",
        "--split-length": "",
        "--kmer-size": "",
        "--skip-predict-frame": "",
        "--prodigal-translation-table": "",
        "threads": 2
    },
    "centrifuge": {
        "threads": 2,
        "run": "",
        "db": ""
    },
    "anvi_run_hmms": {
        "run": true,
        "threads": 2,
        "--also-scan-trnas": true,
        "--installed-hmm-profile": "",
        "--hmm-profile-dir": ""
},
"anvi_run_kegg_kofams": {
  "run": true,
  "threads": 2,
  "--kegg-data-dir": "/home/exacloud/gscratch/BakerLab/charlom/github/anvio/kegg",
  "--hmmer-program": "",
  "--keep-all-hits": true,
  "--log-bitscores": true,
  "--just-do-it": true
},
"anvi_run_ncbi_cogs": {
  "run": true,
  "threads": 2,
  "--cog-data-dir": "",
  "--temporary-dir-path": "",
  "--search-with": ""
},
    "anvi_run_scg_taxonomy": {
        "run": true,
        "threads": 2,
        "--scgs-taxonomy-data-dir": ""
    },
    "anvi_run_trna_scan": {
        "run": false,
        "threads": 2,
        "--trna-cutoff-score": ""
        },
    "anvi_script_reformat_fasta": {
        "run": true,
        "--prefix": "{group}",
        "--simplify-names": true,
        "--keep-ids": "",
        "--exclude-ids": "",
        "--min-len": "",
        "--seq-type": "",
        "threads": ""
    },
    "emapper": {
        "--database": "bact",
        "--usemem": true,
        "--override": true,
        "path_to_emapper_dir": "",
        "threads": ""
    },
    "anvi_script_run_eggnog_mapper": {
    "--use-version": "0.12.6",
        "run": "",
        "--cog-data-dir": "",
        "--drop-previous-annotations": "",
        "threads": ""
    },
    "anvi_get_sequences_for_hmm_hits": {
        "--return-best-hit": true,
        "--align-with": "famsa",
        "--concatenate-genes": true,
        "--get-aa-sequences": true,
        "--hmm-sources": "Bacteria_71",
        "--separator": "",
        "--min-num-bins-gene-occurs": "",
        "--max-num-genes-missing-from-bin": "",
        "--gene-names": "",
        "threads": 10
    },
    "trimal": {
        "-gt": 0.5,
        "additional_params": "",
        "threads": 10
    },
    "iqtree": {
        "threads": 10,
        "-m": "WAG",
        "-bb": 1000,
        "additional_params": ""
    },
    "anvi_pan_genome": {
        "threads": 10,
        "--project-name": "",
        "--genome-names": "",
        "--skip-alignments": "",
        "--align-with": "",
        "--exclude-partial-gene-calls": "",
        "--use-ncbi-blast": "",
        "--minbit": "",
        "--mcl-inflation": 9,
        "--min-occurrence": "",
        "--min-percent-identity": "",
        "--description": "",
        "--overwrite-output-destinations": "",
        "--skip-hierarchical-clustering": "",
        "--enforce-hierarchical-clustering": "",
        "--distance": "",
        "--linkage": ""
         },
    "import_phylogenetic_tree_to_pangenome": {
        "tree_name": "phylogeny",
        "--just-do-it": "",
        "threads": ""
    },
    "anvi_compute_genome_similarity": {
        "run": false,
        "additional_params": "",
        "threads": ""
    },
    "gen_external_genome_file": {
        "threads": ""
    },
    "export_gene_calls_for_centrifuge": {
        "threads": ""
    },
    "anvi_import_taxonomy_for_genes": {
        "threads": ""
    },
    "annotate_contigs_database": {
        "threads": ""
    },
    "anvi_get_sequences_for_gene_calls": {
        "threads": ""
    },
    "gunzip_fasta": {
        "threads": ""
    },
    "reformat_external_gene_calls_table": {
        "threads": ""
        },
    "import_external_functions": {
        "threads": ""
    },
    "anvi_run_pfams": {
        "run": "",
        "--pfam-data-dir": "",
        "threads": ""
    },
    "anvi_gen_genomes_storage": {
        "--gene-caller": "",
        "threads": ""
    },
    "anvi_get_sequences_for_gene_clusters": {
        "--gene-cluster-id": "",
        "--gene-cluster-ids-file": "",
        "--collection-name": "",
        "--bin-id": "",
        "--min-num-genomes-gene-cluster-occurs": "",
        "--max-num-genomes-gene-cluster-occurs": "",
        "--min-num-genes-from-each-genome": "",
        "--max-num-genes-from-each-genome": "",
        "--max-num-gene-clusters-missing-from-genome": "",
        "--min-functional-homogeneity-index": "",
        "--max-functional-homogeneity-index": "",
        "--min-geometric-homogeneity-index": "",
        "--max-geometric-homogeneity-index": "",
        "--add-into-items-additional-data-table": "",
        "--concatenate-gene-clusters": "",
        "--separator": "",
        "--align-with": "",
        "threads": ""
    },
    "project_name": "Prevotella_pangenome",
    "internal_genomes": "",
    "external_genomes": "external_genomes.txt",
    "sequence_source_for_phylogeny": "",
    "output_dirs": {
        "FASTA_DIR": "01_FASTA",
        "CONTIGS_DIR": "02_CONTIGS",
        "PHYLO_DIR": "01_PHYLOGENOMICS",
        "PAN_DIR": "03_PAN",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": 24,
    "config_version": "3",
    "workflow_name": "pangenomics"
}
```

### NOTES: 
* Make sure that the fasta.txt is accesible or add the path to it.
* Use a different project name for your different analysis.
* Do not worry about the external genomes file, it will be created for you during the run.
* Make sure max_threads matches with the cores you asign in the following slurm script.

### 11.4 - Run Anvi'o pangenomics workflow with Slurm

Once you have your .json file and fasta.txt file you can run the  Anvi'o pangenomics workflow.

```sh
anvi-run-workflow -w pangenomics -c pan-configs.json \
  --additional-params --jobs 24 --keep-going --rerun-incomplete \
  --cores 24 \
  --cluster 'sbatch --job-name=Rothia_pangenome --output={log} --error={log} --nodes=1 --ntasks=1 --cpus-per-task=1'
```
The flags --keep-going --rerun-incomplete will allow you to re-run the workflow from where you stopped if you ran into an error and had to repeat the analysis or just if you want to add new genomes to the analysis. Make sure to delete or rename the 03_PAN folder previosuly if you are running the workflow.

Running this workflow will create four main direcotories: *00_LOGS* which includes all errors or information from the run, *01_FASTA* with all the reformatted fasta that Anvi'o uses as an input, *02_CONTIGS* which includes .db files with all the information related to the gene anotations such as keggs, hmms, and cogs that you can explore further using "anvi-db-info". Finally *03_PAN* directory includes two main files *pangenome-PAN.db* and *pangenome-GENOMES.db* which are going to be essential for the following summary and visualization steps.

### 11.4 - Compute genome similarity

Located in the 03_PAN folder run this script to compute genome similarity.

```sh
anvi-compute-genome-similarity --external-genomes ../external_genomes.txt          --program fastANI --output-dir ANI --num-threads 6 --pan-db pangenome-PAN.db
```

### 11.5 - Visualize and summarize the results

Visualize and customize the pangenome from the Anvi'o interactive interface.

```sh
anvi-display-pan -g Prevotella_pangenome-GENOMES.db \
                 -p Prevotella_pangenome-PAN.db
```
*Use a proxy jump from a second terminal to visualize it in your browser if you are working on the ARC*
```sh
ssh -L 8081:localhost:8081 -J user@acc.ohsu.edu user@arc-infra-3.ohsu.edu
```

Add a default collection and summarize the results, this will output a summary.txt file.

```sh
anvi-script-add-default-collection -p pangenome-PAN.db -C DEFAULT
```
```sh
anvi-summarize -p pangenome-PAN.db \
               -g -GENOMES.db \
               -C DEFAULT \
               -o SUMMARY
```

### 12 - Obtain core genes and exclusive gene clusters.

1) Remove any genome you do not want to take into account in the analysis. In my case it will be an external genome I added to root the phylogenetic tree.

Column 4 corresponds to genome_names in the summary.txt file, so we can make a filtered summary.txt file deleting lines with the name of the genome we want to exclude.

```sh
awk -F'\t' '$4 != "GCF_001836735_1_external_genomic"' Rothia_pangenome_gene_clusters_summary.txt > Rothia_pangenome_gene_clusters_summary_filtered.txt
```
2) Identify gene clusters that are in all genomes of the species included in the analysis but NOT in the genome of interest. The output is a .txt file with all the gene cluster IDs. 

```sh
awk -F'\t' '
NR > 1 {
  gc = $2        # gene_cluster_id
  genome = $4    # genome_name
  gcs[gc][genome] = 1
}
END {
  for (gc in gcs) {
    n = 0
    for (genome in gcs[gc]) {
      if (genome == "Rothia") {
        n = -99
        break
      }
      n++
    }
    if (n == 30) {
      print gc
    }
  }
}
' Rothia_pangenome_gene_clusters_summary_filtered.txt > gene_clusters_NOT_in_Rothia.txt
```
Note: Adjust this parameter *(genome == "Rothia")* with the name of your genome of interest as it appears in the summary file. And also adjust *(n == 30)* according to the number of species in the analysis. In my case, once removed the external genome I had 30 genomes apart from my genome of interest (Rothia), 31 in total.

Add all the information of the lines corresponding to the gene clusters IDs as well as the header.

```sh
head -n 1 Rothia_pangenome_gene_clusters_summary.txt > lines_clusters_NOT_in_Rothia.txt
grep -Ff gene_clusters_NOT_in_Rothia.txt Rothia_pangenome_gene_clusters_summary.txt >> gene_clusters_NOT_in_Rothia_final.txt
```

3) Identify gene clusters that are exclusive from the genome of interest.

```sh
awk '$4 ~ /Rothia/ {print $2}' Rothia_pangenome_gene_clusters_summary_filtered.txt | sort | uniq > clusters_in_Rothia.txt
```

```sh
awk '$4 !~ /Rothia/ {print $2}' Rothia_pangenome_gene_clusters_summary_filtered.txt | sort | uniq > clusters_in_others.txt
```

```sh
comm -23 clusters_in_Rothia.txt clusters_in_others.txt > gene_clusters_exclusive_Rothia.txt
```

Add all the information of the lines corresponding to the gene clusters IDs as well as the header.

```sh
head -n 1 Rothia_pangenome_gene_clusters_summary_filtered.txt > lines_exclusive_Rothia.txt
grep -Ff gene_clusters_exclusive_Rothia.txt Rothia_pangenome_gene_clusters_summary_filtered.txt >> gene_clusters_exclusive_Rothia_final.txt 
```

4) Obtain core genes

Obtain the unique combinations of gene_cluster_id and genome_name

```sh
awk '{print $2, $4}' Rothia_pangenome_gene_clusters_summary_filtered.txt | sort -u > cluster_genome_pairs.txt
```

Count how many genomes each gene_cluster_id appears in
```sh
cut -d' ' -f1 cluster_genome_pairs.txt | sort | uniq -c > cluster_counts.txt
```

Filter those that appear in all genomes (CORE genes), replace the number on the $1 == X according to the number of genomes included. 

```sh
awk '$1 == 31 {print $2}' cluster_counts.txt > core_gene_clusters.txt
```

Add all the information of the lines corresponding to the gene clusters IDs as well as the header.

```sh
head -n 1 Rothia_pangenome_gene_clusters_summary_filtered.txt > core_genes_summary.txt
grep -Ff core_gene_clusters.txt Rothia_pangenome_gene_clusters_summary_filtered.txt >> core_genes_Rothia_final.txt
```

## Obteain core genes and build a phylogenomic tree using Anvi'o

### Extract core gene sequences from your Anvi'o pangenome

```sh
anvi-get-sequences-for-gene-clusters -g your_pangenome-GENOMES.db -p your_pangenome-PAN.db -o better.core.fa --min-num-genomes-gene-cluster-occurs 8 --max-num-genes-from-each-genome 1 --min-geometric-homogeneity-index 1 --max-functional-homogeneity-index 0.99 --concatenate-gene-cluster
```
Filters:

--min-num-genomes-gene-cluster-occurs 8: selects gene clusters found in at least 8 genomes (defines your "core").

--max-num-genes-from-each-genome 1: ensures one copy per genome (avoids paralogs).

--min-geometric-homogeneity-index 1: requires perfect geometric homogeneity (genes are evenly distributed).

--max-functional-homogeneity-index 0.99: ensures gene clusters are functionally consistent.

--concatenate-gene-cluster: stitches gene sequences from each genome together to prepare for phylogenetic analysis.

The output file is a fasta file: better.core.fa where core genes from each genome are concatenated into a single sequence per genome.

### Build a maximum likelihood phylogenetic tree using the aligned core gene sequences.

```sh
iqtree -s better-core.fa -nt 7 -m WAG -bb 1000
```

### Visaluze tree in Anvi'o

```sh
anvi-gen-phylogenomic-profile -i better.core.fa -o phylogenomic-profile-Nanosynbacter --prefix better_core
```

```sh
anvi-interactive -p phylogenomic-profile-Nanosynbacter.db -t better.core.fa.contree --manual
```

### More on editing your Anvi'o pangenome figure

1) Export the current layer metadata from the pangenome database. 

This command exports existing "layer" metadata from your pangenome database. The -t layers flag tells Anvi'o to focus on layer-type data, which are annotations associated with genomes (like taxonomy, source, assembly level, others). The output file pan_layers.txt will contain the genome names and the currently available metadata, if any.

```sh
anvi-export-misc-data -p 03_PAN/Nanosynbacter_pangenome_all-PAN.db -t layers -o pan_layers.txt
```

2) Create your custom layer metadata file

Make a .txt file that adds layers of information, like level_of_assembly, or display_name. 

*The genome names in the first column in both files need to be the same The file must be tab-delimited and include a header row.*

```sh
anvi-import-misc-data -p 03_PAN/Nanosynbacter_pangenome_all-PAN.db -t layers names.txt
```

3) Import your custom metadata into the pangenome database

This imports your custom metadata (like genome display names or assembly level) into the .PAN.db database. The -t layers tells Anvi’o that you're modifying layer data. After this step, your metadata will be available in the pangenome interactive interface.

```sh
anvi-display-pan -g Nanosynbacter_pangenome_all-GENOMES.db -p Nanosynbacter_pangenome_all-PAN.db
```

Run anvi-display-pan workflow.
