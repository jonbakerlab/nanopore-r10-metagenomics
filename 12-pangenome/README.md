# 12. Pangenome analysis using anvi'o

Once we identify which species our genome of interest corresponds to (see [step 11](../11-gtdbtk/)), we can run a pangenome analysis.

## Download species-cluster genomes

First, download all genomes from the species cluster in GTDB. This can be done manually or using ncbi-genome-download.

```sh
ncbi-genome-download bacteria -s genbank -A GCA_963515415.1,GCA_927910885.1,GCA_937910925.1,GCA_963531205.1 --formats fasta --output-folder Species_cluster
```
You can also make a `.txt` file with the accession numbers of all the genomes you want to download:

```sh
ACCESSIONS=$(paste -sd, accessions.txt)

ncbi-genome-download bacteria -s genbank -A "$ACCESSIONS" --formats fasta --output-folder .
```
Extract the `.fna.gz` files from each sub-folder created and move them to your current directory:

```sh
mv GCA_* ../../
```
```sh
for d in GCA_*; do
    fna_file=$(find "$d" -name "*.fna.gz" | head -n 1)
    if [ -f "$fna_file" ]; then
        mv "$fna_file" .
    else
        echo "No .fna.gz file found in $d"
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

## Create the fasta.txt file

The `fasta.txt` file required by the anvi'o pangenome workflow has two columns, "name" and "path" (to the `.fasta` files), with one line per genome in the analysis. It is the main input for the pangenome workflow.

Remember to copy or move your genome assembly (`closed_genome_reoriented.fasta`) to the directory containing the species genomes downloaded from NCBI.

```sh
echo -e "name\tpath" > fasta.txt
for f in *.fasta; do
    name=$(basename "$f" .fasta)
    path=$(realpath "$f")
    echo -e "${name}\t${path}" >> fasta.txt
done
```
anvi'o does not accept all characters in the `.txt` file. If you get an error the first time running it, clean the file with:

```sh
awk 'BEGIN{OFS="\t"} NR==1 {print $1, $2; next} {gsub(/[^A-Za-z0-9_]/, "_", $1); print $1, $2}' fasta.txt > fasta.txt.cleaned
mv fasta.txt.cleaned fasta.txt
```

## Create the pan-configs.json file

Copy [`pan-configs.json`](pan-configs.json) from this folder and edit it as needed (e.g. `project_name`, thread counts, `external_genomes`).

Notes:
* Make sure `fasta.txt` is accessible, or add the path to it.
* Use a different project name for each analysis.
* Do not worry about the external genomes file; it is created during the run.
* Make sure `max_threads` matches the cores assigned in the SLURM script.

## Run the anvi'o pangenomics workflow with SLURM

Once you have your `.json` file and `fasta.txt` file, run the workflow with [`run_pangenome_workflow.sh`](run_pangenome_workflow.sh):

```sh
bash run_pangenome_workflow.sh
```
The flags `--keep-going --rerun-incomplete` allow you to re-run the workflow from where you stopped if you hit an error, or to add new genomes to the analysis. Make sure to delete or rename the `03_PAN` folder beforehand if you are re-running the workflow.

Running this workflow creates four main directories: `00_LOGS` (errors and run information), `01_FASTA` (the reformatted FASTA files used as anvi'o input), `02_CONTIGS` (`.db` files with gene annotations such as KEGG, HMMs and COGs, explorable with `anvi-db-info`), and `03_PAN`, which includes the two main files `pangenome-PAN.db` and `pangenome-GENOMES.db` that are essential for the following summary and visualization steps.

## Compute genome similarity

From the `03_PAN` folder, run:

```sh
anvi-compute-genome-similarity --external-genomes ../external_genomes.txt --program fastANI --output-dir ANI --num-threads 6 --pan-db pangenome-PAN.db
```

## Visualize and summarize the results

Visualize and customize the pangenome from the anvi'o interactive interface.

```sh
anvi-display-pan -g Prevotella_pangenome-GENOMES.db \
                 -p Prevotella_pangenome-PAN.db
```
Use a proxy jump from a second terminal to view it in your browser if you are working on the ARC:
```sh
ssh -L 8081:localhost:8081 -J user@acc.ohsu.edu user@arc-infra-3.ohsu.edu
```

Add a default collection and summarize the results; this outputs a `summary.txt` file.

```sh
anvi-script-add-default-collection -p pangenome-PAN.db -C DEFAULT
```
```sh
anvi-summarize -p pangenome-PAN.db \
               -g -GENOMES.db \
               -C DEFAULT \
               -o SUMMARY
```
