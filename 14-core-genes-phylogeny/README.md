# 14. Core genes phylogeny

Obtain core genes and build a phylogenomic tree using anvi'o.

## Extract core gene sequences from your anvi'o pangenome

```sh
anvi-get-sequences-for-gene-clusters -g your_pangenome-GENOMES.db -p your_pangenome-PAN.db -o better-core.fa --min-num-genomes-gene-cluster-occurs 8 --max-num-genes-from-each-genome 1 --min-geometric-homogeneity-index 1 --max-functional-homogeneity-index 0.99 --concatenate-gene-clusters
```
Filters:

* `--min-num-genomes-gene-cluster-occurs 8`: selects gene clusters found in at least 8 genomes (defines your "core").
* `--max-num-genes-from-each-genome 1`: ensures one copy per genome (avoids paralogs).
* `--min-geometric-homogeneity-index 1`: requires perfect geometric homogeneity (genes are evenly distributed).
* `--max-functional-homogeneity-index 0.99`: ensures gene clusters are functionally consistent.
* `--concatenate-gene-clusters`: stitches gene sequences from each genome together to prepare for phylogenetic analysis.

The output is a FASTA file, `better-core.fa`, where the core genes from each genome are concatenated into a single sequence per genome.

## Build a maximum-likelihood phylogenetic tree

Use the aligned core gene sequences:

```sh
iqtree -s better-core.fa -nt 7 -m WAG -bb 1000
```

## Visualize the tree in anvi'o

```sh
anvi-gen-phylogenomic-profile -i better-core.fa -o phylogenomic-profile-Nanosynbacter --prefix better_core
```

```sh
anvi-interactive -p phylogenomic-profile-Nanosynbacter.db -t better-core.fa.contree --manual
```
