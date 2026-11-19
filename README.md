# nanopore-r10-metagenomics

Companion repository for our R10 Nanopore metagenomics study. It contains the code and Supplemental Data for a pipeline that produces complete, high-quality metagenome-assembled genomes (MAGs) from saliva.

## Introduction

High-molecular-weight genomic DNA was extracted from saliva and a metagenomic long-read library was prepared using the Ligation Sequencing gDNA - Native Barcoding Kit 24 V14 (Oxford Nanopore Technologies). It was sequenced on a MinION using an R10 flow cell (Oxford Nanopore Technologies) with the adaptive sampling option.

This pipeline takes the pod5 files output by a MinION sequencer, performs super high accuracy basecalling using Dorado, and aims to obtain complete, high-quality microbial genomes from saliva. This work was done at the Baker Lab with support from Dr. Jonathon Baker and his team, especially Dr. Matthew Barbisan.

## Overview
Each analysis step is documented in its own folder, alongside its associated scripts and files:

1. [Adaptive sampling and pod5 channel split](01-adaptive-sampling/)
2. [Super high accuracy basecalling (Dorado)](02-basecalling/)
3. [Conda environments](03-conda-environments/)
4. [Human read removal (minimap2)](04-human-read-removal/)
5. [Genome assembly (Flye or Myloasm)](05-assembly/)
6. [Taxonomic profiling (Sylph)](06-sylph-profiling/)
7. [Selecting closed contigs](07-closed-contigs/)
8. [Misassembly detection](08-misassemblies/)
9. [Reorienting contigs to dnaA (dnaapler)](09-reorient-dnaapler/)
10. [Genome completeness (CheckM2)](10-checkm2/)
11. [Classification (GTDB-Tk)](11-gtdbtk/)
12. [Pangenome analysis (anvi'o)](12-pangenome/)
13. [Core and exclusive gene clusters](13-core-exclusive-genes/)
14. [Core genes phylogeny](14-core-genes-phylogeny/)
15. [Editing the pangenome figure](15-pangenome-figure/)

## Software

Programs needed: Dorado, minimap2, samtools, parallel, seqkit, Flye and/or Myloasm, Sylph (and sylph-tax), dnaapler, CheckM2, GTDB-Tk, bakta, ncbi-genome-download, anvi'o (developer version), IQ-TREE. See [03-conda-environments](03-conda-environments/).
