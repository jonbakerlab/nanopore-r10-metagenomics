# 03. Install and set up conda environments

Creating conda environments to compartmentalize all packages needed for the analysis is recommended. For more information on installing conda and mamba see the "Introduction to Conda and Mamba" repository (https://github.com/jonbakerlab/conda-microbiome-tutorial). Below is an example of creating a conda environment and installing minimap2. This needs to be repeated for each program needed.

## Create a conda environment

```sh
conda create -n minimap2 -c bioconda minimap2
conda activate minimap2
conda install bioconda::minimap2
```

## Programs needed

minimap2, samtools, parallel, seqkit, Flye, Myloasm, Sylph, dnaapler, bakta, anvi'o (developer version).
