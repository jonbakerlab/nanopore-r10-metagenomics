# 01. Adaptive sampling and pod5 channel split

## A note on adaptive sampling

If you are interested in obtaining a higher ratio of microbial DNA reads to human DNA reads, [adaptive sampling](https://nanoporetech.com/es/document/adaptive-sampling) can help a lot. First, download the human reference genome:

```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

When setting up your run, select the "Deplete Sequences" option in the Adaptive Sampling section and attach the human reference genome file as the alignment reference. This file must be in the same MinKNOW directory to be accessible.

## Splitting pod5 files by channel

This code filters the `.pod5` files from a MinION run into two buckets: one from channels 1-256 (adaptive sampling in depletion mode, to deplete reads mapping to the human genome) and one from channels 257-512. The commands are in [`split_pod5_by_channel.sh`](split_pod5_by_channel.sh).

```bash
bash split_pod5_by_channel.sh
```

The sizes of the two output files were:

```bash
ls -lh summed-reads/
-rw-r-----. 1 bakerjo bakerjo 31G Jul 13 09:39 channels_1_256.pod5
-rw-r-----. 1 bakerjo bakerjo 64G Jul 13 09:41 channels_257_512.pod5
```
