## This code is to filter .pod5 files from a MinION run into two buckets, one from channels 1-256 (adaptive sampling--depletion mode to deplete reads mapping to the human genome) and one from channels 257-512.

```bash
# Generate a table containing read IDs and channels:
pod5 view *.pod5 \
 --include "read_id,channel" \
 --output reads.tsv

# Then create two read ID lists:
awk -F'\t' '$2 >= 1 && $2 <= 256 {print $1}' reads.tsv > channels_1_256.txt

awk -F'\t' '$2 >= 257 && $2 <= 512 {print $1}' reads.tsv > channels_257_512.txt

# Extract the reads into separate POD5 files

pod5 filter *.pod5 \
 --ids channels_1_256.txt \
 --output summed-reads/channels_1_256.pod5

pod5 filter *.pod5 \
 --ids channels_257_512.txt \
 --output summed-reads/channels_257_512.pod5
```

The sizes of the two files was the following:
```bash
ls -lh summed-reads/
-rw-r-----. 1 bakerjo bakerjo 31G Jul 13 09:39 channels_1_256.pod5
-rw-r-----. 1 bakerjo bakerjo 64G Jul 13 09:41 channels_257_512.pod5
```

## Using Myloasm instead of Flye for genome assembly
Myloasm is a recently developed assembler designed specifically for long-read sequencing data. In the authors’ benchmarks of real Oxford Nanopore metagenomes, Myloasm recovered approximately three times as many complete circular contigs as the next-best assembler tested.

Shaw, J., Marin, M.G. & Li, H. High-resolution metagenome assembly for modern long reads with myloasm. _Nat Biotechnol_ (2026). https://doi.org/10.1038/s41587-026-03053-z

## 6-Select circular contigs from the Myloasm assembly
The ```assembly_primary.fa``` file is the main output contig file produced by Myloasm and is the file that should usually be used for downstream analysis.

In Myloasm, confidently circular contigs are labeled directly in the FASTA header as:

_circular-yes_



```bash
# Extract all contigs labeled circular-yes
seqkit grep \
-n \
-r \
-p '_circular-yes_'\
assembly_primary.fa > contigs_closed.fa
```

