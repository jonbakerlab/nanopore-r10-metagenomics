# 04. Filter out human reads using minimap2

Since we are interested in microbial genomic DNA, we filter out all human reads.

## Download the human reference genome

```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

## Index the human genome

```sh
minimap2 -d human.mmi GRCh38.primary_assembly.genome.fa
```

## Align the reads to the human genome

Make sure to unzip your `fastq.gz` files first.

```sh
minimap2 -ax lr:hq human.mmi barcode01/*.fastq > aligned.sam
```

## Filter out human reads

These commands separate mapped (human) and unmapped (non-human) reads. The output is a BAM file that needs to be converted back to FASTQ.

```sh
samtools view -@ 14 -b -F 4 aligned.sam > human_reads.bam
samtools view -@ 14 -b -f 4 aligned.sam > non_human_reads.bam
```
```sh
samtools fastq -@ 14 human_reads.bam > human_reads.fastq
samtools fastq -@ 14 non_human_reads.bam > filtered_metagenomics_reads.fastq
```

## Notes

It is worth keeping the `human_reads.fastq` file to compare file sizes with the `filtered_metagenomics_reads.fastq` file, and to check that the sum of the sizes of both files makes sense relative to the unfiltered file you started with.

You can view file sizes using:
```sh
ls -lh
```
or
```sh
du -sh *yourfile*
```

Concatenate two or more filtered FASTQ files from different runs if needed:

```sh
cat filtered_metagenomics_reads_1.fastq filtered_metagenomics_reads_2.fastq > filtered_metagenomics_reads_merged.fastq
```

Verify the integrity of the merged file:

```sh
parallel -j 14 grep -c "^@" ::: filtered_metagenomics_reads_1.fastq filtered_metagenomics_reads_2.fastq
```
```sh
parallel -j 14 grep -c "^@" ::: filtered_metagenomics_reads_merged.fastq
```
