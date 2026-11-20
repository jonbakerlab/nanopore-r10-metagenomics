# 08. Find misassemblies

## Map reads back to each genome assembly

Use minimap2 to map the filtered reads back to each genome assembly (e.g. `closed_genome.fasta`):

```sh
minimap2 -ax lr:hq -p 1 -o mapping-super.sam -t 12 --secondary-seq closed_genome.fasta filtered_metagenomics_reads.fastq
```

## Convert, sort and index

Convert the SAM file to a BAM file, then sort and index it with anvi'o:

```sh
samtools view -bS -F 4 -@ 12 mapping-super.sam -o mapping-super.bam
```
```sh
anvi-init-bam -T 12 -o mapping-super.bam mapping-super.bam
```

## Run anvi-script-find-misassemblies

To run `anvi-script-find-misassemblies` you need a BAM file of long reads mapped onto an assembly made from those same reads. The script creates two output `.txt` files: one for clipping events and one for contig regions with no coverage.

```sh
anvi-script-find-misassemblies -b mapping-super.bam -o misassemblies-super --just-do-it --clipping-ratio 0.5
```

Note: to learn more about the nature of these errors, see https://merenlab.org/data/benchmarking-long-read-assemblers/
