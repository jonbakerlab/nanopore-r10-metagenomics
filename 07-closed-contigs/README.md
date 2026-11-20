# 07. Select closed contigs that represent bacterial genomes

Closed (circular) contigs are selected from the assembly. Follow the section matching the assembler used in [step 5](../05-assembly/).

## Option A: Flye

Use the `assembly_info.txt` and `assembly.fasta` files from the Flye metagenome assembly.

### Select closed contigs from `assembly_info.txt`

Closed contigs are assigned a "Y" in the 4th column of the file, which corresponds to circularized (Yes/No).

```sh
awk -F'\t' '$4 == "Y" {print $1}' assembly_info.txt > contigs_closed.txt
```

### Extract the circularized contigs from `assembly.fasta`

```sh
seqkit grep -f contigs_closed.txt assembly.fasta > contigs_closed_assembly.fasta
```

### Filter by size

Keep contigs larger than, e.g., 500,000 bp:

```sh
seqkit seq -m 500000 -g contigs_closed_assembly.fasta > closed_genomes_filtered_500k.fasta
```

### Note

You can also extract particular contigs of interest. After running metaFlye, contigs are named `edge_XXX`, where XXX is a unique number for each contig. You can see the name assigned to your contig of interest by selecting the contig in Bandage.

```sh
seqkit grep -p "edge_XXX" assembly.fasta > edge_XXX.fasta
```
```sh
cat edge_XXX.fasta closed_genomes_filtered_500k.fasta > closed_genomes_filtered_with_edge_XXX.fasta
```

### Verify all contigs included in your list

```sh
seqkit seq -n closed_genomes_filtered.fasta
```

### Split each contig into a different FASTA file

```sh
seqkit split -i closed_genomes_filtered.fasta -O output_dir
```

## Option B: Myloasm

In Myloasm, confidently circular contigs are labeled directly in the FASTA header as `_circular-yes_`. Extract them from `assembly_primary.fa`:

```bash
seqkit grep \
  -n \
  -r \
  -p '_circular-yes_' \
  assembly_primary.fa > contigs_closed.fa
```
