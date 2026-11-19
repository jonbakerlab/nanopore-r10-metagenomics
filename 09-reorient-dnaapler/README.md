# 09. Reorient each circular contig to start with the dnaA gene

```sh
dnaapler chromosome -i closed_genome.fasta -o dnaapler/closed_genome_reoriented.fasta --autocomplete nearest
```

The `dnaapler chromosome` command reorients your sequence to begin with the dnaA chromosomal replication initiator gene. The `--autocomplete nearest` option works if the tool cannot find a dnaA gene directly (e.g. if it is not annotated): it then searches for the gene most similar to dnaA and reorients the genome to begin with that gene.
