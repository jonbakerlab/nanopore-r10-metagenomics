# 10. Check genome completeness

CheckM2 evaluates the quality of a genome assembly by predicting completeness (the percentage of expected genes present in the genome) and contamination (the presence of sequences from other organisms).

```sh
checkm2 predict --input closed_genome_reoriented.fasta --output-directory checkm2_output -x fasta
```
