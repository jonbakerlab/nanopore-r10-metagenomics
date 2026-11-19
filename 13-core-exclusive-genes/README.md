# 13. Obtain core genes and exclusive gene clusters

The examples below use a genome of interest named `Rothia` and a pangenome summary file, `Rothia_pangenome_gene_clusters_summary.txt`, from [step 12](../12-pangenome/). Adjust names and genome counts for your own analysis.

## 1) Remove genomes you do not want to include

In our case this was an external genome added to root the phylogenetic tree. Column 4 of the summary file is `genome_name`, so we can make a filtered summary by deleting lines with the name of the genome to exclude.

```sh
awk -F'\t' '$4 != "GCF_001836735_1_external_genomic"' Rothia_pangenome_gene_clusters_summary.txt > Rothia_pangenome_gene_clusters_summary_filtered.txt
```

## 2) Gene clusters in all genomes of the species but NOT in the genome of interest

The output is a `.txt` file with all the gene cluster IDs.

```sh
awk -F'\t' '
NR > 1 {
  gc = $2        # gene_cluster_id
  genome = $4    # genome_name
  gcs[gc][genome] = 1
}
END {
  for (gc in gcs) {
    n = 0
    for (genome in gcs[gc]) {
      if (genome == "Rothia") {
        n = -99
        break
      }
      n++
    }
    if (n == 30) {
      print gc
    }
  }
}
' Rothia_pangenome_gene_clusters_summary_filtered.txt > gene_clusters_NOT_in_Rothia.txt
```
Note: adjust `(genome == "Rothia")` to the name of your genome of interest as it appears in the summary file, and adjust `(n == 30)` to the number of genomes in the analysis. In our case, once the external genome was removed, there were 30 genomes apart from the genome of interest (Rothia), 31 in total.

Add the full lines for these gene cluster IDs, along with the header:

```sh
head -n 1 Rothia_pangenome_gene_clusters_summary.txt > gene_clusters_NOT_in_Rothia_final.txt
grep -Ff gene_clusters_NOT_in_Rothia.txt Rothia_pangenome_gene_clusters_summary.txt >> gene_clusters_NOT_in_Rothia_final.txt
```

## 3) Gene clusters exclusive to the genome of interest

```sh
awk '$4 ~ /Rothia/ {print $2}' Rothia_pangenome_gene_clusters_summary_filtered.txt | sort | uniq > clusters_in_Rothia.txt
```
```sh
awk '$4 !~ /Rothia/ {print $2}' Rothia_pangenome_gene_clusters_summary_filtered.txt | sort | uniq > clusters_in_others.txt
```
```sh
comm -23 clusters_in_Rothia.txt clusters_in_others.txt > gene_clusters_exclusive_Rothia.txt
```

Add the full lines for these gene cluster IDs, along with the header:

```sh
head -n 1 Rothia_pangenome_gene_clusters_summary_filtered.txt > gene_clusters_exclusive_Rothia_final.txt
grep -Ff gene_clusters_exclusive_Rothia.txt Rothia_pangenome_gene_clusters_summary_filtered.txt >> gene_clusters_exclusive_Rothia_final.txt
```

## 4) Obtain core genes

Obtain the unique combinations of `gene_cluster_id` and `genome_name`:

```sh
awk '{print $2, $4}' Rothia_pangenome_gene_clusters_summary_filtered.txt | sort -u > cluster_genome_pairs.txt
```

Count how many genomes each `gene_cluster_id` appears in:
```sh
cut -d' ' -f1 cluster_genome_pairs.txt | sort | uniq -c > cluster_counts.txt
```

Keep those that appear in all genomes (core genes); replace `31` in `$1 == 31` with the number of genomes included:

```sh
awk '$1 == 31 {print $2}' cluster_counts.txt > core_gene_clusters.txt
```

Add the full lines for these gene cluster IDs, along with the header:

```sh
head -n 1 Rothia_pangenome_gene_clusters_summary_filtered.txt > core_genes_Rothia_final.txt
grep -Ff core_gene_clusters.txt Rothia_pangenome_gene_clusters_summary_filtered.txt >> core_genes_Rothia_final.txt
```
