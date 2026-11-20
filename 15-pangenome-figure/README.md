# 15. Editing your anvi'o pangenome figure

## 1) Export the current layer metadata from the pangenome database

This exports existing "layer" metadata from your pangenome database. The `-t layers` flag tells anvi'o to focus on layer-type data, which are annotations associated with genomes (such as taxonomy, source, assembly level). The output file `pan_layers.txt` contains the genome names and any currently available metadata.

```sh
anvi-export-misc-data -p 03_PAN/Nanosynbacter_pangenome_all-PAN.db -t layers -o pan_layers.txt
```

## 2) Create your custom layer metadata file

Make a `.txt` file that adds layers of information, such as `level_of_assembly` or `display_name`.

The genome names in the first column of both files need to be the same. The file must be tab-delimited and include a header row.

## 3) Import your custom metadata into the pangenome database

This imports your custom metadata (such as genome display names or assembly level) into the `PAN.db` database. `-t layers` tells anvi'o that you are modifying layer data. After this step, your metadata is available in the pangenome interactive interface.

```sh
anvi-import-misc-data -p 03_PAN/Nanosynbacter_pangenome_all-PAN.db -t layers names.txt
```

## 4) Display the pangenome

```sh
anvi-display-pan -g Nanosynbacter_pangenome_all-GENOMES.db -p Nanosynbacter_pangenome_all-PAN.db
```
