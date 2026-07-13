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


