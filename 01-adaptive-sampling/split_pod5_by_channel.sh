#!/bin/bash
# Split pod5 files into channels 1-256 (adaptive sampling, depletion mode) and 257-512.

mkdir -p summed-reads

# Generate a table containing read IDs and channels:
pod5 view *.pod5 \
 --include "read_id,channel" \
 --output reads.tsv

# Create two read ID lists:
awk -F'\t' '$2 >= 1 && $2 <= 256 {print $1}' reads.tsv > channels_1_256.txt
awk -F'\t' '$2 >= 257 && $2 <= 512 {print $1}' reads.tsv > channels_257_512.txt

# Extract the reads into separate POD5 files
pod5 filter *.pod5 \
 --ids channels_1_256.txt \
 --output summed-reads/channels_1_256.pod5

pod5 filter *.pod5 \
 --ids channels_257_512.txt \
 --output summed-reads/channels_257_512.pod5
