#!/bin/bash
set -euo pipefail

input_folder=$1
output_folder=$2
sampleID=$3

mkdir -p "${output_folder}"

umi_tools count \
    --per-gene \
    --gene-tag=XT \
    --assigned-status-tag=XS \
    --skip-tags-regex="Unassigned" \
    --per-cell \
    -I "${input_folder}/${sampleID}_nascent_dedup_sorted.bam" \
    -S "${output_folder}/${sampleID}_nascent_counts.tsv.gz" \
    --wide-format-cell-counts
