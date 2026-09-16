#!/bin/bash

# Create a counts matrix


# Define variables

input_folder=$1
output_folder=$2
sampleID=$3

umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --skip-tags-regex="Unassigned" --per-cell -I ${input_folder}/${sampleID}_assigned_sorted.bam -S ${output_folder}/${sampleID}_counts.tsv.gz --wide-format-cell-counts

