#!/bin/bash

# Create deduplicated bam files (dedup on UMI and genomic location)

# Define variables
# input is feature output

input_folder=$1
output_folder=$2
sampleID=$3

umi_tools dedup --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${input_folder}/${sampleID}_assigned_sorted.bam -S ${output_folder}/${sampleID}_dedup.bam -L ${output_folder}/${sampleID}_log.out

# sort and index the bam files

samtools sort -o ${output_folder}/${sampleID}_dedup_sorted.bam ${output_folder}/${sampleID}_dedup.bam 
samtools index ${output_folder}/${sampleID}_dedup_sorted.bam
