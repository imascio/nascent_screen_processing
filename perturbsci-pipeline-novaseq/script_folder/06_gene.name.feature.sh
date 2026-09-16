#!/bin/bash
# Assign reads to genes

# Define variables

gtf_file=$1
sampleID=$2
input_folder=$3
output_folder=$4


featureCounts -a $gtf_file \
              -o ${output_folder}/${sampleID}_gene_assigned \
              -g "gene_name" \
              -t "gene" \
              -R BAM ${input_folder}/${sampleID}Aligned.sortedByCoord.out.bam \
              -T 4;            
samtools sort ${output_folder}/${sampleID}Aligned.sortedByCoord.out.bam.featureCounts.bam -o ${output_folder}/${sampleID}_assigned_sorted.bam;
samtools index ${output_folder}/${sampleID}_assigned_sorted.bam;
