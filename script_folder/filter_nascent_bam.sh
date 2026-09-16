#!/bin/bash
set -euo pipefail

input_folder=$1
sampleID=$2
output_folder=$3
nascent_reads_txt_folder=$4
picard_cmd=$5
java_bin=${6:-java}

input_bam="${input_folder}/${sampleID}_dedup_sorted.bam"
input_read="${nascent_reads_txt_folder}/${sampleID}.txt"
output_bam="${output_folder}/${sampleID}_nascent_dedup.bam"

mkdir -p "${output_folder}"

if [[ ! -s "${input_read}" ]]; then
    echo "No nascent read names for ${sampleID}; writing empty BAM with header." >&2
    samtools view -H "${input_bam}" | samtools view -b -o "${output_bam}" -
else
    if [[ "${picard_cmd}" == *.jar ]]; then
        "${java_bin}" -jar "${picard_cmd}" FilterSamReads \
            I="${input_bam}" \
            O="${output_bam}" \
            READ_LIST_FILE="${input_read}" \
            FILTER=includeReadList
    else
        "${picard_cmd}" FilterSamReads \
            I="${input_bam}" \
            O="${output_bam}" \
            READ_LIST_FILE="${input_read}" \
            FILTER=includeReadList
    fi
fi

samtools sort -o "${output_folder}/${sampleID}_nascent_dedup_sorted.bam" "${output_bam}"
samtools index "${output_folder}/${sampleID}_nascent_dedup_sorted.bam"
