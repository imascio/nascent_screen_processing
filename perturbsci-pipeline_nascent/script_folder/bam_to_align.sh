#!/bin/bash
set -euo pipefail

input_folder=$1
sampleID=$2
output_folder=$3
ref_genome_fa=$4
sam2tsv_cmd=$5
java_bin=${6:-java}

mkdir -p "${output_folder}"

echo "Generating sam2tsv alignment table for ${sampleID}" >&2
if [[ "${sam2tsv_cmd}" == *.jar ]]; then
    "${java_bin}" -jar "${sam2tsv_cmd}" \
        -R "${ref_genome_fa}" \
        "${input_folder}/${sampleID}_dedup_sorted.bam" \
        > "${output_folder}/${sampleID}.align"
else
    "${sam2tsv_cmd}" \
        -R "${ref_genome_fa}" \
        "${input_folder}/${sampleID}_dedup_sorted.bam" \
        > "${output_folder}/${sampleID}.align"
fi
