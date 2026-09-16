#!/bin/bash

sampleID=$1
input_dir=$2
output_dir=$3

echo Checking $sampleID
zcat ${input_dir}/$sampleID*R2*gz | sed 'N;N;N;/\n\n/d' | gzip -c > ${output_dir}/$sampleID.R2.fq.gz

echo "################################"
echo "Finished cleaning reads"
