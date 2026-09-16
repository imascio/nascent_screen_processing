#!/bin/bash

# directory where fastqs are saved
fastq_dir=$1

# names of the files with gex sample names and gdo names to be created
gex_output=$2
gdo_output=$3


find $1 -type f -name '*.fastq.gz' -exec basename {} \; | awk -F'_' '{print $1}' | sort -u | awk -v gex_o="$gex_output" -v gdo_o="$gdo_output" '
/gex/ { print > gex_o }
/gdo/ { print > gdo_o }
'
