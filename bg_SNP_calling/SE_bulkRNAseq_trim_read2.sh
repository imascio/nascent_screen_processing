
input_folder=$1
sample=$2
output_folder=$3

echo Trimming sample: $sample

shopt -s nullglob
read2_files=("$input_folder"/"$sample"*.fastq.gz "$input_folder"/"$sample"*.fq.gz)
shopt -u nullglob

if [ "${#read2_files[@]}" -ne 1 ]; then
    echo "Expected exactly one read 2 FASTQ for sample '$sample' in $input_folder, found ${#read2_files[@]}." >&2
    printf 'Matches:\n' >&2
    printf '  %s\n' "${read2_files[@]}" >&2
    exit 1
fi

trim_galore "${read2_files[0]}" -o "$output_folder"
echo Trimming $sample done.
