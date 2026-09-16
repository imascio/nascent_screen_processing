
input_folder=$1
sample=$2
output_folder=$3

echo Trimming sample: $sample
trim_galore --paired $input_folder/$sample*R1*.gz $input_folder/$sample*R2*.gz -a2 AAAAAAAA --stringency 3 -o $output_folder --length 0
echo Trimming $sample done.
