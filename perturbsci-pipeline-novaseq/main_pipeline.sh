#!/bin/bash

#run as bash main_pipeline.sh > output.log 2> messages.log

## refer to README as this code rarely works - activate the conda environment prior to running the script
## Location of the bashrc file to activate the conda environments.
## activate_location="~/.bashrc"
## source $activate_location
## conda activate EasySci

# Define script directory
script_path="/brahms/shared/EasySci/shell_pipeline/script_folder"

###################### DEFINE INPUT VARIABLES ######################

# define the project folder, this is where all output files will be written
project_folder="/brahms/shared/EasySci/test_rna1"

# define the fastq folder including all fastq files
fastq_folder=$project_folder/fastq

# define the PCR group sample id for each gene experssion fastq file
sample_ID=$project_folder/gex_sampleID.txt

# define the PCR group sample id for each sgRNA fastq file
gdo_sample_ID=$project_folder/gdo_sampleID.txt

# define the output folder where all processing files for gene expression will be written, inside the project folder
gex_processing_folder=$project_folder/gex_processing

# define the output folder where all processing files for sgRNA will be written, inside the project folder
gdo_processing_folder=$project_folder/gdo_processing

# define the length of read 2, which can only be a numeric value of 54 or 55 - the sgRNA counting script will not work if the read 2 length is not 54 or 55
read2_length=55

# feature type will set the feature.sh script to output counts using ENSEMBL IDs or gene names - this parameter must equal "gene_name" or "ENSEMBL" or the script won't run
feature_type="gene_name"

# define the number of gRNA UMI cutoff for keeping cells
cutoff=10

# define the core number for parallel processing
core=16

# define the core number for reads filtering and sorting
samtools_core=16

# define the number of jobs for parallelizing
N_JOBS=16
# STAR uses a lot of memory and cannot run as many jobs when you have large data from novaseq
N_JOBS_STAR=10

# define the location of index files for reads alignment with STAR
# Human genome version 42
index="/brahms/shared/refdata-gex-GRCh38-2024-A/STAR_index_2.7.10b"

# define the gtf file for gene counting
# Human genome version 44 - downloaded from 10X 2024 cellranger build
gtf_file="/brahms/shared/refdata-gex-GRCh38-2024-A/genes/genes.gtf"

# define the gtf file for exon counting
# Human genome version 42
gtf_file_exon="/brahms/shared/EasySci/human/gencode.v42.primary_assembly.annotation.exons.gtf"

# define the location for STAR to make the tmp files - this can be your working directory unless FIFO files cannot be written to it
tmp_folder="/home/mascioi"

# define the location of the ligation barcodes
ligation_barcode=$script_path/Ligation_barcodes_NextSeq.pickle2
# define the location of the RT barcodes
RT_barcode=$script_path/RT_barcodes.pickle2
# define the location of the randomN RT barcodes for the gene counting step
randomN_barcode_file=$script_path/RandomN_RT_barcodes.txt
#define the folder of inner i7 barcode dictionary
inner_i7_bc_file=$script_path/simp_inner_i7_220517.pickle2

#define the folder of gRNA barcode dictionary - this is experiment-specific
gRNA_correction_file=/brahms/mascioi/RBP_screen/RBP_20bp_sgRNA_seq.pickle2
#define the folder containing the gRNA annotation file - this is experiment-specific 
gRNA_annotation_df=/brahms/mascioi/RBP_screen/RBP_20bp_guide_info_table.txt

# table with two columns - "GEX" which is the gex_sampleID and "inneri7" which is the corresponding gdo_sampmleID + inneri7
# this is to rename the gdo matrix to have matching cell names as the gex counts matrix
match_table="/brahms/mascioi/RBP_screen/novaseq/full.gex.gdo.colnames.matching.for.making.seurat.object.csv"
# define the name of the final seuart object - should be .rds extension
seurat_object_name=pipeline_seurat_output.rds


######### NEW STUFF ##########################

############ CREATING SAMPLE SHEET FROM FASTQS ############

# 

input_folder=$fastq_folder
script=$script_path/01_make_sample_ID_files.sh

now=$(date +"%T")
echo "~~~~~~~~~~~~~~~~~~" >&2
echo "01. Collecting samples names from gex and gdo fastq files... $now" >&2

bash $script $input_folder $sample_ID $gdo_sample_ID

now=$(date +"%T")
echo "Sample sheets made. $now" >&2

######### END NEW STUFF ##########################
############ BARCODE EXTRACTION ############

# the script take an input folder, a sample ID list, an output folder, the RT barcode list, the ligation barcode list and core number and a rendom hexamer barcode list. Then it extract the RT barcode from read1, the ligation barocde from read2, correct them to the nearest RT and ligation barcode (with edit distance <= 1), and attach the RT and ligation barcode and UMI sequence to the read name of read1 and read3. Reads with unmatched RT or ligation barcodes are discarded.

input_folder=$fastq_folder
output_folder=$gex_processing_folder/BC_attach
script=$script_path/02_barcode_extraction_parallel.py

now=$(date +"%T")
echo "~~~~~~~~~~~~~~~~~~" >&2
echo "01.5. Changing the name of the gex fastq files... $now" >&2

echo $sample_ID >&2
output_fastq=$project_folder/lane_combined_fastq
mkdir -p $output_fastq

# NextSeq outputs R1 reads as R1 fastq, I1 reads as I1 fastq, I2 reads as R2 fastq, and R2 reads and R3 fastq
# NovaSeq outputs fastqs with the same name as the reads. So I am both combining lanes and changing fastq names to match NextSeq naming scheme and won't have to change later scripts
for sample in $(cat $sample_ID); do 
    echo "changing name $sample"
    cat $input_folder/*$sample*R1*.fastq.gz > $output_fastq/$sample.R1.fastq.gz
    cat $input_folder/*$sample*I2*.fastq.gz > $output_fastq/$sample.R2.fastq.gz
    cat $input_folder/*$sample*R2*.fastq.gz > $output_fastq/$sample.R3.fastq.gz
    cat $input_folder/*$sample*I1*.fastq.gz > $output_fastq/$sample.I1.fastq.gz
done

now=$(date +"%T")
echo "~~~~~~~~~~~~~~~~~~" >&2
echo "02. Attaching barcode and UMI.... $now" >&2

input_folder=$output_fastq
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose python $script $input_folder {} $output_folder $ligation_barcode $RT_barcode $core $randomN_barcode_file :::: $sample_ID

now=$(date +"%T")
echo "Barcode transformed and UMI attached. $now" >&2

################# TRIM READ 2 #################

# This steps trims any poly A stretches present in read 2
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "03. Trimming read 2... $now" >&2

input_folder=$gex_processing_folder/BC_attach
output_folder=$gex_processing_folder/trimmed_fastq
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash $script_path/03_trimming.sh $input_folder {} $output_folder :::: ${sample_ID}

################# REMOVE EMPTY LINES #################

# This step removes any empty fastq lines resulting from trimming reads that are entirely poly A
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "04. Removing empty lines from fastqs... $now" >&2

input_folder=$gex_processing_folder/trimmed_fastq
output_folder=$gex_processing_folder/cleaned_fastq
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash ${script_path}/04_remove_empty_lines.sh {} $input_folder ${output_folder} :::: ${sample_ID}

now=$(date +"%T")
echo "Empty lines removed. $now" >&2

############ ALIGN ############
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "05. Aligning cleaned fastq files... $now" >&2

input_folder=$gex_processing_folder/cleaned_fastq
output_folder=$gex_processing_folder/STAR_alignment
mkdir -p $output_folder

parallel -j ${N_JOBS_STAR} --verbose bash $script_path/05_align.sh $index {} $input_folder $output_folder $tmp_folder :::: ${sample_ID}

now=$(date +"%T")
echo "Alignment complete. $now" >&2

############ FEATURE ############

# This step maps aligned reads to genes using the umi tools feature function
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "06. Mapping reads to genes... $now" >&2

input_folder=$gex_processing_folder/STAR_alignment
output_folder=$gex_processing_folder/feature
mkdir -p $output_folder

# Check the value of $feature_type and set $feature_script accordingly
if [ "$feature_type" == "gene_name" ]; then
    feature_script="06_gene.name.feature.sh"
elif [ "$feature_type" == "ENSEMBL" ]; then
    feature_script="06_ENSEMBL.feature.sh"
else
    echo "Error: feature_type must be either 'gene_name' or 'ENSEMBL'"
    exit 1
fi

# Print the selected script to verify
echo "The selected feature script is: $feature_script"

parallel -j ${N_JOBS} --verbose bash ${script_path}/${feature_script} $gtf_file {} $input_folder $output_folder :::: ${sample_ID}

now=$(date +"%T")
echo "Mapping to genes complete. $now" >&2

############ DEDUP BAMS ############

# Remove UMI duplicates and generate a cell x gene counts matrix
# echo "~~~~~~~~~~~~~~~~~~" >&2
# now=$(date +"%T")
# echo "06.5. Generating a counts matrix! Almost done w gex processing... $now" >&2

# input_folder=$gex_processing_folder/feature
# output_folder=$gex_processing_folder/dedup
# mkdir -p $output_folder

# parallel -j ${N_JOBS} --verbose bash $script_path/06.5_dedup.sh $input_folder $output_folder {} :::: ${sample_ID}

############ COUNT ############

# Remove UMI duplicates and generate a cell x gene counts matrix
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "07. Generating a counts matrix! Almost done w gex processing... $now" >&2

input_folder=$gex_processing_folder/feature
output_folder=$gex_processing_folder/count
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash $script_path/07_count.sh $input_folder $output_folder {} :::: ${sample_ID}


############ GDO ############

mkdir -p $gdo_processing_folder

# Change the file names of raw gdo fastq.gz
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "07.5. Changing the name of the gdo fastq files... $now" >&2

for sample in $(cat $gdo_sample_ID); do 
    echo "changing name $sample"
    cat $fastq_folder/*$sample*R1*.fastq.gz > $output_fastq/$sample.R1.fastq.gz
    cat $fastq_folder/*$sample*I2*.fastq.gz > $output_fastq/$sample.R2.fastq.gz
    cat $fastq_folder/*$sample*R2*.fastq.gz > $output_fastq/$sample.R3.fastq.gz
    cat $fastq_folder/*$sample*I1*.fastq.gz > $output_fastq/$sample.I1.fastq.gz
done
# Run the guide counting script
echo "~~~~~~~~~~~~~~~~~~" >&2
now=$(date +"%T")
echo "08. Processing the gdo reads into single-cell counts matrix (must be reformatted for seurat later)... $now" >&2

## setting the guide counting script depending on the length of read 2
if [ $read2_length -eq 55 ]; then
    guide_script="08_sgrna_20bp_count_lig_in_header.py"
elif [ $read2_length -eq 54 ]; then
    guide_script="08_sgrna_19bp_count.py"
else
    echo "Error: read2_length must be either 54 or 55"
    exit 1
fi
# Print the script variable to verify
echo "The selected script is: $guide_script"

python3 ${script_path}/${guide_script} $output_fastq ${gdo_sample_ID} $gdo_processing_folder $RT_barcode $inner_i7_bc_file $ligation_barcode $gRNA_correction_file $gRNA_annotation_df $cutoff $core


############ SEURAT ############
echo "~~~~~~~~~~~~~~~~~~" >&2
# Make a seurat object from the counts matrices

now=$(date +"%T")
echo "09. Generating a seurat object! Last step... $now" >&2

input_folder=$gex_processing_folder/count
output_folder=$gex_processing_folder/seurat
mkdir -p $output_folder

Rscript $script_path/09_make_seurat_object.R $sample_ID $input_folder $output_folder $seurat_object_name $gdo_processing_folder $script_path $gdo_sample_ID $match_table

now=$(date +"%T")
echo "You're done - good luck! $now" >&2
echo "~~~~~~~~~~~~~~~~~~" >&2
