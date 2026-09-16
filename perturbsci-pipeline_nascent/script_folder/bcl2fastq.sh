#!/bin/bash
#
#SBATCH --job-name=mkfastq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0:30:00
#SBATCH --mem=16G
#SBATCH -o jobs/mkfastq.%j.out 
#SBATCH -e jobs/mkfastq.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=imascio@nygenome.org



#########
## This is an example of the bcl2fastq script you will run on the NYGC cluster
## to get the fastq files demultiplexed based on the PCR barcodes. 
## You run this script on the command line as `sbatch bcl2fastq.sh`.
## you need to change the following based on your own directories/sequencing:
## "dir" -- the directory you are working in aka where the fastqs will be saved
## to (inside the fastq folder that is made).
## "--runfolder-dir" -- the argument in bcl2fastq that points to the sequencing
## run. This directory name can be found in the basespace run info. It also has
## the date the sequencing run started as the beginning.
## You will also need a sample sheet that contains information on which PCR
## barcodes were used (both for rna and guide i7 barcodes). An example is saved
## in this script folder. I always name mine sample_sheet.csv but if you use a 
## different name you need to change the "--sample-sheet" argument in the 
## bcl2fastq command.
#########
module purge
module load bcl2fastq/2.20.0.422
module load multiqc/1.8

dir=/gpfs/commons/groups/satija_lab/imascio/RBP_screen/241022_sequencing_71APA_redo

fastq=${dir}/fastq

mkdir ${fastq}

bcl2fastq --runfolder-dir /gpfs/commons/instruments/nextseq/NB552173/241022_NB552173_0394_AHYNF7BGXW \
	-o ${fastq} \
	--sample-sheet ${dir}/sample_sheet.csv \
	--reports-dir ${fastq}/report \
	--barcode-mismatches 1 \
	--create-fastq-for-index-reads \
	--no-lane-splitting \
	--use-bases-mask Y*,I*,Y*,Y* \
	--minimum-trimmed-read-length 0 \
	--mask-short-adapter-reads 0

multiqc ${fastq}
