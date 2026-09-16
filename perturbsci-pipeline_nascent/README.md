# perturbsci-pipeline
A processing pipeline that will take fastq reads produced from a PerturbSci library preparation with a NextSeq 550 75 cycle kit and output a counts matrix and provides a script to make a seurat object. This is run on Croen due to the parallelization. You need to run `bcl2fastq` yourself on NYGC clusters and move the fastqs over to our servers. There is an example of what the `bcl2fastq.sh` script should look like (with comments explaining how to run it) and an example of the associated sample sheet (`sample_sheet.csv`) required to run `bcl2fastq`.

## Create the conda environment
Create a conda environment with the required packages from the provided YAML file.
```{bash}
conda env create --name EasySci --file=EasySci_conda_env.yaml
```

## Input parameters
These parameters are input at the beginning of the `main_pipeline.sh` file and vary based on your specific experiment. Check out the parameters $feature_type and $read2_length to make sure they match your experiment. They are commented in the `main_pipeline.sh` file.

### Sample ID files
Create two files, one for RNA samples and one for guide samples. These are .txt files where each line is the sample name for the corresponding PCR barcode. Subsequently these are the sample names that are added to the fastq files after demultiplexing as they are provided in the sample sheet mentioned above. An example for 4 PCR barcodes would be:
```
rna01
rna02
rna03
rna04
```
An example of 2 guide enrichment PCR barcodes would be:
```
gdo1
gdo2
```

Additionally, rna01-rna08 are always associated with gdo1, rna09-16 are always associated with gdo2, etc. So make sure if your sample sheet the barcode sequences correspond to this. (Ask me if this is confusing.)

### Genome Reference
This pipeline uses STAR to align your reads. Build or download a reference such as the one from [10X Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest?). If you don't already have STAR indices built from the refernce, you will have to do that with:
```{bash}
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir path/to/star/output/reference --genomeFastaFiles path/to/fasta/refernce.tar.gz
```
Note that if the version of STAR that made the index files is not the same version of STAR that is aligning the data, you will get an error.

### sgRNA Correction File and Info Table
You will need to provide a .pickle2 file of the guide sequences used in your experiment. This can be created with the `gRNA_whitelist_generation.py` script provided. To run this script you will need to provide a file with sgRNA information. This is a tab delimited .txt file with three columns. The file has a header with the columns as `gRNA_seq`, `names`, and `gene_symbol`. The first column is the `gRNA_seq` or sequence for your sgRNAs. Our sequencing structure with a 75 cycle kit only capture the first 19 out of 20 bp which is reflected in this column. If you are using a different kit that captures all 20bp, you will need a new sgRNA information .txt file with 20bp instead of 19bp. There is a parameter called `read2_length` that reflects this information and is explained in the `main_pipeline.sh` file. The second column is the `names` or the name of your sgRNA (i.e. CD55g1). The third column is the `gene_symbol` or the gene symbol of the sgRNA target. Here is an example:
```
gRNA_seq        names   gene_symbol
GCGCGGCCCAGAGGCGGCA     CPSF6g1 CPSF6
GACCTGGGCTCTCACGTAC     CPSF6g2 CPSF6
GGCGGGCAGTAGCCGCTGA     NUDT21g1        NUDT21
GGCTACTGCCCGCCATTAA     NUDT21g2        NUDT21
```
Then change what your sgRNA info table is called in the `gRNA_pilot_screen_file` variable and what you want the .pickle2 file to be called with the `pickle_out` variable in the `gRNA_whitelist_generation.py` script before running it.

Both the info table .txt file and .pickle2 file are inputs for the main script file.

### Gene Name Output
The `feature_type` parameter refers to whether you want the counts matrix to have ENSEMBL IDs or gene names as the row names. This parameter must be equal to "ENSEMBL" or "gene_name", respectively, for the script to work.

## Note about FeatureCounts
The feature step uses featureCounts from subread which defaults to only looking at exon regions. This pipeline has the `-t` option set to `"gene"` to make sure all mapped regions of genes are counted (exon, UTR, intron, etc.). Look to the subread [documentation](https://subread.sourceforge.net/SubreadUsersGuide.pdf) if you wish to change this option.

## Run the Pipeline
After changing all the variables in the beginning of the `main_pipeline.sh` file to fit your experiment, you can run the pipeline. 

First, activate a tmux session (or soemthing similar so when your laptop goes to sleep the pipeline isn't killed.
Then, activate your conda environment, which should be called EasySci if you ran the code from above.

Now you can run the pipeline as:
```{bash}
bash main_pipeline.sh > output.log 2> messages.log
```
Errors and other messages will be saved to messages.log and normal output will be written to output.log.

You will get fastq files and sam files ouput at each processing step in the `gex_processing` directory. The guide counts will be output to the `gdo_processing` directory. You will get a final seurat object (see below) with the guide and gene expression counts in the `gex_processing` directory. Proceed by checking the QC metrics (nCount_RNA, nFeature_RNA, nCount_GDO) and follow Seurat clustering workkflows to analyze your data, including the `MULTIseqDemux` function to assign guides. 

### Making Seurat Object
Since this script has to connect to ENSEMBL servers, which sometimes stall out, it is very likely that it may fail. If you have output in the `gex_processing/count` and `gdo_processing` directories, just open this script up in R and run it line by line on your own, which I often have to do. Additionally, sometimes new package updates may cause dependencies to fail and the script may error out as well. So definitely don't count on this seurat object to be correct or even present when the pipeline finishes. Sorry :(
