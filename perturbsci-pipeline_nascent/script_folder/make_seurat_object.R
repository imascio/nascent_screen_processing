### January 17, 2025 ###
# List of package names you want to install
packages_to_install <- c("Seurat", "SeuratObject", "tidyverse", "data.table", "Azimuth", "genio", "Matrix")

# Check and install packages if not already installed
for (package_name in packages_to_install) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name)
  }
}

# Load the packages
for (package_name in packages_to_install) {
  library(package_name, character.only = TRUE)
}

rm(package_name, packages_to_install)

## read in arguments following terminal command
# order of arguments 1)gex_sampleID.txt 2)input_dir_counts_matrix 3)output_dir 4)output_name.rds 5)gdo_processing_dir 6)script_folder 7)gdo_sampleID.txt
args <- commandArgs(trailingOnly = TRUE)

############################GEX PROCESSING####################################################

### make gex matrix list and rename col names
gex.list <- list()
gex.seuratv5.list <- list()
n_pcr_samples <- count_lines(args[1])
sample_names <- read.table(args[1])
gdo_names <-read.table(args[7])

# if converting to ensemble times out you will have to start the for loop over again
# but can change the starting number from 1 to pick up on the last object made
# also there are lines commented out that can be used instead of ENSEMBL but I have come across issues with it which is why I prefer ENSEMBL
for (i in 1:n_pcr_samples) {
  # read in counts
  name <- sample_names$V1[i]
  print(name)
  c <- fread(file = 
               paste0(args[2],"/",name, "_counts.tsv.gz"))
  # save df without gene column and set gene column as rownames
  rownames(c) <- c$gene
  m <- as.matrix(c[,2:ncol(c)], nrow = nrow(c))
  #genes.ens <- c$gene
  rownames(m) <- c$gene
  #genes.ens.short <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", genes.ens)
  #ids <- mapIds(org.Hs.eg.db,
                #keys=genes.ens,
                #column="SYMBOL",
                #keytype="ENSEMBL",
                #multiVals="first")
  #rownames(m) <-ids
  #n <- Azimuth:::ConvertEnsembleToSymbol(m, species = "human")
  n <- m
  gex.list[[i]] <- n
  # saving seurat objects
  obj <- CreateSeuratObject(counts = n, min.features = 200)
  obj$pcr_barcode <- name
  gex.seuratv5.list[[i]] <- obj
  print("done")
}

# combining all seurat objects together
if (length(gex.seuratv5.list) > 1) {
  objv5 <- merge(gex.seuratv5.list[[1]],
                 gex.seuratv5.list[2:length(gex.seuratv5.list)],
                 add.cell.ids = sample_names$V1)
} else {
  objv5 <- gex.seuratv5.list[[1]]
}

objv5$sample <- gsub("_[0-9]","",objv5$pcr_barcode)


cells <- colnames(objv5)
rtlig.barcodes <- str_split_fixed(cells, "[0-9]_", n=2)[,2]
ligation.barcodes <-  substr(rtlig.barcodes, start = 1, stop = 10)
rt.barcodes <- substr(rtlig.barcodes, start = 11 , stop = 20)

objv5$rt.barcode <- rt.barcodes
objv5$lig.barcode <- ligation.barcodes

# save seurat object that is only gex counts right now
saveRDS(objv5,paste0(args[3],"/gex_processing.rds"))
##############################################################################################

############################GDO PROCESSING####################################################
# prepping guide matrix  - Zihan from Cao Lab gRNA_MM_formatting.R
sciRNAseq_gene_count_summary <- function (gene_count_folder) {
  gene_matrix_dir = paste(gene_count_folder, "/gRNA_mat.count", sep = "")
  df_gene_dir = paste(gene_count_folder, "/gRNA_annotat.report", sep = "")
  df_cell_dir = paste(gene_count_folder, "/cell_gRNA_annotat.report", sep = "")
  
  df_gene = read.csv(df_gene_dir, header = F)
  df_cell = read.csv(df_cell_dir, header = F)
  gene_matrix = read.csv(gene_matrix_dir, header = F)
  colnames(df_gene) = c("gRNA_index", "gRNA_name")
  colnames(df_cell) = c("cell_index", "cell_name")
  colnames(gene_matrix) = c("gRNA_index", "cell_index", "count")
  rownames(df_gene) = df_gene$gene_id
  rownames(df_cell) = df_cell$cell_name
  gene_count = sparseMatrix(i = gene_matrix$gRNA_index, j = gene_matrix$cell_index, x = gene_matrix$count)
  df_gene = df_gene[1:nrow(gene_count), ]
  rownames(gene_count) = df_gene$gRNA_name
  colnames(gene_count) = df_cell$cell_name
  
  return(list(df_cell, df_gene, gene_count))
}

# save location of report folder
report_folder = paste0(args[5],"/gRNA_report")

# format the output into a sparse matrix ready for R analysis
result = sciRNAseq_gene_count_summary(report_folder)
gRNA_df_cell = result[[1]]
gRNA_df_gene = result[[2]]
gRNA_count = result[[3]]
# saving guide data as RData
save(gRNA_df_cell, gRNA_df_gene, gRNA_count, file = paste0(args[3], "/gRNA_Summary.RData"))

#saving sparse matrix to matrix, duplicating the single row and making new row names
gRNA_mat <- as.matrix(gRNA_count)
colnames_gdo <- colnames(gRNA_mat)
# making a table of the combinations of barcodes in order of colnames
gdo_bc_combinations <- data.frame(inneri7 = gsub("\\..*","",colnames_gdo),
                                  sgRNAcapture = substr(gsub(".*\\.","",colnames_gdo) , start = 11 , stop = 20),
                                  lig = substr(gsub(".*\\.","",colnames_gdo) , start = 1 , stop = 10))

# making inner i7 list - need both gex and gdo sample info to correctly combine the pcr rna barcodes
gdo.df <- expand.grid(c("TCGGATTCGG", "CTAAGCCTTG", "CTAACTAGGT", "GCAAGACCGT", "ATGGAACGAA", "TAGAGGCGTT", "GCATCGTATG", "TGGACGACTA"),gdo_names[,1])
sample_inner_combos <- paste0(gdo.df$Var2, gdo.df$Var1)
i7list <- data.frame(i7_ID = sample_names[,1],
                     inneri7 = sample_inner_combos)
gdo_bc_i7ID <- left_join(x = gdo_bc_combinations, y = i7list)

# if your shortdT RT plates do not have matching sgRNA barcodes you can use this commented out code
# to match the sgRNA capture barcode with the barcode of the corresponding shortdT barcode plate (in this example it's plate 2)
# gdo_bc_i7ID$RT <- NA

# all.bc <- read.csv(paste0(args[6],"/rt_lig_ligrc_sgrnacapt_barcodes.csv"))

# for (i in 1: nrow(gdo_bc_i7ID)) {
#   gdo_bc_i7ID$RT[i] <- all.bc$shortdT2[gdo_bc_i7ID$sgRNAcapture[i] == all.bc$sgRNAcapture]
# }

# the new plates have sgRNA capture match the RT shortdT barcodes so RT and sgRNA capture sequences are the same 
gdo_bc_i7ID$RT <- gdo_bc_i7ID$sgRNAcapture

new_gd_colnames <- paste0(gdo_bc_i7ID$i7_ID, "_" ,gdo_bc_i7ID$lig,gdo_bc_i7ID$RT)
colnames(gRNA_mat) <- new_gd_colnames


# save processed guide matrix
saveRDS(gRNA_mat, paste0(args[3],"/gdo_procesed_matrix.rds"))
#################################################################################################

##############################MERGING GEX AND GDO INTO ONE OBJECT################################
### find cells (col names) that intersect with gex and gdo & make new gdo matrix with those
cell.intersect <- intersect(colnames(objv5), colnames(gRNA_mat))

### find cells (col names) that are in gex but missing in gdo
missing_cells <- colnames(objv5)[!(colnames(objv5) %in% colnames(gRNA_mat))]

### make matrix with 0 for gdo counts of missing cells then cbind it to intersected gdo matrix
missing_guides <- matrix(data = 0, nrow = nrow(gRNA_mat), ncol = length(missing_cells))
colnames(missing_guides) <- missing_cells
rownames(missing_guides) <- rownames(gRNA_mat)

gdo.mat <- cbind(gRNA_mat, missing_guides)

# subset gdo.mat to only have cells in gex obj
all(colnames(objv5) %in% colnames(gdo.mat))

gdo.subset <- gdo.mat[,colnames(objv5)]

### merge gex seurat objects and add gdo assay to it
obj <- objv5
obj[["GDO"]] <- CreateAssayObject(counts = gdo.subset)
#################################################################################################

##############################ADD GDO META#############################################
## adding guide metadata - you can skip this section if you want
# making dataframe with cell names, guide counts, and a binning category
guide.meta <- data.frame(cells = names(obj$nCount_GDO), counts = obj$nCount_GDO, guide_group = NA)
# guides less than 10
guide.meta[guide.meta$counts < 10, 3] <- "<10"
# guides 10 to less than 50
guide.meta[guide.meta$counts >= 10 & guide.meta$counts < 50, 3] <- "10-50"
# guides 50 to less than 100
guide.meta[guide.meta$counts >= 50 & guide.meta$counts < 100, 3] <- "50-100"
# guides 100 or more
guide.meta[guide.meta$counts >= 100, 3] <- "100+"


guide.meta2 <- data.frame(sgRNA_count = guide.meta$guide_group, row.names = guide.meta$cells)

obj <- AddMetaData(obj, metadata = guide.meta2)

#######################################################################################

#######################################save the final object!################################################
saveRDS(obj, paste0(args[3],"/",args[4]))
#######################################################################################
