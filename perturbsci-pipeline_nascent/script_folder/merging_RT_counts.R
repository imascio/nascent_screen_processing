# script to add RT type and then merge the separate counts into single cells
# meant to run line by line in RStudio - not written to run as a script from the command line
######
# installing packages
# List of package names you want to install
packages_to_install <- c("Seurat", "tidyverse", "patchwork", "scales", "reshape2", "SeuratData", "mixtools", "readxl", "data.table")

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

######
# load barcode table
all.bc <- read.csv("script_folder/rt_lig_ligrc_sgrnacapt_barcodes.csv")
######
# add the type of RT barcode to a maetadata column
obj$rt.type <- NA

for (i in 1:length(obj$rt.type)) {
  if (obj$rt.barcode[i] %in% c(all.bc$shortdT1,all.bc$shortdT2)) {
    obj$rt.type[i] <- "shortdT"
  } else {
    if (obj$rt.barcode[i] %in% c(all.bc$randN1, all.bc$randN2)) {
      obj$rt.type[i] <- "randHex"
    } else {
      obj$rt.type[i] <- "N/A"
    }
  }
}
######
# subset obj into shortdT and random hexamer only
shortdT_obj <- subset(obj, subset = rt.type == "shortdT")

randHex_obj <- subset(obj, subset = rt.type == "randHex")


sum(!(randHex_obj$rt.barcode %in% c(all.bc$randN1,all.bc$randN2)))
sum(!(shortdT_obj$rt.barcode %in% c(all.bc$shortdT1,all.bc$shortdT2)))

randHex_obj$rt.match <- NA

for (i in 1:length(randHex_obj$rt.barcode)) {
  if(i %% 1000 == 0) {print(i)}
  if (randHex_obj$rt.barcode[i] %in% all.bc$randN1) {
    randHex_obj$rt.match[i] <- all.bc$shortdT1[randHex_obj$rt.barcode[i] == all.bc$randN1]
  }
  if (randHex_obj$rt.barcode[i] %in% all.bc$randN2) {
    randHex_obj$rt.match[i] <- all.bc$shortdT2[randHex_obj$rt.barcode[i] == all.bc$randN2]
  }
}

head(colnames(randHex_obj))
head(randHex_obj$lig.barcode)

new.randHex.colnames <- paste0(randHex_obj$pcr_barcode,"_",randHex_obj$lig.barcode,randHex_obj$rt.match)

randHex_obj_matched <- randHex_obj
colnames(randHex_obj_matched) <- new.randHex.colnames

head(randHex_obj_matched$rt.barcode)
head(randHex_obj_matched$rt.match)
ncol(randHex_obj_matched)
######
# adding metadata column for if the cell has a matching randHex
obj$has.randH <- NA
for(i in 1:length(colnames(obj))) {
  if (colnames(obj)[i] %in% new.randHex.colnames) {
    obj$has.randH[i] <- "Y"
  } else {
    obj$has.randH[i] <- "N"
  }
}

table(obj$has.randH)

shortdT_obj$has.randH <- NA
for(i in 1:length(colnames(shortdT_obj))) {
  if (colnames(shortdT_obj)[i] %in% new.randHex.colnames) {
    shortdT_obj$has.randH[i] <- "Y"
  } else {
    shortdT_obj$has.randH[i] <- "N"
  }
}

table(shortdT_obj$has.randH)




shortdT.RNA <- shortdT_obj[["RNA"]]$counts
shortdT.GDO <- shortdT_obj[["GDO"]]$counts


randHex.RNA <- randHex_obj_matched[["RNA"]]$counts
randHex.GDO <- randHex_obj_matched[["GDO"]]$counts

randhex.unique <- colnames(randHex.RNA)[!(colnames(randHex.RNA) %in% colnames(shortdT.RNA))]
colSums(randHex.RNA[,randhex.unique])

# creating merged rna counts matrix
new.rna <- rbindlist(lapply(list(shortdT.RNA, randHex.RNA), function(x) setDT(as.data.frame(x), 
                                               keep.rownames = TRUE)), fill = TRUE)[, lapply(.SD, sum, na.rm = TRUE), by = rn]
# converting to data frame from data table to remove the rn column and make it the rownames
new.rna.df <- as.data.frame(new.rna)
rownames(new.rna.df) <- new.rna.df$rn
new.rna.df <- dplyr::select(new.rna.df, -c("rn"))

# doing the same with guide counts matrix
new.gdo <- rbindlist(lapply(list(shortdT.GDO, randHex.GDO), function(x) setDT(as.data.frame(x), 
                                               keep.rownames = TRUE)), fill = TRUE)[, lapply(.SD, sum, na.rm = TRUE), by = rn]
new.gdo.df <- as.data.frame(new.gdo)
rownames(new.gdo.df) <- new.gdo.df$rn
new.gdo.df <- dplyr::select(new.gdo.df, -c("rn"))


# make the seurat object
merged.obj <- CreateSeuratObject(counts = new.rna.df, assay = "RNA")
merged.obj[["GDO"]] <- CreateAssayObject(counts = new.gdo.df)



# adding the meta data
merged.obj$pcr.barcode <- gsub("_[A-Z]+","",colnames(merged.obj))

cells <- colnames(merged.obj)
rtlig.barcodes <- str_split_fixed(cells, "[0-9]_", n=2)[,2]
ligation.barcodes <-  substr( rtlig.barcodes , start = 1, stop = 10)
rt.barcodes <- substr( rtlig.barcodes , start = 11 , stop = 20)

merged.obj$rt.barcode <- rt.barcodes
merged.obj$lig.barcode <- ligation.barcodes

merged.obj$rt.type <- NA
for (i in 1:length(merged.obj$rt.type)) {
  if (colnames(merged.obj)[i] %in% colnames(randHex.RNA) & colnames(merged.obj)[i] %in% colnames(shortdT.RNA)){
    merged.obj$rt.type[i] <- "both"
  } else {
    if (colnames(merged.obj)[i] %in% colnames(shortdT.RNA)) {
      merged.obj$rt.type[i] <- "shortdT"
    } else {
      if (colnames(merged.obj)[i] %in% colnames(randHex.RNA)) {
        merged.obj$rt.type[i] <- "randHex"
      }
    }
  }
}
table(merged.obj$rt.type)


merged.obj$percent.mt <- PercentageFeatureSet(merged.obj, pattern = "^MT-")
saveRDS(merged.obj, "seurat_objects/240523_full_object_rhex_dt_merged_noguides_assigned.rds")
