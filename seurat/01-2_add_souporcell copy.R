#/usr/bin/env R
#-----------------------------------------------------------------------
# description : SouporCell labeling
# author      : songnsnow
# date        : 231229
# notes       : 
#-----------------------------------------------------------------------
# INPUT ##########################################################
# batch 1
setwd("/data/project/RCC_PBMC_HWS/workflow/singlecell/seurat/b1")   # set working directory
getwd()     # check
dt <- readRDS('rds/01-1.rds')  # batch SeuratObject file
soup_path <- '/data/project/RCC_PBMC_HWS/workflow/singlecell/souporcell/outs_b1/clusters.tsv'   # SouporCell results PATH
fn <- 'b1'  # filename

# batch 2 (after cb)
setwd("/data/project/RCC_PBMC_HWS/workflow/singlecell/seurat/b2")   # set working directory
getwd()     # check
dt <- readRDS('rds/01-1.rds')  # batch SeuratObject file
soup_path <- '/data/project/RCC_PBMC_HWS/workflow/singlecell/souporcell/outs_b2/clusters.tsv'   # SouporCell results PATH
fn <- 'b2'  # filename

# batch 2 (before cb)
setwd("/data/project/RCC_PBMC_HWS/workflow/singlecell/seurat/b2")   # set working directory
getwd()     # check
dt <- readRDS('rds/01-1.rds')  # batch SeuratObject file
soup_path <- '/data/project/RCC_PBMC_HWS/workflow/singlecell/souporcell/outs_b2/clusters.tsv'   # SouporCell results PATH
fn <- 'b2'  # filename
##################################################################

soup_csv <- read.csv(soup_path, header=TRUE, sep='')
soup_df <- as.data.frame(soup_csv)
soup_ass <- data.frame(soup_df[,1:3])   # only get assigned clusters info by cell barcode

colnames(soup_ass) <- c('cell','status', 'cluster')
rownames(soup_ass) <- soup_ass$cell
soup_ass$cell <- NULL
dt <- AddMetaData(dt, metadata=soup_ass)

df <- as.data.frame(dt@meta.data)
df[ ,c('orig.ident', 'nCount_RNA','nFeature_RNA')] <- list(NULL)
df$cluster <- as.character(df$cluster)

# --------------------------------------------------------------------------------------------------
# Assign RCC patient ids - b1
for (i in 1:nrow(df)){
  if (df[i,"cluster"] == '0'){
    donor.id <- 'RCC.4'
  }
  else if (df[i,"cluster"] == '1'){
    donor.id <- 'RCC.8'
  }
  else if (df[i,"cluster"] == '2'){
    donor.id <-'RCC.1'
  }
  else if (df[i,"cluster"] == '3'){
    donor.id <- 'RCC.3'
  }
  else if (df[i,"cluster"] == '4'){
    donor.id <- 'RCC.7'
  }
  else if (df[i,"cluster"] == '5'){
    donor.id <- 'RCC.6'
  }
  else if (df[i,"cluster"] == '7'){
    donor.id <- 'RCC.5'
  }
  else {
    donor.id <- 'NA'
  }
  df[i,'donor.id'] <- donor.id
}

# Assign RCC patient ids - b2 (after cb)
for (i in 1:nrow(df)){
  if (df[i,"cluster"] == '0'){
    donor.id <- 'RCC.8'
  }
  else if (df[i,"cluster"] == '1'){
    donor.id <- 'RCC.9'
  }
  else if (df[i,"cluster"] == '2'){
    donor.id <-'HC.1'
  }
  else if (df[i,"cluster"] == '4'){
    donor.id <- 'RCC.6'
  }
  else if (df[i,"cluster"] == '5'){
    donor.id <- 'HC.2'
  }
  else if (df[i,"cluster"] == '6'){
    donor.id <- 'RCC.2'
  }
  else {
    donor.id <- 'NA'
  }
  df[i,'donor.id'] <- donor.id
}



# Assign RCC patient ids - b2 (before cb)
for (i in 1:nrow(df)){
  if (df[i,"cluster"] == '0'){
    donor.id <- 'RCC.8'
  }
  else if (df[i,"cluster"] == '1'){
    donor.id <- 'RCC.9'
  }
  else if (df[i,"cluster"] == '2'){
    donor.id <-'HC.1'
  }
  else if (df[i,"cluster"] == '4'){
    donor.id <- 'RCC.6'
  }
  else if (df[i,"cluster"] == '5'){
    donor.id <- 'HC.2'
  }
  else if (df[i,"cluster"] == '6'){
    donor.id <- 'RCC.2'
  }
  else {
    donor.id <- 'NA'
  }
  df[i,'donor.id'] <- donor.id
}
# --------------------------------------------------------------------------------------------------

# Add RCC subtype
NORM <- c("HC.1", "HC.2")
NAG <- c("RCC.1", "RCC.2","RCC.3", "RCC.4")
AG <- c("RCC.5" , "RCC.6", "RCC.7", "RCC.8", "RCC.9")

for (i in 1:nrow(df)) {
  if(df[i, "donor.id"] %in% NAG) {
    subtype <- "NAG"
  }
  else if (df[i, "donor.id"] %in% AG) {
    subtype <- "AG"
  }
  else if (df[i, "donor.id"] %in% NORM) {
    subtype <- "NORM"
  }
  else {
    subtype <- "NA"
  }
  df[i, "subtype"]  <- subtype
}

# add RCC subtype2
BM <- c("RCC.5", "RCC.6")
LM <- c("RCC.7", "RCC.8", "RCC.9")

for (i in 1:nrow(df)) {
  if(df[i, "donor.id"] %in% NAG) {
    subtype <- "NAG"
  }
  else if (df[i, "donor.id"] %in% LM) {
    subtype <- "LM"
  }
  else if (df[i, "donor.id"] %in% BM) {
    subtype <- "BM"
  }
  else if (df[i, "donor.id"] %in% NORM) {
    subtype <- "NORM"
  }
  else {
    subtype <- "NA"
  }
  df[i, "subtype2"]  <- subtype
}

df[ ,c('status', 'cluster')] <- list(NULL)

dt <- AddMetaData(dt, metadata=df)
columns.to.remove <- c("status", "cluster")
for(i in columns.to.remove) {
  dt[[i]] <- NULL
}

saveRDS(dt,'rds/01-2.rds')