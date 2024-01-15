#/usr/bin/env R
#-----------------------------------------------------------------------
# description : souporcell labeling
# author      : rlo
# date        : 231223
# notes       : conda activate sc_pathway_r430
#-----------------------------------------------------------------------
#########################
## Add soup results_b1###
#########################
soup_csv <- read.csv('/data/project/RCC_PBMC_HWS/workflow/singlecell/souporcell/outs_b1/clusters.tsv', header=TRUE, sep='')
soup_df <- as.data.frame(soup_csv)
donor_id <- data.frame(soup_df[,1:3])

colnames(donor_id) <- c('cell','doublet', 'id')
rownames(donor_id) <- donor_id$cell
donor_id$cell <- NULL
dt.1 <- AddMetaData(dt.1, metadata=donor_id)

df <- as.data.frame(dt.1@meta.data)
df$orig.ident <- NULL
df$nCount_RNA <- NULL
df$nFeature_RNA <- NULL
df$id <- as.character(df$id)

# Assign RCC patient ids
for (i in 1:nrow(df)){
  if (df[i,"id"] == '0'){
    donor.id <- 'RCC.4'
  }
  else if (df[i,"id"] == '1'){
    donor.id <- 'RCC.8'
  }
  else if (df[i,"id"] == '2'){
    donor.id <-'RCC.1'
  }
  else if (df[i,"id"] == '3'){
    donor.id <- 'RCC.3'
  }
  else if (df[i,"id"] == '4'){
    donor.id <- 'RCC.7'
  }
  else if (df[i,"id"] == '5'){
    donor.id <- 'RCC.6'
  }
  else if (df[i,"id"] == '7'){
    donor.id <- 'RCC.5'
  }
  else {
    donor.id <- 'unassinged'
  }
  df[i,'donor.id'] <- donor.id
}

# Add RCC subtype
NAG <- c("RCC.1","RCC.3","RCC.4")
AG <- c("RCC.5" , "RCC.6", "RCC.7","RCC.8")

# NAG vs AG
for (i in 1:nrow(df)) {
  if(df[i, "donor.id"] %in% NAG) {
    subtype <- "NAG"
  }
  else if (df[i, "donor.id"] %in% AG) {
    subtype <- "AG"
  }
  else {
    subtype <- "NA"
  }
  df[i, "subtype"]  <- subtype
}

# NAG vs (AG vs BM)
NAG <- c("RCC.1","RCC.3","RCC.4")
BM <- c("RCC.5","RCC.6")
AG <- c("RCC.7" , "RCC.8")

# NAG vs AG
for (i in 1:nrow(df)) {
  if(df[i, "donor.id"] %in% NAG) {
    subtype <- "NAG"
  }
  else if (df[i, "donor.id"] %in% AG) {
    subtype <- "AG"
  }
  else if (df[i, "donor.id"] %in% BM) {
    subtype <- "BM"
  }
  else {
    subtype <- "NA"
  }
  df[i, "subtype2"]  <- subtype
}

df$id <- NULL
df$doublet <- NULL
dt.1 <- AddMetaData(dt.1, metadata=df)


#########################
## Add soup results_b2###
#########################
soup_csv <- read.csv('/data/project/RCC_PBMC_HWS/workflow/singlecell/souporcell/outs_b2/clusters.tsv', header=TRUE, sep='')
soup_df <- as.data.frame(soup_csv)
donor_id <- data.frame(soup_df[,1:3])

colnames(donor_id) <- c('cell','doublet', 'id')
rownames(donor_id) <- donor_id$cell
donor_id$cell <- NULL
dt.2 <- AddMetaData(dt.2, metadata=donor_id)

df <- as.data.frame(dt.2@meta.data)
df$orig.ident <- NULL
df$nCount_RNA <- NULL
df$nFeature_RNA <- NULL
df$id <- as.character(df$id)

# Assign RCC patient ids
for (i in 1:nrow(df)){
  if (df[i,"id"] == '0'){
    donor.id <- 'RCC.8'
  }
  else if (df[i,"id"] == '1'){
    donor.id <- 'RCC.9'
  }
  else if (df[i,"id"] == '2'){
    donor.id <-'HC.1'
  }
  else if (df[i,"id"] == '4'){
    donor.id <- 'RCC.6'
  }
  else if (df[i,"id"] == '5'){
    donor.id <- 'HC.2'
  }
  else if (df[i,"id"] == '6'){
    donor.id <- 'RCC.2'
  }
  else {
    donor.id <- 'unassinged'
  }
  df[i,'donor.id'] <- donor.id
}

# Add RCC subtype
HC <- c("HC.1","HC.2")
NAG <- c("RCC.2")
AG <- c("RCC.6" , "RCC.7", "RCC.8","RCC.9")

# NAG vs AG
for (i in 1:nrow(df)) {
  if(df[i, "donor.id"] %in% NAG) {
    subtype <- "NAG"
  }
  else if (df[i, "donor.id"] %in% AG) {
    subtype <- "AG"
  }
  else if (df[i, "donor.id"] %in% HC) {
    subtype <- "HC"
  }
  else {
    subtype <- "NA"
  }
  df[i, "subtype"]  <- subtype
}

# NAG vs (AG vs BM)
HC <- c("HC.1","HC.2")
NAG <- c("RCC.2")
BM <- c("RCC.6")
AG <- ("RCC.7", "RCC.8","RCC.9")

# NAG vs AG
for (i in 1:nrow(df)) {
  if(df[i, "donor.id"] %in% NAG) {
    subtype <- "NAG"
  }
  else if (df[i, "donor.id"] %in% AG) {
    subtype <- "AG"
  }
  else if (df[i, "donor.id"] %in% BM) {
    subtype <- "BM"
  }
  else if (df[i, "donor.id"] %in% HC) {
    subtype <- "HC"
  }
  else {
    subtype <- "NA"
  }
  df[i, "subtype2"]  <- subtype
}

df$id <- NULL
df$doublet <- NULL
dt.2 <- AddMetaData(dt.2, metadata=df)
