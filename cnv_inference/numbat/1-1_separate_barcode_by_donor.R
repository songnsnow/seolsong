# separate barcode by donor --------------------------------------------------------------------------------------------------------
# load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# load data
dt <- readRDS('/PATH/seuratObject.rds')
dt <- SetIdent(dt, value = dt@meta.data$id)
d.1 <- subset(dt, idents = c('Donor1')) 
d.2 <- subset(dt, idents = c('Donor2')) 
d.3 <- subset(dt, idents = c('Donor3')) #1476

# make barcode file by donor
barcodes <- rownames(d.1@meta.data)
write.table(barcodes, file = '/PATH/numbat/d.1_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(d.2@meta.data)
write.table(barcodes, file = '/PATH/numbat/d.2_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(d.3@meta.data)
write.table(barcodes, file = '/PATH/numbat/d.3_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# make matrix file by donor
matrix <- d.1@assays[["RNA"]]@counts
saveRDS(matrix, '/PATH/numbat/input/d.1_matrix.rds')
matrix <- d.2@assays[["RNA"]]@counts
saveRDS(matrix, '/P/numbat/input/d.2_matrix.rds')
matrix <- d.3@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/d.3_matrix.rds')
matrix <- rcc.4@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.4_matrix.rds')
matrix <- rcc.5@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.5_matrix.rds')
matrix <- rcc.6@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.6_matrix.rds')
matrix <- rcc.7@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.7_matrix.rds')
matrix <- rcc.8@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.8_matrix.rds')
