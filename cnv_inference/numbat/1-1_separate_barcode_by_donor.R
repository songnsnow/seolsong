# separate barcode by donor --------------------------------------------------------------------------------------------------------
# load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# load data
dt <- readRDS('/data/project/RCC_PBMC_HWS/SS/seurat/rds/seurat_4_1.rds')
dt <- SetIdent(dt, value = dt@meta.data$id)
rcc.1 <- subset(dt, idents = c('RCC.1')) #1024
rcc.2 <- subset(dt, idents = c('RCC.2')) #1393
rcc.3 <- subset(dt, idents = c('RCC.3')) #1476
rcc.4 <- subset(dt, idents = c('RCC.4')) #255
rcc.5 <- subset(dt, idents = c('RCC.5')) #1670
rcc.6 <- subset(dt, idents = c('RCC.6')) #948
rcc.7 <- subset(dt, idents = c('RCC.7')) #1105
rcc.8 <- subset(dt, idents = c('RCC.8')) #12418

# make barcode file by donor
barcodes <- rownames(rcc.1@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.1_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.2@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.2_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.3@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.3_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.4@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.4_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.5@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.5_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.6@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.6_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.7@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.7_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
barcodes <- rownames(rcc.8@meta.data)
write.table(barcodes, file = '/data/project/RCC_PBMC_HWS/SS/numbat/rcc.8_barcodes.tsv', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# make matrix file by donor
matrix <- rcc.1@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.1_matrix.rds')
matrix <- rcc.2@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.2_matrix.rds')
matrix <- rcc.3@assays[["RNA"]]@counts
saveRDS(matrix, '/data/project/RCC_PBMC_HWS/SS/numbat/input/rcc.3_matrix.rds')
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
