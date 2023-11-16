# load libraries --------------------------------------------------------------------------------------------------------------------
library(Seurat)


# read seuratObject -----------------------------------------------------------------------------------------------------------------
obj <- readRDS('/data/project/RCC_PBMC_HWS/SS/seurat/rds/seurat_4_2.rds')
obj


# export UMAP embeddings ------------------------------------------------------------------------------------------------------------
umap <- obj[["umap"]]@cell.embeddings
write.csv(umap, "/data/project/RCC_PBMC_HWS/SS/scvelo/rcc_umap.csv")


# export clusters -------------------------------------------------------------------------------------------------------------------
#change csv format...?
cls <- data.frame(clusters=obj$seurat_clusters)
write.csv(cls, "/data/project/RCC_PBMC_HWS/SS/scvelo/rcc_clusters.csv")