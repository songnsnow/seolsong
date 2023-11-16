# install SeuratDisk for as.loom() function
install.packages('Seurat')
library('Seurat')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library('SeuratDisk')  # need SeuratDisk for as.loom() function

# read seuratObject
obj <- readRDS('/data/project/RCC_PBMC_HWS/SS/seurat/rds/seurat_4_2.rds')
obj

# convert seuratObject to loom
rcc_loom <- as.loom(obj, filename="/data/project/RCC_PBMC_HWS/SS/scvelo/rcc.loom", verbose=FALSE)
rcc_loom$close_all()  # always remember to close loom files when done