# installations ------------------------------------------------------------------------------------------
install.packages('Seurat')
install.packages('dplyr')
install.packages('patchwork')
install.packages('ggplot2')
install.packages("sctransform")     # install SCTransform from CRAN
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("glmGamPoi")   # install glmGamPoi from biocmanager