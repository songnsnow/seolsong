# load libraries ------------------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)                # load SCTransform
library(glmGamPoi)                  # load glmGamPoi


# Load PBMC dataset & create SeuratObject -----------------------------------------------------------------
dt <- Read10X(data.dir = "/PATH/filtered_feature_bc_matrix")
dt <- CreateSeuratObject(counts = dt, project = "project_name", min.cells = 3, min.features = 200)


# QC ------------------------------------------------------------------------------------------------------
dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^MT-")
dt <- PercentageFeatureSet(dt, "^RP[SL]", col.name = "percent.rb")
dt <- PercentageFeatureSet(dt, "^HB[^(P)]", col.name = "percent.hb")
dt <- PercentageFeatureSet(dt, "PECAM1|PF4", col.name = "percent.plat")

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb", "percent.plat")
pt <- VlnPlot(dt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()
ggsave('/PATH/qc_vlnplot_bf.png',width=5, height=7,scale=2)

FeatureScatter(dt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.5)
ggsave('/PATH/qc_featurescatter_bf.png',width=7, height=5,scale=2)

save1 <- dt

# subset
dt <- subset(dt, subset = nFeature_RNA > 200 & nFeature_RNA < 7500) 
dt <- subset(dt, subset = percent.mt < 20) 
dt <- subset(dt, subset = percent.hb < 3)

save2 <- dt


# pre-processing (SCTransform) ----------------------------------------------------------------------------------------
dt <- PercentageFeatureSet(dt, pattern = "^MT-", col.name = "percent.mt")
dt <- SCTransform(dt, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
dt <- RunPCA(dt, verbose=FALSE)
dt <- RunUMAP(dt, dims = 1:50, verbose=FALSE)
dt <- FindNeighbors(dt, dims = 1:50, verbose=FALSE)
dt <- FindClusters(dt, resolution = 1.1, verbose=FALSE)

