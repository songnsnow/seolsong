#/usr/bin/env R conda seurat5
#-----------------------------------------------------------------------
# description : Drawing contour plot in R
# author      : rlo
# date        : 230111
# notes       : 
#-----------------------------------------------------------------------
##library loading ------------------------------------------------------
# available libraries
library(Seurat)
library(remotes)
library(BPCells)
library(presto)
library(glmGamPoi)
library(SeuratWrappers)
library(TFBSTools)
library(ggplot2)

# INPUT ##########################################################
setwd("/data/project/RCC_PBMC_HWS/workflow/singlecell/seurat/") # set working directory
# load data (skip if continued)
dt <- readRDS('rds/06-1_2.rds')
#####################################################################

##Run umap
dt <- PercentageFeatureSet(dt, pattern = "^MT-", col.name = "percent.mt")
dt <- SCTransform(dt, vars.to.regress = "percent.mt", verbose = TRUE) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:50) %>%
    FindClusters(resolution = 1.1) %>%
    RunUMAP(dims = 1:50)


# Seurat 객체에서 UMAP dimension을 가져오기
umap_data <- as.data.frame(dt@reductions$umap)

# DimPlot 그리기
dim_plot <- DimPlot(dt, reduction = "umap", group.by = 'integrated.cell.type', label = TRUE) +
  theme_classic()
ggsave('figures/07-1_umap_by.celltypye.png',dim_plot, width=5, height=5,scale=1.5)

# Contour Plot 추가
dim_plot + 
  geom_density_2d(data = umap_data, aes(x = umap_1, y = umap_2), color = "blue", fill = "blue", alpha = 0.5)
ggsave('figures/07-1_umap_density.png',dim_plot, width=5, height=5,scale=1.5)


####
umap_data <- as.data.frame(dt@reductions$umap@cell.embeddings)
colnames(umap_data) <- c("umap_1", "umap_2")

# 클러스터 정보 추가
umap_data$celltype <- as.factor(dt@active.ident)
umap_data$subtype <- as.factor(dt@meta.data$subtype2)

# Contour plot 그리기
plot <- ggplot(umap_data) +
  geom_point(aes(x = umap_1, y = umap_2, color = celltype), show.legend = TRUE) +
  geom_density_2d(aes(x = umap_1, y = umap_2, fill = celltype), alpha = 0.1) +
  stat_density_2d(aes(x = umap_1, y = umap_2, color = subtype, fill = celltype),
                 geom = "polygon", alpha = 0.3) + 
  theme_classic() +
  scale_x_continuous(limits = c(min(umap_data$umap_1) - 0.1 * diff(range(umap_data$umap_1)),
                                 max(umap_data$umap_1) + 0.1 * diff(range(umap_data$umap_1)))) +
  scale_y_continuous(limits = c(min(umap_data$umap_2) - 0.1 * diff(range(umap_data$umap_2)),
                                 max(umap_data$umap_2) + 0.1 * diff(range(umap_data$umap_2)))) +
  facet_wrap(~subtype, scales = "free")

# 그림 저장
ggsave('figures/07-1_contour_plot.png', plot, width = 8, height = 6)
