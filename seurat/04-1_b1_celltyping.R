#/usr/bin/env R conda seu5
#-----------------------------------------------------------------------
# description : Celltyping by canonical markers
# author      : rlo
# date        : 230103
# notes       : 
#-----------------------------------------------------------------------

# INPUT ##########################################################
setwd("/data/project/RCC_PBMC_HWS/workflow/singlecell/seurat/b1") # set working directory
# load data (skip if continued)
dt <- readRDS('rds/03-1.rds')
#####################################################################

# Celltyping ---------------------------------------------------------------------------------------------------------
## big cell types
dt <- SetIdent(dt, value = dt@meta.data$seurat_clusters)

# FeaturePlot
mks <- c('MS4A1','CD38','CD27','JCHAIN','CD14','FCGR3A','SELL','CX3CR1','HLA-DR','ITGB2','ITGAM','NRP1','CD1C','PPBP','CD3E','SLC4A10','FOXP3','CD4','CD8B','LEF1','IL7R','CD44','PRF1','NCAM1') 

FeaturePlot(dt, features = mks, ncol=4)
ggsave('figures/04-1_mks.featureplot.all.png',width=5, height=7, scale=3)

# DotPlot
clsts <- c('0','1','3','8','10','13','16','17','2','12','7','11','14','6','9','18','20','4','5','19','15','26')
Idents(dt) <- factor(Idents(dt), levels= clsts)

dt <- SetIdent(dt, value = dt@meta.data$seurat_clusters)
mks <- c('CD247','CD3E','CD8A','CD8B','CD4','CD14','CD16','CD56','FCGR3A','NCAM1','CD1C','NRP1','FCER1A','MS4A1','CD38','CD27','JCHAIN','PPBP')

dp <- DotPlot(object = dt, features = mks, cluster.idents = TRUE) + theme(axis.text.x=element_text(angle=90)) + coord_flip() + theme_classic()
ggsave('figures/04-1_mks.dotplot.all.png',plot=dp,width=6, height=3, scale=2)



## cell type subsets---------------------------------------------------------------------------------------------------
# DotPlot for T cell states
dt <- SetIdent(dt, value = dt@meta.data$seurat_clusters)
Ts <- subset(dt, idents = c('0','1','2','3','4','6','7','8','10','11','13','14','16','17','21','22','24','25','26'))
#Ts <- c('0','1','2','3','4','6','7','8','10','11','13','14','16','17','21','22','24','25')
# Idents(Ts) <- factor(Idents(Ts), levels= clsts)

mks <- c('PECAM1','CD4','CCR7','CD27','BCL2','LEF1','CD44','CD28','IL7R','GZMB','IL2RB','SPN','KLRG1','PRF1','SLC4A10','FOXP3','TRDC','CD3D','CD45RO','CD62L')
VlnPlot(dt, features = mks, idents= c(0,1,2,3,4,6,7,8,10,11,13,14,16,17,21,22,24,25,26),assay='SCT',ncol=3) & ylim(0,3) & geom_boxplot(width=0.15,fill='white')

DotPlot(object = Ts, features = mks,col.min=0,cluster.idents = TRUE) + theme(axis.text.x=element_text(angle=90)) + coord_flip() + theme_classic()
ggsave('figures/04-1_mks.dotplot.T.png',width=6, height=4, scale=1)

#Vionlin plot for NK states
nk <- c('GNLY','KLRF1', 'TRDC','FCER1G','KLRC1')
VlnPlot(dt, features=nk, idents=c(4,10,14,24,7,17),assay='SCT',ncol=3) & ylim(0,3) & geom_boxplot(width=0.15,fill='white')
ggsave('figures/04-1_mks.vlnplot.nk.png',width=6, height=3,scale=2)

#Violin plot for B states
B <- c()
VlnPlot(dt, features=B, idents=c(12,18,19,20),assay='SCT',ncol=3) & ylim(0,3) & geom_boxplot(width=0.15,fill='white')
ggsave('figures/04-1_mks.vlnplot.B.png',width=6, height=3, scale=2)

# assigning cell id -------------------------------------------------------------------------------------------------
dt <- SetIdent(dt, value = dt@meta.data$seurat_clusters)
#id <- c('T cell','T cell','T cell','T cell','T cell','classical MC','T cell','T cell','T cell','classical MC','T cell',
#        'T cell','B cell','T cell','','classical MC','T cell','T cell','B cell','B cell','B cell',
#        'T cell','T cell','NK','T cell','Platelet','T cell')

id <- c('CD8+ T','CD4+ T','CD4+ T','CD4+ T','CD56dim NK','classical MC','CD4+ T','CD56dim NK','CD8+ T','classical MC','CD56dim NK',
        'CD8+ T','Naive B','CD4+ T','CD56dim NK','classical MC','CD4+ T','CD56dim NK','Plasma','Naive B','Memory B',
        'CD8+ T','MAIT','Non-classical NK','CD56dim NK','platelet','CD4+ T')

lbd <- dt
names(id) <- levels(lbd)
dt <- RenameIdents(lbd, id)
dt$cell.type <- dt@active.ident


# Plot by id
DimPlot(dt, reduction = "umap", group.by='cell.type', label=TRUE) + theme_classic()
ggsave('figures/04-1_processed_umap_by.celltype.png',width=5, height=5,scale=1.5)


#####Azimuth

# run Azimuth for reference ------------------------------------------------------------------------------------------
# load libraries
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

# The RunAzimuth function can take a Seurat object as input
dt <- RunAzimuth(dt, reference = "pbmcref")
dp <- DimPlot(dt, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
ggsave('figures/04-1_azi_l2.png',plot=dp,width=6, height=6, scale=2)

head(dt@meta.data)

saveRDS(dt,'rds/04-1.rds')

###########################
obj <- readRDS('rds/04-1.rds')
obj$percent.mt <- NULL
obj$nCount_SCT <- NULL
obj$nFeature_SCT <- NULL
obj$SCT_snn_res.1.1 <- NULL
obj$predicted.celltype.l1.score <- NULL
obj$predicted.celltype.l2.score <- NULL
obj$predicted.celltype.l3.score <- NULL
obj$mapping.score <- NULL
saveRDS(obj, 'rds/04_2.rds')
