#/usr/bin/env R seu5
#-----------------------------------------------------------------------
# description : load data & create SeuratObject
# author      : songnsnow
# date        : 231229
# notes       : 
#-----------------------------------------------------------------------
# INPUT ##########################################################
setwd("/data/project/RCC_PBMC_HWS/workflow/singlecell/seurat")   # set working directory
getwd()     # check

# batch 1
fn <- 'b1'  # filename
h5_path <- '/data/project/RCC_PBMC_HWS/H372TDSX7/run/outs/filtered_feature_bc_matrix.h5'   # filtered feature bc matrix h5 file PATH

# batch 2 (after cellbender)
fn <- 'b2'  # filename
h5_path <- '/data/project/RCC_PBMC_HWS/workflow/singlecell/cellbender/b2/b2_cellbender_2_filtered.h5'   # filtered feature bc matrix h5 file PATH

# batch 2 (no cellbender)
fn <- 'b2'  # filename
h5_path <- '/data/project/RCC_PBMC_HWS/workflow/singlecell/cellranger/RCC_GEX/rcc_10x/outs/filtered_feature_bc_matrix.h5'   # filtered feature bc matrix h5 file PATH
##################################################################


# set up directories
dir.create(fn)
setwd(fn)
dir.create("rds")
dir.create("figures")


# load libraries
library(Seurat)
library(dplyr)
library(scCustomize)

# Load PBMC dataset & create SeuratObject
# rcc.10x <- Read10X_h5(h5_path)
rcc.10x <- Read_CellBender_h5_Mat(h5_path)  # read from cellbender output
rcc <- CreateSeuratObject(counts = rcc.10x, project = "rcc", min.cells = 3, min.features = 200)
dt <- rcc


saveRDS(dt,'rds/01-1.rds')