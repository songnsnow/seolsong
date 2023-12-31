# command line ----------------------------------------------------------------------------------------------------
ulimit -s unlimited 


# load library ----------------------------------------------------------------------------------------------------
library(infercnv)
library(Seurat)


# create infercnv obj ---------------------------------------------------------------------------------------------
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = as.matrix(obj@active.ident),
                                    delim = "\t",
                                    gene_order_file = "/PATH/inferCNV/hg38_gencode_v27.txt",
                                    ref_group_names = c('ann1','ann2','ann3','ann4')) # error if datasets with different normalization assays are merged


# run inferCNV ----------------------------------------------------------------------------------------------------
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1, 
                             out_dir = '/PATH/inferCNV/outs', 
                             cluster_by_groups = TRUE, 
                             denoise = TRUE,
                             HMM = TRUE)