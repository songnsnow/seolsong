# 0. downsample cells (if needed) -----------------------------------------------------------------------
# downsample cells
obj <- readRDS('/PATH/seuratObject.rds')
obj <- SetIdent(obj, value = obj@meta.data$metadata)

s1 <- subset(obj,idents='metadata_1', downsample=500, seed=1)
s2 <- subset(obj,idents='metadata_2', downsample=500, seed=1)
s3 <- subset(obj,idents='metadata_3', downsample=500, seed=1)

merged <- merge(s1, c(s2, s3))
merged <- SetIdent(merged, value = merged@meta.data$metadata)


# 1. cell type annotations ------------------------------------------------------------------------------
# annotations file 
dt <- readRDS('/PATH/seuratObject.rds')
df <- data.frame(dt$type)
df <- tibble::rownames_to_column(df, "cell")
write.table(df, '/PATH/inferCNV/annotations.txt', quote = FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

# annotations - set (creating inferCNV object):
merged <- SetIdent(merged, value = merged@meta.data$metadata)   # metadata = cell type annotations
annotations_file = as.matrix(merged@active.ident)   # in inferCNV object


# 2. raw counts matrix ----------------------------------------------------------------------------------
# count matrix (in inferCNV)
obj <- readRDS('/PATH/seuratObject.rds')
counts_matrix <- GetAssayData(obj, slot="counts")
counts_matrix <- GetAssayData(merged, slot="counts")
