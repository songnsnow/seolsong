# import packages --------------------------------------------------------------------------------------------------------------
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import sklearn as sk
import matplotlib.pyplot as pl
import scvelo as scv
import ipywidgets as widgets
import seaborn as sns


# Set parameters for plots, including size, color, etc. -------------------------------------------------------------------------
scv.set_figure_params(style="scvelo")
pl.rcParams["figure.figsize"] = (10,10)
Colorss = sns.color_palette("tab10",26)


# Define Path to the matrix output from cellranger ------------------------------------------------------------------------------
Path10x = '/data/project/RCC_PBMC_HWS/H372TDSX7/run/outs/filtered_feature_bc_matrix'
rcc = sc.read_10x_mtx(Path10x,var_names='gene_symbols',cache=True)   # Read matrix
rcc   # Print information on this new object


# import clusters and projection ------------------------------------------------------------------------------------------------
# Read clusters
Clusters = pd.read_csv("/data/project/RCC_PBMC_HWS/SS/scvelo/rcc_clusters.csv", delimiter=',',index_col=0)

# Create a list with only selected RCC Barcodes
RCC_BCs = Clusters.index

# Read UMAP embeddings
UMAP = pd.read_csv("/data/project/RCC_PBMC_HWS/SS/scvelo/rcc_umap.csv", delimiter=',',index_col=0)
UMAP = UMAP.loc[RCC_BCs,]   # Select RCC barcodes in the projection
UMAP = UMAP.to_numpy()   # Transform to Numpy

rcc = rcc[RCC_BCs]   # Filter cells for only selected RCC cells
rcc.obs['Loupe'] = Clusters   # Import cluster information into rcc object
rcc.obsm["X_umap"] = UMAP   # Import UMAP information into rcc object

rcc   # Print information about the object


# import sliced/unsliced counts ------------------------------------------------------------------------------------------------
# Read velocyto output
VelRCC = scv.read('/data/project/RCC_PBMC_HWS/H372TDSX7/run/velocyto/run.loom', cache=True)

# Merge velocyto with cellranger matrix
scv.utils.clean_obs_names(rcc)
scv.utils.clean_obs_names(VelRCC)
rcc = scv.utils.merge(rcc, VelRCC)


# run velocity analysis --------------------------------------------------------------------------------------------------------
# Standard scvelo processing to run Dynamical Mode
scv.pp.filter_and_normalize(rcc, min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(rcc, n_pcs=30, n_neighbors=30)

# Estimate RNA velocity and latent time
scv.tl.recover_dynamics(rcc)
scv.tl.velocity(rcc, mode='dynamical')
scv.tl.velocity_graph(rcc)
scv.tl.recover_latent_time(rcc)


# Visualization ------------------------------------------------------------------------------------------------------------------
# Generate plot with UMAP and cluster
scv.pl.velocity_embedding_stream(rcc,basis="umap",color="Loupe",title='RCC',fontsize=20,legend_fontsize=20,min_mass=2,palette=Colorss,save='/data/project/RCC_PBMC_HWS/SS/scvelo/scVelo-umap-cluster.png')
# AttributeError: property 'categories' of 'Categorical' object has no setter

scv.pl.velocity_embedding_stream(rcc,basis="umap",color="velocity_pseudotime",title='RCC',fontsize=20,legend_fontsize=20,min_mass=2,palette=Colorss,save='/data/project/RCC_PBMC_HWS/SS/scvelo/scVelo-umap-velocity_pseudotime.png')

# Generate plot with UMAP and latent time
scv.pl.velocity_embedding_stream(rcc,basis="umap",color="latent_time",title='RCC',fontsize=20,legend_fontsize=20,min_mass=2,color_map="plasma",save='/data/project/RCC_PBMC_HWS/SS/scvelo/scVelo-umap-latent_time.png')



sc.pl.violin(rcc, keys='latent_time',groupby="Clusters",order=[4,1,5,6,7,3,2],save='/data/project/RCC_PBMC_HWS/SS/scvelo/scVelo-violin-latent_time.png')

Genes=["RETN","LTF","CAMP","ACTB","GCA","LCN2",
"S100A8","MYL6","S100A9","FCGR3B","S100A11","FTH1","IFIT1",
"IFITM3","IFIT3","ISG15","IFIT2","RPS9","NEAT1","MALAT1","NFKBIA","CXCL8"]

scv.pl.heatmap(rcc, var_names=Genes, sortby='latent_time', col_color='Loupe', n_convolve=100,figsize=(16,8),yticklabels=True,sort=True,colorbar=True,show=True,layer="count",save='/data/project/RCC_PBMC_HWS/SS/scvelo/scVelo-heatmap-latent_time.png')



# rank velocity genes
scv.tl.rank_velocity_genes(rcc, groupby='Clusters', min_corr=.3)
df = scv.DataFrame(rcc.uns['rank_velocity_genes']['names'])
df.head()