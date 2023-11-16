# installation (use one of following methods) ----------------------------------------------------------------------------
# command line - *USED*
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
conda install pip
/home/songnsnow/anaconda3/envs/velocyto/bin/pip install velocyto

# R
install.packages('devtools')
library(devtools)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("pcaMethods")
install_github("velocyto-team/velocyto.R")

# docker
cd velocyto.R/dockers/debian9
docker build -t velocyto .
docker run --name velocyto -it velocyto