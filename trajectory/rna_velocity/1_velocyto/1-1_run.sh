# run --------------------------------------------------------------------------------------------------------------------
# command line - *USED*
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

velocyto run10x /data/project/RCC_PBMC_HWS/H372TDSX7/run /data/project/RCC_PBMC_HWS/H372TDSX7/refdata-gex-GRCh38-2020-A/genes/genes.gtf