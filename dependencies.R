# -- Script used to install various dependencies to cluster nodes.

install.packages("tidyverse", "magrittr", "RPMM", "doParallel", repos = "http://cran.cnr.berkeley.edu/")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("impute")