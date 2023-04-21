#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(MASS)
library(SingleCellExperiment)
library(Matrix.utils)
library(reshape2)
library(RNOmni)
library(coloc)

#######################################
# Read in data
#######################################

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/")

files = list.files(pattern="full_results.tsv")

files_in = purrr::map(files,function(x){
  read_tsv(x) %>%
    mutate(CHR = X.CHROM) %>%
    filter(!is.na(P))
})

eqtl_res = do.call("bind_rows",files_in)
saveRDS(eqtl_res,"full_eqtl_results_all.rds")
