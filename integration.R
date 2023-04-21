# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
library(SoupX)
library(DoubletFinder)
library(gridExtra)

# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/basic_qc/individual_post_qc_datasets/")

# bad batches
bad_batches = "S17624-E11"

files = list.files()
files = files[!grepl(bad_batches,files)]

post_qc_datasets = lapply(files,function(x){
  message("reading in ",x)
  df = readRDS(x)
  return(df)
})

all_combo = merge(x = post_qc_datasets[[1]], y = post_qc_datasets[2:length(post_qc_datasets)])


# Integrate datasets
var.features = SelectIntegrationFeatures(post_qc_datasets, nfeatures = 5000)
VariableFeatures(all_combo) = var.features

# Run PCA
all_combo = all_combo %>% RunPCA(assay.use="SCT")

# Run Harmony
all_combo = all_combo %>% RunHarmony(assay.use="SCT",group.by.vars="batch_id")

# save
saveRDS(all_combo, file = "/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/datasets/gex_pbmc_csf_combined.rds")
