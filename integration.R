# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(SoupX)
library(gridExtra)

# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/basic_qc/individual_post_qc_datasets/")

# bad batches
bad_batches = "S17624-E11"

cam_files = list.files(full.names = T)
cam_files = cam_files[!grepl(bad_batches,cam_files)]

# add in TUM samples
tum_files = list.files("/home/hpcjaco1/rds/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/individual_post_qc_datasets/",
                       full.names = T)

# combine CAM & TUM samples
files = append(cam_files,tum_files)

# read in post-QC datasets
post_qc_datasets = list()

read_in_datasets = function(x,cohort){
  message("reading in ",x)
  df = readRDS(x)
  df@meta.data$processing_site = ifelse(cohort=="Cam","CAM","TUM")

  df@meta.data$batch_id = ifelse(
    cohort=="Cam",
    paste0("CAM_",df@meta.data$batch_id),
    paste0("TUM_",stringr::str_remove(pattern=".rds",stringr::str_split(x,pattern="batch_")[[1]][2]))
    )
  df
}

# loop through all files
for(i in c(1:length(cam_files))){
  post_qc_datasets[[length(post_qc_datasets)+1]] = read_in_datasets(x=cam_files[i],cohort="Cam")
}
for(i in c(1:length(tum_files))){
  post_qc_datasets[[length(post_qc_datasets)+1]] = read_in_datasets(x=tum_files[i],cohort="TUM")
}

# check this has worked
message(length(post_qc_datasets))

# combine datasets
all_combo = merge(x = post_qc_datasets[[1]], y = post_qc_datasets[2:length(post_qc_datasets)])

# Integrate datasets by selecting integration features
var.features = SelectIntegrationFeatures(post_qc_datasets, nfeatures = 10000)
VariableFeatures(all_combo) = var.features

# Run PCA
all_combo = all_combo %>% RunPCA(assay.use="SCT")

# Run Harmony
all_combo = all_combo %>% RunHarmony(assay.use="SCT",group.by.vars=c("batch_id","processing_site"))

# save
saveRDS(all_combo, file = "/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/gex_pbmc_csf_combined.rds")
