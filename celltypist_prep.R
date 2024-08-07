# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)

# set WD
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/datasets")

# read in args
args = commandArgs(trailingOnly = T)
x = as.character(as.numeric(args[1]) - 1)
message("Doing cluster ",x)

# read in outputs from deconvolution & integration step
all_combo = readRDS("all_combo_with_pheno.rds")


# export each cluster for celltypist

this_cluster = subset(all_combo, seurat_clusters == x)

# get counts
count_mat = this_cluster@assays$RNA@counts

# transform
count_mat = as.matrix(count_mat)

# transpose
count_mat = t(count_mat)

# write to file
outfile = paste0("all_combo_celltypist_counts_cluster_",x,".tsv")
df = data.frame(count_mat)
write_tsv(df,file = outfile )
