# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
library(gridExtra)

# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/datasets/")

# read in outputs from deconvolution & integration step
all_combo = readRDS("all_combo_with_updated_pheno.rds")

# Hive off B & T cell clusters for VDJ analysis

new_cluster_ids = read_csv("cluster_identities.csv")

# hive off B cells for VDJ analysis
b_cell_clusters =  c("B cells","Plasma cells")
b_cells = subset(all_combo, subset = cell_type %in% b_cell_clusters)
saveRDS(b_cells,"b_cells.rds")

# hive off T cells for VDJ analysis
t_cell_clusters = c("CD4 T cells","CD8 T cells","Tregs","MAIT")
t_cells = subset(all_combo, subset = cell_type %in% t_cell_clusters)
saveRDS(t_cells,"t_cells.rds")
