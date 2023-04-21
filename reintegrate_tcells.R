##############################
###     Load packages      ###
##############################
library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(harmony)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(celldex)
library(SingleR)
library(edgeR)

setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined")


##############################
### Re-process T cell GEX  ###
##############################

# read in T cell GEX data
t_cells = readRDS("./datasets/t_cells.rds")

# Re-run SCTransform & Harmony with new variable genes
DefaultAssay(t_cells)="RNA"
indiv_tcell_data = SplitObject(t_cells, split.by = "batch_id")
indiv_tcell_data_sct = lapply(X = indiv_tcell_data,
                   FUN = SCTransform,
                   method = "glmGamPoi",
                   vars.to.regress = "percent.mt")

var.features = SelectIntegrationFeatures(object.list = indiv_tcell_data_sct, nfeatures = 3000)
t_cells = merge(x = indiv_tcell_data_sct[[1]], y = indiv_tcell_data_sct[2:length(indiv_tcell_data_sct)])
VariableFeatures(t_cells) = var.features

# Run PCA
t_cells = t_cells %>% RunPCA(assay.use="SCT")

# Run Harmony
t_cells = t_cells %>% RunHarmony(assay.use="SCT",group.by.vars="batch_id")
saveRDS(t_cells,"/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/datasets/t_cells_reintegrated.rds")
