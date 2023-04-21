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
### Re-process B cell GEX  ###
##############################

# read in B cell GEX data
b_cells = readRDS("./datasets/b_cells.rds")

# Re-run SCTransform & Harmony with new variable genes
DefaultAssay(b_cells)="RNA"
indiv_bcell_data = SplitObject(b_cells, split.by = "batch_id")
indiv_bcell_data_sct = lapply(X = indiv_bcell_data,
                   FUN = SCTransform,
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE,
                   variable.features.n = 20000,
                   vars.to.regress = "percent.mt")


# remove low count batches
low_count_batches = b_cells@meta.data %>% dplyr::count(batch_id) %>% filter(n<10)
indiv_bcell_data_sct = indiv_bcell_data_sct[!(names(indiv_bcell_data_sct) %in% low_count_batches$batch_id)]

var.features = SelectIntegrationFeatures(object.list = indiv_bcell_data_sct[1:length(indiv_bcell_data_sct)], nfeatures = 5000)
b_cells = merge(x = indiv_bcell_data_sct[[1]], y = indiv_bcell_data_sct[2:length(indiv_bcell_data_sct)])
VariableFeatures(b_cells) = var.features

# Run PCA
b_cells = b_cells %>% RunPCA(assay.use="SCT")

# Run Harmony
b_cells = b_cells %>% RunHarmony(assay.use="SCT",group.by.vars="batch_id")
saveRDS(b_cells,"/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/datasets/b_cells_reintegrated.rds")
