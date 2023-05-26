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

#######################################
# Read in data
#######################################
# get args
args = commandArgs(trailingOnly = TRUE)
message(args)

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/")

# Read in data
all_combo = readRDS("../datasets/all_combo_phenotypes_new_cluster_ids.rds")

# read in SNPs
snps_for_eqtl = read_tsv("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc.bim",col_names=F)
message("SNPs for eQTL testing")
message(nrow(snps_for_eqtl))

# create file with new IDs for plink
codex = read_table("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/data/genotypes/GT_ID.txt",col_names=FALSE)

# sort out long cambridge IDs
new_ids = lapply(codex$X2, function(x){
  if(grepl("C00",x)){
      new_id = str_remove(pattern="C00TU0",x)
      new_id = str_split(pattern="a",new_id)[[1]][1]
      new_id = paste0("TU",new_id)
      new_id = str_remove(pattern="v1",new_id)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else if(grepl("v",x)){
      new_id = str_remove(pattern="v1",x)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else {
      return(x)
  }
})

codex$sc_id = unlist(new_ids)
colnames(codex)[1] = "genotyping_id"
codex = codex %>% distinct(genotyping_id,sc_id)

# get TUM sample data
tum_phenotypes =read_csv("/rds/user/hpcjaco1/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/phenotypes_all.csv") %>%
dplyr::select(GSA_ID,PatID) %>%
filter(!is.na(GSA_ID))
colnames(tum_phenotypes) = c("genotyping_id","sc_id")

# combine
codex = bind_rows(codex,tum_phenotypes) %>%
  distinct(sc_id,.keep_all=T)
colnames(codex)[2] = "iid"

# get TUM sample data for batch 2
tum_phenotypes2 = read_table("/rds/user/hpcjaco1/hpc-work/tum_genos_round2/new_TUM_GenIDs.txt") %>%
filter(!is.na(GenID)) %>%
dplyr::select(2,1)
colnames(tum_phenotypes2) = c("genotyping_id","iid")

# combine
codex = bind_rows(codex,tum_phenotypes2) %>%
  distinct(iid,.keep_all=T)

# subset single-cell data to those in genotyped cohort
eqtl_dataset = subset(all_combo, iid %in% codex$iid)

# add genotyping ID
eqtl_dataset@meta.data = eqtl_dataset@meta.data %>%
  left_join(codex,by="iid") %>%
  mutate(full_cell_id = paste0(source,"_",genotyping_id,"_",rownames(eqtl_dataset@meta.data)))
rownames(eqtl_dataset@meta.data) = colnames(eqtl_dataset)

# rename cells with full ID
eqtl_dataset = RenameCells(eqtl_dataset,new.names = eqtl_dataset@meta.data$full_cell_id)

# set DefaultAssay & remove SCT for size
DefaultAssay(eqtl_dataset)="RNA"
eqtl_dataset[['SCT']] = NULL

# filter to cell type 
cell_type_to_test =  unique(eqtl_dataset@meta.data$cell_type_crude)[as.numeric(args[1])]
cell_type_no_space =  str_replace_all(cell_type_to_test," ","_")
eqtl_dataset = subset(eqtl_dataset,subset = cell_type_crude == cell_type_to_test)

# save file
outfile = paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/eqtl_sc_data_",cell_type_no_space)
saveRDS(eqtl_dataset,outfile)
