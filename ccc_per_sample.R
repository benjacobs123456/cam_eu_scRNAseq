#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(gridExtra)
library(MASS)
library(reshape2)
library(liana)
library(nichenetr)


#######################################
# Read in data
#######################################

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/ccc/")

# get args
args = commandArgs(trailingOnly=T)
index = as.numeric(args[1])

# get sc data
sc_dat = readRDS("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/all_combo_with_updated_pheno.rds")

# subset sc data to just this donor
message("Getting donor ",index," of ",length(unique(sc_dat@meta.data$iid)))
this_donor = unique(sc_dat@meta.data$iid)[index]
message(this_donor)
sc_dat = subset(sc_dat,iid == this_donor & source=="CSF")


#######################################
# LIANA
#######################################

# set SCT as default
DefaultAssay(sc_dat) = "SCT"

# run LIANA
liana_res = liana_wrap(sc_dat, idents_col = "cell_type") %>%
  liana_aggregate()

# add iid
liana_res = liana_res %>%
  mutate(iid = this_donor,
  phenotype = sc_dat@meta.data$phenotype[1])

# save to file
write_csv(liana_res,
paste0("./per_sample/liana_res_",index,".csv")
)
