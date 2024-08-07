
args <- commandArgs(TRUE)

# Load packages
.libPaths("/opt/R/lib/R/library/")
library(reticulate)
library(Seurat)
library(SingleCellExperiment)
library(plyr)
library(Matrix)
library(Matrix.utils)
library(tidyselect)
library(tidyr)
library(dplyr)

outputfolder = args[1]

# read in gene expression data
data = readRDS(args[2]) #/xx/xxx/SC/CAM_TUM/ewtl/CSF/all_combo_renormalizedperbatch.rds 

# define source
selsource = args[3] # CSF or PBMC
data = subset(data, subset = source == selsource)

# convert to sce object
data.sce = as.SingleCellExperiment(data, assay="SCT")

scaledcounts = assay(data.sce,"scaledata")

# average counts per individual and cell type
groups = colData(data.sce)[, c("phenotype","iid")]
groups$BC = rownames(groups)
scaledcounts = t(scaledcounts); scaledcounts = data.frame(scaledcounts)
scaledcounts$BC = rownames(scaledcounts)
scaledcounts = join(scaledcounts,unique(data.frame(groups)),by="BC",type="left")
meanscaledcounts = scaledcounts %>% group_by(phenotype,iid) %>% summarise_at(vars(-BC), mean)
meanscaledcounts$phenotype = NULL
meanscaledcounts = data.frame(meanscaledcounts)
rownames(meanscaledcounts) = meanscaledcounts$iid
meanscaledcounts$iid = NULL

#PC calculation
meanscaledcounts1 = meanscaledcounts[ , which(apply(meanscaledcounts, 2, var) != 0)]
pca = prcomp(meanscaledcounts1,scale=T)
write.csv(pca$x,paste0(outputfolder,"/","PCs.csv"))