# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(gridExtra)

# set WD
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/datasets")

# read in outputs from deconvolution & integration step
all_combo = readRDS("gex_pbmc_csf_combined.rds")

# remove unwanted assays
all_combo[["ADT"]] = NULL
all_combo[["SP"]] = NULL
all_combo[["HTO"]] = NULL

# plots
p1=DimPlot(all_combo,reduction="pca",group.by="batch_id")+NoLegend()
p2=DimPlot(all_combo,reduction="harmony",group.by="batch_id")+NoLegend()
png(res=300,units="in",height=16,width=16,file="harmony_vs_pcs_batch.png")
grid.arrange(p1,p2,nrow=1)
dev.off()

p1=DimPlot(all_combo,reduction="harmony",split.by="processing_site",group.by="source")
p2=DimPlot(all_combo,reduction="harmony",split.by="source",group.by="processing_site")
png(res=300,units="in",height=5,width=8,file="harmony_source_site.png")
grid.arrange(p1,p2,nrow=1)
dev.off()

# Elbow plot
png(res=300,units="in",height=8,width=8,file="elbow.png")
ElbowPlot(all_combo,reduction="harmony")
dev.off()

set.seed(123456)

# remove RBCs
all_combo[['percent.hb']] = PercentageFeatureSet(all_combo, features=c("HBB","HBA1","HBA2"))
table(all_combo@meta.data$percent.hb>0.01)
all_combo = subset(all_combo, subset = percent.hb < 0.01)


# UMAP over a range of PCs and resolutions
# repeat with different parameters to test robustness of clusters
umap_change_pcs = function(x,y){
  message(paste0("Start time:",date()))
  message(paste0("Running UMAP & clustering with ",x," PCs and a resolution of ",y))
  all_combo_test = all_combo %>% RunUMAP(reduction="harmony", dims=1:x) %>%
  FindNeighbors(reduction="harmony",dims=1:x) %>%
  FindClusters(resolution=y,group.singletons = FALSE)
  message(paste0("Finished UMAP & clustering with ",x,"PCs and a resolution of",y,". Moving on."))
  message(paste0("End time:",date()))
  png(res=300,units="in",height=8,width=12,file=paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,".png"))
  print(DimPlot(all_combo_test,label=T))
  dev.off()

  saveRDS(all_combo_test,paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,".rds"))
}

mat = matrix(c(50,50,50,50,50,2.5,3,1,1.5,2),ncol=2)
mapply(umap_change_pcs,x=mat[,1],y=mat[,2]) # takes several hours - not routinely run
