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
all_combo = readRDS("gex_pbmc_csf_combined.rds")

# read in qc metrics for each batch
qc_files = list.files("../basic_qc/",full.name = TRUE)
qc_files = qc_files[grepl("qc_metrics",qc_files)]
qc_metrics_list = lapply(qc_files,function(x){
  df = read_csv(x)
  return(df)
})
qc_metrics = do.call("bind_rows", qc_metrics_list)
write_csv(qc_metrics,"../basic_qc/overall_qc_metrics.csv")

# make and save QC plots
qc_plot = function(y){
  p=ggplot(qc_metrics,aes(x,qc_metrics[[y]],fill=source))+
  geom_col()+
  labs(x="Batch",y=y)+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))
  return(p)
}
plots = lapply(colnames(qc_metrics)[c(5,4,10)],qc_plot)
png(res=300,units="in",height=16,width=16,file="../basic_qc/qc_plots.png")
do.call("grid.arrange",c(plots,ncol=1))
dev.off()

# cells per donor
counts = all_combo@meta.data %>% group_by(donor.id,source,cohort) %>% dplyr::count()
ggplot(counts,aes(donor.id,n,fill=cohort))+facet_wrap(~source)+
geom_col(position=position_dodge())+
labs(x="Donor",y="N")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
theme(axis.text.x=element_text(angle=90))

# harmony vs pcs
p1=DimPlot(all_combo,reduction="pca",group.by="source")
p2=DimPlot(all_combo,reduction="harmony",group.by="source")
p3=DimPlot(all_combo,reduction="pca",group.by="cohort")
p4=DimPlot(all_combo,reduction="harmony",group.by="cohort")

png(res=300,units="in",height=8,width=8,file="harmony_vs_pcs.png")
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

png(res=300,units="in",height=8,width=8,file="harmony_vs_pcs.png")
grid.arrange(p1,p2,nrow=1)
dev.off()

p1=DimPlot(all_combo,reduction="pca",group.by="batch_id")
p2=DimPlot(all_combo,reduction="harmony",group.by="batch_id")
png(res=300,units="in",height=8,width=16,file="harmony_vs_pcs_batch.png")
grid.arrange(p1,p2,nrow=1)
dev.off()

# Elbow plot
png(res=300,units="in",height=8,width=8,file="elbow.png")
ElbowPlot(all_combo,reduction="harmony")
dev.off()

# get var explained
eigenvalues = all_combo@reductions$harmony@stdev^2

centiles = c(50,60,70,80,90,90,99)
sapply(centiles,function(x){
  no_pcs = which(cumsum(eigenvalues) / sum(eigenvalues) * 100 > x)[1]
  return(no_pcs)
})

# see how PCs / Harmony PCs are influenced by batch and source
harmony_embeddings = data.frame(all_combo@reductions$harmony@cell.embeddings)
harmony_embeddings = harmony_embeddings %>% mutate(cell_name = rownames(harmony_embeddings))
pc_embeddings = data.frame(all_combo@reductions$pc@cell.embeddings)
pc_embeddings = pc_embeddings %>% mutate(cell_name = rownames(pc_embeddings))
meta_data = data.frame(all_combo@meta.data)
meta_data = meta_data %>% mutate(cell_name = rownames(meta_data))
joint_df = meta_data %>% left_join(harmony_embeddings,by="cell_name") %>% left_join(pc_embeddings,by="cell_name")

p1=ggplot(joint_df,aes(batch_id,harmony_1,fill=source))+geom_violin(alpha=0.7)+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_brewer(palette="Set2")
p2=ggplot(joint_df,aes(batch_id,PC_1,fill=source))+geom_violin(alpha=0.7)+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_brewer(palette="Set2")
png(res=300,units="in",height=8,width=8,file="harmony_vs_pcs_embeddings.png")
grid.arrange(p1,p2,ncol=1)
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
  png(res=300,units="in",height=8,width=12,file=paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,"blueprint.png"))
  print(DimPlot(all_combo_test,label=T,group.by="ann_blueprint"))
  dev.off()
  png(res=300,units="in",height=8,width=12,file=paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,"hpca.png"))
  print(DimPlot(all_combo_test,label=T,group.by="ann_hpca"))
  dev.off()
  png(res=300,units="in",height=8,width=12,file=paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,"monaco.png"))
  print(DimPlot(all_combo_test,label=T,group.by="ann_monaco"))
  dev.off()
  png(res=300,units="in",height=20,width=20,file=paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,"blueprint.png"))
  print( ggplot(all_combo_test@meta.data,aes(seurat_clusters,fill=ann_blueprint))+geom_bar(position="fill") )
  dev.off()

  saveRDS(all_combo_test,paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,".rds"))
}

mat = matrix(c(50,50,50,50,1.5,2,2.5,1),ncol=2)
mapply(umap_change_pcs,x=mat[,1],y=mat[,2]) # takes several hours - not routinely run
