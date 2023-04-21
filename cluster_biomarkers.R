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
all_combo = readRDS("all_combo_with_UMAP_PCs_50resolution2.5.rds")

# clean and integrated phenotype data
reference = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/EU_ID.txt")
reference$cohort = "EU"
reference$PBMC_PoolSize = as.character(reference$PBMC_PoolSize)
reference_cam = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/GT_ID.txt")
reference_cam$cohort = "Cam"
reference_cam$PBMC_PoolSize = as.character(reference_cam$PBMC_PoolSize)
reference = bind_rows(reference,reference_cam)

# sort out long cambridge IDs
new_ids = lapply(all_combo@meta.data$donor.id, function(x){
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

all_combo[['donor.id']] = unlist(new_ids)

# get rid of donors without an ID
unique_donors = unique(all_combo@meta.data$donor.id)
assigned_donors = unique_donors[!grepl("donor",unique_donors)]
all_combo = subset(all_combo, subset = donor.id %in% assigned_donors )

# merge with phenotype file
cam_pheno = read_csv("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/references/5PrimeCSF_ClinicalData_NoPID_211124.csv")
cam_pheno$shortID = unlist(lapply(cam_pheno$ID, function(x){
  new_id = str_remove(pattern="C00TU0",x)
  new_id = str_split(pattern="a",new_id)[[1]][1]
  new_id = paste0("TU",new_id)
  new_id = str_remove(pattern="v1",new_id)
  new_id = str_remove(pattern="v2",new_id)
  return(new_id)
}))
cam_pheno = cam_pheno %>% select(shortID,Category)
eu_pheno = read_csv("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/references/EU_pheno.csv")
pheno = bind_rows(cam_pheno,eu_pheno)
pheno$donor.id = pheno$shortID

# merge with metadata
donor_phenotypes = lapply(all_combo@meta.data$donor.id,function(x){
  message("processing ",x)
  cc_status = pheno[pheno$donor.id==x,]$Category[1]
  return(cc_status)
}) %>% unlist()

all_combo[['phenotype']] = donor_phenotypes

# remove singletons
cluster_counts = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count() %>% arrange(as.numeric(as.character(seurat_clusters)))
write_csv(cluster_counts,"raw_cluster_counts.csv")
all_combo = subset(all_combo, subset = seurat_clusters != "singleton")

# remove clusters with n<5 cells (1 clusters, 2 cells)
high_count_clusters = cluster_counts %>% filter(n>5)
all_combo = subset(all_combo, subset = seurat_clusters %in% high_count_clusters$seurat_clusters)

# remove patient on tysabri - 917 cells, all PBMC
all_combo = subset(all_combo, subset = donor.id != "TU146")

# make dim plots
p0 = DimPlot(all_combo,label=T, repel=T)
png("../cluster_plots/basic_dimplot_labelled.png",res=300,units="in",width=8,height=8)
p0
dev.off()

p1 = DimPlot(all_combo,split.by="source", repel=T)
png("../cluster_plots/basic_dimplot_by_source.png",res=300,units="in",width=8,height=8)
p1
dev.off()

p2 = DimPlot(all_combo,split.by="phenotype", repel=T)
png("../cluster_plots/basic_dimplot_by_category.png",res=300,units="in",width=8,height=8)
p2
dev.off()

p3 = DimPlot(all_combo,split.by="cohort", repel=T)
png("../cluster_plots/basic_dimplot_by_cohort.png",res=300,units="in",width=8,height=8)
p3
dev.off()

cluster_markers = c("CLEC10A","FCER1A","S100A8", "S100A9", "CD14","CD1C","FCGR3A", "SMIM25","CD3E","CD3D","CD8A", "CD8B","ANXA1","CCL5","CCR7","MKI67","TCF7","IL7R","TRDC","TRGC1","KLRB1", "GNLY", "NKG7","PRSS57","CD79A", "MS4A1","sct_IGKC","LILRA4", "IRF8", "PF4","FOXP3","CTLA4","CD38","CD27","IGHM","IGHG1","CD24")
p4 = FeaturePlot(all_combo,features=cluster_markers, repel=T)
png("../cluster_plots/canonical_markers.png",res=300,units="in",width=16,height=16)
p4
dev.off()

p5 = DotPlot(all_combo,features=cluster_markers) + theme(axis.text.x = element_text(angle=90))
png("../cluster_plots/dotplot_canonical_markers.png",res=300,units="in",width=12,height=8)
p5
dev.off()

short_cluster_markers = c("CLEC10A","S100A8", "CD16","CD3E","CD8A","IL7R", "NKG7","CD79A", "PF4","CD38")
p5 = FeaturePlot(all_combo,features=short_cluster_markers,raster=FALSE)
png("../cluster_plots/few_canonical_markers.png",res=300,units="in",width=12,height=8)
p5
dev.off()

# look at singleR annotations
totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_blueprint) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_blueprint) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_blueprint))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"singler_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_blueprint,label=ann_blueprint))+geom_col()+geom_text(angle=90)+theme_bw()
png("../cluster_plots/singler_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()

totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_monaco) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_monaco) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_monaco))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"monaco_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_monaco,label=ann_monaco))+geom_col()+geom_text(angle=90)+theme_bw()
png("../cluster_plots/monaco_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()

totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_dice) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_dice) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_dice))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"dice_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_dice,label=ann_dice))+geom_col()+geom_text(angle=90)+theme_bw()
png("../cluster_plots/dice_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()

####################################
# export for celltypist
####################################

# export data for celltypist
count_mat = all_combo@assays$RNA@counts %>% t()
for(i in c(1:7)){
  outfile = paste0("all_combo_celltypist_counts_chunk_",i,".tsv")
  start_row = 1 + (i-1) * 50000
  end_row = 50000 * i
  if(i < 7){
    filtered_count_mat = count_mat[c(start_row:end_row),]
    df = data.frame(filtered_count_mat)
    write_tsv(df,file = outfile )
  } else {
    filtered_count_mat = count_mat[c(start_row:nrow(count_mat)),]
    df = data.frame(filtered_count_mat)
    write_tsv(df,file = outfile )
  }
}

# NB must now run sbatch celltypist.sh

# read in celltypist preds
lowres_preds_list = list()
highres_preds_list = list()
for(i in c(1:7)){
  # low res
  filename = paste0("lowres_chunk",i,"predicted_labels.csv")
  preds = read_csv(filename)
  lowres_preds_list[[i]] = preds

  # high res
  filename = paste0("highres_chunk",i,"predicted_labels.csv")
  preds = read_csv(filename)
  highres_preds_list[[i]] = preds

}
lowres_preds_list = do.call("bind_rows",lowres_preds_list)
highres_preds_list = do.call("bind_rows",highres_preds_list)

# add to main Seurat object metadata
all_combo[['ann_celltypist_highres']] = highres_preds_list$predicted_labels


####################################
# cluster biomarkers
####################################
# Find markers
all_combo_markers = FindAllMarkers(all_combo, min.pct=0.25, logfc.threshold = 0.25,only.pos=TRUE,recorrect_umi=FALSE)
write_csv(all_combo_markers,"cluster_biomarkers.csv")
all_combo_markers = read_csv("cluster_biomarkers.csv")


# new manual labelling of clusters IDs
new_cluster_ids = read_csv("cluster_identities.csv")
new_cluster_ids = new_cluster_ids %>% arrange(as.character(cluster))

new_cluster_ids = new_cluster_ids$cell_type
names(new_cluster_ids) = levels(all_combo)
all_combo = RenameIdents(all_combo,new_cluster_ids)

# summarise biomarkers
biomarker_summary = all_combo_markers
biomarkers = lapply(unique(biomarker_summary$cluster),function(x){
  this_cluster = biomarker_summary %>% filter(cluster==x) %>% arrange(desc(avg_log2FC)) %>% head(10)
  genes = paste0(this_cluster$gene,collapse=", ")
  data.frame(x,genes)
})
biomarkers = do.call("bind_rows",biomarkers) %>% tibble::as_tibble()
write_csv(biomarkers,"biomarker_summary.csv")

# save
saveRDS(all_combo,"all_combo_phenotypes_new_cluster_ids.rds")
