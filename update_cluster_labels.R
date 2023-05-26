# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)

# set WD
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/datasets")

# plot consensus celltypist annotations
celltypist_annotations = list()
for(j in c(1:66)){
  i=j-1
  message("Cluster: ",i)
  hires = read_csv(paste0("highres_cluster",i,"predicted_labels.csv"))
  hires$cluster = i
  celltypist_annotations[[j]] = hires
}
celltypist_annotations = do.call("bind_rows",celltypist_annotations)

plot_dat = celltypist_annotations %>%
  group_by(cluster) %>%
  dplyr::count(predicted_labels) %>%
  mutate(prop = n/sum(n)) %>%
  filter(prop >0.05) %>%
  slice_max(prop,n=1)

png("celltypist_annotations.png",res=300,units="in",width=6,height=10)
ggplot(plot_dat,aes(cluster,prop,label=paste0(cluster,":",predicted_labels)))+
geom_col(fill="white",color="black")+
theme_minimal()+
geom_text(hjust=1)+
coord_flip()
dev.off()

write_csv(plot_dat,"celltypist_cluster_calls.csv")

# repeat with lowres


# plot consensus celltypist annotations
celltypist_annotations = list()
for(j in c(1:66)){
  i=j-1
  message("Cluster: ",i)
  lores = read_csv(paste0("lowres_cluster",i,"predicted_labels.csv"))
  lores$cluster = i
  celltypist_annotations[[j]] = lores
}
celltypist_annotations = do.call("bind_rows",celltypist_annotations)


plot_dat = celltypist_annotations %>%
  group_by(cluster) %>%
  dplyr::count(predicted_labels) %>%
  mutate(prop = n/sum(n)) %>%
  filter(prop >0.05) %>%
  slice_max(prop,n=1)

png("celltypist_annotations_lowres.png",res=300,units="in",width=6,height=10)
ggplot(plot_dat,aes(cluster,prop,label=paste0(cluster,":",predicted_labels)))+
geom_col(fill="white",color="black")+
theme_minimal()+
geom_text(hjust=1)+
coord_flip()
dev.off()

write_csv(plot_dat,"celltypist_cluster_calls_lowres.csv")


# read in outputs from deconvolution & integration step
message("Reading in data")
all_combo = readRDS("all_combo_with_pheno.rds")

# check number of clusters
n_clust = length(unique(all_combo@meta.data$seurat_clusters))
message("No. of clusters: ",n_clust)

# new manual labelling of clusters IDs
message("Renaming clusters")
new_cluster_ids = read_csv("updated_cluster_identities.csv")
all_combo@meta.data = all_combo@meta.data %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  left_join(new_cluster_ids %>%
    dplyr::rename("seurat_clusters" = cluster) %>%
    mutate(seurat_clusters = as.character(seurat_clusters)),
    by="seurat_clusters")
rownames(all_combo@meta.data) = colnames(all_combo)

new_cluster_ids = new_cluster_ids %>%
  arrange(as.character(cluster))
new_cluster_ids = new_cluster_ids$cell_type
names(new_cluster_ids) = levels(all_combo)
all_combo = RenameIdents(all_combo,new_cluster_ids)

# read in celltypist annotations
all_clusters = list()
for(j in c(1:66)){
  i=j-1
  message("Doing cluster ",i)
  this_cluster = all_combo@meta.data %>%
    filter(seurat_clusters == i)
  highres = read_csv(paste0("highres_cluster",i,"predicted_labels.csv"))
  lowres = read_csv(paste0("lowres_cluster",i,"predicted_labels.csv"))
  this_cluster$ann_celltypist_highres = highres$predicted_labels
  this_cluster$ann_celltypist_lowres = lowres$predicted_labels
  all_clusters[[j]] = this_cluster
}
all_annotations = do.call("bind_rows",all_clusters)
all_annotations$original_cell_id = rownames(all_annotations)

# merge back together
all_combo@meta.data$original_cell_id = rownames(all_combo@meta.data)

all_combo@meta.data = all_combo@meta.data %>%
  left_join(all_annotations %>%
    dplyr::select(original_cell_id,ann_celltypist_lowres,ann_celltypist_highres),
    by="original_cell_id")
rownames(all_combo@meta.data) = colnames(all_combo)

# process biomarkers for supp table
biomarkers = read_csv("biomarker_summary.csv")
new_cluster_ids = data.frame(new_cluster_ids)
new_cluster_ids$cluster = rownames(new_cluster_ids)

new_cluster_ids = new_cluster_ids %>% left_join(
biomarkers %>%
  dplyr::rename("cluster" = x) %>%
  mutate(cluster = as.character(cluster)), by = "cluster"
)
write_csv(new_cluster_ids,"supp_table_4_clusters.csv")

# save
saveRDS(all_combo,"all_combo_phenotypes_new_cluster_ids.rds")
