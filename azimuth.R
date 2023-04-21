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


# sample 100,000 cells
cells = colnames(all_combo) %>% sample(size=40000)

all_combo[['cell_id']] = colnames(all_combo)
all_combo@meta.data = all_combo@meta.data %>% mutate(keep = ifelse(cell_id %in% cells,"yes","no"))
all_combo = subset(all_combo,subset = keep == "yes")

# remove SCT
DefaultAssay(all_combo) = "RNA"
all_combo[['SCT']] = NULL
all_combo_azimuth = DietSeurat(all_combo)
saveRDS(all_combo_azimuth,"trimmed_for_azimuth.rds")

# process azimuth results
all_combo = readRDS("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/datasets/trimmed_for_azimuth.rds")
preds = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/datasets/azimuth_preds.tsv")
preds = preds %>% filter(mapping.score>0.5)
all_combo@meta.data = all_combo@meta.data %>% filter(cell_id %in% preds$cell) %>% left_join(preds %>% dplyr::rename("cell_id" = "cell"))

# look at azimuth annotations
totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(predicted.celltype.l2) %>% summarise(total = sum(n))
all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(predicted.celltype.l2) %>% left_join(totals,by="seurat_clusters")

azimuth_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(predicted.celltype.l2) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(predicted.celltype.l2))%>% slice_max(n=1,order_by=prop)
write_csv(azimuth_calls,"azimuth_cluster_calls.csv")
p=ggplot(azimuth_calls,aes(seurat_clusters,prop,fill=predicted.celltype.l2,label=predicted.celltype.l2))+geom_col()+geom_text(angle=90)+theme_bw()
png("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/cluster_plots/azimuth_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()
