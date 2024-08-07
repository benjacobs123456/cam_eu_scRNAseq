
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
data = readRDS(args[2]) #/xx/xxx/SC/CAM_TUM/datasets/all_combo_phenotypes_new_cluster_ids.rds

# define source
selsource = args[3] # CSF or PBMC
data = subset(data, subset = source == selsource)

# celltypes
cts = unique(data@meta.data$cell_type_crude)
Idents(data) = "cell_type_crude"

# Rerun SCtransform
batches = SplitObject(data,split.by="batch_id")
for (i in 1:length(batches)){
	print(i)
batches[[i]] = SCTransform(batches[[i]],assay="RNA",vars.to.regress="percent.mt",method = "glmGamPoi",variable.features.n = 10000,return.only.var.genes =F)
}
rm(data)
gc()
datam = merge(x = batches[[1]], y = batches[2:length(batches)])
scaledat = batches[[1]]@assays$SCT@scale.data
for (i in 2:length(batches)){
	print(i)
	scaledat = merge(scaledat,batches[[i]]@assays$SCT@scale.data,by=0,all=T)
	rownames(scaledat) = scaledat$Row.names; scaledat$Row.names = NULL
	
}
rm(batches); gc()
scaledatorig = datam@assays$SCT@data
origrownames = rownames(scaledatorig)
origcolnames = colnames(scaledatorig)
scaledat= scaledat[origrownames,origcolnames]
table(rownames(scaledatorig) == rownames(scaledat))
table(colnames(scaledatorig) == colnames(scaledat))

scaled_counts = datam@assays$SCT@scale.data
datam@assays$SCT@scale.data = as.matrix(scaledat)

saveRDS(datam,paste0(outputfolder,"all_combo_renormalizedperbatch.RDS"))

# Fix sample IDs
datam@meta.data[datam@meta.data$iid == "8A2H3PL1","iid"] = "PatID_5"
datam@meta.data[datam@meta.data$iid == "E9LEH7P8","iid"] = "PatID_56"
datam@meta.data[datam@meta.data$iid == "JEGK54J2","iid"] = "PatID_59"
datam@meta.data[datam@meta.data$iid == "RL4X6288","iid"] = "PatID_3"
datam@meta.data[datam@meta.data$iid == "WAF0FQN5","iid"] = "PatID_52"
datam@meta.data[datam@meta.data$iid == "WQ0EJVR1","iid"] = "PatID_60"

# convert to sce object
data.sce = as.SingleCellExperiment(datam, assay="SCT")

# GENE QC 1 remove genes that are expressed in < 2% of all cells
rawcounts = counts(data.sce)
keep <- rownames(rawcounts)[rowSums(rawcounts != 0) > ncol(rawcounts)*0.02]
scaledcounts = assay(data.sce,"scaledata")
scaledcounts = scaledcounts[rownames(scaledcounts) %in% keep == T,]

# find out which celltype*individual do not have enough cells
min_cells_per_sample = 5

low_counts = datam@meta.data %>%
  group_by(iid,cell_type_crude,phenotype) %>% 
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(group_to_exclude = paste0(cell_type_crude,"_",phenotype,"_",iid))

# average counts per individual and cell type
groups = colData(data.sce)[, c("cell_type_crude", "phenotype","iid")]
groups$BC = rownames(groups)
scaledcounts = t(scaledcounts); scaledcounts = data.frame(scaledcounts)
scaledcounts$BC = rownames(scaledcounts)
scaledcounts = join(scaledcounts,unique(data.frame(groups)),by="BC",type="left")
meanscaledcounts = scaledcounts %>% group_by(cell_type_crude,phenotype,iid) %>% summarise_at(vars(-BC), mean)

# Remove these combinations of celltype*individual
meanscaledcounts$comb = paste0(meanscaledcounts$cell_type_crude,"_",meanscaledcounts$phenotype,"_",meanscaledcounts$iid)
meanscaledcounts = meanscaledcounts[meanscaledcounts$comb %in% low_counts$group_to_exclude == F,]

# Get numbers (for later plots)
test1 = colData(data.sce)
groups1 = test1[paste0(test1$cell_type_crude,"_",test1$phenotype,"_",test1$iid) %in% low_counts$group_to_exclude == F, c("cell_type_crude", "phenotype","iid")]
groups1$comb = paste0(groups1$cell_type_crude,"_",groups1$phenotype,"_",groups1$iid)
nums1 = data.frame(table(groups1$comb))
write.table(nums1,paste0(outputfolder,"numbers_cells.txt"),sep="\t",row.names=F,quote=F)

# Gene QC2 - genes at least in 20% of individuals per cell type
rowsums = rowSums(rawcounts[rownames(rawcounts) %in% gsub("\\.","-",colnames(meanscaledcounts)) == T | rownames(rawcounts) %in% colnames(meanscaledcounts) == T,])
check = (rowsums != 0); keep = names(check[check == TRUE])
rawcountssel = rawcounts[keep,]
rawcountssel = t(rawcountssel); rawcountssel = data.frame(rawcountssel)
rawcountssel$BC = rownames(rawcountssel)
rawcountssel = join(rawcountssel,unique(data.frame(groups)),by="BC",type="left")
sumrawcountssel = rawcountssel %>% group_by(cell_type_crude,phenotype,iid) %>% summarise_at(vars(-BC), sum)
genes_all = c()
for (ct in cts){
	mat = sumrawcountssel[sumrawcountssel$cell_type_crude == ct,]
	mat$cell_type_crude = NULL; mat$iid = NULL; mat$phenotype = NULL
	mat = t(mat)
    genes_ct <- rownames(mat)[rowSums(mat != 0) >= ncol(mat)*0.2]
    genes_all = unique(c(genes_all,genes_ct))
	write.table(genes_ct,paste0(outputfolder,"genes_",ct,".txt"),sep="\t",row.names=F,quote=F,col.names=F)
}
length(genes_all)
meanscaledcountsout = data.frame(meanscaledcounts)
rownames(meanscaledcountsout) = paste0(meanscaledcountsout$cell_type_crude,"_",meanscaledcountsout$phenotype,"_",meanscaledcountsout$iid)
meanscaledcountsout$cell_type_crude = NULL; meanscaledcountsout$iid = NULL; meanscaledcountsout$phenotype = NULL
meanscaledcountsout = data.frame(t(meanscaledcountsout))
meanscaledcountsout = meanscaledcountsout[rownames(meanscaledcountsout) %in% genes_all,]


# Write aggregated gene expression data to file
write.table(meanscaledcountsout,paste0(outputfolder,"data_gex.txt"),sep="\t")



