##############################
###     Load packages      ###
##############################
# NB first run conda activate dandelion
library(reticulate)
use_condaenv("dandelion",required=TRUE)
library(Seurat)
library(dplyr)
library(readr)
ddl = import('dandelion')
library(stringr)
library(harmony)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(celldex)
library(SingleR)
library(edgeR)
library(tidyr)
library(shazam)
library(alakazam)

sessionInfo()
reticulate::py_config()

setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined")

# read in T cell GEX data
t_cells = readRDS("./datasets/t_cells.rds")


##############################
### Combine GEX with VDJ   ###
##############################

# load tcr VDJ files post-processing (by singularity dandelion)
samples = list.files("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/dandelion_inputs/TCR/")
samples = samples[grepl("CSF",samples)|grepl("PBMC",samples)]

files = list()
for (i in 1:length(samples)){
    filename = paste0("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/dandelion_inputs/TCR/",samples[i],"/dandelion/filtered_contig_igblast_db-pass.tsv")
    tcrs = read_tsv(filename)
    tcrs$cell_id = paste0(samples[i],"_",tcrs$cell_id)
    tcrs$sequence_id = paste0(samples[i],"_",tcrs$sequence_id)
    files[[i]] = tcrs
}
combined_tcr = do.call(rbind, files)

# create a column called filter_rna and set it to FALSE
t_cells@meta.data$filter_rna = FALSE

# revert to old barcodes (get rid of extra bits added due to duplicates)
original_barcodes = str_split_fixed(rownames(t_cells@meta.data),pattern="_",n=2)[,1]
t_cells@meta.data$original_barcode = original_barcodes

# add in source and donor id
t_cells@meta.data$original_barcode = paste0(t_cells@meta.data$source,"_",t_cells@meta.data$donor.id,"_",t_cells@meta.data$original_barcode)

# remove duplicate names (n=10)
unique_cells = t_cells@meta.data %>% group_by(original_barcode) %>% dplyr::count() %>% filter(n==1)
t_cells = subset(t_cells, subset = original_barcode %in% unique_cells$original_barcode)

# rename cells with new barcodes
t_cells = RenameCells(t_cells,new.names=t_cells@meta.data$original_barcode)

# filter VDJ and seurat objects so they contain the same cells
combined_tcr = combined_tcr %>% filter(cell_id %in% t_cells@meta.data$original_barcode)
t_cells = subset(t_cells, subset = original_barcode %in% combined_tcr$cell_id)

# change to python object
sc = import("scanpy")
DefaultAssay(t_cells)="RNA"
# convert the meta.data slot to a python friendly object
obs = r_to_py(t_cells@meta.data)
normcounts = r_to_py(Matrix::t(GetAssayData(t_cells)))
adata = sc$AnnData(X = normcounts, obs = obs)
sc$pp$neighbors(adata)

# Combine tcr data and GEX
vdj_results_list = ddl$pp$filter_contigs(combined_tcr,
  adata,
  filter_contig=T,
  filter_rna=T,
  filter_poorqualitycontig=T,
  filter_vj_chains=F,
  filter_missing=T,
  productive_only=T,
  keep_highest_umi=TRUE)

# Print no. of cells passing QC
vdj_results_list[[2]]$obs %>% filter(
has_contig == "True" &
contig_QC_pass == "True" &
filter_contig == "FALSE") %>% nrow()

ann_df = vdj_results_list[[2]]$obs %>% data.frame
qc_stats = ann_df %>% dplyr::count(has_contig,contig_QC_pass,filter_contig,filter_contig_quality,filter_contig_VJ,filter_contig_VDJ)
write_csv(qc_stats,"./tcr/vdj_qc_stats.csv")

#################################
### Clones & network          ###
#################################
# find clones
ddl$tl$find_clones(self=vdj_results_list[[1]],locus='tr')
# find clone size
ddl$tl$clone_size(vdj_results_list[[1]])

# generate tcr network
# ddl$tl$generate_network(vdj_results_list[[1]])

# get back into seurat to visualise
add_col_to_metadata = function(x){
  col_to_add = unlist(vdj_results_list[[1]]$metadata[[x]])
  names(col_to_add) <- row.names(vdj_results_list[[1]]$metadata)
  t_cells <<- AddMetaData(t_cells, col_to_add, x)
}
sapply(colnames(vdj_results_list[[1]]$metadata),add_col_to_metadata)

# get rid of cells without tcr data
cells_in_vdj = rownames(vdj_results_list[[1]]$metadata)
t_cells = subset(t_cells,subset = original_barcode %in% cells_in_vdj)
# rename Idents
t_cells@meta.data$cell_type = Idents(t_cells)


# save progress
saveRDS(t_cells,"t_cells_post_processing.rds")
message("saved progress")
