# conda activate dandelion

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

# read in B cell GEX data
b_cells = readRDS("./datasets/b_cells.rds")


##############################
### Combine GEX with VDJ   ###
##############################

# load BCR VDJ files post-processing (by singularity dandelion)
samples = list.files("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/dandelion_inputs/BCR/")
samples = samples[grepl("CSF",samples)|grepl("PBMC",samples)]

files = list()
for (i in 1:length(samples)){
    filename = paste0("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/dandelion_inputs/BCR/",samples[i],"/dandelion/filtered_contig_igblast_db-pass_genotyped.tsv")
    files[[i]] = read_tsv(filename)
}
combined_bcr = do.call(rbind, files)

# create a column called filter_rna and set it to FALSE
b_cells@meta.data$filter_rna = FALSE

# revert to old barcodes (get rid of extra bits added due to duplicates)
original_barcodes = str_split_fixed(rownames(b_cells@meta.data),pattern="-",n=2)[,1]
b_cells@meta.data$original_barcode = original_barcodes

# add in source and donor id
b_cells@meta.data$original_barcode = paste0(b_cells@meta.data$source,"_",b_cells@meta.data$donor.id,"_",b_cells@meta.data$original_barcode)

# rename cells with new barcodes
b_cells = RenameCells(b_cells,new.names=b_cells@meta.data$original_barcode)

# change to python object
sc = import("scanpy")
DefaultAssay(b_cells)="RNA"
# convert the meta.data slot to a python friendly object
obs = r_to_py(b_cells@meta.data)
normcounts = r_to_py(Matrix::t(GetAssayData(b_cells)))
adata = sc$AnnData(X = normcounts, obs = obs)
sc$pp$neighbors(adata)

# Combine BCR data and GEX
vdj_results_list = ddl$pp$filter_contigs(combined_bcr,
  adata,
  filter_contig=TRUE,
  filter_rna=TRUE,
  filter_poorqualitycontig=TRUE,
  filter_vj_chains=T,
  filter_missing=T,
  filter_poorqualitycontig=T,
  productive_only = T,
  keep_highest_umi=TRUE)

# Print no. of cells passing QC
vdj_results_list[[2]]$obs %>% filter(
has_contig == "True" &
contig_QC_pass == "True" &
filter_contig == "FALSE") %>% nrow()

ann_df = vdj_results_list[[2]]$obs %>% data.frame
qc_stats = ann_df %>% dplyr::count(has_contig,contig_QC_pass,filter_contig,filter_contig_quality,filter_contig_VJ,filter_contig_VDJ)
write_csv(qc_stats,"./bcr/vdj_qc_stats.csv")

# filter cells which don't have 2 contigs exactly
bad_cells = vdj_results_list[[1]]$data %>% group_by(cell_id) %>% dplyr::count() %>% filter(n!=2)
vdj_results_list[[1]]$data = vdj_results_list[[1]]$data %>% filter(!cell_id %in% bad_cells$cell_id)

#################################
### Clones & network          ###
#################################
# find clones
ddl$tl$find_clones(vdj_results_list[[1]])
# find clone size
ddl$tl$clone_size(vdj_results_list[[1]])



# generate BCR network
# ddl$tl$generate_network(vdj_results_list[[1]])

# get back into seurat to visualise
add_col_to_metadata = function(x){
  col_to_add = unlist(vdj_results_list[[1]]$metadata[[x]])
  names(col_to_add) <- row.names(vdj_results_list[[1]]$metadata)
  b_cells <<- AddMetaData(b_cells, col_to_add, x)
}
sapply(colnames(vdj_results_list[[1]]$metadata),add_col_to_metadata)
# get rid of cells without BCR data
cells_in_vdj = rownames(vdj_results_list[[1]]$metadata)
b_cells = subset(b_cells,subset = original_barcode %in% cells_in_vdj)
# rename Idents
b_cells@meta.data$cell_type = Idents(b_cells)

#################################
### SHM                       ###
#################################
# calculate mutations in Shazam
new_airr = shazam::observedMutations(vdj_results_list[[1]]$data,frequency=TRUE,regionDefinition=IMGT_V)
# calculate physicochemical propoerties in alakazam
new_airr = alakazam::aminoAcidProperties(new_airr,seq="junction",trim=TRUE)


# function to add cols from airr/dandelion object into b cell metadata
# first just add col in metadata to match name
b_cells@meta.data$cell_id = b_cells@meta.data$original_barcode

cols_to_include = c("mu_freq_cdr_r","mu_freq_cdr_s","mu_freq_fwr_r",
"mu_freq_fwr_s","junction_aa_length",
"junction_aa_gravy","junction_aa_bulk",
"junction_aa_aliphatic","junction_aa_polarity",
"junction_aa_charge","junction_aa_basic",
"junction_aa_acidic","junction_aa_aromatic")

add_cols_from_airr_to_metadata = function(airr_table,metadata,included_cols=cols_to_include){
  airr_table = airr_table %>% filter(cell_id %in% metadata[['cell_id']])
  airr_table = airr_table %>% select(cell_id,locus,cols_to_include)
  heavy_chains = airr_table %>% filter(locus=="IGH")
  light_chains = airr_table %>% filter(locus=="IGL" | locus =="IGK")

  colnames(heavy_chains)[colnames(heavy_chains)!="cell_id"] = paste0("heavychain_",colnames(heavy_chains)[colnames(heavy_chains)!="cell_id"])
  colnames(light_chains)[colnames(light_chains)!="cell_id"] = paste0("lightchain_",colnames(light_chains)[colnames(light_chains)!="cell_id"])

  metadata = metadata %>% left_join(heavy_chains,by="cell_id") %>%
  left_join(light_chains,by="cell_id")
  metadata
}
b_cells@meta.data = add_cols_from_airr_to_metadata(new_airr,b_cells@meta.data)
# rename Idents
b_cells@meta.data$cell_type = Idents(b_cells)
library(gridExtra)
b_cells@meta.data$r_s_ratio = b_cells@meta.data$heavychain_mu_freq_cdr_r / b_cells@meta.data$heavychain_mu_freq_cdr_s
b_cells@meta.data = b_cells@meta.data  %>% mutate(cdr3_mu_freq_overall = (heavychain_mu_freq_cdr_r+heavychain_mu_freq_cdr_s+lightchain_mu_freq_cdr_r+lightchain_mu_freq_cdr_s)/2)
b_cells@meta.data = b_cells@meta.data  %>% mutate(shm_positive = ifelse(cdr3_mu_freq_overall==0,"No","Yes"))

# save progress
saveRDS(b_cells,"b_cells_post_processing.rds")

# diversity
curve = estimateAbundance(new_airr, group="sample_id", ci=0.95, nboot=100, clone="clone_id", min_n=20)

sample_curve = alphaDiversity(new_airr, group="sample_id", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=10)

unique_donors = b_cells@meta.data %>% distinct(donor.id,.keep_all=TRUE)
cc_status = unique_donors %>% select(donor.id,phenotype)
diversity = sample_curve@diversity %>% separate(sample_id,sep="_",into=c("source","donor.id"))
diversity = cc_status %>% left_join(diversity,by="donor.id")

q2 = diversity %>% filter(q==2) %>% filter(phenotype=="MS")
p2 = ggplot(q2,aes(source,d,fill=source))+geom_boxplot()+ggtitle("Simpson's Diversity Index: MS CSF B cells vs PBMC")+theme_classic()

png(res=300,units="in",height=6,width=6,"./bcr/diversity.png")
p2
dev.off()
