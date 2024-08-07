##############################
###     load packages      ###
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

setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/")

# read in T cell GEX data
t_cells = readRDS("../t_cells.rds")


##############################
### Combine GEX with VDJ   ###
##############################

# load tcr VDJ files post-processing (by singularity dandelion)
meta = read_csv("../dandelion_inputs/TCR/meta_file.csv")

files = list()
for (i in 1:nrow(meta)){
    this_sample = meta$sample[i]
    this_person = meta$individual[i]
    file_path = paste0("../dandelion_inputs/TCR/",this_sample,"/dandelion/filtered_contig_dandelion.tsv")
    if(file.exists(file_path)){
      dat = read_tsv(file_path, col_types = cols_only(
      sequence_id = col_character(),
      sequence = col_character(),
      rev_comp = col_logical(),
      productive = col_logical(),
      v_call = col_character(),
      d_call = col_character(),
      j_call = col_character(),
      sequence_alignment = col_character(),
      germline_alignment = col_character(),
      junction = col_character(),
      junction_aa = col_character(),
      v_cigar = col_character(),
      d_cigar = col_character(),
      j_cigar = col_character(),
      stop_codon = col_logical(),
      vj_in_frame = col_logical(),
      locus = col_character(),
      c_call = col_character(),
      junction_length = col_double(),
      np1_length = col_double(),
      np2_length = col_double(),
      v_sequence_start = col_double(),
      v_sequence_end = col_double(),
      v_germline_start = col_double(),
      v_germline_end = col_double(),
      d_sequence_start = col_double(),
      d_sequence_end = col_double(),
      d_germline_start = col_double(),
      d_germline_end = col_double(),
      j_sequence_start = col_double(),
      j_sequence_end = col_double(),
      j_germline_start = col_double(),
      j_germline_end = col_double(),
      v_score = col_double(),
      v_identity = col_double(),
      v_support = col_double(),
      d_score = col_double(),
      d_identity = col_double(),
      d_support = col_double(),
      j_score = col_double(),
      j_identity = col_double(),
      j_support = col_double(),
      fwr1 = col_character(),
      fwr2 = col_character(),
      fwr3 = col_character(),
      fwr4 = col_character(),
      cdr1 = col_character(),
      cdr2 = col_character(),
      cdr3 = col_character(),
      cell_id = col_character(),
      consensus_count = col_double(),
      duplicate_count = col_double(),
      v_call_10x = col_character(),
      d_call_10x = col_character(),
      j_call_10x = col_character(),
      junction_10x = col_character(),
      junction_10x_aa = col_character(),
      j_support_igblastn = col_double(),
      j_score_igblastn = col_double(),
      j_call_igblastn = col_character(),
      j_call_blastn = col_character(),
      j_identity_blastn = col_double(),
      j_alignment_length_blastn = col_double(),
      j_number_of_mismatches_blastn = col_double(),
      j_number_of_gap_openings_blastn = col_double(),
      j_sequence_start_blastn = col_double(),
      j_sequence_end_blastn = col_double(),
      j_germline_start_blastn = col_double(),
      j_germline_end_blastn = col_double(),
      j_support_blastn = col_double(),
      j_score_blastn = col_double(),
      j_sequence_alignment_blastn = col_character(),
      j_germline_alignment_blastn = col_character(),
      j_source = col_logical(),
      d_support_igblastn = col_double(),
      d_score_igblastn = col_double(),
      d_call_igblastn = col_character(),
      d_call_blastn = col_character(),
      d_identity_blastn = col_double(),
      d_alignment_length_blastn = col_double(),
      d_number_of_mismatches_blastn = col_double(),
      d_number_of_gap_openings_blastn = col_double(),
      d_sequence_start_blastn = col_double(),
      d_sequence_end_blastn = col_double(),
      d_germline_start_blastn = col_double(),
      d_germline_end_blastn = col_double(),
      d_support_blastn = col_double(),
      d_score_blastn = col_double(),
      d_sequence_alignment_blastn = col_character(),
      d_germline_alignment_blastn = col_character(),
      d_source = col_character(),
      junction_aa_length = col_double(),
      fwr1_aa = col_character(),
      fwr2_aa = col_character(),
      fwr3_aa = col_character(),
      fwr4_aa = col_character(),
      cdr1_aa = col_character(),
      cdr2_aa = col_character(),
      cdr3_aa = col_character(),
      sequence_alignment_aa = col_character(),
      v_sequence_alignment_aa = col_character(),
      d_sequence_alignment_aa = col_character(),
      j_sequence_alignment_aa = col_character(),
      j_call_multimappers = col_character(),
      j_call_multiplicity = col_double(),
      j_call_sequence_start_multimappers = col_double(),
      j_call_sequence_end_multimappers = col_double(),
      j_call_support_multimappers = col_double()
    ))

      dat$sequence_id = paste0(this_person,"_",str_remove(dat$sequence_id,paste0(this_sample,"_")))
      dat$cell_id = paste0(this_person,"_",str_remove(dat$cell_id,paste0(this_sample,"_")))
      files[[i]] = dat
    } else {
      message(file_path," does not exist!")
    }
}
combined_tcr = do.call("bind_rows", files)

# find duplicate contigs (keep higher UMI)
combined_tcr = combined_tcr %>%
arrange(desc(consensus_count)) %>%
distinct(sequence_id,.keep_all=T)

# create a column called filter_rna and set it to FALSE
t_cells@meta.data$filter_rna = FALSE

# revert to old barcodes (get rid of extra bits added due to duplicates)
t_cells@meta.data$original_barcode = rownames(t_cells@meta.data)

get_orig_barcode = function(x){
y = str_split_fixed(x,n=2,"_")[1,1]
return(y)
}
t_cells@meta.data$original_barcode = sapply(t_cells@meta.data$original_barcode,get_orig_barcode)

# add in donor id
t_cells@meta.data$original_barcode = paste0(t_cells@meta.data$iid,"_",t_cells@meta.data$original_barcode)

# find duplicates
non_duplicated_barcodes = t_cells@meta.data %>% dplyr::count(original_barcode) %>% filter(n==1)
combined_tcr = combined_tcr %>% filter(cell_id %in% non_duplicated_barcodes$original_barcode)

# filter to cells in VDJ
t_cells = subset(t_cells, original_barcode %in% combined_tcr$cell_id & original_barcode %in% non_duplicated_barcodes$original_barcode)

# filter VDJ to contigs with GEX data
combined_tcr = combined_tcr %>% filter(cell_id %in% t_cells@meta.data$original_barcode)
message("There are", nrow(combined_tcr)," contigs")
message("There are", nrow(combined_tcr %>% distinct(cell_id))," unique cells in TCR data")

# rename cells with new barcodes
t_cells = RenameCells(t_cells,new.names=t_cells@meta.data$original_barcode)
message("There are", nrow(t_cells@meta.data %>% nrow)," unique cells in gex data")


# change to python object
sc = import("scanpy")
DefaultAssay(t_cells)="RNA"
# convert the meta.data slot to a python friendly object
obs = r_to_py(t_cells@meta.data)
normcounts = r_to_py(Matrix::t(GetAssayData(t_cells)))
adata = sc$AnnData(X = normcounts, obs = obs)
sc$pp$neighbors(adata)
message("Finished making python object")

# Combine tcr data and GEX
vdj_results_list = ddl$pp$filter_contigs(combined_tcr,
  adata,
  filter_contig=T,
  filter_rna=T,
  filter_poorqualitycontig=T,
  filter_vj_chains=F,
  filter_missing=T,
  productive_only=T,
  keep_highest_umi=T)

# Print no. of cells passing QC
vdj_results_list[[2]]$obs %>% filter(
has_contig == "True" &
contig_QC_pass == "True" &
filter_contig == "FALSE") %>% nrow()

ann_df = vdj_results_list[[2]]$obs %>% data.frame
qc_stats = ann_df %>% dplyr::count(has_contig,contig_QC_pass,filter_contig,filter_contig_quality,filter_contig_VJ,filter_contig_VDJ)
write_csv(qc_stats,"vdj_qc_stats.csv")

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
