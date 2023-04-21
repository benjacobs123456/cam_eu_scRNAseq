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

setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/BCR")

# read in B cell GEX data
b_cells = readRDS("../b_cells.rds")


##############################
### Combine GEX with VDJ   ###
##############################

# load BCR VDJ files post-processing (by singularity dandelion)
meta = read_csv("../dandelion_inputs/BCR/meta_file.csv")

files = list()
for (i in 1:nrow(meta)){
    this_sample = meta$sample[i]
    this_person = meta$individual[i]
    file_path = paste0("../dandelion_inputs/BCR/",this_sample,"/dandelion/filtered_contig_dandelion.tsv")
    dat = read_tsv(file_path,col_types = cols_only(
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
        d_source = col_logical(),
        v_call_genotyped = col_character(),
        germline_alignment_d_mask = col_character(),
        sample_id = col_character(),
        c_sequence_alignment = col_character(),
        c_germline_alignment = col_character(),
        c_sequence_start = col_double(),
        c_sequence_end = col_double(),
        c_score = col_double(),
        c_identity = col_double(),
        c_call_10x = col_character(),
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
        complete_vdj = col_logical(),
        j_call_multimappers = col_character(),
        j_call_multiplicity = col_double(),
        j_call_sequence_start_multimappers = col_double(),
        j_call_sequence_end_multimappers = col_double(),
        j_call_support_multimappers = col_double(),
        mu_count = col_double(),
        mu_freq = col_double()
      )
    )

    dat$sequence_id = paste0(this_person,"_",str_remove(dat$sequence_id,paste0(this_sample,"_")))
    dat$cell_id = paste0(this_person,"_",str_remove(dat$cell_id,paste0(this_sample,"_")))
    files[[i]] = dat
}
combined_bcr = do.call("bind_rows", files)

# create a column called filter_rna and set it to FALSE
b_cells@meta.data$filter_rna = FALSE

# revert to old barcodes (get rid of extra bits added due to duplicates)
b_cells@meta.data$original_barcode = rownames(b_cells@meta.data)

get_orig_barcode = function(x){
y = str_split_fixed(x,n=2,"_")[1,1]
return(y)
}
b_cells@meta.data$original_barcode = sapply(b_cells@meta.data$original_barcode,get_orig_barcode)

# add in donor id
b_cells@meta.data$original_barcode = paste0(b_cells@meta.data$iid,"_",b_cells@meta.data$original_barcode)

# find duplicates
non_duplicated_barcodes = b_cells@meta.data %>% dplyr::count(original_barcode) %>% filter(n==1)
combined_bcr = combined_bcr %>% filter(cell_id %in% non_duplicated_barcodes$original_barcode)

# filter to cells in VDJ
b_cells = subset(b_cells, original_barcode %in% combined_bcr$cell_id & original_barcode %in% non_duplicated_barcodes$original_barcode)

# filter VDJ to contigs with GEX data
combined_bcr = combined_bcr %>% filter(cell_id %in% b_cells@meta.data$original_barcode)

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
  productive_only = T,
  keep_highest_umi=TRUE)

# Print no. of cells passing QC
vdj_results_list[[2]]$obs %>% filter(
has_contig == "True" &
contig_QC_pass == "True" &
filter_contig == "FALSE") %>% nrow()

ann_df = vdj_results_list[[2]]$obs %>% data.frame
qc_stats = ann_df %>% dplyr::count(has_contig,contig_QC_pass,filter_contig,filter_contig_quality,filter_contig_VJ,filter_contig_VDJ)
write_csv(qc_stats,"vdj_qc_stats.csv")

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

# get back into seurat
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

rownames(b_cells@meta.data) = colnames(b_cells)
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

png(res=300,units="in",height=6,width=6,"diversity.png")
p2
dev.off()
