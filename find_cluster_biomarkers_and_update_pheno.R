# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)

# set WD
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/datasets")

##########################################
#   Phenotype cleaning                   #
##########################################

# read in outputs from deconvolution & integration step
all_combo = readRDS("all_combo_with_UMAP_PCs_50resolution2.5.rds")

# get counts before merge
counts = all_combo@meta.data %>% dplyr::count(donor.id,source)
write_csv(counts,"very_raw_counts_per_donor_before_merge.csv")


# manually correct the ID of one donor from TUM
# add individual_id
all_combo@meta.data = all_combo@meta.data %>%
  mutate(donor.id = ifelse(
  batch_id == "TUM_B2P2" & donor.id == "donor1",
  "B2P2S2",
  donor.id))

# clean and integrate phenotype data
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

# get rid of donors without an ID (i.e failed deconvolution)
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
cam_pheno =  cam_pheno %>%
dplyr::select(shortID,Category,Age,Gender,Category_fine,Oligoclonal,Diagnosis) %>%
dplyr::rename("Sex"= Gender) %>%
mutate(Sex = case_when(
  Sex == "2" ~ "F",
  Sex == "1" ~ "M"
)) %>%
dplyr::rename("ms_subtype" = Category_fine)

# repeat for EU samples
eu_pheno = read_csv("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/references/EU_pheno.csv")
pheno = bind_rows(cam_pheno,eu_pheno)
pheno$donor.id = pheno$shortID
pheno = pheno %>% select(-shortID)

pheno = pheno %>%
mutate(Category = case_when(
Category == "Noninflammatory" ~ "NIND",
Category != "Noninflammatory" ~ Category
))

pheno = pheno %>% mutate(Oligoclonal = case_when(
Oligoclonal == "Pos(3to9)" ~ TRUE,
Oligoclonal == "Pos(>10)" ~ TRUE,
Oligoclonal == "Positive" ~ TRUE,
Oligoclonal == "Negative" ~ FALSE
))

pheno = pheno %>%
dplyr::rename("OCB" = Oligoclonal)

# add in TUM pheno data
tum_pheno = read_csv("/rds/user/hpcjaco1/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/phenotypes_all_new.csv")

tum_pheno = tum_pheno %>%
  select(Sample,PatID,Sex,Age,OCB,GSA_ID) %>%
  dplyr::rename("donor.id" = Sample)

tum_pheno2 = read_csv("/rds/user/hpcjaco1/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/group_assignment_TUM.csv")

tum_pheno = tum_pheno %>%
  left_join(
    tum_pheno2 %>%
      dplyr::select(PatID,CC_Status,MS_subtype) %>%
      distinct(PatID,.keep_all=T),
      by="PatID"
  ) %>%
  dplyr::rename("Category" = CC_Status) %>%
  mutate(Category = case_when(
    Category == "ONIND" ~ "NIND",
    Category == "Infectious" ~ "OINDI",
    Category == "MS" ~ "MS",
    Category == "Inflammatory" ~ "OIND",
    Category == "OIND" ~ "OIND"
    )) %>%
  mutate(
    ms_subtype = case_when(
      MS_subtype == "PPMS" ~ "PPMS",
      MS_subtype == "RRMS" ~ "RMS",
      MS_subtype == "SPMS" ~ "SPMS"
    )
  )

tum_diagnosis = read_csv("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/references/TUM_clinical_info_04_04_24.csv")

tum_pheno = tum_pheno %>%
  left_join(tum_diagnosis %>%
  dplyr::select(PatID,Diagnosis),
  by="PatID")


# combine
pheno = bind_rows(pheno,tum_pheno)

# clean
pheno = pheno %>%
  mutate(ms_subtype = ifelse(Category != "MS",NA,ms_subtype))

# get pheno information from metadata
pheno_data_overall = data.frame()
donor_ids = all_combo@meta.data$donor.id %>% unique()
for(i in c(1:length(donor_ids))){
  message(i)
  this_donor = donor_ids[i]

  pheno_data = NULL

  if(this_donor %in% pheno$donor.id){
    message("Match")
    pheno_data <<- pheno %>% filter(donor.id==this_donor)  %>%
    mutate(new_donor_id = donor.id) %>% distinct(donor.id,.keep_all=T)
  } else {
    message("Not matched, getting GSA ID")
    this_donor_new_id = str_split_fixed(this_donor,"_",2)[2]
    pheno_data <<- pheno %>% filter(GSA_ID==this_donor_new_id) %>%
    mutate(new_donor_id = GSA_ID) %>% distinct(PatID,.keep_all=T)
  }

  pheno_data$original_donor_id = this_donor
  pheno_data_overall <<- bind_rows(pheno_data_overall,pheno_data)
}

# merge
all_combo@meta.data = all_combo@meta.data %>%
  left_join(
    pheno_data_overall %>%
    select(-donor.id) %>%
    dplyr::rename("donor.id" = original_donor_id),
    by = "donor.id")

# create unique individual ID
all_combo@meta.data = all_combo@meta.data %>%
  mutate(iid = ifelse(processing_site == "CAM",donor.id,PatID))

# get counts before merge
counts = all_combo@meta.data %>% dplyr::count(iid,source)
write_csv(counts,"raw_counts_per_donor_before_merge.csv")

#######################################
# Merge overlapping patients
#######################################

sample_overlap = read_csv("/rds/user/hpcjaco1/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/Sample_Overlap.csv")

rownames(all_combo@meta.data) = colnames(all_combo)

# write demographics file

pheno_simple = all_combo@meta.data %>%
  distinct(iid,source,.keep_all=T) %>%
  dplyr::select(iid,Age,Sex,Category,Diagnosis,ms_subtype,OCB,processing_site,source)
rownames(pheno_simple) = NULL

pheno_simple = pheno_simple %>%
  tidyr::pivot_wider(values_from = source,
  names_from = source,
  id_cols = 1:8) %>%
  mutate(CSF = ifelse(!is.na(CSF),"Yes","No")) %>%
  mutate(PBMC = ifelse(!is.na(PBMC),"Yes","No")) %>%
  arrange(Category,processing_site) %>%
  mutate(replicate = ifelse(iid %in% sample_overlap$CAM_SAMPLE_ID | iid %in% sample_overlap$TUM_PatID,"*",""))
write_csv(pheno_simple,"demographics.csv")

all_combo@meta.data = all_combo@meta.data %>%
  left_join(
    sample_overlap %>%
      dplyr::rename("iid" = TUM_PatID),
    by="iid"
  ) %>%
  mutate(iid_prior_to_merging_duplicates = iid) %>%
  mutate(iid = ifelse(iid %in% sample_overlap$TUM_PatID,CAM_SAMPLE_ID,iid))
rownames(all_combo@meta.data) = colnames(all_combo)

message("Donor n after duplicates merged")
tbl1 = all_combo@meta.data %>%
  distinct(iid,source,.keep_all=T) %>%
  dplyr::count(Category,source,processing_site)
tbl1
write_csv(tbl1,"donor_counts_no_duplicates.csv")

rownames(all_combo@meta.data) = colnames(all_combo)

##########################################
#   Cluster biomarkers                   #
##########################################

# remove singletons
all_combo = subset(all_combo, subset = seurat_clusters != "singleton")

# get high count clusters
high_count_clusters = all_combo@meta.data %>% dplyr::count(seurat_clusters) %>% filter(n>=10)

# remove clusters with n<10 cells
all_combo = subset(all_combo, subset = seurat_clusters %in% high_count_clusters$seurat_clusters)

# save
saveRDS(all_combo,"all_combo_with_pheno.rds")

# make dim plots
p0 = DimPlot(all_combo,label=T, repel=T)
png("basic_dimplot_labelled.png",res=300,units="in",width=8,height=8)
p0
dev.off()

p1 = DimPlot(all_combo,split.by="source", repel=T)
png("basic_dimplot_by_source.png",res=300,units="in",width=8,height=8)
p1
dev.off()

cluster_markers = c("CLEC10A","FCER1A","S100A8", "S100A9", "CD14","CD1C","FCGR3A", "SMIM25","CD3E","CD3D","CD8A", "CD8B","ANXA1","CCL5","CCR7","MKI67","TCF7","IL7R","TRDC","TRGC1","KLRB1", "GNLY", "NKG7","PRSS57","CD79A", "MS4A1","sct_IGKC","LILRA4", "IRF8", "PF4","FOXP3","CTLA4","CD38","CD27","IGHM","IGHG1","CD24")
p4 = FeaturePlot(all_combo,features=cluster_markers, repel=T)
png("canonical_markers.png",res=300,units="in",width=16,height=16)
p4
dev.off()

p5 = DotPlot(all_combo,features=cluster_markers) + theme(axis.text.x = element_text(angle=90))
png("dotplot_canonical_markers.png",res=300,units="in",width=12,height=16)
p5
dev.off()

short_cluster_markers = c("CLEC10A","S100A8", "CD16","CD3E","CD8A","IL7R", "NKG7","CD79A", "PF4","CD38")
p5 = FeaturePlot(all_combo,features=short_cluster_markers,raster=FALSE)
png("few_canonical_markers.png",res=300,units="in",width=12,height=8)
p5
dev.off()

# look at singleR annotations
totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_blueprint) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_blueprint) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_blueprint))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"singler_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_blueprint,label=ann_blueprint))+geom_col()+geom_text(angle=90)+theme_bw()
png("singler_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()

totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_monaco) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_monaco) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_monaco))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"monaco_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_monaco,label=ann_monaco))+geom_col()+geom_text(angle=90)+theme_bw()
png("monaco_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()

totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_dice) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_dice) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_dice))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"dice_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_dice,label=ann_dice))+geom_col()+geom_text(angle=90)+theme_bw()
png("dice_calls.png",res=300,units="in",width=12,height=8)
p
dev.off()

####################################
# cluster biomarkers
####################################

# Find markers
all_combo_markers = FindAllMarkers(all_combo, min.pct=0.5, logfc.threshold = 1,only.pos=TRUE,recorrect_umi=FALSE)
write_csv(all_combo_markers,"cluster_biomarkers.csv")

# summarise biomarkers
biomarker_summary = all_combo_markers
biomarkers = lapply(unique(biomarker_summary$cluster),function(x){
  this_cluster = biomarker_summary %>% filter(cluster==x) %>% arrange(desc(avg_log2FC)) %>% head(10)
  genes = paste0(this_cluster$gene,collapse=", ")
  data.frame(x,genes)
})
biomarkers = do.call("bind_rows",biomarkers) %>% tibble::as_tibble()
write_csv(biomarkers,"biomarker_summary.csv")
