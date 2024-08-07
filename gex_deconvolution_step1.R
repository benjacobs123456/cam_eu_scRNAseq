# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
library(SoupX)
library(DoubletFinder)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
print(args)
index = as.numeric(args[1])
message("index is ",index)

# Load references
blueprint = celldex::BlueprintEncodeData()
monaco = celldex::MonacoImmuneData()
hpca = celldex::HumanPrimaryCellAtlasData()
dice = celldex::DatabaseImmuneCellExpressionData()

# set WD
setwd("/rds-d4/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined")

# link file
reference = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/EU_ID.txt")
reference$cohort = "EU"
reference$PBMC_PoolSize = as.character(reference$PBMC_PoolSize)
reference_cam = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/GT_ID.txt")
reference_cam$cohort = "Cam"
reference_cam$PBMC_PoolSize = as.character(reference_cam$PBMC_PoolSize)
reference = bind_rows(reference,reference_cam)

# define list of good batches
good_pbmc_batches = unique(reference$PBMC_GEX)
good_csf_batches = unique(reference$CSF_GEX)

# get sets of files
eu_pbmc_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_EU")
eu_pbmc_samples = eu_pbmc_samples[eu_pbmc_samples %in% good_pbmc_batches]
cam_pbmc_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_output/PBMC/")
cam_pbmc_samples = cam_pbmc_samples[cam_pbmc_samples %in% good_pbmc_batches]
eu_csf_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_EU")
eu_csf_samples = eu_csf_samples[eu_csf_samples %in% good_csf_batches]
cam_csf_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_output/CSF/")
cam_csf_samples = cam_csf_samples[cam_csf_samples %in% good_csf_batches]

unpooled_csf_samples = reference[reference$CSF_PoolSize==1,]
unpooled_pbmc_samples = reference[reference$PBMC_PoolSize==1,]

# get paramaters for this run
sample_list = c(eu_csf_samples, eu_pbmc_samples,
cam_csf_samples, cam_pbmc_samples)
message("getting batch name")
this_sample = sample_list[index]
message("Batch name is ",this_sample)
message("getting batch source")
batch_source = str_split(colnames(reference[grepl(this_sample,reference)])[1],pattern="_")[[1]][1]
message("getting cohort")
this_cohort = reference[reference[[paste0(batch_source,"_GEX")]]==this_sample,]$cohort[1]

message("running qc")

# define function to do the following:
# 1 read in sample for the batch
# 2 calculate basic QC metrics
# 3 get rid of bad cells
# 4 label with donor ids and case status
# 5 get rid of individuals with unpaired samples

# make empty df where we'll store QC stats
start_time=Sys.time()
do_gex_qc = function(x = this_sample,source = batch_source ,cohort = this_cohort){
message("Reading in GEX data for...")
message("Cohort:", cohort)
message("Batch:", x)
message("Source:", source)
message("Starting at ",start_time)

# run SoupX
raw = if(cohort == "Cam"){
  message("Looking in Cambridge folders")
  load10X(dataDir=paste0("/rds-d4/project/sjs1016/rds-sjs1016-msgen/10X_5prime/CRv5_output/",source,"/",x,"/outs/count"))
  } else if(cohort == "EU") {
  message("Looking in EU folders")
  load10X(dataDir=paste0("/rds-d4/project/sjs1016/rds-sjs1016-msgen/10X_5prime/CRv5_EU/",x,"/outs/count"))
  }

# manual method
nonExpressedGeneList = list(HB = c("HBB", "HBA2","HBA1"), IG = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"))
useToEst = estimateNonExpressingCells(raw, nonExpressedGeneList = nonExpressedGeneList)
sc = calculateContaminationFraction(raw,nonExpressedGeneList,useToEst=useToEst)
out = adjustCounts(sc)
df = out

print("Making Seurat object")
df = df %>% CreateSeuratObject()

# stash initial cell number
pre_qc_cells = df@meta.data %>% nrow()
rho = sc$metaData$rho[1]

# stash cohort name in metadata
df[['cohort']] = cohort

# define % MT genes and Hb
df[["percent.mt"]] = PercentageFeatureSet(df,pattern="^MT-")
df[["percent.hb"]] = PercentageFeatureSet(df,pattern="^HB-")
# nFeatures & MT
df = subset(df, subset = nFeature_RNA>100 & percent.mt < 10)

# stash cell number after basic qc
post_mt_nfeature_filters_cells = df@meta.data %>% nrow()

# initial seurat workflow
df = df %>% SCTransform() %>% RunPCA() %>% RunUMAP(dims=1:10) %>% FindNeighbors() %>% FindClusters()

# Run DoubletFinder
# Homotypic Doublet Proportion Estimate
homotypic.prop = modelHomotypic(df@meta.data$seurat_clusters)
nExp_poi = round(0.075*nrow(df@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))

df = doubletFinder_v3(df, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
colnames(df@meta.data)[12] = "dbfinder_call"

no_of_doublets_dbfinder = df@meta.data %>% filter(dbfinder_call=="Doublet") %>% nrow()

# subset
df = subset(df, subset = dbfinder_call =="Singlet")

# set batch id
df[["batch_id"]] = x

# set body fluids source
df[["source"]] = source

# add in donor IDs
# needs to be done separately for unpooled samples and multiplexed samples
if(x %in% unpooled_csf_samples$CSF_GEX){
  # stash in meta data
  this_donor = unpooled_csf_samples[unpooled_csf_samples$CSF_GEX==x,]
  df[["donor.id"]] =  this_donor$ID
} else if(x %in% unpooled_pbmc_samples$PBMC_GEX){
  # stash in meta data
  this_donor = unpooled_pbmc_samples[unpooled_pbmc_samples$PBMC_GEX==x,]
  df[["donor.id"]] =  this_donor$ID
} else {
  # read in decon for this batch
  decon = read_tsv(paste0("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/Donor_Decon/",x,"_donor_ids.tsv"))
  print("Merge with decon files to get per-sample files")

  # get donor IDs
  # remove decon cells where there is >1 call
  # or where the wrong sample was sent (in EU case)

  donor_ids = sapply(rownames(df@meta.data),function(x){
    this_cell = x
    this_donor = decon[decon$cell==this_cell,'donor_id']

    if(cohort == "Cam" & length(this_donor)==1){
      return(this_donor$donor_id)
    } else if(cohort == "EU" & length(this_donor)==1)
        col_of_interest = paste0(source,"_GEX")
        this_batch = reference[reference[[col_of_interest]]==x,]$Batch
        if(all(this_batch == 2)){
          return(this_donor$donor_id)
        } else if(all(this_batch ==1) & grepl("TUM",this_donor)){
          return(this_donor$donor_id)
    } else {
      return(NA)
    }
})

# stash in meta data
df[["donor.id"]] =  unlist(donor_ids)
}

# get some QC metrics pre-QC
no_unassigned_cells = df@meta.data %>% filter(donor.id=="unassigned") %>% nrow()
no_doublets = df@meta.data %>% filter(donor.id=="doublet") %>% nrow()

# now do QC
# remove doublets & unassigned cells
df = subset(df, subset = donor.id  != "doublet")
# unassigned
df = subset(df, subset = donor.id  != "unassigned")

final_cells = df@meta.data %>% nrow()

# store to QC df
qc_metrics = data.frame(x, cohort, source, pre_qc_cells,rho,  post_mt_nfeature_filters_cells, no_of_doublets_dbfinder, no_unassigned_cells, no_doublets, final_cells)

write_csv(qc_metrics,paste0("./basic_qc/qc_metrics_",x,".csv"))
# now normalise for the batch using SC transform
df = df %>% SCTransform(assay="RNA",vars.to.regress="percent.mt",method = "glmGamPoi",variable.features.n = 10000,return.only.var.genes = FALSE)

# plots
png(paste0("./basic_qc/individual_qc_plots/",x,"_qc_plot.png"),height=8,width=8,units="in",res=300)
print(VlnPlot(df,features=c("nCount_RNA","nFeature_RNA","percent.mt"),group.by="donor.id",fill="donor.id")+labs(x="Donor ID"))
dev.off()

message("Running SingleR")
blueprint_annotations = SingleR(test = df@assays$RNA@data, ref = blueprint, labels = blueprint$label.fine)
df[["ann_blueprint"]] = blueprint_annotations$pruned.labels
monaco_annotations = SingleR(test = df@assays$RNA@data, ref = monaco, labels = monaco$label.fine)
df[["ann_monaco"]] = monaco_annotations$pruned.labels
hpca_annotations = SingleR(test = df@assays$RNA@data, ref = hpca, labels = hpca$label.fine)
df[["ann_hpca"]] = hpca_annotations$pruned.labels
dice_annotations = SingleR(test = df@assays$RNA@data, ref = dice, labels = dice$label.fine)
df[["ann_dice"]] = dice_annotations$pruned.labels
message("Finished SingleR")

for(i in unique(df@meta.data$donor.id)){
  png(paste0("./basic_qc/individual_cluster_plots/",x,"individual_cluster_plots_",i,"_qc_plot.png"),height=6,width=6,units="in",res=300)
  plot_df = subset(df, subset = donor.id == i)
  print(DimPlot(plot_df,group.by = "ann_blueprint", label=T, repel=T)+ggtitle(paste0(source,"_",i)))
  dev.off()
}

now_time=Sys.time()
message("Started at: ",start_time,". Currently: ", now_time)

output_file_name = paste0("./basic_qc/individual_post_qc_datasets/post_qc_batch_",x,"source_",source,"cohort_",cohort,".rds")
saveRDS(df,output_file_name)
}

do_gex_qc()
