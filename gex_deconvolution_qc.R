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

# define unmatched samples
unmatched_samples = reference %>% filter(CSF_GEX =="X" | PBMC_GEX =="X")
bad_batches = c("S18303-F11","S17624-E11")

# restrict to samples with both CSF & PBMC
reference = reference %>% filter(CSF_GEX !="X" & PBMC_GEX!="X")

# remove bad batches
reference = reference %>% filter(!(CSF_GEX %in% bad_batches | PBMC_GEX %in% bad_batches))

# define list of good batches
good_pbmc_batches = unique(reference$PBMC_GEX)
good_csf_batches = unique(reference$CSF_GEX)

# get sets of files
eu_pbmc_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_EU")
eu_pbmc_samples = eu_pbmc_samples[!(eu_pbmc_samples %in% bad_batches)]
eu_pbmc_samples = eu_pbmc_samples[!(eu_pbmc_samples %in% unmatched_samples$PBMC_GEX)]
eu_pbmc_samples = eu_pbmc_samples[eu_pbmc_samples %in% good_pbmc_batches]
cam_pbmc_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_output/PBMC/")
cam_pbmc_samples = cam_pbmc_samples[!(cam_pbmc_samples %in% bad_batches)]
cam_pbmc_samples = cam_pbmc_samples[!(cam_pbmc_samples %in% unmatched_samples$PBMC_GEX)]
cam_pbmc_samples = cam_pbmc_samples[cam_pbmc_samples %in% good_pbmc_batches]
eu_csf_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_EU")
eu_csf_samples = eu_csf_samples[!(eu_csf_samples %in% bad_batches)]
eu_csf_samples = eu_csf_samples[!(eu_csf_samples %in% unmatched_samples$CSF_GEX)]
eu_csf_samples = eu_csf_samples[eu_csf_samples %in% good_csf_batches]
cam_csf_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_output/CSF/")
cam_csf_samples = cam_csf_samples[!(cam_csf_samples %in% bad_batches)]
cam_csf_samples = cam_csf_samples[!(cam_csf_samples %in% unmatched_samples$CSF_GEX)]
cam_csf_samples = cam_csf_samples[cam_csf_samples %in% good_csf_batches]

unpooled_csf_samples = reference[reference$CSF_PoolSize==1,]
unpooled_pbmc_samples = reference[reference$PBMC_PoolSize==1,]

cluster_plots = list()
vln_plots = list()

# define function to do the following:
# 1 read in sample for the batch
# 2 calculate basic QC metrics
# 3 get rid of bad cells
# 4 label with donor ids and case status
# 5 get rid of individuals with unpaired samples

# make empty df where we'll store QC stats
qc_df = data.frame()
start_time=Sys.time()
do_gex_qc = function(x,source,cohort){
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
no_cells = df@meta.data %>% nrow()

# stash cohort name in metadata
df[['cohort']] = cohort

# define % MT genes and Hb
df[["percent.mt"]] = PercentageFeatureSet(df,pattern="^MT-")
df[["percent.hb"]] = PercentageFeatureSet(df,pattern="^HB-")
# nFeatures & MT
df = subset(df, subset = nFeature_RNA>100 & percent.mt < 10)

# stash cell number after basic qc
no_cells_filters = df@meta.data %>% nrow()

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
no_donors = df@meta.data %>% distinct(donor.id) %>% nrow()
no_unassigned_cells = df@meta.data %>% filter(donor.id=="unassigned") %>% nrow()
no_doublets = df@meta.data %>% filter(donor.id=="doublet") %>% nrow()
prop_unassigned_cells = no_unassigned_cells / no_cells *100
prop_doublets = no_doublets / no_cells *100
rho = sc$metaData$rho[1]

# make QC plots
violin_qc = VlnPlot(df,features=c("nCount_RNA","nFeature_RNA","percent.mt"),group.by="donor.id",fill="donor.id")+labs(x="Donor ID")+ggtitle(paste0("Batch: ",x,". Cohort: ",cohort,". Source: ",source))
vln_plots[[length(vln_plots)+1]] = violin_qc


# now do QC
# remove doublets & unassigned cells
df = subset(df, subset = donor.id  != "doublet")
# unassigned
df = subset(df, subset = donor.id  != "unassigned")

# store to QC df
qc_metrics = data.frame(source,x,rho,no_of_doublets_dbfinder,no_cells,no_cells_filters,no_donors,no_unassigned_cells,prop_unassigned_cells,prop_doublets,post_QC_n_cells = nrow(df@meta.data))
qc_df <<- bind_rows(qc_df,qc_metrics)

# now normalise for the batch using SC transform
df = df %>% SCTransform(assay="RNA",vars.to.regress="percent.mt",method = "glmGamPoi")

# plots
png(paste0("./basic_qc/individual_cluster_plots/",x,"individual_cellranger_cluster_plots_qc_plot.png"),height=8,width=8,units="in",res=300)
print(DimPlot(df,label=T,split.by="donor.id")+ggtitle(paste0("Batch: ",x,". Source: ",source)))
dev.off()

now_time=Sys.time()
message("Started at: ",start_time,". Currently: ", now_time)
return(df)
}

# cambridge samples
cam_pbmc_output = sapply(cam_pbmc_samples, do_gex_qc,source="PBMC",cohort="Cam")
cam_csf_output = sapply(cam_csf_samples, do_gex_qc,source="CSF",cohort="Cam")

# EU samples
eu_pbmc_output = sapply(eu_pbmc_samples, do_gex_qc,source="PBMC",cohort="EU")
eu_csf_output = sapply(eu_csf_samples, do_gex_qc,source="CSF",cohort="EU")

# merge cam
cam_csf_combo = merge(x=cam_csf_output[[1]],y=cam_csf_output[2:length(cam_csf_output)])
cam_pbmc_combo = merge(x=cam_pbmc_output[[1]],y=cam_pbmc_output[2:length(cam_pbmc_output)])
cam_all_combo = merge(x=cam_csf_combo,y=cam_pbmc_combo)

# merge eu
eu_csf_combo = merge(x=eu_csf_output[[1]],y=eu_csf_output[2:length(eu_csf_output)])
eu_pbmc_combo = merge(x=eu_pbmc_output[[1]],y=eu_pbmc_output[2:length(eu_pbmc_output)])
eu_all_combo = merge(x=eu_csf_combo,y=eu_pbmc_combo)

# merge all samples together
all_combo = merge(x=eu_all_combo,y=cam_all_combo)

# Integrate datasets
var.features = SelectIntegrationFeatures(append(cam_pbmc_output,cam_csf_output,eu_pbmc_output,eu_csf_output), nfeatures = 3000)
VariableFeatures(all_combo) = var.features

# Run PCA
all_combo = all_combo %>% RunPCA(assay.use="SCT")

# Run Harmony
all_combo = all_combo %>% RunHarmony(assay.use="SCT",group.by.vars=c("batch_id","cohort"))

# save
saveRDS(all_combo, file = "./datasets/gex_pbmc_csf_combined.rds")
write_csv(qc_df,"./basic_qc/qc_metrics.csv")

# plots
png(paste0("./basic_qc/violin_qc_plot.png"),height=32,width=32,units="in",res=300)
print(do.call("grid.arrange",c(vln_plots,ncol=6)))
dev.off()

# plots
png(paste0("./basic_qc/individual_cellranger_cluster_plots_qc_plot.png"),height=32,width=32,units="in",res=300)
print(DimPlot(all_combo,label=T,split.by="donor.id",fill="source")+ggtitle(paste0("Donor ID: ",donor.id,". Source: ",source)))
dev.off()
