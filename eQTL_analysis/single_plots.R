library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
#library(Matrix.utils)
library(reshape2)
library(data.table)
library(plyr)
library(gridExtra)
library(ggrepel)
library(grid)


# directories
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
setwd(data.dir)

# Figures 
datall = readRDS(paste0(data.dir,"/","Complete_results_all_celltypes.rds"))
datcsf = datall[datall$source == "CSF",]
datpbmc = datall[datall$source == "PBMC",]

genlist=read.table("/xx/xxx/SC/CAM_TUM/eqtl/gene_list.txt",h=T)
fam_filename = "/xx/xxx/SC/CAM_TUM/eqtl/genotype_data/merged/merged_1to22_SampleQC_VariantQC2.psam"
meta_filename="/xx/xxx/SC/CAM_TUM/datasets/gex_metadata_new.csv"
genIDsTUM_filename="/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs.csv"
genIDsCAM_filename="/xx/xxx/SC/CAM_TUM/sc_for_tum/GT_ID.txt"
genotype_filename="/xx/xxx/SC/CAM_TUM/eqtl/genotype_data/merged/merged_1to22_SampleQC_VariantQC"

fam = fread(fam_filename)
fam = data.frame(fam)
fam$Cohort = "TUM"
fam[grepl("GT_20_",fam[,1]) == T,"Cohort"] = "CAM"

## Pheno data
meta = read.csv2(meta_filename)
tum_ids = read.csv(genIDsTUM_filename)
tum_ids$X = NULL
colnames(tum_ids)[1]="Sample"
cam_ids = read.table(genIDsCAM_filename,h=T,sep="\t")
gen_ids = unique(rbind(tum_ids[,c("Sample","GenID")],setnames(cam_ids[,c("ShortID","GenotypingID")],names(tum_ids[,c("Sample","GenID")]))))
gen_ids = gen_ids[is.na(gen_ids$GenID) == F & gen_ids$GenID != "",]
colnames(gen_ids)[1] = "donor.id"
meta = join(meta,gen_ids,by="donor.id")
meta = meta[is.na(meta$GenID) == F,]
meta = unique(meta[,c("PatID","phenotype","Age","Sex","GenID")])
colnames(meta) = c("Sample","Phenotype","Age","Sex","IID")
meta = meta[meta$IID %in% data.frame(fam)[,1] == T,]
metakeep = meta

# Read gene expression data
gex = read.table("/xx/xxx/SC/CAM_TUM/eqtl/CSF/data_gex.txt",h=T); gexcsf=gex
gex = read.table("/xx/xxx/SC/CAM_TUM/eqtl/PBMC/data_gex.txt",h=T); gexpbmc=gex


args = commandArgs(TRUE)
args = c("rs4810801","PREX1","B.cells","CD4.T.cells","CD8.T.cells","CSF","/xx/xxx/SC/CAM_TUM/eqtl/Targeted/")
#args = c("rs61909096","ETS1","CD4.T.cells","CD8.T.cells","B.cells","CSF","/xx/xxx/SC/CAM_TUM/eqtl/Targeted/")



snp=args[1]
gene=args[2]
ctout1=args[3]
ctout2=args[4]
ctout3=args[5]
cohort="All"
source=args[6]
round=1
data.dir=args[7]

celltypes = c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")



source(paste0(data.dir,"/Scripts/eqtl_plot_function.R"))

plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,cell_types=celltypes,window=500000)

save(plots,file=paste0(data.dir,"/figures/",gene,".RData"))

png(paste0(data.dir,"/figures/eQTL_",gene,"_",snp,"_",ctout1,"_",source,".png"),width=8,height=8,units="in",res=300)
grid.arrange(plots[[1]],plots[[4]],plots[[13]],nrow=2,layout_matrix=cbind(c(1,3),c(2,3)))
grid.text("A",x=unit(0.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(0.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
dev.off()

png(paste0(data.dir,"/figures/eQTL_",gene,"_",snp,"_",ctout1,"_",source,"_V2.png"),width=12,height=8,units="in",res=300)
grid.arrange(plots[[1]],plots[[4]],plots[[13]],nrow=2,layout_matrix=cbind(c(1,2),c(3,3)))
grid.text("A",x=unit(0.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(6.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(0.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
dev.off()

