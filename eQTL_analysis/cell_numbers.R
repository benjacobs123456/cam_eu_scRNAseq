# Import libraries 
library(tidyverse) 
library(broom) 
library(ggplot2) 
#library(qvalue) 
library(matrixStats) 
library(data.table) 
library(dplyr) 
library(magrittr) 
library(tidyr)
library(plyr)
library(gdata)
library(RNOmni)

csf_numsfile ="/xx/xxx/SC/CAM_TUM/eqtl/CSF/numbers_cells.txt"
pbmc_numsfile = "/xx/xxx/SC/CAM_TUM/eqtl/PBMC/numbers_cells.txt"

genIDsTUM_filename="/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs.csv"
genIDsCAM_filename="/xx/xxx/SC/CAM_TUM/sc_for_tum/GT_ID.txt"
meta_filename="/xx/xxx/SC/CAM_TUM/datasets/gex_metadata_new.csv"
fam_filename = "/xx/xxx/SC/CAM_TUM/eqtl/genotype_data/merged/merged_1to22_SampleQC_VariantQC2.psam"

tum_ids = read.csv(genIDsTUM_filename)
tum_ids$X = NULL
colnames(tum_ids)[1]="Sample"
cam_ids = read.table(genIDsCAM_filename,h=T,sep="\t")
gen_ids = unique(rbind(tum_ids[,c("Sample","GenID")],setnames(cam_ids[,c("ShortID","GenotypingID")],names(tum_ids[,c("Sample","GenID")]))))
gen_ids = gen_ids[is.na(gen_ids$GenID) == F & gen_ids$GenID != "",]
colnames(gen_ids)[1] = "donor.id"

meta = fread(meta_filename,select=c("PatID","donor.id"),sep=";")
meta = unique(meta)
meta = join(meta,gen_ids,by="donor.id")
meta = meta[is.na(meta$GenID) == F,]

fam = fread(fam_filename)
fam = data.frame(fam)
meta = meta[meta$GenID %in% fam[,1] == T,]

pbmc_nums = read.table(pbmc_numsfile,h=T,sep="\t"); pbmc_nums$source = "PBMC"
csf_nums = read.table(csf_numsfile,h=T,sep="\t"); csf_nums$source = "CSF"
nums = rbind(pbmc_nums,csf_nums)
nums = separate(nums,col="Var1",sep="_",into=c("celltype","cohort","rem"),remove=F)
nums$CTC = paste0(nums$celltype,"_",nums$cohort)
nums$PatID = NA
for (i in 1:nrow(nums)){
	nums[i,"PatID"] = gsub(paste0(nums[i,"CTC"],"_"),"",nums[i,"Var1"])
}

nums = nums[nums$PatID %in% meta$PatID == T,]

nums[nums$cohort != "MS","cohort"] = "Control"

nums =  nums %>% group_by(source,cohort,celltype) %>% dplyr::summarise(n=sum(Freq),.groups='drop')
nums = data.frame(nums)

write.table(nums,"/xx/xxx/SC/CAM_TUM/eqtl/numbers_cells.txt",sep="\t",row.names=F,quote=F)