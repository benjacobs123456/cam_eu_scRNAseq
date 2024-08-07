library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Matrix.utils)
library(reshape2)
library(data.table)
library(plyr)
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(tidyr)

args = commandArgs(TRUE)

# directories
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted/"
setwd(data.dir)


res = fread(paste0(data.dir,"All_fdr_eSNPs.txt")) # read file with eSNPs
res$source_cohort_ct = paste0(res$source,"_",res$cohort,"_",res$cell_type)
res = res[,c("X.CHROM","POS","SNP","GENE","P","FDR","source_cohort_ct")]
colnames(res)[1] = "CHR"
res = separate(res,col="source_cohort_ct",sep="_",into=c("Source","Cohort","cell_type"),rem=F)
resk = res

results=res

cell_types = unique(res$cell_type)

gnames = read.table(paste0(data.dir,"../references/genes_names.txt"))
colnames(gnames) = c("ENSEMBLID","GENE")
results = join(results,gnames,by="GENE",type="left")
results$Comb = paste0(results$CHR,"_",results$POS,"_",results$ENSEMBLID)
results$CHR_POS = paste0(results$CHR,"_",results$POS)
gtex_tissues = c("Whole_Blood","Adipose_Subcutaneous","Esophagus_Mucosa","Adipose_Visceral_Omentum","Esophagus_Muscularis","Adrenal_Gland","Heart_Atrial_Appendage","Artery_Aorta","Heart_Left_Ventricle","Artery_Coronary","Kidney_Cortex","Artery_Tibial","Liver",
                 "Brain_Amygdala","Lung","Brain_Anterior_cingulate_cortex_BA24","Minor_Salivary_Gland","Brain_Caudate_basal_ganglia","Muscle_Skeletal","Brain_Cerebellar_Hemisphere","Nerve_Tibial","Brain_Cerebellum","Ovary","Brain_Cortex","Pancreas","Brain_Frontal_Cortex_BA9","Pituitary",
                 "Brain_Hippocampus","Prostate","Brain_Hypothalamus","Skin_Not_Sun_Exposed_Suprapubic","Brain_Nucleus_accumbens_basal_ganglia","Skin_Sun_Exposed_Lower_leg","Brain_Putamen_basal_ganglia","Small_Intestine_Terminal_Ileum","Brain_Spinal_cord_cervical_c-1","Spleen",
                 "Brain_Substantia_nigra","Stomach","Breast_Mammary_Tissue","Testis","Cells_Cultured_fibroblasts","Thyroid","Cells_EBV-transformed_lymphocytes","Uterus","Colon_Sigmoid","Vagina","Colon_Transverse","Esophagus_Gastroesophageal_Junction")

results$siginGtex = 0
results$SNP_in_Gtex = 0
results$Gene_in_Gtex = 0
results$window_Gtex = 0
results$window_Gtex_WholeBlood = 0
for (tissue in gtex_tissues){
  print(tissue)
  reftex = fread(paste0(data.dir,"../references/GTEx_Analysis_v8_eQTL/",tissue,".v8.signif_variant_gene_pairs_sub.txt"),h=T)
  colnames(reftex) = c("comb1","ENSEMBLID1")
  reftex = separate(reftex,col="comb1",sep="_",into=c("CHR","POS","A1","A2","build"),rem=T)
  reftex = separate(reftex,col="ENSEMBLID1",sep="\\.",into=c("ENSEMBLID","rem"),rem=T)
  reftex = reftex[,c("CHR","POS","ENSEMBLID")]
  reftex$CHR = gsub("chr","",reftex$CHR)
  reftex$Comb = paste0(reftex$CHR,"_",reftex$POS,"_",reftex$ENSEMBLID)
  reftex$CHR_POS = paste0(reftex$CHR,"_",reftex$POS)
  reftexgenes = fread(paste0(data.dir,"../references/GTEx_Analysis_v8_eQTL/",tissue,".v8.egenes.txt"),h=T)
  reftexgenes = separate(reftexgenes,col="gene_id",sep="\\.",into=c("ENSEMBLID","rem"),rem=T)


  results[results$Comb %in% reftex$Comb == T,"siginGtex"] = 1
  print(table(results$siginGtex))
  results[results$CHR_POS %in% reftex$CHR_POS == T,"SNP_in_Gtex"] = 1
  results[results$ENSEMBLID %in% reftex$ENSEMBLID == T | results$ENSEMBLID %in% reftexgenes$ENSEMBLID == T,"Gene_in_Gtex"] = 1
  results[,paste0("Gtex_",tissue)] = 0
  results[results$Comb %in% reftex$Comb == T,paste0("Gtex_",tissue)] = 1
  
  for (i in 1:nrow(results)){
    #print(i)
    if (results[i,"siginGtex"] == 1 & is.na(results[i,"siginGtex"]) == F){
      results[i,"window_Gtex"] = 1
      if(tissue == "Whole_Blood"){results[i,"window_Gtex_WholeBlood"] = 1}
      next;
    }else{
      chr = results[i,"CHR"]
      start = results[i,"POS"] - 200000
      end = results[i,"POS"] + 200000
      gene = results[i,"ENSEMBLID"]
      reftex = data.frame(reftex)
      if (nrow(reftex[reftex$CHR == chr & reftex$POS >= start & reftex$POS <= end & reftex$ENSEMBLID == gene,]) > 0){
        results[i,"window_Gtex"] = 1
        if(tissue == "Whole_Blood"){results[i,"window_Gtex_WholeBlood"] = 1}
      }
    }
  }
  
  
}


#eQTLGen
results$Comb = paste0(results$SNP,"_",results$GENE)
ref = fread("../references/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")
ref$Comb = paste0(ref$SNP,"_",ref$GeneSymbol)
results$SNP_in_eQTLGen = 0; results$Gene_in_eQTLGen = 0; results$sigeQTLGen = NA
results$window_eQTLGen = 0


refall = fread("../references/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",select=c("SNP","GeneSymbol","Pvalue","SNPChr","SNPPos"))
results[results$SNP %in% refall$SNP == T,"SNP_in_eQTLGen"] = 1
results[results$GENE %in% refall$GeneSymbol == T,"Gene_in_eQTLGen"] = 1
results[results$SNP_in_eQTLGen == 1 & results$Gene_in_eQTLGen == 1,"sigeQTLGen"] = 0
results[results$Comb %in% ref$Comb == T,"sigeQTLGen"] = 1
refnom = refall[refall$Pvalue < 0.05,]
refnom$Comb = paste0(refnom$SNP,"_",refnom$GeneSymbol)
results[results$SNP_in_eQTLGen == 1 & results$Gene_in_eQTLGen == 1,"nomeQTLGen"] = 0
results[results$Comb %in% refnom$Comb == T,"nomeQTLGen"] = 1
for (i in 1:nrow(results)){
  print(i)
  pos = as.numeric(as.character(results[i,"POS"]))
  start = pos -100000
  end = pos +100000
  chr=as.numeric(as.character(results[i,"CHR"]))
  gene=results[i,"GENE"]
  #if(nrow(refall[refall$SNPChr == chr & refall$SNPPos >= start & refall$SNPPos <= end,]) > 0){ results[i,"window_eQTLGen"] = 0}
  if(nrow(ref[ref$SNPChr == chr & ref$SNPPos >= start & ref$SNPPos <= end & ref$GeneSymbol == as.character(gene),]) > 0){ results[i,"window_eQTLGen"] = 1}
}


#Yazar et al. 
yazar = read.csv2("../references/yazar.csv")
results$Comb = paste0(results$SNP,"_",results$ENSEMBLID)
yazar$Comb = paste0(yazar$SNP,"_",yazar$Gene.Ensembl.ID)
results$sigYazar = 0; results$window_Yazar = 0
results[results$Comb %in% yazar$Comb == T,"sigYazar"] = 1
for (i in 1:nrow(results)){
  print(i)
  pos = as.numeric(as.character(results[i,"POS"]))
  start = pos -100000
  end = pos+100000
  chr=as.numeric(as.character(results[i,"CHR"]))
  gene=results[i,"ENSEMBLID"]
  if(nrow(yazar[yazar$Chromosome == chr & yazar$Position >= start & yazar$Position <= end & yazar$Gene.Ensembl.ID == as.character(gene),]) > 0){ results[i,"window_Yazar"] = 1}
}

# Write results to file
fwrite(results,paste0(data.dir,"comparison_with_references.txt"))

