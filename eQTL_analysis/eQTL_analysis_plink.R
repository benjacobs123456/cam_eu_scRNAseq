

# Call the cell type 
args = commandArgs(trailingOnly=TRUE) 
cellLabel  <- args[1]
print(cellLabel) 
chrNumber <- args[2] 
print(chrNumber) 
cohort <- args[3]
source <- args[4]

 
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


# Output directory
output.dir = paste0("/xx/xxx/SC/CAM_TUM/eqtl/Targeted/",source,"/",cohort)
 
# Input filenames 
expression_filename = paste0("/xx/xxx/SC/CAM_TUM/eqtl/",source,"/data_gex.txt")
genotype_filename="/xx/xxx/SC/CAM_TUM/eqtl/Targeted/genotype_data/merged/merged_1to22_SampleQC_VariantQC"
sample_filename="/xx/xxx/SC/CAM_TUM/eqtl/Targeted/genotype_data/merged/merged_1to22_SampleQC_VariantQC.psam"
geneLoc_filename <- "/xx/xxx/SC/CAM_TUM/eqtl/gene_list.txt"
meta_filename="/xx/xxx/SC/CAM_TUM/eqtl/all_combo_with_updated_pheno_metadata.csv"
genIDsTUM_filename="/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs.csv"
genIDsTUM_2023_filename="/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs_2023.csv"
genIDsCAM_filename="/xx/xxx/SC/CAM_TUM/sc_for_tum/GT_ID.txt"
#snpLoc_filename=paste0("/xx/xxx/SC/eqtl/CSF/gen_data/chr",chrNumber,"_qc.snp-stats")
covariate_filename <-  "/xx/xxx/SC/CAM_TUM/eqtl/Targeted/covariates.txt"
pcs_filename <-  paste0("/xx/xxx/SC/CAM_TUM/eqtl/",source,"/PCs.csv")
#pheno_filename <- "/xx/xxx/SC/Seurat/phenotypes_all.csv"
fam_filename = "/xx/xxx/SC/CAM_TUM/eqtl/Targeted/genotype_data/merged/merged_1to22_SampleQC_VariantQC2.psam"
genelist_filename = paste0("/xx/xxx/SC/CAM_TUM/eqtl/",source,"/genes_",gsub("\\."," ",cellLabel),".txt")

# Read in files 
fam = fread(fam_filename)
fam = data.frame(fam)
fam$Cohort = "TUM"
fam[grepl("GT_20_",fam[,1]) == T,"Cohort"] = "CAM"


## Pheno data
meta = read.csv(meta_filename)
meta[meta$iid == "8A2H3PL1","iid"] = "PatID_5"
meta[meta$iid == "E9LEH7P8","iid"] = "PatID_56"
meta[meta$iid == "JEGK54J2","iid"] = "PatID_59"
meta[meta$iid == "RL4X6288","iid"] = "PatID_3"
meta[meta$iid == "WAF0FQN5","iid"] = "PatID_52"
meta[meta$iid == "WQ0EJVR1","iid"] = "PatID_60"
meta = unique(meta[,c("iid","phenotype","Age","Sex")])
tum_ids = read.csv(genIDsTUM_filename)
tum_ids$X = NULL
tum_ids = separate(tum_ids,col="GenID",sep="_TUM_TUM",into=c("IID","rem"),rem=F)
tum_ids$GenID = tum_ids$IID
colnames(tum_ids)[2]="Sample"
tum_ids_2023 = read.csv(genIDsTUM_2023_filename)
tum_ids_2023$X = NULL
colnames(tum_ids_2023)[2]="Sample"
cam_ids = read.table(genIDsCAM_filename,h=T,sep="\t")
gen_ids = unique(rbind(tum_ids[,c("Sample","GenID")],setnames(tum_ids_2023[,c("Sample","GenID_2023")],names(tum_ids[,c("Sample","GenID")])),setnames(cam_ids[,c("ShortID","GenotypingID")],names(tum_ids[,c("Sample","GenID")]))))
gen_ids = gen_ids[is.na(gen_ids$GenID) == F & gen_ids$GenID != "",]
colnames(gen_ids)[1] = "iid"
meta = join(meta,gen_ids,by="iid")
meta = unique(meta[,c("iid","phenotype","Age","Sex","GenID")])
meta = meta[is.na(meta$GenID) == F,]
colnames(meta) = c("Sample","Phenotype","Age","Sex","IID")
meta = meta[meta$IID %in% data.frame(fam)[,1] == T,]
if(cohort == "MS"){
  meta = meta[meta$Phenotype == "MS",]
} else if (cohort == "Control"){
  meta = meta[meta$Phenotype != "MS",]
} else if (cohort == "All"){
  meta = meta
}
pheno = meta

## Covariates 
covariate_df <- fread(covariate_filename) 
dim(covariate_df) 
print(covariate_df[1:5,1:5]) 
## PC factors
pcs_df <- read.csv(pcs_filename)
dim(pcs_df)
print(pcs_df[1:5,1:5]) 
pcs_df[,"Sample"] = pcs_df$X; pcs_df$X = NULL
pcs_df[pcs_df$Sample == "8A2H3PL1","Sample"] = "PatID_5"
pcs_df[pcs_df$Sample == "E9LEH7P8","Sample"] = "PatID_56"
pcs_df[pcs_df$Sample == "JEGK54J2","Sample"] = "PatID_59"
pcs_df[pcs_df$Sample == "RL4X6288","Sample"] = "PatID_3"
pcs_df[pcs_df$Sample == "WAF0FQN5","Sample"] = "PatID_52"
pcs_df[pcs_df$Sample == "WQ0EJVR1","Sample"] = "PatID_60"
pcs_df = join(pcs_df,pheno[,c("Sample","IID")],by="Sample",type="left")
colsel = colnames(pcs_df)[1:20]
covariate_df = join(covariate_df,pcs_df[,c(colsel,"IID")],by="IID",type="left")
covariate_df = covariate_df[covariate_df$IID %in% pcs_df$IID == T,]
colnames(fam)[1] = "IID"
covariate_df = join(covariate_df,fam[,c("IID","Cohort")],by="IID",type="left")

## Count matrix 
expression_df <- fread(expression_filename) 
relgenes = read.table(genelist_filename)
expression_df = expression_df %>% filter(V1 %in% relgenes$V1)
expression_df = data.frame(expression_df); genes = expression_df$V1
expression_df = expression_df[,startsWith(colnames(expression_df),paste0(cellLabel,"_"))]
rownames(expression_df)= genes
dim(expression_df) 
print(expression_df[1:5,1:5])
expression_df=t(expression_df)
expression_df = data.frame(expression_df)
expression_df$ID = rownames(expression_df)
if(cohort == "MS"){
  expression_df = expression_df[grepl(paste0(cohort,"_"),expression_df$ID) == T,]
}
if(cohort == "Control"){
  expression_df = expression_df[grepl(paste0("MS","_"),expression_df$ID) == F,]
}

expression_df = separate(expression_df,col="ID",sep="_",into=c("CT", "cohort","rest","rest2","rest3"),rem=F)
for (i in 1:nrow(expression_df)){expression_df[i,"Sample"] = gsub(paste0(expression_df[i,"CT"],"_",expression_df[i,"cohort"],"_"),"",expression_df[i,"ID"]);}
expression_df$CT = NULL; expression_df$cohort = NULL; expression_df$rest = NULL; expression_df$rest2 = NULL; expression_df$rest3 = NULL 
expression_df = join(expression_df,pheno[,c("Sample","IID")],by="Sample",type="left")
expression_df$Sample = NULL
expression_df = expression_df[expression_df$IID %in% covariate_df$IID == T,]
expression_df = expression_df[,c(ncol(expression_df),1:(ncol(expression_df)-1))]

expression_df = expression_df[expression_df$IID %in% pheno$IID == T & expression_df$IID %in% covariate_df$IID == T,]
covariate_df = covariate_df[covariate_df$IID %in% pheno$IID == T & covariate_df$IID %in% expression_df$IID == T,]

## Gene location file 
geneLoc_df <- fread(geneLoc_filename) 
dim(geneLoc_df) 
print(geneLoc_df[1:5,1:5]) 
geneLoc_df = geneLoc_df[geneLoc_df$chromosome_name == chrNumber,]
expression_df = expression_df[,colnames(expression_df) %in% c("IID",geneLoc_df$hgnc_symbol) == T]

## Remove genes with 0 expression in all samples 
gene_ids <- colnames(expression_df[-1]) 
sample_ids <- expression_df$IID  
 
# Work out which ones have colSums != 0: 
for (i in 2:ncol(expression_df)){
expression_df[,i] = as.numeric(as.character(expression_df[,i]))}
i <- (colSums(expression_df[,-1], na.rm=T) != 0) 
nonzerogenes <- names(i[i==TRUE]) 
expression_df_nonzero <- expression_df[,nonzerogenes] %>%   # !!!! removed ,with=FALSE
  mutate (sampleid = sample_ids) %>%    
  select(sampleid,everything()) 
print(expression_df_nonzero[1:5,1:5]) 
dim(expression_df_nonzero)  
 
# Number of individuals with non-zero expression 
numberOfIndividuals <- nrow(expression_df_nonzero) 
print(sprintf("Total number of individuals = %s",numberOfIndividuals)) 
numberOfIndividualsWithNonZeroExpression <- colSums(expression_df_nonzero[-1]!=0,na.rm=T) 

# Percentage of individuals with non-zero expression 
percentOfIndividualsWithNonZeroExpression <- (numberOfIndividualsWithNonZeroExpression/nrow(expression_df_nonzero))*100 

# Mean Expression 
meanExpression <-  colMeans(expression_df_nonzero[,-1],na.rm=T)

# Expression Variance 
expressionVariance <-  apply(expression_df_nonzero[,-1], MARGIN=2, FUN=var, na.rm=TRUE)  
 
# Combine percentage, mean and variance information 
extraInfo <- data.frame(colnames(expression_df_nonzero[-1]), numberOfIndividuals, 
  numberOfIndividualsWithNonZeroExpression, 
  percentOfIndividualsWithNonZeroExpression, 
  meanExpression, expressionVariance) 
colnames(extraInfo)[1] = "Gene"

fwrite(extraInfo,file=sprintf("%s/%s_%s_extraInfo.tsv",output.dir, cellLabel, chrNumber),quote=F, sep="\t") 

print(sprintf("Number of genes: %s", nrow(extraInfo))) 
 
# Plot gene expression against percent of non-zero individuals 
lessthan20 <- extraInfo[(extraInfo$percentOfIndividualsWithNonZeroExpression<100),] 
#association <- ggplot(lessthan20, aes(x=percentOfIndividualsWithNonZeroExpression,y=meanExpression)) + 
#    geom_point() + 
#    theme_classic() 
#ggsave(association,file=sprintf("%s/%s/%s_%s_expression_lessthan20pctzero.png", output.dir, cohort_out, cellLabelout, chrNumber)) 
 
# Filter genes with less than 10 percent individuals with non-zero expression 
atleast10percent <- extraInfo[(extraInfo$percentOfIndividualsWithNonZeroExpression>10),] 
expression_df_nonzero <- expression_df_nonzero[-1][,colnames(expression_df_nonzero[-1])%in% rownames(atleast10percent)] 
expression_df_nonzero$sampleid <- sample_ids 
gene_ids <- colnames(expression_df_nonzero[-ncol(expression_df_nonzero)]) 
 
print(sprintf("Number of genes after filtering: %s", nrow(atleast10percent))) 

# Prepare covariate file 
covariate_ids <- colnames(covariate_df)[-1] 
covariate_df$sampleid = covariate_df$IID

# Find Gene-SNP pairs for chrNumber 
geneLoc_df <- geneLoc_df[(geneLoc_df$hgnc_symbol %in% gene_ids),] 
geneLoc_df$left <- geneLoc_df$start_position - 500000 
geneLoc_df$right <- geneLoc_df$end_position + 500000 
geneLoc_df$geneid = geneLoc_df$hgnc_symbol

if (length(sample_ids) < 10){
  print("Too few samples")
 # quit()
}
genes = unique(geneLoc_df$geneid)

all_res = data.frame()
start1 = Sys.time()
for (i in 1:length(genes)){
  print(paste0(i,"---",length(genes)))
  
  gen=genes[i]

  # Write pheno file, snps and covariate files
  phenofile = expression_df_nonzero[,c("sampleid",gen)]
  colnames(phenofile)[1] = "IID"
  phenofile = phenofile[is.na(phenofile[,gen]) == F,]
  phenofile$test = RankNorm(phenofile[,gen])
  phenofile[,gen] = phenofile$test; phenofile$test = NULL
  temp_phenofile=paste0(output.dir,"/temp_phenofile_",gen,"_",cellLabel,".txt")
  write.table(phenofile,temp_phenofile,sep="\t",row.names=F,quote=F)
  temp_indfile=paste0(output.dir,"/temp_indfile_",gen,"_",cellLabel,".txt")
  write.table(phenofile[,1],temp_indfile,sep=" ",row.names=F,quote=F,col.names=F)
  covfile = covariate_df[,c("IID","Age","Sex","C1","C2","C3","C4","C5","C6","C7","C8","PC1","PC2","PC3","PC4")]
  covfile[covfile$Sex == "M","Sex"] = 0;  covfile[covfile$Sex == "F","Sex"] = 1
  temp_covfile=paste0(output.dir,"/temp_covfile_",gen,"_",cellLabel,".txt")
  write.table(covfile,temp_covfile,sep="\t",row.names=F,quote=F)

  # Run analysis
  temp_resultsfile=paste0(output.dir,"/temp_resfile_",gen,"_",cellLabel,".txt")
  startbp = geneLoc_df[geneLoc_df$hgnc_symbol == gen,"StarteQTLWindow"]
  endbp = geneLoc_df[geneLoc_df$hgnc_symbol == gen,"EndeQTLWindow"]
  system(paste0("/xx/xxx/bin/plink2 --pfile ",genotype_filename, " --chr ",chrNumber, " --from-bp ",startbp, " --to-bp ", endbp, " --keep ",temp_indfile, " --pheno ",temp_phenofile," --covar ",temp_covfile," --covar-variance-standardize C1 C2 C3 C4 C5 C6 C7 C8 PC1 PC2 PC3 PC4"," --glm hide-covar --vif 50"," --maf 0.05 --out ",temp_resultsfile))
  
  # Read results
  if (file.exists(paste0(temp_resultsfile,".",gen,".glm.linear"))){
  res = read_tsv(paste0(temp_resultsfile,".",gen,".glm.linear"))
  res$GENE = gen
  if (nrow(all_res) == 0){
    all_res = res
  }else{
    all_res = rbind(all_res,res)
  }}


  # Remove temporary files
  file.remove(temp_phenofile)
  file.remove(temp_indfile)
  file.remove(temp_covfile)
  file.remove(paste0(temp_resultsfile,".",gen,".glm.linear"))
  file.remove(paste0(temp_resultsfile,".log"))
}
end1 = Sys.time()
print(end1-start1)
keep = all_res


all_res = separate(all_res,col="ID",sep=";",into=c("SNP1","SNP2","SNP3"),rem=F)
all_res = data.frame(all_res)
all_res$SNP = NA
all_res[grepl("rs",all_res$SNP3) == T & is.na(all_res$SNP3) == F,"SNP"] = all_res[grepl("rs",all_res$SNP3) == T  & is.na(all_res$SNP3) == F,"SNP3"]
all_res[grepl("rs",all_res$SNP2) == T & is.na(all_res$SNP2) == F,"SNP"] = all_res[grepl("rs",all_res$SNP2) == T  & is.na(all_res$SNP2) == F,"SNP2"]
all_res[grepl("rs",all_res$SNP1) == T & is.na(all_res$SNP1) == F,"SNP"] = all_res[grepl("rs",all_res$SNP1) == T  & is.na(all_res$SNP1) == F,"SNP1"]
all_res[,c("SNP1","SNP2","SNP3")] = NULL

# Save 
fwrite(all_res,sprintf("%s/%s_chr%s_correlation_results_plink.tsv", output.dir, cellLabel,chrNumber),sep="\t",quote=F) 

# Save needed tables for next rounds
fwrite(geneLoc_df,sprintf("%s/%s_chr%s_geneLoc_df.tsv", output.dir, cellLabel,chrNumber),sep="\t",quote=F) 
fwrite(covariate_df,sprintf("%s/%s_chr%s_covariate_df.tsv", output.dir, cellLabel,chrNumber),sep="\t",quote=F) 
fwrite(expression_df_nonzero,sprintf("%s/%s_chr%s_exprs_df.tsv", output.dir, cellLabel,chrNumber),sep="\t",quote=F) 

quit()


