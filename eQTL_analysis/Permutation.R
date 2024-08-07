library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Matrix.utils)
library(reshape2)
library(data.table)
library(plyr)
library(gridExtra)
library(ggrepel)

args = commandArgs(TRUE)

rown = args[1]

# directories
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
setwd(data.dir)

# Define function
getPermutationSingle <- function(PRS, modRes, nPerm=100000, p1sided=T, exp_effect="positive") # Permutation test for a single PRS
{

  form <- as.formula(paste0("residuals ~ ",PRS)) # test statistics for the correct model
  coefs <- coef(summary(lm(form, data=modRes)))
  original_p <- coefs[2,4]
  if(p1sided)
  {
    original_p <- original_p/2
    if(exp_effect=="positive" & coefs[2,1]<0) original_p <- 1-original_p else if(exp_effect=="negative" & coefs[2,1]>0) original_p <- 1-original_p  
  }
  
  if(any(is.na(modRes$residuals))) temp <- modRes[-which(is.na(modRes$residuals)),] else temp <- modRes 
  for (i in seq(1,nPerm)) # repetition of analyses using random sampling
  {
    temp$residuals <- sample(temp$residuals,nrow(temp))
    coefs <- coef(summary(lm(form, data=temp)))
    temp_p <- coefs[2,4]
    if(p1sided)
    {
      temp_p <- temp_p/2
      if(exp_effect=="positive" & coefs[2,1]<0) temp_p <- 1-temp_p else if(exp_effect=="negative" & coefs[2,1]>0) temp_p <- 1-temp_p  
    }
    if (i==1) {allP <- temp_p} else {allP <- c(allP, temp_p)}  
  }
  pSmaller <- length(allP[which(allP<=original_p)]) # calculation of permutation p-value
  finalP <- pSmaller/nPerm
  if(finalP==0) finalP <- paste0("<",1/nPerm)
  return(finalP)
}


genotype_filename="/xx/xxx/SC/CAM_TUM/eqtl/genotype_data/merged/merged_1to22_SampleQC_VariantQC"
 

# Get eGenes and eSTNPs
eSNPs = fread(paste0(data.dir,"/All_fdr_eSNPs.txt"))
eSNPs = unique(eSNPs)
eSNPs$Comb = paste0(eSNPs$source,"__",eSNPs$cohort,"__",eSNPs$SNP,"__",eSNPs$GENE)

# Only do this for "novel" hits
novel = fread(paste0(data.dir,"/comparison_with_references.txt"))
novel = data.frame(novel)
novel$Comb = paste0(novel$Source,"__",novel$Cohort,"__",novel$SNP,"__",novel$GENE)
novel = novel[novel$Comb %in% eSNPs$Comb == T,]
novel = novel[(novel$sigeQTLGen != T | is.na(novel$sigeQTLGen) == T) & (novel$siginGtex == 0 | is.na(novel$siginGtex) == T) & (novel$sigYazar == 0 | is.na(novel$sigYazar) == T),]
novel = novel[(novel$nomeQTLGen != 1 | is.na(novel$nomeQTLGen) == T),]
eSNPs = eSNPs[eSNPs$Comb %in% novel$Comb == T,]

eSNPs = data.frame(eSNPs)
# Run Permutation
run_permutation <- function(i){
	snp = as.character(eSNPs[i,"SNP"])
	gene = as.character(eSNPs[i,"GENE"])
	source = as.character(eSNPs[i,"source"])
	ct = as.character(eSNPs[i,"cell_type"])
	chr = as.character(eSNPs[i,"X.CHROM"])
	pos= as.character(eSNPs[i,"POS"])
	cohort="All"

	 # Read relevant files
 	 exprs_df = data.frame(fread(paste0(data.dir,"/../",source,"/",cohort,"/",ct,"_chr",chr,"_exprs_df.tsv")))
 	 exprs_df = exprs_df[,c(gene,"sampleid")]; colnames(exprs_df)[2] = "IID"
  	 temp_covfile = paste0(data.dir,"/../",source,"/",cohort,"/",ct,"_chr",chr,"_covariate_df.tsv")
  	 covariate_df = fread(temp_covfile)

  	 # Read genetic data
  	temp_file=paste0("test_",gene,"_",source,"_",cohort,"_",chr,"_",pos,"_",ct)
    system(paste0("/xx/xxx/bin/plink --bfile ",genotype_filename, " --chr ",chr, " --from-bp ",pos, " --to-bp ", pos, " --recode oxford --out ",temp_file))
	dos = read.table(paste0(temp_file,".gen"))
	dossam = read.table(paste0(temp_file,".sample")); dossam = dossam[-1:-2,]
	file.remove(paste0(temp_file,".gen"))
	file.remove(paste0(temp_file,".sample"))
	file.remove(paste0(temp_file,".nosex"))
	file.remove(paste0(temp_file,".log"))
	dosn = data.frame(dos[,1:5])
	ind = 1
	c = 6
	while(c <= ncol(dos)){
	  individual = dossam[ind,"V2"]
	  for (r in 1:nrow(dos)){
	  if (dos[r,c] == 1){
	    dosn[r,individual] = 2
	  }
	  if (dos[r,c+1] == 1){
	    dosn[r,individual] = 1
	  }
	  if (dos[r,c+2] == 1){
	    dosn[r,individual] = 0
	  }
	}
	  ind = ind + 1 
	  c = c + 3
	}
	dos = dosn[grepl(paste0(snp,";"),dosn$V2) == T,]
	dos = data.frame(t(dos))
	minallele = dos[4,1]
	majallele = dos[5,1]
	dos = dos[-1:-5,]
	dos = data.frame(snp=dos)
	rownames(dos) = colnames(dosn)[6:ncol(dosn)]
	dos$IID = rownames(dos)

	# Prepare data
	dataset = merge(exprs_df,dos,by="IID")
	dataset = merge(dataset,covariate_df,by="IID")
 
	# Get residuals
	dataset = data.frame(dataset)
	dataset$gene = dataset[,gene]
	dataset = dataset[is.na(dataset[,gene]) == F,]
	test = lm(gene ~ as.numeric(snp) + as.factor(Sex) + as.numeric(Age) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(PC1) + as.numeric(PC2) + as.numeric(PC3) + as.numeric(PC4) ,data=dataset)
	pvalorig = summary(test)$coefficients[2,4]
	print(pvalorig)
	dataset$residuals = residuals(lm(gene ~ as.factor(Sex) + as.numeric(Age) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(PC1) + as.numeric(PC2) + as.numeric(PC3) + as.numeric(PC4) ,data=dataset))

	# Run Permutation
	p_perm = NA
	roundUp <- function(x) {10^ceiling(log10(x))}
	nperm = roundUp(1/as.numeric(as.character(eSNPs[i,"P"])))
	print(nperm)
	if(nperm > 1000000){nperm = 1000000}
	#nperm=100000
	p_perm = getPermutationSingle("snp",dataset,nPerm=nperm,p1sided=F)

	out = data.frame(rownumber = i, p_permutation = p_perm)
	out

}

outtab = data.frame(run_permutation(rown))
print(outtab)
write.table(outtab,paste0(data.dir,"/Permutation/perm_",rown,".txt"),sep="\t",row.names=F,quote=F,col.names=F)

