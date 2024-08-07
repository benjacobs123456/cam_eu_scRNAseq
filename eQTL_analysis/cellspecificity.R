library(data.table)
library(plyr)
library(ggplot2)
library(tidyr)
library(gdata)
library(RNOmni)
library(gridExtra)
library(dplyr)



data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"

cell_types=c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")

datall = readRDS(paste0(data.dir,"/","Complete_results_all_celltypes.rds"))
datall = data.frame(datall)
datall$P = as.numeric(as.character(datall$P))
datall = datall[is.na(datall$P) == F,]
dat = datall[datall$source == "CSF",]
dat2 = datall[datall$source == "PBMC",]

eSNPs = fread(paste0(data.dir,"/All_fdr_eSNPs.txt"))



cts = cell_types
# CSF - CSF
all_res = data.frame()
for (i1 in 1:length(cts)){
	ct1 = cts[i1]
	print(ct1)

	res1 = eSNPs[eSNPs$cell_type == ct1 & eSNPs$source == "CSF",]
    res1$SNP_GENE = paste0(res1$SNP,"_",res1$GENE)
    res1 = data.frame(res1)

	

	# one eQTL per gene
	if (nrow(res1) > 0){

	# loop through other cell types
	for (i2 in 1:length(cts)){

		ct2 = cts[i2]
		print(paste0(ct1,"-",ct2))

		if (ct1 != ct2){
			res2 = dat[dat$cell_type == ct2,]
			if (nrow(res2) > 0){
			res2$SNP_GENE = paste0(res2$SNP,"_",res2$GENE)

			colnames(res2) = paste0(colnames(res2),"_","ct2")
			res2 = data.frame(res2)
			res2$SNP_GENE = res2[,paste0("SNP_GENE","_","ct2")]
			res2 = res2[res2$SNP_GENE %in% res1$SNP_GENE == T,]

			res = join(res1,res2,by="SNP_GENE",type="left")
			res = res[is.na(res$A1_ct2) == F,]
			res = data.frame(res)
			res[res$A1 != res$A1_ct2,"BETA_ct2"] = as.numeric(as.character(res[res$A1 != res$A1_ct2,"BETA_ct2"])) * -1
			res$het_z = (as.numeric(as.character(res$BETA))- as.numeric(as.character(res[,paste0("BETA","_","ct2")]))) / ( sqrt(as.numeric(as.character(res$SE))^2  + as.numeric(as.character(res[,paste0("SE","_","ct2")]))^2))
   			res$het_p = 1 - pnorm(abs(res$het_z))
   			res$cts = paste0(ct1,"_",ct2)
   			all_res = rbind(all_res,res)
   		}

		}
	}
}

}
all_res = all_res[is.na(all_res$het_p) == F,]
all_res1 = all_res

#PBMC - PBMC
all_res2 = data.frame()

for (i1 in 1:length(cts)){
	ct1 = cts[i1]
	print(ct1)

	res1 = eSNPs[eSNPs$cell_type == ct1 & eSNPs$source == "PBMC",]
    res1$SNP_GENE = paste0(res1$SNP,"_",res1$GENE)
    res1 = data.frame(res1)
	
	# one eQTL per gene
	if (nrow(res1) > 0){

	# loop through other cell types
	for (i2 in 1:length(cts)){

		ct2 = cts[i2]
		print(paste0(ct1,"-",ct2))

		if (ct1 != ct2){

			res2 = dat2[dat2$cell_type == ct2,]
			if (nrow(res2) > 0){
			res2$SNP_GENE = paste0(res2$SNP,"_",res2$GENE)

			colnames(res2) = paste0(colnames(res2),"_","ct2")
			res2 = data.frame(res2)
			res2$SNP_GENE = res2[,paste0("SNP_GENE","_","ct2")]
			res2 = res2[res2$SNP_GENE %in% res1$SNP_GENE == T,]

			res = join(res1,res2,by="SNP_GENE",type="left")
			res = res[is.na(res$A1_ct2) == F,]
			res = data.frame(res)
			res[res$A1 != res$A1_ct2,"BETA_ct2"] = as.numeric(as.character(res[res$A1 != res$A1_ct2,"BETA_ct2"])) * -1
			res$het_z = (as.numeric(as.character(res$BETA))- as.numeric(as.character(res[,paste0("BETA","_","ct2")]))) / ( sqrt(as.numeric(as.character(res$SE))^2  + as.numeric(as.character(res[,paste0("SE","_","ct2")]))^2))
   			res$het_p = 1 - pnorm(abs(res$het_z))
   			res$cts = paste0(ct1,"_",ct2)
   			all_res2 = rbind(all_res2,res)
   		}

		}
	}
}

}
all_res2 = all_res2[is.na(all_res2$het_p) == F,]


# CSF - PBMC
all_res3 = data.frame()
for (i1 in 1:length(cts)){
	ct1 = cts[i1]
	print(ct1)

	res1 = eSNPs[eSNPs$cell_type == ct1 & eSNPs$source == "CSF",]
    res1$SNP_GENE = paste0(res1$SNP,"_",res1$GENE)
    res1 = data.frame(res1)
	
	

	# one eQTL per gene
	if (nrow(res1) > 0){

	# loop through other cell types
	for (i2 in 1:length(cts)){

		ct2 = cts[i2]
		print(paste0(ct1,"-",ct2))

		if (ct1 != ct2){
			res2 = dat2[dat2$cell_type == ct2,]
			if (nrow(res2) > 0){
			res2$SNP_GENE = paste0(res2$SNP,"_",res2$GENE)

			colnames(res2) = paste0(colnames(res2),"_","ct2")
			res2 = data.frame(res2)
			res2$SNP_GENE = res2[,paste0("SNP_GENE","_","ct2")]
			res2 = res2[res2$SNP_GENE %in% res1$SNP_GENE == T,]

			res = join(res1,res2,by="SNP_GENE",type="left")
			res = res[is.na(res$A1_ct2) == F,]
			res = data.frame(res)
			res[res$A1 != res$A1_ct2,"BETA_ct2"] = as.numeric(as.character(res[res$A1 != res$A1_ct2,"BETA_ct2"])) * -1			
			res$het_z = (as.numeric(as.character(res$BETA))- as.numeric(as.character(res[,paste0("BETA","_","ct2")]))) / ( sqrt(as.numeric(as.character(res$SE))^2  + as.numeric(as.character(res[,paste0("SE","_","ct2")]))^2))
   			res$het_p = 1 - pnorm(abs(res$het_z))
   			res$cts = paste0(ct1,"_",ct2)
   			all_res3 = rbind(all_res3,res)
   		}

		}
	}
}

}
all_res3 = all_res3[is.na(all_res3$het_p) == F,]

#PBMC - CSF
all_res4 = data.frame()
for (i1 in 1:length(cts)){
	ct1 = cts[i1]
	print(ct1)

	res1 = eSNPs[eSNPs$cell_type == ct1 & eSNPs$source == "PBMC",]
    res1$SNP_GENE = paste0(res1$SNP,"_",res1$GENE)
    res1 = data.frame(res1)
	
	
	# one eQTL per gene
	if (nrow(res1) > 0){
	
	# loop through other cell types
	for (i2 in 1:length(cts)){

		ct2 = cts[i2]
		print(paste0(ct1,"-",ct2))

		if (ct1 != ct2){

			res2 = dat[dat$cell_type == ct2,]
			if (nrow(res2) > 0){
			res2$SNP_GENE = paste0(res2$SNP,"_",res2$GENE)

			colnames(res2) = paste0(colnames(res2),"_","ct2")
			res2 = data.frame(res2)
			res2$SNP_GENE = res2[,paste0("SNP_GENE","_","ct2")]
			res2 = res2[res2$SNP_GENE %in% res1$SNP_GENE == T,]

			res = join(res1,res2,by="SNP_GENE",type="left")
			res = res[is.na(res$A1_ct2) == F,]
			res = data.frame(res)
			res[res$A1 != res$A1_ct2,"BETA_ct2"] = as.numeric(as.character(res[res$A1 != res$A1_ct2,"BETA_ct2"])) * -1		
			res$het_z = (as.numeric(as.character(res$BETA))- as.numeric(as.character(res[,paste0("BETA","_","ct2")]))) / ( sqrt(as.numeric(as.character(res$SE))^2  + as.numeric(as.character(res[,paste0("SE","_","ct2")]))^2))
   			res$het_p = 1 - pnorm(abs(res$het_z))
   			res$cts = paste0(ct1,"_",ct2)
   			all_res4 = rbind(all_res4,res)
   		}

		}
	}
}

}
all_res4 = all_res4[is.na(all_res4$het_p) == F,]



all_res = rbind(all_res1,all_res2,all_res3,all_res4)

all_res$P = as.numeric(as.character(all_res$P))
all_res$P_ct2 = as.numeric(as.character(all_res$P_ct2))
all_res = all_res[order(all_res$het_p),]

fwrite(all_res,paste0(data.dir,"/cellspecificeffects.tsv"))



# Correlation of Z scores
test = all_res
test$BETA = as.numeric(as.character(test$BETA))
test$BETA_ct2 = as.numeric(as.character(test$BETA_ct2))
test$SE = as.numeric(as.character(test$SE))
test$SE_ct2 = as.numeric(as.character(test$SE_ct2))
test$Z = test$BETA/test$SE
test$Z_ct2 = test$BETA_ct2/test$SE_ct2
print("Correlation of Z scores:")
print(cor.test(test$Z,test$Z_ct2)$estimate)

# eQTLs with evidence for cell type specificity 
# Defitionion: select eQTLs with FDR < 0.05 & P > 0.05 in all other cell types  and with effect in same cell type but other source
all_res$het_p_fdr = p.adjust(all_res$het_p,method="fdr")
test = all_res
test$keep = 0
compall = datall[datall$GENE %in% test$GENE == T & datall$P < 0.05,]
compall = compall[compall$ID %in% test$ID == T,]
test = data.frame(test)
for (i in 1:nrow(test)){
	if(test[i,"het_p_fdr"] < 0.05){
		ct1 = test[i,"cell_type"]
		source1 = test[i,"source"]
		comp = test[test$GENE == test[i,"GENE"] & test$ID == test[i,"ID"] & test$cell_type_ct2 != ct1,]
		if (nrow(comp) > 0){
		if (ct1 %in% c("B.cells","Plasma.cells") == T) {comp = comp[comp$cell_type_ct2 %in% c("Plasma.cells","B.cells") == F,]}
		if (ct1 %in% c("CD4.T.cells","Tregs","CD8.T.cells") == T){comp = comp[comp$cell_type_ct2 %in% c("CD4.T.cells","Tregs","CD8.T.cells") == F,]}
	}
		if(nrow(comp[comp$P_ct2 < 0.01,]) == 0){
			
			#if(nrow(compall[compall$cell_type== ct1 & compall$source != source1 & compall$GENE == test[i,"GENE"] & compall$ID == test[i,"ID"],]) > 0){
				test[i,"keep"] = 1
			#}
		}
	}
}
keep1 = test
test = test[test$keep == 1,]



# Select some interesting results
novel = fread(paste0(data.dir,"/comparison_with_references.txt"))
test$Comb = paste0(test$SNP,"_",test$GENE)
novel$Comb = paste0(novel$SNP,"_",novel$GENE)
test = join(test,novel,by="Comb",type="left")
test = unique(test)
testn = test[test$siginGtex == 0 & (test$sigeQTLGen == 0 | is.na(test$sigeQTLGen) == T) & is.na(test$SNP) == F,]
refall = fread("/xx/xxx/SC/CAM_TUM/eqtl/references/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",select=c("SNP","GeneSymbol","Pvalue","SNPChr","SNPPos"))
ref = fread("/xx/xxx/SC/CAM_TUM/eqtl/references/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",select=c("SNP","GeneSymbol","Pvalue","SNPChr","SNPPos"))
yazar = read.csv2("/xx/xxx/SC/CAM_TUM/eqtl/references/yazar.csv")


refallk = refall
refall = refall[refall$P < 0.05,]
refalln= refallk[refallk$P >= 0.05,]
refall$Comb = paste0(refall$SNP,"_",refall$GeneSymbol)
refalln$Comb = paste0(refalln$SNP,"_",refalln$GeneSymbol)
testna = testn[testn$Comb %in% refall$Comb == F & is.na(testn$SNP) == F,]
refalln[refalln$Comb %in% testna$Comb == T,]
#testna = unique(testna[,1:24])

test = unique(test)

#####

eSNPs = fread(paste0(data.dir,"/All_fdr_eSNPs_permutation_results.txt"))
eSNPs$Comb = paste0(eSNPs$SNP,"_",eSNPs$GENE)
eSNPsperm = eSNPs[eSNPs$perm_sig == 1,]
test = test[test$Comb %in% eSNPsperm$Comb == T | (test$sigeQTLGen == T & is.na(test$sigeQTLGen) == F) | (test$siginGtex == T & is.na(test$siginGtex) == F) | (test$sigYazar == T & is.na(test$sigYazar) == F),]

# Interesing findings
test = test[test$cell_type %in% c("CD4.T.cells","CD8.T.cells") == T | (test$Comb %in% testn$Comb == T & test$Comb %in% testna$Comb == T),  ]

print("Number of unique combinations of celltype, snp, gene with evidence for cell type specificity: ")
test$comb = paste0(test$cell_type,"_",test$ID,"_",test$GENE)
print(length(unique(test$comb)))

inf = unique(test[,c("cell_type","ID","GENE"),])
inftab = data.frame(table(inf$cell_type)); inftab = inftab[order(inftab$Freq,decreasing=T),]
inftab$total = nrow(inf); inftab$per = round(inftab$Freq/inftab$total*100,1)
print("Cell type distribution:")
print(inftab)

# Prepare supplementary table
out = test
out = join(out,eSNPs[,c("Comb","perm_p")],by="Comb",type="left")
out = out[,c("X.CHROM","POS","SNP","GENE","source","cell_type","BETA","SE","P","FDR","perm_p")]
out = unique(out)
out[is.na(out$perm_p) == T,"perm_p"] = "previously reported"
colnames(out) = c("Chr","Pos","SNP","Gene","Source","Cell_type","Beta","SE","p","FDR p","permutation p")
write.csv(out,paste0(data.dir,"/tables/supplementary_table_cellspecific_eqtls.csv"))


#Figures

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

# Read control summary statistics
datcsf = datall[datall$source == "CSF",]
datpbmc = datall[datall$source == "PBMC",]


celltypes = c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")


source("/xx/xxx/SC/CAM_TUM/eqtl/Scripts/eqtl_plot_function_All.R")


out = out[out$Gene != "XCL1",] # can't plot this one

myplots = list()

for (i in 1:nrow(out)){
  snp=as.character(out[i,"SNP"])
  gene=as.character(out[i,"Gene"])
  
  ctout1=as.character(out[i,"Cell_type"])
  if (ctout1 == "CD4.T.cells"){ctout2 = "CD8.T.cells";ctout3 = "B.cells"}
  if (ctout1 != "CD4.T.cells"){
    if (ctout1 == "B.cells"){
      ctout2 = "CD4.T.cells";ctout3 = "CD14.Mono"
    }else{
      ctout2 = "CD4.T.cells";ctout3 = "B.cells"
    }
    }
  source = as.character(out[i,"Source"])
  cohort = "All"
  round=1
  plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,cell_types=celltypes,window=500000)

  myplots[[length(myplots)+1]] = plots

}


dir.create(paste0(data.dir,"/figures/cellspecificity"))

for (i in 1:nrow(out)){

   print(i)

  snp=as.character(out[i,"SNP"])
  gene=as.character(out[i,"Gene"])
  ctout1=as.character(out[i,"Cell_type"])
  source = as.character(out[i,"Source"])
  cohort = "All"

  png(paste0(data.dir,"/figures/cellspecificity/",gene,"_",snp,"_",ctout1,"_",source,"_",cohort,".png"),width=12,height=15,units="in",res=300)
  grid.arrange(myplots[[i]][[1]],myplots[[i]][[4]],myplots[[i]][[2]],myplots[[i]][[5]],myplots[[i]][[3]],myplots[[i]][[6]],myplots[[i]][[13]],ncol=4,layout_matrix=cbind(c(1,5),c(2,6),c(3,7),c(4,7)))
  dev.off()


}
