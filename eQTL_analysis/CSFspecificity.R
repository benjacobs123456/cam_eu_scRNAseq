library(data.table)
library(plyr)


data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"

cell_types=c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")


dat = readRDS(paste0(data.dir,"/","Complete_results_all_celltypes.rds"))
eSNPs = fread(paste0(data.dir,"/All_fdr_eSNPs.txt"))

dat = data.frame(dat)
dat$P = as.numeric(as.character(dat$P))
dat = dat[is.na(dat$P) == F,]
dat$SNP_GENE = paste0(dat$SNP,"_",dat$GENE)




cts = unique(dat$cell_type)

all_res = data.frame()
for (i1 in 1:length(cell_types)){
	ct1 = cell_types[i1]
	print(ct1)

  res1 = eSNPs[eSNPs$cell_type == ct1 & eSNPs$source == "CSF",]
  res1$SNP_GENE = paste0(res1$SNP,"_",res1$GENE)
  res1 = data.frame(res1)

  if (nrow(res1) > 0){

	# PBMC
	res2 = dat[dat$SNP_GENE %in% res1$SNP_GENE == T & dat$source == "PBMC",]
	if (nrow(res2) > 0){

	colnames(res2) = paste0(colnames(res2),"_","PBMC")
	res2 = data.frame(res2)
	res2$SNP_GENE = res2[,paste0("SNP_GENE","_","PBMC")]

	res = join(res1,res2,by="SNP_GENE",type="left")
	res = res[is.na(res$A1_PBMC) == F,]
  res = data.frame(res)
	res[res$A1 != res$A1_PBMC,"BETA_PBMC"] = as.numeric(as.character(res[res$A1 != res$A1_PBMC,"BETA_PBMC"])) * -1
	res = data.frame(res)
	res$het_z = (as.numeric(as.character(res$BETA))- as.numeric(as.character(res[,paste0("BETA","_","PBMC")]))) / ( sqrt(as.numeric(as.character(res$SE))^2  + as.numeric(as.character(res[,paste0("SE","_","PBMC")]))^2))
   	res$het_p = 1 - pnorm(abs(res$het_z))
   	all_res = rbind(all_res,res)
  
   		}
	
}

}


all_res2 = data.frame()
for (i1 in 1:length(cell_types)){
	ct1 = cell_types[i1]
	print(ct1)

	res1 = eSNPs[eSNPs$cell_type == ct1 & eSNPs$source == "PBMC",]
  res1$SNP_GENE = paste0(res1$SNP,"_",res1$GENE)
  res1 = data.frame(res1)


  if (nrow(res1) > 0){

	# CSF
	res2 = dat[dat$SNP_GENE %in% res1$SNP_GENE == T & dat$source == "CSF",]
	if (nrow(res2) > 0){

	colnames(res2) = paste0(colnames(res2),"_","CSF")
	res2$SNP_GENE = res2[,paste0("SNP_GENE","_","CSF")]

	res = join(res1,res2,by="SNP_GENE",type="left")
	res = res[is.na(res$A1_CSF) == F,]
  res = data.frame(res)
	res[res$A1 != res$A1_CSF,"BETA_CSF"] = as.numeric(as.character(res[res$A1 != res$A1_CSF,"BETA_CSF"])) * -1
	res = data.frame(res)
	res$het_z = (as.numeric(as.character(res$BETA))- as.numeric(as.character(res[,paste0("BETA","_","CSF")]))) / ( sqrt(as.numeric(as.character(res$SE))^2  + as.numeric(as.character(res[,paste0("SE","_","CSF")]))^2))
   	res$het_p = 1 - pnorm(abs(res$het_z))
   	all_res2 = rbind(all_res2,res)
   		}
	
}

}

all_res = all_res[is.na(all_res$het_p) == F,]
all_res2 = all_res2[is.na(all_res2$het_p) == F,]
all_res$source1 = "CSF"; all_res$source2 = "PBMC"
all_res2$source1 = "PBMC"; all_res2$source2 = "CSF"
colnames(all_res) = gsub("_PBMC","_source2",colnames(all_res))
colnames(all_res2) = gsub("_CSF","_source2",colnames(all_res2))
all_res = rbind(all_res,all_res2)

all_res$P = as.numeric(as.character(all_res$P))
all_res$P_source2 = as.numeric(as.character(all_res$P_source2))
all_res = all_res[order(all_res$het_p),]

fwrite(all_res,paste0(data.dir,"/","sourcespecificeffects.tsv"))



# Correlation of Z scores
# Do this only test with more than 30 observations - and only comparing the same cell types
#test = all_res[all_res$OBS_CT > 30 & all_res$OBS_CT_source2 > 30,]
test = all_res
test$BETA = as.numeric(as.character(test$BETA))
test$BETA_source2 = as.numeric(as.character(test$BETA_source2))
test$SE = as.numeric(as.character(test$SE))
test$SE_source2 = as.numeric(as.character(test$SE_source2))
test$Z = test$BETA/test$SE
test$Z_source2 = test$BETA_source2/test$SE_source2
test = test[test$cell_type == test$cell_type_source2,]
print("Correlation of Z scores:")
print(cor.test(test$Z,test$Z_source2)$estimate)
plot(test$Z,test$Z_source2)

# eQTLs with evidence for source specificity 
# Defitionion: select eQTLs with FDR_ct < 0.05 & P > 0.01 in all cell types in the other source
all_res$het_p_fdr = p.adjust(all_res$het_p,method="fdr")
test = all_res
test = data.frame(test)
test$keep = 0
for (i in 1:nrow(test)){
  print(i)
	if(test[i,"het_p_fdr"] < 0.05){
		source1 = test[i,"source"]
		comp = test[test$GENE == test[i,"GENE"] & test$ID == test[i,"ID"] & test$source_source2 != source1,]
    comp2 = dat[dat$GENE == test[i,"GENE"] & dat$ID == test[i,"ID"] & dat$source != source1,]
		if(nrow(comp[comp$P_source2 < 0.01,]) == 0 & nrow(comp2[comp2$P < 0.01,]) == 0 ){
			test[i,"keep"] = 1
		}
	}
}
keep1 = test
test = test[test$keep == 1,]
test = test[test$cell_type == test$cell_type_source2,]
test = test[order(test$het_p),]



print("Number of unique combinations of celltype, snp, gene with evidence for source type specificity: ")
test$comb = paste0(test$cell_type,"_",test$ID,"_",test$GENE)
print(length(unique(test$comb)))

inf = unique(test[,c("cell_type","ID","GENE"),])
inftab = data.frame(table(inf$cell_type)); inftab = inftab[order(inftab$Freq,decreasing=T),]
inftab$total = nrow(inf); inftab$per = round(inftab$Freq/inftab$total*100,1)
print("Cell type distribution:")
print(inftab)
inf = unique(test[,c("cell_type","ID","GENE","source"),])
inftab2 = data.frame(table(inf$source)); inftab2 = inftab2[order(inftab2$Freq,decreasing=T),]
inftab2$total = nrow(inf); inftab2$per = round(inftab2$Freq/inftab2$total*100,1)
print("Source distribution:")
print(inftab2)



# Select some interesting results
novel = fread(paste0(data.dir,"/comparison_with_references.txt"))
novel = data.frame(novel)
test1 = test[test$cell_type == test$cell_type_source2,]
novel$Comb = paste0(novel$SNP,"_",novel$GENE)
test1$Comb = paste0(test1$SNP,"_",test1$GENE)
novel[novel$sigeQTLGen == T & is.na(novel$sigeQTLGen) == F,"sigeQTLGen"] = 1
novel[novel$sigeQTLGen == F & is.na(novel$sigeQTLGen) == F,"sigeQTLGen"] = 0
test1 = join(test1,novel,by="Comb",type="left")
test1_pbmc = test1
test1 = test1[test1$siginGtex == 0 & (test1$sigeQTLGen == 0 | is.na(test1$sigeQTLGen) == T) & test1$sigYazar == 0 & is.na(test1$SNP) == F,]
test1 = test1[is.na(test1$Source) == F,]
refall = fread("/xx/xxx/SC/CAM_TUM/eqtl/references/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",select=c("SNP","GeneSymbol","Pvalue","SNPChr","SNPPos"))
ref = fread("/xx/xxx/SC/CAM_TUM/eqtl/references/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",select=c("SNP","GeneSymbol","Pvalue","SNPChr","SNPPos"))
yazar = read.csv2("/xx/xxx/SC/CAM_TUM/eqtl/references/yazar.csv")

refallk = refall
refall = refall[refall$P < 0.05,]
refalln= refallk[refallk$P >= 0.05,]
refall$Comb = paste0(refall$SNP,"_",refall$GeneSymbol)
refalln$Comb = paste0(refalln$SNP,"_",refalln$GeneSymbol)
test1 = test1[test1$Comb %in% refall$Comb == F & is.na(test1$SNP) == F,]
refalln[refalln$Comb %in% test1$Comb == T,]

test1 = unique(test1)
test1_csf = test1[test1$source == "CSF",]

test1_pbmc = test1_pbmc[test1_pbmc$source == "PBMC",]
test1_pbmc = unique(test1_pbmc)

#Figures

genlist=read.table("/xx/xxx/SC/CAM_TUM/eqtl/gene_list.txt",h=T)
fam_filename = "/xx/xxx/SC/CAM_TUM/eqtl/genotype_data/merged/merged_1to22_SampleQC_VariantQC2.psam"
meta_filename="/xx/xxx/SC/CAM_TUM/datasets/gex_metadata_new.csv"
genIDsTUM_filename="/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs.csv"
genIDsCAM_filename="/xx/xxx/SC/CAM_TUM/sc_for_tum/GT_ID.txt"
genotype_filename="/xx/xxx/SC/CAM_TUM/eqtl/genotype_data/merged/merged_1to22_SampleQC_VariantQC"

datcsf = dat[dat$source == "CSF",]
datpbmc = dat[dat$source == "PBMC",]

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

celltypes = c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")


# Check which CSF specific results survive permutation
csf = test1_csf
eSNPs = fread(paste0(data.dir,"/All_fdr_eSNPs_permutation_results.txt"))
csf$Comb = paste0(csf$SNP,"_",csf$GENE)
eSNPs$Comb = paste0(eSNPs$SNP,"_",eSNPs$GENE)
eSNPsperm = eSNPs[eSNPs$perm_sig == 1,]
csf = csf[csf$Comb %in% eSNPsperm$Comb == T,]
#csf = csf[,c(1:79,128:134)]

# Prepare supplementary table
out = csf
out = join(out,eSNPs[,c("Comb","perm_p")],by="Comb",type="left")
out = out[,c("X.CHROM","POS","SNP","GENE","cell_type","BETA","SE","P","FDR","perm_p")]
colnames(out) = c("Chr","Pos","SNP","Gene","Cell_type","Beta","SE","p","FDR p","permutation p")
write.csv(out,paste0(data.dir,"/tables/supplementary_table_CSFspecific_eqtls.csv"))


# Manually check eQTL window for those
ref = fread(paste0(data.dir,"/../references/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"))
ref = ref[ref$GeneSymbol %in% csf[csf$window_eQTLGen == 1,"GENE"] == T,]
gtex = data.frame()
for (i in 1:nrow(csf)){
  if (csf[i,"window_eQTLGen"] == 1){
    out = c(csf[i,"SNP"],ref[ref$GeneSymbol == csf[i,"GENE"],"SNP"])
    out = c(out,as.list(yazar[yazar$Gene.ID == csf[i,"GENE"],"SNP"]))

    write.table(out,paste0(data.dir,"/remove_",csf[i,"GENE"],"_eQTLSNPs.txt"),sep="\t",quote=F,row.names=F,col.names=F)
  }
}


# Figures CSF
source(paste0(data.dir,"/Scripts/eqtl_plot_function.R"))

myplots_csf = list()

for (i in 1:nrow(csf)){
  print(i)
  snp=as.character(csf[i,"SNP"])
  gene=as.character(csf[i,"GENE"])
  
  ctout1=as.character(csf[i,"cell_type"])
  if (ctout1 == "CD4.T.cells"){ctout2 = "CD8.T.cell";ctout3 = "B.cells"}
  if (ctout1 != "CD4.T.cells"){
    if (ctout1 == "B.cells"){
      ctout2 = "CD4.T.cells";ctout3 = "CD8.T.cells"
    }else{
      ctout2 = "CD4.T.cells";ctout3 = "B.cells"
    }
    }
  source = as.character(csf[i,"source"])
  cohort = as.character(csf[i,"cohort"])
  round=1
  plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,cell_types=celltypes,window=500000)

  myplots_csf[[length(myplots_csf)+1]] = plots

}

dir.create(paste0(data.dir,"/figures/CSFspecificity"))
dir.create(paste0(data.dir,"/figures/CSFspecificity/CSF"))
dir.create(paste0(data.dir,"/figures/CSFspecificity/PBMC"))

for (i in 1:nrow(csf)){

  snp=as.character(csf[i,"SNP"])
  gene=as.character(csf[i,"GENE"])
  ctout1=as.character(csf[i,"cell_type"])
  source = as.character(csf[i,"source"])

  png(paste0(data.dir,"/figures/CSFspecificity/CSF/",gene,"_",snp,"_",ctout1,"_",source,".png"),width=8,height=8,units="in",res=300)
  grid.arrange(myplots_csf[[i]][[1]],myplots_csf[[i]][[4]],myplots_csf[[i]][[7]],myplots_csf[[i]][[10]],myplots_csf[[i]][[13]],nrow=3)
  dev.off()

}

