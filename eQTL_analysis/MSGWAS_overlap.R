library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Matrix.utils)
library(reshape2)
library(data.table)
library(plyr)
library(gridExtra)
library(coloc)
library(ggrepel)

args = commandArgs(TRUE)

# directories
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
setwd(data.dir)
dir.create(paste0(data.dir,"/figures/colocMS/"))

cell_types=c("CD4.T.cells","CD8.T.cells","Tregs","B.cells","Plasma.cells")


# read in MS GWAS
ms_gwas = fread("/xx/shared_data/IMSGC_projects_summary_statistics/Science_2019/Discovery_phase_GWAS/meta-analysis/discovery_metav3.0.meta2")
ms_gwas = ms_gwas[is.na(ms_gwas$P) == F,]
ms_gwas = ms_gwas %>% mutate(chr_pos = paste0("chr",CHR,":",BP))

ms_gwas_hg38 = read.table("/xx/shared_data/IMSGC_projects_summary_statistics/Science_2019/Discovery_phase_GWAS/meta-analysis/discovery_metav3.0.meta2_hg38.bed",sep="\t")
ms_gwas_hg38$V4 = gsub("_",":",ms_gwas_hg38$V4)

# read in hg38 coordinates
ms_gwas_hg38 = ms_gwas %>%
  left_join(ms_gwas_hg38 %>% dplyr::rename(SNP = V4),by="SNP") %>%
  filter(!is.na(V3)) %>%
  filter(!is.na(P) & !is.na(OR)) %>%
  mutate(BP = V3) %>%
  dplyr::select(-V1,-V2,-V3,-V5) %>%
  mutate(SNP = ifelse(grepl("rs",SNP),SNP,paste0("chr",CHR,":",BP)))


ms_gwas_hg38 = ms_gwas_hg38 %>% mutate(chr_pos = paste0("chr",CHR,":",BP))
ms_gwas = ms_gwas_hg38


# Get all eGenes (in MS)
dat = readRDS(paste0(data.dir,"/","Complete_results.rds"))



all_genes = fread(paste0(data.dir,"/All_fdr_eSNPs.txt"))

# Select relevant eQTL results
dat$GENE_CT = paste0(dat$GENE,"_",dat$cell_type)
all_genes$GENE_CT = paste0(all_genes$GENE,"_",all_genes$cell_type)
dat = dat[dat$GENE_CT %in% all_genes$GENE_CT == T,]



all_genes$leadSNP = all_genes$SNP
all_genes_csf = all_genes[all_genes$source == "CSF",]
all_genes_pbmc = all_genes[all_genes$source == "PBMC",]
dat1 = dat[dat$source == "CSF",]
dat2 = dat[dat$source == "PBMC",]
dat1 = join(dat1,unique(all_genes_csf[,c("GENE_CT","leadSNP")]),by="GENE_CT",type="left")
dat2 = join(dat2,unique(all_genes_pbmc[,c("GENE_CT","leadSNP")]),by="GENE_CT",type="left")
dat1 = dat1[is.na(dat1$leadSNP) == F,]
dat2 = dat2[is.na(dat2$leadSNP) == F,]
snps1 = unique(dat1[,c("SNP","POS")]); snps1$leadSNP = snps1$SNP; snps1$leadPOS = snps1$POS
snps2 = unique(dat2[,c("SNP","POS")]); snps2$leadSNP = snps2$SNP; snps2$leadPOS = snps2$POS
dat1 = join(dat1,snps1[,c("leadSNP","leadPOS")],by="leadSNP",type="left")
dat2 = join(dat2,snps2[,c("leadSNP","leadPOS")],by="leadSNP",type="left")
dat1$CHR_POS = paste0("chr",dat1$X.CHROM,":",dat1$POS)
dat2$CHR_POS = paste0("chr",dat2$X.CHROM,":",dat2$POS)
dat1 = dat1[dat1$CHR_POS %in% ms_gwas$chr_pos == T,]
dat2 = dat2[dat2$CHR_POS %in% ms_gwas$chr_pos == T,]
dat1 = data.frame(dat1)
dat2 = data.frame(dat2)

all_genes = data.frame(all_genes)



# Colocalization
coloc_function <- function (x) {
	
  ct = all_genes[x,"cell_type"]
	gene = all_genes[x,"GENE"]
  source = all_genes[x,"source"]
  cohort = all_genes[x,"cohort"]
  leadSNP = all_genes[x,"leadSNP"]
 
	if (source == "CSF"){res = dat1[dat1$GENE == gene & dat1$cell_type == ct & dat1$leadSNP == as.character(leadSNP),]}
  if (source == "PBMC"){res = dat2[dat2$GENE == gene & dat2$cell_type == ct & dat2$leadSNP == as.character(leadSNP),]}
 
  colnames(res)[1] = "CHR"
	res$P = as.numeric(as.character(res$P))
	res$BETA = as.numeric(as.character(res$BETA))
	res$SE = as.numeric(as.character(res$SE))
  res = res[is.na(res$BETA) == F,]

  window=200000
  res$start_position = res$leadPOS - window
  res$end_position = res$leadPOS + window
  res = res[res$POS >= res$start_position & res$POS <= res$end_position,]
   
  res$leadSNP = NULL
  res$leadPOS = NULL
  res$start_position = NULL
  res$end_position = NULL
  res_coloc = res
  res_coloc = unique(res_coloc)


  if(nrow(res_coloc) > 1){

  
	ms_gwas_coloc = ms_gwas %>% filter(chr_pos %in% res_coloc$CHR_POS)
	ms_gwas_coloc$CHR_POS = ms_gwas_coloc$chr_pos
  ms_gwas_coloc = ms_gwas_coloc %>% distinct(CHR_POS,.keep_all=T)

  # Not needed - does not make a difference
  res_coloc$comb = paste0(res_coloc$CHR_POS,"_",res_coloc$REF,"_",res_coloc$ALT)
  ms_gwas_coloc$comb1 = paste0(ms_gwas_coloc$CHR_POS,"_",ms_gwas_coloc$A1,"_",ms_gwas_coloc$A2)
  ms_gwas_coloc$comb2 = paste0(ms_gwas_coloc$CHR_POS,"_",ms_gwas_coloc$A2,"_",ms_gwas_coloc$A1)
  res_coloc = res_coloc[res_coloc$comb %in% ms_gwas_coloc$comb1 == T | res_coloc$comb %in% ms_gwas_coloc$comb2 == T, ]
  #res_coloc[res_coloc$comb %in% ms_gwas_coloc$comb2 == T, "BETA"] = res_coloc[res_coloc$comb %in% ms_gwas_coloc$comb2 == T, "BETA"] * -1

  if (nrow(ms_gwas_coloc) > 1){

	# format
        ms_gwas_coloc = ms_gwas_coloc %>%
          mutate(beta = log(OR)) %>%
          mutate(z = (qnorm(1-P/2))) %>%
          mutate( se = abs(beta / z)) %>%
          mutate( se = ifelse(z == 0,0,se)) %>%
          mutate(varbeta = se ^ 2) #%>% 
          ms_gwas_coloc = ms_gwas_coloc %>%
          mutate(varbeta = ifelse(varbeta == 0,min(ms_gwas_coloc$varbeta[ms_gwas_coloc$varbeta!=0]),varbeta))
        
        coloc_ms_gene_input = list(
          beta = ms_gwas_coloc$beta,
          varbeta = ms_gwas_coloc$varbeta,
          snp = ms_gwas_coloc$CHR_POS,
          position = ms_gwas_coloc$BP,
          type="cc")
        
        coloc_eqtl_gene_input = list(
          beta = res_coloc$BETA,
          varbeta = res_coloc$SE^2,
          snp = res_coloc$CHR_POS,
          position = res_coloc$POS,
          type="quant",
          sdY = 1)
        
        coloc_abf = coloc.abf(coloc_ms_gene_input,coloc_eqtl_gene_input)
        print(coloc_abf)
        coloc_res = coloc_abf$summary %>% t %>% data.frame %>% mutate(gene = gene, cell = ct, source = source, cohort=cohort)
        
        if(coloc_abf$summary["PP.H4.abf"]>0.7){
          p=ggplot(coloc_abf$results,aes(position,SNP.PP.H4))+
            geom_point()+
            labs(y="Posterior probability of shared causal variant",x="Genomic position")+
            theme_minimal()+
            geom_text_repel(data = coloc_abf$results %>% filter(SNP.PP.H4>0.5),mapping = aes(position,SNP.PP.H4,label = snp))+
            ggtitle(paste0(coloc_res$gene,"\n",coloc_res$cell))
          plot_title = paste0(data.dir,"/figures/colocMS/coloc_",coloc_res$gene,"_",coloc_res$cell,"_",coloc_res$source,".png")
          png(plot_title,res=300,units="in",width=4,height=4)
          print(p)
          dev.off()
          
        }
       
        coloc_restable = coloc_abf$results
        bestSNP = coloc_restable[coloc_restable$SNP.PP.H4 == max(coloc_restable$SNP.PP.H4),"snp"][1]
        coloc_res$SNPhighestPP = res_coloc[res_coloc$CHR_POS == bestSNP,"SNP"][1]
        coloc_res

      }}
} 


all_tests = data.frame()
for (i in 1:nrow(all_genes)){
  print(i)
  test = coloc_function(i)
  all_tests = rbind(all_tests,test)
}

fwrite(all_tests,paste0(data.dir,"/Colocalization_results.tsv"))


# Relevant results
sig = all_tests[all_tests$PP.H4.abf > 0.7,]
all_genes$gene = all_genes$GENE; all_genes$cell = all_genes$cell_type
sig = join(sig,all_genes[,c("gene","SNP","cell","source")],by=c("gene","cell","source"),type="left")
sig$mineQTLSNP = sig$SNP; sig$SNP = NULL
sig$mineQTLp = NA


sig = data.frame(sig)

# Add some relevant information
for (i in 1:nrow(sig)){
  this_source = sig[i,"source"]
  this_ct = sig[i,"cell"]
  this_gene = sig[i,"gene"]
  this_SNP = sig[i,"mineQTLSNP"]
  this_cohort = sig[i,"cohort"]
  this_SNPPP = sig[i,"SNPhighestPP"]
  if (this_source == "CSF"){res = dat1[dat1$GENE == this_gene & dat1$cell_type == this_ct & dat1$leadSNP == as.character(this_SNP),]}
  if (this_source == "PBMC"){res = dat2[dat2$GENE == this_gene & dat2$cell_type == this_ct & dat2$leadSNP == as.character(this_SNP),]}

  res  = data.frame(res)
  sig[i,"mineQTLp"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNP),"P"][1]
  sig[i,"mineQTLp_FDR"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNP),"FDR"][1]
  sig[i,"mineQTLp_POS"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNP),"POS"][1]
  sig[i,"mineQTLp_CHR"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNP),"X.CHROM"][1]
  sig[i,"mineQTLp_chr_pos"] = paste0("chr",sig[i,"mineQTLp_CHR"],":",sig[i,"mineQTLp_POS"])
  sig[i,"mineQTLp_MSp"] = ms_gwas[ms_gwas$chr_pos == sig[i,"mineQTLp_chr_pos"],"P"][1]
  sig[i,"SNPhighestPP_eQTLp"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNPPP),"P"][1]
  sig[i,"SNPhighestPP_eQTLp_FDR"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNPPP),"FDR"][1]
  sig[i,"SNPhighestPP_POS"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNPPP),"POS"][1]
  sig[i,"SNPhighestPP_CHR"] = res[res$GENE == as.character(this_gene) & res$SNP == as.character(this_SNPPP),"X.CHROM"][1]
  sig[i,"SNPhighestPP_chr_pos"] = paste0("chr",sig[i,"SNPhighestPP_CHR"],":",sig[i,"SNPhighestPP_POS"])
  sig[i,"SNPhighestPP_MSp"] = ms_gwas[ms_gwas$chr_pos == sig[i,"SNPhighestPP_chr_pos"],"P"][1]
  
}
sig = sig[is.na(sig$gene) == F,]
sig = sig[is.na(sig$mineQTLp) == F,]



prevrep = fread(paste0(data.dir,"/comparison_with_references.txt"))
all_genes$Comb2 = paste0(all_genes$GENE,"_",all_genes$cell_type,"_",all_genes$source)
sig$Comb2 = paste0(sig$gene,"_",sig$cell,"_",sig$source)
sig = join(sig,all_genes[,c("Comb2","leadSNP")],by="Comb2",type="left")
sig$Comb2 = NULL; all_genes$Comb2 = NULL
sig$Comb2 = paste0(sig$leadSNP,"_",sig$gene,"_",sig$cell,"_",sig$source)
prevrep$Comb2 = paste0(prevrep$SNP,"_",prevrep$GENE,"_",prevrep$cell_type,"_",prevrep$Source)
sig = join(sig,prevrep[,c("siginGtex","sigeQTLGen","nomeQTLGen","sigYazar","Comb2")],by="Comb2",type="left")
sig[(sig$sigeQTLGen == T & is.na(sig$sigeQTLGen) == F) | (sig$siginGtex == 1 & is.na(sig$siginGtex) == F) | (sig$sigYazar == 1 & is.na(sig$sigYazar) == F),"prevrep" ] = 1
sig[(sig$nomeQTLGen == 1 & is.na(sig$nomeQTLGen) == F) ,"prevrepnom" ] = 1
prevrep$Comb2 = NULL

perm = fread(paste0(data.dir,"/All_fdr_eSNPs_permutation_results.txt"))
perm$Comb2 = paste0(perm$SNP,"_",perm$GENE,"_",perm$cell_type,"_",perm$source)
sig = join(sig,perm[,c("Comb2","permutation_p","perm_sig")],by="Comb2",type="left")


# Write significant results to file
fwrite(sig,paste0(data.dir,"/Colocalization_results_significant.tsv"))

# Create Table
out = sig
out = out[,c("gene","cell","source","PP.H4.abf","SNPhighestPP","SNPhighestPP_CHR","SNPhighestPP_POS","SNPhighestPP_eQTLp","SNPhighestPP_MSp","mineQTLSNP","mineQTLp","prevrep","permutation_p")]
out = unique(out)
out[is.na(out$permutation_p) == T,"permutation_p"] = "previously reported"
out$prevrep = NULL
colnames(out) = c("Gene","Cell_type","Source","PP_coloc","SNP_highestPP","Chr","Pos","eQTP_pvalue","MSGWAS_pvalue","eQTL_leadSNP","eQTL_leadSNP_pvalue","permutation_pvalue")
write.csv(out,paste0(data.dir,"/tables/supplementary_table_coloc_results.csv"))

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


# Read MS gwas data 
ms_gwas_plot = fread("/xx/shared_data/IMSGC_projects_summary_statistics/Science_2019/Discovery_phase_GWAS/meta-analysis/discovery_metav3.0.meta2")
ms_gwas_plot = ms_gwas_plot[is.na(ms_gwas_plot$P) == F,]
ms_gwas_plot = ms_gwas_plot %>% mutate(chr_pos = paste0("chr",CHR,":",BP))

ms_gwas_hg38_plot = read.table("/xx/shared_data/IMSGC_projects_summary_statistics/Science_2019/Discovery_phase_GWAS/meta-analysis/discovery_metav3.0.meta2_hg38.bed",sep="\t")
ms_gwas_hg38_plot$V4 = gsub("_",":",ms_gwas_hg38_plot$V4)

ms_gwas_hg38_plot = ms_gwas_plot %>%
  left_join(ms_gwas_hg38_plot %>% dplyr::rename(SNP = V4),by="SNP") %>%
  filter(!is.na(V3)) %>%
  filter(!is.na(P) & !is.na(OR)) %>%
  mutate(BP = V3) %>%
  dplyr::select(-V1,-V2,-V3,-V5) %>%
  mutate(SNP = ifelse(grepl("rs",SNP),SNP,paste0("chr",CHR,":",BP)))
ms_gwas_hg38_plot = ms_gwas_hg38_plot %>% mutate(chr_pos = paste0("chr",CHR,":",BP))
ms_gwas_plot = ms_gwas_hg38_plot
ms_gwas_plot$P = -log10(ms_gwas_plot$P)
ms_gwas_plot$POS = ms_gwas_plot$BP
ms_gwas_plot$POS = ms_gwas_plot$POS/1000000

celltypes=c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")


source(paste0(data.dir,"/Scripts/eqtl_plot_function_colocMSrisk.R"))

myplots = list()

for (i in 1:nrow(sig)){

  snp=as.character(sig[i,"SNPhighestPP"])
  gene=as.character(sig[i,"gene"])
  cohort=as.character(sig[i,"cohort"])
  
  ctout1=as.character(sig[i,"cell"])
  if (ctout1 == "CD4.T.cells"){ctout2 = "CD8.T.cell";ctout3 = "B.cells"}
  if (ctout1 != "CD4.T.cells"){
    if (ctout1 == "B.cells"){
      ctout2 = "CD4.T.cells";ctout3 = "CD8.T.cells"
    }else{
      ctout2 = "CD8.T.cells";ctout3 = "B.cells"
    }
    }
  source = as.character(sig[i,"source"])
  round=1
  plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,window=500000,cell_types = celltypes)

  myplots[[length(myplots)+1]] = plots

}


for (i in 1:nrow(sig)){


  snp=as.character(sig[i,"SNPhighestPP"])
  gene=as.character(sig[i,"gene"])
  ctout1=as.character(sig[i,"cell"])
  source = as.character(sig[i,"source"])
  cohort=as.character(sig[i,"cohort"])

  png(paste0(data.dir,"/figures/colocMS/",gene,"_",snp,"_",ctout1,"_",source,"_",cohort,".png"),width=8,height=8,units="in",res=300)
  grid.arrange(myplots[[i]][[1]],myplots[[i]][[4]],myplots[[i]][[14]],myplots[[i]][[13]],nrow=2)
  dev.off()

}


source(paste0(data.dir,"/Scripts/eqtl_plot_function_colocMSrisk.R"))

# Figure ZC2HC1A
i=11
snp=as.character(sig[i,"SNPhighestPP"])
gene=as.character(sig[i,"gene"])
cohort=as.character(sig[i,"cohort"])
ctout1=as.character(sig[i,"cell"])
ctout2 = "CD4.T.cells";ctout3 = "CD8.T.cells"
source = as.character(sig[i,"source"])
round=1
plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,window=500000,cell_types = celltypes)
zc2hc1a = plots
save(zc2hc1a,file=paste0(data.dir,"/figures/zc2hc1a.RData"))
png(paste0(data.dir,"/figures/coloc_ZC2HC1A.png"),width=8,height=8,units="in",res=300)
grid.arrange(plots[[1]],plots[[4]],plots[[15]],plots[[13]],nrow=2)
grid.text("A",x=unit(0.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(0.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("D",x=unit(4.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
dev.off()


source(paste0(data.dir,"/Scripts/eqtl_plot_function_colocMSrisk.R"))


# Figure AHI1
i=6
snp=as.character(sig[i,"SNPhighestPP"])
gene=as.character(sig[i,"gene"])
cohort=as.character(sig[i,"cohort"])
ctout1=as.character(sig[i,"cell"])
ctout2 = "CD4.T.cells";ctout3 = "CD8.T.cells"
source = as.character(sig[i,"source"])
round=1
plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,window=500000,cell_types = celltypes)
ahi1 = plots
save(ahi1,file=paste0(data.dir,"/figures/ahi1.RData"))
png(paste0(data.dir,"/figures/coloc_AHI1.png"),width=8,height=8,units="in",res=300)
grid.arrange(plots[[1]],plots[[4]],plots[[15]],plots[[13]],nrow=2)
grid.text("A",x=unit(0.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(0.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("D",x=unit(4.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
dev.off()

# Figure EAF2
i=7
snp=as.character(sig[i,"SNPhighestPP"])
gene=as.character(sig[i,"gene"])
cohort=as.character(sig[i,"cohort"])
ctout1=as.character(sig[i,"cell"])
ctout2 = "CD4.T.cells";ctout3 = "CD8.T.cells"
source = as.character(sig[i,"source"])
round=1
plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,window=500000,cell_types = celltypes)
eaf2 = plots
save(eaf2,file=paste0(data.dir,"/figures/eaf2.RData"))
png(paste0(data.dir,"/figures/coloc_EAF2.png"),width=8,height=8,units="in",res=300)
grid.arrange(plots[[1]],plots[[4]],plots[[15]],plots[[13]],nrow=2)
grid.text("A",x=unit(0.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(0.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("D",x=unit(4.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
dev.off()




source(paste0(data.dir,"/Scripts/eqtl_plot_function_colocMSrisk.R"))


# Figure ETS1
all_res = fread(paste0(data.dir,"/Colocalization_results.tsv"))
all_res = all_res[all_res$gene == "ETS1",]
snp=as.character(all_res[1,"SNPhighestPP"])
gene=as.character(all_res[1,"gene"])
cohort=as.character(all_res[1,"cohort"])
ctout1=as.character(all_res[1,"cell"])
ctout2 = "CD4.T.cells";ctout3 = "CD8.T.cells"
source = as.character(all_res[1,"source"])
round=1
plots = eQTL_plot(snp=snp,gene=gene,ctout1 =ctout1,ctout2 =ctout2,ctout3=ctout3,cohort=cohort,source=source,round=round,window=500000,cell_types = celltypes)
ets1 = plots
save(ets1,file=paste0(data.dir,"/figures/ets1.RData"))
png(paste0(data.dir,"/figures/coloc_ETS1.png"),width=8,height=8,units="in",res=300)
grid.arrange(plots[[1]],plots[[4]],plots[[15]],plots[[13]],nrow=2)
grid.text("A",x=unit(0.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.1,"in"),y=unit(7.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(0.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("D",x=unit(4.1,"in"),y=unit(3.9,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
dev.off()

####