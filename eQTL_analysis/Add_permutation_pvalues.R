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

# Get eGenes and eSTNPs
eSNPs = fread(paste0(data.dir,"/All_fdr_eSNPs.txt"))
eSNPs = unique(eSNPs)
eSNPsall = eSNPs
eSNPs$Comb = paste0(eSNPs$source,"__",eSNPs$cohort,"__",eSNPs$SNP,"__",eSNPs$GENE)


# Only do this for "novel" hits
novel = fread(paste0(data.dir,"/comparison_with_references.txt"))
novel = data.frame(novel)
novel$Comb = paste0(novel$Source,"__",novel$Cohort,"__",novel$SNP,"__",novel$GENE)
novel = novel[novel$Comb %in% eSNPs$Comb == T,]
novel = novel[(novel$sigeQTLGen != T | is.na(novel$sigeQTLGen) == T) & (novel$siginGtex == 0 | is.na(novel$siginGtex) == T) & (novel$sigYazar == 0 | is.na(novel$sigYazar) == T),]
novel = novel[(novel$nomeQTLGen != T | is.na(novel$nomeQTLGen) == T),]

eSNPs = data.frame(eSNPs)
eSNPs$rowname = NA
eSNPs[eSNPs$Comb %in% novel$Comb == T,"rowname"] = 1:nrow(eSNPs[eSNPs$Comb %in% novel$Comb == T,])

# Read permutation results
res = read.table(paste0(data.dir,"/permutation_results.txt"),sep="\t")
colnames(res) = c("rowname","perm_p")
eSNPs = join(eSNPs,res,by="rowname",type="left")
eSNPs$permutation_p = as.numeric(as.character(eSNPs$perm_p))

# Select results with significant perm_p
sw_fdr = max(eSNPs$P)
eSNPs$perm_sig = NA
eSNPs[is.na(eSNPs$perm_p) == F,"perm_sig"] = 0
eSNPs[eSNPs$permutation_p < sw_fdr & is.na(eSNPs$permutation_p) == F,"perm_sig"] = 1
eSNPs[eSNPs$perm_p %in% c("<1e-05","<1e-06","<1e-07","<1e-08","<1e-09") & is.na(eSNPs$perm_p) == F,"perm_sig"] = 1
fwrite(eSNPs,paste0(data.dir,"/All_fdr_eSNPs_permutation_results.txt"))

# Prepare figures
results = fread(paste0(data.dir,"/comparison_with_references.txt"))
eSNPsall$Comb2 = paste0(eSNPsall$SNP,"_",eSNPsall$GENE)
results$Comb2 = paste0(results$SNP,"_",results$GENE)
comp = join(eSNPsall,results[,c("Comb2","siginGtex","sigeQTLGen","sigYazar","nomeQTLGen","window_Yazar","window_Gtex","window_eQTLGen")],by="Comb2",type="left")
eSNPs$Comb2 = paste0(eSNPs$SNP,"_",eSNPs$GENE)
comp = join(comp,eSNPs[,c("Comb2","perm_p","permutation_p","perm_sig")],by="Comb2",type="left")
comp = unique(comp)
comp$prev = "novel"
comp[(comp$nomeQTLGen == 1 & is.na(comp$nomeQTLGen) == F) ,"prev"] = "previously reported at \nnominal significance"
comp[(comp$sigeQTLGen == T & is.na(comp$sigeQTLGen) == F) | (comp$siginGtex == 1 & is.na(comp$siginGtex) == F) | (comp$sigYazar == 1 & is.na(comp$sigYazar) == F),"prev"] = "previously \nreported"
plotdat = data.frame(table(comp[,c("source","prev")]))
plotdat$prev = factor(plotdat$prev,ordered=T,levels=c("novel","previously reported at \nnominal significance","previously \nreported"))

p=ggplot(plotdat,aes(x=source,y=Freq,fill=prev))+
  geom_col() + theme_minimal() + 
  scale_fill_manual(values=c("goldenrod","grey70","grey40")) + 
  ylab("Number of \nassociation signals") + xlab("") + theme(legend.title = element_blank(), legend.position="bottom",legend.margin=margin(0,0,0,0), legend.box.margin=margin(-15,-10,-5,-10))+
   theme(plot.title = element_text(family="sans",size=9,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9,angle=90), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))+
  theme(strip.text.x = element_text(family="sans",size=9,face = "bold"),strip.text.y = element_text(family="sans",size=9,face = "bold"))

plotdat2 = data.frame(table(comp[comp$prev == "novel",c("source","perm_sig")]))
plotdat2[plotdat2$perm_sig == 1 & is.na(plotdat2$perm_sig) == F,"Perm_Sig"] = "Yes"
plotdat2[plotdat2$perm_sig == 0,"Perm_Sig"] = "No"
plotdat2$Perm_Sig = factor(plotdat2$Perm_Sig,ordered=T,levels=c("Yes","No"))
p2=ggplot(plotdat2,aes(x=source,y=Freq,fill=Perm_Sig))+
  geom_col() + theme_minimal() + 
  #scale_fill_discrete(name="significant permutation p-value")+
  scale_fill_manual(values=c("firebrick3","darkgrey"),name="significant permutation p-value") +  theme(legend.position="bottom",legend.margin=margin(0,0,0,0), legend.box.margin=margin(-15,-10,-5,-10)) + 
  ylab("Number of 'novel' \nassociation signals") + xlab("") +
   theme(plot.title = element_text(family="sans",size=9,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9,angle=90), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))+
  theme(strip.text.x = element_text(family="sans",size=9,face = "bold"),strip.text.y = element_text(family="sans",size=9,face = "bold"))

f1_3 = p
f1_4 = p2
save(f1_3,file=paste0(data.dir,"/figures/Fig1_3.RData"))
save(f1_4,file=paste0(data.dir,"/figures/Fig1_4.RData"))

# Closer look at significant "novel" results
sig = comp[comp$perm_sig == 1,]
nrow(sig[sig$window_Gtex == 0 & sig$window_Yazar == 0 & sig$window_eQTLGen == 0,])

# Prepare supplementary table
out = comp[comp$prev == "novel"]
out = out[,c("X.CHROM","POS","SNP","GENE","source","cell_type","BETA","SE","P","FDR","perm_p")]
colnames(out) = c("Chr","Pos","SNP","Gene","Source","Cell_type","Beta","SE","p","FDR p","permutation p")
write.csv(out,paste0(data.dir,"/tables/supplementary_table_novel_eqtls.csv"))
