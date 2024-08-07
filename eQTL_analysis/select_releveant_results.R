
library(readxl)
library(tidyr)
library(plyr)
library(dplyr)
library(data.table)

# set directory
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
setwd(data.dir)

# Read eQTL results (complete)
dat = readRDS(paste0(data.dir,"/","All_complete_results.rds"))

# Select relevant genes
genesimsgc = read.csv("References/Table_IMSGC.csv")[,2]
genesde = read.csv("References/POOLED_UNPAIRED_all_de_res.csv")
genesmigration = read.table("References/genes_TcellMigration.txt",h=T)
genesde = genesde[genesde$cell_type %in% c("B cells","CD4 T cells","CD8 T cells") == T,]
genesde = genesde[genesde$P_adj < 0.05,]
genesde = genesde[abs(genesde$logFC) > 0.5,]
csfgenes1 = unique(genesde$gene)
csfgenes = unique(c(csfgenes1,genesimsgc,genesmigration$GenID))
length(csfgenes)

reldat = dat[dat$cell_type %in% c("B.cells","CD4.T.cells","CD8.T.cells"),]
reldat = reldat[reldat$GENE %in% csfgenes == T,]

reldat$P = as.numeric(as.character(reldat$P))
reldat = data.frame(reldat)
reldat = reldat[is.na(reldat$P) == F,]
reldat$FDR = NA
reldat$FDR = p.adjust(reldat$P,method="fdr")

reldat$Comb = paste0(reldat$SNP,"_",reldat$GENE,"_",reldat$cell_type,"_",reldat$source)
reldat$P = as.numeric(as.character(reldat$P))
datfdr = reldat %>% filter(FDR< 0.1 & !is.na(FDR))  %>% group_by(cell_type,GENE,source) %>% arrange(P) %>%  filter(row_number()==1)


fwrite(datfdr,"All_fdr_eSNPs.txt")
saveRDS(reldat,paste0(data.dir,"/Complete_results.rds"))

# Get list of relevant genes
genes = unique(reldat$GENE)
write.table(genes,paste0(data.dir,"/selected_genes.txt"),sep="\t",row.names=F,quote=F,col.names=F)

# Save data for other cell types (all)
reldat$Comb2 = paste0(reldat$SNP,"_",reldat$GENE)
dat$Comb2 = paste0(dat$SNP,"_",dat$GENE)
outdat = dat[dat$Comb2 %in% reldat$Comb2 == T,]
saveRDS(outdat,paste0(data.dir,"/Complete_results_all_celltypes.rds"))

