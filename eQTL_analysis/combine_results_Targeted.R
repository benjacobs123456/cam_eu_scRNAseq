library(data.table)
library(plyr)

args = commandArgs(TRUE)

data.dir=args[1]
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"

cell_types=c("CD4.T.cells","CD8.T.cells","Tregs", "B.cells","Plasma.cells","NK.cells","mDCs","pDCs","CD14.Mono","CD16.Mono","MAIT.cells")

window = 500000
geneLoc = read.table("/xx/xxx/SC/CAM_TUM/eqtl/gene_list.txt",h=T)
geneLoc$left = geneLoc$start_position - as.numeric(as.character(window))
geneLoc$right = geneLoc$end_position + as.numeric(as.character(window))
geneLoc$GENE = geneLoc$hgnc_symbol

dat_all = data.frame()
dat_all_all = data.frame()

cohort = args[2]

for (source in c("CSF","PBMC")){
		dat_all = data.frame()
		for (i in 1:length(cell_types)){
			ct = cell_types[i]
			print(ct)
			file = paste0(data.dir,"/",source,"/",cohort,"/",ct,"_results1.tsv")
			if (file.exists(file)){
  		    dat = fread(file)
  		    if (nrow(dat) > 0){
			dat$POS = as.numeric(as.character(dat$POS))
			dat = join(dat,geneLoc[,c("GENE","left","right")],by="GENE",type="left")
			dat = dat[is.na(dat$POS) == F & dat$POS >= dat$left & dat$POS <= dat$right,]
			dat$cell_type = ct
			dat$source=source
			dat$cohort=cohort
			dat_all = rbind(dat_all,dat)}
			}

		}
 
 		dat_all = dat_all[is.na(dat_all$P) == F,]
		dat_all_all = rbind(dat_all_all,dat_all)
		dat_all$FDR = p.adjust(dat_all$P, method="fdr")
        saveRDS(dat_all,paste0(data.dir,"/",cohort,"_",source,"_complete_results.rds"))
}


dat_all_all$FDR = p.adjust(dat_all_all$P, method="fdr")
saveRDS(dat_all_all,paste0(data.dir,"/",cohort,"_complete_results.rds"))

