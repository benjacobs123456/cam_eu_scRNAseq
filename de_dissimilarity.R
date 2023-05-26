#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(gridExtra)
library(MASS)
library(reshape2)

#######################################
# Read in data
#######################################


# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/de_plots/")

# read in files
de_files = list.files(path="../de_results/",pattern="edgeR_de_tests_",full.names=T)
de_res = purrr::map(de_files, function(x){

short_name = str_remove(x,"../de_results//edgeR_de_tests_")
short_name = str_split(short_name," - ")[[1]]
short_name = str_split(short_name,"_")

read_csv(x) %>%
mutate(pheno1 = short_name[[1]][1],
pheno2 = short_name[[2]][1],
source1 = short_name[[1]][2],
source2 = short_name[[2]][2],
cell_type = str_remove(short_name[[2]][3],".csv"))
})

de_res = do.call("bind_rows",de_res)

# bonf
de_res = de_res %>%
  mutate(fdr = p.adjust(PValue,method="bonf"))

#######################################
# CSF v PBMC
#######################################


find_het_genes = function(pheno){

  # find MS significant tests vs control
  ms_sig = de_res %>%
    filter(fdr < 0.05 &
      pheno1=="MS" &
      pheno2 =="MS" &
      source1 == "CSF" &
      source2 == "PBMC")

  # find OIND vs control CSF
  csf_v_pbmc = de_res %>%
    filter(pheno1==pheno &
      pheno2 ==pheno &
      source1 == "CSF" &
      source2 == "PBMC")

  # get list of clusters
  clusters = unique(ms_sig$cell_type)

  plots = list()
  all_de_genes = list()
  # loop through each cluster
  for(cluster in clusters){

    # filter to this cluster only
    ms = ms_sig %>% filter(cell_type==cluster)
    control = csf_v_pbmc %>% filter(cell_type==cluster)

    # filter to intersection
    ms = ms %>% filter(gene %in% control$gene)
    control = control %>% filter(gene %in% ms$gene)

    # combine
    combo_dat = ms %>% left_join(control,by="gene")

    # find discordant genes
    combo_dat = combo_dat %>%
    mutate(discordant = ifelse(
      fdr.y < 0.05 & sign(logFC.x)!=sign(logFC.y),
      "Discordant","Concordant"))

    # add to overall
    all_de_genes[[length(all_de_genes)+1]] = combo_dat

    # plot
    rho = cor.test(combo_dat$logFC.x,combo_dat$logFC.y)$estimate
    n_genes = nrow(combo_dat)

    p = ggplot(combo_dat,aes(logFC.x,logFC.y,col=discordant))+
    geom_point()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    theme_minimal()+
    scale_color_manual(values = c("black","red"))+
    labs(x="LogFC (MS)",y=paste0("LogFC (",pheno,")"))+
    geom_label_repel(data = combo_dat %>% filter(discordant=="Discordant"),mapping = aes(logFC.x,logFC.y,label=gene),max.overlaps = 10, max.time=5,max.iter=1e6,force=1000)+
    ggtitle(paste0(cluster,"\nrho=",round(rho,2),"\nN genes=",n_genes))+
    theme(legend.position="none")

    # add to plot list
    plots[[length(plots)+1]] = p
  }

  # combine plots


  png(paste0("de_dissimilarity_plots_",pheno,".png"),res=600,units="in",width=10,height=10)
  do.call("grid.arrange",c(grobs = plots,top=pheno))
  dev.off()

  # joint plot
  all_de_genes = do.call("bind_rows",all_de_genes)

  n_col = length(unique(all_de_genes$cell_type.x))
  colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
  colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

  all_de_genes$cell_type.x = factor(
  all_de_genes$cell_type.x,
  ordered=TRUE,
  levels = c("B cells",
  "Plasma cells",
  "mDCs",
  "HSPCs",
  "CD14 Mono",
  "CD16 Mono",
  "Macrophages",
  "CD4 T cells",
  "CD8 T cells",
  "Tregs",
  "MAIT cells",
  "NK cells",
  "pDCs"
  ))


p = ggplot(all_de_genes,aes(logFC.x,logFC.y,col=cell_type.x))+
    geom_point()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    theme_minimal()+
    scale_color_brewer(palette="Paired")+
    labs(x="LogFC (MS)",y=paste0("LogFC (",pheno,")"),color="Cell type")+
    geom_abline(intercept=0,slope=1,linetype="dashed",alpha=0.5)
    png(paste0("de_dissimilarity_plots_all_together",pheno,".png"),res=600,units="in",width=5,height=4)
    print(p)
      dev.off()
  
}
find_het_genes("OIND")
find_het_genes("Control")
find_het_genes("OINDI")

