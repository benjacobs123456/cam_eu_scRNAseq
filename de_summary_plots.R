library(tidyverse)


#########################
# DA SUMMARY PLOTS
########################
# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/da_plots/")

# find files
files = list.files(full.names = T,pattern = "crude_labels_edgeR_da_tests")

# read in
da_res = purrr::map(files,function(x){
  read_csv(x) %>%
    mutate(contrast = str_remove_all(str_remove_all(x,"./crude_labels_edgeR_da_tests_"),".csv"))

})

da_res = do.call("bind_rows",da_res)

# get FDR value
da_res$P_adj = p.adjust(da_res$PValue,method="fdr")

# plot
da_res = da_res %>% filter(contrast %in% c(
"MS_CSF - Control_CSF",
"MS_CSF - OIND_CSF",
"MS_CSF - OINDI_CSF",
"MS_PBMC - Control_PBMC",
"MS_PBMC - OIND_PBMC",
"MS_PBMC - OINDI_PBMC",
"MS_CSF - MS_PBMC",
"OIND_CSF - OIND_PBMC",
"OINDI_CSF - OINDI_PBMC",
"Control_CSF - Control_PBMC"))

da_res$contrast = factor(da_res$contrast, levels = c(
"MS_CSF - Control_CSF",
"MS_CSF - OIND_CSF",
"MS_CSF - OINDI_CSF",
"MS_PBMC - Control_PBMC",
"MS_PBMC - OIND_PBMC",
"MS_PBMC - OINDI_PBMC",
"MS_CSF - MS_PBMC",
"OIND_CSF - OIND_PBMC",
"OINDI_CSF - OINDI_PBMC",
"Control_CSF - Control_PBMC"),ordered=T)

da_res$cell = factor(da_res$cell, levels = c(
"B cells",
"Plasma cells",
"Macrophages",
"CD16 Mono",
"CD14 Mono",
"mDCs",
"HSPCs",
"pDCs",
"NK cells",
"CD8 T cells",
"CD4 T cells",
"MAIT cells",
"Tregs"),ordered=T)

# plot
p = ggplot(da_res %>%
  mutate(p_label = case_when(
    P_adj < 0.0005 ~ "***",
    P_adj < 0.005 ~ "**",
    P_adj < 0.05 ~ "*",
    P_adj >= 0.05 ~ ""
  )) %>%
  mutate(logFC = ifelse(P_adj > 0.05, NA, logFC)),
  aes(cell,contrast,fill=logFC,label = p_label))+
  geom_tile(color="black")+
  geom_text()+
  theme_minimal()+
  scale_fill_gradient2(low="purple",high="orange",midpoint = 0)+
  labs(x="Cell type",y = "Comparison",fill="Log fold\nchange")

png("da_summary_plot_all_cells.png",res=600,units="in",width=13,height=3)
p
dev.off()

# repeat for celltypist

# find files
files = list.files(full.names = T,pattern = "celltypist_edgeR_da_tests")

# read in
da_res = purrr::map(files,function(x){
  read_csv(x) %>%
    mutate(contrast = str_remove_all(str_remove_all(x,"./celltypist_edgeR_da_tests_"),".csv"))

})

da_res = do.call("bind_rows",da_res)

# get FDR
da_res$P_adj = p.adjust(da_res$PValue,method="fdr")

# plot
da_res = da_res %>% filter(contrast %in% c(
"MS_CSF - Control_CSF",
"MS_CSF - OIND_CSF",
"MS_CSF - OINDI_CSF",
"MS_PBMC - Control_PBMC",
"MS_PBMC - OIND_PBMC",
"MS_PBMC - OINDI_PBMC",
"MS_CSF - MS_PBMC",
"OIND_CSF - OIND_PBMC",
"OINDI_CSF - OINDI_PBMC",
"Control_CSF - Control_PBMC"))

da_res$contrast = factor(da_res$contrast, levels = c(
"MS_CSF - Control_CSF",
"MS_CSF - OIND_CSF",
"MS_CSF - OINDI_CSF",
"MS_PBMC - Control_PBMC",
"MS_PBMC - OIND_PBMC",
"MS_PBMC - OINDI_PBMC",
"MS_CSF - MS_PBMC",
"OIND_CSF - OIND_PBMC",
"OINDI_CSF - OINDI_PBMC",
"Control_CSF - Control_PBMC"),ordered=T)

p = ggplot(da_res %>%
  mutate(p_label = case_when(
    P_adj < 0.0005 ~ "***",
    P_adj < 0.005 ~ "**",
    P_adj < 0.05 ~ "*",
    P_adj >= 0.05 ~ ""
  )) %>%
  mutate(logFC = ifelse(P_adj > 0.05, NA, logFC)),
  aes(cell,contrast,fill=logFC,label = p_label))+
  geom_tile(color="black")+
  geom_text()+
  theme_minimal()+
  scale_fill_gradient2(low="purple",high="orange",midpoint = 0)+
  labs(x="Cell type",y = "Comparison",fill="Log fold\nchange")+
  theme(axis.text.x = element_text(angle=90))


png("da_summary_plot_all_cells_celltypist.png",res=600,units="in",width=13,height=6)
p
dev.off()

##############################
# DE SUMMARY PLOTS
##############################

# set WD
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/de_plots")

# list files
files = list.files("../de_results/",pattern="edgeR_de",full.names=T)

# read in files
de_res = purrr::map(files,function(x){
  y = str_remove(str_remove(x,"../de_results//edgeR_de_tests_"),".csv")
  cell = str_split(y,"_")[[1]][4]
  y = str_remove(y,paste0("_",cell))
  read_csv(x) %>%
    mutate(contrast = y, cell_type = cell)
})

# combine
de_res = do.call("bind_rows",de_res)

# recalculate P_adj with global Bonf
de_res = de_res %>%
  mutate(P_adj = p.adjust(PValue,method="bonf"))


sig_effects0 =  de_res %>% filter(P_adj < 0.05) %>%
  mutate(gene_cell = paste0(gene,"_",cell_type))


# summary numbers
plot_dat = de_res %>%
  group_by(contrast,cell_type) %>%
  mutate(sig = ifelse(P_adj < 0.05, "Y","N")) %>%
  dplyr::count(sig) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(total = sum(n)) %>%
  filter(sig=="Y") %>%
  mutate(orig_contrast = contrast) %>%
  tidyr::separate(contrast,sep=" - ",into=c("pheno1","pheno2")) %>%
  filter(grepl("CSF",pheno1) & grepl("PBMC",pheno2)) %>%
  tidyr::separate(pheno1,"_",into=c("pheno","other"))

plot_dat %>% filter(pheno=="MS") %>% arrange(prop)
plot_dat %>% filter(pheno=="MS") %>% arrange(prop) %>% ungroup %>%summarise(median(prop))


# define function to make plots
make_plot = function(comparison = "MS_CSF - OIND_CSF"){

  # initialise plot list
  plots = list()
  res_overall = list()
  # loop through cell types
  for(this_cell_type in unique(de_res$cell_type)){


    # plot
    de_res = de_res %>%
      mutate(
        direction = case_when(
          P_adj < 0.05 & logFC>0 ~ "Up",
          P_adj < 0.05 & logFC<0 ~ "Down",
          P_adj >=0.05 ~ "nonsig"
        ))
    # de plot
    colours = c("Up" = "red", "Down" = "blue", "nonsig" = "grey")
    plot_dat = de_res %>% filter(contrast == comparison & cell_type == this_cell_type)
    p=ggplot(plot_dat,aes(logFC,-log10(PValue),color=direction,label=gene))+
      theme_bw()+
      geom_point()+
      ggrepel::geom_text_repel(data = plot_dat %>% filter(P_adj<0.05) %>% arrange(PValue) %>% filter(logFC>1) %>% head(5))+
      scale_color_manual(values = colours)+
      ggtitle(str_remove_all(this_cell_type,"_"))+
      theme(legend.position = "none")

    plots[[length(plots)+1]] = p
    res_overall[[length(res_overall)+1]] = plot_dat
  }


  png(paste0("de_",comparison,".png"),res=600,units="in",height=10,width=10)
  print(gridExtra::grid.arrange(grobs = plots, top = comparison))
  dev.off()
  res_overall = do.call("bind_rows",res_overall)
  write_csv(res_overall,paste0("de_",comparison,".csv"))
}

# save overall results
write_csv(de_res,"all_de_res.csv")
make_plot("MS_CSF - OIND_CSF")
make_plot("MS_CSF - Control_CSF")
make_plot("MS_CSF - OINDI_CSF")
make_plot("MS_CSF - MS_PBMC")
make_plot("OIND_CSF - OIND_PBMC")
make_plot("OINDI_CSF - OINDI_PBMC")
make_plot("Control_CSF - Control_PBMC")
make_plot("MS_PBMC - OIND_PBMC")
make_plot("MS_PBMC - OINDI_PBMC")

de_res = de_res %>% filter(cell_type != "CD16 Mono")
make_plot("MS_PBMC - Control_PBMC")




de_res %>%
  filter(contrast %in% c(
  "MS_CSF - OINDI_CSF",
  "MS_CSF - OIND_CSF"
  )) %>%
  filter(P_adj)



  ##############################
  # POOLED DE SUMMARY PLOTS
  ##############################

  # set WD
  setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/de_plots")

  # list files
  files = list.files("../de_results/",pattern="POOLED_UNPAIRED_DE_edgeR_de_tests_",full.names=T)

  # read in files
  de_res = purrr::map(files,function(x){
    y = str_remove(str_remove(x,"../de_results//POOLED_UNPAIRED_DE_edgeR_de_tests_"),".csv")
    read_csv(x) %>%
      mutate(cell_type = y)
  })

  # combine
  de_res = do.call("bind_rows",de_res)

  # recalculate P_adj with global Bonf
  de_res = de_res %>%
    mutate(P_adj = p.adjust(PValue,method="bonf"))


  # summary numbers
  plot_dat = de_res %>%
    group_by(cell_type) %>%
    mutate(sig = ifelse(P_adj < 0.05, "Y","N")) %>%
    dplyr::count(sig) %>%
    mutate(prop = n/sum(n)) %>%
    mutate(total = sum(n)) %>%
    filter(sig=="Y")

# get cell counts
cell_counts = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/da_plots/all_cell_counts.csv") %>%
  pivot_wider(names_from = source,values_from = n, id_cols = cell_type_crude) %>%
  dplyr::rename("cell_type" = cell_type_crude)

plot_dat = plot_dat %>%
left_join(cell_counts,by="cell_type")

# set up palette
n_col = plot_dat$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)
plot_dat$cell_type = factor(plot_dat$cell_type,
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
ggplot(plot_dat,aes(log10(CSF),n,col=cell_type,label=cell_type))+
  geom_point(size=3)+
  theme_minimal()+
  scale_color_manual(values = colour_pal)+
  labs(x="Log10(N total cells in CSF)",y="N significant DE genes")+
  geom_label_repel(fill="black",show.legend=F)+
  theme(legend.position="none")

ggplot(plot_dat,aes(log10(PBMC),n,col=cell_type,label=cell_type))+
  geom_point(size=3)+
  theme_minimal()+
  scale_color_manual(values = colour_pal)+
  labs(x="Log10(N total cells in PBMC)",y="N significant DE genes")+
  geom_label_repel(fill="black",show.legend=F)+
  theme(legend.position="none")


de_res$contrast = "CSF - PBMC"
comparison = "CSF - PBMC"

# save overall results
write_csv(de_res,"POOLED_UNPAIRED_all_de_res.csv")
make_plot("CSF - PBMC")


# correlation
pooled = read_csv("POOLED_UNPAIRED_all_de_res.csv")
unpooled = read_csv("all_de_res.csv") %>% filter(contrast == "MS_CSF - MS_PBMC")

cor.test(combo$logFC.x,combo$logFC.y)
combo = pooled %>%
left_join(unpooled,by=c("cell_type","gene"))
ggplot(combo,aes(logFC.x,logFC.y,col=cell_type))+geom_point()+
labs(x="logFC (All CSF - All PBMC)",
y="logFC (MS CSF - MS PBMC)")+
theme_minimal()+
geom_abline(intercept=0, slope =1, alpha=0.5)




# define function to make plots
make_plot_pooled_csf_v_pbmc = function(){

  # initialise plot list
  plots = list()
  res_overall = list()
  # loop through cell types
  for(this_cell_type in unique(de_res$cell_type)){


    # plot
    de_res = de_res %>%
      mutate(
        direction = case_when(
          P_adj < 0.05 & logFC>0 ~ "Up",
          P_adj < 0.05 & logFC<0 ~ "Down",
          P_adj >=0.05 ~ "nonsig"
        ))
    # de plot
    colours = c("Up" = "red", "Down" = "blue", "nonsig" = "grey")
    plot_dat = de_res %>% filter(cell_type == this_cell_type)
    p=ggplot(plot_dat,aes(logFC,-log10(PValue),color=direction,label=gene))+
      theme_bw()+
      geom_point(size=1,alpha=0.8)+
      ggrepel::geom_text_repel(data = plot_dat %>% filter(P_adj<0.05) %>% arrange(PValue) %>% filter(logFC>1) %>% head(10),
        size=3,
        max.overlaps=Inf,
        min.segment.length=0,
        force=10,
        force_pull=0,
        max.iter=1e6,
        segment.color="black")+
      scale_color_manual(values = colours)+
      ggtitle(str_remove_all(this_cell_type,"_"))+
      theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

    plots[[length(plots)+1]] = p
    res_overall[[length(res_overall)+1]] = plot_dat
  }


  png(paste0("POOLED_de_",comparison,".png"),res=600,units="in",height=8,width=8)
  print(gridExtra::grid.arrange(grobs = plots, top = "CSF vs PBMC (pooled)"))
  dev.off()
  res_overall = do.call("bind_rows",res_overall)
}


  # filter to cell types of interest
de_res = de_res %>% filter(cell_type %in%
c("Plasma cells","CD8 T cells","CD4 T cells","B cells"))
  sig_effects =  de_res %>% filter(P_adj < 0.05) %>%
    mutate(gene_cell = paste0(gene,"_",cell_type))

make_plot_pooled_csf_v_pbmc()
