library(tidyverse)


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

# recalculate P_adj with strict Bonferroni 
de_res = de_res %>%
  mutate(P_adj = p.adjust(PValue,method="bonf"))

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

ggplot(plot_dat,aes(cell_type,prop*100,fill=pheno))+
  geom_col(color="black",position=position_dodge())+
  theme_minimal()+
  scale_fill_brewer(palette="Set1")

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


make_plot("MS_CSF - OIND_CSF")
make_plot("MS_CSF - Control_CSF")

make_plot("MS_CSF - MS_PBMC")
make_plot("OIND_CSF - OIND_PBMC")
make_plot("Control_CSF - Control_PBMC")
make_plot("MS_PBMC - Control_PBMC")
make_plot("MS_PBMC - Control_PBMC")


# get gene sets for plotting 
library(msigdbr)
library(fgsea)

# retrieve gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
hallmark_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)


chol_genes = hallmark_list$HALLMARK_CHOLESTEROL_HOMEOSTASIS
mtor_genes = hallmark_list$HALLMARK_MTORC1_SIGNALING
oxphos_genes = hallmark_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION


# pathway = chol_genes 
pathway = oxphos_genes
comparison = "MS_CSF - MS_PBMC"

# plot
this_pathway = de_res %>%
  filter(gene %in% pathway & contrast == comparison) %>%
  mutate(
    direction = case_when(
      P_adj < 0.05 & logFC>0 ~ "Up",
      P_adj < 0.05 & logFC<0 ~ "Down",
      P_adj >=0.05 ~ "nonsig"
    ))
top_genes = this_pathway %>%
  arrange(desc(logFC)) %>%
  distinct(gene) %>%
  head(n=30)
this_pathway = this_pathway %>% filter(gene %in% top_genes$gene)

colours = c("Up" = "orange", "Down" = "purple", "nonsig" = "grey")

ggplot(this_pathway,aes(cell_type,gene,fill=direction))+
  geom_tile()+
  scale_fill_manual(values = colours)+
  theme_minimal()+
  theme(legend.position = "none")



# de plot
