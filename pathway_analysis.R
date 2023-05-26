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
library(decoupleR)

#######################################
# Read in data
#######################################

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/pathway_analysis/")


#######################################
# PROGENY
#######################################

# define testing comparisons
comparisons = c(
  "MS_PBMC - Control_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_CSF - Control_CSF",
  "MS_CSF - OIND_CSF",
  "MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
"OINDI_CSF - OINDI_PBMC",
"MS_CSF - OINDI_CSF",
"MS_PBMC - OINDI_PBMC")

# define cluster names
clusters = c("B cells","CD14 Mono","CD16 Mono","CD4 T cells","CD8 T cells","MAIT","mDCs","NK cells","pDCs","Plasma cells","Tregs")


# initialise results lists
res_overall = list()


# set permutation number
do_progeny = function(n_perm=10000){

  # grab progeny data
  net = get_progeny(organism = 'human', top = 100)

  for(x in 1:length(levels(factor(clusters)))){
    for(y in 1:length(comparisons)){
      # get name of DE file to read in
      comparison = comparisons[y]

      this_comparison = stringr::str_replace_all(comparison,"_"," ")
      this_comparison = stringr::str_replace_all(this_comparison,"-","vs")

      cluster = levels(factor(clusters))[x]
      message("Doing pathway analysis for cluster ",cluster)
      message("Comparison:", comparison)

      # read in DE file
      file = paste0("../de_results/edgeR_de_tests_",comparison,"_",cluster,".csv")
      if(!file.exists(file)){
        message("File does not exist. Moving on")
        next
      }
      de = read_csv(file)

      # format for decoupleR
      de_mat = matrix(de$logFC)
      rownames(de_mat) = de$gene

      # run decoupleR
      acts = run_wmean(mat=de_mat, net=net, .source='source', .target='target',
                        .mor='weight', times = n_perm, minsize = 10)

      # filter to just norm_wmean
      acts = acts %>% filter(statistic=="norm_wmean")

      # add in data
      acts$comparison = comparison
      acts$cluster = cluster

      # add to main results
      res_overall[[length(res_overall)+1]] <<- acts
    }
  }
}

do_progeny()
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"progeny_res_overall.csv")
res_overall$fdr = p.adjust(res_overall$p_value,method="fdr")

# plot
for(x in unique(res_overall$comparison)){
  top_paths = res_overall %>%
    filter(comparison == x) %>%
    arrange(fdr) %>% distinct(source,.keep_all=TRUE) %>%
    head(n=20)
  plot_data = res_overall %>%
    filter(source %in% top_paths$source) %>%
    filter(comparison == x) %>%
    arrange(score)
  plot_data$pathway = factor(plot_data$source,levels=unique(plot_data$source))

  # fill in 0s
  plot_data = plot_data %>%
    mutate(combo = paste0(cluster,"_",pathway))
  missing_combinations = expand.grid(plot_data$cluster,plot_data$pathway) %>%
    distinct() %>%
    tibble() %>%
    mutate(combo = paste0(Var1,"_",Var2)) %>%
    filter(!combo %in% plot_data$combo) %>%
    dplyr::rename("cluster" = Var1,"pathway"=Var2) %>%
    mutate(score = NA)
  plot_data = plot_data %>% bind_rows(plot_data,missing_combinations)


  # set labels
  plot_data = plot_data %>%
    mutate(p_label = case_when(
      fdr >= 0.01 ~ "",
      fdr < 0.01 & fdr >= 0.001 ~"*",
      fdr < 0.001 & fdr >= 0.0001 ~ "**",
      fdr < 0.0001 ~"***"
    ))
  png(paste0("progeny_summary_",x,".png"),res=300,units="in",width=8,height=4)


  p=ggplot(plot_data,aes(cluster,pathway,fill=score,label=p_label))+
    geom_tile(col="black")+
    geom_text(size=4)+
    scale_fill_gradient(low="purple",high="orange")+
    theme_classic()+
    labs(x="Cell type",y="Pathway")+
    ggtitle(x)+
    theme(axis.text.x=element_text(angle=90,vjust=-0.05))
  print(p)
  dev.off()
}


#######################################
# DOROTHEA
#######################################

# initialise results lists
res_overall = list()

# set permutation number
do_dorothea = function(n_perm=10000){

  # grab progeny data
  net = get_dorothea(organism='human', levels=c('A', 'B', 'C'))

  for(x in 1:length(levels(factor(clusters)))){
    for(y in 1:length(comparisons)){
      # get name of DE file to read in
      comparison = comparisons[y]
      this_comparison = stringr::str_replace_all(comparison,"_"," ")
      this_comparison = stringr::str_replace_all(this_comparison,"-","vs")

      cluster = levels(factor(clusters))[x]
      message("Doing pathway analysis for cluster ",cluster)
      message("Comparison:", comparison)

      # read in DE file
      file = paste0("../de_results/edgeR_de_tests_",comparison,"_",cluster,".csv")
      if(!file.exists(file)){
        message("File does not exist. Moving on")
        next
      }
      de = read_csv(file)

      # format for decoupleR
      de_mat = matrix(de$logFC)
      rownames(de_mat) = de$gene

      # run decoupleR
      acts = run_wmean(mat=de_mat, net=net, .source='source', .target='target',
                        .mor='mor', times = n_perm, minsize = 10)


      # filter to just norm_wmean
      acts = acts %>% filter(statistic=="norm_wmean")

      # add in data
      acts$comparison = comparison
      acts$cluster = cluster

      # add to main results
      res_overall[[length(res_overall)+1]] <<- acts
    }
  }
}

do_dorothea()
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"dorothea_res_overall.csv")
res_overall$fdr = p.adjust(res_overall$p_value,method="fdr")

# plot
for(x in unique(res_overall$comparison)){
  top_paths = res_overall %>%
    filter(comparison == x) %>%
    arrange(fdr) %>% distinct(source,.keep_all=TRUE) %>%
    head(n=20)
  plot_data = res_overall %>%
    filter(source %in% top_paths$source) %>%
    filter(comparison == x) %>%
    arrange(score)
  plot_data$pathway = factor(plot_data$source,levels=unique(plot_data$source))

  # fill in 0s
  plot_data = plot_data %>%
    mutate(combo = paste0(cluster,"_",pathway))
  missing_combinations = expand.grid(plot_data$cluster,plot_data$pathway) %>%
    distinct() %>%
    tibble() %>%
    mutate(combo = paste0(Var1,"_",Var2)) %>%
    filter(!combo %in% plot_data$combo) %>%
    dplyr::rename("cluster" = Var1,"pathway"=Var2) %>%
    mutate(score = NA)
  plot_data = plot_data %>% bind_rows(plot_data,missing_combinations)


  # set labels
  plot_data = plot_data %>%
    mutate(p_label = case_when(
      fdr >= 0.01 ~ "",
      fdr < 0.01 & fdr >= 0.001 ~"*",
      fdr < 0.001 & fdr >= 0.0001 ~ "**",
      fdr < 0.0001 ~"***"
    ))
  png(paste0("dorothea_summary_",x,".png"),res=300,units="in",width=8,height=4)


  p=ggplot(plot_data,aes(cluster,pathway,fill=score,label=p_label))+
    geom_tile(col="black")+
    geom_text(size=4)+
    scale_fill_gradient(low="purple",high="orange")+
    theme_classic()+
    labs(x="Cell type",y="Pathway")+
    ggtitle(x)+
    theme(axis.text.x=element_text(angle=90,vjust=-0.05))
  print(p)
  dev.off()
}
