#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(MASS)
library(SingleCellExperiment)
library(Matrix.utils)
library(reshape2)

#######################################
# Read in data
#######################################


# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/gsea/")

#######################################
# GSEA
#######################################

library(msigdbr)
library(fgsea)

# retrieve gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
kegg_gene_sets = msigdbr(species = "Homo sapiens", category = "C2") %>% filter(gs_subcat=="CP:KEGG")
reactome_gene_sets = msigdbr(species = "Homo sapiens", category = "C2") %>% filter(gs_subcat=="CP:REACTOME")
go_gene_sets = msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_subcat=="GO:BP")

hallmark_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
kegg_list = split(x = kegg_gene_sets$gene_symbol, f = kegg_gene_sets$gs_name)
reactome_list = split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
go_list = split(x = go_gene_sets$gene_symbol, f = go_gene_sets$gs_name)

# set seed
set.seed(1)


comparisons = c(
  "MS_CSF - Control_CSF",
  "MS_PBMC - Control_PBMC",
  "MS_CSF - MS_PBMC",
  "Control_CSF - Control_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_PBMC - OIND_PBMC")

res_overall = list()
plots = list()

clusters = c("B cells","CD14 Mono","CD16 Mono","CD4 T cells","CD8 T cells","MAIT","mDCs","NK cells","pDCs","Plasma cells","Tregs")


nperm = 1000
do_gsea = function(geneset = "hallmark"){
  genelist = eval(parse(text = paste0(geneset,"_list")))
  for(x in 1:length(levels(factor(clusters)))){
    for(y in 1:length(comparisons)){

      comparison = comparisons[y]
      this_comparison = if(comparison == "MS_CSF - MS_PBMC"){
        "MS CSF vs MS PBMC"
      } else if(comparison=="Control_CSF - Control_PBMC"){
        "Control CSF vs Control PBMC"
      } else if(comparison=="MS_PBMC - Control_PBMC"){
        "MS PBMC vs Control PBMC"
      } else if(comparison=="MS_CSF - Control_CSF"){
        "MS CSF vs Control CSF"
      } else if(comparison=="MS_CSF - OIND_CSF"){
        "MS CSF vs OIND CSF"
      } else if(comparison=="MS_PBMC - OIND_PBMC"){
        "MS PBMC vs OIND PBMC"
      }


      cluster = levels(factor(clusters))[x]
      message("Doing GSEA for cluster ",cluster)
      message("Comparison:", comparison)

      # read in and rank
      file = paste0("../de/de_results/edgeR_de_tests_",comparison,"_",cluster,".csv")
      if(!file.exists(file)){
        message("File does not exist. Moving on")
        next
      }
      de = read_csv(file)
      ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
      ranked_genes_vector = ranked_genes$logFC
      names(ranked_genes_vector) = ranked_genes$gene
      message("There are ",nrow(ranked_genes)," genes in this analysis")

      # do gsea
      res = fgsea(genelist, stats = ranked_genes_vector, minSize=10,eps=0, nPermSimple = nperm)
      sig_pathways = res %>% filter(padj<0.05)

      if(nrow(sig_pathways)>0){
        leading_edges = list()
        for(i in c(1:nrow(sig_pathways))){
          leading_edge = sig_pathways$leadingEdge[[i]] %>% paste0(collapse=", ")
          leading_edges[[i]] = leading_edge
        }
        leading_edges = unlist(leading_edges)
        sig_pathways$leadingEdge = leading_edges

        write_csv(sig_pathways,file=paste0("../gsea/leading_edge_plots/sig_pathways_",geneset,"_",cluster,"_",comparison,".csv"))

        # plot leading edge enrichment
        png(paste0("../gsea/leading_edge_plots/gsea_enrichment_",geneset,"_",cluster,"_",comparison,".png"),res=300,units="in",height=8,width=8)
        print(plotGseaTable(genelist[sig_pathways$pathway], ranked_genes_vector, sig_pathways %>% arrange(desc(abs(NES))) %>% head(n=10) ))
        dev.off()
      }

      res = res %>% arrange(padj) %>%  dplyr::select(1,2,3,5,6) %>%  mutate(cell_type = cluster) %>%  mutate(comparison = this_comparison) %>% data.frame()
      res_overall[[(length(res_overall)+1)]] <<- res
      res = res %>% arrange(NES)
      res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
      topres = res %>% arrange(padj) %>% head(n=20) %>% arrange(desc(NES))
      topres$direction = ifelse(topres$NES>0,"up","down")
      colours = c("up" = "red","down"="blue")
      plot = ggplot(topres,aes(NES,pathway,label=ifelse(padj<0.05,"*"," "),alpha=-log10(pval),fill=direction))+
        geom_col()+
        theme_bw()+
        ggtitle(paste0("Comparison of ",cluster," in ",this_comparison)) +
        labs(x="Normalised enrichment score",fill=paste0("Up or down-regulated in \n ",this_comparison))+
        geom_text()+
        scale_fill_manual(values = colours, labels = c(paste0("Up in ",this_comparison),paste0("Down in ",this_comparison)))+
        scale_alpha(guide="none")

      png(paste0("../gsea/gsea_results_",geneset,"_",cluster,"_",comparison,".png"),res=300,units="in",height=8,width=8)
      print(plot)
      dev.off()
    }
  }
}

do_gsea("hallmark")
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"../gsea/res_overall.csv")
message("Finished GSEA")
res_overall$fdr = p.adjust(res_overall$pval,method="fdr")


for(x in c( "MS CSF vs MS PBMC","MS CSF vs Control CSF","Control CSF vs Control PBMC","MS PBMC vs Control PBMC", "MS CSF vs OIND CSF", "MS PBMC vs OIND PBMC" )){
  top_paths = res_overall %>% filter(comparison == x) %>% arrange(fdr) %>% distinct(pathway,.keep_all=TRUE) %>% head(n=10)
  plot_data = res_overall %>% filter(pathway %in% top_paths$pathway) %>% filter(comparison == x) %>% filter(grepl("HALLMARK",pathway)) %>% arrange(NES)
  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  plot_data$pathway = str_remove(plot_data$pathway,"HALLMARK_")
  png(paste0("../gsea/summary_",x,".png"),res=300,units="in",width=8,height=4)
  p=ggplot(plot_data,aes(cell_type,pathway,fill=NES,label=ifelse(fdr<0.01,"*","")))+
    geom_tile(col="black")+
    geom_text(size=4)+
    scale_fill_gradient(low="purple",high="orange")+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)+
    theme(axis.text.x=element_text(angle=90,vjust=-0.05))
  print(p)
  dev.off()
}

res_overall = list()
do_gsea("go")
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"../gsea/go_res_overall.csv")

res_overall = list()
do_gsea("reactome")
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"../gsea/reactome_res_overall.csv")

res_overall = list()
do_gsea("kegg")
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"../gsea/kegg_res_overall.csv")
