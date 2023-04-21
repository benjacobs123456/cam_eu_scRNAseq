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
library(reshape2)

#######################################
# Read in data
#######################################
# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/tcr/")

# Read in data
t_cells = readRDS("../t_cells_post_processing.rds")
rownames(t_cells@meta.data) = colnames(t_cells)


#######################################
# Big picture descriptive stats
#######################################
message("No. cells:")
nrow(t_cells@meta.data)

#######################################
# Recluster and reannotate
#######################################


t_subsets_blueprint = c(
"CD4+ T-cells",
"CD4+ Tcm",
"CD4+ Tem",
"CD8+ T-cells",
"CD8+ Tcm",
"CD8+ Tem",
"Tregs")

t_subsets_monaco = c("T regulatory cells",
"Th1/Th17 cells",
"Effector memory CD8 T cells",
"Th2 cells",
"Central memory CD8 T cells",
"Follicular helper T cells",
"Th17 cells",
"Terminal effector CD8 T cells",
"Th1 cells",
"Naive CD4 T cells",
"Naive CD8 T cells",
"Terminal effector CD4 T cells")


t_cell_markers = c("CD3G","IL7R","CD8A","TCF7","FOXP3","CCL5","CCR7","GZMH","CD27")


# filter by annotation
filtered_tcells = subset(t_cells, subset = ann_monaco %in% t_subsets_monaco)
discarded_cells = t_cells@meta.data %>% filter(!cell_id %in% filtered_tcells@meta.data$cell_id)
table(discarded_cells$ann_monaco)

# recluster
set.seed(1)
DefaultAssay(filtered_tcells)="SCT"
filtered_tcells = RunUMAP(filtered_tcells,reduction="harmony",dims=1:50)
filtered_tcells = FindNeighbors(filtered_tcells)
filtered_tcells = FindClusters(filtered_tcells,resolution=0.5)
filtered_tcells = SetIdent(filtered_tcells,value = filtered_tcells[['ann_monaco']])

p1=DimPlot(filtered_tcells)
p2=FeaturePlot(filtered_tcells,features=t_cell_markers)
p3=DotPlot(filtered_tcells,features=t_cell_markers)

png("dimplot.png",res=300,units="in",width=6,height=4)
p1
dev.off()

png("featureplot.png",res=300,units="in",width=8,height=8)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=10,height=4)
p3
dev.off()

# recode blueprint annotations for plotting (simplify)
non_t_cell_annotations = c("CLP","Erythrocytes","Megakaryocytes","Monocytes","NK cells")
filtered_tcells@meta.data = filtered_tcells@meta.data %>%
mutate(ann_blueprint_simple  = ifelse(ann_blueprint %in% non_t_cell_annotations,"Non-T cell",ann_blueprint)) %>%
mutate(ann_blueprint_simple  = ifelse(is.na(ann_blueprint_simple),"Non-T cell",ann_blueprint_simple))
filtered_tcells@meta.data$ann_blueprint_simple = as.factor(filtered_tcells@meta.data$ann_blueprint_simple)

p1=ggplot(filtered_tcells@meta.data,aes(ann_monaco,fill=ann_blueprint_simple))+geom_bar(position="fill")+labs(x="Monaco annotation",y="Proportion",fill="Blueprint/ENCODE annotation")+theme_bw()+scale_fill_brewer(palette="Set2")+theme(axis.text.x=element_text(angle=90))

png("monaco_v_blueprint_annotations.png",res=300,height=4,width=12,units="in")
p1
dev.off()

t_cells = filtered_tcells


# cluster biomarkers
biomarkers = FindAllMarkers(t_cells,recorrect_umi=FALSE,logfc.threshold=0.25,min.pct=0.25,only.pos=TRUE)
topbiomarkers = biomarkers %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=5)
write_csv(topbiomarkers,"biomarkers.csv")

#######################################
# DA & composition
#######################################

colour_pal <- RColorBrewer::brewer.pal(12, "Set2")
colour_pal <- grDevices::colorRampPalette(colour_pal)(12)

png("proportions.png",res=300,width=8,height=4,units="in")
ggplot(t_cells@meta.data,aes(phenotype,fill=ann_monaco))+
geom_bar(position="fill")+
facet_wrap(~source)+
scale_fill_manual(values=colour_pal)+
theme_bw()+
labs(fill="B cell subset",x="Phenotype",y="Proportion of T cell pool")
dev.off()

# DA testing
# create new unique ID with donor and source
t_cells@meta.data$cell_type = Idents(t_cells)
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
t_cells@meta.data = t_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# stash sample info
sample_info = t_cells@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(t_cells@meta.data$cell_type,t_cells@meta.data$donor_source)

# filter out clusters with 0 counts
abundances = abundances[rowSums(abundances)>0,]

sample_info = sample_info %>% filter(donor_source %in% colnames(abundances))
sample_info = sample_info[match(colnames(abundances),sample_info$donor_source),]
sample_info$grouping = paste0(sample_info$phenotype,"_",sample_info$source)
sample_info$grouping = sapply(sample_info$grouping, function(x){
    if(x == "Noninflammatory_CSF"){
      return("Control_CSF")
    } else if(x == "Noninflammatory_PBMC"){
      return("Control_PBMC")
    } else {
      return(x)
    }
    })

y.ab = DGEList(abundances, samples=sample_info)

# update sample info
y.ab$samples =  y.ab$samples %>% left_join(t_cells@meta.data %>% distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Gender,drb1_1501_dose,donor_source),by="donor_source")

# now loop
da_overall_list = list()
results_df = data.frame()
design <- model.matrix(~ 0 + grouping + Age + Gender,y.ab$sample)
colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Gender")

y.ab <- estimateDisp(y.ab, design,trend="none")
fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# define contrast for testing

contrast =  makeContrasts(
  MS_CSF - MS_PBMC,
  OIND_CSF - OIND_PBMC,
  Control_CSF - Control_PBMC,
  MS_CSF - OIND_CSF,
  MS_CSF - Control_CSF,
  OIND_CSF - Control_CSF,
  MS_PBMC - OIND_PBMC,
  MS_PBMC - Control_PBMC,
  OIND_PBMC - Control_PBMC,
  levels = design
)

# do da tests
for(i in 1:length(colnames(contrast))){
  contrast_name = colnames(contrast)[i]
  res = glmQLFTest(fit, contrast = contrast[,i])

  # write to file
  print(summary(decideTests(res)))
  res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
  res$cell = rownames(res)
  write_csv(res,paste0("edgeR_da_tests_",contrast_name,".csv"))
  res$significant = ifelse(res$P_adj<0.01,"yes","no")
  res$direction = ifelse(res$logFC>0,"Up","Down")
  comparison_label = contrast_name

  # de plot
  colours = c("Up" = "red","Down" = "blue")
  plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,label=cell))+
  theme_bw()+
  geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
  geom_vline(xintercept = 0,alpha=0.2)+
  geom_point()+
  NoLegend()+
  geom_text_repel(data=res %>% arrange(PValue) %>% head(n=30),mapping=aes(),max.overlaps=100)+
  ggtitle(comparison_label)+
  scale_color_manual(values = colours)

  png(paste0("da_plot_",contrast_name,"_.png"),res=300,units="in",height=4,width=4)
  print(plot)
  dev.off()
}

#######################################
# DE
#######################################


DefaultAssay(t_cells)="RNA"
t_cells@meta.data$cell_type = t_cells@meta.data$ann_monaco
t_cells@meta.data = t_cells@meta.data %>% mutate(cell_type = ifelse(cell_type=="Th1/Th17 cells","Th1_Th17 cells",cell_type))

# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
t_cells@meta.data = t_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# convert to sce object
t_cells.sce = as.SingleCellExperiment(t_cells)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = t_cells@meta.data %>%
  group_by(donor.id,source,phenotype,cell_type) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# aggregate counts
groups = colData(t_cells.sce)[, c("cell_type", "phenotype","source","donor.id")]
aggregated_counts  = aggregate.Matrix(t(counts(t_cells.sce)),
groupings = groups, fun = "sum")

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]


# get names for clusters
clusters = levels(factor(t_cells@meta.data$cell_type))

# split by cluster
filtered_matrices = lapply(clusters,function(x){
  filtered_matrix = aggregated_counts[grepl(x,rownames(aggregated_counts)),] %>% t()
})
names(filtered_matrices)=clusters

results_df = data.frame()
de_results = lapply(clusters,function(cell_type){

  de_input = filtered_matrices[[cell_type]]

  group_vector = lapply(colnames(de_input),function(y){
      if(grepl("MS_PBMC",y)){
        "MS_PBMC"
      } else if(grepl("MS_CSF",y)){
        "MS_CSF"
      } else if(grepl("Noninflammatory_CSF",y)){
        "Control_CSF"
      } else if(grepl("Noninflammatory_PBMC",y)){
        "Control_PBMC"
      } else if(grepl("OIND_CSF",y)){
        "OIND_CSF"
      } else if(grepl("OIND_PBMC",y)){
        "OIND_PBMC"
      }
    }) %>% unlist %>% factor()

    # make the DGE object
    y=DGEList(de_input,group=group_vector,remove.zeros=TRUE)

    # update sample info
    y$samples =  y$samples %>% mutate(full_cell_id = rownames(y$samples)) %>% left_join(t_cells@meta.data %>% dplyr::select(Age,Gender,drb1_1501_dose,full_cell_id) %>% distinct(full_cell_id,.keep_all=TRUE),by="full_cell_id")

    # make the design matrix
    message("Doing DE for ",cell_type)
    design = model.matrix(~0+group_vector+Age+Gender,y$samples)
    colnames(design) = c(levels(group_vector),"Age","Gender")


    keep = filterByExpr(
      y,
      design = design,
      group = group_vector,
      min.count = 10,
      min.total.count = 15,
      large.n = 10,
      min.prop = 0.7)
    y = y[keep, , keep.lib.sizes=FALSE]
    y = calcNormFactors(y)
    y = estimateDisp(y,design,robust=TRUE)

    if(is.na(y$common.dispersion)){
      message("Couldn't fit dispersion model. Skipping")
    } else {

    fit = glmQLFit(y, design, robust=TRUE)

    # define contrast for testing

    contrast =  if(length((unique(y$samples$group)))==6){
      makeContrasts(
        MS_CSF - MS_PBMC,
        OIND_CSF - OIND_PBMC,
        Control_CSF - Control_PBMC,
        MS_CSF - OIND_CSF,
        MS_CSF - Control_CSF,
        OIND_CSF - Control_CSF,
        MS_PBMC - OIND_PBMC,
        MS_PBMC - Control_PBMC,
        OIND_PBMC - Control_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==5 & (!("Control_CSF" %in% unique(y$samples$group)))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        OIND_CSF - OIND_PBMC,
        MS_CSF - OIND_CSF,
        MS_PBMC - OIND_PBMC,
        MS_PBMC - Control_PBMC,
        OIND_PBMC - Control_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==4 & !("Control_PBMC" %in% unique(y$samples$group)) & !("Control_CSF" %in% unique(y$samples$group))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        OIND_CSF - OIND_PBMC,
        MS_CSF - OIND_CSF,
        MS_PBMC - OIND_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==5 & (!("OIND_CSF" %in% unique(y$samples$group)))){
        makeContrasts(
        MS_CSF - MS_PBMC,
        Control_CSF - Control_PBMC,
        MS_CSF - Control_CSF,
        MS_PBMC - OIND_PBMC,
        MS_PBMC - Control_PBMC,
        OIND_PBMC - Control_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==4 & !("OIND_PBMC" %in% unique(y$samples$group)) & !("OIND_CSF" %in% unique(y$samples$group))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        Control_CSF - Control_PBMC,
        MS_CSF - Control_CSF,
        MS_PBMC - Control_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==4 & !("Control_CSF" %in% unique(y$samples$group)) & !("OIND_CSF" %in% unique(y$samples$group))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        MS_PBMC - OIND_PBMC,
        MS_PBMC - Control_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==3 & !("OIND_PBMC" %in% unique(y$samples$group)) & !("Control_PBMC" %in% unique(y$samples$group))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        MS_CSF - OIND_CSF,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==3 & !("OIND_CSF" %in% unique(y$samples$group)) & !("Control_PBMC" %in% unique(y$samples$group))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        MS_PBMC - OIND_PBMC,
        levels = design
      )
    } else if(length((unique(y$samples$group)))==2){
      makeContrasts(
        MS_CSF - MS_PBMC,
        levels = design
      )
    }

    # do de tests
    for(i in 1:length(colnames(contrast))){
    contrast_name = colnames(contrast)[i]
    res = glmQLFTest(fit, contrast = contrast[,i])

    # write to file
    print(summary(decideTests(res)))
    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$gene = rownames(res)
    write_csv(res,paste0("edgeR_de_tests_",contrast_name,"_",cell_type,".csv"))
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
    total_genes = nrow(res)
    non_sig = sum(res$significant=="no")
    sig_up = sum(res$significant=="yes" & res$direction=="Up")
    sig_down = sum(res$significant=="yes" & res$direction=="Down")
    comparison_label = contrast_name
    # de plot
    colours = c("Up" = "red", "Down" = "blue")
    plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,alpha=ifelse(significant=="yes",1,0.1),label=gene))+
    theme_bw()+
    geom_point()+
    NoLegend()+
    geom_label_repel(data=res %>% arrange(PValue) %>% head(n=30),mapping=aes())+
    ggtitle(paste0(cell_type,":", comparison_label))+
    scale_color_manual(values = colours)

    png(paste0("de_plot_",contrast_name,"_",cell_type,".png"),res=300,units="in",height=8,width=8)
    print(plot)
    dev.off()

    qlf_estimates = res %>%  dplyr::select(gene,PValue,logFC)


    # stash results
    res_list = data.frame(
      "Cell type" = cell_type,
      "Comparison" = contrast_name,
      "QLF - N genes tested" = total_genes,
      "QLF - Non-significant (FDR >0.01)" = non_sig,
      "QLF - Up-regulated (FDR<0.01)" = sig_up,
      "QLF - Down-regulated (FDR<0.01)" = sig_down)

    results_df <<- bind_rows(results_df,res_list)
  }
  }
})


results_df[,3:ncol(results_df)] = sapply(results_df[,3:ncol(results_df)], as.numeric)


results_file = results_df %>% mutate(Proportion_significant = (`QLF...Up.regulated..FDR.0.01.` + `QLF...Up.regulated..FDR.0.01.`) / `QLF...N.genes.tested`)
write_csv(results_file,"overall_de_edger_results.csv")
results_df = results_df %>% melt(id.vars=c("Cell.type","Comparison"))
results_df$variable = str_remove(results_df$variable,pattern="QLF...")
results_df = results_df %>% filter(!grepl("N.genes.tested",variable))

p1=ggplot(results_df,aes(Cell.type,value,fill=variable,label=value))+
facet_wrap(~Comparison,ncol=3)+
geom_col()+theme_classic() + theme(axis.text.x=element_text(angle=90))+labs(y="n genes",x="Cell type", fill="Key")+scale_fill_brewer(palette="Set2")
png(paste0("overall_de_gene_count_plot.png"),res=300,units="in",height=8,width=8)
p1
dev.off()

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

clusters = t_cells$cell_type %>% unique

nperm = 100000
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
      file = paste0("edgeR_de_tests_",comparison,"_",cluster,".csv")
      if(!file.exists(file)){
        next
      }
      de = read_csv(file)
      ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
      ranked_genes_vector = ranked_genes$logFC
      names(ranked_genes_vector) = ranked_genes$gene
      message("There are ",nrow(ranked_genes)," genes in this analysis")

      # do gsea
      res = fgsea(genelist, stats = ranked_genes_vector, minSize=10,eps=0, nPermSimple = nperm)
      sig_pathways = res %>% filter(padj<0.1)

      if(nrow(sig_pathways)>0){
        leading_edges = list()
        for(i in c(1:nrow(sig_pathways))){
          leading_edge = sig_pathways$leadingEdge[[i]] %>% paste0(collapse=", ")
          leading_edges[[i]] = leading_edge
        }
        leading_edges = unlist(leading_edges)
        sig_pathways$leadingEdge = leading_edges

        write_csv(sig_pathways,file=paste0("sig_pathways_",geneset,"_",cluster,"_",comparison,".csv"))
      }

      res = res %>% arrange(padj) %>%  dplyr::select(1,2,3,5,6) %>%  mutate(ann_monaco = cluster) %>%  mutate(comparison = this_comparison) %>% data.frame()
      res_overall[[(length(res_overall)+1)]] <<- res
      res = res %>% arrange(NES)
      res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
      topres = res %>% arrange(padj) %>% head(n=10) %>% arrange(desc(NES))
      plot = ggplot(topres,aes(NES,pathway,label=ifelse(padj<0.05,"*"," "),alpha=-log10(pval),fill=ifelse(NES>0,paste0("Up in ",this_comparison),paste0("Down in ",this_comparison))))+
        geom_col()+
        theme_bw()+
        ggtitle(paste0("Comparison of ",cluster," in ",this_comparison)) +
        labs(x="Normalised enrichment score",fill=paste0("Up or down-regulated in \n ",this_comparison))+
        geom_text()

      png(paste0("gsea_results_",geneset,"_",cluster,"_",comparison,".png"),res=300,units="in",height=8,width=8)
      print(plot)
      dev.off()
    }
  }
}

do_gsea("hallmark")
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"res_overall.csv")
res_overall$fdr = p.adjust(res_overall$pval,method="fdr")
for(x in c( "MS CSF vs MS PBMC","MS CSF vs Control CSF","Control CSF vs Control PBMC","MS PBMC vs Control PBMC","MS PBMC vs OIND PBMC","MS CSF vs OIND CSF" )){
  plot_data = res_overall %>% filter(comparison == x) %>% filter(grepl("HALLMARK",pathway)) %>% arrange(fdr)
  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  png(paste0("summary_",x,".png"),res=300,units="in",width=6,height=8)
  p=ggplot(plot_data,aes(ann_monaco,pathway,fill=NES,label=ifelse(fdr<0.1,"*","")))+
    geom_tile()+
    geom_text()+
    scale_fill_viridis_c()+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)+
    theme(axis.text.x=element_text(angle=90))
  print(p)
  dev.off()
}

res_overall = list()
do_gsea("go")
res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"go_res_overall.csv")
res_overall$fdr = p.adjust(res_overall$pval,method="fdr")
for(x in c( "MS CSF vs MS PBMC","MS CSF vs Control CSF","Control CSF vs Control PBMC","MS PBMC vs Control PBMC","MS PBMC vs OIND PBMC","MS CSF vs OIND CSF" )){
  plot_data = res_overall %>% filter(comparison == x) %>% filter(grepl("GO",pathway)) %>% arrange(fdr)
  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  png(paste0("go_summary_",x,".png"),res=300,units="in",width=6,height=8)
  p=ggplot(plot_data,aes(ann_monaco,pathway,fill=NES,label=ifelse(fdr<0.1,"*","")))+
    geom_tile()+
    geom_text()+
    scale_fill_viridis_c()+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)+
    theme(axis.text.x=element_text(angle=90))
  print(p)
  dev.off()
}

#do_gsea("kegg")
#do_gsea("reactome")



#######################################
# Clonal analysis
#######################################

# redefine clones ID using just TRB clone ID

create_trb_clone_id = function(x){
  full_id = x
  split_id = str_split(full_id,pattern="_")
  new_id = paste0(split_id[[1]][c(1:3)],collapse="_")
  return(new_id)
}
t_cells@meta.data$trb_clone_id = sapply(t_cells@meta.data$clone_id, create_trb_clone_id)
t_cells@meta.data$clone_id = t_cells@meta.data$trb_clone_id

# define clones
expanded_clones = t_cells@meta.data %>% group_by(donor.id,clone_id) %>% dplyr::count() %>% filter(n>1) %>% mutate(donor_clone = paste0(donor.id,"_",clone_id))

t_cells@meta.data = t_cells@meta.data %>% mutate(donor_clone = paste0(donor.id,"_",clone_id)) %>%
mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

# big picture numbers
table(t_cells@meta.data$expanded_clone)
t_cells@meta.data %>% distinct(donor_clone) %>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n>1)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n==2)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n>5)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n==1)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% arrange(desc(n)) %>% head

# recalculate clonal size per donor
clonal_size = t_cells@meta.data %>% dplyr::count(donor_clone) %>% dplyr::rename("clonal_size" = "n")
t_cells@meta.data = t_cells@meta.data %>% left_join(clonal_size,by="donor_clone")


# donors with expanded clones
donors_with_expanded_clones = t_cells@meta.data  %>% group_by(donor.id) %>% dplyr::count(donor_clone) %>% filter(n>1) %>% distinct(donor.id)
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_with_expanded_clones = ifelse(donor.id %in% donors_with_expanded_clones$donor.id,"Yes","No"))

# histogram
p=ggplot(t_cells@meta.data,aes(clonal_size,fill=phenotype))+geom_histogram()+theme_bw()+scale_fill_brewer(palette="Set2")+labs(x="Clone size")
png("clonal_histogram.png",res=300,height=2,width=4,units="in")
p
dev.off()

# clonality of csf and pbmc
p1=ggplot(t_cells@meta.data,aes(phenotype,fill=expanded_clone))+
geom_bar(position="fill")+
facet_wrap(~source)+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Expanded clone?",x="Phenotype",y="Proportion of T cell pool")

png("clonal_proportions_source.png",res=300,height=4,width=7,units="in")
p1
dev.off()

# quantify clonality of csf
counts = t_cells@meta.data %>% group_by(source,phenotype) %>% dplyr::count(expanded_clone)
totals = t_cells@meta.data %>% group_by(source,phenotype) %>% dplyr::count(expanded_clone) %>% summarise(totals = sum(n))
counts %>% left_join(totals,by=c("source","phenotype")) %>% mutate(prop = n/totals)

#######################################
# Clonal phenotypes
#######################################

colour_pal = RColorBrewer::brewer.pal(12, "Set2")
colour_pal = grDevices::colorRampPalette(colour_pal)(12)

p1=ggplot(t_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,fill=cell_type))+
geom_bar(position="fill")+
scale_fill_manual(values=colour_pal)+
theme_bw()+
labs(fill="Cell type",x="Expanded clone",y="Proportion of T cell pool")

p2=ggplot(t_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,fill=source))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Source",x="Expanded clone",y="Proportion of T cell pool")

# phenotypes of clonal cells - TRBV types
t_cells@meta.data$trbv_family = sapply(t_cells@meta.data$v_call_VDJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  return(y)
})


no_trbv_families = t_cells@meta.data$trbv_family %>% unique %>% length
colour_pal = RColorBrewer::brewer.pal(12, "Set2")
colour_pal = grDevices::colorRampPalette(colour_pal)(no_trbv_families)

p3=ggplot(t_cells@meta.data %>% filter(phenotype=="MS") %>% distinct(donor_clone,.keep_all=TRUE),aes(expanded_clone,fill=trbv_family))+
geom_bar(position="fill")+
scale_fill_manual(values = colour_pal)+
theme_bw()+
labs(fill="TRBV family",x="Expanded clone",y="Proportion of T cell pool")

p4=ggplot(t_cells@meta.data %>% filter(phenotype=="MS") %>% distinct(donor_clone,.keep_all=TRUE),aes(expanded_clone,fill=status_summary))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="TR chain pairing",x="Expanded clone",y="Proportion of T cell pool")

# cdr3 length
t_cells@meta.data = t_cells@meta.data %>% mutate(cdr3_length = nchar(junction_aa_VDJ))
p5 = ggplot(t_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,cdr3_length,fill=source))+
geom_boxplot()+
theme_bw()+
labs(fill="Source",x="Expanded clone",y="CDR3 length")

# dim plot highlighting expanded clones
rownames(t_cells@meta.data) = colnames(t_cells)
p6=DimPlot(subset(t_cells,subset = phenotype == "MS"),group.by="expanded_clone")
png("clonal_phenotypes.png",res=300,height=10,width=14,units="in")
grid.arrange(p1,p2,p3,p4,p5,p6+ggtitle(""))
dev.off()

#######################################
# Clonal relationships
#######################################

# find public clones in the dataset

public_clones = t_cells@meta.data %>% group_by(clone_id,donor.id) %>% dplyr::count() %>% group_by(clone_id) %>% dplyr::count(clone_id) %>% arrange(desc(n)) %>% filter(n>1)
t_cells@meta.data = t_cells@meta.data %>% mutate(private_status = ifelse(clone_id %in% public_clones$clone_id,"Public","Private"))

ms_specific_clones = sapply(public_clones$clone_id, function(clone){
  this_clone = t_cells@meta.data %>% filter(clone_id == clone)
  ms_specific = all(this_clone$phenotype == "MS")
  ms_specific
})
public_clones$ms_specific = ms_specific_clones

# lookup in database
specific_clones = public_clones %>% filter(ms_specific==TRUE)
cdr3s = t_cells@meta.data %>%
filter(clone_id %in% specific_clones$clone_id) %>%
distinct(clone_id,.keep_all=TRUE) %>%
dplyr::select(contains("junction"),clone_id,v_call_VDJ,j_call_VDJ)
specific_clones = specific_clones %>% left_join(cdr3s,by="clone_id")

# vdj DB
# downloaded from https://github.com/antigenomics/vdjdb-db/releases/tag/2021-09-05
vdj_db = read_tsv("vdjdb.txt")

specific_clones %>% left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope)) %>%
distinct(junction_aa_VDJ,.keep_all=TRUE) %>%
dplyr::count(antigen.species)

specific_clones_overlap_with_known_cdr3 = specific_clones %>% left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope)) %>%
distinct(junction_aa_VDJ,.keep_all=TRUE)
write_csv(specific_clones_overlap_with_known_cdr3,"specific_clones_overlap_with_known_cdr3.csv")

# whole dataset
tcr_db_gex = t_cells@meta.data %>% left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope)) %>%
dplyr::count(antigen.species,phenotype)

totals = t_cells@meta.data %>% group_by(phenotype) %>% dplyr::count() %>% dplyr::rename("total_cells" = n)
tcr_db_gex = tcr_db_gex %>% left_join(totals,by="phenotype")

p1=ggplot(tcr_db_gex,aes(n,antigen.species,fill=phenotype))+geom_col(position=position_dodge())
p2=ggplot(tcr_db_gex,aes(n/total_cells,antigen.species,fill=phenotype))+geom_col(position=position_dodge())


png("pathogen_tcrs.png",res=300,units="in",height=6,width=12)
grid.arrange(p1,p2,nrow=1)
dev.off()


# repeat for clonal vs nonclonal in ms
# whole dataset
tcr_db_gex = t_cells@meta.data %>% filter(phenotype=="MS") %>%
left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope)) %>%
dplyr::count(antigen.species,expanded_clone)

totals = t_cells@meta.data %>% group_by(expanded_clone) %>% dplyr::count() %>% dplyr::rename("total_cells" = n)
tcr_db_gex = tcr_db_gex %>% left_join(totals,by="expanded_clone")

p1=ggplot(tcr_db_gex,aes(n,antigen.species,fill=expanded_clone))+geom_col(position=position_dodge())
p2=ggplot(tcr_db_gex,aes(n/total_cells,antigen.species,fill=expanded_clone))+geom_col(position=position_dodge())

png("pathogen_tcrs_by_expanded_just_ms.png",res=300,units="in",height=6,width=12)
grid.arrange(p1,p2,nrow=1)
dev.off()

# repeat for clonal vs nonclonal in controls
# whole dataset
tcr_db_gex = t_cells@meta.data %>% filter(phenotype=="Noninflammatory") %>%
left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope)) %>%
dplyr::count(antigen.species,expanded_clone)

totals = t_cells@meta.data %>% group_by(expanded_clone) %>% dplyr::count() %>% dplyr::rename("total_cells" = n)
tcr_db_gex = tcr_db_gex %>% left_join(totals,by="expanded_clone")

p1=ggplot(tcr_db_gex,aes(n,antigen.species,fill=expanded_clone))+geom_col(position=position_dodge())
p2=ggplot(tcr_db_gex,aes(n/total_cells,antigen.species,fill=expanded_clone))+geom_col(position=position_dodge())

png("pathogen_tcrs_by_expanded_just_control.png",res=300,units="in",height=6,width=12)
grid.arrange(p1,p2,nrow=1)
dev.off()

# by source
tcr_db_gex = t_cells@meta.data %>% left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope))

ebv_counts = tcr_db_gex %>% filter(antigen.species=="EBV") %>% dplyr::count(source,phenotype,donor.id,cell_type)

# permutation
# whole dataset


bootstrap_ebv_enrichment = function(pheno1 = "Noninflammatory",
clonal1 =  "Expanded",
pheno2 = "Noninflammatory",
clonal2 = "Not expanded")
{
  reference_data = t_cells@meta.data %>%
  filter(phenotype==pheno1 & expanded_clone==clonal1)
  n_cells = reference_data %>% nrow

  reference_data = reference_data %>%
  left_join(
  vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
  by = c("junction_aa_VDJ")) %>%
  filter(gene=="TRB") %>%
  filter(!is.na(antigen.epitope)) %>%
  dplyr::count(antigen.species=="EBV")

  n_ebv_cells_ref = reference_data[2,2]

  message(n_ebv_cells_ref, " EBV+ cells ")
  message(n_cells, " total cells ")

  results = list()
  for(i in c(1:1000)){
  message("iteration ",i)
  test_data = t_cells@meta.data %>%
  filter(phenotype==pheno2 & expanded_clone==clonal2) %>%
  sample_n(size = n_cells,replace=FALSE)
  test_data = test_data %>%
  left_join(
  vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
  by = c("junction_aa_VDJ")) %>%
  filter(gene=="TRB") %>%
  filter(!is.na(antigen.epitope)) %>%
  dplyr::count(antigen.species=="EBV")
  n_ebv_cells = test_data[2,2]
  message(n_ebv_cells, " + EBV cells in this run")

  if(!is.na(n_ebv_cells)){
    results[[i]] = n_ebv_cells
  } else {
    results[[i]] = 0
  }
  }
  results = unlist(results)
  pval = sum(results>n_ebv_cells_ref)/1000
  message(pval)
}

bootstrap_ebv_enrichment(pheno1 = "Noninflammatory",clonal1 =  "Expanded",pheno2 = "Noninflammatory",clonal2 = "Not expanded")
bootstrap_ebv_enrichment(pheno1 = "MS",clonal1 =  "Expanded",pheno2 = "MS",clonal2 = "Not expanded")
bootstrap_ebv_enrichment(pheno1 = "OIND",clonal1 =  "Expanded",pheno2 = "OIND",clonal2 = "Not expanded")








################################
# de clonal vs not
################################

DefaultAssay(t_cells) = "RNA"
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_expanded = paste0(donor.id,"_",expanded_clone))
rownames(t_cells@meta.data) = colnames(t_cells)



for(cluster in unique(t_cells@meta.data$cell_type)){
message("Cluster: ",cluster)

cells_for_de = subset(t_cells, subset = phenotype=="MS" & source == "CSF" & cell_type == cluster)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = cells_for_de@meta.data %>%
  group_by(donor.id,expanded_clone) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(donor.id,"_",expanded_clone))


# convert to sce object
cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

# aggregate counts
groups = colData(cells_for_de.sce)[, c("donor.id","expanded_clone")]
aggregated_counts  = aggregate.Matrix(t(counts(cells_for_de.sce)),
groupings = groups, fun = "sum") %>% t()

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]

group_vector = lapply(colnames(aggregated_counts),function(y){
      if(grepl("Not expanded",y)){
        "Not_expanded"
      } else if(grepl("Expanded",y)){
        "Expanded"
      }
    }) %>% unlist %>% factor()

# make the DGE object
y=DGEList(aggregated_counts,group=group_vector,remove.zeros=TRUE)

# update sample info
y$samples =  y$samples %>%
mutate(donor_expanded = rownames(y$samples)) %>%
left_join(t_cells@meta.data %>%
dplyr::select(Age,Gender,drb1_1501_dose,donor_expanded) %>%
distinct(donor_expanded,.keep_all=TRUE),by="donor_expanded")

design = model.matrix(~0+group_vector+Age+Gender,y$samples)
colnames(design) = c(levels(group_vector),"Age","Gender")

keep = filterByExpr(
  y,
  design = design,
  group = group_vector,
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7)
y = y[keep, , keep.lib.sizes=FALSE]
y = calcNormFactors(y)
y = estimateDisp(y,design,robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
contrast_to_test = makeContrasts("Expanded - Not_expanded",levels = design)
res = glmQLFTest(fit, contrast = contrast_to_test)

print(summary(decideTests(res)))
res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
res$gene = rownames(res)
res$significant = ifelse(res$P_adj<0.05,"yes","no")
res$direction = ifelse(res$logFC>0,"Up","Down")
total_genes = nrow(res)
non_sig = sum(res$significant=="no")
sig_up = sum(res$significant=="yes" & res$direction=="Up")
sig_down = sum(res$significant=="yes" & res$direction=="Down")


# de plot
colours = c("Up" = "red","Down" = "blue")
plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,alpha=ifelse(significant=="yes",1,0.1),label=gene))+
theme_bw()+
geom_point()+
NoLegend()+
scale_color_manual(values = colours)+
geom_label_repel(data=res %>% arrange(PValue) %>% head(n=50),mapping=aes())+
guides(alpha=FALSE)

png(paste0(cluster,"_expanded_vs_not_csf_ms.png"),res=300,units="in",height=4,width=4)
print(plot)
dev.off()


de = res
ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
ranked_genes_vector = ranked_genes$logFC
names(ranked_genes_vector) = ranked_genes$gene
message("There are ",nrow(ranked_genes)," genes in this analysis")
write_csv(res,file=paste0(cluster,"_expanded_vs_not_csf_ms.csv"))

# do gsea
nperm=10000
geneset="hallmark"
genelist = eval(parse(text = paste0(geneset,"_list")))
res = fgsea(genelist, stats = ranked_genes_vector, minSize=10,eps=0, nPermSimple = nperm)

res = res %>% arrange(NES)
res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
topres = res %>% arrange(padj) %>% head(n=10) %>% arrange(desc(NES))
plot = ggplot(topres,aes(NES,pathway,label=ifelse(padj<0.1,"*"," "),alpha=-log10(pval),fill=ifelse(NES>0,"Up in expanded clones","Down in expanded clones")))+
  geom_col()+
  theme_bw()+
  labs(x="Normalised enrichment score",fill=paste0("Up or down-regulated in expanded clones"))+
  geom_text()+
  scale_fill_manual(values=c("blue","red"))+
  guides(alpha=FALSE)

png(paste0(cluster,"_gsea_results_hallmark_expanded.png"),res=300,units="in",height=6,width=8)
print(plot)
dev.off()
}


################################
# de DRB1*15:01
################################

DefaultAssay(t_cells) = "RNA"
t_cells@meta.data = t_cells@meta.data %>% mutate(drb_status = ifelse(drb1_1501_dose >= 1,"Positive","Negative"))
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_hla = paste0(donor.id,"_",drb_status))

for(cluster in unique(t_cells@meta.data$cell_type)){
message("Cluster: ",cluster)

cells_for_de = subset(t_cells, subset = phenotype=="MS" & source == "CSF" & cell_type == cluster)

# get rid of NAs
donors_to_keep = cells_for_de@meta.data %>% filter(!is.na(drb_status))
cells_for_de = subset(cells_for_de, subset = donor.id %in% donors_to_keep$donor.id)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = cells_for_de@meta.data %>%
  group_by(donor.id,drb_status) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(donor.id,"_",drb_status))


# convert to sce object
cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

# aggregate counts
groups = colData(cells_for_de.sce)[, c("donor.id","drb_status")]
aggregated_counts  = aggregate.Matrix(t(counts(cells_for_de.sce)),
groupings = groups, fun = "sum") %>% t()

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]

group_vector = lapply(colnames(aggregated_counts),function(y){
      if(grepl("Negative",y)){
        "Negative"
      } else if(grepl("Positive",y)){
        "Positive"
      }
    }) %>% unlist %>% factor()

# make the DGE object
y=DGEList(aggregated_counts,group=group_vector,remove.zeros=TRUE)

# update sample info
y$samples =  y$samples %>%
mutate(donor_hla = rownames(y$samples)) %>%
left_join(t_cells@meta.data %>%
dplyr::select(Age,Gender,drb1_1501_dose,donor_hla) %>%
distinct(donor_hla,.keep_all=TRUE),by="donor_hla")

design = model.matrix(~0+group_vector+Age+Gender,y$samples)
colnames(design) = c(levels(group_vector),"Age","Gender")

keep = filterByExpr(
  y,
  design = design,
  group = group_vector,
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7)
y = y[keep, , keep.lib.sizes=FALSE]
y = calcNormFactors(y)
y = estimateDisp(y,design,robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
contrast_to_test = makeContrasts("Positive - Negative",levels = design)
res = glmQLFTest(fit, contrast = contrast_to_test)

print(summary(decideTests(res)))
res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
res$gene = rownames(res)
res$significant = ifelse(res$P_adj<0.05,"yes","no")
res$direction = ifelse(res$logFC>0,"Up","Down")
total_genes = nrow(res)
non_sig = sum(res$significant=="no")
sig_up = sum(res$significant=="yes" & res$direction=="Up")
sig_down = sum(res$significant=="yes" & res$direction=="Down")


# de plot
colours = c("Up" = "red","Down" = "blue")
plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,alpha=ifelse(significant=="yes",1,0.1),label=gene))+
theme_bw()+
geom_point()+
NoLegend()+
scale_color_manual(values = colours)+
geom_label_repel(data=res %>% arrange(PValue) %>% head(n=50),mapping=aes())+
guides(alpha=FALSE)

png(paste0(cluster,"_hla_drb1_15_pos_v_neg_csf_ms.png"),res=300,units="in",height=6,width=6)
print(plot)
dev.off()


de = res
ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
ranked_genes_vector = ranked_genes$logFC
names(ranked_genes_vector) = ranked_genes$gene
message("There are ",nrow(ranked_genes)," genes in this analysis")
write_csv(res,file=paste0(cluster,"_expanded_vs_not_csf_ms.csv"))
}
