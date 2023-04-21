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
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/bcr/")

# Read in data
b_cells = readRDS("../b_cells_post_processing.rds")
rownames(b_cells@meta.data) = colnames(b_cells)


#######################################
# Big picture descriptive stats
#######################################
message("No. cells:")
nrow(b_cells@meta.data)

#######################################
# Recluster and reannotate
#######################################

b_subsets_monaco = c("Exhausted B cells","Naive B cells","Non-switched memory B cells","Plasmablasts","Switched memory B cells")
b_subsets_blueprint = c("B-cells","Class-switched memory B-cells","Memory B-cells","naive B-cells","Plasma cells")
b_cell_markers = c("IGHD","IGHG1","IGHM","CD27","CD38","MKI67")

# filter by annotation
filtered_bcells = subset(b_cells, subset = ann_monaco %in% b_subsets_monaco)

# recluster
set.seed(1)
DefaultAssay(filtered_bcells)="SCT"
filtered_bcells = RunUMAP(filtered_bcells,reduction="harmony",dims=1:50)
filtered_bcells = FindNeighbors(filtered_bcells)
filtered_bcells = FindClusters(filtered_bcells,resolution=0.5)
filtered_bcells = SetIdent(filtered_bcells,value = filtered_bcells[['ann_monaco']])

p1=DimPlot(filtered_bcells)
p2=FeaturePlot(filtered_bcells,features=b_cell_markers)
p3=DotPlot(filtered_bcells,features=b_cell_markers)

png("dimplot.png",res=300,units="in",width=6,height=4)
p1
dev.off()

png("featureplot.png",res=300,units="in",width=6,height=6)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=8,height=4)
p3
dev.off()

# recode blueprint annotations for plotting (simplify)
non_b_cell_annotations = c("CD4+ T-cells","CD4+ Tem","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","CLP","Erythrocytes","Megakaryocytes","Monocytes","NK cells","Tregs")
filtered_bcells@meta.data = filtered_bcells@meta.data %>%
mutate(ann_blueprint_simple  = ifelse(ann_blueprint %in% non_b_cell_annotations,"Non-B cell",ann_blueprint)) %>%
mutate(ann_blueprint_simple  = ifelse(is.na(ann_blueprint_simple),"Non-B cell",ann_blueprint_simple))
filtered_bcells@meta.data$ann_blueprint_simple = as.factor(filtered_bcells@meta.data$ann_blueprint_simple)

p1=ggplot(filtered_bcells@meta.data,aes(ann_monaco,fill=ann_blueprint_simple))+geom_bar(position="fill")+labs(x="Monaco annotation",y="Proportion",fill="Blueprint/ENCODE annotation")+theme_bw()+scale_fill_brewer(palette="Set2")

png("monaco_v_blueprint_annotations.png",res=300,height=4,width=10,units="in")
p1
dev.off()

b_cells = filtered_bcells

# export for azimuth
b_cells_azimuth = DietSeurat(b_cells)
DefaultAssay(b_cells_azimuth) = "RNA"
b_cells_azimuth[['SCT']] = NULL
saveRDS(b_cells_azimuth,"bcells_for_azimuth.rds")

# load azimuth annotations back in
preds = read_tsv("azimuth_pred.tsv")
preds = preds %>% filter(mapping.score>0.5)
b_cells_azimuth_preds = b_cells@meta.data %>% filter(cell_id %in% preds$cell) %>% left_join(preds %>% dplyr::rename("cell_id" = "cell"))

# look at azimuth annotations
totals = b_cells_azimuth_preds %>% group_by(ann_monaco) %>% dplyr::count(predicted.celltype.l2) %>% summarise(total = sum(n))

azimuth_calls = b_cells_azimuth_preds %>% group_by(ann_monaco) %>%
dplyr::count(predicted.celltype.l3) %>%
left_join(totals,by="ann_monaco") %>%
mutate(prop = n/total*100) %>% filter(prop > 5) %>%
filter(!is.na(predicted.celltype.l3))

write_csv(azimuth_calls,"azimuth_cluster_calls.csv")
p=ggplot(azimuth_calls,aes(ann_monaco,prop,fill=predicted.celltype.l3))+geom_col()+theme_bw()
png("azimuth_calls.png",res=300,units="in",width=10,height=6)
p
dev.off()

# cluster biomarkers
biomarkers = FindAllMarkers(b_cells,recorrect_umi=FALSE,logfc.threshold=0,min.pct=0,only.pos=TRUE)
topbiomarkers = biomarkers %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=5)
write_csv(topbiomarkers,"biomarkers.csv")
write_csv(biomarkers,"all_biomarkers.csv")

#######################################
# DA & composition
#######################################

n_col = b_cells@meta.data$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Set2")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

b_cells@meta.data$cell_type = factor(b_cells@meta.data$cell_type,ordered=TRUE)
p=DimPlot(b_cells,raster=F,group.by="cell_type")+scale_color_manual(values = colour_pal)+ggtitle("")
png("dim_plot_simple_labels.png",res=300,units="in",width=8,height=6)
p
dev.off()

p=DimPlot(subset(b_cells, subset = phenotype=="MS"),raster=F,group.by="cell_type",split.by="source")+
scale_color_manual(values = colour_pal)+ggtitle("")
png("dim_plot_simple_labels_just_ms_by_source.png",res=300,units="in",width=6,height=4)
p
dev.off()


b_cells@meta.data$phenotype = factor(b_cells@meta.data$phenotype,levels=c("Noninflammatory","OIND","MS"),ordered=T)

png("proportions.png",res=300,width=8,height=4,units="in")
ggplot(b_cells@meta.data,aes(phenotype,fill=cell_type))+
scale_fill_manual(values = colour_pal)+
geom_bar(position="fill")+
facet_wrap(~source)+
theme_bw()+
labs(fill="B cell subset",x="Phenotype",y="Proportion of B cell pool")
dev.off()


png("proportions_just_ms.png",res=300,width=8,height=4,units="in")
ggplot(b_cells@meta.data %>% filter(phenotype=="MS"),aes(source,fill=cell_type))+
scale_fill_manual(values = colour_pal)+
geom_bar(position="fill")+
theme_bw()+
labs(fill="B cell subset",x="Phenotype",y="Proportion of B cell pool")
dev.off()

# DA testing
# create new unique ID with donor and source
b_cells@meta.data$cell_type = Idents(b_cells)
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
b_cells@meta.data = b_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# stash sample info
sample_info = b_cells@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(b_cells@meta.data$cell_type,b_cells@meta.data$donor_source)

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
y.ab$samples =  y.ab$samples %>% left_join(b_cells@meta.data %>% distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Gender,drb1_1501_dose,donor_source),by="donor_source")

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

  # refactor cell types (for plotting)
  res$cell = factor(res$cell,ordered=TRUE,levels = levels(factor(b_cells@meta.data$cell_type %>% unique,ordered=TRUE)))

  # da plot
  plot = ggplot(res,aes(logFC,-log10(PValue),label=cell,fill=cell))+
  theme_classic()+
  scale_fill_manual(values = colour_pal)+
  geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
  geom_vline(xintercept = 0,alpha=0.2)+
  NoLegend()+
  geom_label_repel(data=res,mapping=aes(),max.overlaps=100)+
  ggtitle(comparison_label)+
  scale_y_continuous(limits=c(0,30))+
  scale_x_continuous(limits=c(-5,5))+
  geom_point(shape=16)

  png(paste0("da_plot_",contrast_name,"_.png"),res=300,units="in",height=4,width=4)
  print(plot)
  dev.off()
}

#######################################
# DE
#######################################


DefaultAssay(b_cells)="RNA"

# create new unique ID with donor and source
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
b_cells@meta.data = b_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# convert to sce object
b_cells.sce = as.SingleCellExperiment(b_cells)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = b_cells@meta.data %>%
  group_by(donor.id,source,phenotype,cell_type) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# aggregate counts
groups = colData(b_cells.sce)[, c("ident", "phenotype","source","donor.id")]
aggregated_counts  = aggregate.Matrix(t(counts(b_cells.sce)),
groupings = groups, fun = "sum")

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]


# get names for clusters
clusters = levels(factor(b_cells@meta.data$cell_type))

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
    y$samples =  y$samples %>% mutate(full_cell_id = rownames(y$samples)) %>% left_join(b_cells@meta.data %>% dplyr::select(Age,Gender,drb1_1501_dose,full_cell_id) %>% distinct(full_cell_id,.keep_all=TRUE),by="full_cell_id")

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



    bonf = 0.05 / length(res$PValue)
    res$significant = ifelse(res$PValue<bonf,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
    res = res %>% mutate(direction = ifelse(PValue>bonf,"nonsig",direction))

    total_genes = nrow(res)
    non_sig = sum(res$significant=="no")
    sig_up = sum(res$significant=="yes" & res$direction=="Up")
    sig_down = sum(res$significant=="yes" & res$direction=="Down")
    comparison_label = contrast_name
    # de plot
    colours = c("Up" = "red", "Down" = "blue", "nonsig" = "grey")
    plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,label=gene))+
    theme_bw()+
    geom_point()+
    NoLegend()+
    geom_label_repel(data=res %>% arrange(PValue) %>% head(n=10),mapping=aes(),max.overlaps=500)+
    ggtitle(paste0(cell_type,":", comparison_label))+
    scale_color_manual(values = colours)+
    scale_x_continuous(limits=c(-12,12))+
    scale_y_continuous(limits=c(0,20))

    png(paste0("de_plot_",contrast_name,"_",cell_type,".png"),res=300,units="in",height=6,width=6)
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

clusters = b_cells$ann_monaco %>% unique

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
  top_pathways = plot_data %>% arrange(desc(abs(NES))) %>% distinct(pathway) %>% head(8)
  plot_data = plot_data %>% filter(pathway %in% top_pathways$pathway)
  plot_data$pathway = str_remove(plot_data$pathway,"HALLMARK_")
  plot_data$ann_monaco = factor(plot_data$ann_monaco,levels = c("Naive B cells","Non-switched memory B cells","Switched memory B cells","Plasmablasts","Exhausted B cells"))

  plot_data$ann_monaco_simple = recode(plot_data$ann_monaco,
  "Naive B cells" = "N",
  "Non-switched memory B cells" = "USM",
  "Switched memory B cells" = "SM",
  "Plasmablasts" = "PB",
  "Exhausted B cells" = "EX") %>% factor()


  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  png(paste0("summary_",x,".png"),res=300,units="in",width=6,height=4)
  p=ggplot(plot_data,aes(ann_monaco_simple,pathway,fill=NES,label=ifelse(padj<0.1,ifelse(padj<0.01,ifelse(padj<0.001,"***","**"),"*"),"")))+
    geom_tile()+
    geom_text()+
    scale_fill_viridis_c()+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)
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


# MTORC1 heatmap

mtor_genes = hallmark_list$HALLMARK_MTORC1_SIGNALING
pb_mtor = read_csv("edgeR_de_tests_MS_CSF - MS_PBMC_Plasmablasts.csv")
pb_mtor = pb_mtor %>% filter(gene %in% mtor_genes) %>% filter(P_adj<0.05 & logFC > 0)


#######################################
# Clonal analysis
#######################################

# define clones
expanded_clones = b_cells@meta.data %>% group_by(donor.id,clone_id) %>% dplyr::count() %>% filter(n>1) %>% mutate(donor_clone = paste0(donor.id,"_",clone_id))

b_cells@meta.data = b_cells@meta.data %>% mutate(donor_clone = paste0(donor.id,"_",clone_id)) %>%
mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

# big picture numbers
table(b_cells@meta.data$expanded_clone)
b_cells@meta.data %>% distinct(donor_clone) %>% nrow
b_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n>1)%>% nrow
b_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n==2)%>% nrow
b_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n>5)%>% nrow
b_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n==1)%>% nrow
b_cells@meta.data %>% dplyr::count(donor_clone) %>% arrange(desc(n)) %>% head

# recalculate clonal size per donor
clonal_size = b_cells@meta.data %>% dplyr::count(donor_clone) %>% dplyr::rename("clonal_size" = "n")
b_cells@meta.data = b_cells@meta.data %>% left_join(clonal_size,by="donor_clone")

# donors with expanded clones
donors_with_expanded_clones = b_cells@meta.data  %>% group_by(donor.id) %>% dplyr::count(donor_clone) %>% filter(n>1) %>% distinct(donor.id)
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_with_expanded_clones = ifelse(donor.id %in% donors_with_expanded_clones$donor.id,"Yes","No"))

# histogram
p=ggplot(b_cells@meta.data,aes(clonal_size,fill=phenotype))+geom_histogram()+theme_bw()+scale_fill_brewer(palette="Set2")+labs(x="Clone size")
png("clonal_histogram.png",res=300,height=2,width=4,units="in")
p
dev.off()

# clonality of csf and pbmc
p1=ggplot(b_cells@meta.data %>% filter(source=="CSF"),aes(y=phenotype,fill=expanded_clone))+
geom_bar(position="fill")+
facet_wrap(~source)+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Expanded clone?",y="Phenotype",x="Proportion of B cell pool")+
scale_x_continuous(breaks = c(0,0.5,1))


png("clonal_proportions_source.png",res=300,height=2,width=4,units="in")
p1
dev.off()


# quantify clonality of csf
counts = b_cells@meta.data %>% group_by(source,phenotype) %>% dplyr::count(expanded_clone)
totals = b_cells@meta.data %>% group_by(source,phenotype) %>% dplyr::count(expanded_clone) %>% summarise(totals = sum(n))
counts %>% left_join(totals,by=c("source","phenotype")) %>% mutate(prop = n/totals)

# per individual
counts = b_cells@meta.data %>% group_by(source,phenotype,donor.id) %>% dplyr::count(expanded_clone)
totals = b_cells@meta.data %>% group_by(source,phenotype,donor.id) %>% dplyr::count() %>% mutate(totals = n) %>% dplyr::select(-n)
props = counts %>% left_join(totals,by=c("source","phenotype","donor.id")) %>% mutate(prop = n/totals)

ggplot(props %>% filter(expanded_clone=="Expanded"),aes(phenotype,prop))+geom_violin()+facet_wrap(~source)


#######################################
# Clonal phenotypes
#######################################

rownames(b_cells@meta.data) = colnames(b_cells)


png("clones_celltypes.png",res=300,height=3,width=6,units="in")
DimPlot(subset(b_cells, subset = phenotype=="MS"),split.by="expanded_clone")
dev.off()
png("clones_source.png",res=300,height=3,width=5,units="in")
DimPlot(subset(b_cells, subset = phenotype=="MS"),split.by="expanded_clone",group.by="source")
dev.off()
png("clones_isotype.png",res=300,height=3,width=5,units="in")
DimPlot(subset(b_cells, subset = phenotype=="MS"),split.by="expanded_clone",group.by="isotype")
dev.off()

p1=ggplot(b_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,fill=cell_type))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Cell type",x="Expanded clone",y="Proportion of B cell pool")

p2=ggplot(b_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,fill=source))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Source",x="Expanded clone",y="Proportion of B cell pool")

p3=ggplot(b_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,fill=shm_positive))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="SHM",x="Expanded clone",y="Proportion of B cell pool")

p4=ggplot(b_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,fill=isotype))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Isotype",x="Expanded clone",y="Proportion of B cell pool")

# phenotypes of clonal cells - IGHV types
b_cells@meta.data$ighv_family = sapply(b_cells@meta.data$v_call_genotyped_VDJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  return(y)
})

p5=ggplot(b_cells@meta.data %>% filter(phenotype=="MS") %>% distinct(donor_clone,.keep_all=TRUE),aes(expanded_clone,fill=ighv_family))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="IGHV family",x="Expanded clone",y="Proportion of B cell pool")

p6=ggplot(b_cells@meta.data %>% filter(phenotype=="MS") %>% distinct(donor_clone,.keep_all=TRUE),aes(expanded_clone,fill=status_summary))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Ig chain pairing",x="Expanded clone",y="Proportion of B cell pool")

color_pal = RColorBrewer::brewer.pal(12, "Set2")
color_pal = grDevices::colorRampPalette(color_pal)(15)
p7=ggplot(b_cells@meta.data %>% filter(phenotype=="MS") %>% distinct(donor_clone,.keep_all=TRUE),aes(expanded_clone,fill=c_call_VDJ))+
geom_bar(position="fill")+
scale_fill_manual(values=color_pal)+
theme_bw()+
labs(fill="IGHC call",x="Expanded clone",y="Proportion of B cell pool")

# cdr3 length
b_cells@meta.data = b_cells@meta.data %>% mutate(cdr3_length = nchar(junction_aa_VDJ))
p8 = ggplot(b_cells@meta.data %>% filter(phenotype=="MS"),aes(expanded_clone,cdr3_length,fill=source))+
geom_boxplot()+
theme_bw()+
labs(fill="Source",x="Expanded clone",y="CDR3 length")

# dim plot highlighting expanded clones
rownames(b_cells@meta.data) = colnames(b_cells)
p9=DimPlot(subset(b_cells,subset = phenotype == "MS"),group.by="expanded_clone")
png("clonal_phenotypes.png",res=300,height=10,width=14,units="in")
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9+ggtitle(""))
dev.off()

#######################################
# Clonal relationships
#######################################

# find public clones in the dataset

private_status_df = lapply(unique(b_cells@meta.data$clone_id),function(x){
  this_clone_data = b_cells@meta.data %>% filter(clone_id == x)
  this_clone_data = this_clone_data %>%
  mutate(private_status = ifelse(all(this_clone_data$donor.id == this_clone_data$donor.id[1]),"Private","Public"))
  this_clone_data
})
private_status_df = do.call("bind_rows",private_status_df) %>% distinct(clone_id,.keep_all=TRUE)

public_clone = private_status_df %>% filter(private_status=="Public") %>% dplyr::select(clone_id)
b_cells@meta.data %>% filter(clone_id %in% public_clone$clone_id)


bb_status_df = lapply(unique(b_cells@meta.data$donor_clone),function(x){
this_clone_data = b_cells@meta.data %>% filter(donor_clone == x)
this_clone_data = this_clone_data %>%
mutate(bb_status = ifelse(all(this_clone_data$source == "CSF"),"CSF only",NA)) %>%
mutate(bb_status = ifelse(all(this_clone_data$source == "PBMC"),"PBMC only",bb_status)) %>%
mutate(bb_status = ifelse("CSF" %in% this_clone_data$source & "PBMC" %in% this_clone_data$source ,"CSF & PBMC", bb_status))
this_clone_data
})
bb_status_df = do.call("bind_rows",bb_status_df) %>% distinct(donor_clone,.keep_all=TRUE)


p1=ggplot(bb_status_df %>% group_by(bb_status) %>% dplyr::count(),aes(bb_status,n,label=n,fill=bb_status))+geom_col()+geom_text(vjust=-1)+labs(x="Clone private to CSF or PBMC?",y="N clones")+NoLegend()
png("bb_status.png",res=300,height=6,width=6,units="in")
p1
dev.off()

# make clonal connection plot
shared = bb_status_df %>% filter(bb_status=="CSF & PBMC")
shared_cells = b_cells@meta.data %>% filter(donor_clone %in% shared$donor_clone)
counts = shared_cells %>% group_by(donor.id,source,ann_monaco,donor_clone) %>% dplyr::count()

clonal_plot = ggplot(counts,aes(source,ann_monaco,size=n,group=donor_clone,col=donor_clone))+geom_point()+geom_line()+facet_wrap(~donor.id)+theme_bw()+NoLegend()+labs(y="Cell type")
png("clonal_plot.png",res=300,height=6,width=4,units="in")
clonal_plot
dev.off()


#################################
# check for public clones
################################

bcr_db = read_csv("bcr_db.csv")
b_cells@meta.data %>% filter(junction_aa_VDJ %in% bcr_db$CDR3.heavy.aa)

################################
# de clonal vs not
################################

DefaultAssay(b_cells) = "RNA"
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_expanded = paste0(donor.id,"_",expanded_clone))
rownames(b_cells@meta.data) = colnames(b_cells)


leading_edges_mtorc1 = c()
gsea_results_overall = data.frame()
for(cluster in c("Exhausted B cells","Plasmablasts")){
message("Cluster: ",cluster)

cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == "CSF" & ann_monaco == cluster)

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
left_join(b_cells@meta.data %>%
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

bonf = 0.05 / length(res$PValue)
res$significant = ifelse(res$PValue<bonf,"yes","no")
res$direction = ifelse(res$logFC>0,"Up","Down")
res = res %>% mutate(direction = ifelse(PValue>bonf,"nonsig",direction))

total_genes = nrow(res)
non_sig = sum(res$significant=="no")
sig_up = sum(res$significant=="yes" & res$direction=="Up")
sig_down = sum(res$significant=="yes" & res$direction=="Down")
# de plot
colours = c("Up" = "red", "Down" = "blue", "nonsig" = "grey")
plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,label=gene))+
theme_bw()+
geom_point()+
NoLegend()+
geom_label_repel(data=res %>% arrange(PValue) %>% head(n=10),mapping=aes(),max.overlaps=500)+
ggtitle(cluster)+
scale_color_manual(values = colours)+
scale_x_continuous(limits=c(-7,7))+
scale_y_continuous(limits=c(0,10))

png(paste0(cluster,"_expanded_vs_not_csf_pcs_ms.png"),res=300,units="in",height=6,width=6)
print(plot)
dev.off()


de = res
ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
ranked_genes_vector = ranked_genes$logFC
names(ranked_genes_vector) = ranked_genes$gene
message("There are ",nrow(ranked_genes)," genes in this analysis")
write_csv(res,file=paste0(cluster,"_expanded_vs_not_csf_pcs_ms.csv"))

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

leading_edge_mtorc1 = res[res$pathway=="HALLMARK_MTORC1_SIGNALING",'leadingEdge']$leadingEdge %>% unlist()
leading_edges_mtorc1 <<- c(leading_edges_mtorc1,leading_edge_mtorc1)

res = data.frame(res)
res$cell_type = cluster
gsea_results_overall <<- rbind(gsea_results_overall,res)

png(paste0(cluster,"_gsea_results_hallmark_expanded_pcs.png"),res=300,units="in",height=6,width=8)
print(plot)
dev.off()
}

# summarise gsea results
res_overall = gsea_results_overall
res_overall$fdr = p.adjust(res_overall$pval,method="fdr")
plot_data = res_overall %>% filter(grepl("HALLMARK",pathway)) %>% arrange(fdr)
top_pathways = plot_data %>% arrange(desc(abs(NES))) %>% distinct(pathway) %>% head(8)
plot_data = plot_data %>% filter(pathway %in% top_pathways$pathway)
plot_data$pathway = str_remove(plot_data$pathway,"HALLMARK_")
plot_data$cell_type = factor(plot_data$cell_type,levels = c("Plasmablasts","Exhausted B cells"))
plot_data$cell_type = recode(plot_data$cell_type,
"Plasmablasts" = "PB",
"Exhausted B cells" = "EX") %>% factor()


  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  png("summary_gsea_clonal_ms_csf_gsea.png",res=300,units="in",width=5,height=3)
  p=ggplot(plot_data,aes(cell_type,pathway,fill=NES,label=ifelse(padj<0.1,ifelse(padj<0.01,ifelse(padj<0.001,"***","**"),"*"),"")))+
    geom_tile()+
    geom_text()+
    scale_fill_viridis_c()+
    theme_classic()+
    labs(x="Cell type",y="Gene set")
    print(p)
  dev.off()





unique_leading_edges_mtorc1 = unique(leading_edges_mtorc1)
pc_top_de = read_csv("Plasmablasts_expanded_vs_not_csf_pcs_ms.csv")
pc_genes = pc_top_de %>% arrange(desc(logFC),PValue) %>% filter(gene %in% leading_edges_mtorc1) %>% head(10)
DoHeatmap(subset(b_cells,subset = ann_monaco == "Plasmablasts"),features=pc_genes$gene,slot="data",group.by="expanded_clone")
ex_top_de = read_csv("Exhausted B cells_expanded_vs_not_csf_pcs_ms.csv")
ex_genes = ex_top_de %>% arrange(desc(logFC),PValue) %>% filter(gene %in% leading_edges_mtorc1) %>% head(10)

DoHeatmap(subset(b_cells,subset = ann_monaco == "Exhausted B cells"),features=ex_genes$gene,slot="data",group.by="expanded_clone")

#################################
# de clonal vs not broad clusters
#################################

pcs = subset(b_cells,cell_type_crude=="Plasma cells")

DefaultAssay(pcs) = "RNA"
pcs@meta.data = pcs@meta.data %>% mutate(donor_expanded = paste0(donor.id,"_",expanded_clone))
rownames(pcs@meta.data) = colnames(pcs)

cells_for_de = subset(pcs, subset = phenotype=="MS" & source == "CSF")

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
left_join(pcs@meta.data %>%
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

res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
res$gene = rownames(res)

bonf = 0.05 / length(res$PValue)
res$significant = ifelse(res$PValue<bonf,"yes","no")
res$direction = ifelse(res$logFC>0,"Up","Down")
res = res %>% mutate(direction = ifelse(PValue>bonf,"nonsig",direction))

total_genes = nrow(res)
non_sig = sum(res$significant=="no")
sig_up = sum(res$significant=="yes" & res$direction=="Up")
sig_down = sum(res$significant=="yes" & res$direction=="Down")
# de plot
colours = c("Up" = "red", "Down" = "blue", "nonsig" = "grey")
plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,label=gene))+
theme_bw()+
geom_point()+
NoLegend()+
geom_label_repel(data=res %>% arrange(PValue) %>% head(n=10),mapping=aes(),max.overlaps=500)+
ggtitle(cluster)+
scale_color_manual(values = colours)+
scale_x_continuous(limits=c(-7,7))+
scale_y_continuous(limits=c(0,10))

png("broad_clusters_pcs_expanded_vs_not_csf_pcs_ms.png",res=300,units="in",height=6,width=6)
print(plot)
dev.off()


de = res
ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
ranked_genes_vector = ranked_genes$logFC
names(ranked_genes_vector) = ranked_genes$gene
message("There are ",nrow(ranked_genes)," genes in this analysis")
write_csv(res,file="broad_clusters_pcs_expanded_vs_not_csf_pcs_ms.csv")

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

res = data.frame(res)
res$cell_type = cluster

png(paste0(cluster,"_gsea_results_hallmark_expanded_pcs.png"),res=300,units="in",height=6,width=8)
print(plot)
dev.off()

png("broad_clusters_pcs_expanded_vs_not_csf_pcs_ms_gsea_results_hallmark.png",res=300,units="in",height=6,width=8)
print(plot)
dev.off()



################################
# pathway scores
################################
mtorc1 = hallmark_list$HALLMARK_MTORC1_SIGNALING


add_hallmark_score = function(hallmark_pathway){
message(hallmark_pathway)
features = list(c(eval(parse(text=hallmark_pathway))))
# module score
DefaultAssay(b_cells) = "RNA"
b_cells = AddModuleScore(b_cells,features = features,nbin=10,ctrl=5,name=eval(hallmark_pathway))
b_cells[[eval(hallmark_pathway)]] = RNOmni::RankNorm(b_cells@meta.data[[paste0(eval(hallmark_pathway),"1")]])
b_cells
}

b_cells = add_hallmark_score("mtorc1")

png("vln_gene_scores_b_cells.png",res=300,units="in",width=4,height=4)
VlnPlot(subset(b_cells, subset = source=="CSF" & phenotype=="MS" & cell_type %in% c("Exhausted B cells","Plasmablasts")),features=c("mtorc1"),split.by="expanded_clone")+scale_fill_brewer(palette="Set3")+ggtitle("")+labs(y="MTORC1 gene score",x="")
dev.off()

compare_mtor_score = function(cell){
a = b_cells@meta.data %>% filter(cell_type == cell & source=="CSF" & phenotype=="MS") %>% filter(expanded_clone=="Expanded")
b = b_cells@meta.data %>% filter(cell_type == cell & source=="CSF" & phenotype=="MS") %>% filter(expanded_clone=="Not expanded")
t.test(a$mtorc1,b$mtorc1)$p.value
}

FeatureScatter(b_cells,feature1 = "mtorc1", feature2 = "MKI67",group.by="ann_monaco")


pcs = subset(b_cells,subset = ann_monaco %in% c("Plasmablasts"))
DefaultAssay(pcs)="SCT"
pcs = RunUMAP(pcs,reduction="harmony",dims=1:50)
pcs = FindNeighbors(pcs)
pcs = FindClusters(pcs,resolution=0.1)
DefaultAssay(pcs) = "RNA"
FeaturePlot(pcs,features="mtorc1")

png("eaf2.png",res=300,units="in",width=7,height=5)
FeaturePlot(subset(b_cells, subset = source == "CSF" & phenotype=="MS"),features="EAF2",split.by="expanded_clone")
dev.off()

VlnPlot(subset(b_cells,subset = phenotype == "MS" & source == "CSF" & cell_type %in% c("Exhausted B cells","Plasmablasts")),features="EAF2",split.by="expanded_clone")

VlnPlot(subset(b_cells,subset = phenotype == "MS" & source == "CSF"),features="EAF2",split.by="expanded_clone")
VlnPlot(subset(b_cells,subset = phenotype == "MS" & source == "CSF"),features="EAF2",group.by="expanded_clone")
VlnPlot(subset(b_cells,subset = phenotype == "MS" & source == "CSF" & ann_monaco =="Plasmablasts"),features="EAF2",group.by="expanded_clone")
VlnPlot(subset(b_cells,subset = phenotype == "MS" & source == "CSF" & ann_monaco =="Exhausted B cells"),features="EAF2",group.by="expanded_clone")


FeaturePlot(subset(b_cells,subset = phenotype == "OIND" & source == "CSF"),features="EAF2",split.by="expanded_clone")

################################
# de clonal vs not in single donors
################################


overall_results_df = data.frame()
# n=1 clonal vs not
DefaultAssay(b_cells) = "RNA"
# get names for clusters
clusters = levels(factor(b_cells@meta.data$ann_monaco))


donor_list = unique(b_cells@meta.data$donor.id)

do_de_per_donor = function(donor, cell_type = "Plasmablasts", source = "CSF"){
  message("Donor: ",donor)
  message("Source: ",source)
  message("Cell type:", cell_type)
  cells = b_cells@meta.data %>% filter(phenotype=="MS" & source == source & donor.id ==donor & ann_monaco == cell_type)
  no_expanded = cells %>% filter(expanded_clone=="Expanded")
  non_expanded = cells %>% filter(expanded_clone=="Not expanded")
  n_clones = no_expanded %>% distinct(clone_id) %>% nrow
  if(nrow(cells)==0 | nrow(no_expanded)<10 | nrow(non_expanded)==0){
    return(NA)
  }
  cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == source & donor.id ==donor & ann_monaco == cell_type)

  # distinct expanded clones
  expanded_clones = cells_for_de@meta.data %>% filter(expanded_clone=="Expanded") %>% distinct(donor_clone)

  for(clone in expanded_clones$donor_clone){

    # filter out other expanded cells
    cells_to_keep = cells_for_de@meta.data %>% filter(!(expanded_clone == "Expanded" & donor_clone!=clone))
    cells_for_de_filtered = subset(cells_for_de, subset = full_cell_id %in% cells_to_keep$full_cell_id)
    no_expanded = cells_to_keep %>% filter(expanded_clone=="Expanded") %>% nrow

    if(no_expanded<10){
      next
    } else {
    # convert to sce object
    cells_for_de.sce = as.SingleCellExperiment(cells_for_de_filtered)

    # aggregate counts
    groups = colData(cells_for_de.sce)[, c("donor.id","expanded_clone")]
    aggregated_counts  = aggregate.Matrix(t(counts(cells_for_de.sce)),
    groupings = groups, fun = "sum") %>% t()

    group_vector = lapply(colnames(aggregated_counts),function(y){
          if(grepl("Not expanded",y)){
            "Not_expanded"
          } else if(grepl("Expanded",y)){
            "Expanded"
          }
        }) %>% unlist %>% factor(levels=c("Not_expanded","Expanded"))
    # make the DGE object
    y=DGEList(aggregated_counts,group=group_vector,remove.zeros=TRUE)
    design = model.matrix(~0+group_vector)
    colnames(design) = levels(group_vector)

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
    res = exactTest(y,dispersion=0.4^2)
    print(summary(decideTests(res)))

    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$gene = rownames(res)
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res = res %>% mutate(direction = ifelse(logFC>0,"Up","Down")) %>% mutate(direction = ifelse(P_adj<0.05,direction,"neither"))

    total_genes = nrow(res)
    non_sig = sum(res$significant=="no")
    sig_up = sum(res$significant=="yes" & res$direction=="Up")
    sig_down = sum(res$significant=="yes" & res$direction=="Down")

    res$donor = donor
    res$source = source
    res$cell_type = cell_type
    res$donor_clone = clone
    overall_results_df <<- bind_rows(overall_results_df,res)

    # de plot
    colours = c("Up" = "red","Down" = "blue", "neither" = "grey")
    plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,alpha=ifelse(significant=="yes",1,0.1),label=gene))+
    theme_bw()+
    geom_point()+
    NoLegend()+
    geom_text_repel(data=res %>% arrange(PValue) %>% head(n=10),max_overlaps=50,mapping=aes())+
    ggtitle(paste0("Clone ID: ",clone,"\nExpanded n: ",no_expanded,"\nDonor: ",donor))+scale_color_manual(values = colours)
    return(plot)
    }
  }
}

plots = purrr::map(donor_list,do_de_per_donor)

plots = plots[!is.na(plots)]
null_vec = lapply(plots,function(x){
  is.null(x)
}) %>% unlist()
plots = plots[!null_vec]

png("pcs_clonal_vs_not_csf_single_donors.png",res=300,units="in",height=16,width=16)
grid.arrange(grobs=plots)
dev.off()

sig_results = overall_results_df %>% filter(P_adj < 0.05) %>% filter(!grepl("IG",gene))

write_csv(overall_results_df,"overall_de_expanded_v_not_single_donors_ms_csf.csv")
write_csv(sig_results,"sig_results_no_ig_overall_de_expanded_v_not_single_donors_ms_csf.csv")

colours = c("Up" = "red","Down" = "blue", "neither" = "grey")
plot = ggplot(sig_results,aes(logFC,-log10(PValue),color=direction,label=gene))+
theme_bw()+
facet_wrap(~donor_clone)+
geom_point()+
NoLegend()+
geom_text_repel(data=sig_results,mapping=aes())+
scale_color_manual(values = colours)

png("pcs_clonal_vs_not_csf_single_donors_no_ig.png",res=300,units="in",height=16,width=16)
plot
dev.off()

#######################################
# correlation with clinical phenotypes
#######################################


# merge with phenotype file
cam_pheno = read_csv("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/references/5PrimeCSF_ClinicalData_NoPID_211124.csv")

# sort out long cambridge IDs
new_ids = lapply(cam_pheno$ID, function(x){
  if(grepl("C00",x)){
      new_id = str_remove(pattern="C00TU0",x)
      new_id = str_split(pattern="a",new_id)[[1]][1]
      new_id = paste0("TU",new_id)
      new_id = str_remove(pattern="v1",new_id)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else if(grepl("v",x)){
      new_id = str_remove(pattern="v1",x)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else {
      return(x)
  }
})

cam_pheno$shortID = unlist(new_ids)
cam_hla = read_tsv("../../references/drb1_status.tsv")
cam_pheno = cam_pheno %>% left_join(cam_hla %>% dplyr::rename(shortID = donor_id) %>% dplyr::select(shortID,drb1_1501_dose)) %>% dplyr::select(-ID) %>% dplyr::rename(ID = shortID)

cam_pheno = b_cells@meta.data %>% filter(cohort=="Cam") %>% left_join(cam_pheno %>% dplyr::rename("donor.id" = "ID"),by="donor.id")

cam_pheno$OCB_status = recode(cam_pheno$Oligoclonal,
"Negative" = "Negative",
"Pos(>10)" = "Positive",
"Pos(3to9)" = "Negative")

counts = cam_pheno %>% group_by(donor.id,phenotype,source,OCB_status,Lym,Age.y,Gender.y,drb1_1501_dose.y,Meds) %>% dplyr::count()
expanded = cam_pheno %>% group_by(donor.id,phenotype,source,expanded_clone) %>% dplyr::count() %>% filter(expanded_clone=="Expanded")
combo = expanded %>% left_join(counts,by=c("donor.id","source","phenotype"))
combo = combo %>% mutate(proportion_expanded = n.x/n.y)

p1=ggplot(combo %>% filter(!is.na(OCB_status)) %>% filter(phenotype=="MS"),aes(OCB_status,proportion_expanded,fill=source))+geom_boxplot(position=position_dodge(width=1))+geom_point(alpha=0.1,position=position_dodge(width=1))+scale_fill_brewer(palette="Set2")


p2=ggplot(combo %>% filter(!is.na(drb1_1501_dose.y)),aes(drb1_1501_dose.y,proportion_expanded,fill=source))+geom_boxplot(position=position_dodge(width=1))+geom_point(alpha=0.1,position=position_dodge(width=1))+scale_fill_brewer(palette="Set2")
p3=ggplot(combo %>% filter(!is.na(drb1_1501_dose.y)),aes(factor(drb1_1501_dose.y),proportion_expanded,fill=source))+geom_boxplot(position=position_dodge(width=1))+geom_point(alpha=0.1,position=position_dodge(width=1))+scale_fill_brewer(palette="Set2")
p4=ggplot(combo ,aes(Age.y,proportion_expanded))+geom_point()+facet_wrap(~source)


##############################################
# drb1 eqtl
##############################################

# export plasma cells from ms csf only
csf_pcs = subset(b_cells,subset = source=="CSF" & phenotype=="MS" & ann_monaco=="Plasmablasts")

################################
# de drb1
################################

DefaultAssay(b_cells) = "RNA"
b_cells@meta.data = b_cells@meta.data %>% mutate(drb_status = ifelse(drb1_1501_dose >= 1,"Positive","Negative"))
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_hla = paste0(donor.id,"_",drb_status))

for(cluster in unique(b_cells@meta.data$cell_type)){
message("Cluster: ",cluster)

cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == "CSF" & ann_monaco == cluster)

# get rid of NAs
donors_to_keep = cells_for_de@meta.data %>% filter(!is.na(drb_status))
cells_for_de = subset(cells_for_de, subset = donor.id %in% donors_to_keep$donor.id)

#cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == "CSF" & ann_monaco %in% c("Plasmablasts","Exhausted B cells"))

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
left_join(b_cells@meta.data %>%
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
geom_label_repel(data=res %>% arrange(PValue) %>% head(n=50),mapping=aes())

png(paste0(cluster,"_hla_drb1_15_pos_v_neg_csf_ms.png"),res=300,units="in",height=8,width=8)
print(plot)
dev.off()


de = res
ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
ranked_genes_vector = ranked_genes$logFC
names(ranked_genes_vector) = ranked_genes$gene
message("There are ",nrow(ranked_genes)," genes in this analysis")
write_csv(res,file=paste0(cluster,"_hla_pos_vs_not_csf_pcs_ms.csv"))
}


#################################
# EAF2
#################################
genos = read_table("~/eaf_genotypes_ad.raw")
codex = read_table("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/data/genotypes/GT_ID.txt",col_names=FALSE)

genos = genos %>% left_join(codex %>% rename("X1" = "IID"),by="IID") %>% dplyr::select(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE)

# sort out long cambridge IDs
new_ids = lapply(genos$X2, function(x){
  if(grepl("C00",x)){
      new_id = str_remove(pattern="C00TU0",x)
      new_id = str_split(pattern="a",new_id)[[1]][1]
      new_id = paste0("TU",new_id)
      new_id = str_remove(pattern="v1",new_id)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else if(grepl("v",x)){
      new_id = str_remove(pattern="v1",x)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else {
      return(x)
  }
})

genos$donor.id = unlist(new_ids)

b_cells@meta.data = b_cells@meta.data %>% left_join(genos,by="donor.id")
rownames(b_cells@meta.data) = colnames(b_cells)
genotyped_cohort = subset(b_cells,subset = donor.id %in% genos$donor.id)
donors_to_keep = genotyped_cohort@meta.data %>% filter(!is.na(rs2331964_T))
genotyped_cohort = subset(genotyped_cohort,subset = donor.id %in% donors_to_keep$donor.id)

genotyped_cohort@meta.data$rs2331964_T = recode(genotyped_cohort@meta.data$rs2331964_T,
"0" = "C/C",
"1" = "T/C",
"2" = "T/T",
)

genotyped_cohort@meta.data %>% distinct(donor.id,.keep_all=TRUE) %>% dplyr::count(phenotype,rs2331964_T)
VlnPlot(subset(genotyped_cohort,subset = phenotype=="MS"),features="EAF2",group.by="rs2331964_T")
VlnPlot(subset(genotyped_cohort,subset = phenotype=="MS" & source=="CSF" & ann_monaco=="Exhausted B cells"),features="EAF2",group.by="rs2331964_T")


################################
# de eaf2
################################

DefaultAssay(genotyped_cohort) = "RNA"
genotyped_cohort@meta.data = genotyped_cohort@meta.data %>% mutate(eaf2_status = ifelse(rs2331964_HET == 1,"Het","Hom"))
genotyped_cohort@meta.data = genotyped_cohort@meta.data %>% mutate(donor_snp = paste0(donor.id,"_",eaf2_status))

for(source_of_interest in c("CSF","PBMC")){
for(cluster in c("Exhausted B cells","Plasmablasts")){
cells_for_de = subset(genotyped_cohort, subset = source == source_of_interest & ann_monaco == cluster)

# get rid of NAs
donors_to_keep = cells_for_de@meta.data %>% filter(!is.na(eaf2_status))
cells_for_de = subset(cells_for_de, subset = donor.id %in% donors_to_keep$donor.id)


# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 1

low_counts = cells_for_de@meta.data %>%
  group_by(donor.id,eaf2_status) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(donor.id,"_",eaf2_status))


# convert to sce object
cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

# aggregate counts
groups = colData(cells_for_de.sce)[, c("donor.id","eaf2_status")]
aggregated_counts  = aggregate.Matrix(t(counts(cells_for_de.sce)),
groupings = groups, fun = "sum") %>% t()

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]

group_vector = lapply(colnames(aggregated_counts),function(y){
      if(grepl("Het",y)){
        "Het"
      } else if(grepl("Hom",y)){
        "Hom"
      }
    }) %>% unlist %>% factor()



# make the DGE object
y=DGEList(aggregated_counts,group=group_vector,remove.zeros=TRUE)

# update sample info
y$samples =  y$samples %>%
mutate(donor_snp = rownames(y$samples)) %>%
left_join(genotyped_cohort@meta.data %>%
dplyr::select(Age,Gender,rs2331964_T,eaf2_status,donor_snp) %>%
distinct(donor_snp,.keep_all=TRUE),by="donor_snp")

design = model.matrix(~0+group_vector+Age+Gender,y$samples)
colnames(design) = c(levels(group_vector),"Age","Gender")

y = calcNormFactors(y)

df = data.frame(y$samples,"eaf2" = y$counts[rownames(y$counts)=="EAF2",]) %>% mutate(eaf2_norm = eaf2 * norm.factors)
df = df %>% filter(eaf2_norm != 0 )

k = kruskal.test(data=df,eaf2_norm ~ rs2331964_T)

p=ggplot(df,aes(rs2331964_T,eaf2_norm,fill=rs2331964_T))+geom_point()+geom_boxplot(alpha=0.5)+labs(y="EAF2 expression")+ggtitle(cluster)+theme_bw()+
annotate("label",label = paste0("P = ",round(k$p.value,2)),x=1.5,y=100)
png(paste0(cluster,"_",source_of_interest,"eaf2_csf_ms.png"),res=300,units="in",height=4,width=4)
print(p)
dev.off()
}
}

plotdat = subset(genotyped_cohort,subset = ann_monaco=="Plasmablasts" & source == "CSF")
DefaultAssay(plotdat="SCT")
plotdat = RunUMAP(plotdat,dims=1:50,reduction="harmony")
DefaultAssay(plotdat)="RNA"
FeaturePlot(plotdat,features="EAF2",split.by="rs2331964_T")



#######################################
# exhaustion markers
#######################################

biomarkers = FindMarkers(b_cells,recorrect_umi=FALSE,logfc.threshold=0.25,min.pct=0.25,only.pos=TRUE,ident.1="Exhausted B cells")
exhaustion_score = biomarkers %>% arrange(p_val_adj) %>% head(100)
exhaustion_score$gene = rownames(exhaustion_score)

features = exhaustion_score$gene
# module score
DefaultAssay(b_cells) = "RNA"
b_cells = AddModuleScore(b_cells,features = features,nbin=10,ctrl=5,name="Exhaustion")
b_cells[['exhaustion_score_z']] = RNOmni::RankNorm(b_cells@meta.data[[paste0("Exhaustion","1")]])
