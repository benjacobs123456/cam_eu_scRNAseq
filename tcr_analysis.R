#######################################
# Load packages
#######################################
library(tidyr)
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
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR")

# Read in data
t_cells = readRDS("t_cells_post_processing.rds")

rownames(t_cells@meta.data) = colnames(t_cells)

# define plot theme
theme_umap = function(){
    theme_minimal() %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
}
#######################################
# Big picture descriptive stats
#######################################

message("No. cells:")
message(nrow(t_cells@meta.data))

#######################################
# Clustering and reannotation
#######################################

# recluster
DefaultAssay(t_cells)="SCT"
t_cells = RunUMAP(t_cells,reduction="harmony",dims=1:50)
t_cells = FindNeighbors(t_cells)
t_cells = FindClusters(t_cells,resolution=1)

t_cell_names = c("Tcm/Naive helper T cells",
"Tem/Trm cytotoxic T cells",
"Tcm/Naive cytotoxic T cells",
"Regulatory T cells",
"MAIT cells",
"Tem/Effector helper T cells",
"Tem/Temra cytotoxic T cells",
"Cycling T cells",
"Early lymphoid/T lymphoid",
"Memory CD4+ cytotoxic T cells",
"CD8a/a",
"Follicular helper T cells",
"Type 1 helper T cells",
"Double-positive thymocytes",
"Treg(diff)",
"Type 17 helper T cells",
"Double-negative thymocytes",
"NKT cells",
"Trm cytotoxic T cells",
"CD8a/b(entry)",
"ILC",
"Tem/Effector helper T cells PD1+")

# filter to just celltypist T cell annotations
t_cells = subset(t_cells,subset = ann_celltypist_highres %in% t_cell_names)
t_cells = SetIdent(t_cells,value="ann_celltypist_highres")


t_cell_markers = c("CD3G","IL7R","CD8A","TCF7","FOXP3","CCL5","CCR7","MKI67","CD27","TRGV9","TRDV2")

# remove outliers (doublets)
outliers = subset(t_cells,UMAP_1 > 7)
t_cells = subset(t_cells,UMAP_1 < 7)

# recluster
DefaultAssay(t_cells)="SCT"
t_cells = RunUMAP(t_cells,reduction="harmony",dims=1:50)
t_cells = FindNeighbors(t_cells)
t_cells = FindClusters(t_cells,resolution=1)
t_cells = SetIdent(t_cells,value="ann_celltypist_highres")

##############################
# cluster biomarkers
##############################
#biomarkers = FindAllMarkers(t_cells,recorrect_umi=FALSE,logfc.threshold=0,min.pct=0,only.pos=TRUE)
#topbiomarkers = biomarkers %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=5)
#write_csv(topbiomarkers,"biomarkers.csv")
#write_csv(biomarkers,"all_biomarkers.csv")



#######################################
# DA & composition
#######################################
# restrict to abundant cells
abundant_cells = t_cells@meta.data %>% dplyr::count(ann_celltypist_highres) %>% filter(n>100)
t_cells = subset(t_cells, ann_celltypist_highres %in% abundant_cells$ann_celltypist_highres)

# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(iid,"_",source))

t_cells@meta.data$cell_type = t_cells@meta.data$ann_celltypist_highres
n_col = t_cells@meta.data$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

t_cells@meta.data$cell_type = factor(
  t_cells@meta.data$cell_type,
  ordered=TRUE,
  levels = c("Tcm/Naive helper T cells",
"Tem/Trm cytotoxic T cells",
"Regulatory T cells",
"Tem/Temra cytotoxic T cells",
"Tem/Effector helper T cells",
"Memory CD4+ cytotoxic T cells",
"MAIT cells",
"Tcm/Naive cytotoxic T cells",
"Type 1 helper T cells",
"Cycling T cells",
"Follicular helper T cells",
"Early lymphoid/T lymphoid"))


# plots to check annotations
p1=DimPlot(t_cells)
p2=FeaturePlot(t_cells,features=t_cell_markers)
p3=DotPlot(t_cells,features=t_cell_markers)

png("dimplot.png",res=600,units="in",width=6,height=6)
p1 + theme_umap() +  scale_color_manual(values = colour_pal)
dev.off()

png("featureplot.png",res=600,units="in",width=6,height=6)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=10,height=12)
p3
dev.off()

p=DimPlot(t_cells,label=F,raster=F,group.by="cell_type")+
  scale_color_manual(values = colour_pal)+
  ggtitle("")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="",y="")

png("dim_plot_simple_labels.png",res=600,units="in",width=5,height=4)
p
dev.off()

p=DimPlot(t_cells,split.by="source",label=F,raster=F,group.by="cell_type")+
  scale_color_manual(values = colour_pal)+
  ggtitle("")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="",y="")

png("dim_plot_simple_labels_source.png",res=600,units="in",width=5,height=4)
p
dev.off()

# stash sample info
sample_info = t_cells@meta.data %>%
  dplyr::select(iid,source,phenotype,donor_source) %>%
  distinct(donor_source,.keep_all=TRUE)

abundances = table(t_cells@meta.data$cell_type,t_cells@meta.data$donor_source)

# filter out clusters with 5 counts
abundances = abundances[rowSums(abundances)>5,]

sample_info = sample_info %>% filter(donor_source %in% colnames(abundances))
sample_info = sample_info[match(colnames(abundances),sample_info$donor_source),]
sample_info$grouping = paste0(sample_info$phenotype,"_",sample_info$source)
sample_info$grouping = sapply(sample_info$grouping, function(x){
  if(x == "NIND_CSF"){
    return("Control_CSF")
  } else if(x == "NIND_PBMC"){
    return("Control_PBMC")
  } else {
    return(x)
  }
})

y.ab = DGEList(abundances, samples=sample_info)

# update sample info
y.ab$samples =  y.ab$samples %>% left_join(t_cells@meta.data %>%
                                             distinct(donor_source,.keep_all=TRUE) %>%
                                             dplyr::select(Age,Sex,donor_source),by="donor_source")

# now loop
da_overall_list = list()
results_df = data.frame()
design <- model.matrix(~ 0 + grouping + Age + Sex,y.ab$sample)
colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Sex")

y.ab <- estimateDisp(y.ab, design,trend="none")
fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# define contrast for testing

contrast =  makeContrasts(
  MS_CSF - OIND_CSF,
  MS_CSF - Control_CSF,
  OIND_CSF - Control_CSF,
  MS_PBMC - OIND_PBMC,
  MS_PBMC - Control_PBMC,
  OIND_PBMC - Control_PBMC,
  MS_CSF - MS_PBMC,
  OIND_CSF - OIND_PBMC,
  OINDI_CSF - OINDI_PBMC,
  MS_CSF - OINDI_CSF,
  Control_CSF - Control_PBMC,
  levels = design
)

da_plots = list()
# do da tests
for(i in 1:length(colnames(contrast))){
  contrast_name = colnames(contrast)[i]
  res = glmQLFTest(fit, contrast = contrast[,i])

  # write to file
  print(summary(decideTests(res)))
  res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
  res$cell = rownames(res)
  write_csv(res,paste0("./crude_labels_edgeR_da_tests_",contrast_name,".csv"))
  res$significant = ifelse(res$P_adj<0.01,"yes","no")
  res$direction = ifelse(res$logFC>0,"Up","Down")
  comparison_label = contrast_name

  # refactor cell types (for plotting)
  res$cell = factor(res$cell,ordered=TRUE,
  levels = c("Tcm/Naive helper T cells",
"Tem/Trm cytotoxic T cells",
"Regulatory T cells",
"Tem/Temra cytotoxic T cells",
"Tem/Effector helper T cells",
"Memory CD4+ cytotoxic T cells",
"MAIT cells",
"Tcm/Naive cytotoxic T cells",
"Type 1 helper T cells",
"Cycling T cells",
"Follicular helper T cells",
"Early lymphoid/T lymphoid"))
  # da plot

  plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
    theme_classic()+
    scale_color_manual(values = colour_pal)+
    geom_text_repel(size=3,max.overlaps=100)+
    geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
    geom_vline(xintercept = 0,alpha=0.2)+
    NoLegend()+
    ggtitle(comparison_label)+
    scale_y_continuous(limits=c(0,50))+
    scale_x_continuous(limits=c(-9,9))+
    geom_point(shape=16,size=3)


  png(paste0("./crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
  print(plot)
  dev.off()
}

# DA summary plots
# find files
files = list.files(full.names = T,pattern = "crude_labels_edgeR_da_tests")

# read in
da_res = purrr::map(files,function(x){
  read_csv(x) %>%
    mutate(contrast = str_remove_all(str_remove_all(x,"./crude_labels_edgeR_da_tests_"),".csv"))

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

p = ggplot(da_res %>% mutate(logFC = ifelse(P_adj > 0.05, NA, logFC)),aes(cell,contrast,fill=logFC,label = ifelse(P_adj<0.05,"*","")))+
  geom_tile(color="black")+
  geom_text()+
  theme_minimal()+
  scale_fill_gradient2(low="purple",high="orange",midpoint = 0)+
  labs(x="Cell type",y = "Comparison",fill="Log fold\nchange")+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0))

png("da_summary_plot_all_cells_celltypist.png",res=600,units="in",width=10,height=4)
p
dev.off()


#######################################
# Phenotypes & DA
#######################################

# reorder phenotypes
t_cells@meta.data$phenotype = factor(t_cells@meta.data$phenotype,levels=c("NIND","OIND","OINDI","MS"),ordered=TRUE)

# proportion plots
png("./crude_labels_proportion_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(t_cells@meta.data,aes(phenotype,fill=cell_type))+
  geom_bar(position="fill",color="black")+
  facet_wrap(~source)+
  scale_fill_manual(values = colour_pal)+
  theme_classic()+
  labs(x="Phenotype",y="Proportion",fill="Cell type")
dev.off()

plots = list()
for(pheno in c("MS","OINDI","OIND","NIND")){
  p=ggplot(t_cells@meta.data %>%
             filter(phenotype==pheno),
           aes(donor.id,fill=cell_type))+
    geom_bar(position="fill",color="black")+
    facet_wrap(~source)+
    scale_fill_manual(values = colour_pal)+
    theme_classic()+
    labs(x="Donor",y="Proportion",fill="Cell type")+
    ggtitle(pheno)+
    theme(legend.position="none",axis.text.x=element_blank())
  plots[[length(plots)+1]] = p
}

png("./crude_labels_proportion_cell_counts_barplot_per_individual.png",res=300,width=12,height=4,units="in")
do.call("grid.arrange",plots)
dev.off()

#######################################
# DE
#######################################

DefaultAssay(t_cells)="RNA"

# change cell type names to get rid of / characters

# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(iid,"_",source))
t_cells@meta.data = t_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

# convert to sce object
t_cells.sce = as.SingleCellExperiment(t_cells)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = t_cells@meta.data %>%
  group_by(iid,source,phenotype,cell_type) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

# aggregate counts
groups = colData(t_cells.sce)[, c("ident", "phenotype","source","iid")]
aggregated_counts  = aggregate.Matrix(t(counts(t_cells.sce)),
                                      groupings = groups, fun = "sum")

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]


# get names for clusters
clusters = unique(t_cells@meta.data$ann_celltypist_highres)

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
    } else if(grepl("NIND_CSF",y)){
      "Control_CSF"
    } else if(grepl("NIND_PBMC",y)){
      "Control_PBMC"
    } else if(grepl("OIND_CSF",y)){
      "OIND_CSF"
    } else if(grepl("OIND_PBMC",y)){
      "OIND_PBMC"
    } else if(grepl("OINDI_CSF",y)){
      "OINDI_CSF"
    } else if(grepl("OINDI_PBMC",y)){
      "OINDI_PBMC"
    }
  }) %>% unlist %>% factor()

  # make the DGE object
  y=DGEList(de_input,group=group_vector,remove.zeros=T)

  # update sample info
  y$samples =  y$samples %>% mutate(full_cell_id = rownames(y$samples)) %>%
    left_join(t_cells@meta.data %>%
                dplyr::select(Age,Sex,full_cell_id) %>%
                distinct(full_cell_id,.keep_all=TRUE),by="full_cell_id")

                # make the design matrix
                message("Doing DE for ",cell_type)
                design = model.matrix(~0+group_vector+Age+Sex,y$samples)
                colnames(design) = c(levels(group_vector),"Age","Sex")


                keep = filterByExpr(
                  y,
                  design = design,
                  group = group_vector,
                  min.count = 10,
                  min.total.count = 500,
                  large.n = 100,
                  min.prop = 0.9)
                y = y[keep, , keep.lib.sizes=FALSE]
                y = calcNormFactors(y)
                y = estimateDisp(y,design,robust=TRUE)
                fit = glmQLFit(y, design, robust=TRUE)

                # define contrast for testing
                contrast = if(ncol(design)==10){
                makeContrasts(
                    MS_CSF - MS_PBMC,
                    OIND_CSF - OIND_PBMC,
                    Control_CSF - Control_PBMC,
                    MS_CSF - OINDI_CSF,
                    OINDI_CSF - OINDI_PBMC,
                    MS_PBMC - OINDI_PBMC,
                    MS_CSF - OIND_CSF,
                    MS_CSF - Control_CSF,
                    OIND_CSF - Control_CSF,
                    MS_PBMC - OIND_PBMC,
                    MS_PBMC - Control_PBMC,
                    OIND_PBMC - Control_PBMC,
                    levels = design
                  )} else if(
                    ncol(design)!=10 &
                    (!("Control_PBMC" %in% colnames(design)) |
                    !("Control_CSF" %in% colnames(design)) ) &
                    "OINDI_CSF" %in% colnames(design) & "OIND_PBMC" %in% colnames(design) & ("OIND_CSF" %in% colnames(design)))
                  {
                  makeContrasts(
                  MS_CSF - MS_PBMC,
                  OIND_CSF - OIND_PBMC,
                  MS_CSF - OINDI_CSF,
                  OINDI_CSF - OINDI_PBMC,
                  MS_PBMC - OINDI_PBMC,
                  MS_CSF - OIND_CSF,
                  MS_PBMC - OIND_PBMC,
                  levels = design
                  )} else if(
                    ncol(design)!=10 &
                    "Control_PBMC" %in% colnames(design) &
                    "Control_CSF" %in% colnames(design) &   "OINDI_CSF" %in% colnames(design) &
                    (!("OIND_PBMC" %in% colnames(design)) |
                    !("OIND_CSF" %in% colnames(design)) )
                    ){
                  makeContrasts(
                  MS_CSF - MS_PBMC,
                  Control_CSF - Control_PBMC,
                  MS_CSF - OINDI_CSF,
                  OINDI_CSF - OINDI_PBMC,
                  MS_PBMC - OINDI_PBMC,
                  MS_CSF - Control_CSF,
                  OIND_CSF - Control_CSF,
                  MS_PBMC - Control_PBMC,
                  levels = design
                  )}  else if(
                    ncol(design)!=10 &
                    (!("Control_PBMC" %in% colnames(design)) |
                    !("Control_CSF" %in% colnames(design)) |
                    !("OIND_PBMC" %in% colnames(design)) |
                    !("OIND_CSF" %in% colnames(design)) ) &
                    "OINDI_CSF" %in% colnames(design)
                   ){
                  makeContrasts(
                  MS_CSF - MS_PBMC,
                  MS_CSF - OINDI_CSF,
                  OINDI_CSF - OINDI_PBMC,
                  MS_PBMC - OINDI_PBMC,
                  levels = design
                  )} else if(
                    ncol(design)!=10 &
                    !("Control_PBMC" %in% colnames(design)) |
                    !("Control_CSF" %in% colnames(design)) |
                    !("OIND_PBMC" %in% colnames(design)) |
                    !("OIND_CSF" %in% colnames(design)) |
                    !("OINDI_CSF" %in% colnames(design))
                   ){
                  makeContrasts(
                  MS_CSF - MS_PBMC,
                  MS_PBMC - OINDI_PBMC,
                  levels = design
                  )}



  # do de tests
  for(i in 1:length(colnames(contrast))){
    contrast_name = colnames(contrast)[i]
    res = glmQLFTest(fit, contrast = contrast[,i])

    # write to file
    print(summary(decideTests(res)))
    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$gene = rownames(res)
    cell_type = str_replace_all(cell_type,"/","_")

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

    png(paste0("de_plot_",contrast_name,"_",cell_type,".png"),res=300,units="in",height=4,width=4)
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


# summary plots
# list files
files = list.files(pattern="edgeR_de",full.names=T)

# read in files
de_res = purrr::map(files,function(x){
  y = str_remove(str_remove(x,"./edgeR_de_tests_"),".csv")

  z = str_locate_all(y,"_")[[1]][3,1]

  cell = substr(y,z+1,nchar(y))
  y = substr(y,0,z-1)

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


  png(paste0("de_",comparison,".png"),res=600,units="in",height=6,width=6)
  print(gridExtra::grid.arrange(grobs = plots, top = comparison))
  dev.off()
  res_overall = do.call("bind_rows",res_overall)
  write_csv(res_overall,paste0("de_",comparison,".csv"))
}
make_plot("MS_CSF - OIND_CSF")
make_plot("MS_CSF - MS_PBMC")
make_plot("MS_CSF - OINDI_CSF")

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
  "MS_PBMC - OIND_PBMC",
  "MS_CSF - OINDI_CSF",
  "MS_PBMC - OINDI_PBMC")

res_overall = list()

clusters = t_cells$cell_type %>% unique

nperm = 10000
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
      } else if(comparison=="MS_CSF - OINDI_CSF"){
        "MS CSF - OINDI CSF"
      } else if(comparison=="MS_PBMC - OINDI_PBMC"){
        "MS PBMC vs OINDI PBMC"
      }



      cluster = levels(factor(clusters))[x]
      cluster = str_replace_all(cluster,"/","_")

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


      res = res %>% arrange(padj) %>%  dplyr::select(1,2,3,5,6) %>% mutate(cell_type = cluster) %>% mutate(comparison = this_comparison) %>% data.frame()
      res_overall[[(length(res_overall)+1)]] <<- res
      res = res %>% arrange(NES)
      res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
    }
  }
}

do_gsea("hallmark")

res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"res_overall.csv")
res_overall$fdr = p.adjust(res_overall$pval,method="fdr")
for(x in c( "MS CSF vs MS PBMC","MS CSF vs Control CSF","Control CSF vs Control PBMC","MS PBMC vs Control PBMC","MS PBMC vs OIND PBMC","MS CSF vs OIND CSF" )){


  plot_data = res_overall %>%
    filter(comparison == x) %>%
    filter(grepl("HALLMARK",pathway)) %>%
    arrange(fdr)
  plot_data$pathway = str_remove(plot_data$pathway,"HALLMARK_")
  top_pathways = plot_data %>% distinct(pathway,.keep_all = TRUE) %>% arrange(pval) %>% head(n=10)
  plot_data = plot_data %>% filter(pathway %in% top_pathways$pathway)
  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))

  png(paste0("summary_",x,".png"),res=300,units="in",width=7,height=4)
  p=ggplot(plot_data,aes(cell_type,pathway,fill=NES,label=ifelse(padj<0.01,ifelse(padj<0.001,ifelse(padj<0.0001,"***","**"),"*"),"")))+
    geom_tile(color="black")+
    geom_text()+
    scale_fill_gradient(low="purple",high="orange")+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)
  print(p)
  dev.off()
}

# run for other sets (go, kegg, reactome)
#res_overall = list()
#do_gsea("go")
#do_gsea("kegg")
#do_gsea("reactome")

#######################################
# Repertoire analysis CSF v periphery
#######################################
rownames(t_cells@meta.data) = colnames(t_cells)

# define clones
expanded_clones = t_cells@meta.data %>%
  group_by(iid,clone_id) %>%
  dplyr::count() %>%
  filter(n>1) %>%
  mutate(donor_clone = paste0(iid,"_",clone_id))

t_cells@meta.data = t_cells@meta.data %>%
  mutate(donor_clone = paste0(iid,"_",clone_id)) %>%
  mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))


# plot cell types in expanded vs not
png("./expanded_plot_celltypes_barplot.png",res=300,width=7,height=4,units="in")
plot_dat = t_cells@meta.data %>%
  mutate(phenotype_new = ifelse(as.character(phenotype) == "OINDI","ID",as.character(phenotype))) %>%
  mutate(expanded_clone = ifelse(expanded_clone == "Expanded","Y","N"))
plot_dat$phenotype_new = factor(plot_dat$phenotype_new,levels = c("NIND","OIND","ID","MS"),ordered=T)
ggplot(plot_dat,
  aes(expanded_clone,fill=cell_type))+
  geom_bar(position="fill",color="black")+
  facet_grid(source~phenotype_new)+
  scale_fill_manual(values = colour_pal)+
  theme_minimal()+
  labs(x="Expanded clone?",y="Proportion of\nB/plasma cells",fill="Cell type")
dev.off()



# define TRBV types
t_cells@meta.data$trbv_family = sapply(t_cells@meta.data$v_call_VDJ,function(x){
  y = str_split(x,"\\|",2)[[1]][1]
  y=str_split(y,pattern="-",n=2)[[1]][1]
  y = str_remove_all(y,"D")
  y = str_split(y,"\\/",2)[[1]][1]

  return(y)
})

t_cells@meta.data$trav_family = sapply(t_cells@meta.data$v_call_VJ,function(x){
  y = str_split(x,"\\|",2)[[1]][1]
  y=str_split(y,pattern="-",n=2)[[1]][1]
  y = str_remove_all(y,"D")
  y = str_split(y,"\\/",2)[[1]][1]

  return(y)
})


# simplify v call
all_calls = list()
for(i in c(1:length(t_cells@meta.data$v_call_VDJ))){
  message(i)
  this_call = t_cells@meta.data$v_call_VDJ[i]

  # take first call if ambiguous
  this_call = str_split(this_call,"\\|",2)[[1]][1]

  # split to simple call
  this_call = ifelse(str_count(this_call,"-")>1,
    paste0(str_split(this_call,"-")[[1]][1],"-",str_split(this_call,"-")[[1]][2]),
    this_call
  )
  this_call = str_remove_all(this_call,"D")
  this_call = str_remove_all(this_call,"/OR9-2")
  all_calls[[i]] = this_call
}
# add to metadata
t_cells@meta.data$v_call_simple = unlist(all_calls)


# repeat for tra
# simplify v call
all_calls = list()
for(i in c(1:length(t_cells@meta.data$v_call_VJ))){
  message(i)
  this_call = t_cells@meta.data$v_call_VJ[i]

  # take first call if ambiguous
  this_call = str_split(this_call,"\\|",2)[[1]][1]
  # get rid of pseudogenes
  this_call = str_split(this_call,"\\/",2)[[1]][1]

  # split to simple call
  this_call = ifelse(str_count(this_call,"-")>1,
    paste0(str_split(this_call,"-")[[1]][1],"-",str_split(this_call,"-")[[1]][2]),
    this_call
  )
  this_call = str_remove_all(this_call,"D")
  all_calls[[i]] = this_call
}

# add to metadata
t_cells@meta.data$trav_call_simple = unlist(all_calls)


## plots
# define plot theme
theme_umap = function(){
  theme_minimal() %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
}

# cell types
# compare csf vs periphery
png("csf_v_pbmc_celltypes.png",res=300,height=3,width=6,units="in")
DimPlot(t_cells,split.by="source")+scale_color_manual(values = colour_pal)+
  theme_umap()
dev.off()


## trbv
png("csf_v_pbmc_trbv.png",res=300,height=3,width=6,units="in")
DimPlot(t_cells,split.by="source",group.by="trbv_family")+
  theme_umap()

  dev.off()



#######################################
#  TRBV usage
#######################################

do_da_var = function(x, contrasts_to_test,variable,ylim=10, data = t_cells@meta.data,plot_title,var_label){

  # get rid of donors with low counts
  donors_to_keep = data %>% dplyr::count(donor_source) %>% filter(n>1)
  da_dat = data %>% filter(donor_source %in% donors_to_keep$donor_source)

  # stash sample info
  sample_info = da_dat %>%
    dplyr::select(iid,source,x,donor_source) %>%
    filter(!is.na(.data[[x]])) %>%
    distinct(donor_source,.keep_all=TRUE)

  data_for_abundances = da_dat %>%
    dplyr::select(iid,source,x,donor_source,.data[[variable]]) %>%
    filter(!is.na(.data[[x]]))
  abundances = table(data_for_abundances[[variable]],data_for_abundances$donor_source)

  # filter out clusters with <10 counts
  abundances = abundances[rowSums(abundances)>10,]

  sample_info = sample_info %>% filter(donor_source %in% colnames(abundances))
  sample_info = sample_info[match(colnames(abundances),sample_info$donor_source),]
  sample_info$grouping = paste0(sample_info[[x]],"_",sample_info$source)


  y.ab = DGEList(abundances, samples=sample_info)

  # update sample info
  y.ab$samples =  y.ab$samples %>% left_join(da_dat %>%
                                               filter(!is.na(.data[[x]])) %>%
                                               distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Sex,donor_source),by="donor_source")

  da_overall_list = list()
  results_df = data.frame()
  design <- model.matrix(~ 0 + grouping + Age + Sex,y.ab$sample)
  colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Sex")

  y.ab <- estimateDisp(y.ab, design,trend="none")
  fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

  # define contrast for testing

  contrast =  makeContrasts(
    contrasts = contrasts_to_test,
    levels = design
  )

  da_overall_res = list()
  # do da tests
  for(i in c(1:length(contrasts_to_test))){
    contrast_name = colnames(contrast)[i]
    res = glmQLFTest(fit, contrast = contrast[,i])

    # write to file
    print(summary(decideTests(res)))
    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$cell = rownames(res)
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res = res %>%
      mutate(direction = case_when(
        P_adj > 0.05 ~ "nonsig",
        P_adj <= 0.05 & logFC>0 ~ "Up",
        P_adj <= 0.05 & logFC<0 ~ "Down"
        ))

    comparison_label = contrast_name

    colours = c("Up" = "red", "Down" = "blue","nonsig"="grey")

    # da plot
    plot = ggplot(res,aes(logFC,-log10(PValue),label=cell, colour = direction))+
      theme_classic()+
      geom_point(shape=16,size=2)+
      geom_text_repel(data = res %>% filter(P_adj<0.05),
                      mapping = aes(logFC,-log10(PValue)),
                      size=3,max.overlaps=100,max.time = 10,max.iter = 100000)+
      geom_vline(xintercept = 0,alpha=0.2)+
      NoLegend()+
      ggtitle(comparison_label)+
      scale_y_continuous(limits=c(0,ylim))+
      scale_x_continuous(limits=c(-5,5))+
      scale_color_manual(values = colours)


    png(paste0(plot_title,"_pheno_comparisons_da_plot_",contrast_name,"_",variable,".png"),res=600,units="in",height=3,width=3)
    print(plot)
    dev.off()
    write_csv(res,
              file=paste0(plot_title,"_pheno_comparisons_da_plot_",contrast_name,"_",variable,".csv"))
    da_overall_res[[i]] = res %>% mutate(contrast = comparison_label)
  }
  da_overall_res = do.call("bind_rows",da_overall_res)

# make overall plot

# get FDR value
da_overall_res$P_adj = p.adjust(da_overall_res$PValue,method="fdr")

# plot
da_overall_res = da_overall_res %>% filter(contrast %in% c(
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

da_overall_res$contrast = factor(da_overall_res$contrast, levels = c(
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

da_overall_res = da_overall_res %>% arrange(cell)
da_overall_res$cell = factor(da_overall_res$cell, levels = unique(da_overall_res$cell) ,ordered=T)


p = ggplot(da_overall_res %>%
  mutate(logFC = ifelse(P_adj > 0.05, NA, logFC)),
  aes(cell,contrast,fill=logFC,label = ifelse(P_adj<0.05,"*","")))+
  geom_tile(color="black")+
  geom_text()+
  theme_minimal()+
  scale_fill_gradient2(low="purple",high="orange",midpoint = 0)+
  labs(x=var_label,y = "Comparison",fill="Log fold\nchange")+
  theme(axis.text.x=element_text(angle=90,vjust=0))

fileout = paste0("da_summary_plot_",plot_title,"_",variable,".png")
png(fileout,res=600,units="in",width=14,height=4)
print(p)
dev.off()

}

default_contrasts = c("MS_CSF - MS_PBMC",
                      "OIND_CSF - OIND_PBMC",
                      "OINDI_CSF - OINDI_PBMC",
                      "NIND_CSF - NIND_PBMC",
                      "MS_PBMC - OIND_PBMC",
                      "MS_PBMC - OINDI_PBMC",
                      "MS_PBMC - NIND_PBMC",
                      "MS_CSF - OIND_CSF",
                      "MS_CSF - OINDI_CSF",
                      "MS_CSF - NIND_CSF")


# run DA
do_da_var(x = "phenotype",
          contrasts_to_test = default_contrasts,
          variable="trbv_family",
          var_label = "TRBV family",
          plot_title = "all",
          ylim=20)

do_da_var(data = t_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
          ,
          x = "phenotype",
          contrasts_to_test = default_contrasts,
          variable="trbv_family",
          var_label = "TRBV family",
          plot_title = "eachclone",
          ylim=20)

do_da_var(data = t_cells@meta.data %>%
  mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>%
  distinct(donor_clone_source,.keep_all=TRUE) %>%
  filter(cell_type != "MAIT cells")
          ,
          x = "phenotype",
          contrasts_to_test = default_contrasts,
          variable="trbv_family",
          var_label = "TRBV family",
          plot_title = "eachclone_no_mait",
          ylim=20)

# individual gene segments
do_da_var(data = t_cells@meta.data,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                "MS_CSF - OINDI_CSF"),
          variable="v_call_simple",
          plot_title = "trbv_segments",
          var_label = "TRBV gene",
          ylim=20)

# individual gene segments
do_da_var(data = t_cells@meta.data %>%
            mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>%
            distinct(donor_clone_source,.keep_all=TRUE),
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                "MS_CSF - OINDI_CSF"),
          variable="v_call_simple",
          plot_title = "trbv_segments_each_clone",
          var_label = "TRBV gene",
          ylim=20)


# TRA  family
do_da_var(data = t_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
          ,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                "MS_CSF - OINDI_CSF"),
          variable="trav_family",
          plot_title = "trav_call_simple",
          var_label = "TRAV gene",
          ylim=20)


do_da_var(data = t_cells@meta.data  %>%
filter(grepl("^T",trav_call_simple)) %>%
mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>%
distinct(donor_clone_source,.keep_all=TRUE)
          ,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                  "MS_CSF - OINDI_CSF"),
          variable="trav_call_simple",
          plot_title = "trav_eachclone",
          var_label = "TRAV gene",
          ylim=20)


#######################################
# Clonal analysis
#######################################

simplify_pval = function(x){
  ifelse(
    is.na(x),
    "NA",
    ifelse(x<0.0001,
           "***",
           ifelse(x<0.001,
                  "**",
                  ifelse(x<0.01,
                         "*",
                         as.character(round(x,2))
                  ))))
}

# quantify clonality of csf
binary_comparison = function(variable,plot_title,level="Yes"){
  props = t_cells@meta.data %>%
    group_by(iid,source,phenotype) %>%
    dplyr::count(.data[[variable]],.drop=F) %>%
    mutate(total = sum(n), prop = n/total) %>%
    filter(!is.na(.data[[variable]])) %>%
    pivot_wider(id_cols = c(iid,source,phenotype,total),
    values_from = c(n,prop),
    names_from=variable) %>%
    mutate(new_n = ifelse(
      is.na(.data[[paste0("n_",level)]]),
      0,
      .data[[paste0("n_",level)]]
      )) %>%
    mutate(prop = new_n / total)


    # compare MS with NIND
    a = props[props$source=="CSF" & props$phenotype=="NIND",][['prop']]
    b = props[props$source=="CSF" & props$phenotype=="MS",][['prop']]
    message("MS CSF vs NIND CSF")
    message(round(median(a),2)," ",round(median(b),2))
    message(round(t.test(a,b)$p.value,2))

    a = props[props$source=="CSF" & props$phenotype=="OIND",][['prop']]
    b = props[props$source=="CSF" & props$phenotype=="MS",][['prop']]
    message("MS CSF vs OINDI CSF")
    message(round(median(a),2)," ",round(median(b),2))
    message(round(t.test(a,b)$p.value,2))

    a = props[props$source=="CSF" & props$phenotype=="OINDI",][['prop']]
    b = props[props$source=="CSF" & props$phenotype=="MS",][['prop']]
    message("MS CSF vs OINDI CSF")
    message(round(median(a),2)," ",round(median(b),2))
    message(round(t.test(a,b)$p.value,2))


  phenos = unique(props$phenotype)
  pvals = list()
  for(i in c(1:length(phenos))){
    a = props[props$source=="CSF" & props$phenotype==phenos[i],][['prop']]
    b = props[props$source=="PBMC" & props$phenotype==phenos[i],][['prop']]

    p=if(length(a)<=1 | length(b) <=1){
      NA
    } else {
      t.test(a,b)$p.value
    }
    message(p)
    pvals[[length(pvals)+1]] = p
  }
  pvals = data.frame(phenotype = phenos,P = unlist(pvals))

  props = props %>%
    left_join(pvals,by="phenotype") %>%
    mutate(P = simplify_pval(P))

  p=ggplot(props,aes(phenotype,prop,fill=source))+
    geom_boxplot(color="black")+
    scale_fill_brewer(palette="Set1")+
    ggtitle(plot_title)+
    geom_text(mapping = aes(x = phenotype,y = 1.05,label = P))+
    theme_minimal()+
    labs(x="Phenotype",y="Proportion")


  png(file=paste0(plot_title,"_comparison.png"),res=300,units="in",width=4,height=4)
  print(p)
  dev.off()
}

# big picture numbers

# define clones
expanded_clones = t_cells@meta.data %>%
  group_by(iid,clone_id) %>%
  dplyr::count() %>%
  filter(n>1) %>%
  mutate(donor_clone = paste0(iid,"_",clone_id))

t_cells@meta.data = t_cells@meta.data %>%
  mutate(donor_clone = paste0(iid,"_",clone_id)) %>%
  mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

binary_comparison(variable = "expanded_clone",plot_title="Clonal expansion",level="Expanded")

t_cells@meta.data %>% nrow
table(t_cells@meta.data$expanded_clone)
t_cells@meta.data %>% distinct(donor_clone) %>% nrow
t_cells@meta.data %>% dplyr::count(expanded_clone) %>% mutate (n/sum(n))
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n>1)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n==2)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n>5)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% filter(n==1)%>% nrow
t_cells@meta.data %>% dplyr::count(donor_clone) %>% arrange(desc(n)) %>% head

# recalculate clonal size per donor
clonal_size = t_cells@meta.data %>% dplyr::count(donor_clone) %>% dplyr::rename("clonal_size" = "n")
t_cells@meta.data = t_cells@meta.data %>% left_join(clonal_size,by="donor_clone")

# donors with expanded clones
donors_with_expanded_clones = t_cells@meta.data  %>%
  group_by(iid) %>%
  dplyr::count(donor_clone) %>%
  filter(n>1) %>%
  distinct(iid)
t_cells@meta.data = t_cells@meta.data %>%
  mutate(donor_with_expanded_clones = ifelse(iid %in% donors_with_expanded_clones$iid,"Yes","No"))

donors_with_expanded_clones_csf = t_cells@meta.data  %>%
  filter(source=="CSF") %>%
  group_by(iid) %>%
  dplyr::count(donor_clone) %>%
  filter(n>1) %>%
  distinct(iid)
t_cells@meta.data = t_cells@meta.data %>%
  mutate(donor_with_expanded_clones_csf = ifelse(iid %in% donors_with_expanded_clones_csf$iid,"Yes","No"))

# get some counts
t_cells@meta.data %>%
  group_by(phenotype) %>%
  distinct(iid,.keep_all=T) %>%
  dplyr::count(donor_with_expanded_clones) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(donor_with_expanded_clones=="Yes")

# histogram
p=ggplot(t_cells@meta.data,
         aes(clonal_size,fill=phenotype))+
  geom_histogram()+
  theme_bw()+
  scale_fill_brewer(palette="Set2")+
  labs(x="Clone size")
png("clonal_histogram.png",res=300,height=2,width=4,units="in")
p
dev.off()



#######################################
# look at relation to phenotypes
#######################################

cells = unique(t_cells@meta.data$ann_celltypist_highres)

dat = t_cells@meta.data %>%
  filter(cell_type %in% cells) %>%
  group_by(cell_type,iid,source,phenotype,Age,Sex,OCB,ms_subtype) %>%
  dplyr::count(expanded_clone) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(total>1 & expanded_clone=="Expanded")

dat2 = t_cells@meta.data %>%
  group_by(iid,source,phenotype,Age,Sex,OCB,ms_subtype) %>%
  dplyr::count(expanded_clone) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(total>1 & expanded_clone=="Expanded")

dat %>% group_by(phenotype,source) %>% summarise(median(prop),range(prop))
dat %>% filter(phenotype=="MS") %>% group_by(OCB,source) %>% summarise(median(total),range(total))

p0 = ggplot(
  dat,
  aes(phenotype,
      prop,
      fill=source)
)+
  facet_wrap(~cell_type)+
  geom_boxplot(alpha=0.8)+
  geom_point(alpha=0.7,position=position_dodge(width=0.75))+
  theme_minimal()+
  scale_fill_brewer(palette="Set1")+
  labs(x="Phenotype",y="Proportion of \nexpanded B cells",fill="Compartment")
png("clonal_proportions_individual_level.png",res=300,height=4,width=6,units="in")
p0
dev.off()

# plot per individual
rownames(t_cells@meta.data) = colnames(t_cells)
p=DimPlot(subset(t_cells,source=="CSF" & phenotype=="MS"),group.by="expanded_clone",split.by="iid",ncol=9)+
  theme_umap()
png("./indiv_clonal_dim_plots_ms.png",res=300,units="in",height=8,width=8)
p
dev.off()


# recode ms subtype
dat = dat %>%
  mutate(ms_subtype_binary = ifelse(ms_subtype == "SPMS","RMS",ms_subtype))
dat2 = dat2 %>%
  mutate(ms_subtype_binary = ifelse(ms_subtype == "SPMS","RMS",ms_subtype))

# counts
pheno_corr = function(variable,plot_title){
  plot_dat = dat2 %>% filter(phenotype=="MS" & !is.na(.data[[variable]]))

  counts = plot_dat %>% filter(source=="CSF") %>% distinct(iid,.keep_all=T) %>% ungroup %>% dplyr::count(.data[[variable]])
  print(counts)

  pvals = list()
  phenos = unique(plot_dat$source)

  for(i in c(1:length(phenos))){
    a = plot_dat[plot_dat$source==phenos[i] & plot_dat[[variable]]==levels(factor(plot_dat[[variable]]))[1],][['prop']]
    b = plot_dat[plot_dat$source==phenos[i] & plot_dat[[variable]]==levels(factor(plot_dat[[variable]]))[2],][['prop']]

    p=if(length(a)<=1 | length(b) <=1){
      NA
    } else {
      t.test(a,b)$p.value
    }
    message(p)
    pvals[[length(pvals)+1]] = p
  }
  pvals = data.frame(source = phenos,P = unlist(pvals))

  # combine
  plot_dat = plot_dat %>%
    left_join(pvals,by="source") %>%
    mutate(P = simplify_pval(P))

  p1=ggplot(
    plot_dat,
    aes(.data[[variable]],
        prop))+
    geom_boxplot(alpha=0.8)+
    theme_minimal()+
    facet_wrap(~source)+
    scale_fill_brewer(palette="Set3")+
    labs(x=plot_title,y="Proportion of \nexpanded B cells")+
    geom_text(mapping = aes(.data[[variable]],y=1.05,label = P))

  png(paste0("clonal_proportions_vs_",plot_title,".png"),res=300,height=4,width=8,units="in")
  print(p1)
  dev.off()
}

pheno_corr_celltype = function(variable,plot_title){
  plot_dat = dat %>% filter(phenotype=="MS" & !is.na(.data[[variable]]))


  pvals = list()
  plot_dat = plot_dat %>% mutate(cell_source = paste0(cell_type,"_",source))
  phenos = unique(plot_dat$cell_source)

  for(i in c(1:length(phenos))){
    a = plot_dat[plot_dat$cell_source==phenos[i] & plot_dat[[variable]]==levels(factor(plot_dat[[variable]]))[1],][['prop']]
    b = plot_dat[plot_dat$cell_source==phenos[i] & plot_dat[[variable]]==levels(factor(plot_dat[[variable]]))[2],][['prop']]

    p=if(length(a)<=1 | length(b) <=1){
      NA
    } else {
      t.test(a,b)$p.value
    }
    message(p)
    pvals[[length(pvals)+1]] = p
  }
  pvals = data.frame(cell_source = phenos,P = unlist(pvals))

  # combine
  plot_dat = plot_dat %>%
    left_join(pvals,by="cell_source") %>%
    mutate(P = simplify_pval(P))

  p1=ggplot(
    plot_dat,
    aes(cell_type,
        prop,
        fill=.data[[variable]]))+
    geom_boxplot(alpha=0.8)+
    theme_minimal()+
    facet_wrap(~source)+
    scale_fill_brewer(palette="Set3")+
    labs(x="Cell type",y="Proportion of \nexpanded B cells",fill=plot_title)+
    geom_text(mapping = aes(cell_type,y=1.05,label = P))

  png(paste0("clonal_proportions_vs_",plot_title,".png"),res=300,height=4,width=8,units="in")
  print(p1)
  dev.off()
}
pheno_corr("OCB","OCB status")
pheno_corr("ms_subtype","MS subtype")


#######################################
# Clonal phenotypes
#######################################

rownames(t_cells@meta.data) = colnames(t_cells)
# compare expanded vs non-expanded in MS
expanded_v_not_plot = function(x){
  DimPlot(subset(t_cells,phenotype=="MS"),
          split.by="expanded_clone",group.by=x)+
    theme_umap()
}

rownames(t_cells@meta.data) = colnames(t_cells)
p0=expanded_v_not_plot("cell_type")+ggtitle("Cell type")
p1=expanded_v_not_plot("trbv_family")+ggtitle("TRBV family")
p2=expanded_v_not_plot("v_call_simple")+ggtitle("TRBV")
p3=expanded_v_not_plot("trav_call_simple")+ggtitle("TRAV")
png("expanded_v_not_tcrs_just_ms.png",res=300,units="in",height=16,width=16)
grid.arrange(p0,p1,p2,p3)
dev.off()



#######################################
# Clonal analysis
#######################################


# define DA function for clonal vs non-clonal
do_da_var_clonal = function(x,variable,ylim=10, data = t_cells@meta.data,plot_title){

  data = data %>% mutate(donor_expanded = paste0(iid,"_",expanded_clone))

  # get rid of donors with low counts
  donors_to_keep = data %>% dplyr::count(donor_expanded) %>% filter(n>1)
  da_dat = data %>% filter(donor_expanded %in% donors_to_keep$donor_expanded)

  # stash sample info
  sample_info = da_dat %>%
    dplyr::select(iid,expanded_clone,x,donor_expanded) %>%
    filter(!is.na(.data[[x]])) %>%
    distinct(donor_expanded,.keep_all=TRUE)

  data_for_abundances = da_dat %>%
    dplyr::select(iid,source,x,donor_expanded,.data[[variable]]) %>%
    filter(!is.na(.data[[x]]))
  abundances = table(data_for_abundances[[variable]],data_for_abundances$donor_expanded)

  # filter out clusters with <10 counts
  abundances = abundances[rowSums(abundances)>10,]

  sample_info = sample_info %>% filter(donor_expanded %in% colnames(abundances))
  sample_info = sample_info[match(colnames(abundances),sample_info$donor_expanded),]
  sample_info$grouping = paste0(sample_info[[x]],"_",sample_info$expanded_clone)
  sample_info = sample_info %>%
    separate(grouping,sep="_",into=c("var","grouping")) %>%
    dplyr::select(-var) %>%
    mutate(grouping = str_replace_all(grouping," ","_"))

  y.ab = DGEList(abundances, samples=sample_info)

  # update sample info
  y.ab$samples =  y.ab$samples %>% left_join(da_dat %>%
                                               filter(!is.na(.data[[x]])) %>%
                                               distinct(donor_expanded,.keep_all=TRUE) %>% dplyr::select(Age,Sex,donor_expanded),by="donor_expanded")

  da_overall_list = list()
  results_df = data.frame()
  design <- model.matrix(~ 0 + grouping + Age + Sex,y.ab$sample)
  colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Sex")

  y.ab <- estimateDisp(y.ab, design,trend="none")
  fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

  # define contrast for testing

  contrast =  makeContrasts(
    contrasts = "Expanded - Not_expanded",
    levels = design
  )

  res = glmQLFTest(fit, contrast = contrast)

    # write to file
    print(summary(decideTests(res)))
    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$cell = rownames(res)
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res = res %>%
      mutate(direction = case_when(
        P_adj > 0.05 ~ "nonsig",
        P_adj <= 0.05 & logFC>0 ~ "Up",
        P_adj <= 0.05 & logFC<0 ~ "Down"
        ))

    comparison_label = "Expanded vs Not expanded"

    colours = c("Up" = "red", "Down" = "blue","nonsig"="grey")

    # da plot
    plot = ggplot(res,aes(logFC,-log10(PValue),label=cell, colour = direction))+
      theme_classic()+
      geom_point(shape=16,size=2)+
      geom_text_repel(data = res %>% filter(P_adj<0.05),
                      mapping = aes(logFC,-log10(PValue)),
                      size=3,max.overlaps=100,max.time = 10,max.iter = 100000)+
      geom_vline(xintercept = 0,alpha=0.2)+
      NoLegend()+
      ggtitle(comparison_label)+
      scale_y_continuous(limits=c(0,ylim))+
      scale_x_continuous(limits=c(-5,5))+
      scale_color_manual(values = colours)

      write_csv(res,
                file=paste0(plot_title,"_pheno_comparisons_da_plot_",contrast_name,"_",variable,".csv"))

    png(paste0(plot_title,"_pheno_comparisons_da_plot_expanded_vs_not_expanded_",variable,".png"),res=600,units="in",height=3,width=3)
    print(plot)
    dev.off()
}

do_da_var_clonal(x="cell_type",
variable="cell_type",
data = t_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_celltypes_ms",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = t_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_trbv_ms",
ylim=50)

do_da_var_clonal(x="trav_call_simple",
variable="trav_call_simple",
data = t_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_trav_ms",
ylim=50)


do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = t_cells@meta.data %>% filter(phenotype=="MS") %>% filter(cell_type != "MAIT cells"),
plot_title="clonal_trbv_ms_no_mait",
ylim=50)

do_da_var_clonal(x="trav_call_simple",
variable="trav_call_simple",
data = t_cells@meta.data %>% filter(phenotype=="MS") %>% filter(cell_type != "MAIT cells"),
plot_title="clonal_trav_ms",
ylim=50)



################################
# de clonal vs not
################################

t_cells@meta.data = t_cells@meta.data %>% mutate(donor_expanded = paste0(iid,"_",expanded_clone))

do_de = function(dat,plot_title,mincells=2){
  DefaultAssay(dat) = "RNA"
  cells_for_de = dat

  # print out params
  message("Cells:",cells_for_de@meta.data %>% nrow)
  message("Expanded cells:",cells_for_de@meta.data %>% filter(expanded_clone=="Expanded") %>% nrow)
  message("Non-expanded cells:",cells_for_de@meta.data %>% filter(expanded_clone=="Not expanded") %>% nrow)
  message("Donors expanded:",cells_for_de@meta.data %>% filter(expanded_clone=="Expanded") %>% distinct(iid) %>% nrow)
  message("Donors non-expanded:",cells_for_de@meta.data %>% filter(expanded_clone=="Not expanded") %>% distinct(iid) %>% nrow)

  rownames(cells_for_de@meta.data) = colnames(cells_for_de)

  # tabulate to find out which 'groups' will have insufficient cells for DE
  min_cells_per_sample = mincells

  low_counts = cells_for_de@meta.data %>%
    group_by(iid,expanded_clone) %>%
    dplyr::count() %>%
    arrange(n) %>%
    filter(n<min_cells_per_sample) %>%
    mutate(donor_to_exclude = paste0(iid,"_",expanded_clone))

  # convert to sce object
  cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

  # aggregate counts
  groups = colData(cells_for_de.sce)[, c("iid","expanded_clone")]
  aggregated_counts  = aggregate.Matrix(t(counts(cells_for_de.sce)),
                                        groupings = groups, fun = "sum") %>% t()

  # remove groups with low cell counts for DE (<n cells)
  aggregated_counts = aggregated_counts[!colnames(aggregated_counts) %in% low_counts$donor_to_exclude,]

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
                dplyr::select(Age,Sex,donor_expanded) %>%
                distinct(donor_expanded,.keep_all=TRUE),by="donor_expanded")

  design = model.matrix(~0+group_vector+Age+Sex,y$samples)
  colnames(design) = c(levels(group_vector),"Age","Sex")

  keep = filterByExpr(
    y,
    design = design,
    group = group_vector,
    min.count = 10,
    min.total.count = 500,
    large.n = 100,
    min.prop = 0.9)

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
    scale_color_manual(values = colours)+
    theme(panel.grid = element_blank())

  png(paste0(plot_title,"_expanded_vs_not.png"),res=300,units="in",height=4,width=4)
  print(plot)
  dev.off()

  # gsea
  de = res
  ranked_genes = de %>% arrange(desc(logFC)) %>% dplyr::select(gene,logFC)
  ranked_genes_vector = ranked_genes$logFC
  names(ranked_genes_vector) = ranked_genes$gene
  message("There are ",nrow(ranked_genes)," genes in this analysis")
  write_csv(res,file=paste0(plot_title,"_de_expanded_vs_not.csv"))

  # do gsea
  nperm=10000
  geneset="hallmark"
  genelist = eval(parse(text = paste0(geneset,"_list")))
  res = fgsea(genelist, stats = ranked_genes_vector, minSize=10,eps=0, nPermSimple = nperm)
  res = res %>% arrange(NES)
  res$fdr = p.adjust(res$pval,method="fdr")
  res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
  topres = res %>% arrange(pval) %>% head(n=20) %>% arrange(desc(NES))
  topres$pathway = str_remove(topres$pathway,"HALLMARK_")
  topres$pathway = factor(topres$pathway,ordered=TRUE,levels=topres$pathway)
  topres$direction = ifelse(topres$NES>0,"Up","Down")
  write_csv(res,file=paste0(plot_title,"_gsea_expanded_vs_not.csv"))

  colours = c("Up" = "red", "Down" = "blue")
  # NB labeled with nominal significance
  plot = ggplot(topres,aes(NES,pathway,label=ifelse(fdr<0.01,"*"," "),fill=direction))+
    theme_bw()+
    labs(x="Normalised enrichment score",fill="Pathway enrichment")+
    geom_col(color="black")+
    geom_text(size=5,position = position_dodge(0.5))+
    guides(alpha=FALSE)+
    theme_minimal()

  png(paste0(plot_title,"_gsea_expanded_vs_not.png"),res=300,units="in",height=3,width=6)
  print(plot)
  dev.off()
}


t_cell_types = c("Tcm/Naive helper T cells",
"Tem/Trm cytotoxic T cells",
"Regulatory T cells",
"Tem/Temra cytotoxic T cells",
"Tem/Effector helper T cells",
"Memory CD4+ cytotoxic T cells",
"MAIT cells",
"Tcm/Naive cytotoxic T cells",
"Type 1 helper T cells",
"Cycling T cells",
"Follicular helper T cells",
"Early lymphoid/T lymphoid")

rownames(t_cells@meta.data) = colnames(t_cells)
for(cell in t_cell_types){
  new_name = str_replace_all(cell,"\\/","_")
  do_de(
    dat = subset(t_cells,phenotype=="MS" & source == "CSF" & ann_celltypist_highres == cell),
    plot_title = paste0("MS_CSF_",new_name)
  )

}


# cdr3 seqs
cdr3_sequences = vdjdb$cdr3

logo_plot = function(cdr3_sequences, length=15){

# filter to specified length
cdr3_sequences = cdr3_sequences[nchar(cdr3_sequences)==length]

# split string
plot_df = list()
for(i in c(1:length)){
  plot_df[[i]] = substr(cdr3_sequences,i,i)
}
# combine
plot_df = do.call("cbind",plot_df)

# make df
plot_df = data.frame(plot_df)
colnames(plot_df) = c(1:length)

# lengthen
plot_df = tidyr::pivot_longer(plot_df,cols = c(1:length))


# counts & proportions
plot_df = plot_df %>%
  dplyr::group_by(name) %>%
  dplyr::count(value) %>%
  dplyr::mutate(prop = n/sum(n))

# get top letter for each pos
plot_df = plot_df %>%
dplyr::group_by(name) %>%
dplyr::slice_max(prop,with_ties=F) %>%
dplyr::arrange(name) %>%
dplyr::mutate(name = as.numeric(name))

# summarise
ggplot(plot_df,aes(name,prop,fill=value,label=value))+
geom_col(color="black")+
geom_text(vjust=0)+
theme_minimal()+
labs(x="Position",y="Proportion")
}


#################################
# check for public clones
################################

library(circlize)

public_clones = t_cells@meta.data %>%
  group_by(clone_id) %>%
  dplyr::count(iid)  %>%
  dplyr::count(clone_id) %>%
  filter(n>1)

public_clones_with_tra = t_cells@meta.data %>%
filter(status == "TRB + TRA") %>%
  group_by(clone_id) %>%
  dplyr::count(iid)  %>%
  dplyr::count(clone_id) %>%
  filter(n>1)

overall_dat = data.frame()
ms_specific_clones = list()
for(i in c(1:nrow(public_clones))){
  message(i)
    dat = t_cells@meta.data %>%
    mutate(iid = paste0(iid,"_",phenotype)) %>%
    filter(clone_id %in% public_clones$clone_id[i])
  ms_specific =  data.frame(public_clones$clone_id[i],all(dat$phenotype=="MS"))
  ms_specific_clones[[i]] = ms_specific
  dat = dat %>%
    dplyr::count(iid,clone_id)
  dat = expand.grid(dat$iid,dat$iid) %>%
    filter(Var1!=Var2)
  colnames(dat)= c("from","to")
  overall_dat <<- bind_rows(overall_dat,dat)
}
overall_dat$value = 1
overall_dat$clone_id = NULL
ms_specific_clones = do.call("bind_rows",ms_specific_clones)

# ms-specific
colnames(ms_specific_clones) = c("clone_id","ms_specific")
ms_specific_clones %>% dplyr::count(ms_specific) %>%
mutate(prop = n/sum(n))

ms_specific_clones = ms_specific_clones %>% filter(ms_specific==T)
grid.cols = t_cells@meta.data %>%
  distinct(iid,phenotype) %>%
  mutate(col = case_when(
    phenotype == "MS" ~ "purple",
    phenotype=="NIND" ~ "blue",
    phenotype=="OIND"~ "red",
    phenotype=="OINDI"~"orange"
  ))
grid_col_vec = grid.cols$col
names(grid_col_vec) = paste0(grid.cols$iid,"_",grid.cols$phenotype)

set.seed(123)

png("chord_diag.png",res=600,units="in",width=4,height=4)
chordDiagram(overall_dat,grid.col = grid_col_vec,annotationTrack = c("name","grid"))
dev.off()


############################################
# lookup in database
############################################

# save metadate (so can be read back in quickly)
# write_csv(t_cells@meta.data,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/metadata.csv")

t_cells@meta.data = t_cells@meta.data %>%
mutate(ms_specific = ifelse(clone_id %in% ms_specific_clones$clone_id,"MS specific","Not MS specific"))
rownames(t_cells@meta.data) = colnames(t_cells)
# DimPlot(subset(t_cells,phenotype=="MS"),split.by="ms_specific")

vdjdb = immunarch::dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb",.species = "HomoSapiens")
# mcpas = immunarch::dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human")

# define combination of trbv & cdr3
vdjdb = vdjdb %>%
  filter(gene == "TRB") %>%
  mutate(trbv_cdr3 = paste0(`v.segm`,"_",cdr3))

# define pathogen-specific t cells
dat = t_cells@meta.data %>%
  mutate(trbv_cdr3 = paste0(v_call_VDJ,"_",junction_aa_VDJ)) %>%
  left_join(vdjdb,by="trbv_cdr3")

# count
denominators = dat  %>%
  distinct(iid,source,phenotype) %>%
  dplyr::count(source,phenotype)

ebv = dat  %>%
  group_by(iid,source,phenotype) %>%
  dplyr::count(antigen.species) %>%
  mutate(total_cells = sum(n)) %>%
  filter(antigen.species == "EBV") %>%
  mutate(prop_ebv_specific = n/total_cells) %>%
  mutate(any_ebv_specific = ifelse(prop_ebv_specific>0,"yes","no")) %>%
  group_by(source,phenotype) %>%
  dplyr::count(any_ebv_specific) %>%
  left_join(denominators,by=c("source","phenotype")) %>%
  mutate(pct = n.x/n.y*100) %>%
  mutate("EBV" = paste0(n.x," / ",n.y," (",round(pct,1),"%)")) %>%
  dplyr::select(1,2,7)

cmv =   dat %>%
    group_by(iid,source,phenotype) %>%
    dplyr::count(antigen.species) %>%
    mutate(total_cells = sum(n)) %>%
    filter(antigen.species == "CMV") %>%
    mutate(prop_cmv_specific = n/total_cells) %>%
    mutate(any_cmv_specific = ifelse(prop_cmv_specific>0,"yes","no")) %>%
    group_by(source,phenotype) %>%
    dplyr::count(any_cmv_specific) %>%
    left_join(denominators,by=c("source","phenotype")) %>%
    mutate(pct = n.x/n.y*100) %>%
    mutate("CMV" = paste0(n.x," / ",n.y," (",round(pct,1),"%)")) %>%
    dplyr::select(1,2,7)



viral_dat = ebv %>%
  left_join(cmv,by=c("source","phenotype")) %>%
  pivot_wider(id_cols = phenotype,names_from = source,values_from = c("EBV","CMV"))
write_csv(viral_dat,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/viral_specificity.csv")


# DE EBV-specific in MS vs other

# subset to ebv-specific cells
ebv_specific = vdjdb %>%
  filter(gene == "TRB" & antigen.species=="EBV") %>%
  mutate(trbv_cdr3 = paste0(`v.segm`,"_",cdr3))
cmv_specific = vdjdb %>%
  filter(gene == "TRB" & antigen.species=="CMV") %>%
  mutate(trbv_cdr3 = paste0(`v.segm`,"_",cdr3))

# define pathogen-specific t cells
t_cells@meta.data = t_cells@meta.data %>%
  mutate(trbv_cdr3 = paste0(v_call_VDJ,"_",junction_aa_VDJ))

t_cells@meta.data$ebv_specific = ifelse(t_cells@meta.data$trbv_cdr3 %in% ebv_specific$trbv_cdr3,"yes","no")
t_cells@meta.data$cmv_specific = ifelse(t_cells@meta.data$trbv_cdr3 %in% cmv_specific$trbv_cdr3,"yes","no")


do_de_viral = function(dat,plot_title,mincells=2){
  DefaultAssay(dat) = "RNA"
  cells_for_de = dat

  rownames(cells_for_de@meta.data) = colnames(cells_for_de)
  print(cells_for_de@meta.data %>% dplyr::count(virus_specific))
  # tabulate to find out which 'groups' will have insufficient cells for DE
  min_cells_per_sample = mincells

  low_counts = cells_for_de@meta.data %>%
    group_by(iid,virus_specific) %>%
    dplyr::count() %>%
    arrange(n) %>%
    filter(n<min_cells_per_sample) %>%
    mutate(donor_to_exclude = paste0(iid,"_",virus_specific))

  # convert to sce object
  cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

  # aggregate counts
  groups = colData(cells_for_de.sce)[, c("iid","virus_specific")]
  aggregated_counts  = aggregate.Matrix(t(counts(cells_for_de.sce)),
                                        groupings = groups, fun = "sum") %>% t()

  # remove groups with low cell counts for DE (<n cells)
  aggregated_counts = aggregated_counts[!colnames(aggregated_counts) %in% low_counts$donor_to_exclude,]

  group_vector = lapply(colnames(aggregated_counts),function(y){
    if(grepl("no",y)){
      "Not_Virus_specific"
    } else if(grepl("yes",y)){
      "Virus_specific"
    }
  }) %>% unlist %>% factor()

  # make the DGE object
  y=DGEList(aggregated_counts,group=group_vector,remove.zeros=TRUE)

  # update sample info
  y$samples =  y$samples %>%
    mutate(donor_virus = rownames(y$samples)) %>%
    left_join(t_cells@meta.data %>%
                dplyr::select(Age,Sex,donor_virus) %>%
                distinct(donor_virus,.keep_all=TRUE),by="donor_virus")

  design = model.matrix(~0+group_vector+Age+Sex,y$samples)
  colnames(design) = c(levels(group_vector),"Age","Sex")

  keep = filterByExpr(
    y,
    design = design,
    group = group_vector,
    min.count = 10,
    min.total.count = 500,
    large.n = 100,
    min.prop = 0.9)

  y = y[keep, , keep.lib.sizes=FALSE]
  y = calcNormFactors(y)
  y = estimateDisp(y,design,robust=TRUE)
  fit = glmQLFit(y, design, robust=TRUE)
  contrast_to_test = makeContrasts("Virus_specific - Not_Virus_specific",levels = design)
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
    scale_color_manual(values = colours)+
    theme(panel.grid = element_blank())

  png(paste0(plot_title,"_virus_specifc_vs_not.png"),res=300,units="in",height=4,width=4)
  print(plot)
  dev.off()

  # save
  outfile_csv = paste0(plot_title,"_virus_specifc_vs_not.csv")
  write_csv(res,outfile_csv)

}

# EBV
t_cells@meta.data = t_cells@meta.data %>%
  mutate(virus_specific = ebv_specific) %>%
  mutate(donor_virus = paste0(iid,"_",virus_specific))
rownames(t_cells@meta.data) = colnames(t_cells)
do_de_viral(subset(t_cells,ann_celltypist_highres=="Tem/Trm cytotoxic T cells" & phenotype=="MS" & source=="CSF"),"EBV_specific_ms_csf_Tem_rm",mincells=2)
do_de_viral(subset(t_cells,ann_celltypist_highres=="Tem/Trm cytotoxic T cells" & phenotype=="OIND" & source=="CSF"),"EBV_specific_oind_csf_Tem_rm",mincells=2)

# CMV
t_cells@meta.data = t_cells@meta.data %>%
  mutate(virus_specific = cmv_specific) %>%
  mutate(donor_virus = paste0(iid,"_",virus_specific))
rownames(t_cells@meta.data) = colnames(t_cells)
do_de_viral(subset(t_cells,ann_celltypist_highres=="Tem/Trm cytotoxic T cells" & phenotype=="MS" & source=="CSF"),"CMV_specific_ms_csf_Tem_rm",mincells=2)
do_de_viral(subset(t_cells,ann_celltypist_highres=="Tem/Trm cytotoxic T cells" & phenotype=="OIND" & source=="CSF"),"CMV_specific_oind_csf_Tem_rm",mincells=2)

# compare
ebv_ms = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/EBV_specific_ms_csf_Tem_rm_virus_specifc_vs_not.csv")
cmv_ms = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/CMV_specific_ms_csf_Tem_rm_virus_specifc_vs_not.csv")
ebv_oind = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/EBV_specific_oind_csf_Tem_rm_virus_specifc_vs_not.csv")

# join
combo_dat = ebv_ms %>%
  left_join(cmv_ms,by="gene")


p1 = ggplot(combo_dat %>% filter(P_adj.x < 0.05),aes(-log10(PValue.x),-log10(PValue.y),label=gene))+
  geom_point()+
  labs(x="-log10 P(EBV)",y="-log10 P(CMV)")+
  ggrepel::geom_text_repel(data = combo_dat %>% filter(P_adj.x < 0.05))+
  theme_bw()

p2 = ggplot(combo_dat %>% filter(P_adj.x < 0.05),aes(logFC.x,logFC.y,label=gene))+
  geom_point()+
  labs(x="logFC(EBV)",y="logFC(CMV)")+
  ggrepel::geom_text_repel(data = combo_dat %>% filter(P_adj.x < 0.05))+
  theme_bw()


# find matches in VDJDB
trb_matches = t_cells@meta.data %>%
  dplyr::rename("cdr3" = junction_aa_VDJ) %>%
  left_join(vdjdb,by="cdr3") %>%
  mutate(antigen.species = ifelse(is.na(antigen.species),"No match",antigen.species))

trb_matches %>%
  dplyr::count(antigen.species != "No match")

expanded_plot_dat = trb_matches %>%
group_by(expanded_clone,phenotype) %>%
  dplyr::count(antigen.species) %>%
    mutate(prop = n/sum(n)) %>%
    print(n=100)

p = ggplot(expanded_plot_dat %>% filter(antigen.species != "No match"),
aes(antigen.species,prop,fill=expanded_clone))+
geom_col(position= position_dodge(),color="black")+
facet_wrap(~phenotype,nrow=1)+
scale_fill_brewer(palette="Set1")+
theme_minimal()+
labs(x="Antigen target", y="Proportion",fill="Expanded clone?")+
coord_flip()+
scale_y_continuous(breaks = c(0,0.02))

png("clonal_epitopes.png",res=600,units="in",width=8,height=3)
p
dev.off()

# repeat just CSF
expanded_plot_dat_csf = trb_matches %>%
filter(source=="CSF") %>%
group_by(expanded_clone,phenotype) %>%
  dplyr::count(antigen.species) %>%
    mutate(prop = n/sum(n))

p = ggplot(expanded_plot_dat_csf %>% filter(antigen.species != "No match"),
aes(antigen.species,prop,fill=expanded_clone))+
geom_col(position= position_dodge(),color="black")+
facet_wrap(~phenotype,nrow=1)+
scale_fill_brewer(palette="Set1")+
theme_minimal()+
labs(x="Antigen target", y="Proportion",fill="Expanded clone?")+
coord_flip()+
scale_y_continuous(breaks = c(0,0.02))

png("clonal_epitopes_just_csf.png",res=600,units="in",width=8,height=3)
p
dev.off()

# repeat with epitopes
expanded_plot_dat_csf = trb_matches %>%
filter(source=="CSF") %>%
group_by(expanded_clone,phenotype) %>%
  dplyr::count(antigen.species,antigen.gene) %>%
    mutate(prop = n/sum(n))


p = ggplot(expanded_plot_dat_csf %>% filter(antigen.species != "No match") %>% mutate(antigen = paste0(antigen.species,"_",antigen.gene)),
aes(antigen,prop,fill=expanded_clone))+
geom_col(position= position_dodge(),color="black")+
facet_wrap(~phenotype,nrow=1)+
scale_fill_brewer(palette="Set1")+
theme_minimal()+
labs(x="Antigen target", y="Proportion",fill="Expanded clone?")+
coord_flip()+
scale_y_continuous(breaks = c(0,0.02))

png("clonal_epitopes_just_csf_full_epitope.png",res=600,units="in",width=8,height=6)
p
dev.off()




# find exact matches in VDJDB
trb_exact_matches = t_cells@meta.data %>%
  dplyr::rename("cdr3" = junction_aa_VDJ) %>%
  left_join(vdjdb %>% filter(Chain=="TRB"),by="cdr3") %>%
  filter(v_call_simple == v.segm & j_call_VDJ == j.segm) %>%
  filter(!is.na(antigen.species) &  antigen.species != "No match")

# find exact matches in VDJDB
tra_exact_matches = t_cells@meta.data %>%
  dplyr::rename("cdr3" = junction_aa_VJ) %>%
  left_join(vdjdb %>% filter(Chain=="TRA"),by="cdr3") %>%
  filter(trav_call_simple == v.segm & j_call_VJ == j.segm) %>%
  filter(!is.na(antigen.species) &  antigen.species != "No match")

matches = t_cells@meta.data  %>%
filter(
  original_barcode %in% tra_exact_matches$original_barcode &
  original_barcode %in% trb_exact_matches$original_barcode
  )

# check same target
full_matches = list()
for(i in c(1:nrow(matches))){

  this_cell = matches[i,]
  tra = tra_exact_matches  %>% filter(original_barcode == this_cell$original_barcode)
  trb = trb_exact_matches  %>% filter(original_barcode == this_cell$original_barcode)

  # loop through trb matches
  for(j in c(1:nrow(trb))){
    this_trb = trb[j,]

    matched_tras = tra %>%
      filter(antigen.species == this_trb$antigen.species &
      antigen.epitope == this_trb$antigen.epitope &
      mhc.class == this_trb$mhc.class &
      mhc.a == this_trb$mhc.a
      )

    full_matches[[length(full_matches)+1]] = this_trb %>%
    dplyr::select(original_barcode,v.segm,antigen.species,antigen.epitope,mhc.class) %>%
    left_join(matched_tras  %>%
    dplyr::select(original_barcode,v.segm,antigen.species,antigen.epitope,mhc.class),
      by = c("original_barcode","antigen.species","antigen.epitope","mhc.class"))

}
}

full_matches = do.call("bind_rows",full_matches)

full_matches = full_matches %>% filter(!is.na(v.segm.y)) %>%
left_join(t_cells@meta.data,by="original_barcode")

# LOOK AT SPECIFIC MATCHES
ebv_matches = full_matches %>% filter(antigen.species == "EBV")
# make indicator in main dataset
ebv = t_cells@meta.data %>%
  mutate(has_specific_match = ifelse(
  original_barcode %in% ebv_matches$original_barcode,
  "ebv_match",
  "not_ebv_match"
  ))

plot_dat = ebv %>%
group_by(phenotype,source) %>%
dplyr::count(has_specific_match) %>%
mutate(prop = n/sum(n)) %>%
filter(has_specific_match=="ebv_match")

ggplot(plot_dat ,
aes(phenotype,prop,label=n,fill=source))+
geom_col(color="black",position=position_dodge())+
scale_fill_brewer(palette="Set1")+
theme_minimal()

# REPEAT FOR CMV
ebv_matches = full_matches %>% filter(antigen.species == "CMV")
# make indicator in main dataset
ebv = t_cells@meta.data %>%
  mutate(has_specific_match = ifelse(
  original_barcode %in% ebv_matches$original_barcode,
  "ebv_match",
  "not_ebv_match"
  ))

plot_dat = ebv %>%
group_by(phenotype,source) %>%
dplyr::count(has_specific_match) %>%
mutate(prop = n/sum(n)) %>%
filter(has_specific_match=="ebv_match")

ggplot(plot_dat ,
aes(phenotype,prop,label=n,fill=source))+
geom_col(color="black",position=position_dodge())+
scale_fill_brewer(palette="Set1")+
theme_minimal()


# filter out low counts
lowcounts = trb_matches %>%
  dplyr::count(antigen.species) %>%
  arrange(n) %>%
  filter(n<50)
trb_matches = trb_matches %>%
  filter(!antigen.species %in% lowcounts$antigen.species)

# plot
plot_dat = trb_matches %>%
dplyr::count(source,phenotype,antigen.species) %>%
group_by(source,phenotype) %>%
mutate(prop = n/sum(n))

# permute
n_perm = 1000
overall_res = list()
trb_matches = trb_matches %>%
dplyr::select(source,phenotype,antigen.species)
for(i in c(1:n_perm)){
  message(i)

  sampled_dat = sample_n(trb_matches,size = nrow(trb_matches),replace=T)
  plot_dat = sampled_dat %>%
  dplyr::count(source,phenotype,antigen.species) %>%
  group_by(source,phenotype) %>%
  mutate(prop = n/sum(n)) %>%
  pivot_wider(
  id_cols = c(source,antigen.species),
  names_from = phenotype,
  values_from = prop) %>%
  mutate(ms_max = ifelse(
    (MS > OINDI | is.na(OINDI) ) &
    (MS > OIND | is.na(OIND) ) &
    (MS > NIND | is.na(NIND) ),
    "MS",
    "Other"
  ))
  overall_res[[i]] = plot_dat
}
overall_res = do.call("bind_rows",overall_res)

p = ggplot(plot_dat %>%
mutate(prop = as.numeric(prop)) %>%
filter(antigen.species != "No match"),aes(antigen.species,prop,fill=phenotype))+
geom_col(color="black",position = position_dodge())+facet_wrap(~source)+
coord_flip()+
scale_fill_brewer(palette="Set1")+
theme_minimal()+
labs(y="Proportion",x="Antigen target")+
scale_y_continuous(breaks=c(0,0.005,0.01))
png("epitope_specificity.png",res=600,units="in",width=5,height=3)
p
dev.off()

# pval vs nind
pval_vs_nind = overall_res %>%
mutate(MS_bigger_than_NIND = ifelse(MS>NIND,"yes","no")) %>%
group_by(source,antigen.species) %>%
dplyr::count(MS_bigger_than_NIND) %>%
mutate(prop = (n)/(sum(n)+1) ) %>%
filter(MS_bigger_than_NIND == "yes") %>%
mutate(pval = 1-prop) %>%
arrange(pval)

# pval vs oind
pval_vs_oind = overall_res %>%
mutate(MS_bigger_than_OIND = ifelse(MS>OIND,"yes","no")) %>%
group_by(source,antigen.species) %>%
dplyr::count(MS_bigger_than_OIND) %>%
mutate(prop = (n)/(sum(n)+1) ) %>%
filter(MS_bigger_than_OIND == "yes") %>%
mutate(pval = 1-prop) %>%
arrange(pval)

pval_vs_oindi = overall_res %>%
mutate(MS_bigger_than_OINDI = ifelse(MS>OINDI,"yes","no")) %>%
group_by(source,antigen.species) %>%
dplyr::count(MS_bigger_than_OINDI) %>%
mutate(prop = (n)/(sum(n)+1) ) %>%
filter(MS_bigger_than_OINDI == "yes") %>%
mutate(pval = 1-prop) %>%
arrange(pval)

overall_epitope_spec = bind_rows(pval_vs_nind,pval_vs_oind,pval_vs_oindi) %>%
mutate(p_adj = p.adjust(pval,method="fdr"))



# DE EBV-specific vs other
