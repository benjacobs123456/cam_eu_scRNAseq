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

t_cell_markers = c("CD3G","IL7R","CD8A","TCF7","FOXP3","CCL5","CCR7","GZMH","CD27")

t_subsets_celltypist = c("Tcm/Naive helper T cells",
"Tem/Trm cytotoxic T cells",
"Regulatory T cells",
"Tem/Temra cytotoxic T cells",
"Tem/Effector helper T cells",
"Memory CD4+ cytotoxic T cells",
"Type 17 helper T cells",
"MAIT cells",
"Tcm/Naive cytotoxic T cells",
"Type 1 helper T cells",
"Cycling T cells",
"CD8a/a",
"Follicular helper T cells",
"CD8a/b(entry)",
"Early lymphoid/T lymphoid",
"NKT cells",
"Treg(diff)",
"Trm cytotoxic T cells",
"Tem/Effector helper T cells PD1+")

# filter by annotation
t_cells = subset(t_cells, subset = ann_celltypist_highres %in% t_subsets_celltypist)


# recluster
set.seed(1)
DefaultAssay(t_cells)="SCT"
t_cells = RunUMAP(t_cells,reduction="harmony",dims=1:50)
t_cells = FindNeighbors(t_cells)
t_cells = FindClusters(t_cells,resolution=0.5)
t_cells = SetIdent(t_cells,value = t_cells[['ann_celltypist_highres']])

# plots to check annotations
p1=DimPlot(t_cells)
p2=FeaturePlot(t_cells,features=t_cell_markers)
p3=DotPlot(t_cells,features=t_cell_markers)

png("dimplot.png",res=300,units="in",width=6,height=4)
p1
dev.off()

png("featureplot.png",res=300,units="in",width=8,height=8)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=10,height=4)
p3
dev.off()


##############################
# cluster biomarkers
##############################
#biomarkers = FindAllMarkers(t_cells,recorrect_umi=FALSE,logfc.threshold=0.25,min.pct=0.25,only.pos=TRUE)
#topbiomarkers = biomarkers %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=5)
#write_csv(topbiomarkers,"biomarkers.csv")
#write_csv(biomarkers,"all_biomarkers.csv")

#######################################
# DA & composition
#######################################

common_cell_types = t_cells@meta.data %>% dplyr::count(cell_type) %>% filter(n>20)

# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))

t_cells@meta.data$cell_type = t_cells@meta.data$ann_celltypist_highres
n_col = t_cells@meta.data$cell_type %>% unique %>% length
colour_pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_col)

t_cells@meta.data$cell_type = factor(
t_cells@meta.data$cell_type,
ordered=TRUE,
levels = c("Tcm/Naive helper T cells",
"Tem/Trm cytotoxic T cells",
"Regulatory T cells",
"Tem/Temra cytotoxic T cells",
"Tem/Effector helper T cells",
"Memory CD4+ cytotoxic T cells",
"Type 17 helper T cells",
"MAIT cells",
"Tcm/Naive cytotoxic T cells",
"Type 1 helper T cells",
"Cycling T cells",
"CD8a/a",
"Follicular helper T cells",
"CD8a/b(entry)",
"Early lymphoid/T lymphoid",
"NKT cells",
"Treg(diff)",
"Trm cytotoxic T cells",
"Tem/Effector helper T cells PD1+"))


p=DimPlot(t_cells,label=F,raster=F,group.by="cell_type")+
scale_color_manual(values = colour_pal)+
ggtitle("")+
theme_minimal()+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("dim_plot_simple_labels.png",res=300,units="in",width=5,height=4)
p
dev.off()

# split by source
p=DimPlot(subset(t_cells,ann_celltypist_highres %in% common_cell_types$cell_type),label=F,raster=F,group.by="cell_type",split.by="source")+
scale_color_manual(values = colour_pal)+
ggtitle("")+
theme_minimal()+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("dim_plot_simple_labels_source.png",res=300,units="in",width=5,height=5)
p
dev.off()

# stash sample info
sample_info = t_cells@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(t_cells@meta.data$cell_type,t_cells@meta.data$donor_source)

# filter out clusters with 0 counts
abundances = abundances[rowSums(abundances)>20,]

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
  MS_CSF - OIND_CSF,
  MS_CSF - Control_CSF,
  OIND_CSF - Control_CSF,
  MS_PBMC - OIND_PBMC,
  MS_PBMC - Control_PBMC,
  OIND_PBMC - Control_PBMC,
  MS_CSF - MS_PBMC,
  OIND_CSF - OIND_PBMC,
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
  "Type 17 helper T cells",
  "MAIT cells",
  "Tcm/Naive cytotoxic T cells",
  "Type 1 helper T cells",
  "Cycling T cells",
  "CD8a/a",
  "Follicular helper T cells",
  "CD8a/b(entry)",
  "Early lymphoid/T lymphoid",
  "NKT cells",
  "Treg(diff)",
  "Trm cytotoxic T cells",
  "Tem/Effector helper T cells PD1+"))


  # da plot

  plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
  theme_classic()+
  scale_color_manual(values = colour_pal)+
  geom_text_repel(data = res %>% filter(PValue < 0.005 / nrow(res)),
  mapping = aes(logFC,-log10(PValue)),
  size=2,force=100,max.time=5,max.overlaps=100)+
  geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
  geom_vline(xintercept = 0,alpha=0.2)+
  NoLegend()+
  ggtitle(comparison_label)+
  scale_y_continuous(limits=c(0,30))+
  scale_x_continuous(limits=c(-9,9))+
  geom_point(shape=16,size=3)


  png(paste0("./crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
  print(plot)
  dev.off()
}

#######################################
# Phenotypes & DA
#######################################

# reorder phenotypes
t_cells@meta.data$phenotype = factor(t_cells@meta.data$phenotype,levels=c("Noninflammatory","OIND","MS"),ordered=TRUE)

# proportion plots
png("./crude_labels_proportion_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(t_cells@meta.data %>% filter(cell_type %in% common_cell_types$cell_type),aes(phenotype,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Phenotype",y="Proportion",fill="Cell type")
dev.off()

plots = list()
for(pheno in c("MS","OIND","Noninflammatory")){
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

# relationship between cell type proportions and phenotypes
# define ocb positivity
t_cells@meta.data = t_cells@meta.data %>%
mutate(oligo_pos = ifelse(Oligoclonal %in% c("Pos(>10)","Pos(3to9)"),"OCB+","OCB-")) %>%
mutate(oligo_pos = ifelse(cohort=="EU" & is.na(Oligoclonal),"OCB+",oligo_pos)) %>%
mutate(drb_pos = ifelse(drb1_1501_dose > 0 ,"DRB1*15+","DRB1*15-"))

plots = list()
for(pheno in c("RMS","PPMS")){
p=ggplot(t_cells@meta.data %>%
filter(Category_fine==pheno),
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

for(oligo_status in c("OCB+","OCB-")){
p=ggplot(t_cells@meta.data %>%
filter(phenotype=="MS" & oligo_pos == oligo_status),
aes(donor.id,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Donor",y="Proportion",fill="Cell type")+
ggtitle(oligo_status)+
theme(legend.position="none",axis.text.x=element_blank())
plots[[length(plots)+1]] = p
}

for(drb_status in c("DRB1*15+","DRB1*15-")){
p=ggplot(t_cells@meta.data %>%
filter(phenotype=="MS" & drb_pos == drb_status & !is.na(drb_pos)),
aes(donor.id,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Donor",y="Proportion",fill="Cell type")+
ggtitle(drb_status)+
theme(legend.position="none",axis.text.x=element_blank())
plots[[length(plots)+1]] = p
}

png("./crude_labels_proportion_cell_counts_barplot_per_individual_by_ms_subtype.png",res=300,width=12,height=4,units="in")
do.call("grid.arrange",plots)
dev.off()

# repeat DA pairwise between OCB +/-, DRB, and PPMS vs RMS

# first clean variables
t_cells@meta.data = t_cells@meta.data %>%
mutate(oligo_pos = ifelse(
oligo_pos == "OCB+",
"OCBpos",
"OCBneg"
))

t_cells@meta.data = t_cells@meta.data %>%
mutate(drb_pos = ifelse(
drb_pos == "DRB1*15+",
"DRBpos",
ifelse(!is.na(drb_pos),
"DRBneg",
NA
)))

do_da = function(x, contrasts_to_test){
  # stash sample info
  sample_info = t_cells@meta.data %>%
  filter(phenotype=="MS") %>%
  dplyr::select(donor.id,source,x,donor_source) %>%
  filter(!is.na(.data[[x]])) %>%
  distinct(donor_source,.keep_all=TRUE)

  data_for_abundances = t_cells@meta.data %>%
  filter(phenotype=="MS") %>%
  dplyr::select(donor.id,source,x,donor_source,cell_type) %>%
  filter(!is.na(.data[[x]]))
  abundances = table(data_for_abundances$cell_type,data_for_abundances$donor_source)

  # filter out clusters with <10 counts
  abundances = abundances[rowSums(abundances)>10,]

  sample_info = sample_info %>% filter(donor_source %in% colnames(abundances))
  sample_info = sample_info[match(colnames(abundances),sample_info$donor_source),]
  sample_info$grouping = paste0(sample_info[[x]],"_",sample_info$source)


  y.ab = DGEList(abundances, samples=sample_info)

  # update sample info
  y.ab$samples =  y.ab$samples %>% left_join(t_cells@meta.data %>%
  filter(!is.na(.data[[x]])) %>%
  distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Gender,donor_source),by="donor_source")

  # now loop
  da_overall_list = list()
  results_df = data.frame()
  design <- model.matrix(~ 0 + grouping + Age + Gender,y.ab$sample)
  colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Gender")

  y.ab <- estimateDisp(y.ab, design,trend="none")
  fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

  # define contrast for testing

  contrast =  makeContrasts(
    contrasts = contrasts_to_test,
    levels = design
  )

  da_plots = list()
  # do da tests
  for(i in c(1:2)){
    contrast_name = colnames(contrast)[i]
    res = glmQLFTest(fit, contrast = contrast[,i])

    # write to file
    print(summary(decideTests(res)))
    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$cell = rownames(res)
    res$significant = ifelse(res$P_adj<0.01,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
    comparison_label = contrast_name
    # refactor cell types (for plotting)
    res$cell = factor(res$cell,ordered=TRUE,
    levels = c(
    "Tcm/Naive helper T cells",
    "Tem/Trm cytotosix T cells",
    "Tcm/Naive cytotoxic T cells",
    "Regulatory T cells",
    "MAIT T cells",
    "Tem/Effector helper T cells",
    "Tem/Temra cytotoxic T cells",
    "Cycling T cells",
    "Early lymphoid/T lymphoid",
    "Memory CD4+ cytotoxic T cells",
    "CD8a/a",
    "Follicular helper T cells",
    "Type 1 helper T cells",
    "Type 17 helper T cells",
    "Trm cytotoxic T cells"))

    # da plot
    plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
    theme_classic()+
    scale_color_manual(values = colour_pal)+
    geom_text_repel(size=3,max.overlaps=100)+
    geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
    geom_vline(xintercept = 0,alpha=0.2)+
    NoLegend()+
    ggtitle(comparison_label)+
    scale_y_continuous(limits=c(0,10))+
    scale_x_continuous(limits=c(-9,9))+
    geom_point(shape=16,size=3)


    png(paste0("./pheno_comparisons_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
    print(plot)
    dev.off()
  }
}

# run DA
do_da(x = "Category_fine",
contrasts_to_test = c("RMS_CSF - PPMS_CSF","RMS_PBMC - PPMS_PBMC"))
do_da(x = "oligo_pos",
contrasts_to_test = c("OCBpos_CSF - OCBneg_CSF","OCBpos_PBMC - OCBneg_PBMC"))
do_da(x = "drb_pos",
contrasts_to_test = c("DRBpos_CSF - DRBneg_CSF","DRBpos_PBMC - DRBneg_PBMC"))


p=DimPlot(subset(t_cells,phenotype=="MS" & source=="CSF"),split.by="donor.id",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),plot.title = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())
png("./indiv_dimplots_ms.png",res=300,units="in",height=8,width=8)
p
dev.off()

#######################################
# DE
#######################################

DefaultAssay(t_cells)="RNA"

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
groups = colData(t_cells.sce)[, c("ident", "phenotype","source","donor.id")]
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

    cell_type_plot = str_replace_all(cell_type,"/","_")

    bonf = 0.05 / length(res$PValue)
    res$significant = ifelse(res$PValue<bonf,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
    res = res %>% mutate(direction = ifelse(PValue>bonf,"nonsig",direction))
    outfile = paste0("edgeR_de_tests_",contrast_name,"_",cell_type_plot,".csv")
    write_csv(
      res,
      file = outfile)

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

    png(paste0("de_plot_",contrast_name,"_",cell_type_plot,".png"),res=300,units="in",height=4,width=4)
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

clusters = t_cells$cell_type %>% unique

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

      cluster = str_replace_all(levels(factor(clusters))[x],"/","_")
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

  png(paste0("summary_",x,".png"),res=300,units="in",width=6,height=6)
  p=ggplot(plot_data,aes(cell_type,pathway,fill=NES,label=ifelse(padj<0.01,ifelse(padj<0.001,ifelse(padj<0.0001,"***","**"),"*"),"")))+
    geom_tile(color="black")+
    geom_text()+
    scale_fill_gradient(low="purple",high="orange")+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)+
    theme(axis.text.x = element_text(angle=90))
  print(p)
  dev.off()
}

# run for other sets (go, kegg, reactome)
#res_overall = list()
#do_gsea("go")
#do_gsea("kegg")
#do_gsea("reactome")

#######################################
# TCR analysis
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


#######################################
# Repertoire analysis CSF v periphery
#######################################
rownames(t_cells@meta.data) = colnames(t_cells)

## define variables
# rename NIND (for graphs)
t_cells@meta.data$phenotype = recode(t_cells@meta.data$phenotype, "Noninflammatory" = "NIND")

# define clones
expanded_clones = t_cells@meta.data %>% group_by(donor.id,clone_id) %>% dplyr::count() %>% filter(n>1) %>% mutate(donor_clone = paste0(donor.id,"_",clone_id))
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_clone = paste0(donor.id,"_",clone_id)) %>%
mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

# define TRBV types
t_cells@meta.data$trbv_family = sapply(t_cells@meta.data$v_call_VDJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  return(y)
})


## plots

# cell types
# compare csf vs periphery
png("csf_v_pbmc_celltypes.png",res=300,height=3,width=6,units="in")
DimPlot(t_cells,split.by="source")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

## trbv
png("csf_v_pbmc_trbv.png",res=300,height=3,width=6,units="in")
DimPlot(t_cells,split.by="source",group.by="trbv_family")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


# compare csf vs periphery - just ms
csf_v_pbmc_plot = function(x){
DimPlot(subset(t_cells,phenotype=="MS"),split.by="source",group.by=x)+
scale_color_manual(values = colour_pal)+
theme_minimal()+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
var_list = list("cell_type","expanded_clone","trbv_family","status_summary")
plot_list = lapply(var_list,csf_v_pbmc_plot)

png("csf_v_pbmc_just_ms_tcr_repertoire.png",res=300,units="in",height=8,width=8)
do.call("grid.arrange",plot_list)
dev.off()


# bar plots
make_categorical_plot = function(x){
  ggplot(t_cells@meta.data,aes(phenotype,fill=t_cells@meta.data[[x]]))+
  geom_bar(position="fill",color="black")+
  facet_wrap(~source)+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  labs(x="Phenotype",fill=x,y="Proportion")
}

make_continuous_plot = function(x){
  ggplot(t_cells@meta.data,aes(phenotype,fill=cell_type,y=t_cells@meta.data[[x]]))+
  geom_boxplot(color="black")+
  facet_wrap(~source)+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  labs(x="Phenotype",y=x)
}

var_list = list("expanded_clone","trbv_family","status_summary","c_call_VDJ")
plot_list = lapply(var_list,make_categorical_plot)

png("csf_v_pbmc_tcr_repertoire.png",res=300,units="in",height=10,width=10)
do.call("grid.arrange",plot_list)
dev.off()

#######################################
#  TRBV usage
#######################################

do_da_var = function(x, contrasts_to_test,variable,ylim=10, data = t_cells@meta.data,plot_title){

  # get rid of donors with low counts
  donors_to_keep = data %>% dplyr::count(donor_source) %>% filter(n>1)
  da_dat = data %>% filter(donor_source %in% donors_to_keep$donor_source)

  # stash sample info
  sample_info = da_dat %>%
  dplyr::select(donor.id,source,x,donor_source) %>%
  filter(!is.na(.data[[x]])) %>%
  distinct(donor_source,.keep_all=TRUE)

  data_for_abundances = da_dat %>%
  dplyr::select(donor.id,source,x,donor_source,.data[[variable]]) %>%
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
  distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Gender,donor_source),by="donor_source")

  da_overall_list = list()
  results_df = data.frame()
  design <- model.matrix(~ 0 + grouping + Age + Gender,y.ab$sample)
  colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Gender")

  y.ab <- estimateDisp(y.ab, design,trend="none")
  fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

  # define contrast for testing

  contrast =  makeContrasts(
    contrasts = contrasts_to_test,
    levels = design
  )

  da_plots = list()
  # do da tests
  for(i in c(1:length(contrasts_to_test))){
    contrast_name = colnames(contrast)[i]
    res = glmQLFTest(fit, contrast = contrast[,i])

    # write to file
    print(summary(decideTests(res)))
    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$cell = rownames(res)
    res$significant = ifelse(res$P_adj<0.01,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
    comparison_label = contrast_name

    # da plot
    plot = ggplot(res,aes(logFC,-log10(PValue),label=cell))+
    theme_classic()+
    geom_text_repel(data = res %>% filter(P_adj<0.05),
      mapping = aes(logFC,-log10(PValue)),
      size=3,max.overlaps=100)+
    geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
    geom_vline(xintercept = 0,alpha=0.2)+
    NoLegend()+
    ggtitle(comparison_label)+
    scale_y_continuous(limits=c(0,ylim))+
    scale_x_continuous(limits=c(-5,5))+
    geom_point(shape=16,size=2)


    png(paste0(plot_title,"_pheno_comparisons_da_plot_",contrast_name,"_",variable,".png"),res=300,units="in",height=3,width=3)
    print(plot)
    dev.off()
    write_csv(res,
    file=paste0(plot_title,"_pheno_comparisons_da_plot_",contrast_name,"_",variable,".csv"))
  }
}

# run DA
do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "NIND_CSF - NIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_CSF - NIND_CSF"),
variable="trbv_family",
plot_title = "trb_family",
ylim=20)

# individual gene segments
do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF"),
variable="v_call_VDJ",
plot_title = "trb_allgenes",
ylim=20)

# alpha chain
do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF"),
variable="v_call_VJ",
plot_title = "tra_allgenes",
ylim=20)


#######################################
# Clonal analysis
#######################################

main_clonal_cells = t_cells@meta.data %>% dplyr::group_by(expanded_clone) %>% dplyr::count(cell_type) %>% mutate(prop = n / sum(n)) %>% filter(expanded_clone=="Expanded") %>% arrange(desc(prop))%>% filter(prop>0.01)

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

binary_comparison = function(variable,plot_title,level="Yes"){
  props = t_cells@meta.data %>%
    group_by(donor.id,source,phenotype) %>%
    dplyr::count(.data[[variable]]) %>%
    mutate(prop = n/sum(n)) %>%
    filter(.data[[variable]]==level)

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
binary_comparison_celltypes = function(variable,plot_title,level="Yes"){
  props = t_cells@meta.data %>%
  filter(cell_type %in% main_clonal_cells$cell_type) %>%
    group_by(donor.id,source,phenotype,cell_type) %>%
    dplyr::count(.data[[variable]]) %>%
    mutate(prop = n/sum(n)) %>%
    filter(.data[[variable]]==level) %>%
    mutate(cell_pheno  = paste0(phenotype,"_",cell_type))

  phenos = unique(props$cell_pheno)

  pvals = list()
  for(i in c(1:length(phenos))){
  message(i)
    a = props[props$source=="CSF" & props$cell_pheno==phenos[i],][['prop']]
    b = props[props$source=="PBMC" & props$cell_pheno==phenos[i],][['prop']]

    p=if(length(a)<=1 | length(b) <=1){
      NA
    } else if (max(a) == min(a) & max(a) == max(b) & max(b) == min(b)){
      1
    } else {
      t.test(a,b)$p.value
    }
    message(p)
    pvals[[length(pvals)+1]] = p
  }
  pvals = data.frame(cell_pheno = phenos,P = unlist(pvals))

  props = props %>%
    left_join(pvals,by="cell_pheno") %>%
    mutate(P = simplify_pval(P))

  p=ggplot(props,aes(phenotype,prop,fill=source))+
    geom_boxplot(color="black")+
    facet_wrap(~cell_type)+
    scale_fill_brewer(palette="Set1")+
    ggtitle(plot_title)+
    geom_text(mapping = aes(x = phenotype,y = 1.05,label = P))+
    theme_minimal()+
    labs(x="Phenotype",y="Proportion")


  png(file=paste0(plot_title,"_comparison.png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
}

binary_comparison(variable = "expanded_clone",plot_title="Clonal expansion",level="Expanded")
binary_comparison_celltypes(variable = "expanded_clone",plot_title="Clonal expansion (by cell type)",level="Expanded")

# compare csf in ms vs controls

common_cell_types = t_cells@meta.data %>% dplyr::count(cell_type) %>% filter(n>200)
ms_v_cont_csf_plot = function(x){
DimPlot(subset(t_cells,source=="CSF" & cell_type %in% common_cell_types$cell_type),split.by="phenotype",group.by=x)+
scale_color_brewer(palette="Set1")+
theme_minimal()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

p0=ms_v_cont_csf_plot("cell_type")+ggtitle("Cell type")
p1=ms_v_cont_csf_plot("expanded_clone")+ggtitle("Clonal expansion")
png("ms_v_cont_csf_tcr_repertoire.png",res=300,units="in",height=3,width=3)
grid.arrange(p0,p1)
dev.off()

# big picture numbers
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
donors_with_expanded_clones = t_cells@meta.data  %>% group_by(donor.id) %>% dplyr::count(donor_clone) %>% filter(n>1) %>% distinct(donor.id)
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_with_expanded_clones = ifelse(donor.id %in% donors_with_expanded_clones$donor.id,"Yes","No"))

donors_with_expanded_clones_csf = t_cells@meta.data  %>% filter(source=="CSF") %>% group_by(donor.id) %>% dplyr::count(donor_clone) %>% filter(n>1) %>% distinct(donor.id)
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_with_expanded_clones_csf = ifelse(donor.id %in% donors_with_expanded_clones_csf$donor.id,"Yes","No"))

# get some counts
t_cells@meta.data %>%
  group_by(phenotype) %>%
  distinct(donor.id,.keep_all=T) %>%
  dplyr::count(donor_with_expanded_clones) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(donor_with_expanded_clones=="Yes")

t_cells@meta.data %>%
  group_by(phenotype) %>%
  distinct(donor.id,.keep_all=T) %>%
  dplyr::count(donor_with_expanded_clones_csf) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(donor_with_expanded_clones_csf=="Yes")

# histogram
p=ggplot(t_cells@meta.data,aes(clonal_size,fill=phenotype))+geom_histogram()+theme_bw()+scale_fill_brewer(palette="Set2")+labs(x="Clone size")
png("clonal_histogram.png",res=300,height=2,width=4,units="in")
p
dev.off()

# clonality of csf and pbmc
t_cells@meta.data$expanded_clone = factor(t_cells@meta.data$expanded_clone,levels=c("Not expanded","Expanded"),ordered=TRUE)

p1=ggplot(t_cells@meta.data %>% filter(source=="CSF"),aes(phenotype,fill=expanded_clone))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Expanded clone?",x="Phenotype",y="Proportion of B cell pool")+
png("clonal_proportions_source.png",res=300,height=4,width=4,units="in")
p1
dev.off()



#######################################
# Clonal relationships
#######################################

# find public clones in the dataset

public_clones = t_cells@meta.data %>% group_by(clone_id,donor.id) %>% dplyr::count() %>% group_by(clone_id) %>% dplyr::count(donor.id) %>% arrange(desc(n)) %>% filter(n>1)
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

t_cells@meta.data$expanded_clone = factor(t_cells@meta.data$expanded_clone,ordered=T)
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

p1=ggplot(tcr_db_gex,aes(n,antigen.species,fill=expanded_clone))+geom_col(color="black",position=position_dodge())
p2=ggplot(tcr_db_gex,aes(n/total_cells,antigen.species,fill=expanded_clone))+geom_col(color="black",position=position_dodge())

png("pathogen_tcrs_by_expanded_just_ms.png",res=300,units="in",height=6,width=12)
grid.arrange(p1,p2,nrow=1)
dev.off()

# repeat for clonal vs nonclonal in controls
# whole dataset
tcr_db_gex = t_cells@meta.data %>% filter(phenotype=="NIND") %>%
left_join(
vdj_db %>% dplyr::rename("junction_aa_VDJ" = cdr3),
by = c("junction_aa_VDJ")) %>%
filter(gene=="TRB") %>%
filter(!is.na(antigen.epitope)) %>%
dplyr::count(antigen.species,expanded_clone)

totals = t_cells@meta.data %>% group_by(expanded_clone) %>% dplyr::count() %>% dplyr::rename("total_cells" = n)
tcr_db_gex = tcr_db_gex %>% left_join(totals,by="expanded_clone")

p1=ggplot(tcr_db_gex,aes(n,antigen.species,fill=expanded_clone))+geom_col(color="black",position=position_dodge())
p2=ggplot(tcr_db_gex,aes(n/total_cells,antigen.species,fill=expanded_clone))+geom_col(color="black",position=position_dodge())

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


bootstrap_ebv_enrichment = function(pheno1 = "NIND",
clonal1 =  "Expanded",
pheno2 = "NIND",
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

bootstrap_ebv_enrichment(pheno1 = "NIND",clonal1 =  "Expanded",pheno2 = "NIND",clonal2 = "Not expanded")
bootstrap_ebv_enrichment(pheno1 = "MS",clonal1 =  "Expanded",pheno2 = "MS",clonal2 = "Not expanded")
bootstrap_ebv_enrichment(pheno1 = "OIND",clonal1 =  "Expanded",pheno2 = "OIND",clonal2 = "Not expanded")


# combine with bcr data
