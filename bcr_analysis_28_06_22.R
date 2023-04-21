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

b_cell_markers = c("IGHD","IGHG1","IGHM","CD27","CD38","MKI67")

# recluster
set.seed(1)
DefaultAssay(b_cells)="SCT"
b_cells = RunUMAP(b_cells,reduction="harmony",dims=1:50)
b_cells = FindNeighbors(b_cells)
b_cells = FindClusters(b_cells,resolution=0.2)

# export for azimuth
# b_cells_azimuth = DietSeurat(b_cells)
# DefaultAssay(b_cells_azimuth) = "RNA"
# b_cells_azimuth[['SCT']] = NULL
# saveRDS(b_cells_azimuth,"bcells_for_azimuth.rds")
# load azimuth annotations back in
preds = read_tsv("azimuth_pred.tsv")
b_cells@meta.data = b_cells@meta.data %>% left_join(preds %>% dplyr::rename("cell_id" = "cell"))
rownames(b_cells@meta.data) = colnames(b_cells)

#filter to just celltypist B cell annotations
b_cells = subset(b_cells,subset = ann_celltypist_highres %in% c("B cells","Cycling B cells","Germinal center B cells","Large pre-B cells","Memory B cells","Naive B cells","Small pre-B cells","Plasma cells","Transitional B cells"))
b_cells = SetIdent(b_cells,value="ann_celltypist_highres")

########################
# basic plots
########################

p1=DimPlot(b_cells)
p2=FeaturePlot(b_cells,features=b_cell_markers)
p3=DotPlot(b_cells,features=b_cell_markers)

png("dimplot.png",res=300,units="in",width=6,height=4)
p1
dev.off()

png("featureplot.png",res=300,units="in",width=5,height=5)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=12,height=4)
p3
dev.off()

##############################
# cluster biomarkers
##############################
biomarkers = FindAllMarkers(b_cells,recorrect_umi=FALSE,logfc.threshold=0,min.pct=0,only.pos=TRUE)
topbiomarkers = biomarkers %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=5)
write_csv(topbiomarkers,"biomarkers.csv")
write_csv(biomarkers,"all_biomarkers.csv")

#######################################
# DA & composition
#######################################
# create new unique ID with donor and source
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))

b_cells@meta.data$cell_type = b_cells@meta.data$ann_celltypist_highres
n_col = b_cells@meta.data$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

b_cells@meta.data$cell_type = factor(
b_cells@meta.data$cell_type,
ordered=TRUE,
levels = c("Naive B cells","Memory B cells","Plasma cells",
"Transitional B cells","B cells","Germinal center B cells",
"Large pre-B cells","Cycling B cells","Small pre-B cells"))

p=DimPlot(b_cells,label=F,raster=F,group.by="cell_type")+
scale_color_manual(values = colour_pal)+
ggtitle("")+
theme_minimal()+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("dim_plot_simple_labels.png",res=300,units="in",width=5,height=4)
p
dev.off()

# stash sample info
sample_info = b_cells@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(b_cells@meta.data$cell_type,b_cells@meta.data$donor_source)

# filter out clusters with 0 counts
abundances = abundances[rowSums(abundances)>10,]

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
  levels = c("Naive B cells","Memory B cells","Plasma cells",
  "Transitional B cells","B cells","Germinal center B cells",
  "Large pre-B cells","Cycling B cells","Small pre-B cells"))

  # da plot

  plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
  theme_classic()+
  scale_color_manual(values = colour_pal)+
  geom_text_repel(size=3,max.overlaps=100)+
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
b_cells@meta.data$phenotype = factor(b_cells@meta.data$phenotype,levels=c("Noninflammatory","OIND","MS"),ordered=TRUE)

# proportion plots
png("./crude_labels_proportion_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(b_cells@meta.data,aes(phenotype,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Phenotype",y="Proportion",fill="Cell type")
dev.off()

plots = list()
for(pheno in c("MS","OIND","Noninflammatory")){
p=ggplot(b_cells@meta.data %>%
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
b_cells@meta.data = b_cells@meta.data %>%
mutate(oligo_pos = ifelse(Oligoclonal %in% c("Pos(>10)","Pos(3to9)"),"OCB+","OCB-")) %>%
mutate(oligo_pos = ifelse(cohort=="EU" & is.na(Oligoclonal),"OCB+",oligo_pos)) %>%
mutate(drb_pos = ifelse(drb1_1501_dose > 0 ,"DRB1*15+","DRB1*15-"))

plots = list()
for(pheno in c("RMS","PPMS")){
p=ggplot(b_cells@meta.data %>%
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
p=ggplot(b_cells@meta.data %>%
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
p=ggplot(b_cells@meta.data %>%
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
b_cells@meta.data = b_cells@meta.data %>%
mutate(oligo_pos = ifelse(
oligo_pos == "OCB+",
"OCBpos",
"OCBneg"
))

b_cells@meta.data = b_cells@meta.data %>%
mutate(drb_pos = ifelse(
drb_pos == "DRB1*15+",
"DRBpos",
ifelse(!is.na(drb_pos),
"DRBneg",
NA
)))

do_da = function(x, contrasts_to_test){
  # stash sample info
  sample_info = b_cells@meta.data %>%
  filter(phenotype=="MS") %>%
  dplyr::select(donor.id,source,x,donor_source) %>%
  filter(!is.na(.data[[x]])) %>%
  distinct(donor_source,.keep_all=TRUE)

  data_for_abundances = b_cells@meta.data %>%
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
  y.ab$samples =  y.ab$samples %>% left_join(b_cells@meta.data %>%
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
    levels = c("Naive B cells","Memory B cells","Plasma cells",
    "Transitional B cells","B cells","Germinal center B cells",
    "Large pre-B cells","Cycling B cells","Small pre-B cells"))

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


p=DimPlot(subset(b_cells,phenotype=="MS" & source=="CSF"),split.by="donor.id",ncol=9)+
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

clusters = b_cells$cell_type %>% unique

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
        for(i in c(1:nrow(sig_pathways))){
          leading_edge = sig_pathways$leadingEdge[[i]]
          path_name = sig_pathways$pathway[i] %>% str_remove("HALLMARK_")
          png(paste0(path_name,"_",cluster,"_",comparison,"gsea_heatmap.png"),res=300,width=8,height=8,units="in")
          DoHeatmap(subset(b_cells, cell_type==cluster & phenotype=="MS"),
          features=leading_edge,
          slot="data",
          group.by="source")
          dev.off()

        }

        write_csv(sig_pathways,file=paste0("sig_pathways_",geneset,"_",cluster,"_",comparison,".csv"))
      }

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
  plot_data$cell_type = factor(plot_data$cell_type,levels=c("Naive B cells","Memory B cells","Plasma cells"),ordered=TRUE)

  png(paste0("summary_",x,".png"),res=300,units="in",width=6,height=6)
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
rownames(b_cells@meta.data) = colnames(b_cells)

# rename NIND (for graphs)
b_cells@meta.data$phenotype = recode(b_cells@meta.data$phenotype, "Noninflammatory" = "NIND")

# define clones
expanded_clones = b_cells@meta.data %>% group_by(donor.id,clone_id) %>% dplyr::count() %>% filter(n>1) %>% mutate(donor_clone = paste0(donor.id,"_",clone_id))
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_clone = paste0(donor.id,"_",clone_id)) %>%
mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

# compare csf vs periphery
png("csf_v_pbmc_celltypes.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
png("csf_v_pbmc_isotypes.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source",group.by="isotype")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

# compare csf in ms vs controls
common_cell_types = c("Naive B cells","Memory B cells","Plasma cells")
ms_v_cont_csf_plot = function(x){
DimPlot(subset(b_cells,source=="CSF" & cell_type %in% common_cell_types),split.by="phenotype",group.by=x)+
scale_color_brewer(palette="Set1")+
theme_minimal()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

p0=ms_v_cont_csf_plot("cell_type")+ggtitle("Cell type")
p1=ms_v_cont_csf_plot("expanded_clone")+ggtitle("Clonal expansion")
png("ms_v_cont_csf_bcr_repertoire.png",res=300,units="in",height=3,width=3)
grid.arrange(p0,p1)
dev.off()

# compare csf vs periphery - just ms
csf_v_pbmc_plot = function(x){
DimPlot(subset(b_cells,phenotype=="MS"),split.by="source",group.by=x)+
scale_color_manual(values = colour_pal)+
theme_minimal()+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
var_list = list("cell_type","expanded_clone","shm_positive","isotype","ighv_family","status_summary")
plot_list = lapply(var_list,csf_v_pbmc_plot)

png("csf_v_pbmc_just_ms_bcr_repertoire.png",res=300,units="in",height=8,width=8)
do.call("grid.arrange",plot_list)
dev.off()


# bar plots
make_categorical_plot = function(x){
  ggplot(b_cells@meta.data,aes(phenotype,fill=b_cells@meta.data[[x]]))+
  geom_bar(position="fill")+
  facet_wrap(~source)+
  theme_bw()+
  labs(x="Phenotype",fill=x,y="Proportion")
}

make_continuous_plot = function(x){
  ggplot(b_cells@meta.data,aes(phenotype,fill=cell_type,y=b_cells@meta.data[[x]]))+
  geom_boxplot()+
  facet_wrap(~source)+
  theme_bw()+
  scale_fill_manual(values = colour_pal)+
  labs(x="Phenotype",y=x)
}

# phenotypes of clonal cells - IGHV types
b_cells@meta.data$ighv_family = sapply(b_cells@meta.data$v_call_genotyped_VDJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  return(y)
})

var_list = list("expanded_clone","shm_positive","isotype","ighv_family","status_summary","c_call_VDJ")
plot_list = lapply(var_list,make_categorical_plot)

png("csf_v_pbmc_bcr_repertoire.png",res=300,units="in",height=10,width=10)
do.call("grid.arrange",plot_list)
dev.off()


# cdr3 length
b_cells@meta.data = b_cells@meta.data %>% mutate(cdr3_length = nchar(junction_aa_VDJ))


#######################################
# Clonal analysis
#######################################


# big picture numbers
table(b_cells@meta.data$expanded_clone)
b_cells@meta.data %>% distinct(donor_clone) %>% nrow
b_cells@meta.data %>% dplyr::count(expanded_clone) %>% mutate (n/sum(n))
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

donors_with_expanded_clones_csf = b_cells@meta.data  %>% filter(source=="CSF") %>% group_by(donor.id) %>% dplyr::count(donor_clone) %>% filter(n>1) %>% distinct(donor.id)
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_with_expanded_clones_csf = ifelse(donor.id %in% donors_with_expanded_clones_csf$donor.id,"Yes","No"))

b_cells@meta.data %>%
  group_by(phenotype) %>%
  distinct(donor.id,.keep_all=T) %>%
  dplyr::count(donor_with_expanded_clones) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(donor_with_expanded_clones=="Yes")

b_cells@meta.data %>%
  group_by(phenotype) %>%
  distinct(donor.id,.keep_all=T) %>%
  dplyr::count(donor_with_expanded_clones_csf) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(donor_with_expanded_clones_csf=="Yes")

# histogram
p=ggplot(b_cells@meta.data,aes(clonal_size,fill=phenotype))+geom_histogram()+theme_bw()+scale_fill_brewer(palette="Set2")+labs(x="Clone size")
png("clonal_histogram.png",res=300,height=2,width=4,units="in")
p
dev.off()

# clonality of csf and pbmc
b_cells@meta.data$expanded_clone = factor(b_cells@meta.data$expanded_clone,levels=c("Not expanded","Expanded"),ordered=TRUE)

p1=ggplot(b_cells@meta.data %>% filter(source=="CSF"),aes(phenotype,fill=expanded_clone))+
geom_bar(position="fill")+
scale_fill_brewer(palette="Set2")+
theme_bw()+
labs(fill="Expanded clone?",x="Phenotype",y="Proportion of B cell pool")+
png("clonal_proportions_source.png",res=300,height=4,width=4,units="in")
p1
dev.off()


# quantify clonality of csf
counts = b_cells@meta.data %>% group_by(source,phenotype) %>% dplyr::count(expanded_clone)
totals = b_cells@meta.data %>% group_by(source,phenotype) %>% dplyr::count(expanded_clone) %>% summarise(totals = sum(n))
counts %>% left_join(totals,by=c("source","phenotype")) %>% mutate(prop = n/totals)

# per individual
counts = b_cells@meta.data %>% group_by(source,phenotype,donor.id) %>% dplyr::count(expanded_clone)
totals = b_cells@meta.data %>% group_by(source,phenotype,donor.id) %>% dplyr::count() %>% mutate(totals = n) %>% dplyr::select(-n)
props = counts %>% left_join(totals,by=c("source","phenotype","donor.id")) %>% mutate(prop = n/totals) %>%
filter(expanded_clone == "Expanded") %>%
mutate(donor_source = paste0(donor.id,"_",source)) %>%
filter(source == "CSF")

# add in 0 counts
all_donors = b_cells@meta.data %>% distinct(source,phenotype,donor.id) %>%
mutate(donor_source = paste0(donor.id,"_",source))
props = props %>%
bind_rows(all_donors %>% filter(!donor_source %in% props$donor_source))
props$prop = tidyr::replace_na(props$prop,0)

props %>% group_by(source,phenotype) %>% summarise(mean(prop))

#######################################
# look at relation to phenotypes
#######################################
phenotypes_clones_df = props %>% left_join(b_cells@meta.data %>% distinct(donor.id,.keep_all=TRUE),by="donor.id")
phenotypes_clones_df$phenotype.x = factor(
phenotypes_clones_df$phenotype.x,
levels = c("NIND","OIND","MS"),
ordered=T
)

get_t_test = function(pheno1,pheno2){
a = phenotypes_clones_df[phenotypes_clones_df$source.x=="CSF" & phenotypes_clones_df$phenotype.x==pheno1,][['prop']]
b = phenotypes_clones_df[phenotypes_clones_df$source.x=="CSF" & phenotypes_clones_df$phenotype.x==pheno2,][['prop']]
t.test(a,b)$p.value
}
phenos = c("MS","NIND","OIND")
pheno_df = expand.grid(phenos,phenos) %>% filter(!Var1 == Var2)
pheno_df$pvals = mapply(get_t_test,as.character(pheno_df$Var1),as.character(pheno_df$Var2))

p0 = ggplot(phenotypes_clones_df,aes(phenotype.x,prop,fill=source.x,col))+geom_boxplot(alpha=0.8)+geom_point(alpha=0.7,position=position_dodge(width=0.75))+theme_minimal()+scale_fill_brewer(palette="Set1")+labs(x="Phenotype",y="Proportion of \nexpanded B cells",fill="Compartment")
png("clonal_proportions_individual_level.png",res=300,height=4,width=4,units="in")
p0
dev.off()

# plot per individual
rownames(b_cells@meta.data) = colnames(b_cells)

p=DimPlot(subset(b_cells,source=="CSF" & phenotype=="MS"),group.by="expanded_clone",split.by="donor.id",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())
png("./indiv_clonal_dim_plots_ms.png",res=300,units="in",height=8,width=8)
p
dev.off()

# plot vs age

model = glm(data = phenotypes_clones_df,
prop ~ factor(phenotype.x,ordered=F) + Age) %>% summary
model$coefficients[4,4]
pval = round(model$coefficients[4,4],2)

p0=ggplot(phenotypes_clones_df,aes(Age,prop,col=phenotype.x))+
geom_point(alpha=0.8)+
geom_smooth(method="lm",se=F)+
theme_minimal()+
scale_color_brewer(palette="Set1")+
labs(x="Age",y="Proportion of \nexpanded B cells",color="Phenotype")+
annotate("label",label = paste0("P=",pval),x=Inf,y=Inf,vjust=1,hjust=1)

png("clonal_proportions_individual_level_vs_age.png",res=300,height=4,width=4,units="in")
p0
dev.off()

# repeat with periphery (sense check )

# per individual
counts = b_cells@meta.data %>% group_by(source,phenotype,donor.id) %>% dplyr::count(expanded_clone)
totals = b_cells@meta.data %>% group_by(source,phenotype,donor.id) %>% dplyr::count() %>% mutate(totals = n) %>% dplyr::select(-n)
props = counts %>% left_join(totals,by=c("source","phenotype","donor.id")) %>% mutate(prop = n/totals) %>%
filter(expanded_clone == "Expanded") %>%
mutate(donor_source = paste0(donor.id,"_",source)) %>%
filter(source == "PBMC")

# add in 0 counts
all_donors = b_cells@meta.data %>% distinct(source,phenotype,donor.id) %>%
mutate(donor_source = paste0(donor.id,"_",source))
props = props %>%
bind_rows(all_donors %>% filter(!donor_source %in% props$donor_source))
props$prop = tidyr::replace_na(props$prop,0)

props %>% group_by(source,phenotype) %>% summarise(mean(prop))





phenotypes_clones_df = phenotypes_clones_df %>% mutate(oligo_pos = ifelse(Oligoclonal %in% c("Pos(>10)","Pos(3to9)"),"OCB+","OCB-"))

ggplot(phenotypes_clones_df %>%
filter(Category_fine %in% c("RMS","PPMS")),
aes(Category_fine,prop))+
facet_wrap(~oligo_pos)+
geom_jitter()+
geom_boxplot()+
theme_minimal()

ggplot(phenotypes_clones_df %>%
filter(Category_fine %in% c("RMS","PPMS","OIND","Noninflammatory")),
aes(Category_fine,fill=donor_with_expanded_clones_csf))+
facet_wrap(~oligo_pos)+
geom_bar(position="fill")+
theme_minimal()


ggplot(phenotypes_clones_df %>%
filter(Category_fine %in% c("RMS","PPMS")),
aes(oligo_pos,prop))+
geom_jitter()+
geom_boxplot()+
theme_minimal()

get_t_test = function(pheno1,pheno2){
a = phenotypes_clones_df[phenotypes_clones_df$source.x=="CSF" & phenotypes_clones_df$phenotype.x=="MS" & phenotypes_clones_df$oligo_pos=="OCB-",][['prop']]
b = phenotypes_clones_df[phenotypes_clones_df$source.x=="CSF" & phenotypes_clones_df$phenotype.x=="MS"& phenotypes_clones_df$oligo_pos=="OCB+",][['prop']]
t.test(a,b)$p.value
}


ggplot(phenotypes_clones_df %>%
filter(Category_fine %in% c("RMS","PPMS")),
aes(oligo_pos,fill=donor_with_expanded_clones_csf))+
geom_bar(position="fill")+
theme_minimal()

phenotypes_clones_df = phenotypes_clones_df %>% mutate(drb_pos = ifelse(drb1_1501_dose > 0 ,"DRB1*15+","DRB1*15-"))
ggplot(phenotypes_clones_df %>%
filter(!is.na(drb_pos) && phenotype.x=="MS"),
aes(factor(drb_pos),prop))+
geom_jitter()+
geom_boxplot()+
theme_minimal()

ggplot(phenotypes_clones_df %>%
filter(!is.na(drb_pos) && phenotype.x=="MS"),
aes(factor(drb_pos),fill=donor_with_expanded_clones_csf))+
geom_bar(position="fill")+
theme_minimal()

get_t_test = function(pheno1,pheno2){
a = phenotypes_clones_df[phenotypes_clones_df$source.x=="CSF" & phenotypes_clones_df$phenotype.x=="MS" & phenotypes_clones_df$drb_pos=="DRB1*15-",][['prop']]
b = phenotypes_clones_df[phenotypes_clones_df$source.x=="CSF" & phenotypes_clones_df$phenotype.x=="MS"& phenotypes_clones_df$drb_pos=="DRB1*15+",][['prop']]
t.test(a,b)$p.value
}

ggplot(phenotypes_clones_df %>%
filter(Category_fine %in% c("RMS","PPMS")),
aes(Age,prop,color=Category_fine))+
geom_point()+
theme_minimal()

ggplot(phenotypes_clones_df %>%
filter(Category_fine %in% c("RMS","PPMS")),
aes(drb_pos,fill=donor_with_expanded_clones_csf))+
geom_bar(position="fill")+
theme_minimal()

#######################################
# Clonal phenotypes
#######################################

rownames(b_cells@meta.data) = colnames(b_cells)
# compare expanded vs non-expanded in MS
common_cell_types = c("Naive B cells","Memory B cells","Plasma cells")
expanded_v_not_plot = function(x){
DimPlot(subset(b_cells,phenotype=="MS" & cell_type %in% common_cell_types),
split.by="expanded_clone",group.by=x)+
scale_color_brewer(palette="Set1")+
theme_minimal()+
theme(axis.text.y=element_blank(),
axis.text.x=element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank())
}

p0=expanded_v_not_plot("cell_type")+ggtitle("Cell type")
p1=expanded_v_not_plot("isotype")+ggtitle("Isotype")
p2=expanded_v_not_plot("shm_positive")+ggtitle("SHM")
p3=expanded_v_not_plot("ighv_family")+ggtitle("IGHV")
png("expanded_v_not_bcrs_just_ms.png",res=300,units="in",height=5,width=7)
grid.arrange(p0,p1,p2,p3)
dev.off()


png("clones_celltypes.png",res=300,height=3,width=6,units="in")
DimPlot(subset(b_cells, subset = phenotype=="MS"),split.by="expanded_clone")
dev.off()
good_cells =  b_cells@meta.data %>% dplyr::count(predicted.celltype.l3) %>% filter(n>20)
png("clones_celltypes2.png",res=300,height=3,width=6,units="in")
DimPlot(subset(b_cells, subset = phenotype=="MS" & predicted.celltype.l3 %in% good_cells$predicted.celltype.l3 ),group.by="predicted.celltype.l3",split.by="expanded_clone")
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

# fidelity within each clone
expanded = b_cells@meta.data %>% filter(phenotype=="MS" & expanded_clone=="Expanded")

same_cell_lgl = list()
same_isotype_lgl = list()
same_source_lgl = list()
same_shm_lgl = list()

for(i in c(1:length(unique(expanded$donor_clone)))){
  this_clone = unique(expanded$donor_clone)[i]
  this_clone_data = expanded %>% filter(donor_clone == this_clone)
  same_cell = all(this_clone_data$cell_type[1]==this_clone_data$cell_type)
  same_isotype = all(this_clone_data$isotype[1]==this_clone_data$isotype)
  same_source = all(this_clone_data$source[1]==this_clone_data$source)
  same_shm = all(this_clone_data$shm_positive[1]==this_clone_data$shm_positive)
  same_cell_lgl[i] = same_cell
  same_isotype_lgl[i] = same_isotype
  same_source_lgl[i] = same_source
  same_shm_lgl[i] = same_shm
}

clone_df = data.frame(
donor_clone = unique(expanded$donor_clone),
same_cell = unlist(same_cell_lgl),
same_isotype = unlist(same_isotype_lgl),
same_source = unlist(same_source_lgl),
same_shm = unlist(same_shm_lgl))
table(clone_df$same_cell)[1]/nrow(clone_df)*100
table(clone_df$same_isotype)[1]/nrow(clone_df)*100
table(clone_df$same_source)[1]/nrow(clone_df)*100

discordant_clones = clone_df %>% filter(same_cell==FALSE)
p0=ggplot(expanded %>% filter(donor_clone %in% discordant_clones$donor_clone),aes(donor_clone,fill=cell_type))+geom_bar(position="fill")+coord_flip()+theme_minimal()+scale_fill_brewer(palette="Paired")+labs(y="Proportion of clone",x="Clone ID",fill="Cell type")
discordant_clones2 = clone_df %>% filter(same_isotype==FALSE)
p1=ggplot(expanded %>% filter(donor_clone %in% discordant_clones2$donor_clone),
aes(donor_clone,fill=isotype))+geom_bar(position="fill")+coord_flip()+theme_minimal()+scale_fill_brewer(palette="Paired")+labs(y="Proportion of clone",x="Clone ID",fill="Isotype")
discordant_clones3 = clone_df %>% filter(same_source==FALSE)
p2=ggplot(expanded %>% filter(donor_clone %in% discordant_clones3$donor_clone),
aes(donor_clone,fill=source))+geom_bar(position="fill")+coord_flip()+theme_minimal()+scale_fill_brewer(palette="Paired")+labs(y="Proportion of clone",x="Clone ID",fill="Compartment")

png("clonal_phenotypes.png",res=300,height=8,width=8,units="in")
lay = rbind(
c(1,1,1,1,1,1,1,1,1,1),
c(1,1,1,1,1,1,1,1,1,1),
c(2,2,2,2,2,2,2,2,2,NA),
c(3,3,3,3,3,3,3,3,3,NA)
)
grid.arrange(p0,p1,p2,layout_matrix = lay)

dev.off()

# exemplar clone
plot_clone = expanded %>% filter(donor_clone %in% discordant_clones$donor_clone)


clone_plot = function(clone_to_plot){
plot_clone_data = plot_clone %>%
dplyr::count(donor_clone,phenotype,source,cell_type,isotype,shm_positive) %>%
filter(donor_clone == clone_to_plot)

png(paste0(clone_to_plot,"clonal_heterogeneity_cell_types.png"),res=300,units="in",width=4,height=3)
p=ggplot(plot_clone_data,aes(source,cell_type,color=isotype,size=n,alpha=0.8))+
geom_point(position=position_dodge(width=1))+
theme_minimal()+
guides(alpha=FALSE)+
geom_vline(xintercept=1.5,alpha=0.1)+
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
ggtitle(plot_clone_data$donor_clone[1])
print(p)
dev.off()
}

clone_plot("9X14GT24_77_9_3_449")
clone_plot("TU695_450_11_3_28")
clone_plot("TU682_117_16_1_149")




#################################
# check for public clones
################################
public_clones = b_cells@meta.data %>% group_by(clone_id) %>% dplyr::count(donor.id)  %>% dplyr::count(clone_id) %>% arrange(desc(n)) %>% filter(n>1)

b_cells@meta.data %>% filter(clone_id %in% public_clones$clone_id)

bcr_db = read_csv("bcr_db.csv")
b_cells@meta.data %>% filter(junction_aa_VDJ %in% bcr_db$CDR3.heavy.aa)

# lanz
lanz_cdr3 = read_csv("../Lanz_AllCsfCdr3ForBenJacobs.csv")
b_cells@meta.data %>% filter(junction_aa_VDJ %in% lanz_cdr3$`JUNCTION | HC`)
b_cells@meta.data %>% filter(junction_aa_VJ %in% lanz_cdr3$`JUNCTION | LC`)

# gasperi
gasperi = read_delim("bcr_data_02052022.csv",delim=";")

b_cells@meta.data %>% filter(junction_aa_VDJ %in% gasperi$junction_aa_VDJ)

################################
# de clonal vs not
################################

DefaultAssay(b_cells) = "RNA"
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_expanded = paste0(donor.id,"_",expanded_clone))
rownames(b_cells@meta.data) = colnames(b_cells)

cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == "CSF" & cell_type == "Plasma cells")
# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 2

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
ggtitle("Plasma cells")+
scale_color_manual(values = colours)+
scale_x_continuous(limits=c(-7,7))+
scale_y_continuous(limits=c(0,10))

png("pcs_expanded_vs_not_csf_pcs_ms.png",res=300,units="in",height=6,width=6)
print(plot)
dev.off()

# gsea
de = res
ranked_genes = de %>% arrange(desc(logFC)) %>% dplyr::select(gene,logFC)
ranked_genes_vector = ranked_genes$logFC
names(ranked_genes_vector) = ranked_genes$gene
message("There are ",nrow(ranked_genes)," genes in this analysis")
write_csv(res,file="pcs_expanded_vs_not_csf_ms.csv")

# do gsea
nperm=1000000
geneset="hallmark"
genelist = eval(parse(text = paste0(geneset,"_list")))
res = fgsea(genelist, stats = ranked_genes_vector, minSize=10,eps=0, nPermSimple = nperm)
res = res %>% arrange(NES)
res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
topres = res %>% arrange(pval) %>% head(n=20) %>% arrange(desc(NES))
topres$pathway = str_remove(topres$pathway,"HALLMARK_")
topres$pathway = factor(topres$pathway,ordered=TRUE,levels=topres$pathway)
topres$direction = ifelse(topres$NES>0,"Up","Down")

colours = c("Up" = "red", "Down" = "blue")
# NB labeled with nominal significance
plot = ggplot(topres,aes(NES,pathway,label=ifelse(pval<0.1,"*"," "),fill=direction))+
  theme_bw()+
  labs(x="Normalised enrichment score",fill="Pathway enrichment")+
  geom_col(color="black")+
  geom_text(size=5,position = position_dodge(0.5))+
  scale_fill_manual(values = colours,labels = c("Up in expanded clones","Down in expanded clones"))+
  guides(alpha=FALSE)+
  theme_minimal()

leading_edge_protein = res[res$pathway=="HALLMARK_PROTEIN_SECRETION",'leadingEdge']$leadingEdge %>% unlist()

png("protein_secretion_heatmap.png",res=300,width=8,height=8,units="in")
DoHeatmap(subset(b_cells, cell_type=="Plasma cells" & phenotype=="MS" & source=="CSF"),
features=leading_edge_protein,
slot="data",
group.by="expanded_clone")
dev.off()

png("pcs_csf_gsea_results_hallmark_expanded.png",res=300,units="in",height=4,width=7)
print(plot)
dev.off()

################################
# pathway scores
################################
protein_sec = res[res$pathway=="HALLMARK_PROTEIN_SECRETION",'leadingEdge']$leadingEdge %>% unlist()
glyc = res[res$pathway=="HALLMARK_GLYCOLYSIS",'leadingEdge']$leadingEdge %>% unlist()
oxphos = res[res$pathway=="HALLMARK_OXIDATIVE_PHOSPHORYLATION",'leadingEdge']$leadingEdge %>% unlist()
xeno = res[res$pathway=="HALLMARK_XENOBIOTIC_METABOLISM",'leadingEdge']$leadingEdge %>% unlist()
upr = res[res$pathway=="HALLMARK_UNFOLDED_PROTEIN_RESPONSE",'leadingEdge']$leadingEdge %>% unlist()
mtorc1 = res[res$pathway=="HALLMARK_MTORC1_SIGNALING",'leadingEdge']$leadingEdge %>% unlist()
clonal_sig = c(protein_sec,glyc,oxphos,xeno,upr,mtorc1)


# just look at PCs
pcs = subset(b_cells,subset = cell_type =="Plasma cells")
# recluster
set.seed(1)
DefaultAssay(pcs)="SCT"
pcs = RunUMAP(pcs,reduction="harmony",dims=1:50)
pcs = FindNeighbors(pcs)
pcs = FindClusters(pcs,resolution=0.2)

add_score = function(genelist,name = "gene_score"){
features = list(genelist)
# module score
DefaultAssay(pcs) = "RNA"
pcs = AddModuleScore(pcs,features = features,nbin=10,ctrl=1000,name=name)
pcs[[paste0(name,"_z")]] = RNOmni::RankNorm(pcs@meta.data[[paste0(eval(name),"1")]])
pcs
}

pcs = add_score(clonal_sig)

# stats
a = pcs@meta.data %>% filter(phenotype=="MS") %>% filter(source == "CSF"  & cell_type=="Plasma cells" & expanded_clone=="Not expanded")
b = pcs@meta.data %>% filter(phenotype=="MS") %>% filter(source == "PBMC" & cell_type=="Plasma cells" & expanded_clone=="Not expanded")
t.test(a$gene_score_z,b$gene_score_z)$p.value

a = pcs@meta.data %>% filter(phenotype=="MS") %>% filter(source == "CSF"  & cell_type=="Plasma cells" & expanded_clone=="Not expanded")
b = pcs@meta.data %>% filter(phenotype=="MS") %>% filter(source == "CSF" & cell_type=="Plasma cells" & expanded_clone=="Expanded")
t.test(a$gene_score_z,b$gene_score_z)$p.value

png("bcells_annotations_CSF.png",res=300,units="in",width=4,height=4)
FeaturePlot(subset(pcs, source == "CSF"),features="gene_score_z")+ggtitle("CSF")+scale_colour_gradient(low="darkblue",high="orange",limits=c(-4,4))
dev.off()
png("bcells_annotations_PBMC.png",res=300,units="in",width=4,height=4)
FeaturePlot(subset(pcs, source == "PBMC"),features="gene_score_z")+ggtitle("PBMC")+scale_colour_gradient(low="darkblue",high="orange",limits=c(-4,4))
dev.off()

png("clonal_sig_gene_score_vln.png",res=300,units="in",width=4,height=4)
VlnPlot(subset(pcs,cell_type=="Plasma cells"),features="gene_score_z",group.by="source",split.by="expanded_clone")+ggtitle("Clonal \ngene score")+scale_fill_brewer(palette="Set3")+
stat_summary(fun="median",col="red",position=position_dodge(width=0.9))
dev.off()

png("mtorc1_sig_gene_score_vln.png",res=300,units="in",width=4,height=4)
VlnPlot(subset(pcs, subset = expanded_clone=="Not expanded"),features="mtorc1_z",group.by="source")+ggtitle("MTORC1 \ngene score")+scale_fill_brewer(palette="Set3")+labs(x="Compartment",y="MTORC1 gene score")+
stat_summary(fun="median",col="red",position=position_dodge(width=0.9))+
annotate(geom="segment",x=1,xend=2,y=3.5,yend=3.5)+
annotate(geom="text",x=1.5,y=3.8,label = "P < 0.0001")+
scale_y_continuous(limits=c(-4,4))
dev.off()


png("mtorc1_sig_gene_score_vln_all_cells.png",res=300,units="in",width=4,height=4)
VlnPlot(subset(pcs,phenotype=="MS"),features="mtorc1_z",split.by="source",group.by="cell_type")+
ggtitle("MTORC1 \ngene score")+
scale_fill_brewer(palette="Set3")+
labs(x="Cell type",fill="Compartment",y="MTORC1 gene score")+
stat_summary(fun="median",col="red",position=position_dodge(width=0.9))+
scale_y_continuous(limits=c(-4,4))
dev.off()



FeatureScatter(pcs,feature1 = "mtorc1_z", feature2 = "clonal_size",group.by="cell_type")

################################
# de clonal vs not in single donors
################################


overall_results_df = data.frame()
# n=1 clonal vs not
DefaultAssay(b_cells) = "RNA"
# get names for clusters
clusters = "Plasmablast"


donor_list = unique(b_cells@meta.data$donor.id)

do_de_per_donor = function(donor, cell_type = "Plasmablast", source = "CSF"){
  message("Donor: ",donor)
  message("Source: ",source)
  message("Cell type:", cell_type)
  cells = b_cells@meta.data %>% filter(phenotype=="MS" & source == source & donor.id ==donor & cell_type == cell_type)
  no_expanded = cells %>% filter(expanded_clone=="Expanded")
  non_expanded = cells %>% filter(expanded_clone=="Not expanded")
  n_clones = no_expanded %>% distinct(clone_id) %>% nrow
  if(nrow(cells)==0 | nrow(no_expanded)<10 | nrow(non_expanded)==0){
    return(NA)
  }
  cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == source & donor.id ==donor & cell_type == cell_type)

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

sig_results = overall_results_df %>% filter(!grepl("IG",gene))

colours = c("Up" = "red","Down" = "blue", "neither" = "grey")
plot = ggplot(sig_results,aes(logFC,-log10(PValue),color=direction,label=gene))+
theme_bw()+
facet_wrap(~donor_clone)+
geom_point()+
NoLegend()+
geom_text_repel(data=sig_results %>% filter(P_adj<0.05),mapping=aes())+
scale_color_manual(values = colours)

png("pcs_clonal_vs_not_csf_single_donors_no_ig.png",res=300,units="in",height=16,width=16)
plot
dev.off()

# look at biggest donor
donor = subset(pcs,subset = donor.id == "TUM_SC_13" & source=="CSF")
sig_results = overall_results_df %>% filter(P_adj < 0.05) %>% filter(!grepl("IG",gene))

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



drb = b_cells@meta.data %>% filter(!is.na(drb1_1501_dose))
counts = drb %>% group_by(donor.id,phenotype,source,drb1_1501_dose) %>% dplyr::count()
expanded = drb %>% group_by(donor.id,phenotype,source,expanded_clone,drb1_1501_dose) %>% dplyr::count() %>% filter(expanded_clone=="Expanded")
combo = expanded %>% left_join(counts,by=c("donor.id","source","phenotype"))
combo = combo %>% mutate(proportion_expanded = n.x/n.y)

ggplot(combo %>% filter(phenotype=="MS" & source=="CSF"),aes(factor(drb1_1501_dose.x),proportion_expanded))+geom_boxplot()
ggplot(expanded %>% filter(phenotype=="MS" & source=="CSF"),aes(drb1_1501_dose,n))+geom_point()

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

cells_for_de = subset(b_cells, subset = phenotype=="MS" & source == "CSF" & cell_type == cluster)

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
VlnPlot(subset(genotyped_cohort,subset = phenotype=="MS" & source=="CSF" & cell_type=="Plasmablast"),features=c("EAF2","CD86"),group.by="rs2331964_T")


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
