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

# filter to just celltypist B cell annotations
b_cells = subset(b_cells,subset = ann_celltypist_highres %in% c("B cells","Cycling B cells","Germinal center B cells","Large pre-B cells","Memory B cells","Naive B cells","Small pre-B cells","Plasma cells","Transitional B cells"))
b_cells = SetIdent(b_cells,value="ann_celltypist_highres")

# plots to check annotations
p1=DimPlot(b_cells)
p2=FeaturePlot(b_cells,features=b_cell_markers)
p3=DotPlot(b_cells,features=b_cell_markers)

png("dimplot.png",res=300,units="in",width=12,height=12)
p1
dev.off()

png("featureplot.png",res=300,units="in",width=12,height=12)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=10,height=12)
p3
dev.off()

##############################
# cluster biomarkers
##############################
#biomarkers = FindAllMarkers(b_cells,recorrect_umi=FALSE,logfc.threshold=0,min.pct=0,only.pos=TRUE)
#topbiomarkers = biomarkers %>% group_by(cluster) %>% slice_min(order_by=p_val_adj,n=5)
#write_csv(topbiomarkers,"biomarkers.csv")
#write_csv(biomarkers,"all_biomarkers.csv")

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
clusters = c("Plasma cells","Memory B cells","Naive B cells")

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

## define variables
# rename NIND (for graphs)
b_cells@meta.data$phenotype = recode(b_cells@meta.data$phenotype, "Noninflammatory" = "NIND")

# define clones
expanded_clones = b_cells@meta.data %>% group_by(donor.id,clone_id) %>% dplyr::count() %>% filter(n>1) %>% mutate(donor_clone = paste0(donor.id,"_",clone_id))
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_clone = paste0(donor.id,"_",clone_id)) %>%
mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

# define IGHV types
b_cells@meta.data$ighv_family = sapply(b_cells@meta.data$v_call_genotyped_VDJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  return(y)
})

# cdr3 length
b_cells@meta.data = b_cells@meta.data %>% mutate(cdr3_length = nchar(junction_aa_VDJ))

## plots

# cell types
# compare csf vs periphery
png("csf_v_pbmc_celltypes.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

## isotypes
png("csf_v_pbmc_isotypes.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source",group.by="isotype")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

## ighv
png("csf_v_pbmc_ighv.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source",group.by="ighv_family")+scale_color_manual(values = colour_pal)+
theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
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
  geom_bar(position="fill",color="black")+
  facet_wrap(~source)+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  labs(x="Phenotype",fill=x,y="Proportion")
}

make_continuous_plot = function(x){
  ggplot(b_cells@meta.data,aes(phenotype,fill=cell_type,y=b_cells@meta.data[[x]]))+
  geom_boxplot(color="black")+
  facet_wrap(~source)+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  labs(x="Phenotype",y=x)
}

var_list = list("expanded_clone","shm_positive","isotype","ighv_family","status_summary","c_call_VDJ")
plot_list = lapply(var_list,make_categorical_plot)

png("csf_v_pbmc_bcr_repertoire.png",res=300,units="in",height=10,width=10)
do.call("grid.arrange",plot_list)
dev.off()

#######################################
#  IGHV usage
#######################################
cells = c("Naive B cells","Memory B cells","Plasma cells")

do_da_var = function(x, contrasts_to_test,variable,ylim=10, data = b_cells@meta.data,plot_title){

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
variable="ighv_family",
plot_title = "all",
ylim=20)

# repeat, sampling each clone only once
do_da_var(data = b_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
,
x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "NIND_CSF - NIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_CSF - NIND_CSF"),
variable="ighv_family",
plot_title = "eachclone",
ylim=20)

# individual gene segments
do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF"),
variable="v_call_genotyped_VDJ",
plot_title = "all",
ylim=20)

# light chain
do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF"),
variable="v_call_genotyped_VJ",
plot_title = "lightchain",
ylim=20)

do_da_var(data = b_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
,
x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF"),
variable="v_call_genotyped_VJ",
plot_title = "lightchain_eachclone",
ylim=20)


# cell-type specific
do_da_var(x = "phenotype",
data = b_cells@meta.data %>% filter(cell_type == "Plasma cells")
,
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_PBMC - NIND_PBMC"),
variable="ighv_family",
plot_title = "plasmacells",
ylim=20)

# cell-type specific
do_da_var(x = "phenotype",
data = b_cells@meta.data %>% filter(cell_type == "Memory B cells")
,
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF",
    "MS_PBMC - NIND_PBMC"),
variable="ighv_family",
plot_title = "memb",
ylim=20)

# cell-type specific
do_da_var(x = "phenotype",
data = b_cells@meta.data %>% filter(cell_type == "Naive B cells")
,
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_PBMC - NIND_PBMC"),
variable="ighv_family",
plot_title = "naive",
ylim=20)

#######################################
#  Isotypes
#######################################

do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "NIND_CSF - NIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_CSF - NIND_CSF"),
variable="c_call_VDJ",
plot_title = "isotypes",
ylim=50)

do_da_var(x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "NIND_CSF - NIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_CSF - NIND_CSF"),
variable="c_call_VJ",
plot_title = "isotypes",
ylim=50)

# cell types
do_da_var(x = "phenotype",
data = b_cells@meta.data %>% filter(cell_type == "Naive B cells")
,
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "NIND_CSF - NIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF",
  "MS_CSF - NIND_CSF"),
variable="c_call_VDJ",
plot_title = "naive",
ylim=50)

# cell types
do_da_var(x = "phenotype",
data = b_cells@meta.data %>% filter(cell_type == "Plasma cells")
,
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF"),
variable="c_call_VDJ",
plot_title = "plasmacells",
ylim=50)

# cell types
do_da_var(x = "phenotype",
data = b_cells@meta.data %>% filter(cell_type == "Memory B cells")
,
contrasts_to_test = c("MS_CSF - MS_PBMC",
  "OIND_CSF - OIND_PBMC",
  "MS_PBMC - OIND_PBMC",
  "MS_PBMC - NIND_PBMC",
  "MS_CSF - OIND_CSF"),
variable="c_call_VDJ",
plot_title = "memb",
ylim=50)

#####################################
# shm
#####################################

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
  props = b_cells@meta.data %>%
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

binary_comparison(variable = "shm_positive",plot_title="SHM")

# repeat, cell-type specific
binary_comparison_celltypes = function(variable,plot_title,level="Yes"){
  props = b_cells@meta.data %>%
  filter(cell_type %in% c("Naive B cells","Plasma cells","Memory B cells")) %>%
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

# counts
b_cells@meta.data %>% group_by(cell_type,source) %>% dplyr::count(shm_positive) %>% mutate(prop = n/sum(n)*100) %>% filter(shm_positive=="Yes")
binary_comparison_celltypes(variable = "shm_positive",plot_title="SHM")

##############################
# mutational load
##############################

do_cont_comparison = function(variable, plot_title){
  pvals = list()

  dat = b_cells@meta.data %>%
    filter(phenotype=="MS" & cell_type %in% cells)
  for(i in c(1:length(cells))){
    a = dat[dat$source=="CSF" & dat$cell_type==cells[i],][[variable]]
    b = dat[dat$source=="PBMC" & dat$cell_type==cells[i],][[variable]]
    pvals[[i]] = t.test(a,b)$p.value
  }
  pvals = data.frame(cell_type = cells,P = unlist(pvals))

  dat = dat %>%
    left_join(pvals,by="cell_type") %>%
    mutate(P = simplify_pval(P))

  p=ggplot(dat,aes(cell_type,.data[[variable]],fill=source))+
    geom_boxplot(color="black")+
    scale_fill_brewer(palette="Set1")+
    ggtitle(plot_title)+
    geom_text(mapping = aes(x = cell_type,y = 1,label = P))+
    theme_minimal()+
    labs(x="Cell type",y="Mutational load")


  png(file=paste0(plot_title,"_comparison.png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
}
do_cont_comparison(variable = "heavychain_mu_freq_cdr_r",plot_title = "Replacement mutations")
do_cont_comparison(variable = "heavychain_mu_freq_cdr_s",plot_title = "Silent mutations")


do_cont_comparison(variable = "cdr3_length",plot_title = "CDR3 length")

# repeat by isotype
do_cont_comparison_by_isotype = function(variable, plot_title){
  pvals = list()
  dat = b_cells@meta.data %>%
    filter(phenotype=="MS")
  isotypes = unique(dat$isotype)

  for(i in c(1:length(isotypes))){
    a = dat[dat$source=="CSF" & dat$isotype==isotypes[i],][[variable]]
    b = dat[dat$source=="PBMC" & dat$isotype==isotypes[i],][[variable]]
    pvals[[i]] = t.test(a,b)$p.value
  }
  pvals = data.frame(isotype = isotypes,P = unlist(pvals))

  dat = dat %>%
    left_join(pvals,by="isotype") %>%
    mutate(P = simplify_pval(P))

  p=ggplot(dat,aes(isotype,.data[[variable]],fill=source))+
    geom_boxplot(color="black")+
    scale_fill_brewer(palette="Set1")+
    ggtitle(plot_title)+
    geom_text(mapping = aes(x = isotype,y = 1,label = P))+
    theme_minimal()+
    labs(x="Cell type",y="Mutational load")


  png(file=paste0(plot_title,"_comparison.png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
}
do_cont_comparison_by_isotype(variable = "heavychain_mu_freq_cdr_r",plot_title = "Replacement mutations (by isotype)")
do_cont_comparison_by_isotype(variable = "heavychain_mu_freq_cdr_s",plot_title = "Silent mutations (by isotype)")
do_cont_comparison_by_isotype(variable = "cdr3_length",plot_title = "CDR3 length (by isotype)")


#######################################
# Clonal analysis
#######################################

binary_comparison(variable = "expanded_clone",plot_title="Clonal expansion",level="Expanded")
binary_comparison_celltypes(variable = "expanded_clone",plot_title="Clonal expansion (by cell type)",level="Expanded")

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

# get some counts
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


#######################################
# look at relation to phenotypes
#######################################
dat = b_cells@meta.data %>%
  filter(cell_type %in% cells) %>%
  group_by(cell_type,donor.id,source,phenotype,drb_pos,oligo_pos,Category_fine,Age,Gender) %>%
  dplyr::count(expanded_clone) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(total>1 & expanded_clone=="Expanded")

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
rownames(b_cells@meta.data) = colnames(b_cells)
p=DimPlot(subset(b_cells,source=="CSF" & phenotype=="MS"),group.by="expanded_clone",split.by="donor.id",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())
png("./indiv_clonal_dim_plots_ms.png",res=300,units="in",height=8,width=8)
p
dev.off()

# models
model_csf = glm(data = dat %>% filter(source=="CSF"),
prop ~ Age + Gender) %>% summary
model_pbmc = glm(data = dat %>% filter(phenotype=="MS" & source=="PBMC"),
prop ~ Age + Gender) %>% summary

p0=ggplot(dat %>% filter(phenotype=="MS"),aes(Age,prop,col=source))+
geom_point(alpha=0.8)+
theme_minimal()+
facet_wrap(~source)+
geom_smooth(method="lm",alpha=0.1,width=0.1,se=F)+
scale_color_brewer(palette="Set1")+
labs(x="Age",y="Proportion of \nexpanded B cells",color="Source")

p1=ggplot(dat %>% filter(phenotype=="MS") %>%
mutate(Gender = ifelse(Gender == 1,"Male","Female")),
aes(Gender,prop,fill=source))+
geom_boxplot(alpha=0.8)+
theme_minimal()+
scale_color_brewer(palette="Set1")+
labs(x="Gender",y="Proportion of \nexpanded B cells",fill="Source")

png("clonal_proportions_individual_level_vs_age.png",res=300,height=4,width=8,units="in")
grid.arrange(p0,p1)
dev.off()



variable = "oligo_pos"
plot_title = "DRB status"
pheno_corr = function(variable,plot_title){
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

pheno_corr("drb_pos","DRB status")
pheno_corr("oligo_pos","OCB status")
pheno_corr("Category_fine","MS subtype")

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

expanded_v_not_plot_csf = function(x){
DimPlot(subset(b_cells,phenotype=="MS" & cell_type %in% common_cell_types & source=="CSF"),
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

p0=expanded_v_not_plot_csf("cell_type")+ggtitle("Cell type")
p1=expanded_v_not_plot_csf("isotype")+ggtitle("Isotype")
p2=expanded_v_not_plot_csf("shm_positive")+ggtitle("SHM")
p3=expanded_v_not_plot_csf("ighv_family")+ggtitle("IGHV")
png("expanded_v_not_bcrs_just_ms_just_csf.png",res=300,units="in",height=5,width=7)
grid.arrange(p0,p1,p2,p3)
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
clone_plot = function(clone_to_plot){
  plot_clone_data = expanded %>%
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
  clone_dat = b_cells@meta.data %>% filter(donor_clone==clone_to_plot)
  write_csv(clone_dat,paste0(clone_to_plot,"_clone_dat.csv"))
}

# look at "9X14GT24_96_9_17_3325"
clone_plot("9X14GT24_96_9_17_3325")

clone_plot("TU740_22_8_17_3141")
clone_plot("FW1LNFL5_96_10_22_2890")
clone_plot("JEGK54J2_188_7_11_752")


expanded %>% filter(shm_positive=="No") %>% dplyr::count(donor_clone)
expanded %>% filter(shm_positive=="No") %>%
  dplyr::count(donor_clone,donor.id) %>%
  dplyr::count(donor.id) %>% nrow()

expanded_unmutated = expanded %>% filter(shm_positive=="No")
expanded_unmutated %>%  dplyr::count(isotype)
expanded_unmutated %>%  dplyr::count(cell_type)
expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(v_call_genotyped_VDJ)
expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(v_call_genotyped_VJ)
expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(c_call_VDJ)

expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(c_call_VJ)
expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(oligo_pos)
expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(drb_pos)
expanded_unmutated %>%  distinct(donor_clone,.keep_all=T) %>% dplyr::count(Category_fine)

# see if present in others
shm_neg_clones = b_cells@meta.data %>%
  filter(expanded_clone=="Expanded" & shm_positive=="No") %>%
  dplyr::count(donor_clone,donor.id,phenotype)

b_cells@meta.data %>%
  filter(donor_clone %in% shm_neg_clones$donor_clone) %>%
  group_by(donor_clone) %>%
  dplyr::count(shm_positive) %>%
  mutate(prop = n/sum(n))


b_cells@meta.data %>%
  filter(donor_clone %in% shm_neg_clones$donor_clone) %>%
  group_by(donor_clone) %>%
  dplyr::count(isotype) %>%
  mutate(prop = n/sum(n))

b_cells@meta.data %>%
  filter(donor_clone %in% shm_neg_clones$donor_clone) %>%
  group_by(donor_clone) %>%
  dplyr::count(cell_type) %>%
  mutate(prop = n/sum(n))

# look at in more detail
expanded_shmneg = subset(b_cells,
  cell_id %in% expanded_unmutated$cell_id)

# compare cdr3 in csf vs pbmc
res = list()
for(x in unique(expanded$donor_clone)){
res[[length(res)+1]] = expanded %>%
  filter(donor_clone == x) %>%
  group_by(source) %>%
  summarise(cdr3mu = mean(cdr3_mu_freq_overall)) %>%
  mutate(donor_clone = x)
}
res = do.call("bind_rows",res)
donors_to_keep = res %>% dplyr::count(donor_clone) %>% filter(n==2)
res %>% filter(donor_clone %in% donors_to_keep$donor_clone)

#################################
# check for public clones
################################
public_clones = b_cells@meta.data %>%
  group_by(clone_id) %>%
  dplyr::count(donor.id)  %>%
  dplyr::count(clone_id) %>%
  arrange(desc(n)) %>% filter(n>1)

b_cells@meta.data %>% filter(clone_id %in% public_clones$clone_id)

bcr_db = read_csv("bcr_db.csv")
b_cells@meta.data %>% filter(junction_aa_VDJ %in% bcr_db$CDR3.heavy.aa)
bcr_db %>% filter(CDR3.heavy.aa %in% b_cells@meta.data$junction_aa_VDJ)

# lanz
lanz_cdr3 = read_csv("../Lanz_AllCsfCdr3ForBenJacobs.csv")
b_cells@meta.data %>% filter(junction_aa_VDJ %in% lanz_cdr3$`JUNCTION | HC`)
b_cells@meta.data %>% filter(junction_aa_VJ %in% lanz_cdr3$`JUNCTION | LC`)


# hamming distance from expanded clones to lanz
expanded_cdr3 = b_cells@meta.data %>% distinct(donor_clone,.keep_all=T)
expanded_cdr3 = expanded_cdr3$junction_aa_VDJ

overlaps = list()
for(i in c(1:length(expanded_cdr3))){
  message(i)
  this_cdr3 = expanded_cdr3[i]
  # check length
  length_cdr3 = nchar(this_cdr3)

  # filter lanz to same length
  matching_length_cdr3_lanz = lanz_cdr3 %>% filter(nchar(`JUNCTION | HC`)==length_cdr3)
  matching_length_cdr3_lanz = matching_length_cdr3_lanz$`JUNCTION | HC`

  # iterate through each of these
  for(j in c(1:length(matching_length_cdr3_lanz))){
    cdr3_to_match = matching_length_cdr3_lanz[j]
    if(length(matching_length_cdr3_lanz)==0){
      next
    } else {
    # iterate along each letter to see if different
    hamming = c()
      for(k in c(1:length_cdr3)){
        res = substr(cdr3_to_match,k,k) == substr(this_cdr3,k,k)
        hamming <<- c(hamming,res)
      }
    lnh = sum(hamming) / length_cdr3
    if(lnh>0.7){
      message("LNH = ", round(lnh,2))
      overlaps[[length(overlaps)+1]] = data.frame(this_cdr3,cdr3_to_match,lnh)
    }
  }
}
}

overlaps = do.call("bind_rows",overlaps)

overlap_detail = b_cells@meta.data %>%
  filter(junction_aa_VDJ %in% overlaps$this_cdr3) %>%
  left_join(overlaps %>% distinct(this_cdr3,.keep_all=T) %>% dplyr::rename("junction_aa_VDJ"=this_cdr3),by="junction_aa_VDJ") %>%
  left_join(lanz_cdr3 %>% dplyr::rename("cdr3_to_match" = `JUNCTION | HC`),by="cdr3_to_match")

write_csv(overlap_detail,"overlaps.csv")

################################
# de clonal vs not
################################

b_cells@meta.data = b_cells@meta.data %>% mutate(donor_expanded = paste0(donor.id,"_",expanded_clone))

do_de = function(dat,plot_title){
  DefaultAssay(dat) = "RNA"
  cells_for_de = dat

  # print out params
  message("Cells:",cells_for_de@meta.data %>% nrow)
  message("Expanded cells:",cells_for_de@meta.data %>% filter(expanded_clone=="Expanded") %>% nrow)
  message("Non-expanded cells:",cells_for_de@meta.data %>% filter(expanded_clone=="Not expanded") %>% nrow)
  message("Donors expanded:",cells_for_de@meta.data %>% filter(expanded_clone=="Expanded") %>% distinct(donor.id) %>% nrow)
  message("Donors non-expanded:",cells_for_de@meta.data %>% filter(expanded_clone=="Not expanded") %>% distinct(donor.id) %>% nrow)

  rownames(cells_for_de@meta.data) = colnames(cells_for_de)

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
  dplyr::select(Age,Gender,donor_expanded) %>%
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
  scale_color_manual(values = colours)+
  scale_x_continuous(limits=c(-7,7))+
  scale_y_continuous(limits=c(0,10))

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
  nperm=100000
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
  plot = ggplot(topres,aes(NES,pathway,label=ifelse(fdr<0.05,"*"," "),fill=direction))+
    theme_bw()+
    labs(x="Normalised enrichment score",fill="Pathway enrichment")+
    geom_col(color="black")+
    geom_text(size=5,position = position_dodge(0.5))+
    scale_fill_manual(values = colours,labels = c("Up in expanded clones","Down in expanded clones"))+
    guides(alpha=FALSE)+
    theme_minimal()

  png(paste0(plot_title,"_gsea_expanded_vs_not.png"),res=300,units="in",height=3,width=6)
  print(plot)
  dev.off()
}

do_de(
  dat = subset(b_cells,phenotype=="MS" & source == "CSF" & cell_type == "Plasma cells"),
  plot_title = "MS_CSF_PCs"
  )

  do_de(
    dat = subset(b_cells,phenotype=="MS" & source == "PBMC" & cell_type == "Plasma cells"),
    plot_title = "MS_PBMC_PCs"
    )

do_de(
  dat = subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type == "Memory B cells"),
  plot_title = "MS_CSF_MemB"
  )

# heatmap
features =
c("PPIB",
"LDHA",
"SUB1",
"ARHGDIB",
"ARPC1B",
"VOPP1")

dat = subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type == "Memory B cells")
png("featureplot_clonal_markers.png",res=300,units="in",width=10,height=10)
FeaturePlot(b_cells,features=features,split.by="expanded_clone")
dev.off()

# gene score
clonal_sig = features
add_score = function(genelist,name = "gene_score"){
features = list(genelist)
# module score
DefaultAssay(b_cells) = "RNA"
b_cells = AddModuleScore(b_cells,features = features,nbin=10,ctrl=1000,name=name)
b_cells[[paste0(name,"_z")]] = RNOmni::RankNorm(b_cells@meta.data[[paste0(eval(name),"1")]])
b_cells
}

b_cells = add_score(clonal_sig)


png("clonal_sig_genes_all_cells.png",res=300,units="in",width=4,height=4)
FeaturePlot(b_cells,features = "gene_score_z")+
scale_colour_gradient(low="darkblue",high="orange",limits=c(-4,4))+
theme_minimal()+
theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
                ggtitle(" ")
dev.off()

png("clonal_sig_genes_ms_csf.png",res=300,units="in",width=4,height=4)
FeaturePlot(subset(b_cells,phenotype=="MS" & source=="CSF"),features = "gene_score_z")+
scale_colour_gradient(low="darkblue",high="orange",limits=c(-4,4))+
theme_minimal()+
theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
                ggtitle("CSF")
dev.off()

png("clonal_sig_genes_ms_pbmc.png",res=300,units="in",width=4,height=4)
FeaturePlot(subset(b_cells,phenotype=="MS" & source=="PBMC"),features = "gene_score_z")+
scale_colour_gradient(low="darkblue",high="orange",limits=c(-4,4))+
theme_minimal()+
theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
                ggtitle("PBMC")
dev.off()

cell_types = c("Naive B cells","Memory B cells","Plasma cells")

outputs = list()
for(cell in cell_types){
  for(this_source in c("CSF","PBMC")){
    a = b_cells@meta.data %>% filter(phenotype=="MS") %>% filter(source == this_source  & cell_type==cell)
    b = b_cells@meta.data %>% filter(phenotype=="OIND") %>% filter(source == this_source & cell_type==cell)
    p = ifelse(nrow(a)==0 | nrow(b)==0,NA,t.test(a$gene_score_z,b$gene_score_z)$p.value)
    df = data.frame(this_source,cell,p,comparison = "MS_OIND",mean_a = mean(a$gene_score_z),mean_b = mean(b$gene_score_z))
    a = b_cells@meta.data %>% filter(phenotype=="MS") %>% filter(source == this_source  & cell_type==cell)
    b = b_cells@meta.data %>% filter(phenotype=="NIND") %>% filter(source == this_source & cell_type==cell)
    p = ifelse(nrow(a)==0 | nrow(b)==0,NA,t.test(a$gene_score_z,b$gene_score_z)$p.value)
    df2 = data.frame(this_source,cell,p,comparison = "MS_NIND",mean_a = mean(a$gene_score_z),mean_b = mean(b$gene_score_z))
    a = b_cells@meta.data %>% filter(phenotype=="OIND") %>% filter(source == this_source  & cell_type==cell)
    b = b_cells@meta.data %>% filter(phenotype=="NIND") %>% filter(source == this_source & cell_type==cell)
    p = ifelse(nrow(a)==0 | nrow(b)==0,NA,t.test(a$gene_score_z,b$gene_score_z)$p.value)
    df3 = data.frame(this_source,cell,p,comparison = "OIND_NIND",mean_a = mean(a$gene_score_z),mean_b = mean(b$gene_score_z))

    outputs[[length(outputs)+1]] = data.frame(bind_rows(df,df2,df3))
  }
}
outputs = do.call("bind_rows",outputs)
outputs %>% mutate(fdr = p.adjust(p,method="fdr")) %>% mutate(sig = ifelse(fdr<0.1,"yes","no")) %>% mutate(up_in_ms = ifelse(mean_a>mean_b,"up in MS"," ")) %>% dplyr::select(-mean_a,-mean_b,-p)

p = ggplot(b_cells@meta.data %>% filter(cell_type %in% cell_types & phenotype %in% c("OIND","MS")),
  aes(cell_type,gene_score_z,fill=phenotype))+
  facet_wrap(~source)+
  geom_boxplot()+
  theme_minimal()+
  labs(x="Cell type",y="Clonal gene signature \nZ score")

png("clonal_gene_score.png",res=300,units="in",width=8,height=4)
p
dev.off()


FeatureScatter(b_cells,feature1 = "gene_score_z", feature2 = "MKI68",group.by="cell_type")



png("bcells_annotations_PBMC.png",res=300,units="in",width=4,height=4)
FeaturePlot(subset(pcs, source == "PBMC"),features="gene_score_z")+ggtitle("PBMC")+scale_colour_gradient(low="darkblue",high="orange",limits=c(-4,4))
dev.off()


  do_de(
    dat = subset(b_cells,phenotype=="MS" & source=="PBMC" & cell_type == "Memory B cells"),
    plot_title = "MS_PBMC_MemB"
    )




################################
# de clonal vs not in single donors
################################


overall_results_df = data.frame()
# n=1 clonal vs not
DefaultAssay(b_cells) = "RNA"
# get names for clusters
clusters = "Plasma cells"


donor_list = unique(b_cells@meta.data$donor.id)

do_de_per_donor = function(donor, cell_type = "Plasma cells", source = "CSF"){
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

sig_results = overall_results_df %>% filter(P_adj < 0.05) %>% filter(!grepl("IG",gene))
genes_to_keep =  sig_results %>% dplyr::count(gene) %>% filter(n>=3)
plot_dat = overall_results_df %>% filter(gene %in% genes_to_keep$gene)
png("indiv_de_heatmap.png",res=300,units="in",height=6,width=6)
ggplot(plot_dat,aes(donor_clone,gene,fill=direction))+
  geom_tile(color="black")+
  theme_minimal()+
  scale_fill_brewer()
dev.off()


##################################
# repeat for memory
plot = do_de_per_donor(cell_type="Memory B cells",donor="9X14GT24")

png("big_clone.png",res=300,units="in",height=6,width=6)
plot
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

