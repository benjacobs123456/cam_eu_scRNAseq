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
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/de/")

#######################################
# Updated phenotypes and hla
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
cam_pheno = cam_pheno %>% dplyr::select(shortID,Age,Gender,Poly,Lym,Prot,CSFGlu,BloodGlu,GluRatio,OligoBands,Oligoclonal,Category_fine) %>% left_join(cam_hla %>% dplyr::rename(shortID = donor_id) %>% dplyr::select(shortID,drb1_1501_dose)) %>% dplyr::rename(ID = shortID)

eu_pheno = read_csv("../../references/Info_SC_samples_with_HLA.csv")
eu_pheno_hla = read_csv("../../references/HLA_part2.csv")

eu_pheno = eu_pheno %>%
  filter(is.na(HLA_DRB1_1501)) %>%
  dplyr::select(-contains("HLA")) %>%
  left_join(eu_pheno_hla,by="ID") %>%
  dplyr::select(-5) %>%
  bind_rows(eu_pheno %>% filter(!is.na(HLA_DRB1_1501))) %>%
  dplyr::select(-7)

eu_pheno = eu_pheno %>% dplyr::select(ID,Age,Sex,HLA_DRB1_1501) %>% dplyr::rename(Gender = Sex, drb1_1501_dose = HLA_DRB1_1501)
eu_pheno = eu_pheno %>% mutate(Gender = ifelse(Gender == "female",2,1))
pheno = bind_rows(cam_pheno,eu_pheno)
write_tsv(pheno,"../../references/overall_phenotypes_with_hla.tsv")

#######################################
# Read in data
#######################################

# Read in data
all_combo = readRDS("../datasets/all_combo_phenotypes_new_cluster_ids.rds")

# Set default assay
DefaultAssay(all_combo)="RNA"

# combine w pheno data
pheno = pheno %>% dplyr::rename(donor.id = ID) %>% filter(donor.id %in% all_combo@meta.data$donor.id) %>% distinct(donor.id,.keep_all=TRUE)

all_combo@meta.data = all_combo@meta.data %>% left_join(pheno,by="donor.id")
rownames(all_combo@meta.data) = colnames(all_combo)

# stash new names
all_combo[['cell_type']] = Idents(all_combo)

# remove platelets
all_combo = subset(all_combo,subset = cell_type != "Platelets")

# merge with crude cell labels
new_cluster_ids = read_csv("../datasets/cluster_identities.csv")
new_cluster_ids = new_cluster_ids %>% dplyr::select(cell_type,cell_type_crude)

all_combo@meta.data = all_combo@meta.data %>% left_join(new_cluster_ids,by="cell_type")
rownames(all_combo@meta.data) = colnames(all_combo)

# big picture numbers
table(all_combo@meta.data$source,all_combo@meta.data$phenotype)

# quick split by pheno
ms_subtypes_csf = DimPlot(subset(all_combo,subset = source == "CSF" & Category_fine %in% c("RMS","PPMS")),split.by="Category_fine",group.by="cell_type_crude")
ms_subtypes_pbmc = DimPlot(subset(all_combo,subset = source == "PBMC" & Category_fine %in% c("RMS","PPMS")),split.by="Category_fine",group.by="cell_type_crude")

#######################################
# Extract B and T cells (& myeloid)
#######################################

# hive off B cells for VDJ analysis
b_cell_clusters =  c("B cells","Plasma cells")
b_cells = subset(all_combo, subset = cell_type_crude %in% b_cell_clusters)
saveRDS(b_cells,"../datasets/b_cells.rds")

# hive off T cells for VDJ analysis
t_cell_clusters = c("CD4 T cells","CD8 T cells","Tregs","MAIT")
t_cells = subset(all_combo, subset = cell_type_crude %in% t_cell_clusters)
saveRDS(t_cells,"../datasets/t_cells.rds")

# myeloid cells
myeloid_clusters =  c("CD14 Mono","CD16 Mono","mDCs","Macrophages","pDCs")
myeloid = subset(all_combo, subset = cell_type_crude %in% myeloid_clusters)
saveRDS(myeloid,"../datasets/myeloid.rds")


#######################################
# Diff abundance with edgeR
#######################################
# create new unique ID with donor and source
all_combo@meta.data = all_combo@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
all_combo@meta.data = all_combo@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))


# stash sample info
sample_info = all_combo@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

# filter out clusters with <10 counts
abundances = abundances[rowSums(abundances)>=10,]

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
y.ab$samples =  y.ab$samples %>%
  left_join(all_combo@meta.data %>% distinct(donor_source,.keep_all=TRUE) %>%
  dplyr::select(Age,Gender,donor_source),by="donor_source")

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
  write_csv(res,paste0("./da_plots/edgeR_da_tests_",contrast_name,".csv"))
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

  png(paste0("./da_plots/da_plot_",contrast_name,"_.png"),res=300,units="in",height=6,width=6)
  print(plot)
  dev.off()
}

counts = all_combo@meta.data %>% group_by(cell_type,donor.id,phenotype,source) %>% dplyr::count() %>% arrange(cell_type,source)
cell_medians = counts %>% group_by(phenotype,source,cell_type) %>% summarise(median = median(n),max = max(n))

png("./da_plots/absolute_cell_counts.png",res=300,width=20,height=8,units="in")
ggplot(counts,aes(cell_type,n,fill=phenotype))+
facet_wrap(~source)+
geom_boxplot(position=position_dodge(width=1))+
scale_y_log10()+
geom_text(data=cell_medians,mapping=aes(cell_type,1.2*max,label=median),position=position_dodge(width=1))+
theme_bw()+
labs(x="Cell type",y="Log10(cell count)")+
scale_fill_brewer(palette="Set2")
dev.off()

totals = counts %>% group_by(donor.id,source) %>% summarise(total = sum(n))
counts = counts %>% left_join(totals,by=c("donor.id","source"))
counts = counts %>% mutate(prop = n / total)

png("./da_plots/proportion_cell_counts.png",res=300,width=20,height=8,units="in")
ggplot(counts,aes(cell_type,prop,fill=phenotype))+facet_wrap(~source)+geom_boxplot()+theme_bw()+scale_y_log10()+labs(x="Cell type",y="Cell proportion")+scale_fill_brewer(palette="Set2")
dev.off()

png("./da_plots/proportion_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(all_combo@meta.data,aes(phenotype,fill=cell_type))+geom_bar(position="fill")+facet_wrap(~source)+scale_fill_brewer(palette="Set3")+theme_classic()+labs(x="Phenotype",y="Proportion",fill="Cell type")
dev.off()

png("./da_plots/absolute_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(all_combo@meta.data,aes(phenotype,fill=cell_type))+geom_bar()+facet_wrap(~source)+scale_fill_brewer(palette="Set3")+theme_classic()+labs(x="Phenotype",y="Absolute cell count",fill="Cell type")
dev.off()


###############################################
# Diff abundance with edgeR - using CellTypist
###############################################
# NOW SWITCH LABELS TO CELLTYPIS CELL TYPE LABELS FOR DA
all_combo = SetIdent(all_combo,value=all_combo[['ann_celltypist_highres']])
all_combo[['cell_type']] = all_combo[['ann_celltypist_highres']]


# create new unique ID with donor and source
all_combo@meta.data = all_combo@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
all_combo@meta.data = all_combo@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))


# stash sample info
sample_info = all_combo@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

# filter out clusters with <50 counts
abundances = abundances[rowSums(abundances)>50,]

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
y.ab$samples =  y.ab$samples %>% left_join(all_combo@meta.data %>% distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Gender,drb1_1501_dose,donor_source),by="donor_source")

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
  write_csv(res,paste0("./da_plots/celltypist_edgeR_da_tests_",contrast_name,".csv"))
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

  png(paste0("./da_plots/celltypist_da_plot_",contrast_name,"_.png"),res=300,units="in",height=6,width=6)
  print(plot)
  dev.off()
}



# NOW SWITCH LABELS TO CRUDE CELL TYPE LABELS FOR CRUDE DA AND DE
all_combo = SetIdent(all_combo,value=all_combo[['cell_type_crude']])
all_combo[['cell_type_fine']] = all_combo[['cell_type']]
all_combo[['cell_type']] = all_combo[['cell_type_crude']]




#######################################
# Diff abundance with edgeR
#######################################
# create new unique ID with donor and source
all_combo@meta.data = all_combo@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
all_combo@meta.data = all_combo@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

n_col = all_combo@meta.data$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

all_combo@meta.data$cell_type = factor(
all_combo@meta.data$cell_type,
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
"MAIT",
"NK cells",
"pDCs"
))

p=DimPlot(all_combo,label=F,raster=F,group.by="cell_type")+
scale_color_manual(values = colour_pal)+
ggtitle("")
png("dim_plot_simple_labels.png",res=300,units="in",width=5,height=4)
p
dev.off()


# stash sample info
sample_info = all_combo@meta.data %>%
dplyr::select(donor.id,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

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
y.ab$samples =  y.ab$samples %>% left_join(all_combo@meta.data %>% distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Gender,drb1_1501_dose,donor_source),by="donor_source")

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
  write_csv(res,paste0("./da_plots/crude_labels_edgeR_da_tests_",contrast_name,".csv"))
  res$significant = ifelse(res$P_adj<0.01,"yes","no")
  res$direction = ifelse(res$logFC>0,"Up","Down")
  comparison_label = contrast_name

  # refactor cell types (for plotting)
  res$cell = factor(res$cell,ordered=TRUE,levels = c("B cells",
"Plasma cells",
"mDCs",
"HSPCs",
"CD14 Mono",
"CD16 Mono",
"Macrophages",
"CD4 T cells",
"CD8 T cells",
"Tregs",
"MAIT",
"NK cells",
"pDCs"
))

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
  scale_x_continuous(limits=c(-4,4))+
  geom_point(shape=16,size=3)


  png(paste0("./da_plots/crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
  print(plot)
  dev.off()
}


# NB repeat for CSF vs PBMC - different scales
contrast =  makeContrasts(
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
  write_csv(res,paste0("./da_plots/crude_labels_edgeR_da_tests_",contrast_name,".csv"))
  res$significant = ifelse(res$P_adj<0.01,"yes","no")
  res$direction = ifelse(res$logFC>0,"Up","Down")
  comparison_label = contrast_name

  # refactor cell types (for plotting)
  res$cell = factor(res$cell,ordered=TRUE,levels = levels(factor(all_combo@meta.data$cell_type %>% unique,ordered=TRUE)))

  # de plot



  plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
  theme_classic()+
  scale_color_manual(values = colour_pal)+
  geom_text_repel(size=3,max.overlaps=100)+
  geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
  geom_vline(xintercept = 0,alpha=0.2)+
  NoLegend()+
  ggtitle(comparison_label)+
  scale_y_continuous(limits=c(0,40))+
  scale_x_continuous(limits=c(-6,6))+
  geom_point(shape=16,size=3)


  png(paste0("./da_plots/crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
  print(plot)
  dev.off()
}



# reorder phenotypes
all_combo@meta.data$phenotype = factor(all_combo@meta.data$phenotype,levels=c("Noninflammatory","OIND","MS"),ordered=TRUE)

counts = all_combo@meta.data %>% group_by(cell_type,donor.id,phenotype,source) %>% dplyr::count() %>% arrange(cell_type,source)
cell_medians = counts %>% group_by(phenotype,source,cell_type) %>% summarise(median = median(n),max = max(n))

png("./da_plots/crude_labels_absolute_cell_counts.png",res=300,width=20,height=8,units="in")
ggplot(counts,aes(cell_type,n,fill=phenotype))+
facet_wrap(~source)+
geom_boxplot(position=position_dodge(width=1))+
scale_y_log10()+
geom_text(data=cell_medians,mapping=aes(cell_type,1.2*max,label=median),position=position_dodge(width=1))+
theme_bw()+
labs(x="Cell type",y="Log10(cell count)")+
scale_fill_brewer(palette="Set2")
dev.off()

totals = counts %>% group_by(donor.id,source) %>% summarise(total = sum(n))
counts = counts %>% left_join(totals,by=c("donor.id","source"))
counts = counts %>% mutate(prop = n / total)

png("./da_plots/crude_labels_proportion_cell_counts.png",res=300,width=20,height=8,units="in")
ggplot(counts,aes(cell_type,prop,fill=phenotype))+facet_wrap(~source)+geom_boxplot()+theme_bw()+scale_y_log10()+labs(x="Cell type",y="Cell proportion")+scale_fill_brewer(palette="Set2")
dev.off()

png("./da_plots/crude_labels_proportion_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(all_combo@meta.data,aes(phenotype,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Phenotype",y="Proportion",fill="Cell type")
dev.off()

png("./da_plots/crude_labels_proportion_cell_counts_barplot_just_csf.png",res=300,width=4,height=4,units="in")
ggplot(all_combo@meta.data %>% filter(source=="CSF"),aes(phenotype,fill=cell_type))+
geom_bar(position="fill",color="black")+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Phenotype",y="Proportion",fill="Cell type")
dev.off()

png("./da_plots/crude_labels_absolute_cell_counts_barplot.png",res=300,width=7,height=4,units="in")
ggplot(all_combo@meta.data,aes(phenotype,fill=cell_type))+geom_bar()+facet_wrap(~source)+scale_fill_brewer(palette="Set3")+theme_classic()+labs(x="Phenotype",y="Absolute cell count",fill="Cell type")
dev.off()


##############################
###    Compare celltypist
##############################
lowres_celltypes = all_combo@meta.data %>% dplyr::count(ann_celltypist_lowres) %>% filter(n>50)
highres_celltypes = all_combo@meta.data %>% dplyr::count(ann_celltypist_highres) %>% filter(n>50)

png("celltypist_anno_lowres.png",res=300,units="in",height=6,width=8)
ggplot(all_combo@meta.data %>% filter(ann_celltypist_lowres %in% lowres_celltypes$ann_celltypist_lowres),
aes(cell_type,fill=ann_celltypist_lowres))+
geom_bar(position="fill",color="black")+
theme_bw()+
scale_fill_brewer(palette="Paired")
dev.off()

png("celltypist_anno.png",res=300,units="in",height=6,width=12)
ggplot(all_combo@meta.data %>% filter(ann_celltypist_highres %in% highres_celltypes$ann_celltypist_highres),
aes(cell_type,fill=ann_celltypist_highres))+
geom_bar(position="fill",color="black")+
theme_bw()
dev.off()

##############################
###    DE with edgeR       ###
##############################

# convert to sce object
all_combo.sce = as.SingleCellExperiment(all_combo)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = all_combo@meta.data %>%
  group_by(donor.id,source,phenotype,cell_type) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# aggregate counts
groups = colData(all_combo.sce)[, c("ident", "phenotype","source","donor.id")]
aggregated_counts  = aggregate.Matrix(t(counts(all_combo.sce)),
groupings = groups, fun = "sum")

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]


# get names for clusters
clusters = levels(factor(all_combo@meta.data$cell_type))
clusters = clusters[!clusters %in% c("HSPCs","Macrophages")]

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
    y$samples =  y$samples %>% mutate(full_cell_id = rownames(y$samples)) %>% left_join(all_combo@meta.data %>% dplyr::select(Age,Gender,drb1_1501_dose,full_cell_id) %>% distinct(full_cell_id,.keep_all=TRUE),by="full_cell_id")


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
    } else if(length((unique(y$samples$group)))==3 & !("OIND_PBMC" %in% unique(y$samples$group)) & !("Control_PBMC" %in% unique(y$samples$group))){
      makeContrasts(
        MS_CSF - MS_PBMC,
        MS_CSF - OIND_CSF,
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
    write_csv(res,paste0("./de_results/edgeR_de_tests_",contrast_name,"_",cell_type,".csv"))
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
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
    scale_color_manual(values = colours)

    png(paste0("./de_plots/de_plot_",contrast_name,"_",cell_type,".png"),res=300,units="in",height=4,width=4)
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

######################################
# individual heterogeneity
######################################


plots = list()
for(pheno in c("MS","OIND","Noninflammatory")){
p=ggplot(all_combo@meta.data %>%
filter(phenotype==pheno),
aes(donor.id,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Donor",y="Proportion",fill="Cell type")+
ggtitle(pheno)+
theme(axis.text.x=element_blank())
plots[[length(plots)+1]] = p
}

png("./crude_labels_proportion_cell_counts_barplot_per_individual.png",res=300,width=12,height=12,units="in")
do.call("grid.arrange",plots)
dev.off()

# relationship between cell type proportions and phenotypes
# define ocb positivity
all_combo@meta.data = all_combo@meta.data %>%
mutate(oligo_pos = ifelse(Oligoclonal %in% c("Pos(>10)","Pos(3to9)"),"OCB+","OCB-")) %>%
mutate(oligo_pos = ifelse(cohort=="EU" & is.na(Oligoclonal),"OCB+",oligo_pos)) %>%
mutate(drb_pos = ifelse(drb1_1501_dose > 0 ,"DRB1*15+","DRB1*15-"))

plots = list()
for(pheno in c("RMS","PPMS")){
p=ggplot(all_combo@meta.data %>%
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
p=ggplot(all_combo@meta.data %>%
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
p=ggplot(all_combo@meta.data %>%
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

png("./crude_labels_proportion_cell_counts_barplot_per_individual_by_ms_subtype.png",res=300,width=12,height=12,units="in")
do.call("grid.arrange",plots)
dev.off()

# repeat DA pairwise between OCB +/-, DRB, and PPMS vs RMS

# first clean variables
all_combo@meta.data = all_combo@meta.data %>%
  mutate(oligo_pos = ifelse(
  oligo_pos == "OCB+",
  "OCBpos",
  "OCBneg"
  ))

all_combo@meta.data = all_combo@meta.data %>%
mutate(drb_pos = ifelse(
drb_pos == "DRB1*15+",
"DRBpos",
ifelse(!is.na(drb_pos),
"DRBneg",
NA
)))

do_da = function(x, contrasts_to_test){

  # stash sample info
  sample_info = all_combo@meta.data %>%
  filter(phenotype=="MS") %>%
  dplyr::select(donor.id,source,x,donor_source) %>%
  filter(!is.na(.data[[x]])) %>%
  distinct(donor_source,.keep_all=TRUE)

  data_for_abundances = all_combo@meta.data %>%
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
  y.ab$samples =  y.ab$samples %>% left_join(all_combo@meta.data %>%
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
    res$cell = factor(res$cell,ordered=TRUE,levels = c("B cells",
  "Plasma cells",
  "mDCs",
  "HSPCs",
  "CD14 Mono",
  "CD16 Mono",
  "Macrophages",
  "CD4 T cells",
  "CD8 T cells",
  "Tregs",
  "MAIT",
  "NK cells",
  "pDCs"
  ))

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


p=DimPlot(subset(all_combo,phenotype=="MS" & source=="CSF"),split.by="donor.id",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),plot.title = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_color_manual(values = colour_pal)+
ggtitle("")
png("./indiv_dimplots_ms_csf.png",res=300,units="in",height=8,width=8)
p
dev.off()

p=DimPlot(subset(all_combo,phenotype=="MS" & source=="PBMC"),split.by="donor.id",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),plot.title = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_color_manual(values = colour_pal)+
ggtitle("")
png("./indiv_dimplots_ms_csf.png",res=300,units="in",height=8,width=8)
p
dev.off()


saveRDS(all_combo,"../datasets/all_combo_with_updated_pheno.rds")
