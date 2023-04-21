#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(MASS)
library(SingleCellExperiment)
library(Matrix.utils)
library(reshape2)

#######################################
# functions etc
#######################################

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

# define plotting function
plot_fx = function(path,x,width=4,height=4){
  png(path,res=600,units="in",width=width,height=height)
  print(x)
  dev.off()
}

#######################################
# Read in data
#######################################

# set WD
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/datasets")

# Read in data
all_combo = readRDS("all_combo_phenotypes_new_cluster_ids.rds")

# Set default assay
DefaultAssay(all_combo)="RNA"

# stash new names
all_combo[['cell_type']] = Idents(all_combo)

# remove platelets
all_combo = subset(all_combo,subset = cell_type != "Platelets")

# big picture numbers
message("Cell numbers")
tbl0 = all_combo@meta.data %>%
  dplyr::count(source,Category)
tbl0

message("Donor n")
tbl1 = all_combo@meta.data %>%
  distinct(iid,source,.keep_all=T) %>%
  dplyr::count(Category,source,processing_site)
tbl1
write_csv(tbl0,"cell_counts.csv")
write_csv(tbl1,"donor_counts.csv")

# get metadata for dandelion preprocessing
metadata_for_tum_dandelion = all_combo@meta.data %>%
  mutate(original_barcode = rownames(all_combo@meta.data)) %>%
  filter(processing_site == "TUM") %>%
  dplyr::select(batch_id,source,donor.id,iid,original_barcode)

write_csv(metadata_for_tum_dandelion,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/metadata_for_tum_dandelion.csv")

#######################################
# Merge overlapping patients
#######################################

sample_overlap = read_csv("/rds/user/hpcjaco1/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/Sample_Overlap.csv")

# save replicates
replicates = subset(all_combo, iid %in% sample_overlap$TUM_PatID | iid %in% sample_overlap$CAM_SAMPLE_ID)

# examine replicates
for(i in c(1:nrow(sample_overlap))){
  this_source = "CSF"
  dat_to_plot = subset(replicates, iid == sample_overlap$TUM_PatID[i] | iid == sample_overlap$CAM_SAMPLE_ID[i] & source == this_source)
  p = DimPlot(dat_to_plot,split.by="processing_site",group.by="cell_type_crude")+NoLegend() +theme_umap()
  plot_fx(paste0("subject_",i,"replicate_csf.png"),p)
}

# look at major batch effects
p = DimPlot(subset(all_combo, source=="CSF" & Category == "MS"),split.by="processing_site",group.by="cell_type_crude")+NoLegend() +theme_umap()
plot_fx("batch_effects_umap_csf_ms.png",p)


# write demographics file

pheno_simple = all_combo@meta.data %>%
  distinct(iid,source,.keep_all=T) %>%
  dplyr::select(iid,Age,Sex,Category,OCB,processing_site,source)
rownames(pheno_simple) = NULL

pheno_simple = pheno_simple %>%
tidyr::pivot_wider(values_from = source,
names_from = source,
id_cols = 1:6) %>%
mutate(CSF = ifelse(!is.na(CSF),"Yes","No")) %>%
mutate(PBMC = ifelse(!is.na(PBMC),"Yes","No")) %>%
arrange(Category,processing_site) %>%
mutate(replicate = ifelse(iid %in% sample_overlap$CAM_SAMPLE_ID | iid %in% sample_overlap$TUM_PatID,"*",""))
write_csv(pheno_simple,"demographics.csv")

all_combo@meta.data = all_combo@meta.data %>%
  left_join(
    sample_overlap %>%
      dplyr::rename("iid" = TUM_PatID),
    by="iid"
  ) %>%
  mutate(iid_prior_to_merging_duplicates = iid) %>%
  mutate(iid = ifelse(iid %in% sample_overlap$TUM_PatID,CAM_SAMPLE_ID,iid))
rownames(all_combo@meta.data) = colnames(all_combo)



message("Donor n after duplicates merged")
tbl1 = all_combo@meta.data %>%
  distinct(iid,source,.keep_all=T) %>%
  dplyr::count(Category,source,processing_site)
tbl1
write_csv(tbl1,"donor_counts_no_duplicates.csv")

#######################################
# Diff abundance with edgeR
#######################################
# create new unique ID with donor and source
all_combo@meta.data = all_combo@meta.data %>% mutate(donor_source = paste0(iid,"_",source))
all_combo@meta.data$phenotype = all_combo@meta.data$Category
all_combo@meta.data = all_combo@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

# stash sample info
sample_info = all_combo@meta.data %>%
  dplyr::select(iid,source,phenotype,donor_source) %>%
  distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

# filter out clusters with <10 counts
abundances = abundances[rowSums(abundances)>=10,]

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
y.ab$samples =  y.ab$samples %>%
  left_join(
    all_combo@meta.data %>%
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
setwd("../")
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

###############################################
# Diff abundance with edgeR - using CellTypist
###############################################
# NOW SWITCH LABELS TO CELLTYPIST CELL TYPE LABELS FOR DA
all_combo = SetIdent(all_combo,value=all_combo[['ann_celltypist_highres']])
all_combo[['cell_type']] = all_combo[['ann_celltypist_highres']]


# create new unique ID with donor and source
all_combo@meta.data = all_combo@meta.data %>% mutate(donor_source = paste0(iid,"_",source))
all_combo@meta.data = all_combo@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",iid))


# stash sample info
sample_info = all_combo@meta.data %>%
dplyr::select(iid,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

# filter out clusters with <50 counts
abundances = abundances[rowSums(abundances)>50,]

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
y.ab$samples =  y.ab$samples %>%
    left_join(all_combo@meta.data %>%
                distinct(donor_source,.keep_all=TRUE) %>%
                dplyr::select(Age,Sex,donor_source),
              by="donor_source")

# now loop
da_overall_list = list()
results_df = data.frame()
design <- model.matrix(~ 0 + grouping + Age + Sex,y.ab$sample)
colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Sex")

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
all_combo@meta.data = all_combo@meta.data %>% mutate(donor_source = paste0(iid,"_",source))
all_combo@meta.data = all_combo@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

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
"MAIT cells",
"NK cells",
"pDCs"
))

# plot
p=DimPlot(all_combo,label=F,raster=F,group.by="cell_type")+
scale_color_manual(values = colour_pal)+
theme_umap()
plot_fx("./da_plots/dim_plot_simple_labels.png",p)

# stash sample info
sample_info = all_combo@meta.data %>%
dplyr::select(iid,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

# filter out clusters with 0 counts
abundances = abundances[rowSums(abundances)>0,]

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
y.ab$samples =  y.ab$samples %>%
  left_join(all_combo@meta.data %>%
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
"MAIT cells",
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
  scale_y_continuous(limits=c(0,20))+
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
all_combo@meta.data$phenotype = factor(all_combo@meta.data$phenotype,levels=c("NIND","OIND","MS"),ordered=TRUE)

counts = all_combo@meta.data %>% group_by(cell_type,iid,phenotype,source) %>% dplyr::count() %>% arrange(cell_type,source)
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

totals = counts %>% group_by(iid,source) %>% summarise(total = sum(n))
counts = counts %>% left_join(totals,by=c("iid","source"))
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


######################################
# individual heterogeneity
######################################

# just MS - restrict to samples with CSF & PBMC cells
plot_dat = all_combo@meta.data %>%
filter(Category=="MS")
counts = plot_dat %>%
dplyr::count(iid,source) %>%
tidyr::pivot_wider(id_cols = iid,names_from=source,values_from = n) %>%
filter(CSF > 0 & PBMC > 0)
plot_dat = plot_dat %>% filter(iid %in% counts$iid)

# order by PC %
pc_pcts = plot_dat %>%
filter(source=="CSF") %>%
group_by(iid) %>%
dplyr::count(cell_type) %>%
mutate(plasma_cell_prop = n/sum(n)) %>%
filter(cell_type == "Plasma cells") %>%
dplyr::select(iid,plasma_cell_prop)

plot_dat = plot_dat %>% left_join(pc_pcts,by="iid") %>%
mutate(plasma_cell_prop = ifelse(is.na(plasma_cell_prop),0,plasma_cell_prop)) %>%
arrange(desc(plasma_cell_prop))

plot_dat$iid = factor(plot_dat$iid,levels = unique(plot_dat$iid),ordered=T)
p=ggplot(plot_dat,
aes(iid,fill=cell_type))+
geom_bar(position="fill",color="black")+
facet_wrap(~source)+
scale_fill_manual(values = colour_pal)+
theme_classic()+
labs(x="Donor",y="Proportion",fill="Cell type")+
theme(axis.text.x=element_blank())

png("./da_plots/per_individual_ms_celltypes.png",res=600,width=10,height=4,units="in")
p
dev.off()

# plot pc % vs OCB status
plot_dat = all_combo@meta.data %>%
filter(Category=="MS")
counts = plot_dat %>%
dplyr::count(iid,source) %>%
tidyr::pivot_wider(id_cols = iid,names_from=source,values_from = n) %>%
filter(CSF > 0)
plot_dat = plot_dat %>% filter(iid %in% counts$iid)

pc_pcts = plot_dat %>%
filter(source=="CSF") %>%
group_by(iid,OCB) %>%
dplyr::count(cell_type) %>%
mutate(plasma_cell_prop = n/sum(n)) %>%
filter(cell_type == "Plasma cells") %>%
dplyr::select(iid,plasma_cell_prop,OCB)

p1=ggplot(pc_pcts,aes(OCB,plasma_cell_prop))+
geom_boxplot()+
theme_minimal()+
labs(y="Plasma cell proportion",x="CSF oligoclonal bands")

# define RMS v PMS
all_combo@meta.data = all_combo@meta.data %>%
mutate(Category_fine = case_when(
  Category != "MS" ~ "Not_MS",
  Category == "MS" & is.na(Category_fine) & processing_site=="TUM" & iid == "TUM_SC_19" ~ "PMS",
  Category == "MS" & is.na(Category_fine) & processing_site=="TUM" & iid != "TUM_SC_19" ~ "RMS",
  Category == "MS" & !is.na(Category_fine) ~ Category_fine,
))
plot_dat = all_combo@meta.data %>%
filter(Category=="MS")

counts = plot_dat %>%
dplyr::count(iid,source) %>%
tidyr::pivot_wider(id_cols = iid,names_from=source,values_from = n) %>%
filter(CSF > 0)
plot_dat = plot_dat %>% filter(iid %in% counts$iid)

pc_pcts = plot_dat %>%
filter(source=="CSF") %>%
group_by(iid,Category_fine) %>%
dplyr::count(cell_type) %>%
mutate(plasma_cell_prop = n/sum(n)) %>%
filter(cell_type == "Plasma cells") %>%
dplyr::select(iid,plasma_cell_prop,Category_fine)

p2=ggplot(pc_pcts,aes(Category_fine,plasma_cell_prop))+
geom_boxplot()+
theme_minimal()+
labs(y="Plasma cell proportion",x="MS subtype")

png("./da_plots/phenotypes_het.png",res=600,width=4,height=4,units="in")
grid.arrange(p1,p2)
dev.off()

# repeat DA pairwise between OCB +/-, DRB, and PPMS vs RMS

do_da = function(x, contrasts_to_test){

  all_combo@meta.data = all_combo@meta.data %>%
    mutate(donor_source = paste0(iid,"_",source))

  # stash sample info
  sample_info = all_combo@meta.data %>%
  filter(Category=="MS") %>%
  dplyr::select(iid,source,x,donor_source) %>%
  filter(!is.na(.data[[x]])) %>%
  distinct(donor_source,.keep_all=TRUE)

  data_for_abundances = all_combo@meta.data %>%
  filter(Category=="MS") %>%
  dplyr::select(iid,source,x,donor_source,cell_type) %>%
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
  distinct(donor_source,.keep_all=TRUE) %>% dplyr::select(Age,Sex,donor_source),by="donor_source")

  # now loop
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
  "MAIT cells",
  "NK cells",
  "pDCs"
  ))

    # da plot
    plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
    theme_classic()+
    scale_color_manual(values = colour_pal)+
    geom_text_repel(size=3,max.overlaps=1)+
    geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
    geom_vline(xintercept = 0,alpha=0.2)+
    NoLegend()+
    ggtitle(comparison_label)+
    scale_y_continuous(limits=c(0,7))+
    scale_x_continuous(limits=c(-5,5))+
    geom_point(shape=16,size=3)


    png(paste0("./da_plots/pheno_comparisons_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
    print(plot)
    dev.off()
    write_csv(res,paste0("./da_plots/pheno_comparisons_da_plot_",contrast_name,"_.csv"))

  }
}

# run DA
do_da(x = "Category_fine",
contrasts_to_test = c("RMS_CSF - PPMS_CSF","RMS_PBMC - PPMS_PBMC"))
do_da(x = "OCB",
contrasts_to_test = c("TRUE_CSF - FALSE_CSF","TRUE_PBMC - FALSE_PBMC"))

# get counts
count_pheno = function(x){
dat = all_combo@meta.data %>%
filter(Category=="MS")
counts = dat %>%
  dplyr::count(iid,source) %>%
  filter(source=="CSF")
dat = dat %>% filter(iid %in% counts$iid)
dat %>%
distinct(iid,.keep_all=T) %>%
dplyr::count(.data[[x]]) %>%
mutate(prop = n/sum(n))
}
count_pheno("OCB")
count_pheno("Category_fine")

# individual dim plots
p=DimPlot(subset(all_combo,Category=="MS" & source=="CSF"),split.by="iid",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),plot.title = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_color_manual(values = colour_pal)+
ggtitle("")
png("./indiv_dimplots_ms_csf.png",res=300,units="in",height=8,width=8)
p
dev.off()

p=DimPlot(subset(all_combo,Category=="MS" & source=="PBMC"),split.by="iid",ncol=9)+
theme_minimal()+
NoLegend()+
theme(axis.text.x = element_blank(),axis.text.y=element_blank(),plot.title = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_color_manual(values = colour_pal)+
ggtitle("")
png("./indiv_dimplots_ms_csf.png",res=300,units="in",height=8,width=8)
p
dev.off()

##############################
###    DE with edgeR       ###
##############################

# convert to sce object
all_combo.sce = as.SingleCellExperiment(all_combo)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = all_combo@meta.data %>%
  group_by(iid,source,phenotype,cell_type) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

# aggregate counts
groups = colData(all_combo.sce)[, c("ident", "phenotype","source","iid")]
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
      } else if(grepl("NIND_CSF",y)){
        "Control_CSF"
      } else if(grepl("NIND_PBMC",y)){
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
    y$samples =  y$samples %>% mutate(full_cell_id = rownames(y$samples)) %>%
      left_join(all_combo@meta.data %>% dplyr::select(Age,Sex,full_cell_id) %>%
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


#######################################
# Extract B and T cells
#######################################

# hive off B cells for VDJ analysis
b_cell_clusters =  c("B cells","Plasma cells")
b_cells = subset(all_combo, subset = cell_type_crude %in% b_cell_clusters)
saveRDS(b_cells,"b_cells.rds")

# hive off T cells for VDJ analysis
t_cell_clusters = c("CD4 T cells","CD8 T cells","Tregs","MAIT cells")
t_cells = subset(all_combo, subset = cell_type_crude %in% t_cell_clusters)
saveRDS(t_cells,"t_cells.rds")

# save whole dataset
saveRDS(all_combo,"all_combo_with_updated_pheno.rds")

#######################################
# Diff DA
#######################################

# read
ms = read_csv("edgeR_da_tests_MS_CSF - MS_PBMC.csv") %>% mutate(pheno = "MS")
oind = read_csv("edgeR_da_tests_OIND_CSF - OIND_PBMC.csv") %>% mutate(pheno = "OIND")
nind = read_csv("edgeR_da_tests_Control_CSF - Control_PBMC.csv") %>% mutate(pheno = "NIND")

# combine
combo_dat = bind_rows(ms,oind,nind)

# plot
ggplot(combo_dat,aes(logFC,cell,color=pheno,size=-log10(PValue)))+
geom_point()+
theme_minimal()+
geom_vline(xintercept=0,linetype="dashed")


plot_dat = combo_dat %>%
  filter(pheno %in% c("MS","NIND")) %>%
  pivot_wider(id_cols = cell,
    names_from = pheno,
    values_from = c(logFC, PValue))

plot_dat = plot_dat %>%
mutate(discordant = ifelse(
PValue_MS < 0.05 & sign(logFC_MS) != sign(logFC_NIND) & abs(logFC_MS)>0.5,
"yes","no"
))

ggplot(plot_dat,aes(logFC_NIND,logFC_MS,label=cell,color=discordant))+
geom_point()+
theme_minimal()+
geom_vline(xintercept=0)+
geom_hline(yintercept=0)+
ggrepel::geom_text_repel(data = plot_dat %>% filter(discordant=="yes"),force=50)+
scale_color_manual(values = c("grey","red"))+
theme(legend.position="none")

# repeat with OIND
plot_dat = combo_dat %>%
  filter(pheno %in% c("OIND","NIND")) %>%
  pivot_wider(id_cols = cell,
    names_from = pheno,
    values_from = c(logFC, PValue))

plot_dat = plot_dat %>%
mutate(discordant = ifelse(
PValue_OIND < 0.05 & sign(logFC_OIND) != sign(logFC_NIND) & abs(logFC_OIND)>0.5,
"yes","no"
))

ggplot(plot_dat,aes(logFC_NIND,logFC_OIND,label=cell,color=discordant))+
geom_point()+
theme_minimal()+
geom_vline(xintercept=0)+
geom_hline(yintercept=0)+
ggrepel::geom_text_repel(data = plot_dat %>% filter(discordant=="yes"),force=50)+
scale_color_manual(values = c("grey","red"))+
theme(legend.position="none")
