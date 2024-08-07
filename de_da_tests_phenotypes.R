#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
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
# Get cell numbers per person
#######################################

# set palette and get rid of platelets
all_combo@meta.data$cell_type_crude = factor(
all_combo@meta.data$cell_type_crude,
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

n_col = all_combo@meta.data$cell_type_crude %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

# add pseudo-pseudo id
pats_codex = all_combo@meta.data %>%
  distinct(iid,Category) %>%
  arrange(Category,iid) %>%
  mutate(donor_num = row_number()) %>%
  mutate(fully_anonymous_pseudoid = paste0("sc_",donor_num,"_",Category)) %>%
  dplyr::select(iid,fully_anonymous_pseudoid)

# add to main dataset
all_combo@meta.data = all_combo@meta.data %>%
  left_join(pats_codex,by="iid")
cell_counts_per_person = all_combo@meta.data %>%
  group_by(fully_anonymous_pseudoid,source,Category,Diagnosis,OCB) %>%
  dplyr::count(cell_type_crude) %>%
  mutate(total_cells = sum(n)) %>%
  mutate(pct = 100* n / total_cells) %>%
  tidyr::pivot_wider(id_cols = c(fully_anonymous_pseudoid,Category,Diagnosis,OCB),names_from = c(source,cell_type_crude), values_from=c(pct,total_cells))

abundant_cell_type_per_person = all_combo@meta.data %>%
  group_by(fully_anonymous_pseudoid,source,Category,Diagnosis,OCB) %>%
  dplyr::count(cell_type_crude) %>%
  slice_max(n,with_ties=F) %>%
  ungroup() %>%
  pivot_wider(id_cols = 1,values_from = cell_type_crude,names_from=source)
cell_counts_per_person = cell_counts_per_person %>%
  left_join(abundant_cell_type_per_person,by="fully_anonymous_pseudoid")

write_csv(cell_counts_per_person,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/cell_counts_per_person.csv")

# save codex
write_csv(pats_codex,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/pat_codex.csv")

# replace NAs with 0s

# get average CSF cell proportions per group
cell_types_to_count = all_combo@meta.data %>% distinct(cell_type_crude,fully_anonymous_pseudoid,Category)
summary_csf_props = all_combo@meta.data %>%
 filter(!is.na(cell_type_crude) & source == "CSF") %>%
 group_by(fully_anonymous_pseudoid,Category) %>%
 dplyr::count(cell_type_crude)
cell_types_to_count = cell_types_to_count %>%
  left_join(summary_csf_props,by=c("fully_anonymous_pseudoid","Category","cell_type_crude")) %>%
  filter(!is.na(cell_type_crude)) %>%
  mutate(n = ifelse(is.na(n),0,n)) %>%
  group_by(fully_anonymous_pseudoid,Category) %>%
  mutate(total_cells = sum(n)) %>%
 filter(total_cells > 10) %>%
 mutate(pct = 100* n / total_cells) %>%
 group_by(Category,cell_type_crude) %>%
 summarise(median = median(pct),iqr_pct = IQR(pct,na.rm=T)) %>%
 mutate(med_iqr = paste0(round(median,1),"% (",round(iqr_pct,1),")")) %>%
 pivot_wider(id_cols = cell_type_crude,names_from = Category,values_from = med_iqr)

csf_counts = cell_types_to_count %>% mutate(source = "CSF")

# repeat for pbmc
# get average pbmc cell proportions per group
cell_types_to_count = all_combo@meta.data %>% distinct(cell_type_crude,fully_anonymous_pseudoid,Category)
summary_pbmc_props = all_combo@meta.data %>%
 filter(!is.na(cell_type_crude) & source == "PBMC") %>%
 group_by(fully_anonymous_pseudoid,Category) %>%
 dplyr::count(cell_type_crude)
cell_types_to_count = cell_types_to_count %>%
  left_join(summary_pbmc_props,by=c("fully_anonymous_pseudoid","Category","cell_type_crude")) %>%
  filter(!is.na(cell_type_crude)) %>%
  mutate(n = ifelse(is.na(n),0,n)) %>%
  group_by(fully_anonymous_pseudoid,Category) %>%
  mutate(total_cells = sum(n)) %>%
 filter(total_cells > 10) %>%
 mutate(pct = 100* n / total_cells) %>%
 group_by(Category,cell_type_crude) %>%
 summarise(median = median(pct),iqr_pct = IQR(pct,na.rm=T)) %>%
 mutate(med_iqr = paste0(round(median,1),"% (",round(iqr_pct,1),")")) %>%
 pivot_wider(id_cols = cell_type_crude,names_from = Category,values_from = med_iqr)

pbmc_counts = cell_types_to_count %>% mutate(source = "pbmc")

all_counts_csf_pbmc_cohort = bind_rows(csf_counts,pbmc_counts)
write_csv(all_counts_csf_pbmc_cohort,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/da_plots/pbmc_and_csf_cell_counts_summary_cohort.csv")


# simple per individual plots
cell_proportions_per_person = all_combo@meta.data %>%
  filter(!is.na(cell_type_crude)) %>%
  group_by(fully_anonymous_pseudoid,source,Category,Diagnosis,OCB) %>%
  dplyr::count(cell_type_crude)

# make NAs into 0s
cell_types_to_count = all_combo@meta.data %>% distinct(cell_type_crude,fully_anonymous_pseudoid,Category)

cell_proportions_per_person = cell_types_to_count %>%
  left_join(cell_proportions_per_person,by=c("fully_anonymous_pseudoid","Category","cell_type_crude")) %>%
  filter(!is.na(cell_type_crude)) %>%
  mutate(n = ifelse(is.na(n),0,n)) %>%
  group_by(fully_anonymous_pseudoid,source,Category,Diagnosis,OCB) %>%
  mutate(total_cells = sum(n)) %>%
  mutate(pct = 100* n / total_cells)

# plot just pbmc
pbmc_cell_props = cell_proportions_per_person %>%
  filter(source=="PBMC")
pbmc_cell_props$Category = factor(pbmc_cell_props$Category,levels = c("MS","OINDI","OIND","NIND"),ordered=T)

# order by pc %
ordered_iids = pbmc_cell_props %>%
  filter(cell_type_crude=="Plasma cells") %>%
  mutate(pc_pct = pct) %>%
  ungroup() %>%
  dplyr::select(fully_anonymous_pseudoid,pc_pct)
pbmc_cell_props = pbmc_cell_props %>% left_join(ordered_iids,by="fully_anonymous_pseudoid") %>%
mutate(pc_pct = ifelse(is.na(pc_pct),0,pc_pct)) %>%
arrange(Category,desc(pc_pct))

pbmc_cell_props$fully_anonymous_pseudoid = factor(pbmc_cell_props$fully_anonymous_pseudoid,levels = unique(pbmc_cell_props$fully_anonymous_pseudoid),ordered=T)

divides = pbmc_cell_props %>% ungroup %>%  distinct(fully_anonymous_pseudoid,Category) %>% dplyr::count(Category) %>% mutate(cumsum = cumsum(n))

p = ggplot(pbmc_cell_props,
  aes(fully_anonymous_pseudoid,pct,fill=cell_type_crude,color=Category))+
  geom_col(position=position_stack(),color="black")+
  scale_fill_manual(values = colour_pal)+
  geom_vline(xintercept=divides$cumsum[1]+0.5,linetype="dashed")+
  geom_vline(xintercept=divides$cumsum[2]+0.5,linetype="dashed")+
  geom_vline(xintercept=divides$cumsum[3]+0.5,linetype="dashed")+
  annotate("text",label="MS",x=38,y=105)+
  annotate("text",label="ID",x=79,y=105)+
  annotate("text",label="OIND",x=86,y=105)+
  annotate("text",label="NIND",x=106,y=105)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,20,40,60,80,100))+
  labs(x="Individual",y="% of PBMC cells",fill="Cell type")

png("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/da_plots/individual_level_proportion_plots_just_pbmc.png",res=900,units="in",width=10,height=6)
p
dev.off()


# plot just csf
csf_cell_props = cell_proportions_per_person %>%
  filter(source=="CSF")
csf_cell_props$Category = factor(csf_cell_props$Category,levels = c("MS","OINDI","OIND","NIND"),ordered=T)

# order by pc %
ordered_fully_anonymous_pseudoids = csf_cell_props %>%
  filter(cell_type_crude=="Plasma cells") %>%
  mutate(pc_pct = pct) %>%
  ungroup() %>%
  dplyr::select(fully_anonymous_pseudoid,pc_pct)
csf_cell_props = csf_cell_props %>% left_join(ordered_fully_anonymous_pseudoids,by="fully_anonymous_pseudoid") %>%
mutate(pc_pct = ifelse(is.na(pc_pct),0,pc_pct)) %>%
arrange(Category,desc(pc_pct))

csf_cell_props$fully_anonymous_pseudoid = factor(csf_cell_props$fully_anonymous_pseudoid,levels = unique(csf_cell_props$fully_anonymous_pseudoid),ordered=T)

divides = csf_cell_props %>% ungroup %>%  distinct(fully_anonymous_pseudoid,Category) %>% dplyr::count(Category) %>% mutate(cumsum = cumsum(n))

p = ggplot(csf_cell_props,
  aes(fully_anonymous_pseudoid,pct,fill=cell_type_crude,color=Category))+
  geom_col(position=position_stack(),color="black")+
  scale_fill_manual(values = colour_pal)+
  geom_vline(xintercept=divides$cumsum[1]+0.5,linetype="dashed")+
  geom_vline(xintercept=divides$cumsum[2]+0.5,linetype="dashed")+
  geom_vline(xintercept=divides$cumsum[3]+0.5,linetype="dashed")+
  annotate("text",label="MS",x=60,y=105)+
  annotate("text",label="ID",x=134,y=105)+
  annotate("text",label="OIND",x=155,y=105)+
  annotate("text",label="NIND",x=180,y=105)+
  theme_minimal()+
  theme(axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,20,40,60,80,100))+
  labs(x="Individual",y="% of CSF cells",fill="Cell type")

png("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/da_plots/individual_level_proportion_plots_just_csf.png",res=900,units="in",width=10,height=6)
p
dev.off()


# make demographics table
n = all_combo@meta.data %>%
  distinct(fully_anonymous_pseudoid,.keep_all=T) %>%
  dplyr::count(Category) %>%
  mutate(n = as.character(n)) %>%
  mutate("var" = "total") %>%
  pivot_wider(id_cols = var,names_from = Category,values_from=n)

sex = all_combo@meta.data %>%
  distinct(fully_anonymous_pseudoid,.keep_all=T) %>%
  group_by(Category) %>%
  dplyr::count(Sex) %>%
  mutate(pct = 100*n/sum(n)) %>%
  mutate(n_pct = paste0(n," (",round(pct,1),"%)")) %>%
  pivot_wider(id_cols = Sex,names_from = Category,values_from=n_pct)

ocb = all_combo@meta.data %>%
  distinct(fully_anonymous_pseudoid,source,.keep_all=T) %>%
  filter(!is.na(OCB) & source=="CSF") %>%
  group_by(Category) %>%
  dplyr::count(OCB) %>%
  mutate(pct = 100*n/sum(n)) %>%
  mutate(n_pct = paste0(n," (",round(pct,1),"%)")) %>%
  pivot_wider(id_cols = OCB,names_from = Category,values_from=n_pct)

age = all_combo@meta.data %>%
    distinct(fully_anonymous_pseudoid,.keep_all=T) %>%
    group_by(Category) %>%
    summarise(median_age = median(Age),iqr_age = IQR(Age)) %>%
    mutate(med_iqr = paste0(median_age," (",round(iqr_age,1),")"),"var" = "Age") %>%
    pivot_wider(id_cols = var,names_from = Category,values_from=med_iqr)

dx = all_combo@meta.data %>%
  distinct(fully_anonymous_pseudoid,.keep_all=T) %>%
  group_by(Category) %>%
  dplyr::count(Diagnosis) %>%
  mutate(pct = 100*n/sum(n)) %>%
  mutate(n_pct = paste0(n," (",round(pct,1),"%)")) %>%
  pivot_wider(id_cols = Diagnosis,names_from = Category,values_from=n_pct)

# combine
demog = bind_rows(n,sex,age,ocb,dx)
write_csv(demog,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/summary_demographics.csv")

# restore row names of seurat objects
rownames(all_combo@meta.data) = colnames(all_combo)

# remove platelets
all_combo = subset(all_combo,subset = cell_type != "Platelets")

# Look at but then exclude tysabri donor PBMC

tysabri_donor = subset(all_combo, iid == "TU146")
all_combo = subset(all_combo, iid != "TU146")



#######################################
# Diff abundance with edgeR
#######################################
setwd("../")
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
design <- model.matrix(~ 0 + grouping + Age + Sex,y.ab$samples)
colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Sex")

y.ab <- estimateDisp(y.ab, design,trend="none")
fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# define contrast for testing

contrast =  makeContrasts(
  MS_CSF - MS_PBMC,
  OIND_CSF - OIND_PBMC,
  OINDI_CSF - OINDI_PBMC,
  Control_CSF - Control_PBMC,
  MS_CSF - OIND_CSF,
  MS_CSF - OINDI_CSF,
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
  OINDI_CSF - OINDI_PBMC,
  Control_CSF - Control_PBMC,
  MS_CSF - OIND_CSF,
  MS_CSF - OINDI_CSF,
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

# set up palette
n_col = all_combo@meta.data$cell_type_crude %>% unique %>% length
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


# plot
p=DimPlot(all_combo,label=F,raster=F,group.by="cell_type",split.by="source")+
scale_color_manual(values = colour_pal)+
theme_umap()
plot_fx("./da_plots/dim_plot_simple_labels_by_source.png",p,width=6)

# stash sample info
sample_info = all_combo@meta.data %>%
dplyr::select(iid,source,phenotype,donor_source) %>%
distinct(donor_source,.keep_all=TRUE)

abundances = table(all_combo@meta.data$cell_type,all_combo@meta.data$donor_source)

# filter out clusters with 0 counts
abundances = abundances[rowSums(abundances)>10,]

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
  MS_CSF - OINDI_CSF,
  MS_CSF - Control_CSF,
  OIND_CSF - Control_CSF,
  OINDI_CSF - Control_CSF,
  MS_PBMC - OIND_PBMC,
  MS_PBMC - Control_PBMC,
  MS_PBMC - OINDI_PBMC,
  OIND_PBMC - Control_PBMC,
  OINDI_CSF - OIND_CSF,
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
  res_df = res %>% mutate(contrast = contrast_name,cell_type = cell)
  da_overall_list[[length(da_overall_list)+1]] = res_df
}

da_overall_list = do.call("bind_rows",da_overall_list)
rownames(da_overall_list) = NULL
da_overall_list = da_overall_list %>%
  mutate(P_adj = p.adjust(PValue,method="fdr")) %>%
  dplyr::select(contrast, cell, logFC, PValue, P_adj)
write_csv(da_overall_list,"./da_plots/all_da_results.csv")


# NB repeat for CSF vs PBMC - different scales
contrast =  makeContrasts(
  MS_CSF - MS_PBMC,
  OIND_CSF - OIND_PBMC,
  OINDI_CSF - OINDI_PBMC,
  Control_CSF - Control_PBMC,
  levels = design
)
da_overall_list = list()




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
  scale_y_continuous(limits=c(0,100))+
  scale_x_continuous(limits=c(-10,10))+
  geom_point(shape=16,size=3)


  png(paste0("./da_plots/crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
  print(plot)
  dev.off()
  res_df = res %>% mutate(contrast = contrast_name,cell_type = cell)
  da_overall_list[[length(da_overall_list)+1]] = res_df
}

da_overall_list = do.call("bind_rows",da_overall_list)

rownames(da_overall_list) = NULL
da_overall_list = da_overall_list %>%
  mutate(P_adj = p.adjust(PValue,method="fdr")) %>%
  dplyr::select(contrast, cell, logFC, PValue, P_adj)
write_csv(da_overall_list,"./da_plots/all_da_results_csf_v_pbmc.csv")

# plot differences
da_overall_list$cell = factor(da_overall_list$cell,ordered=TRUE,levels = c("B cells",
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

da_overall_list = da_overall_list %>%
tidyr::separate(contrast,"_",into=c("pheno","other")) %>%
  tidyr::pivot_wider(id_cols = cell,
    names_from = pheno,
    values_from = c(logFC,PValue))


make_comparison_plot = function(pheno1,pheno2){
plot_dat = da_overall_list
plot_dat = plot_dat %>% dplyr::select(1,contains(pheno1),contains(pheno2))
col_name1 = paste0("PValue_",pheno1)
col_name2 = paste0("logFC_",pheno1)
col_name3 = paste0("PValue_",pheno2)
col_name4 = paste0("logFC_",pheno2)

plot_dat = plot_dat %>%
  mutate(specific = ifelse(.data[[col_name1]] < 0.05/nrow(plot_dat) & .data[[col_name3]] > 0.05, "yes","no"))

p = ggplot(plot_dat,
aes(.data[[col_name4]],.data[[col_name2]],colour=specific, label = cell))+
geom_point()+
geom_abline(linetype="dashed",slope=1,intercept=0)+
theme_minimal()+
geom_vline(xintercept=0,alpha=0.5)+
geom_hline(yintercept=0,alpha=0.5)+
scale_color_manual(values = c("yes" = "red", "no" = "grey"))+
geom_text_repel(show.legend=F)+
ggtitle(paste0(pheno1," vs ",pheno2))

outfile = paste0("diff_da_",pheno1,"_",pheno2,".png")
png(outfile,res=600,units="in",width=6,height=6)
print(p)
dev.off()
}
make_comparison_plot("MS","Control")
make_comparison_plot("OIND","Control")
make_comparison_plot("OINDI","Control")




# reorder phenotypes
all_combo@meta.data$phenotype = factor(all_combo@meta.data$phenotype,levels=c("NIND","OIND","OINDI","MS"),ordered=TRUE)

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

#######################################
# Unpaired DA en masse
#######################################

sample_info$grouping = sample_info$source
y.ab = DGEList(abundances, samples=sample_info)
# update sample info
y.ab$samples =  y.ab$samples %>%
  left_join(all_combo@meta.data %>%
              distinct(donor_source,.keep_all=TRUE) %>%
              dplyr::select(Age,Sex,donor_source),by="donor_source")

# now loop
design <- model.matrix(~ 0 + grouping + Age + Sex,y.ab$sample)
colnames(design) = c(levels(factor(y.ab$samples$grouping)),"Age","Sex")

y.ab <- estimateDisp(y.ab, design,trend="none")
fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# define contrast for testing
contrast =  makeContrasts(
  CSF - PBMC,
  levels = design
)

# do da tests
contrast_name = "CSF vs PBMC"
res = glmQLFTest(fit, contrast = contrast)

# write to file
print(summary(decideTests(res)))
res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
res$cell = rownames(res)
write_csv(res,paste0("./da_plots/UNPAIRED_EN_MASSE_crude_labels_edgeR_da_tests_",contrast_name,".csv"))
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
ymax = max(-log10(res$PValue)) *1.1
xmax = max(abs(res$logFC)) * 1.1

  plot = ggplot(res,aes(logFC,-log10(PValue),color=cell,label=cell))+
  theme_classic()+
  scale_color_manual(values = colour_pal)+
  geom_text_repel(size=3,max.overlaps=100)+
  geom_hline(yintercept= -log10(0.05/length(res$logFC)),alpha=0.2)+
  geom_vline(xintercept = 0,alpha=0.2)+
  NoLegend()+
  ggtitle(comparison_label)+
  scale_y_continuous(limits=c(0,ymax))+
  scale_x_continuous(limits=c(-xmax,xmax))+
  geom_point(shape=16,size=3)


  png(paste0("./da_plots/UNPAIRED_EN_MASSE_crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=4)
  print(plot)
  dev.off()

######################################
# PCs
######################################

pcs = prcomp(abundances %>% t, center = T, scale = T)
pcs = pcs$x %>%
  data.frame() %>%
  mutate(donor_source = rownames(pcs$x) ) %>%
  left_join(sample_info,by="donor_source")

p = ggplot(pcs,aes(PC1,PC2,col=phenotype))+
  geom_point()+
  facet_wrap(~source)+
  theme_minimal()+
  scale_color_brewer(palette="Set1")


png("./da_plots/pc_plots_by_phenotype.png"),res=900,units="in",height=6,width=6)
print(p)
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
pc_pcts %>%
  mutate(plasma_cell_prop = plasma_cell_prop * 100) %>%
  ungroup %>%
  summarise(median(plasma_cell_prop),IQR(plasma_cell_prop),
  min(plasma_cell_prop),
  max(plasma_cell_prop)
  )
p1=ggplot(pc_pcts,aes(OCB,plasma_cell_prop))+
geom_boxplot()+
theme_minimal()+
labs(y="Plasma cell proportion",x="CSF oligoclonal bands")

# define RMS v PMS
all_combo@meta.data = all_combo@meta.data %>%
mutate(ms_subtype = case_when(
  Diagnosis == "PPMS" ~ "PPMS"
  Diagnosis == "Relapsing MS" ~ "RMS"
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
filter(!is.na(ms_subtype)) %>%
group_by(iid,ms_subtype) %>%
dplyr::count(cell_type) %>%
mutate(plasma_cell_prop = n/sum(n)) %>%
filter(cell_type == "Plasma cells") %>%
dplyr::select(iid,plasma_cell_prop,ms_subtype)

p2=ggplot(pc_pcts,aes(ms_subtype,plasma_cell_prop))+
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
do_da(x = "ms_subtype",
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
count_pheno("ms_subtype")

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
table(factor(all_combo@meta.data$cell_type))

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
      } else if(grepl("OINDI_PBMC",y)){
        "OINDI_PBMC"
      } else if(grepl("OINDI_CSF",y)){
        "OINDI_CSF"
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


# write summary of all cell counts
cell_counts_all = all_combo@meta.data %>%
dplyr::count(source,cell_type_crude)
write_csv(cell_counts_all,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/da_plots/all_cell_counts.csv")

##############################
###    UNPAIRED POOLED DE  ###
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
groups = colData(all_combo.sce)[, c("ident","source","iid")]
aggregated_counts  = aggregate.Matrix(t(counts(all_combo.sce)),
groupings = groups, fun = "sum")

# remove groups with low cell counts for DE (<n cells)
aggregated_counts = aggregated_counts[!rownames(aggregated_counts) %in% low_counts$donor_to_exclude,]

# get names for clusters
clusters = levels(factor(all_combo@meta.data$cell_type))
clusters = clusters[!clusters %in% c("HSPCs","Macrophages")]
table(factor(all_combo@meta.data$cell_type))

# split by cluster
filtered_matrices = lapply(clusters,function(x){
  filtered_matrix = aggregated_counts[grepl(x,rownames(aggregated_counts)),] %>% t()
})
names(filtered_matrices)=clusters

results_df = data.frame()
de_results = lapply(clusters,function(cell_type){

  de_input = filtered_matrices[[cell_type]]

  group_vector = lapply(colnames(de_input),function(y){
      if(grepl("_PBMC",y)){
        "PBMC"
      } else if(grepl("_CSF",y)){
        "CSF"
      }
      }) %>% unlist %>% factor()

    # make the DGE object
    y=DGEList(de_input,group=group_vector,remove.zeros=TRUE)

    # update sample info
    y$samples =  y$samples %>% mutate(full_cell_id = rownames(y$samples)) %>%
      left_join(all_combo@meta.data %>%
      mutate(full_cell_id = paste0(cell_type,"_",source,"_",iid)) %>%
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
    contrast1 = makeContrasts(CSF - PBMC,levels = design)

    # do de tests
    res = glmQLFTest(fit, contrast = contrast1)

    # write to file
    print(summary(decideTests(res)))

    res = res$table %>% mutate(P_adj = p.adjust(PValue,method="fdr")) %>% arrange(PValue)
    res$gene = rownames(res)
    write_csv(res,paste0("./de_results/POOLED_UNPAIRED_DE_edgeR_de_tests_",cell_type,".csv"))
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
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
    ggtitle(paste0(cell_type," unpaired CSF v PBMC"))+
    scale_color_manual(values = colours)

    png(paste0("./de_plots/de_plot_POOLED_UNPAIRED_",cell_type,".png"),
    res=300,units="in",height=4,width=4)
    print(plot)
    dev.off()


})



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
message("N cells:",nrow(all_combo@meta.data))
message("N people:",nrow(all_combo@meta.data %>% distinct(iid)))
saveRDS(all_combo,"all_combo_with_updated_pheno.rds")
