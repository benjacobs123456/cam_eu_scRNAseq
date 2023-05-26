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

# filter to just celltypist B cell annotations
t_cells = subset(t_cells,subset = ann_celltypist_highres %in% t_cell_names)
t_cells = SetIdent(t_cells,value="ann_celltypist_highres")


t_cell_markers = c("CD3G","IL7R","CD8A","TCF7","FOXP3","CCL5","CCR7","MKI67","CD27","TRGV9","TRDV2")


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
abundant_cells = t_cells@meta.data %>% dplyr::count(ann_celltypist_highres) %>% filter(n>200)
t_cells = subset(t_cells, ann_celltypist_highres %in% abundant_cells$ann_celltypist_highres)

# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(iid,"_",source))

# PICK UP HERE
t_cells@meta.data$cell_type = t_cells@meta.data$ann_celltypist_highres
n_col = t_cells@meta.data$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

t_cells@meta.data$cell_type = factor(
  t_cells@meta.data$cell_type,
  ordered=TRUE,
  levels = c("Naive B cells","Memory B cells","Plasma cells",
             "Transitional B cells","B cells","Germinal center B cells",
             "Large pre-B cells","Cycling B cells","Small pre-B cells","Pre-pro-B cells","Follicular B cells"))



# plots to check annotations
p1=DimPlot(t_cells)
p2=FeaturePlot(t_cells,features=b_cell_markers)
p3=DotPlot(t_cells,features=b_cell_markers)

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

png("dim_plot_simple_labels_source.png",res=600,units="in",width=5,height=3)
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
    scale_y_continuous(limits=c(0,50))+
    scale_x_continuous(limits=c(-9,9))+
    geom_point(shape=16,size=3)


  png(paste0("./crude_labels_da_plot_",contrast_name,"_.png"),res=300,units="in",height=3,width=3)
  print(plot)
  dev.off()
}






#######################################
# Basic plots
#######################################

filtered_tcells = SetIdent(object = filtered_tcells,value=filtered_tcells[['ann_monaco']])

p1=DimPlot(filtered_tcells,group.by="ann_monaco")
p2=FeaturePlot(filtered_tcells,features=t_cell_markers)

png("feature_plot.png",res=300,units="in",width=6,height=4)
p1
dev.off()

png("feature_plot2.png",res=300,units="in",width=6,height=6)
p2
dev.off()

color_pal <- RColorBrewer::brewer.pal(12, "Set2")
color_pal <- grDevices::colorRampPalette(color_pal)(12)
p1= ggplot(filtered_tcells@meta.data,aes(phenotype,fill=ann_monaco))+
geom_bar(position="fill")+
facet_wrap(~source)+
scale_fill_manual(values=color_pal)+
theme_bw()+
labs(fill="T cell subset",x="Phenotype",y="Proportion of T cell pool")

png("t_cell_proportions.png",res=300,units="in",width=8,height=5)
p1
dev.off()

annotations = filtered_tcells@meta.data %>% group_by(ann_monaco) %>% dplyr::count(ann_monaco)
p1=ggplot(annotations,aes(ann_monaco,n,fill=ann_monaco))+geom_col()+theme_bw()+labs(x="Blueprint annotation",fill="Blueprint/ENCODE annotation")
p2=DimPlot(filtered_tcells)
p3=DotPlot(filtered_tcells,features=t_cell_markers)
p4=FeaturePlot(filtered_tcells,features=t_cell_markers)

png("annotated_tcells.png",res=300,height=6,width=10,units="in")
grid.arrange(p1,p2)
dev.off()
png("annotated_tcells_dotplot.png",res=300,height=3,width=8,units="in")
p3
dev.off()

png("annotated_tcells_features.png",res=300,height=6,width=6,units="in")
p4
dev.off()

t_cells = filtered_tcells

#######################################
# Plots
#######################################

# clonal %
expanded_clones = t_cells@meta.data %>% group_by(donor.id,clone_id) %>% dplyr::count() %>% filter(n>1) %>% mutate(donor_clone = paste0(donor.id,"_",clone_id))

t_cells@meta.data = t_cells@meta.data %>% mutate(donor_clone = paste0(donor.id,"_",clone_id)) %>%
mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))

p1=ggplot(t_cells@meta.data,aes(phenotype,fill=expanded_clone))+
geom_bar(position="fill")+
facet_wrap(~source)+
scale_fill_manual(values=color_pal)+
theme_bw()+
labs(fill="T cell subset",x="Phenotype",y="Proportion of T cell pool")

png("clonal_proportions.png",res=300,height=6,width=6,units="in")
p1
dev.off()


clones_per_donor = t_cells@meta.data %>% group_by(donor.id,phenotype) %>% dplyr::count(expanded_clone) %>% filter(expanded_clone=="Expanded")
totals =  t_cells@meta.data %>% group_by(donor.id,phenotype) %>% dplyr::count(expanded_clone) %>% summarise(total = sum(n))
clones_per_donor = clones_per_donor %>% left_join(totals,by=c("donor.id","phenotype"))  %>% mutate(prop = n/total)
p1=ggplot(clones_per_donor,aes(phenotype,prop,fill=phenotype))+geom_boxplot()+theme_classic()+labs(x="Phenotype",y="Proportion of repertoire which is part of a clone")

# clonal % per source
clones_per_donor = t_cells@meta.data %>% group_by(donor.id,phenotype,source,cohort) %>% dplyr::count(expanded_clone) %>% filter(expanded_clone=="Expanded")
totals =  t_cells@meta.data %>% group_by(donor.id,phenotype,source) %>% dplyr::count(expanded_clone) %>% summarise(total = sum(n))
clones_per_donor = clones_per_donor %>% left_join(totals,by=c("donor.id","phenotype","source"))  %>% mutate(prop = n/total)
p2=ggplot(clones_per_donor,aes(phenotype,prop,fill=phenotype))+geom_boxplot()+theme_classic()+labs(x="Phenotype",y="Proportion of repertoire which is part of a clone")+facet_wrap(~source)

ggplot(clones_per_donor,aes(n,fill=phenotype))+geom_histogram()+facet_wrap(~source)+labs(x="Expanded cells per donor")

png("expanded_plots.png",res=300,height=4,width=8,units="in")
p2
dev.off()

ms = clones_per_donor[clones_per_donor$source=="CSF" & clones_per_donor$phenotype=="MS",]$prop
oind = clones_per_donor[clones_per_donor$source=="CSF" & clones_per_donor$phenotype=="OIND",]$prop
cont = clones_per_donor[clones_per_donor$source=="CSF" & clones_per_donor$phenotype=="Noninflammatory",]$prop
t.test(ms,oind)$p.value
t.test(ms,cont)$p.value
t.test(oind,cont)$p.value


ms = clones_per_donor[clones_per_donor$source=="PBMC" & clones_per_donor$phenotype=="MS",]$prop
oind = clones_per_donor[clones_per_donor$source=="PBMC" & clones_per_donor$phenotype=="OIND",]$prop
cont = clones_per_donor[clones_per_donor$source=="PBMC" & clones_per_donor$phenotype=="Noninflammatory",]$prop
t.test(ms,oind)$p.value
t.test(ms,cont)$p.value
t.test(oind,cont)$p.value


# phenotypes of clonal cells - cell types
t_subsets = t_subsets_blueprint
identities = t_cells@meta.data %>% filter(ann_blueprint %in% t_subsets & expanded_clone=="Expanded") %>% group_by(donor.id,phenotype,ann_blueprint) %>% dplyr::count(ann_blueprint)
totals = identities %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
identities_expanded = identities %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Expanded")

identities = t_cells@meta.data %>% filter(ann_blueprint %in% t_subsets & expanded_clone=="Not expanded") %>% group_by(donor.id,phenotype,ann_blueprint) %>% dplyr::count(ann_blueprint)
totals = identities %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
identities_nonexpanded = identities %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Not expanded")

identities_overall = bind_rows(identities_nonexpanded,identities_expanded)

p1=ggplot(identities_overall,aes(ann_blueprint,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=-0.25))+
  labs(x="Phenotype",y="Cell type proportions")+facet_wrap(~phenotype)+
  ggtitle("Expanded vs non-expanded repertoire: cell type proportions")

png("expanded_singler.png",res=300,height=4,width=8,units="in")
p1
dev.off()

p1=ggplot(identities_overall %>% filter(phenotype=="MS"),aes(ann_blueprint,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  labs(y="Cell type proportions",x="")+
  theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.6))

png("expanded_singler_just_ms.png",res=300,height=3,width=4,units="in")
p1
dev.off()


# phenotypes of clonal cells - TRBV types
t_cells@meta.data$trbv_family = sapply(t_cells@meta.data$v_call_VDJ,function(x){
y=str_split(x,pattern="-",n=2)[[1]][1]
return(y)
})
trbv = t_cells@meta.data %>% filter(expanded_clone=="Expanded") %>% group_by(donor.id,phenotype,trbv_family) %>% dplyr::count(trbv_family)
totals = trbv %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
trbv_expanded = trbv %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Expanded")

trbv = t_cells@meta.data %>% filter(expanded_clone=="Not expanded") %>% group_by(donor.id,phenotype,trbv_family) %>% dplyr::count(trbv_family)
totals = trbv %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
trbv_nonexpanded = trbv %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Not expanded")

trbv_overall = bind_rows(trbv_nonexpanded,trbv_expanded)

p1=ggplot(trbv_overall %>% filter(phenotype=="MS"),aes(trbv_family,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=-0.05))+
  labs(x="Phenotype",y="trbv usage proportions")+facet_wrap(~phenotype)+
  ggtitle("Expanded vs non-expanded repertoire: trbv proportions")

png("expanded_trbv.png",res=300,height=4,width=8,units="in")
p1
dev.off()

# counts
count_tbl = t_cells@meta.data %>% group_by(source,donor.id,cohort,phenotype) %>% dplyr::count()
p1=ggplot(count_tbl %>% filter(source=="CSF"),
aes(donor.id,n,fill=cohort))+
geom_col()+
theme_classic()+
ggtitle("CSF")+
scale_y_continuous(limits=c(0,1200))+
theme(axis.text.x=element_text(angle=90))

p2=ggplot(count_tbl %>% filter(source=="PBMC"),
  aes(donor.id,n,fill=cohort))+
  geom_col()+
  theme_classic()+
  ggtitle("PBMC")+
  scale_y_continuous(limits=c(0,1200))+
  theme(axis.text.x=element_text(angle=90))


# clonal relationship plot
expanded = subset(t_cells, subset = expanded_clone == "Expanded")
expanded@meta.data$clone_id %>% unique %>% length



bb_status_df = lapply(unique(expanded@meta.data$clone_id),function(x){
this_clone_data = expanded@meta.data %>% filter(clone_id == x)
this_clone_data = this_clone_data %>%
mutate(bb_status = ifelse(all(this_clone_data$source == "CSF"),"CSF only",NA)) %>%
mutate(bb_status = ifelse(all(this_clone_data$source == "PBMC"),"PBMC only",bb_status)) %>%
mutate(bb_status = ifelse("CSF" %in% this_clone_data$source & "PBMC" %in% this_clone_data$source ,"CSF & PBMC", bb_status))
this_clone_data
})
bb_status_df = do.call("bind_rows",bb_status_df) %>% distinct(clone_id,.keep_all=TRUE)


private_status_df = lapply(unique(expanded@meta.data$clone_id),function(x){
this_clone_data = expanded@meta.data %>% filter(clone_id == x)
this_clone_data = this_clone_data %>%
mutate(private_status = ifelse(all(this_clone_data$donor.id == this_clone_data$donor.id[1]),"Private","Public"))
this_clone_data
})
private_status_df = do.call("bind_rows",private_status_df) %>% distinct(clone_id,.keep_all=TRUE)

p1=ggplot(private_status_df %>% group_by(private_status) %>% dplyr::count(),aes(private_status,n,label=n,fill=private_status))+geom_col()+geom_text(vjust=-1)+labs(x="Clone private to 1 donor?",y="N clones")+NoLegend()
p2=ggplot(bb_status_df %>% group_by(bb_status) %>% dplyr::count(),aes(bb_status,n,label=n,fill=bb_status))+geom_col()+geom_text(vjust=-1)+labs(x="Clone private to CSF or PBMC?",y="N clones")+NoLegend()
png("clones_private.png",res=300,height=6,width=6,units="in")
grid.arrange(p1,p2,ncol=2)
dev.off()

# make clonal connection plot
shared = bb_status_df %>% filter(bb_status=="CSF & PBMC")
shared_cells = t_cells@meta.data %>% filter(clone_id %in% shared$clone_id)
counts = shared_cells %>% group_by(donor.id,source,ann_blueprint,clone_id) %>% dplyr::count()

clonal_plot = ggplot(counts,aes(source,ann_blueprint,size=n,group=clone_id,col=clone_id))+geom_point()+geom_line()+facet_wrap(~donor.id)+theme_bw()+NoLegend()+labs(y="Cell type")
png("clonal_plot.png",res=300,height=12,width=12,units="in")
clonal_plot
dev.off()




# plot private vs public clones
# make clonal connection plot (public clones)
shared = private_status_df %>% filter(private_status =="Public")
shared_cells = t_cells@meta.data %>% filter(clone_id %in% shared$clone_id)
counts = shared_cells %>% group_by(donor.id,source,ann_blueprint,clone_id) %>% dplyr::count()

overall_public_clone_df = data.frame()
for(i in 1:length(unique(shared_cells$donor.id))){
  message("Doing donor ",i, " of ",length(unique(shared_cells$donor.id)))
  donor = unique(shared_cells$donor.id)[i]
  this_donor_clones = shared_cells %>% filter(donor.id==donor)
  donor1_pheno = this_donor_clones$phenotype[1]
  this_donor_clones = unique(this_donor_clones$clone_id)

  for(j in this_donor_clones){
    message("doing clone ",j)
    this_donor_partners = shared_cells %>% filter(donor.id!=donor & clone_id == j)
    this_donor_partners = this_donor_partners %>% distinct(clone_id,.keep_all=TRUE)
    for(k in 1:length(this_donor_partners$donor.id)){
      donor1 = donor
      donor2 = this_donor_partners$donor.id
      donor2_pheno = this_donor_partners$phenotype[1]

      clone.id = j
      message("adding to main df")
      overall_public_clone_df <<- bind_rows(overall_public_clone_df,data.frame(clone.id,donor1,donor1_pheno,donor2, donor2_pheno))
    }
  }
}


ms_specific_clones = sapply(overall_public_clone_df$clone.id, function(clone){
  this_clone = overall_public_clone_df %>% filter(clone.id == clone)
  ms_specific = all(this_clone$donor1_pheno == "MS") & all(this_clone$donor2_pheno == "MS")
  ms_specific
})
overall_public_clone_df$ms_specific = ms_specific_clones

clone_heatmap = ggplot(overall_public_clone_df,aes(donor1,donor2,fill=ms_specific))+geom_tile()+theme(axis.text.x = element_text(angle=90))

png("clonal_heatmap.png",res=300,height=12,width=12,units="in")
clone_heatmap
dev.off()


# lookup in database
specific_clones = overall_public_clone_df %>% filter(ms_specific==TRUE)
cdr3s = t_cells@meta.data %>%
filter(clone_id %in% specific_clones$clone.id) %>%
distinct(clone_id,.keep_all=TRUE) %>%
dplyr::select(contains("junction"),clone_id,v_call_VDJ,j_call_VDJ)  %>%
dplyr::rename("clone.id" = clone_id)
specific_clones = specific_clones %>% left_join(cdr3s,by="clone.id")

tcr_db = read_csv("tcr_db.csv")
specific_clones %>% left_join(
tcr_db %>% dplyr::rename("junction_aa_VDJ" = CDR3.beta.aa, "junction_aa_VJ" = CDR3.alpha.aa),
by = c("junction_aa_VDJ")) %>% filter(!is.na(ICDname)) %>% head


# export for TREX


ms = specific_clones %>% dplyr::select(v_call_VDJ, junction_aa_VDJ, j_call_VDJ)
colnames(ms) = c("TRBV_gene","CDR3_beta","TRBJ_gene")
ms = ms %>% filter(grepl("TRBV",TRBV_gene))
ms = ms %>% distinct(CDR3_beta,.keep_all=TRUE)
write_tsv(ms,"trex.tsv")

#######################################
# DE
#######################################

DefaultAssay(t_cells)="RNA"
t_cells@meta.data$cell_type = t_cells@meta.data$ann_blueprint
t_cells@meta.data$cell_type = factor(t_cells@meta.data$cell_type)
levels(t_cells@meta.data$cell_type) = c("CD4_Naive","CD4_Tcm","CD4_Tem","CD8_Naive","CD8_Tcm","CD8_Tem","Tregs")


# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))
t_cells@meta.data = t_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",donor.id))

# convert to sce object
t_cells.sce = as.SingleCellExperiment(t_cells)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 2

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
    write_csv(res,paste0("edgeR_de_tests_",contrast_name,"_",cell_type,".csv"))
    res$significant = ifelse(res$P_adj<0.01,"yes","no")
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
# Diff abundance with edgeR
#######################################


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
y.ab = calcNormFactors(y.ab)
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

counts = t_cells@meta.data %>% group_by(cell_type,donor.id,phenotype,source) %>% dplyr::count() %>% arrange(cell_type,source)
cell_medians = counts %>% group_by(phenotype,source,cell_type,) %>% summarise(median = median(n),max = max(n))

png("absolute_cell_counts.png",res=300,width=4,height=4,units="in")
p1=ggplot(counts %>% filter(source=="CSF"),aes(cell_type,n,fill=phenotype))+
geom_boxplot(outlier.shape=NA)+
geom_jitter(alpha=0.1)+
geom_point(data=cell_medians %>% filter(source=="CSF"),mapping=aes(cell_type,median),alpha=0.5,position=position_dodge(width=1))+
theme_bw()+
labs(x="Cell type",y="Cell count per donor")+
scale_fill_brewer(palette="Set2")+
ggtitle("CSF")

p2=ggplot(counts %>% filter(source=="PBMC"),aes(cell_type,n,fill=phenotype))+
geom_boxplot(outlier.shape=NA)+
geom_jitter(alpha=0.1)+
geom_point(data=cell_medians %>% filter(source=="PBMC"),mapping=aes(cell_type,median),alpha=0.5,position=position_dodge(width=1))+
theme_bw()+
labs(x="Cell type",y="Cell count per donor")+
scale_fill_brewer(palette="Set2")+
ggtitle("PBMC")
grid.arrange(p1,p2,ncol=1)
dev.off()

totals = counts %>% group_by(donor.id,source) %>% summarise(total = sum(n))
counts = counts %>% left_join(totals,by=c("donor.id","source"))
counts = counts %>% mutate(prop = n / total)

png("proportion_cell_counts.png",res=300,width=20,height=8,units="in")
ggplot(counts,aes(cell_type,prop,fill=phenotype))+facet_wrap(~source)+geom_boxplot()+theme_bw()+scale_y_log10()+labs(x="Cell type",y="Cell proportion")+scale_fill_brewer(palette="Set2")
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
plots = list()

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

      res = res %>% arrange(padj) %>%  dplyr::select(1,2,3,5,6) %>%  mutate(cell_type = cluster) %>%  mutate(comparison = this_comparison) %>% data.frame()
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
#do_gsea("reactome")
#do_gsea("kegg")
#do_gsea("go")


res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"gsea_res_overall.csv")

res_overall$fdr = p.adjust(res_overall$pval,method="fdr")


for(x in c( "MS CSF vs MS PBMC","MS CSF vs Control CSF","Control CSF vs Control PBMC","MS PBMC vs Control PBMC" )){
  plot_data = res_overall %>% filter(comparison == x) %>% filter(grepl("HALLMARK",pathway)) %>% arrange(NES)
  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  png(paste0("summary_",x,".png"),res=300,units="in",width=12,height=16)
  p=ggplot(plot_data,aes(cell_type,pathway,fill=NES,label=ifelse(fdr<0.05,"*","")))+
    geom_tile()+
    geom_text(size=10)+
    scale_fill_viridis_c()+
    theme_classic()+
    labs(x="Cell type",y="Gene set")+
    ggtitle(x)
  print(p)
  dev.off()
}


################################
# de clonal vs not
################################

DefaultAssay(t_cells) = "RNA"


de_plots = list()
gsea_plots = list()
for(i in 1:length(clusters)){
  message("Doing clonal vs non-clonal analysis for ",clusters[i])
  cells_for_de = subset(t_cells, subset = phenotype=="MS" & source == "CSF" & cell_type == clusters[i])
  # convert to sce object
  cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

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
      }) %>% unlist %>% factor()
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
  plot = ggplot(res,aes(logFC,-log10(PValue),color=ifelse(direction=="Up","red","blue"),alpha=ifelse(significant=="yes",1,0.1),label=gene))+
  theme_bw()+
  geom_point()+
  NoLegend()+
  geom_label_repel(data=res %>% arrange(PValue) %>% head(n=50),mapping=aes())+
  ggtitle(clusters[i])


  de = res
  ranked_genes = de %>% arrange(logFC) %>% dplyr::select(gene,logFC)
  ranked_genes_vector = ranked_genes$logFC
  names(ranked_genes_vector) = ranked_genes$gene
  message("There are ",nrow(ranked_genes)," genes in this analysis")

  # do gsea
  nperm=10000
  geneset="hallmark"
  genelist = eval(parse(text = paste0(geneset,"_list")))
  res = fgsea(genelist, stats = ranked_genes_vector, minSize=10,eps=0, nPermSimple = nperm)

  res = res %>% arrange(NES)
  res$pathway = factor(res$pathway,ordered=TRUE,levels=res$pathway)
  topres = res %>% arrange(padj) %>% head(n=40) %>% arrange(desc(NES))
  plot2 = ggplot(topres,aes(NES,pathway,label=ifelse(padj<0.1,"*"," "),alpha=-log10(pval),fill=ifelse(NES>0,"Up in expanded clones","Down in expanded clones")))+
    geom_col()+
    theme_bw()+
    labs(x="Normalised enrichment score",fill=paste0("Up or down-regulated in expanded clones"))+
    geom_text()+
    ggtitle(clusters[i])

    de_plots[[i]] = plot
    gsea_plots[[i]] = plot2
  }

png("gsea_plot_clonal_vs_not.png",res=300,units="in",height=16,width=16)
do.call("gridExtra::grid.arrange",gsea_plots)
dev.off()




################################
# de clonal vs not in single donors
################################


# n=1 clonal vs not
DefaultAssay(t_cells) = "RNA"
# get names for clusters
clusters = levels(factor(t_cells@meta.data$cell_type))


donor_list = unique(t_cells@meta.data$donor.id)

do_de_per_donor = function(donor, cell_type = "CD8 Tem", source = "CSF"){
  message("Donor: ",donor)
  message("Source: ",source)
  message("Cell type:", cell_type)
  cells = t_cells@meta.data %>% filter(phenotype=="MS" & source == source & donor.id ==donor & cell_type == cell_type)
  no_expanded = cells %>% filter(expanded_clone=="Expanded")
  non_expanded = cells %>% filter(expanded_clone=="Not expanded")
  n_clones = no_expanded %>% distinct(clone_id) %>% nrow
  if(nrow(cells)==0 | nrow(no_expanded)<50 | nrow(non_expanded)==0){
    return(NA)
  }
  cells_for_de = subset(t_cells, subset = phenotype=="MS" & source == source & donor.id ==donor & cell_type == cell_type)

  # convert to sce object
  cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

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

  # de plot
  colours = c("Up" = "red","Down" = "blue", "neither" = "grey")
  plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,alpha=ifelse(significant=="yes",1,0.1),label=gene))+
  theme_bw()+
  geom_point()+
  NoLegend()+
  geom_text_repel(data=res %>% arrange(PValue) %>% head(n=10),max_overlaps=50,mapping=aes())+
  ggtitle(paste0("Donor: ",donor,"\nExpanded n: ",nrow(no_expanded),"\nN clones: ",n_clones))+scale_color_manual(values = colours)
  return(plot)
}

plots = purrr::map(donor_list,do_de_per_donor)
png("cd8_tem_clonal_vs_not_csf.png",res=300,units="in",height=16,width=16)
print(grid.arrange(grobs=plots[!is.na(plots)]))
dev.off()

do_de_per_donor = function(donor, cell_type = "CD8 Tcm", source = "CSF"){
  message("Donor: ",donor)
  message("Source: ",source)
  message("Cell type:", cell_type)
  cells = t_cells@meta.data %>% filter(phenotype=="MS" & source == source & donor.id ==donor & cell_type == cell_type)
  no_expanded = cells %>% filter(expanded_clone=="Expanded")
  non_expanded = cells %>% filter(expanded_clone=="Not expanded")
  n_clones = no_expanded %>% distinct(clone_id) %>% nrow
  if(nrow(cells)==0 | nrow(no_expanded)<50 | nrow(non_expanded)==0){
    return(NA)
  }
  cells_for_de = subset(t_cells, subset = phenotype=="MS" & source == source & donor.id ==donor & cell_type == cell_type)

  # convert to sce object
  cells_for_de.sce = as.SingleCellExperiment(cells_for_de)

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

  # de plot
  colours = c("Up" = "red","Down" = "blue", "neither" = "grey")
  plot = ggplot(res,aes(logFC,-log10(PValue),color=direction,alpha=ifelse(significant=="yes",1,0.1),label=gene))+
  theme_bw()+
  geom_point()+
  NoLegend()+
  geom_text_repel(data=res %>% arrange(PValue) %>% head(n=10),max_overlaps=50,mapping=aes())+
  ggtitle(paste0("Donor: ",donor,"\nExpanded n: ",nrow(no_expanded),"\nN clones: ",n_clones))+scale_color_manual(values = colours)
  return(plot)
}

plots = purrr::map(donor_list,do_de_per_donor)
png("cd8_tcm_clonal_vs_not_csf.png",res=300,units="in",height=16,width=16)
print(grid.arrange(grobs=plots[!is.na(plots)]))
dev.off()
