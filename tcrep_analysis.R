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
message(nrow(t_cells@meta.data))

# recluster
DefaultAssay(t_cells)="SCT"
t_cells = RunUMAP(t_cells,reduction="harmony",dims=1:50)
t_cells = FindNeighbors(t_cells)
t_cells = FindClusters(t_cells,resolution=0.6)

# Annotate with SingleR
monaco = celldex::MonacoImmuneData()
monaco_annotations = SingleR::SingleR(test = t_cells@assays$SCT@data,ref=monaco,labels=monaco$label.fine)
t_cells[["ann_monaco"]] = monaco_annotations$pruned.labels

t_subsets_monaco = c("Naive CD8 T cells",
"Central memory CD8 T cells",
"Effector memory CD8 T cells",
"Terminal effector CD8 T cells",
"MAIT cells",
"Vd2 gd T cells",
"Non-Vd2 gd T cells",
"Follicular helper T cells",
"T regulatory cells",
"Th1 cells",
"Th1/Th17 cells",
"Th17 cells",
"Th2 cells",
"Naive CD4 T cells",
"Terminal effector CD4 T cells")

t_subsets_blueprint = c(
"CD4+ T-cells",
"CD4+ Tcm",
"CD4+ Tem",
"CD8+ T-cells",
"CD8+ Tcm",
"CD8+ Tem",
"Tregs")


t_cell_markers = c("IL7R","CD8A","TCF7","FOXP3","CCL5")

# filter by annotation
filtered_tcells = subset(t_cells, subset = ann_monaco %in% t_subsets_monaco & ann_blueprint %in% t_subsets_blueprint)
annotations = filtered_tcells@meta.data %>% group_by(ann_monaco) %>% dplyr::count(ann_blueprint)


  p1=DimPlot(filtered_tcells,group.by="ann_monaco")
  p2=FeaturePlot(filtered_tcells,features=t_cell_markers)

  png("feature_plot.png",res=300,units="in",width=6,height=4)
  p1
  dev.off()

  png("feature_plot2.png",res=300,units="in",width=6,height=6)
  p2
  dev.off()


new_idents = filtered_tcells@meta.data$ann_monaco
names(new_idents) = colnames(filtered_tcells)
filtered_tcells = SetIdent(object = filtered_tcells,value=new_idents)

annotations = filtered_tcells@meta.data %>% group_by(ann_monaco) %>% dplyr::count(ann_blueprint)
p1=ggplot(annotations,aes(ann_monaco,n,fill=ann_blueprint))+geom_col()+theme_bw()+labs(x="Monaco annotation",fill="Blueprint/ENCODE annotation")
p2=DimPlot(filtered_tcells)
p3=DotPlot(filtered_tcells,features=t_cell_markers)
p4=FeaturePlot(filtered_tcells,features=t_cell_markers)

png("annotated_tcells.png",res=300,height=6,width=10,units="in")
grid.arrange(p1,p2)
dev.off()
png("annotated_tcells_dotplot.png",res=300,height=3,width=6,units="in")
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
t_cells@meta.data$expanded_clone = ifelse(t_cells@meta.data$clone_id_size>1,"Expanded","Not expanded")

clones_per_donor = t_cells@meta.data %>% group_by(donor.id,phenotype) %>% dplyr::count(expanded_clone) %>% filter(expanded_clone=="Expanded")
totals =  t_cells@meta.data %>% group_by(donor.id,phenotype) %>% dplyr::count(expanded_clone) %>% summarise(total = sum(n))
clones_per_donor = clones_per_donor %>% left_join(totals,by=c("donor.id","phenotype"))  %>% mutate(prop = n/total)
p1=ggplot(clones_per_donor,aes(phenotype,prop,fill=phenotype))+geom_boxplot()+theme_classic()+labs(x="Phenotype",y="Proportion of repertoire which is part of a clone")

# clonal % per source
clones_per_donor = t_cells@meta.data %>% group_by(donor.id,phenotype,source) %>% dplyr::count(expanded_clone) %>% filter(expanded_clone=="Expanded")
totals =  t_cells@meta.data %>% group_by(donor.id,phenotype,source) %>% dplyr::count(expanded_clone) %>% summarise(total = sum(n))
clones_per_donor = clones_per_donor %>% left_join(totals,by=c("donor.id","phenotype","source"))  %>% mutate(prop = n/total)
p2=ggplot(clones_per_donor,aes(phenotype,prop,fill=phenotype))+geom_boxplot()+theme_classic()+labs(x="Phenotype",y="Proportion of repertoire which is part of a clone")+facet_wrap(~source)

png("expanded_plots.png",res=300,height=4,width=14,units="in")
grid.arrange(p1,p2,ncol=2)
dev.off()

# phenotypes of clonal cells - isotypes
isotypes = t_cells@meta.data %>% filter(expanded_clone=="Expanded") %>% group_by(donor.id,phenotype,isotype) %>% dplyr::count(isotype)
totals = isotypes %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
isotypes_expanded = isotypes %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Expanded")


isotypes = t_cells@meta.data %>% filter(expanded_clone=="Not expanded") %>% group_by(donor.id,phenotype,isotype) %>% dplyr::count(isotype)
totals = isotypes %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
isotypes_nonexpanded = isotypes %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total)%>% mutate(expanded = "Not expanded")

isotypes_overall = bind_rows(isotypes_nonexpanded,isotypes_expanded)

p1=ggplot(isotypes_overall,aes(isotype,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  labs(x="Phenotype",y="Isotype proportion")+facet_wrap(~phenotype)+
  ggtitle("Expanded vs non-expanded repertoire: isotype proportions")

png("expanded_isotypes.png",res=300,height=3,width=8,units="in")
p1
dev.off()

p1=ggplot(isotypes_overall %>% filter(phenotype=="MS"),aes(isotype,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  labs(x="Phenotype",y="Isotype proportion")

png("expanded_isotypes_just_ms.png",res=300,height=3,width=4,units="in")
p1
dev.off()


# phenotypes of clonal cells - cell types
b_subsets = c("Plasma cells","Naive B cells","Memory B cells")
identities = t_cells@meta.data %>% filter(ann_monaco %in% b_subsets & expanded_clone=="Expanded") %>% group_by(donor.id,phenotype,ann_monaco) %>% dplyr::count(ann_monaco)
totals = identities %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
identities_expanded = identities %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Expanded")

identities = t_cells@meta.data %>% filter(ann_monaco %in% b_subsets & expanded_clone=="Not expanded") %>% group_by(donor.id,phenotype,ann_monaco) %>% dplyr::count(ann_monaco)
totals = identities %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
identities_nonexpanded = identities %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Not expanded")

identities_overall = bind_rows(identities_nonexpanded,identities_expanded)

p1=ggplot(identities_overall,aes(ann_monaco,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=-0.25))+
  labs(x="Phenotype",y="Cell type proportions")+facet_wrap(~phenotype)+
  ggtitle("Expanded vs non-expanded repertoire: cell type proportions")

png("expanded_singler.png",res=300,height=4,width=8,units="in")
p1
dev.off()

p1=ggplot(identities_overall %>% filter(phenotype=="MS"),aes(ann_monaco,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  labs(y="Cell type proportions",x="")+
  theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.6))

png("expanded_singler_just_ms.png",res=300,height=3,width=4,units="in")
p1
dev.off()


# phenotypes of clonal cells - TRBV types
t_cells@meta.data$ighv_family = sapply(t_cells@meta.data$v_call_genotyped_VDJ,function(x){
y=str_split(x,pattern="-",n=2)[[1]][1]
return(y)
})
ighv = t_cells@meta.data %>% filter(expanded_clone=="Expanded") %>% group_by(donor.id,phenotype,ighv_family) %>% dplyr::count(ighv_family)
totals = ighv %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
ighv_expanded = ighv %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Expanded")

ighv = t_cells@meta.data %>% filter(expanded_clone=="Not expanded") %>% group_by(donor.id,phenotype,ighv_family) %>% dplyr::count(ighv_family)
totals = ighv %>% group_by(donor.id,phenotype) %>%  summarise(total = sum(n))
ighv_nonexpanded = ighv %>% left_join(totals,by=c("donor.id","phenotype")) %>% mutate(prop = n/total) %>% mutate(expanded = "Not expanded")

ighv_overall = bind_rows(ighv_nonexpanded,ighv_expanded)

p1=ggplot(ighv_overall,aes(ighv_family,prop,fill=expanded))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=-0.05))+
  labs(x="Phenotype",y="IGHV usage proportions")+facet_wrap(~phenotype)+
  ggtitle("Expanded vs non-expanded repertoire: IGHV proportions")

png("expanded_ighv.png",res=300,height=4,width=8,units="in")
p1
dev.off()

# counts
count_tbl = t_cells@meta.data %>% group_by(source,donor.id,cohort) %>% dplyr::count()
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

tcr_db = read_csv("tcr_db.csv")

# make clonal connection plot
shared = bb_status_df %>% filter(bb_status=="CSF & PBMC")
shared_cells = t_cells@meta.data %>% filter(clone_id %in% shared$clone_id)
counts = shared_cells %>% group_by(donor.id,source,ann_monaco,clone_id) %>% dplyr::count()

clonal_plot = ggplot(counts,aes(source,ann_monaco,size=n,group=clone_id,col=clone_id))+geom_point()+geom_line()+facet_wrap(~donor.id)+theme_bw()+NoLegend()+labs(y="Cell type")
png("clonal_plot.png",res=300,height=4,width=4,units="in")
clonal_plot
dev.off()

shared_cells = shared_cells %>% mutate(cdr3_length = nchar(junction_aa_VDJ))
ggplot(shared_cells,aes(source,cdr3_mu_freq_overall,group=clone_id,col=clone_id))+
geom_point()+geom_line()+facet_wrap(~donor.id)+theme_bw()+NoLegend()+labs(y="Cell type")




# cluster profiler
# DE with B cell subsets
# correlation with genetics
# correlation with phenotype

#######################################
# Plots
#######################################

DefaultAssay(t_cells)="RNA"
t_cells@meta.data$cell_type = t_cells@meta.data$ann_monaco

# convert to sce object
t_cells.sce = as.SingleCellExperiment(t_cells)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 5

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

    # make the design matrix
    message("Doing DE for ",cell_type)
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
    res$significant = ifelse(res$P_adj<0.05,"yes","no")
    res$direction = ifelse(res$logFC>0,"Up","Down")
    total_genes = nrow(res)
    non_sig = sum(res$significant=="no")
    sig_up = sum(res$significant=="yes" & res$direction=="Up")
    sig_down = sum(res$significant=="yes" & res$direction=="Down")
    comparison_label = contrast_name
    # de plot
    plot = ggplot(res,aes(logFC,-log10(PValue),color=ifelse(direction=="Up","red","blue"),alpha=ifelse(significant=="yes",1,0.5),label=gene))+
    theme_bw()+
    geom_point()+
    NoLegend()+
    geom_label_repel(data=res %>% arrange(PValue) %>% head(n=50),mapping=aes(),max_overlaps=50)+
    ggtitle(paste0(cell_type,":", comparison_label))

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

# create new unique ID with donor and source
t_cells@meta.data = t_cells@meta.data %>% mutate(donor_source = paste0(donor.id,"_",source))

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

# now loop
da_overall_list = list()


results_df = data.frame()

design <- model.matrix(~ 0 + grouping,y.ab$sample)
colnames(design) = levels(factor((y.ab$sample$grouping)))
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
  plot = ggplot(res,aes(logFC,-log10(PValue),color=ifelse(direction=="Up","red","blue"),label=cell))+
  theme_bw()+
  geom_point()+
  NoLegend()+
  geom_label_repel(data=res %>% arrange(PValue) %>% head(n=30),mapping=aes())+
  ggtitle(comparison_label)+
  geom_hline(yintercept = -log10(0.05/3))

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

clusters = t_cells$ann_monaco %>% unique

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
do_gsea("kegg")
do_gsea("go")


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)

# make nice kegg plots
file = paste0("edgeR_de_tests_MS_CSF - MS_PBMC_Plasma cells.csv")
de = read_csv(file)

df = de
entrez = bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df = df %>% left_join(entrez %>% rename("gene" = SYMBOL),by="gene")
df = df %>% filter(!is.na(ENTREZID)) %>% dplyr::select(-gene)
gene_vector = df$logFC
names(gene_vector) = df$ENTREZID
# Overlay the expression data onto this pathway
pathview(gene.data=gene_vector, species="hsa", pathway.id="hsa04064")
pathview(gene.data=gene_vector, species="hsa", pathway.id="hsa04668")
pathview(gene.data=gene_vector, species="hsa", pathway.id="hsa00759")
pathview(gene.data=gene_vector, species="hsa", pathway.id="hsa04062")
pathview(gene.data=gene_vector, species="hsa", pathway.id="hsa04150")

# heatmap
nfkb_genes = hallmark_list$HALLMARK_TNFA_SIGNALING_VIA_NFKB
de_genes = de %>% arrange(logFC) %>% filter(!grepl("^MT-",gene)) %>% filter(!grepl("^RP",gene)) %>% filter(gene %in% nfkb_genes) %>% head(30)
leading_edge = read_csv("sig_pathways_hallmark_Plasma cells_MS_CSF - MS_PBMC.csv")
leading_edge_genes = leading_edge %>% filter(pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
leading_edge_genes = c(str_split(leading_edge_genes$leadingEdge,pattern=",")[[1]]) %>% str_replace(" ","")

png("heatmap_memB_tnf.png",res=300,units="in",height=8,width=8)
DoHeatmap(subset(t_cells, subset = phenotype == "MS" & ann_monaco == "Plasma cells"),features=leading_edge_genes,slot="data",group.by="source")
dev.off()

res_overall = do.call("rbind",res_overall)
write_csv(res_overall,"../gsea/res_overall.csv")

res_overall$fdr = p.adjust(res_overall$pval,method="fdr")


for(x in c( "MS CSF vs MS PBMC","MS CSF vs Control CSF","Control CSF vs Control PBMC","MS PBMC vs Control PBMC","MS PBMC vs OIND PBMC","MS CSF vs OIND CSF" )){
  plot_data = res_overall %>% filter(comparison == x) %>% filter(grepl("HALLMARK",pathway)) %>% arrange(fdr)
  plot_data$pathway = factor(plot_data$pathway,levels=unique(plot_data$pathway))
  png(paste0("summary_",x,".png"),res=300,units="in",width=6,height=8)
  p=ggplot(plot_data,aes(cell_type,pathway,fill=NES,label=ifelse(fdr<0.1,"*","")))+
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

# check for public clones

tcr_db = read_csv("tcr_db.csv")

# csf vs pbmc

p1=ggplot(t_cells@meta.data,
aes(cell_type,fill=isotype))+
geom_bar()+
ggtitle("Isotypes")+
facet_wrap(~source)+
theme_classic()+
labs(x="Cell type",y="N Isotype")

p2=ggplot(t_cells@meta.data,
aes(cell_type,fill=shm_positive))+
geom_bar()+
ggtitle("SHM")+
facet_wrap(~source)+
theme_classic()+
labs(x="Cell type",y="N SHM+")

p3=ggplot(t_cells@meta.data,
aes(cell_type,nchar(junction_aa_VDJ),fill=source))+
geom_boxplot()+
ggtitle("CDR3 length")+
theme_classic()+
labs(x="Cell type",y="CDR3 length")

p4=ggplot(t_cells@meta.data,
aes(ann_monaco,fill=source))+
geom_bar()+
ggtitle("Cell type (finer annotations)")+
theme_classic()+
facet_wrap(~source)+
labs(x="Cell type",y="N")+
theme(axis.text.x=element_text(angle=90))

png("csf_v_pbmc.png",res=300,height=8,width=14,units="in")
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()

################################
# de clonal vs not
################################

DefaultAssay(t_cells) = "RNA"
cells_for_de = subset(t_cells, subset = phenotype=="MS" & source == "CSF" & ann_monaco == "Plasma cells")
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
geom_label_repel(data=res %>% arrange(PValue) %>% head(n=50),mapping=aes())

png("expanded_vs_not_csf_pcs_ms.png",res=300,units="in",height=8,width=8)
print(plot)
dev.off()

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
plot = ggplot(topres,aes(NES,pathway,label=ifelse(padj<0.1,"*"," "),alpha=-log10(pval),fill=ifelse(NES>0,"Up in expanded clones","Down in expanded clones")))+
  geom_col()+
  theme_bw()+
  labs(x="Normalised enrichment score",fill=paste0("Up or down-regulated in expanded clones"))+
  geom_text()

png("gsea_results_hallmark_expanded_pcs.png",res=300,units="in",height=8,width=16)
print(plot)
dev.off()

