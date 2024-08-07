#######################################
# Load packages
#######################################

library(tidyverse)
library(Seurat)
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
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/BCR")

# Read in data
b_cells = readRDS("b_cells_post_processing.rds")
rownames(b_cells@meta.data) = colnames(b_cells)

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
nrow(b_cells@meta.data)

#######################################
# Clonal definition sens check
#######################################

sensitivity_analysis = b_cells@meta.data %>%
  mutate(cdr3_length_vdj = nchar(junction_aa_VDJ)) %>%
  mutate(cdr3_length_vj = nchar(junction_aa_VJ))


sensitivity_analysis %>%
  group_by(iid,v_call_genotyped_VJ,v_call_genotyped_VDJ) %>%
  dplyr::count() %>%
  filter(n>1)

sensitivity_analysis %>%
  group_by(iid,v_call_genotyped_VJ,v_call_genotyped_VDJ) %>%
  dplyr::count() %>%
  filter(n>1) %>%
  ungroup %>%
  distinct(iid)

clonal_groups = sensitivity_analysis %>%
  group_by(iid,v_call_genotyped_VJ,v_call_genotyped_VDJ, cdr3_length_vdj,cdr3_length_vj) %>%
  dplyr::count() %>%
  filter(n>1)


# calculate lnh
lnh_fx = function(a,b){

  running_total = list()
  str_length = nchar(as.character(a))
  for(n in c(1:str_length)){
    running_total[[n]] = substr(a,n,n) == substr(b,n,n)
  }
  lnh = sum(unlist(running_total)) / str_length
  lnh
}

all_cells = list()
for(i in c(1:nrow(clonal_groups))){

matches = sensitivity_analysis %>%
filter(
iid == clonal_groups$iid[i] &
v_call_genotyped_VJ == clonal_groups$v_call_genotyped_VJ[i] &
v_call_genotyped_VDJ == clonal_groups$v_call_genotyped_VDJ[i] &
cdr3_length_vj == clonal_groups$cdr3_length_vj[i] &
cdr3_length_vdj == clonal_groups$cdr3_length_vdj[i]
)

# make comparison df
grid = expand.grid(unique(matches$junction_aa_VDJ),unique(matches$junction_aa_VDJ))
LNH = list()
for(j in c(1:nrow(grid))){
  LNH[[j]] = lnh_fx(grid$Var1[j],grid$Var2[j])
}
grid$LNH = unlist(LNH)
grid$clonotype_id = paste0("cloneid_",i)
all_cells[[i]] = grid
}

all_cells = do.call("bind_rows",all_cells)


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
b_cells = subset(b_cells,subset = ann_celltypist_highres %in% c("B cells","Cycling B cells","Germinal center B cells","Large pre-B cells","Memory B cells","Naive B cells","Small pre-B cells","Plasma cells","Transitional B cells","Follicular B cells","Pre-pro-B cells"))
b_cells = SetIdent(b_cells,value="ann_celltypist_highres")

cell_umap_embeddings1 = b_cells@reductions$umap@cell.embeddings %>%
data.frame() %>%
filter(UMAP_1 < -5)
cell_umap_embeddings2 = b_cells@reductions$umap@cell.embeddings %>%
data.frame() %>%
filter(UMAP_1 > -5)

# filter
b_cells = subset(b_cells,
 ( ann_celltypist_highres == "Plasma cells" & cell_id %in% rownames(cell_umap_embeddings1) ) |
 ( ann_celltypist_highres != "Plasma cells" & cell_id %in% rownames(cell_umap_embeddings2) )
   )

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
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_source = paste0(iid,"_",source))

b_cells@meta.data$cell_type = b_cells@meta.data$ann_celltypist_highres
n_col = b_cells@meta.data$cell_type %>% unique %>% length
colour_pal <- RColorBrewer::brewer.pal(n_col, "Paired")
colour_pal <- grDevices::colorRampPalette(colour_pal)(n_col)

b_cells@meta.data$cell_type = factor(
  b_cells@meta.data$cell_type,
  ordered=TRUE,
  levels = c("Naive B cells","Memory B cells","Plasma cells",
             "Transitional B cells","B cells","Germinal center B cells",
             "Large pre-B cells","Cycling B cells","Small pre-B cells","Pre-pro-B cells","Follicular B cells"))



# plots to check annotations
p1=DimPlot(b_cells)
p2=FeaturePlot(b_cells,features=b_cell_markers)
p3=DotPlot(b_cells,features=b_cell_markers)

png("dimplot.png",res=600,units="in",width=6,height=6)
p1 + theme_umap() +  scale_color_manual(values = colour_pal)
dev.off()

png("featureplot.png",res=600,units="in",width=6,height=6)
p2
dev.off()

png("dotplot.png",res=300,units="in",width=10,height=12)
p3
dev.off()

p=DimPlot(b_cells,label=F,raster=F,group.by="cell_type")+
  scale_color_manual(values = colour_pal)+
  ggtitle("")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="",y="")

png("dim_plot_simple_labels.png",res=600,units="in",width=5,height=4)
p
dev.off()

p=DimPlot(b_cells,split.by="source",label=F,raster=F,group.by="cell_type")+
  scale_color_manual(values = colour_pal)+
  ggtitle("")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="",y="")

png("dim_plot_simple_labels_source.png",res=600,units="in",width=5,height=3)
p
dev.off()

# stash sample info
sample_info = b_cells@meta.data %>%
  dplyr::select(iid,source,phenotype,donor_source) %>%
  distinct(donor_source,.keep_all=TRUE)

abundances = table(b_cells@meta.data$cell_type,b_cells@meta.data$donor_source)

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
y.ab$samples =  y.ab$samples %>% left_join(b_cells@meta.data %>%
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
# Phenotypes & DA
#######################################

# reorder phenotypes
b_cells@meta.data$phenotype = factor(b_cells@meta.data$phenotype,levels=c("NIND","OIND","OINDI","MS"),ordered=TRUE)

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
for(pheno in c("MS","OINDI","OIND","NIND")){
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

#######################################
# DE
#######################################

DefaultAssay(b_cells)="RNA"

# create new unique ID with donor and source
b_cells@meta.data = b_cells@meta.data %>% mutate(donor_source = paste0(iid,"_",source))
b_cells@meta.data = b_cells@meta.data %>% mutate(full_cell_id = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

# convert to sce object
b_cells.sce = as.SingleCellExperiment(b_cells)

# tabulate to find out which 'groups' will have insufficient cells for DE
min_cells_per_sample = 10

low_counts = b_cells@meta.data %>%
  group_by(iid,source,phenotype,cell_type) %>%
  dplyr::count() %>%
  arrange(n) %>%
  filter(n<min_cells_per_sample) %>%
  mutate(donor_to_exclude = paste0(cell_type,"_",phenotype,"_",source,"_",iid))

# aggregate counts
groups = colData(b_cells.sce)[, c("ident", "phenotype","source","iid")]
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
    left_join(b_cells@meta.data %>%
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
  cell = str_split(y,"_")[[1]][4]
  y = str_remove(y,paste0("_",cell))
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
# CCL22
#######################################
DefaultAssay(b_cells) = "SCT"
FeaturePlot(b_cells,features=c("CCL22"))

ccl22 = subset(b_cells,CCL22 > 0)


b_cells@meta.data = b_cells@meta.data %>%
mutate(ccl22 = ifelse(
  cell_id %in% ccl22@meta.data$cell_id,
  "yes",
  "no"
  ))
b_cells = SetIdent(b_cells,value="ccl22")
rownames(b_cells@meta.data) = colnames(b_cells)
markers = FindMarkers(b_cells,ident.1="yes",recorrect_umi=F,logfc.threshold=1,min.pct=0.5)



#######################################
# Repertoire analysis CSF v periphery
#######################################
rownames(b_cells@meta.data) = colnames(b_cells)

# define clones
expanded_clones = b_cells@meta.data %>%
  group_by(iid,clone_id) %>%
  dplyr::count() %>%
  filter(n>1) %>%
  mutate(donor_clone = paste0(iid,"_",clone_id))

b_cells@meta.data = b_cells@meta.data %>%
  mutate(donor_clone = paste0(iid,"_",clone_id)) %>%
  mutate(expanded_clone = ifelse(donor_clone %in% expanded_clones$donor_clone,"Expanded","Not expanded"))


# define IGHV types
b_cells@meta.data$ighv_family = sapply(b_cells@meta.data$v_call_genotyped_VDJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  y = str_remove_all(y,"D")
  return(y)
})

b_cells@meta.data$iglightchain_family = sapply(b_cells@meta.data$v_call_genotyped_VJ,function(x){
  y=str_split(x,pattern="-",n=2)[[1]][1]
  y = str_remove_all(y,"D")
  return(y)
})


# cdr3 length
b_cells@meta.data = b_cells@meta.data %>% mutate(cdr3_length = nchar(junction_aa_VDJ))

# simplify v call
all_calls = list()
for(i in c(1:length(b_cells@meta.data$v_call_genotyped_VDJ))){
message(i)
this_call = b_cells@meta.data$v_call_genotyped_VDJ[i]

# take first call if ambiguous
this_call = str_split(this_call,"\\|",2)[[1]][1]

# split to simple call
this_call = ifelse(str_count(this_call,"-")>1,
  paste0(str_split(this_call,"-")[[1]][1],"-",str_split(this_call,"-")[[1]][2]),
  this_call
)
this_call = str_remove_all(this_call,"D")
all_calls[[i]] = this_call
}

# add to metadata
b_cells@meta.data$v_call_simple = unlist(all_calls)

# repeat for light chain
# simplify v call
all_calls = list()
for(i in c(1:length(b_cells@meta.data$v_call_genotyped_VJ))){
message(i)
this_call = b_cells@meta.data$v_call_genotyped_VJ[i]

# take first call if ambiguous
this_call = str_split(this_call,"\\|",2)[[1]][1]

# split to simple call
this_call = ifelse(str_count(this_call,"-")>1,
  paste0(str_split(this_call,"-")[[1]][1],"-",str_split(this_call,"-")[[1]][2]),
  this_call
)
this_call = str_remove_all(this_call,"D")
all_calls[[i]] = this_call
}

# add to metadata
b_cells@meta.data$lightchain_call_simple = unlist(all_calls)

# clean isotype data
# simplify v call
all_calls = list()
for(i in c(1:length(b_cells@meta.data$c_call_VDJ))){
message(i)
this_call = b_cells@meta.data$c_call_VDJ[i]

# take first call if ambiguous
this_call = str_split(this_call,"\\,",2)[[1]][1]

# split to simple call
this_call = ifelse(str_count(this_call,"-")>1,
  paste0(str_split(this_call,"-")[[1]][1],"-",str_split(this_call,"-")[[1]][2]),
  this_call
)
all_calls[[i]] = this_call
}

# add to metadata
b_cells@meta.data$clean_isotype = unlist(all_calls)

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
DimPlot(b_cells,split.by="source")+scale_color_manual(values = colour_pal)+
  theme_umap()
dev.off()

## isotypes
png("csf_v_pbmc_isotypes.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source",group.by="clean_isotype")+scale_color_manual(values = colour_pal)+
theme_umap()
  dev.off()

## ighv
png("csf_v_pbmc_ighv.png",res=300,height=3,width=6,units="in")
DimPlot(b_cells,split.by="source",group.by="ighv_family")+scale_color_manual(values = colour_pal)+
  theme_umap()

  dev.off()


# compare csf vs periphery - just ms
csf_v_pbmc_plot = function(x){
  DimPlot(subset(b_cells,phenotype=="MS"),split.by="source",group.by=x)+
    scale_color_manual(values = colour_pal)+
    theme_umap()
}
var_list = list("cell_type","expanded_clone","shm_positive","isotype","ighv_family","status_summary")
plot_list = lapply(var_list,csf_v_pbmc_plot)

png("csf_v_pbmc_just_ms_bcr_repertoire.png",res=300,units="in",height=8,width=8)
do.call("grid.arrange",plot_list)
dev.off()

# repeat for all phenos
# compare csf vs periphery - just ms
csf_v_pbmc_plot_all = function(x){
  DimPlot(b_cells,split.by="source",group.by=x)+
    scale_color_manual(values = colour_pal)+
    theme_umap()
}
plot_list = lapply(var_list,csf_v_pbmc_plot_all)

png("csf_v_pbmc_all_phenos_bcr_repertoire.png",res=300,units="in",height=8,width=8)
do.call("grid.arrange",plot_list)
dev.off()



# bar plots
make_categorical_plot = function(x){
  ggplot(b_cells@meta.data,aes(phenotype,fill=b_cells@meta.data[[x]]))+
    geom_bar(position="fill",color="black")+
    facet_wrap(~source)+
    theme_bw()+
    labs(x="Phenotype",fill=x,y="Proportion")+
    scale_fill_brewer(palette="Paired")
}

make_continuous_plot = function(x){
  ggplot(b_cells@meta.data,aes(phenotype,fill=cell_type,y=b_cells@meta.data[[x]]))+
    geom_boxplot(color="black")+
    facet_wrap(~source)+
    theme_bw()+
    labs(x="Phenotype",y=x)+
    scale_fill_brewer(palette="Paired")
}

var_list = list("expanded_clone","shm_positive","isotype","ighv_family","status_summary","cell_type")
plot_list = lapply(var_list,make_categorical_plot)

png("csf_v_pbmc_bcr_repertoire.png",res=300,units="in",height=14,width=14)
do.call("grid.arrange",plot_list)
dev.off()

#######################################
#  IGHV usage
#######################################
cells = c("Naive B cells","Memory B cells","Plasma cells")

do_da_var = function(x, contrasts_to_test,variable,ylim=10, data = b_cells@meta.data,plot_title,var_label){

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
        P_adj > 0.1 ~ "nonsig",
        P_adj <= 0.1 & logFC>0 ~ "Up",
        P_adj <= 0.1 & logFC<0 ~ "Down"
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
  mutate(logFC = ifelse(P_adj > 0.1, NA, logFC)),
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

fileout_csv = paste0("da_summary_plot_",plot_title,"_",variable,".csv")
write_csv(da_overall_res,fileout_csv)
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
          variable="ighv_family",
          var_label = "IGHV family",
          plot_title = "all",
          ylim=20)

# repeat, sampling each clone only once
do_da_var(data = b_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
          ,
          x = "phenotype",
          contrasts_to_test = default_contrasts,
          variable="ighv_family",
          var_label = "IGHV family",
          plot_title = "eachclone",
          ylim=20)



# individual gene segments
do_da_var(data = b_cells@meta.data,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "OIND_CSF - OIND_PBMC",
                                "OINDI_CSF - OINDI_PBMC",
                                "NIND_CSF - NIND_PBMC",
                                "MS_CSF - OIND_CSF",
                                "MS_CSF - OINDI_CSF",
                                "MS_CSF - NIND_CSF"),
          variable="v_call_simple",
          var_label = "IGHV allele",
          plot_title = "ighv_segments",
          ylim=40)


# by cell type
do_da_var(data = b_cells@meta.data %>% filter(cell_type=="Plasma cells"),
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                "MS_CSF - OINDI_CSF"),
          variable="v_call_simple",
          var_label = "IGHV allele",
          plot_title = "ighv_segments_pcs",
          ylim=20)

# by cell type
do_da_var(data = b_cells@meta.data %>% filter(cell_type=="Memory B cells"),
x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
            "MS_CSF - OIND_CSF",
            "MS_CSF - OINDI_CSF"),
variable="v_call_simple",
var_label = "IGHV allele",
plot_title = "ighv_segments_memb",
ylim=20)
# by cell type
do_da_var(data = b_cells@meta.data %>% filter(cell_type=="Naive B cells"),
x = "phenotype",
contrasts_to_test = c("MS_CSF - MS_PBMC",
                      "MS_CSF - OIND_CSF",
                      "MS_CSF - OINDI_CSF"),
variable="v_call_simple",
var_label = "IGHV allele",
plot_title = "ighv_segments_naive",
ylim=20)

# individual gene segments per clone
do_da_var(data = b_cells@meta.data %>%
            mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>%
            distinct(donor_clone_source,.keep_all=TRUE),
          x = "phenotype",
          contrasts_to_test =      c("MS_CSF - MS_PBMC",
                                      "OIND_CSF - OIND_PBMC",
                                      "OINDI_CSF - OINDI_PBMC",
                                      "NIND_CSF - NIND_PBMC",
                                      "MS_CSF - OIND_CSF",
                                      "MS_CSF - OINDI_CSF",
                                      "MS_CSF - NIND_CSF"),
          variable="v_call_simple",
          var_label = "IGHV allele",
          plot_title = "eachclone",
          ylim=20)




# light chain family
do_da_var(data = b_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
          ,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "OIND_CSF - OIND_PBMC",
                                "OINDI_CSF - OINDI_PBMC",
                                "NIND_CSF - NIND_PBMC",
                                "MS_CSF - OIND_CSF",
                                "MS_CSF - OINDI_CSF",
                                "MS_CSF - NIND_CSF"),
          variable="iglightchain_family",
          var_label = "Ig light chain",
          plot_title = "lightchain_family_eachclone",
          ylim=10)


# make summary table for paper
ighv_family_res = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/BCR/da_summary_plot_eachclone_ighv_family.csv") %>%
  mutate(fdr = p.adjust(PValue,method="fdr")) %>%
  mutate(direction = case_when(
    logFC > 0 & fdr <= 0.1 ~ "Up",
    logFC < 0 & fdr <= 0.1 ~ "Down",
    fdr > 0.1 ~ "NS")) %>%
    mutate(p_simple = ifelse(PValue < 0.001,"<0.001",round(PValue,3))) %>%
    pivot_wider(id_cols = cell,names_from = contrast,values_from = c(direction,p_simple,logFC))

iglkv_family_res = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/BCR/da_summary_plot_lightchain_family_eachclone_iglightchain_family.csv") %>%
  mutate(fdr = p.adjust(PValue,method="fdr")) %>%
  mutate(direction = case_when(
    logFC > 0 & fdr <= 0.1 ~ "Up",
    logFC < 0 & fdr <= 0.1 ~ "Down",
    fdr > 0.1 ~ "NS")) %>%
    mutate(p_simple = ifelse(PValue < 0.001,"<0.001",round(PValue,3))) %>%
    pivot_wider(id_cols = cell,names_from = contrast,values_from = c(direction,p_simple,logFC))

all_iglkhv_res = bind_rows(ighv_family_res,iglkv_family_res)
write_csv(all_iglkhv_res,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/BCR/overall_ighv_iglkv_res.csv")

# light chain family
do_da_var(data = b_cells@meta.data,
        x = "phenotype",
        contrasts_to_test = c("MS_CSF - MS_PBMC",
                              "MS_CSF - OIND_CSF",
                              "MS_CSF - OINDI_CSF"),
        variable="iglightchain_family",
        var_label = "Ig light chain",
        plot_title = "lightchain_family",
        ylim=20)


do_da_var(data = b_cells@meta.data %>% mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE)
          ,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                  "MS_CSF - OINDI_CSF"),
          variable="lightchain_call_simple",
          plot_title = "lightchain_eachclone",
          var_label = "Ig light chain",
          ylim=20)


do_da_var(data = b_cells@meta.data
          ,
          x = "phenotype",
          contrasts_to_test = c("MS_CSF - MS_PBMC",
                                "MS_CSF - OIND_CSF",
                                  "MS_CSF - OINDI_CSF"),
          variable="lightchain_call_simple",
          plot_title = "lightchain",
          var_label = "Ig light chain",
          ylim=20)


#######################################
#  Isotypes
#######################################


do_da_var(data = b_cells@meta.data,
x = "phenotype",
contrasts_to_test = default_contrasts,
variable="clean_isotype",
plot_title = "isotypes",
var_label = "Isotype",
ylim=50)


do_da_var(data = b_cells@meta.data %>%
  mutate(donor_clone_source = paste0(donor_clone,"_",source)) %>% distinct(donor_clone_source,.keep_all=TRUE),
          x = "phenotype",
contrasts_to_test = default_contrasts,
variable="clean_isotype",
plot_title = "isotypes_eachclone",
var_label = "Isotype",
ylim=50)


# cell types
do_da_var(data = b_cells@meta.data %>% filter(cell_type == "Naive B cells"),
x = "phenotype",
contrasts_to_test = default_contrasts,
variable="clean_isotype",
plot_title = "naive_isotypes",
var_label = "Isotype",
ylim=50)

# cell types
do_da_var(data = b_cells@meta.data %>% filter(cell_type == "Memory B cells"),
x = "phenotype",
contrasts_to_test = default_contrasts,
variable="clean_isotype",
plot_title = "mem_isotypes",
var_label = "Isotype",
ylim=50)

do_da_var(data = b_cells@meta.data %>% filter(cell_type == "Plasma cells"),
x = "phenotype",
contrasts_to_test = default_contrasts,
variable="clean_isotype",
plot_title = "pc_isotypes",
var_label = "Isotype",
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

binary_comparison(variable = "shm_positive",plot_title="SHM")


##############################
# mutational load
##############################

do_cont_comparison = function(variable, plot_title,axislab){
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
    labs(x="Cell type",y=axislab)


  png(file=paste0(plot_title,"_comparison.png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
}

do_cont_comparison2 = function(variable, plot_title,axislab,pheno){
  pvals = list()

  dat = b_cells@meta.data %>%
    filter(phenotype==pheno & cell_type %in% cells)
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
    labs(x="Cell type",y=axislab)


  png(file=paste0(plot_title,"_comparison.png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
}

do_cont_comparison(variable = "heavychain_mu_freq_cdr_r",plot_title = "Replacement mutations",axislab = "Mutational load")
do_cont_comparison(variable = "heavychain_mu_freq_cdr_s",plot_title = "Silent mutations",axislab = "Mutational load")

do_cont_comparison(variable = "cdr3_length",plot_title = "CDR3 length",axislab = "CDR3 length")

do_cont_comparison2(variable = "cdr3_length",plot_title = "CDR3 length (MS)",axislab = "CDR3 length",pheno="MS")
do_cont_comparison2(variable = "cdr3_length",plot_title = "CDR3 length (OINDI)",axislab = "CDR3 length",pheno="OINDI")
do_cont_comparison2(variable = "cdr3_length",plot_title = "CDR3 length (OIND)",axislab = "CDR3 length",pheno="OIND")
do_cont_comparison2(variable = "cdr3_length",plot_title = "CDR3 length (NIND)",axislab = "CDR3 length",pheno="NIND")

p=ggplot(b_cells@meta.data %>% filter(cell_type %in% cells),
aes(phenotype,cdr3_length,fill=cell_type))+
facet_wrap(~source)+
geom_boxplot()+
scale_fill_brewer(palette = "Set1")+
theme_minimal()+
labs(x="Phenotype",y="CDR3 length (AA)")
png("cdr3_all_cohorts.png",res=900,units="in",height=3,width=6)
p
dev.off()


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

# define DA function for clonal vs non-clonal
do_da_var_clonal = function(x,variable,ylim=10, data = b_cells@meta.data,plot_title){

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
    contrast_name = comparison_label
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
data = b_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_celltypes_ms",
ylim=50)

do_da_var_clonal(x="cell_type",
variable="cell_type",
data = b_cells@meta.data,
plot_title="clonal_celltypes_allpheno",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_ighv_ms",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data,
plot_title="clonal_ighv_allpheno",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="ighv_family",
data = b_cells@meta.data,
plot_title="clonal_ighv_fam_allpheno",
ylim=50)

# by cell type
do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="MS" & cell_type=="Plasma cells"),
plot_title="clonal_ighv_ms_pcs",
ylim=50)
do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="MS" & cell_type=="Memory B cells"),
plot_title="clonal_ighv_ms_memb",
ylim=50)
do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="MS" & cell_type =="Naive B cells"),
plot_title="clonal_ighv_ms_naive",
ylim=50)


do_da_var_clonal(x="clean_isotype",
variable="clean_isotype",
data = b_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_isotypes_ms",
ylim=50)

do_da_var_clonal(x="lightchain_call_simple",
variable="lightchain_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="MS"),
plot_title="clonal_iglight_ms",
ylim=20)

do_da_var_clonal(x="lightchain_call_simple",
variable="lightchain_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="OINDI"),
plot_title="clonal_iglight_oindi",
ylim=20)

do_da_var_clonal(x="lightchain_call_simple",
variable="lightchain_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="OIND"),
plot_title="clonal_iglight_oind",
ylim=20)


# repeat for OINDI
do_da_var_clonal(x="cell_type",
variable="cell_type",
data = b_cells@meta.data %>% filter(phenotype=="OINDI"),
plot_title="clonal_celltypes_OINDI",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="OINDI"),
plot_title="clonal_ighv_OINDI",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="OIND"),
plot_title="clonal_ighv_OIND",
ylim=50)

do_da_var_clonal(x="v_call_simple",
variable="v_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="NIND"),
plot_title="clonal_ighv_NIND",
ylim=50)

do_da_var_clonal(x="clean_isotype",
variable="clean_isotype",
data = b_cells@meta.data %>% filter(phenotype=="OINDI"),
plot_title="clonal_isotypes_OINDI",
ylim=50)

do_da_var_clonal(x="lightchain_call_simple",
variable="lightchain_call_simple",
data = b_cells@meta.data %>% filter(phenotype=="OINDI"),
plot_title="clonal_iglight_OINDI",
ylim=50)



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
b_cells@meta.data %>% nrow
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
donors_with_expanded_clones = b_cells@meta.data  %>%
  group_by(iid) %>%
  dplyr::count(donor_clone) %>%
  filter(n>1) %>%
  distinct(iid)
b_cells@meta.data = b_cells@meta.data %>%
  mutate(donor_with_expanded_clones = ifelse(iid %in% donors_with_expanded_clones$iid,"Yes","No"))

donors_with_expanded_clones_csf = b_cells@meta.data  %>%
  filter(source=="CSF") %>%
  group_by(iid) %>%
  dplyr::count(donor_clone) %>%
  filter(n>1) %>%
  distinct(iid)
b_cells@meta.data = b_cells@meta.data %>%
  mutate(donor_with_expanded_clones_csf = ifelse(iid %in% donors_with_expanded_clones_csf$iid,"Yes","No"))

# get some counts
b_cells@meta.data %>%
  group_by(phenotype) %>%
  distinct(iid,.keep_all=T) %>%
  dplyr::count(donor_with_expanded_clones) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(donor_with_expanded_clones=="Yes")

# get some counts
# read in main dataset
counts_per_person = read.csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/cell_counts_per_person.csv")
codex = read.csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/pat_codex.csv")
counts_per_person = counts_per_person %>%
  left_join(codex,by="fully_anonymous_pseudoid") %>%
  filter(iid %in% b_cells@meta.data$iid)
# get denominator (i.e. no. of people with any csf cells)
summ_clonal = counts_per_person %>%
  filter(!is.na(total_cells_CSF_B.cells)) %>%
  tibble() %>%
  mutate(donor_with_expanded_clones_csf = ifelse(iid %in% donors_with_expanded_clones_csf$iid,"Yes","No")) %>%
  group_by(Category) %>%
  dplyr::count(donor_with_expanded_clones_csf) %>%
  mutate(total = sum(n), pct = 100*n/total) %>%
  filter(donor_with_expanded_clones_csf == "Yes")
write_csv(summ_clonal,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/clonal_props.csv")


# histogram
p=ggplot(b_cells@meta.data,
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



dat = b_cells@meta.data %>%
  filter(cell_type %in% cells) %>%
  group_by(cell_type,iid,source,phenotype,Age,Sex,OCB,ms_subtype) %>%
  dplyr::count(expanded_clone) %>%
  mutate(total = sum(n), prop = n/sum(n)) %>%
  filter(total>1 & expanded_clone=="Expanded")

dat2 = b_cells@meta.data %>%
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
rownames(b_cells@meta.data) = colnames(b_cells)
p=DimPlot(subset(b_cells,source=="CSF" & phenotype=="MS"),group.by="expanded_clone",split.by="iid",ncol=9)+
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

rownames(b_cells@meta.data) = colnames(b_cells)
# compare expanded vs non-expanded in MS
common_cell_types = c("Naive B cells","Memory B cells","Plasma cells")
expanded_v_not_plot = function(x){
  DimPlot(subset(b_cells,phenotype=="MS" & cell_type %in% common_cell_types),
          split.by="expanded_clone",group.by=x)+
    scale_color_brewer(palette="Set1")+
    theme_umap()
}

p0=expanded_v_not_plot("cell_type")+ggtitle("Cell type")
p1=expanded_v_not_plot("clean_isotype")+ggtitle("Isotype")
p2=expanded_v_not_plot("shm_positive")+ggtitle("SHM")
p3=expanded_v_not_plot("ighv_family")+ggtitle("IGHV")
png("expanded_v_not_bcrs_just_ms.png",res=300,units="in",height=5,width=7)
grid.arrange(p0,p1,p2,p3)
dev.off()

expanded_v_not_plot_csf = function(x){
  DimPlot(subset(b_cells,phenotype=="MS" & cell_type %in% common_cell_types & source=="CSF"),
          split.by="expanded_clone",group.by=x)+
    scale_color_brewer(palette="Set1")+
    theme_umap()
}

p0=expanded_v_not_plot_csf("cell_type")+ggtitle("Cell type")
p1=expanded_v_not_plot_csf("isotype")+ggtitle("Isotype")
p2=expanded_v_not_plot_csf("shm_positive")+ggtitle("SHM")
p3=expanded_v_not_plot_csf("ighv_family")+ggtitle("IGHV")
png("expanded_v_not_bcrs_just_ms_just_csf.png",res=300,units="in",height=5,width=7)
grid.arrange(p0,p1,p2,p3)
dev.off()

# bar plots

png("./expanded_plot_celltypes_barplot.png",res=300,width=7,height=4,units="in")
plot_dat = b_cells@meta.data %>%
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

png("./expanded_plot_isotypes_barplot.png",res=300,width=7,height=4,units="in")
plot_dat = b_cells@meta.data %>%
  mutate(phenotype_new = ifelse(as.character(phenotype) == "OINDI","ID",as.character(phenotype))) %>%
  mutate(expanded_clone = ifelse(expanded_clone == "Expanded","Y","N"))
plot_dat$phenotype_new = factor(plot_dat$phenotype_new,levels = c("NIND","OIND","ID","MS"),ordered=T)
ggplot(plot_dat,
  aes(expanded_clone,fill=clean_isotype))+
  geom_bar(position="fill",color="black")+
  facet_grid(source~phenotype_new)+
  scale_fill_manual(values = colour_pal)+
  theme_minimal()+
  labs(x="Expanded clone?",y="Proportion of\nB/plasma cells",fill="Isotype")
dev.off()

png("./expanded_plot_ighv_fam_barplot.png",res=300,width=7,height=4,units="in")
plot_dat = b_cells@meta.data %>%
  mutate(phenotype_new = ifelse(as.character(phenotype) == "OINDI","ID",as.character(phenotype))) %>%
  mutate(expanded_clone = ifelse(expanded_clone == "Expanded","Y","N"))
plot_dat$phenotype_new = factor(plot_dat$phenotype_new,levels = c("NIND","OIND","ID","MS"),ordered=T)
ggplot(plot_dat,
  aes(expanded_clone,fill=ighv_family))+
  geom_bar(position="fill",color="black")+
  facet_grid(source~phenotype_new)+
  scale_fill_manual(values = colour_pal)+
  theme_minimal()+
  labs(x="Expanded clone?",y="Proportion of\nB/plasma cells",fill="IGHV family")
dev.off()


png("./expanded_plot_ighlkv_fam_barplot.png",res=300,width=7,height=4,units="in")
plot_dat = b_cells@meta.data %>%
  mutate(phenotype_new = ifelse(as.character(phenotype) == "OINDI","ID",as.character(phenotype))) %>%
  mutate(expanded_clone = ifelse(expanded_clone == "Expanded","Y","N"))
plot_dat$phenotype_new = factor(plot_dat$phenotype_new,levels = c("NIND","OIND","ID","MS"),ordered=T)
ggplot(plot_dat,
  aes(expanded_clone,fill=iglightchain_family))+
  geom_bar(position="fill",color="black")+
  facet_grid(source~phenotype_new)+
  scale_fill_manual(values = colour_pal)+
  theme_minimal()+
  labs(x="Expanded clone?",y="Proportion of\nB/plasma cells",fill="IGHV family")
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

expanded %>%
  dplyr::count(donor_clone,phenotype,source,cell_type,isotype,shm_positive) %>%
  arrange(desc(n)) %>% head(5)
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
clone_plot(clone_to_plot = "PatID_65_97_10_2_2215")

# find clones with no expansion
naive_type_clones = b_cells@meta.data %>%
  filter(expanded_clone=="Expanded") %>%
  group_by(donor_clone,phenotype,iid) %>%
  dplyr::count(shm_positive,isotype,source) %>%
  mutate(prop = n/sum(n)) %>%
  filter(shm_positive=="No" & isotype %in% c("IgD","IgM")) %>%
  filter(prop == 1) %>%
  ungroup %>%
  distinct(donor_clone,phenotype,iid,source,n)
naive_type_clones %>% distinct(iid,.keep_all=T) %>% dplyr::count(phenotype)

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

library(circlize)

public_clones = b_cells@meta.data %>%
  group_by(clone_id) %>%
  dplyr::count(iid)  %>%
  dplyr::count(clone_id) %>%
  filter(n>1)

overall_dat = data.frame()
for(i in c(1:nrow(public_clones))){
  dat = b_cells@meta.data %>%
    mutate(iid = paste0(iid,"_",phenotype)) %>%
    filter(clone_id %in% public_clones$clone_id[i]) %>%
    dplyr::count(iid,clone_id) %>%
    tidyr::pivot_wider(id_cols = clone_id,
                       values_from=iid,
                       names_from=iid)
  colnames(dat)[c(2,3)] = c("from","to")
  overall_dat <<- bind_rows(overall_dat,dat)
}
overall_dat$value = 1
overall_dat$clone_id = NULL

grid.cols = b_cells@meta.data %>%
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

# lanz
lanz_cdr3 = read_csv("~/rds/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/Lanz_AllCsfCdr3ForBenJacobs.csv")
b_cells@meta.data %>% filter(junction_aa_VDJ %in% lanz_cdr3$`JUNCTION | HC`)

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
  left_join(overlaps %>% arrange(desc(lnh)) %>% distinct(this_cdr3,.keep_all=T) %>% dplyr::rename("junction_aa_VDJ"=this_cdr3),by="junction_aa_VDJ") %>%
  left_join(lanz_cdr3 %>% dplyr::rename("cdr3_to_match" = `JUNCTION | HC`),by="cdr3_to_match") %>%
  distinct(cell_id,.keep_all=T)

write_csv(overlap_detail,"overlaps.csv")

# write function to clean column in Lanz data
clean_col = function(x){
  y = str_remove_all(str_remove_all(x,"Homsap ")," F")
  z = stringr::str_split(y,"\\*")[[1]][1]
  z = stringr::str_remove_all(z,"D")
  z = paste0(stringr::str_split(z,"\\-")[[1]][1],
  "-",
  stringr::str_split(z,"\\-")[[1]][2])
  z
}
clean_col2 = function(x){
  y = str_remove_all(str_remove_all(x,"Homsap ")," F")
  z = stringr::str_split(y,"\\*")[[1]][1]
  z = paste0(stringr::str_split(z,"\\-")[[1]][1],
  "-",
  stringr::str_split(z,"\\-")[[1]][2])
  z
}
clean_col3 = function(x){
  y = str_remove_all(str_remove_all(x,"Homsap ")," F")
  z = stringr::str_split(y,"\\*")[[1]][1]
  z
}

overlap_detail$ighv_lanz = sapply(overlap_detail$`V.GENE.and.allele | HC`,clean_col)
overlap_detail$igv_light_lanz = sapply(overlap_detail$`V.GENE.and.allele | LC`,clean_col)
overlap_detail$ighd_lanz = sapply( overlap_detail$`D.GENE.and.allele | HC`,clean_col2)
overlap_detail$ighj_lanz = sapply( overlap_detail$`J.GENE.and.allele | HC`,clean_col3)
overlap_detail$igj_light_lanz = sapply(overlap_detail$`J.GENE.and.allele | LC`,clean_col3)


overlap_detail %>%
  filter(v_call_simple==ighv_lanz &
  lightchain_call_simple == igv_light_lanz &
  d_call_VDJ == ighd_lanz &
  j_call_VDJ == ighj_lanz &
  j_call_VJ == igj_light_lanz
   )

################################
# de clonal vs not
################################

b_cells@meta.data = b_cells@meta.data %>% mutate(donor_expanded = paste0(iid,"_",expanded_clone))

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
    left_join(b_cells@meta.data %>%
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


do_de(
  dat = subset(b_cells,phenotype=="MS" & source == "CSF" & cell_type == "Plasma cells"),
  plot_title = "MS_CSF_PCs"
)


do_de(
  dat = subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type == "Memory B cells"),
  plot_title = "MS_CSF_MemB"
  )

do_de(
  dat = subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type == "Memory B cells"),
  plot_title = "MS_CSF_MemB",
  mincells = 10
  )

do_de(
  dat = subset(b_cells,phenotype=="MS" & source == "PBMC" & cell_type == "Plasma cells"),
  plot_title = "MS_PBMC_PCs",
  mincells = 10
)


do_de(
  dat = subset(b_cells,phenotype=="OIND" & source=="CSF" & cell_type == "Memory B cells"),
  plot_title = "OIND_CSF_MemB"
)

do_de(
  dat = subset(b_cells,phenotype=="OINDI" & source=="CSF" & cell_type == "Memory B cells"),
  plot_title = "OINDI_CSF_MemB"
)
do_de(
  dat = subset(b_cells,phenotype=="OIND" & source=="CSF" & cell_type == "Plasma cells"),
  plot_title = "OIND_CSF_PCs"
)

do_de(
  dat = subset(b_cells,source=="CSF" & cell_type == "Plasma cells"),
  plot_title = "all_pheno_CSF_PCs"
)

do_de(
  dat = subset(b_cells,source=="CSF" & cell_type == "Memory B cells"),
  plot_title = "all_pheno_CSF_MemB"
)

# sense check
png("clonal_memb_umap.png",res=300,units="in",width=4,height=4)
DimPlot(subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type=="Memory B cells"),group.by="expanded_clone")+
  theme_umap()
dev.off()

png("clonal_memb_umap_crude_clusters.png",res=300,units="in",width=4,height=4)
DimPlot(subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type=="Memory B cells"),group.by="cell_type_crude",split.by="expanded_clone")+
  theme_umap()
dev.off()


png("sub1_ms_csf.png",res=300,units="in",width=6,height=4)
FeaturePlot(subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type=="Memory B cells"),split.by="expanded_clone",features="SUB1")
dev.off()

png("sub1_oindi_csf.png",res=300,units="in",width=6,height=4)
FeaturePlot(subset(b_cells,phenotype=="OINDI" & source=="CSF" & cell_type=="Memory B cells"),split.by="expanded_clone",features="SUB1")
dev.off()

# heatmap

top_features = read.csv("MS_CSF_MemB_de_expanded_vs_not.csv") %>%
  filter(logFC >0.5) %>%
  filter(P_adj<0.05) %>%
  filter(!grepl("^IG",gene))

selected_features = top_features$gene
DefaultAssay(b_cells)="SCT"

dat = subset(b_cells,phenotype=="MS" & source=="CSF" & cell_type == "Memory B cells")

png("selected_featureplot_clonal_markers.png",res=300,units="in",width=8,height=16)
FeaturePlot(b_cells,features=selected_features,split.by="expanded_clone")
dev.off()

png("selected_featureplot_clonal_markers_just_expanded_memb.png",res=300,units="in",width=8,height=8)
FeaturePlot(dat,features=selected_features,split.by="expanded_clone")
dev.off()



# plot sub1 vs other markers
markers_dat = b_cells@assays$SCT[rownames(b_cells@assays$SCT) %in% c("SUB1","CAPZB","ARPC5","CD27","IGHD","IGHM","IGHG1","CD38","TCL1A","MS4A1"),] %>%
t() %>%
data.frame()


library(corrplot)
corr_plot_sub1 = corrplot(cor(markers_dat,method="spearman"), order = 'hclust', tl.col = 'black',
cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10),diag=F,addCoef.col = 'black')
png("corrplot_sub1_clonal_markers.png",res=300,units="in",width=6,height=6)
print(corr_plot_sub1)
dev.off()

# plot vs shm
sub1_expression =  data.frame(b_cells@assays$RNA["SUB1",] %>% t())
sub1_expression$cell_id = rownames(sub1_expression)
sub1_dat = b_cells@meta.data %>% left_join(sub1_expression,by="cell_id")

ggplot(sub1_dat %>% filter(Category=="MS"),aes(expanded_clone,SUB1))+geom_point()+facet_grid(source~cell_type)




# gene score
clonal_sig = top_features$gene
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
FeaturePlot(b_cells,features = "gene_score_z",split.by="expanded_clone")+
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
    a = b_cells@meta.data %>% filter(expanded_clone=="Not expanded")  %>% filter(phenotype=="MS") %>% filter(source == this_source  & cell_type==cell)
    b = b_cells@meta.data %>% filter(expanded_clone=="Not expanded")  %>% filter(phenotype=="OIND") %>% filter(source == this_source & cell_type==cell)
    p = ifelse(nrow(a)==0 | nrow(b)==0,NA,t.test(a$gene_score_z,b$gene_score_z)$p.value)
    df = data.frame(this_source,cell,p,comparison = "MS_OIND",mean_a = mean(a$gene_score_z),mean_b = mean(b$gene_score_z))
    a = b_cells@meta.data %>% filter(expanded_clone=="Not expanded")  %>% filter(phenotype=="MS") %>% filter(source == this_source  & cell_type==cell)
    b = b_cells@meta.data %>% filter(expanded_clone=="Not expanded")  %>% filter(phenotype=="OINDI") %>% filter(source == this_source & cell_type==cell)
    p = ifelse(nrow(a)==0 | nrow(b)==0,NA,t.test(a$gene_score_z,b$gene_score_z)$p.value)
    df2 = data.frame(this_source,cell,p,comparison = "MS_OINDI",mean_a = mean(a$gene_score_z),mean_b = mean(b$gene_score_z))
    a = b_cells@meta.data %>% filter(expanded_clone=="Not expanded")  %>% filter(phenotype=="OIND") %>% filter(source == this_source  & cell_type==cell)
    b = b_cells@meta.data %>% filter(expanded_clone=="Not expanded")  %>% filter(phenotype=="OINDI") %>% filter(source == this_source & cell_type==cell)
    p = ifelse(nrow(a)==0 | nrow(b)==0,NA,t.test(a$gene_score_z,b$gene_score_z)$p.value)
    df3 = data.frame(this_source,cell,p,comparison = "OIND_OINDI",mean_a = mean(a$gene_score_z),mean_b = mean(b$gene_score_z))

    outputs[[length(outputs)+1]] = data.frame(bind_rows(df,df2,df3))
  }
}
outputs = do.call("bind_rows",outputs)
outputs %>%
  mutate(fdr = p.adjust(p,method="fdr")) %>%
  mutate(sig = ifelse(fdr<0.1,"yes","no")) %>%
  mutate(up_in_ms = ifelse(mean_a>mean_b,"up in MS"," ")) %>%
  dplyr::select(-mean_a,-mean_b,-p) %>%
  filter(cell=="Memory B cells")

p = ggplot(b_cells@meta.data %>%
             filter(expanded_clone=="Not expanded" & phenotype %in% c("OIND","OINDI","MS")) %>%
             filter(cell_type == "Memory B cells"),
           aes(phenotype,gene_score_z,fill=source))+
  geom_boxplot()+
  theme_minimal()+
  labs(x="Phenotype",y="Clonal gene signature \nZ score")

png("clonal_gene_score.png",res=300,units="in",width=8,height=4)
p
dev.off()

png("clonal_gene_score_vs_size.png",res=300,units="in",width=8,height=4)
FeatureScatter(b_cells,feature1 = "gene_score_z", feature2 = "clonal_size",group.by="cell_type")
dev.off()
png("clonal_gene_score_vs_mki67.png",res=300,units="in",width=8,height=4)
FeatureScatter(b_cells,feature1 = "gene_score_z", feature2 = "MKI67",group.by="cell_type")
dev.off()



#######################################
# DOROTHEA
#######################################

library(decoupleR)

# initialise results lists
n_perm = 10000

# grab progeny data
net = get_dorothea(organism='human', levels=c('A', 'B', 'C'))

# read in DE file
de = read_csv("MS_CSF_MemB_de_expanded_vs_not.csv")

# format for decoupleR
de_mat = matrix(de$logFC)
rownames(de_mat) = de$gene

# run decoupleR
acts = run_wmean(mat=de_mat, net=net, .source='source', .target='target',
                  .mor='mor', times = n_perm, minsize = 10)

# filter to just norm_wmean
acts = acts %>% filter(statistic=="norm_wmean")

rfx5 = de %>%
left_join(net %>%
filter(source=="RFX5") %>%
dplyr::rename("gene" = target),
by="gene") %>% filter(!is.na(mor))

png("rfx5.png",res=600,units="in",height=4,width=4)
ggplot(rfx5,aes(logFC,mor,label = gene,size=-log(PValue)))+geom_point()+
geom_text_repel()+
theme_minimal()
dev.off()

# run progeny
net = get_progeny(organism = 'human', top = 1000)

acts = run_wmean(mat=de_mat, net=net, .source='source', .target='target',
                  .mor='weight', times = n_perm, minsize = 10)

# filter to just norm_wmean
acts = acts %>% filter(statistic=="norm_wmean")




################################
# individual genes
################################

png("hcst.png",res=600,units="in",width=4,height=4)
FeaturePlot(dat,features="HCST")+
  scale_colour_gradient(low="darkblue",high="orange")+
  theme_umap()
dev.off()

png("sub1.png",res=600,units="in",width=4,height=4)
FeaturePlot(dat,features="SUB1")+
  scale_colour_gradient(low="darkblue",high="orange")+
  theme_umap()
dev.off()
