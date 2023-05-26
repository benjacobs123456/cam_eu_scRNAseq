#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(MASS)
library(reshape2)
library(liana)
library(corrplot)


#######################################
# Read in data
#######################################

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/ccc/")

# find files
files = list.files("./per_sample/",full.names = T)

# read in
liana_res_overall = read_csv(files)


# get sc data
sc_dat = readRDS("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/all_combo_with_updated_pheno.rds")

# get metadata
sc_dat_metadata = sc_dat@meta.data
write_csv(sc_dat_metadata,"sc_dat_metadata.csv")
# sc_dat_metadata = read_csv("sc_dat_metadata.csv")

#######################################
# Correlation between cell numbers
#######################################

# subset sc data to just MS CSF
sc_dat_metadata = read_csv("sc_dat_metadata.csv") %>%
  filter(source=="CSF" & Category=="MS")

counts_per_person = sc_dat_metadata %>%
  filter(!cell_type_crude %in% c("HSPCs","Platelets")) %>%
  group_by(iid) %>%
  dplyr::count(cell_type_crude) %>%
  mutate(prop = n/sum(n)) %>%
  tidyr::pivot_wider(id_cols = iid,values_from = prop,names_from=cell_type_crude) %>%
  ungroup()

cor_mat = cor(counts_per_person[,-1],use = "complete.obs")
cor_pvals = corrplot::cor.mtest(counts_per_person[,-1],use = "complete.obs")

# get no. of pairwise comparisons
j=0
for(i in c(11:1)){
  j <<- j + (i-1)
}
pbonf = 0.05/j
library(corrplot)
png("corrplot_celltypes.png",res=600,units="in",width=6,height=6)
corrplot::corrplot(cor_mat,p.mat = cor_pvals$p,
                   number.cex = 0.8, order = 'hclust', diag=FALSE, insig='label_sig',
                   tl.srt = 45, col = COL2('PuOr', 10), type="lower",
                   sig.level = pbonf)
dev.off()

p1=ggplot(counts_per_person,aes(`CD14 Mono`,`B cells`))+
  geom_point()+
  theme_minimal()+
  geom_smooth(method="lm",linetype="dashed",se=F)

p2=ggplot(counts_per_person,aes(`CD14 Mono`,`Plasma cells`))+
  geom_point()+
  theme_minimal()+
  geom_smooth(method="lm",linetype="dashed",se=F)

p3=ggplot(counts_per_person,aes(`B cells`,`Plasma cells`))+
  geom_point()+
  theme_minimal()+
  geom_smooth(method="lm",linetype="dashed",se=F)

png("corrplot_celltypes_2.png",res=600,units="in",width=6,height=4)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

#######################################
# LIANA
#######################################

# summarise interaction
liana_res_overall = liana_res_overall %>%
  mutate(interaction = paste0(source,"_",target,"_",ligand.complex,"_",receptor.complex)) %>%
  mutate(interaction_nocell = paste0(ligand.complex,"_",receptor.complex)) %>%
  mutate(sig = ifelse(aggregate_rank < 0.05, "yes","no"))

#
b_pcs = liana_res_overall %>%
  filter(target %in% c("B cells","Plasma cells"))

# find sig pathways
find_sig_interactions = function(
  pheno = "MS",
  cell_type = "Plasma cells"
){

  # filter
  int_dat = b_pcs %>%
    filter(target==cell_type & phenotype==pheno)
  n_samples = int_dat %>% distinct(iid) %>% nrow
  message("N=",n_samples)

  # find sig interactions
  top_interactions = int_dat %>%
    group_by(interaction) %>%
    dplyr::count(sig) %>%
    mutate(prop_sig = n/sum(n), total = sum(n)) %>%
    filter(sig=="yes") %>%
    arrange(desc(prop_sig))


  int_dat %>%
    filter(interaction %in% top_interactions$interaction) %>%
    group_by(interaction,source,target,ligand.complex,receptor.complex,phenotype) %>%
    summarise(median_LR = median(sca.LRscore)) %>%
    left_join(top_interactions,by="interaction") %>%
    arrange(desc(prop_sig)) %>%
    filter(prop_sig>0.9 & total>0.5*n_samples)


}

res = bind_rows(
  find_sig_interactions("MS","B cells"),
  find_sig_interactions("MS","Plasma cells")
)
print(res,n=30) %>% filter(ligand.complex=="MIF")
find_sig_interactions("OIND","B cells")
find_sig_interactions("OIND","Plasma cells")
find_sig_interactions("NIND","B cells")
find_sig_interactions("NIND","Plasma cells")


# filter to chemokines & receptors
chemokine_list = read_tsv("chemokines.txt")
liana_res_overall = liana_res_overall %>%
  filter(receptor.complex %in% chemokine_list$`Approved symbol`)

# find sig pathways
find_sig_interactions = function(
  pheno = "MS",
  cell_type = "Plasma cells"
){

  # filter
  int_dat = liana_res_overall %>%
    filter(target==cell_type & phenotype==pheno)
  n_samples = int_dat %>% distinct(iid) %>% nrow
  message("N=",n_samples)

  # find sig interactions
  top_interactions = int_dat %>%
    group_by(interaction) %>%
    dplyr::count(sig) %>%
    mutate(prop_sig = n/sum(n), total = sum(n)) %>%
    filter(sig=="yes") %>%
    arrange(desc(prop_sig))


  int_dat %>%
    filter(interaction %in% top_interactions$interaction) %>%
    group_by(interaction,source,target,ligand.complex,receptor.complex,phenotype) %>%
    summarise(median_LR = median(sca.LRscore)) %>%
    left_join(top_interactions,by="interaction") %>%
    arrange(desc(prop_sig))
}


res = bind_rows(
  find_sig_interactions("MS","B cells"),
  find_sig_interactions("MS","Plasma cells"),
  find_sig_interactions("OIND","B cells"),
  find_sig_interactions("OIND","Plasma cells"),
  find_sig_interactions("OINDI","B cells"),
  find_sig_interactions("OINDI","Plasma cells"),
  find_sig_interactions("NIND","B cells"),
  find_sig_interactions("NIND","Plasma cells")
)

res$phenotype = factor(res$phenotype,levels=c("NIND","OINDI","OIND","MS"),ordered = T)
write_csv(res %>% filter(total>10),"liana_res_chemokines.csv")

liana_res_overall %>% filter(receptor.complex=="CCR4" & ligand.complex=="CCL22") %>%
  arrange(aggregate_rank) %>%
  filter(target=="Plasma cells")

just_ms = res %>%
  filter(phenotype=="MS") %>%
  arrange(desc(prop_sig)) %>%
  filter(total>10)
plot_dat = just_ms %>% mutate(Interaction = paste0(ligand.complex," > ", receptor.complex))
p1=ggplot(plot_dat,
          aes(Interaction,source,label=paste0(round(prop_sig*100,1),"%"),fill=100*prop_sig))+
  geom_tile()+
  facet_wrap(~target)+
  geom_text(size=3,color="black")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  scale_fill_gradient(low="purple",high="orange")+
  labs(y="Source of ligand",fill="% of samples\nwith P<0.05")

png("liana_res_test_ms_bcells.png",res=600,units="in",height=6,width=8)
p1
dev.off()


#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(gridExtra)
library(MASS)
library(reshape2)


#######################################
# Read in data
#######################################
library(Seurat)

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/ccc/")

# Read in data
b_cells = readRDS("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/BCR/b_cells_post_processing.rds")
rownames(b_cells@meta.data) = colnames(b_cells)

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

# filter doublets (based on UMAP)
png("test.png")
DimPlot(b_cells,split.by="cell_type_crude")
dev.off()

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


p1=FeaturePlot(b_cells,split.by="source",features=c("CXCR3","CXCR4","IGHD","CD27"))

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


b_cells = subset(b_cells,subset = ann_celltypist_highres %in% c("Memory B cells","Naive B cells","Plasma cells"))
p=DimPlot(subset(b_cells, phenotype=="MS"),group.by="ann_celltypist_highres")+
  scale_color_brewer(palette="Paired")+
  theme_umap()
png("b_cells.png",res=600,units="in",width=3,height=1)
p
dev.off()
png("b_cells_cc_receptors.png",res=600,units="in",width=8,height=8)
p1
dev.off()




# plot
DefaultAssay(b_cells)="RNA"
FeaturePlot(b_cells,features="CCL22")
