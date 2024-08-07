
# Load packages
library(plyr)
library(tidyselect)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ragg)
library(grid)
library(gridExtra)


met = fread("/xx/xxx/SC/CAM_TUM/eqtl/all_combo_with_updated_pheno_metadata.tsv")
met$cell_type = as.character(met$cell_type)
met[met$cell_type == "Plasma cells","cell_type"] = "ASCs"
met = data.frame(met)

cell_types=c("CD4 T cells","CD8 T cells","Tregs","MAIT cells","NK cells", "B cells","ASCs","CD14 Mono","CD16 Mono","mDCs","pDCs","Macrophages","HSPCs")
my_colors=c("#F7BF13","#E57601","#944E04","#5B5B5B","#CC0000","#ABD3F8","#3879B4","#6A329F","#C90076","#B6D7A8","#016D45","#B4A7D6","#C4C4C4","#000000","#07DA63")

met$cell_type = factor(met$cell_type,ordered=T,levels=cell_types)

met[met$iid == "8A2H3PL1","iid"] = "PatID_5"
met[met$iid == "E9LEH7P8","iid"] = "PatID_56"
met[met$iid == "JEGK54J2","iid"] = "PatID_59"
met[met$iid == "RL4X6288","iid"] = "PatID_3"
met[met$iid == "WAF0FQN5","iid"] = "PatID_52"
met[met$iid == "WQ0EJVR1","iid"] = "PatID_60"

plot_data = data.frame(table(met[,c("iid","source","cell_type")]))
plot_data = join(plot_data,unique(met[,c("iid","phenotype")]),by="iid",type="left")
plot_data$phenotype = as.character(plot_data$phenotype)
plot_data[plot_data$phenotype == "OINDI","phenotype"] = "ID"

plot_data$source = factor(plot_data$source,ordered=T,levels=c("CSF","PBMC"))
plot_data$phenotype = factor(plot_data$phenotype,ordered=T,levels=c("MS","OIND","ID","NIND"))
plot_data$cell_type = factor(plot_data$cell_type,ordered=T,levels=cell_types)

legend_title = "Cell type"

test = data.frame(table(met[,c("iid","source")]))
test= test[test$Freq != 0,]
ind_with_csf = test[test$source == "CSF",]
ind_with_pbmc = test[test$source == "PBMC",]

plot_data_csf  = plot_data[plot_data$source == "CSF",]
plot_data_csf = plot_data_csf[plot_data_csf$iid %in% ind_with_csf$iid == T,]
p1=ggplot(plot_data_csf,aes(iid,Freq,fill=cell_type))+
  geom_col(position=position_dodge())+
  #coord_flip()+
  facet_grid(cell_type ~ phenotype,scales="free")+
  xlab("Study samples")+ylab("N. of cells")+
  scale_fill_manual(legend_title,values=my_colors)+
  ggtitle("CSF")+
  theme(text = element_text(family="sans",size=7))+
   theme(plot.title = element_text(family="sans",size=15,face = "bold",hjust=0.5), text = element_text(size=15,family="sans"), axis.text.x = element_text(family="sans",size=15,angle=90), axis.text.y = element_text(family="sans",size=10), axis.title.x = element_text(family="sans",size=15,face = "bold"), axis.title.y = element_text(family="sans",size=15,face = "bold"), legend.text=element_text(family="sans",size=15))+
  theme(strip.text.x = element_text(family="sans",size=15,face = "bold"),strip.text.y = element_blank())+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank())


  plot_data_pbmc  = plot_data[plot_data$source == "PBMC",]
  plot_data_pbmc = plot_data_pbmc[plot_data_pbmc$iid %in% ind_with_pbmc$iid == T,]
  p2=ggplot(plot_data_pbmc,aes(iid,Freq,fill=cell_type))+
  geom_col(position=position_dodge())+
  #coord_flip()+
  facet_grid(cell_type ~ phenotype,scales="free")+
  xlab("Study samples")+ylab("N. of cells")+
  scale_fill_manual(legend_title,values=my_colors)+
  theme(text = element_text(family="sans",size=7))+
  ggtitle("PBMC")+
  theme(plot.title = element_text(family="sans",size=15,face = "bold",hjust=0.5), text = element_text(size=15,family="sans"), axis.text.x = element_text(family="sans",size=15,angle=90), axis.text.y = element_text(family="sans",size=10), axis.title.x = element_text(family="sans",size=15,face = "bold"), axis.title.y = element_text(family="sans",size=15,face = "bold"), legend.text=element_text(family="sans",size=15))+
  theme(strip.text.x = element_text(family="sans",size=15,face = "bold"),strip.text.y = element_blank())+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank())

plots = list()
plots[[length(plots)+1]] = p1
plots[[length(plots)+1]] = p2

data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
setwd(data.dir)


cnum = plots
ragg::agg_png(paste0(data.dir,"/figures/Supplementary_Figure_cellnumbers.png"),units="in",height=9,width=7.5,res=300,scaling=0.55)
grid.arrange(cnum[[1]],cnum[[2]],nrow=2)
grid.text("A",x=unit(0.15,"in"),y=unit(16.2,"in"),gp = gpar(fontsize=15, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(0.15,"in"),y=unit(8.1,"in"),gp = gpar(fontsize=15, fontfamily="sans",fontface="bold"))
dev.off()