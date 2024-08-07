
library(grid)
library(gridExtra)
library(ggplot2)
library(plotscale)
library(ragg)
#library(RNOmni)

# directories
data.dir="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
setwd(data.dir)


# Supplementary Figure 1 
load(paste0(data.dir,"/figures/Fig1_3.RData"))
load(paste0(data.dir,"/figures/Fig1_4.RData"))

ragg::agg_png(paste0(data.dir,"/figures/Supplementary_Figure1.png"),units="in",height=3,width=6,res=300,scaling=1)
f1_3 = f1_3  + theme(plot.title = element_text(family="sans",size=8,face = "bold",hjust=0.5), text = element_text(size=8,family="sans"), strip.text.x = element_text(size=8,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=8), axis.text.y = element_text(family="sans",size=8), axis.title.x = element_text(family="sans",size=8,face = "bold"), axis.title.y = element_text(family="sans",size=8,face = "bold"), legend.text=element_text(family="sans",size=8))
f1_4 = f1_4  + theme(plot.title = element_text(family="sans",size=8,face = "bold",hjust=0.5), text = element_text(size=8,family="sans"), strip.text.x = element_text(size=8,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=8), axis.text.y = element_text(family="sans",size=8), axis.title.x = element_text(family="sans",size=8,face = "bold"), axis.title.y = element_text(family="sans",size=8,face = "bold"), legend.text=element_text(family="sans",size=8))
grid.arrange(f1_3,f1_4,nrow=1)
grid.text("A",x=unit(0.1,"in"),y=unit(2.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(3.1,"in"),y=unit(2.9,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
dev.off()

# Main eQTL figure
load(paste0(data.dir,"/figures/zc2hc1a.RData"))
load(paste0(data.dir,"/figures/ets1.RData"))
ragg::agg_png(paste0(data.dir,"/figures/maineQTLfigure_v2.png"),units="in",height=4,width=6,res=300,scaling=0.6)
grid.arrange(ets1[[1]],ets1[[15]],ets1[[13]],zc2hc1a[[1]],zc2hc1a[[15]],zc2hc1a[[13]],nrow=2,ncol=3)
grid.text("A",x=unit(0.25,"in"),y=unit(6.5,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(3.58,"in"),y=unit(6.5,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(6.95,"in"),y=unit(6.5,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("D",x=unit(0.25,"in"),y=unit(3.25,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("E",x=unit(3.58,"in"),y=unit(3.25,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))
grid.text("F",x=unit(6.95,"in"),y=unit(3.25,"in"),gp = gpar(fontsize=10, fontfamily="sans",fontface="bold"))

dev.off()

# Supplementary Figure 2
load(paste0(data.dir,"/figures/PREX1.RData"))
prex1 = plots
ragg::agg_png(paste0(data.dir,"/figures/Supplementary_Figure2.png"),units="in",height=2.1,width=4.2,res=300,scaling=0.55)
grid.arrange(prex1[[1]],prex1[[13]],nrow=1)
grid.text("A",x=unit(0.15,"in"),y=unit(3.7,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.1,"in"),y=unit(3.7,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
dev.off()

# Supplementary Figure 3
library(grid)
load(paste0(data.dir,"/figures/ahi1.RData"))
load(paste0(data.dir,"/figures/eaf2.RData"))
ragg::agg_png(paste0(data.dir,"/figures/Supplementary_Figure3.png"),units="in",height=4.2,width=6.5,res=300,scaling=0.55)
grid.arrange(ahi1[[1]],ahi1[[15]],ahi1[[13]],eaf2[[1]],eaf2[[15]],eaf2[[13]],nrow=2,widths=c(1,1,0.8))
grid.text("A",x=unit(0.15,"in"),y=unit(7.5,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("B",x=unit(4.38,"in"),y=unit(7.5,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("C",x=unit(8.66,"in"),y=unit(7.5,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("D",x=unit(0.15,"in"),y=unit(3.75,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("E",x=unit(4.38,"in"),y=unit(3.75,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
grid.text("F",x=unit(8.66,"in"),y=unit(3.75,"in"),gp = gpar(fontsize=12, fontfamily="sans",fontface="bold"))
dev.off()

# Supp figure cell numbers

# see /xx/xxx/SC/CAM_TUM/eqtl/Targeted/Scripts/numbers_cellspersample.R
