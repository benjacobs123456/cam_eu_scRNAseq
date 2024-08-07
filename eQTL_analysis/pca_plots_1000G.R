args = commandArgs(TRUE)

require(ggplot2)
require(RColorBrewer)
library(gridExtra)
library(plyr)

dat = read.table(args[1])
colnames(dat)[1:11] = c("sample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
ref = read.table(args[2],h=T)
dataset = args[3]
outfile = args[4]

# Cohort assignment
for (i in 1:nrow(dat)){
	if(grepl("TUM_",dat[i,"sample"]) == F && grepl("GT_",dat[i,"sample"]) == F){
	ids = unlist(strsplit(dat[i,"sample"],split="_"))
	dat[i,"sample"] = ids[1]
}
}
dat = join(dat,ref,by="sample",type="left")
if(dataset == "ALL"){
	dat$pop = dat$super_pop
}
dat[is.na(dat$pop) == T & grepl("TUM_",dat$sample) == F,"pop"] = "CAM"
dat[is.na(dat$pop) == T & grepl("TUM_",dat$sample) == T,"pop"] = "TUM"
if(dataset == "EUR"){
	dat$pop = factor(dat$pop, ordered=TRUE,levels=c("CEU","FIN","GBR","IBS","TSI","CAM","TUM"))
}
if(dataset == "ALL"){
	dat$pop = factor(dat$pop, ordered=TRUE,levels=c("AFR","AMR","EAS","EUR","SAS","CAM","TUM"))
}


# Scaling
for (j in 2:11){
	dat[,j] = scale(dat[,j])
}

# Colors
colours <- length(levels(dat[,which(colnames(dat)=="pop")]))
if (colours<3) colours <- 3
if(colours>9)
  {
  enh_palette <- colorRampPalette(brewer.pal(9, "Set1"))
  c_values <- rev(enh_palette(colours+1))[2:(colours+1)]
  } else {c_values <- brewer.pal(colours, "Set1")}

if (dataset == "EUR"){
	group.colors <- c(CEU = c_values[1], FIN = c_values[2], GBR =c_values[3], IBS =c_values[4], TSI = c_values[5], CAM="grey70", TUM="grey30")
}
if (dataset == "ALL"){
	group.colors <- c(AFR = c_values[1], AMR = c_values[2], EAS =c_values[3], EUR =c_values[4], SAS = c_values[5], CAM="grey70", TUM="grey30")
}


dat = dat[order(dat$pop),]
p = ggplot(data=dat, aes(x=PC1, y=PC2, colour=pop)) +geom_point(size=I(2),alpha=0.6) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_colour_manual(values=group.colors) +geom_text(aes(label=pop),check_overlap=T) + xlim(min(dat$PC1),max(dat$PC1)) +ylim(min(dat$PC2),max(dat$PC2))
dat1 = dat[dat$pop %in% c("TUM","CAM") == F,]
p1 = ggplot(data=dat1, aes(x=PC1, y=PC2, colour=pop)) +geom_point(size=I(2),alpha=0.6) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_colour_manual(values=group.colors) +geom_text(aes(label=pop),check_overlap=T) + xlim(min(dat$PC1),max(dat$PC1)) +ylim(min(dat$PC2),max(dat$PC2))
dat2 = dat[dat$pop %in% c("TUM","CAM") == T,]
p2 = ggplot(data=dat2, aes(x=PC1, y=PC2, colour=pop)) +geom_point(size=I(2),alpha=0.6) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_colour_manual(values=group.colors) +geom_text(aes(label=pop),check_overlap=T) + xlim(min(dat$PC1),max(dat$PC1)) +ylim(min(dat$PC2),max(dat$PC2))

png(outfile,height=14,width=25,res=300,units="cm")
grid.arrange(p,p1,p2,nrow=1)
dev.off()

