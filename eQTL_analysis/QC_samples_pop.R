args <- commandArgs(TRUE)


source("/xx/xxx/SC/CAM_TUM/eqtl/Scripts/clusterAnalysis.r")
source("/xx/xxx/SC/CAM_TUM/eqtl/Scripts/MDS_plots.r")
require(ggplot2)
require(RColorBrewer)
library(gridExtra)
library(grid)


mds <- read.table(args[1], h=F)
colnames(mds) = c("IID","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
FAM <- read.table(args[2])
outfiles = args[3]
outplot = args[4]


colnames(FAM) <- c("IID", "SEX")
mergedMDS <- merge(FAM,mds,by="IID")

threshold <- 4
outList1 <- removeOutlier(mergedMDS, threshold, save=F, colour="sex")
print("Number of population outliers (> 4 SD from the mean in C1 and C2):")
print(nrow(outList1$data[,1:2]))
write.table(outList1$data[,1],paste0(outfiles,"1.txt",collapse=""),c=F,r=F,qu=F)
print(paste0(outfiles,"1.txt",collapse=""))

outList2 <- removeOutlier(mergedMDS, threshold, save=F, colour="sex", C1="C3", C2="C4")
print("Number of population outliers (> 4 SD from the mean in C3 and C4):")
print(nrow(outList2$data[outList2$data[,2] %in% outList1$data[,2] == FALSE,1:2]))
write.table(outList2$data[,1],paste0(outfiles,"2.txt",collapse=""),c=F,r=F,qu=F)

outList3 <- removeOutlier(mergedMDS, threshold, save=F, colour="sex", C1="C5", C2="C6")
print("Number of population outliers (> 4 SD from the mean in C5 and C6):")
print(nrow(outList3$data[outList3$data[,2] %in% outList1$data[,2] == FALSE & outList3$data[,2] %in% outList2$data[,2] == FALSE,1:2]))
write.table(outList3$data[,1],paste0(outfiles,"3.txt",collapse=""),c=F,r=F,qu=F)

out1 <- outList1$data[,1:2]
out2 <- outList2$data[,1:2]
out3 <- outList3$data[,1:2]
outpop <- unique(rbind(out1,out2,out3))

mergedMDS$cohort = "CAM"
mergedMDS[grepl("TUM",mergedMDS$IID) == T,"cohort"] = "TUM"


#print("MDS plot before removal of outliers")
colour="cohort"

# scale MDS components
mergedMDS[,"C1"] <- scale(mergedMDS[,"C1"])
mergedMDS[,"C2"] <- scale(mergedMDS[,"C2"])
mergedMDS[,"C3"] <- scale(mergedMDS[,"C3"])
mergedMDS[,"C4"] <- scale(mergedMDS[,"C4"])
mergedMDS[,"C5"] <- scale(mergedMDS[,"C5"])
mergedMDS[,"C6"] <- scale(mergedMDS[,"C6"])
mergedMDS[,"C7"] <- scale(mergedMDS[,"C7"])
mergedMDS[,"C8"] <- scale(mergedMDS[,"C8"])

# plot
colours <- length(levels(mergedMDS[,which(colnames(mergedMDS)==colour)]))
if (colours<3) colours <- 3
if(colours>9)
{
  enh_palette <- colorRampPalette(brewer.pal(9, "Set1"))
  c_values <- rev(enh_palette(colours+1))[2:(colours+1)]
} else {c_values <- brewer.pal(colours, "Set1")}

# plot
leg_rows=1
pl1 <- ggplot(data=mergedMDS, aes_string(x="C1", y="C2", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC1") +ylab("PC2") +guides(col=guide_legend(nrow=leg_rows,override.aes=list(size=I(2))))
pl2 <- ggplot(data=mergedMDS, aes_string(x="C3", y="C4", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC3") +ylab("PC4")
pl3 <- ggplot(data=mergedMDS, aes_string(x="C5", y="C6", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC5") +ylab("PC6")
pl4 <- ggplot(data=mergedMDS, aes_string(x="C7", y="C8", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC7") +ylab("PC8")


grid_arrange_shared_legend <- function(...)
  {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] +theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x +theme(legend.position="none"))),legend,ncol=1,heights=unit.c(unit(1,"npc")-lheight,lheight))
  }

png(outplot)
grid_arrange_shared_legend(pl1, pl2, pl3, pl4)
dev.off()
