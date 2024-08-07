require(ggplot2)
require(RColorBrewer)

plotMDS <- function(mergedMDS,Cx="C1",Cy="C2",colour="sex",save=T)
  {
  nCx <- which(colnames(mergedMDS)==Cx)
  nCy <- which(colnames(mergedMDS)==Cy)

  # scale MDS components
  mergedMDS[,nCx] <- scale(mergedMDS[,nCx])
  mergedMDS[,nCy] <- scale(mergedMDS[,nCy])

  # plot
  colours <- length(levels(mergedMDS[,which(colnames(mergedMDS)==colour)]))
  if (colours<3) colours <- 3
  if(colours>9)
    {
    enh_palette <- colorRampPalette(brewer.pal(9, "Set1"))
    c_values <- rev(enh_palette(colours+1))[2:(colours+1)]
    } else {c_values <- brewer.pal(colours, "Set1")}

  print(ggplot(data=mergedMDS, aes_string(x=Cx, y=Cy, colour=colour)) +geom_point(size=I(3)) +theme_bw() +
    theme(panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_x_continuous(breaks=seq(-100,100,0.5)) +
    scale_y_continuous(breaks=seq(-100,100,0.5)) +scale_colour_manual(values=c_values))
  if(save) ggsave(paste("MDS_",Cx,"_",Cy,"_",colour,".pdf",sep=""))
  }


plotMDS_Group <- function(mergedMDS,group_ids,group_seq,group_names,Cx="C1",Cy="C2",save=T)
  {
  nCx <- which(colnames(mergedMDS)==Cx)
  nCy <- which(colnames(mergedMDS)==Cy)

  # scale MDS components
  mergedMDS[,nCx] <- scale(mergedMDS[,nCx])
  mergedMDS[,nCy] <- scale(mergedMDS[,nCy])

  # plot according to groups
  if(length(group_ids)>9)
    {
    enh_palette <- colorRampPalette(brewer.pal(9, "Set1"))
    # c_values <- enh_palette(length(group_ids))
    # c_values <- rev(enh_palette(length(group_ids)))
    c_values <- rev(enh_palette(length(group_ids)+1))[2:(length(group_ids)+1)]
    } else {c_values <- brewer.pal(length(group_ids), "Set1")}

  mergedMDS["group"] <- NA
  for(i in 1:length(group_ids))
    {
    mergedMDS[which(grepl(group_ids[i],mergedMDS$FID)),]$group <- group_ids[i]
    }

  mergedMDS$group <- factor(mergedMDS$group, levels=group_seq, labels=group_names)
  mergedMDS <- mergedMDS[with(mergedMDS,order(group)),]

  # plot
  print(ggplot(data=mergedMDS, aes_string(x=Cx, y=Cy, colour="group")) +geom_point(size=I(3)) +theme_bw() +
    theme(panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_x_continuous(breaks=seq(-100,100,0.5)) +
    scale_y_continuous(breaks=seq(-100,100,0.5)) +scale_colour_manual(values=c_values, name="Cohort"))
  if(save) ggsave(paste("MDS_",Cx,"_",Cy,".pdf",sep=""))
  }


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


plotMDS_all <- function(mergedMDS,colour="sex",leg_rows=1,save=T)
  {
  require(grid)
  require(gridExtra)

  # scale MDS components
  nC <- grep("^C[0-9]",colnames(mergedMDS))

  for(i in nC)
    {
    mergedMDS[,i] <- scale(mergedMDS[,i])
    }

  colours <- length(levels(mergedMDS[,which(colnames(mergedMDS)==colour)]))
  if (colours<3) colours <- 3
  if(colours>9)
    {
    enh_palette <- colorRampPalette(brewer.pal(9, "Set1"))
    c_values <- rev(enh_palette(colours+1))[2:(colours+1)]
    } else {c_values <- brewer.pal(colours, "Set1")}

  # plot
  pl1 <- ggplot(data=mergedMDS, aes_string(x="C1", y="C2", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 1") +ylab("MDS component 2") +guides(col=guide_legend(nrow=leg_rows,override.aes=list(size=I(2))))
  pl2 <- ggplot(data=mergedMDS, aes_string(x="C3", y="C4", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 3") +ylab("MDS component 4")
  pl3 <- ggplot(data=mergedMDS, aes_string(x="C5", y="C6", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 5") +ylab("MDS component 6")
  pl4 <- ggplot(data=mergedMDS, aes_string(x="C7", y="C8", colour=colour)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 7") +ylab("MDS component 8")

  if(save) png(paste("MDS_C1-C8","_",colour,".png",sep=""), type="cairo", width=3200, height=3200, units="px", pointsize=12, res=300)
  grid_arrange_shared_legend(pl1, pl2, pl3, pl4)
  if(save) dev.off()
  }


plotMDS_all_Group <- function(mergedMDS,group_ids,group_seq,group_names,leg_rows=1,save=T) # leg_rows specifies the number of rows for the legend
  {
  require(grid)
  require(gridExtra)

  # scale MDS components
  nC <- grep("^C[0-9]",colnames(mergedMDS))

  for(i in nC)
    {
    mergedMDS[,i] <- scale(mergedMDS[,i])
    }

  # plot according to groups
  if(length(group_ids)>9)
    {
    enh_palette <- colorRampPalette(brewer.pal(9, "Set1"))
    # c_values <- rev(enh_palette(length(group_ids)))
    # c_values <- enh_palette(length(group_ids))
    c_values <- rev(enh_palette(length(group_ids)+1))[2:(length(group_ids)+1)]
    } else {c_values <- brewer.pal(length(group_ids), "Set1")}

  mergedMDS["group"] <- NA
  for(i in 1:length(group_ids))
    {
    mergedMDS[which(grepl(group_ids[i],mergedMDS$FID)),]$group <- group_ids[i]
    }

  mergedMDS$group <- factor(mergedMDS$group, levels=group_seq, labels=group_names)
  mergedMDS <- mergedMDS[with(mergedMDS,order(group)),]

  pl1 <- ggplot(data=mergedMDS, aes(x=C1, y=C2, colour=group)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values, name="Cohort") +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 1") +ylab("MDS component 2") +guides(col=guide_legend(nrow=leg_rows,override.aes=list(size=I(2))))
  pl2 <- ggplot(data=mergedMDS, aes(x=C3, y=C4, colour=group)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values, name="Cohort") +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 3") +ylab("MDS component 4")
  pl3 <- ggplot(data=mergedMDS, aes(x=C5, y=C6, colour=group)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values, name="Cohort") +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 5") +ylab("MDS component 6")
  pl4 <- ggplot(data=mergedMDS, aes(x=C7, y=C8, colour=group)) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="blank",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values, name="Cohort") +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("MDS component 7") +ylab("MDS component 8")

  if(save) png("MDS_C1-C8.png", width=3200, height=3200, units="px", pointsize=12, res=300)
  grid_arrange_shared_legend(pl1, pl2, pl3, pl4)
  if(save) dev.off()
  }
