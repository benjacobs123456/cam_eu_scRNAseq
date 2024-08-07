require(ggplot2)

removeOutlier <- function(mergedData, threshold, save=T, C1="C1", C2="C2", colour="sex", base="")
  {
  remove <- mergedData

  nC1 <- which(colnames(remove)==C1)
  nC2 <- which(colnames(remove)==C2)

  remove[,nC1] <- scale(remove[,nC1])
  remove[,nC2] <- scale(remove[,nC2])

  cluster1 <- ggplot(data=remove, aes_string(x=C1, y=C2, colour=colour)) +geom_point(size=I(3)) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +scale_colour_manual(values=c("#377EB8","#E41A1C","#4DAF4A"))
  if(save)
    {
    if(base=="") {ggsave(paste("cluster_",C1,"-",C2,"_s0.pdf",sep=""))} else {ggsave(paste("cluster_",C1,"-",C2,"_s",base,"_s0.pdf",sep=""))}
    }

  if(length(remove[which(abs(remove[,nC1])>threshold),]$IID)!=0)
    {
    out <- as.character(remove[which(abs(remove[,nC1])>threshold),]$IID)
    remove <- remove[-which(abs(remove[,nC1])>threshold),]
    remove[,nC1] <- scale(remove[,nC1])
    remove[,nC2] <- scale(remove[,nC2])
    }
  else
    {
    out <- ""
    }

  while(length(which(abs(remove[,nC1])>threshold))!=0)
    {
    out <- append(out, as.character(remove[which(abs(remove[,nC1])>threshold),]$IID))
    remove <- remove[-which(abs(remove[,nC1])>threshold),]
    remove[,nC1] <- scale(remove[,nC1])
    remove[,nC2] <- scale(remove[,nC2])
    }

  if(length(remove[[1]])==0)
    {
    remove <- mergedData
    remove[,nC1] <- scale(remove[,nC1])
    remove[,nC2] <- scale(remove[,nC2])
    }

  if(length(remove[which(abs(remove[,nC2])>threshold),]$IID)!=0)
    {
    out <- append(out, as.character(remove[which(abs(remove[,nC2])>threshold),]$IID))
    remove <- remove[-which(abs(remove[,nC2])>threshold),]
    remove[,nC1] <- scale(remove[,nC1])
    remove[,nC2] <- scale(remove[,nC2])
    }

  while(length(which(abs(remove[,nC2])>threshold))!=0)
    {
    out <- append(out, as.character(remove[which(abs(remove[,nC2])>threshold),]$IID))
    remove <- remove[-which(abs(remove[,nC2])>threshold),]
    remove[,nC1] <- scale(remove[,nC1])
    remove[,nC2] <- scale(remove[,nC2])
    }

  remove <- mergedData[!(mergedData$IID %in% out),]
  remove[,nC1] <- scale(remove[,nC1])
  remove[,nC2] <- scale(remove[,nC2])

  cluster2 <- ggplot(data=remove, aes_string(x=C1, y=C2, colour=colour)) +geom_point(size=I(3)) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +scale_colour_manual(values=c("#377EB8","#E41A1C","#4DAF4A"))
  if(save)
    {
    if(base=="") {ggsave(paste("cluster_",C1,"-",C2,"_s",threshold,".pdf",sep=""))} else {ggsave(paste("cluster_",C1,"-",C2,"_s",base,"_s",threshold,".pdf",sep=""))}
    }

  outData <- mergedData[which(mergedData$IID %in% out),]
  outList <- list(data=outData,plot1=cluster1,plot2=cluster2)
  return(outList)
  }

removeOutlierC1 <- function(mergedData, threshold1, threshold2, save=TRUE, colour="sex", base="")
  {
  remove <- mergedData
  remove$C1 <- scale(remove$C1)
  remove$C2 <- scale(remove$C2)

  cluster1 <- ggplot(data=remove, aes_string(x="C1", y="C2", colour=colour)) +geom_point(size=I(3)) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +scale_colour_manual(values=c("#377EB8","#E41A1C","#4DAF4A"))
  if(save)
    {
    if(base=="") {ggsave("cluster_s0.pdf")} else {ggsave(paste("cluster_s",base,"_s0.pdf",sep=""))}
    }

  if(length(remove[which(remove$C1<threshold1),]$IID)!=0)
    {
    out <- as.character(remove[which(remove$C1<threshold1),]$IID)
    remove <- remove[-which(remove$C1<threshold1),]
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)
    }
  else
    {
    out <- ""
    }

  while(length(which(remove$C1>threshold2))!=0)
    {
    out <- append(out, as.character(remove[which(remove$C1>threshold2),]$IID))
    remove <- remove[-which(remove$C1>threshold2),]
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)
    }

  if(length(remove[[1]])==0)
    {
    remove <- mergedData
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)
    }

  if(length(remove[which(remove$C1>threshold2),]$IID)!=0)
    {
    out <- append(out, as.character(remove[which(remove$C1<threshold1),]$IID))
    remove <- remove[-which(remove$C1>threshold2),]
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)
    }

  while(length(which(remove$C1<threshold1))!=0)
    {
    out <- append(out, as.character(remove[which(remove$C1<threshold1),]$IID))
    remove <- remove[-which(remove$C1<threshold1),]
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)
    }

  if(length(remove[[1]])==0)
    {
    remove <- mergedData
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)
    }

  remove <- mergedData[!(mergedData$IID %in% out),]
  remove$C1 <- scale(remove$C1)
  remove$C2 <- scale(remove$C2)

  cluster2 <- ggplot(data=remove, aes_string(x="C1", y="C2", colour=colour)) +geom_point(size=I(3)) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +scale_colour_manual(values=c("#377EB8","#E41A1C","#4DAF4A"))
  if(save)
    {
    if(base=="") {ggsave(paste("cluster_s",threshold,".pdf",sep=""))} else {ggsave(paste("cluster_s",base,"_C1_s",threshold,".pdf",sep=""))}
    }

  outData <- mergedData[which(mergedData$IID %in% out),]
  outList <- list(data=outData,plot1=cluster1,plot2=cluster2)
  return(outList)
  }

  removeOutlierC2 <- function(mergedData, threshold1, threshold2, save=TRUE, colour="sex", base="")
    {
    remove <- mergedData
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)

    cluster1 <- ggplot(data=remove, aes_string(x="C1", y="C2", colour=colour)) +geom_point(size=I(3)) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +scale_colour_manual(values=c("#377EB8","#E41A1C","#4DAF4A"))
    if(save)
      {
      if(base=="") {ggsave("cluster_s0.pdf")} else {ggsave(paste("cluster_s",base,"_s0.pdf",sep=""))}
      }

    if(length(remove[which(remove$C2<threshold1),]$IID)!=0)
      {
      out <- as.character(remove[which(remove$C2<threshold1),]$IID)
      remove <- remove[-which(remove$C2<threshold1),]
      remove$C1 <- scale(remove$C1)
      remove$C2 <- scale(remove$C2)
      }
    else
      {
      out <- ""
      }

    while(length(which(remove$C2<threshold1))!=0)
      {
      out <- append(out, as.character(remove[which(remove$C2<threshold1),]$IID))
      remove <- remove[-which(remove$C2<threshold1),]
      remove$C1 <- scale(remove$C1)
      remove$C2 <- scale(remove$C2)
      }

    if(length(remove[[1]])==0)
      {
      remove <- mergedData
      remove$C1 <- scale(remove$C1)
      remove$C2 <- scale(remove$C2)
      }

    if(length(remove[which(remove$C2>threshold2),]$IID)!=0)
      {
      out <- append(out, as.character(remove[which(remove$C2>threshold2),]$IID))
      remove <- remove[-which(remove$C2>threshold2),]
      remove$C1 <- scale(remove$C1)
      remove$C2 <- scale(remove$C2)
      }

    while(length(which(remove$C2>threshold2))!=0)
      {
      out <- append(out, as.character(remove[which(remove$C2>threshold2),]$IID))
      remove <- remove[-which(remove$C2>threshold2),]
      remove$C1 <- scale(remove$C1)
      remove$C2 <- scale(remove$C2)
      }

    remove <- mergedData[!(mergedData$IID %in% out),]
    remove$C1 <- scale(remove$C1)
    remove$C2 <- scale(remove$C2)

    cluster2 <- ggplot(data=remove, aes_string(x="C1", y="C2", colour=colour)) +geom_point(size=I(3)) +theme_bw() +theme(panel.grid.major=element_line(colour="grey60")) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +scale_colour_manual(values=c("#377EB8","#E41A1C","#4DAF4A"))
    if(save)
      {
      if(base=="") {ggsave(paste("cluster_s",threshold,".pdf",sep=""))} else {ggsave(paste("cluster_s",base,"_C2_s",threshold,".pdf",sep=""))}
      }

    outData <- mergedData[which(mergedData$IID %in% out),]
    outList <- list(data=outData,plot1=cluster1,plot2=cluster2)
    return(outList)
    }
