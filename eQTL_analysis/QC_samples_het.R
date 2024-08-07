args <- commandArgs(TRUE)



het <- read.table(args[1], h=T)
colnames(het)[5] = "F"
het["Fscaled"] <- scale(het$F)
remove <- het[which(abs(het$Fscaled)>5),]
write.table(remove[,1],args[2],sep="\t",row.names=F,quote=F,col.names=F)

