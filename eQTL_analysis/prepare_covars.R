args = commandArgs(TRUE)

library(plyr)
library(tidyr)


pcaFile = args[1]
metaFile = args[2]
outputfile = args[3]

pca = read.table(pcaFile,h=F)
meta = read.csv(metaFile)
tumIDs = read.csv("/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs.csv"); tumIDs = tumIDs[is.na(tumIDs$GenID) == F,]
tumIDs2 = read.csv("/xx/xxx/SC/CAM_TUM/TUM_genetic_IDs_2023.csv"); tumIDs2 = tumIDs2[is.na(tumIDs2$GenID) == F,]
camIDs = read.table("/xx/xxx/SC/CAM_TUM/sc_for_tum/GT_ID.txt",h=T,sep="\t")

out = pca[,1:11]
colnames(out) = c("IID","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")

meta[meta$iid == "8A2H3PL1","iid"] = "PatID_5"
meta[meta$iid == "E9LEH7P8","iid"] = "PatID_56"
meta[meta$iid == "JEGK54J2","iid"] = "PatID_59"
meta[meta$iid == "RL4X6288","iid"] = "PatID_3"
meta[meta$iid == "WAF0FQN5","iid"] = "PatID_52"
meta[meta$iid == "WQ0EJVR1","iid"] = "PatID_60"

tumIDs = separate(tumIDs,col="GenID",sep="_TUM_TUM",into=c("GenID2","rem"),rem=F)
tumIDs$GenID = tumIDs$GenID2; tumIDs$GenID2 = NULL; tumIDs$rem = NULL
tumIDs = tumIDs[tumIDs$GenID %in% pca$V1 == T,]
tumIDs2 = tumIDs2[tumIDs2$GenID_2023 %in% pca$V1 == T,]


meta = unique(meta[,c("iid","Age","Sex")])
colnames(meta)[1] = "PatID"; meta = join(meta,tumIDs[,c("PatID","GenID")],by="PatID",type="left")
meta = join(meta,tumIDs2[,c("PatID","GenID_2023")],by="PatID",type="left")
colnames(meta)[1] = "ShortID"; meta = join(meta,camIDs[,c("ShortID","GenotypingID")],by="ShortID",type="left")
meta[is.na(meta$GenID) == T,"GenID"] = meta[is.na(meta$GenID) == T,"GenID_2023"] 
meta[is.na(meta$GenID) == T,"GenID"] = meta[is.na(meta$GenID) == T,"GenotypingID"] 
meta = unique(meta[is.na(meta$GenID) == F & meta$GenID != "",c("GenID","Age","Sex")])
meta$IID = meta$GenID



out = join(out,meta[,c("IID","Age","Sex")],by="IID",type="left")



write.table(out,outputfile,sep="\t",row.names=F,quote=F)