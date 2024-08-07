args = commandArgs(TRUE)

library(ggplot2)


dat1 = read.table(args[1])
dat2 = read.table(args[2])
dat3 = read.table(args[3])

inf = read.table(args[4])
outfile = args[5]
outpng1 = args[6]
outpng = args[7]


colnames(dat1)[5] = "CAM_freq"
colnames(dat2)[5] = "TUM_freq"
colnames(dat3)[5] = "t1000G_freq"
dat = merge(dat1,dat2[,c(2,5)],by="V2")
library(plyr)
dat = join(dat,dat3[,c(2,5)],by="V2")

dat$CAM_freq = as.numeric(dat$CAM_freq)
dat$TUM_freq = as.numeric(dat$TUM_freq)
dat$t1000G_freq = as.numeric(dat$t1000G_freq)
dat= dat[is.na(dat$V2) == F,]

library(tidyr)
inf = tidyr::separate(inf,col="V2",sep=";",into=c("ID1","ID2","ID3"),rem=F)
inf$ID2 = gsub("chr","",inf$ID2)
inf$ID3 = gsub("chr","",inf$ID3)

# Find variants to keep
inf$ID = inf$ID2
inf[inf$ID2 %in% dat1$V2 == F,"ID"] = inf[inf$ID2 %in% dat1$V2 == F,"ID3"]
colnames(dat3)[2] = "ID"
inf = join(inf,dat3[,c(2,5)],by="ID",type="left") 

# Plots before QC
# CAM vs. TUM
corr = cor.test(dat$CAM_freq,dat$TUM_freq)$estimate
corr = round(corr,3)
p1 = ggplot(dat,aes(x=CAM_freq,y=TUM_freq)) + geom_point(color="dodgerblue4",alpha=0.6)+
	theme_minimal()+
	xlab("MAF in CAM dataset") + ylab("MAF in TUM dataset")+
	ggtitle(paste0("corr. coeff. = ",corr))

corr = cor.test(inf$V5,inf$t1000G_freq)$estimate
corr = round(corr,3)
p2 = ggplot(inf[is.na(inf$t1000G_freq) == F,],aes(x=V5,y=t1000G_freq)) + geom_point(color="dodgerblue4",alpha=0.6)+
	theme_minimal()+
	xlab("MAF in CAM/TUM dataset") + ylab("MAF in 1000G")+
	ggtitle(paste0("corr. coeff. = ",corr))

library(gridExtra)
png(outpng1,width=15,height=15,units="in",res=300)
grid.arrange(p1,p2,ncol=2)
dev.off()


# Remove rare variants
dat = dat[dat$TUM_freq >= 0.01 & dat$TUM_freq <= 0.99 & dat$CAM_freq >= 0.01 & dat$CAM_freq <= 0.99,]
dat = dat[dat$t1000G_freq >= 0.01 & dat$t1000G_freq <= 0.99 ,]

corr = cor.test(dat$CAM_freq,dat$TUM_freq)$estimate
print(corr)

# Remove variants with AF differing significantly from reference AF (1000G EUR)
removedat = dat[(abs(dat$CAM_freq - dat$t1000G_freq) > 0.2 | abs(dat$TUM_freq - dat$t1000G_freq) > 0.2 | is.na(dat$t1000G_freq) == T) & is.na(dat$V2) == F,]
dat = dat[dat$V2 %in% removedat$V2 == F,]

# Remove ambiguous SNPs with a MAF between 0.4 and 0.6
amb = dat[grepl(":C:G",dat$V2) == T | grepl(":G:C",dat$V2) == T | grepl(":T:A",dat$V2) == T | grepl(":A:T",dat$V2) == T,]
amb = amb[(amb$CAM_freq > 0.4 & amb$CAM_freq < 0.6) | (amb$TUM_freq > 0.4 & amb$TUM_freq < 0.6),]
dat = dat[dat$V2 %in% amb$V2 == F,]
dat = unique(dat)

# Find variants to keep
inf = inf[inf$ID %in% dat$V2 == T,]
write.table(inf$V2,outfile,sep="\t",row.names=F,quote=F)

# Plots MAF vs MAF
# CAM vs. TUM
corr = cor.test(dat$CAM_freq,dat$TUM_freq)$estimate
corr = round(corr,3)
p1 = ggplot(dat,aes(x=CAM_freq,y=TUM_freq)) + geom_point(color="dodgerblue4",alpha=0.6)+
	theme_minimal()+
	xlab("MAF in CAM dataset") + ylab("MAF in TUM dataset")+
	ggtitle(paste0("corr. coeff. = ",corr))

corr = cor.test(inf$V5,inf$t1000G_freq)$estimate
corr = round(corr,3)
p2 = ggplot(inf[is.na(inf$t1000G_freq) == F,],aes(x=V5,y=t1000G_freq)) + geom_point(color="dodgerblue4",alpha=0.6)+
	theme_minimal()+
	xlab("MAF in CAM/TUM dataset") + ylab("MAF in 1000G")+
	ggtitle(paste0("corr. coeff. = ",corr))

library(gridExtra)
png(outpng,width=15,height=15,units="in",res=300)
grid.arrange(p1,p2,ncol=2)
dev.off()
