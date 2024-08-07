library(data.table)
library(plyr)
library(ggplot2)
library(tidyr)
library(gdata)
library(RNOmni)
library(gridExtra)
library(dplyr)
library(ggrepel)


eQTL_plot = function(snp="rs4676756",gene="EAF2",ctout1 = "B.cells",ctout2 = "CD4.T.cells",ctout3="CD14.Mono",cohort="MS",source="CSF",round=1,window=500000,cell_types = c("CD4.T.cells","CD8.T.cells")){

	if(source == "CSF"){
		source2 = "PBMC"
		gexsource1 = gexcsf
		gexsource2 = gexpbmc
		plotcolor1 = "dodgerblue4"
        plotcolor2 = "firebrick"
        ressource1 = datcsf[datcsf$GENE == gene,]
        ressource2 = datpbmc[datpbmc$GENE == gene,]
         
	}
	if(source == "PBMC"){
		source2 = "CSF"
		gexsource1 = gexpbmc
		gexsource2 = gexcsf
		plotcolor1 = "firebrick"
        plotcolor2 = "dodgerblue4"
        ressource1 = datpbmc[datpbmc$GENE == gene,]
        ressource2 = datcsf[datcsf$GENE == gene,]

	}
	pheno = metakeep

	ct1=gsub("\\."," ",ctout1)
	ct2=gsub("\\."," ",ctout2)
	ct3=gsub("\\."," ",ctout3)

	gene_ens = genlist[genlist$hgnc_symbol == gene,"ensembl_gene_id"]
	chr = genlist[genlist$hgnc_symbol == gene,"chromosome_name"]

	# read genetic data
	start = genlist[genlist$hgnc_symbol == gene,"start_position"] - window
	if(start < 0){start = 0}
	end = genlist[genlist$hgnc_symbol == gene,"end_position"] + window
	chrNumber = genlist[genlist$hgnc_symbol == gene,"chromosome_name"]
	temp_file=paste0("test_",gene,"_",source,"_",cohort,"_")
	system(paste0("/xx/xxx/bin/plink --bfile ",genotype_filename, " --chr ",chrNumber, " --from-bp ",start, " --to-bp ", end, " --recode oxford --out ",temp_file))
	dos = read.table(paste0(temp_file,".gen"))
	dossam = read.table(paste0(temp_file,".sample")); dossam = dossam[-1:-2,]
	snpname = dos[grepl(paste0(snp,";"),dos[,2]) == T,2]
	temp_file2 = paste0(temp_file,"r2")
	system(paste0("/xx/xxx/bin/plink --bfile ",genotype_filename, " --chr ",chrNumber, " --from-bp ",start, " --to-bp ", end, " --r2 --ld-snp '",snpname,"' --ld-window 99999 --ld-window-r2 0 "," --out ",temp_file2)) # for LD
	ld = read.table(paste0(temp_file2,".ld"),h=T)
	file.remove(paste0(temp_file,".gen"))
	file.remove(paste0(temp_file,".sample"))
	file.remove(paste0(temp_file,".nosex"))
	file.remove(paste0(temp_file,".log"))
	file.remove(paste0(temp_file2,".log"))
	file.remove(paste0(temp_file2,".ld"))
	file.remove(paste0(temp_file2,".nosex"))
	dosn = data.frame(dos[,1:5])
	ind = 1
	c = 6
	while(c <= ncol(dos)){
	  individual = dossam[ind,"V2"]
	  for (r in 1:nrow(dos)){
	  if (dos[r,c] == 1){
	    dosn[r,individual] = 2
	  }
	  if (dos[r,c+1] == 1){
	    dosn[r,individual] = 1
	  }
	  if (dos[r,c+2] == 1){
	    dosn[r,individual] = 0
	  }
	}
	  ind = ind + 1 
	  c = c + 3
	}
	dos = dosn[grepl(paste0(snp,";"),dosn$V2) == T,]
	dos = data.frame(t(dos))
	minallele = dos[4,1]
	majallele = dos[5,1]
	dos = dos[-1:-5,]
	dos = data.frame(snp=dos)
	rownames(dos) = colnames(dosn)[6:ncol(dosn)]
	dos$IID = rownames(dos)
	dos = dos[dos$IID %in% pheno$IID == T,]


	# Prepare gene experession data
	gexsource1 = gexsource1[gene,]
	gex1 = gexsource1[,startsWith(colnames(gexsource1),paste0(ctout1,"_")) == T]
	gex1 = data.frame(t(gex1))
	gex1$sam = rownames(gex1)
	gex1$sam = gsub(paste0(ctout1,"_"),"",gex1$sam)
	gex1$sam = gsub("OIND_","",gex1$sam)
	gex1$sam = gsub("ONIND_","",gex1$sam)
	gex1$sam = gsub("MS_","",gex1$sam)
	gex1$sam = gsub("OINDI_","",gex1$sam)
	gex1$Sample = gex1$sam; gex1$sam = NULL
	gex1 = join(gex1,pheno[,c("IID","Sample")],by="Sample",type="left")
	gex1 = gex1[gex1$IID %in% rownames(dos) == T,]
	gex1[,gene]= RankNorm(gex1[,gene])
	gex2 = gexsource1[,startsWith(colnames(gexsource1),paste0(ctout2,"_")) == T]
	gex2 = data.frame(t(gex2))
	gex2$sam = rownames(gex2)
	gex2$sam = gsub(paste0(ctout2,"_"),"",gex2$sam)
	gex2$sam = gsub("OIND_","",gex2$sam)
	gex2$sam = gsub("ONIND_","",gex2$sam)
	gex2$sam = gsub("MS_","",gex2$sam)
	gex2$Sample = gex2$sam; gex2$sam = NULL
	gex2 = join(gex2,pheno[,c("IID","Sample")],by="Sample",type="left")
	gex2 = gex2[gex2$IID %in% rownames(dos) == T,]
	gex2[,gene]= RankNorm(gex2[,gene])
	gex3 = gexsource1[,startsWith(colnames(gexsource1),paste0(ctout3,"_")) == T]
	gex3 = data.frame(t(gex3))
	gex3$sam = rownames(gex3)
	gex3$sam = gsub(paste0(ctout3,"_"),"",gex3$sam)
	gex3$sam = gsub("OIND_","",gex3$sam)
	gex3$sam = gsub("ONIND_","",gex3$sam)
	gex3$sam = gsub("MS_","",gex3$sam)
	gex3$Sample = gex3$sam; gex3$sam = NULL
	gex3 = join(gex3,pheno[,c("IID","Sample")],by="Sample",type="left")
	gex3 = gex3[gex3$IID %in% rownames(dos) == T,]
	gex3[,gene]= RankNorm(gex3[,gene])

	dat1 = join(gex1,dos,by="IID",type="left")
	dat1$gene = dat1[,gene]
	dat2 = join(gex2,dos,by="IID",type="left")
	dat2$gene = dat2[,gene]
	dat3 = join(gex3,dos,by="IID",type="left")
	dat3$gene = dat3[,gene]


	dat1 = dat1[is.na(dat1$snp) == F & dat1$snp != "NA",]
	dat2 = dat2[is.na(dat2$snp) == F & dat2$snp != "NA",]
	dat3 = dat3[is.na(dat3$snp) == F & dat3$snp != "NA",]

    genefoundinsource2 = F
	if (gene %in% rownames(gexsource2) == T){
		genefoundinsource2 = T
	}else{genefoundinsource2 == F}
	if (genefoundinsource2 == T){
	gexsource2 = gexsource2[gene,]
	gexs21 = gexsource2[,startsWith(colnames(gexsource2),paste0(ctout1,"_")) == T]
	gexs21 = data.frame(t(gexs21))
	gexs21$sam = rownames(gexs21)
	gexs21$sam = gsub(paste0(ctout1,"_"),"",gexs21$sam)
	gexs21$sam = gsub("OIND_","",gexs21$sam)
	gexs21$sam = gsub("ONIND_","",gexs21$sam)
	gexs21$sam = gsub("MS_","",gexs21$sam)
	gexs21$Sample = gexs21$sam; gexs21$sam = NULL
	gexs21 = join(gexs21,pheno[,c("IID","Sample")],by="Sample",type="left")
	gexs21 = gexs21[gexs21$IID %in% rownames(dos) == T,]
	gexs21[,gene]= RankNorm(gexs21[,gene])
	gexs22 = gexsource2[,startsWith(colnames(gexsource2),paste0(ctout2,"_")) == T]
	gexs22 = data.frame(t(gexs22))
	gexs22$sam = rownames(gexs22)
	gexs22$sam = gsub(paste0(ctout2,"_"),"",gexs22$sam)
	gexs22$sam = gsub("OIND_","",gexs22$sam)
	gexs22$sam = gsub("ONIND_","",gexs22$sam)
	gexs22$sam = gsub("MS_","",gexs22$sam)
	gexs22$Sample = gexs22$sam; gexs22$sam = NULL
	gexs22 = join(gexs22,pheno[,c("IID","Sample")],by="Sample",type="left")
	gexs22 = gexs22[gexs22$IID %in% rownames(dos) == T,]
	gexs22[,gene]= RankNorm(gexs22[,gene])
	gexs23 = gexsource2[,startsWith(colnames(gexsource2),paste0(ctout3,"_")) == T]
	gexs23 = data.frame(t(gexs23))
	gexs23$sam = rownames(gexs23)
	gexs23$sam = gsub(paste0(ctout3,"_"),"",gexs23$sam)
	gexs23$sam = gsub("OIND_","",gexs23$sam)
	gexs23$sam = gsub("ONIND_","",gexs23$sam)
	gexs23$sam = gsub("MS_","",gexs23$sam)
	gexs23$Sample = gexs23$sam; gexs23$sam = NULL
	gexs23 = join(gexs23,pheno[,c("IID","Sample")],by="Sample",type="left")
	gexs23 = gexs23[gexs23$IID %in% rownames(dos) == T,]
	gexs23[,gene]= RankNorm(gexs23[,gene])

	dats21 = join(gexs21,dos,by="IID",type="left")
	dats21$gene = dats21[,gene]
	dats22 = join(gexs22,dos,by="IID",type="left")
	dats22$gene = dats22[,gene]
	dats23 = join(gexs23,dos,by="IID",type="left")
	dats23$gene = dats23[,gene]

	dats21 = dats21[is.na(dats21$snp) == F & dats21$snp != "NA",]
	dats22 = dats22[is.na(dats22$snp) == F & dats22$snp != "NA",]
	dats23 = dats23[is.na(dats23$snp) == F & dats23$snp != "NA",]

}	

	

	


	# Prepare eQTL results of the locus for cohort and source
	res = ressource1[ressource1$cohort == cohort,]
	res = res[res$GENE == gene,]
	res = data.frame(res)
	res$POS = as.numeric(as.character(res$POS))
	res$P = as.numeric(as.character(res$P))
	pos0 = res[grepl(paste0(snp,";"),res$ID) == T,"POS"][1]
	pos0 = as.numeric(as.character(pos0))
	pleft = pos0-200000
	if (pleft < 0){pleft = 0}
	pright = pos0+200000
	res = res[res$POS >= pleft & res$POS <= pright,]
	res$p = -log10(res$P)
	res$BP = res$POS
	res$POS = res$POS/1000000
	res1 = res[res$cell_type == ctout1,]
	res2 = res[res$cell_type == ctout2,]
	res3 = res[res$cell_type == ctout3,]

	

	res = ressource2[ressource2$cohort == cohort,]
	res = res[res$GENE == gene,]
	res = data.frame(res)
	print(nrow(res))
	if (nrow(res) > 0){
	res$POS = as.numeric(as.character(res$POS))
	res$P = as.numeric(as.character(res$P))
	pos0 = res[grepl(paste0(snp,";"),res$ID) == T,"POS"][1]
	pos0 = as.numeric(as.character(pos0))
	pleft = pos0-200000
	if (pleft < 0){pleft = 0}
	pright = pos0+200000
	res = res[res$POS >= pleft & res$POS <= pright,]
	res$p = -log10(res$P)
	res$POS = res$POS/1000000
	ress21 = res[res$cell_type == ctout1,]
	ress22 = res[res$cell_type == ctout2,]
	ress23 = res[res$cell_type == ctout3,]
    }

    if(nrow(res) > 0){
    	maxy = max(res1$p,res2$p,res3$p,ress21$p,ress22$p,ress23$p)+0.2
    }else{
    	maxy = max(res1$p,res2$p,res3$p)+0.2
    }
    # Association plots
 	ld$ID = ld$SNP_B
    res1 = join(res1,ld[,c("ID","R2")],by="ID",type="left")
    res1 = res1[is.na(res1$R2) == F,]
    res1[res1$SNP == snp,"R2"] = 1.2
    res1$R2d = NA; res1[res1$R2 <= 0.2,"R2d"] = 0.2
    res1[res1$R2 > 0.2 & res1$R2 <=0.4,"R2d"] = 0.4
    res1[res1$R2 > 0.4 & res1$R2 <= 0.6,"R2d"] = 0.6
    res1[res1$R2 > 0.6 & res1$R2 <= 0.8,"R2d"] = 0.8
    res1[res1$R2 > 0.8 & res1$R2 <= 1,"R2d"] = 1.0
    res1[res1$SNP == snp,"R2d"] = 1.2
    res1$R2d = factor(res1$R2d,levels=c("0.2","0.4","0.6","0.8","1","1.2"))
    locuszoomcolors=c("#01178B","#87CEEB","#006300","#FFA500","#FF2500","white")
    res1 = res1[order(res1$R2,decreasing=T),]
	p1=ggplot(res1)+
	            geom_point(aes(POS,p,fill=R2d),shape=21,color="black",size=2.2,key_glyph='rect')+ 
				geom_point(size=3,aes(POS[R2==1.2],p[R2==1.2]),shape=23,fill="#9F34F0",color="black")+
	            scale_fill_manual(values=locuszoomcolors,breaks=c("0.2","0.4","0.6","0.8","1"),labels=c(expression(phantom(x)<=0.2),expression(phantom(x)<=0.4),expression(phantom(x)<=0.6),expression(phantom(x)<=0.7),expression(phantom(x)<="1.0")),drop=FALSE,name=bquote(r^2))+
	            geom_text_repel(data=res1 %>% head(n=1),mapping=aes(x=POS,y=p,label=SNP),nudge_y=0.1,size=3)+
	            #scale_color_gradient(low="grey50",high="darkorange",name=bquote(r^2))+
	            labs(y="-log10(eQTL p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
	            theme_minimal()+
	#            ggtitle(paste0("Association plot for ",gene," expression in\n",source," ",ct1))+
	            ylim(0,maxy)+
	            theme(legend.position=c(0.95,0.8),legend.key=element_rect(color="black",size=1.2))+
	            # ZCD2HC1A and ETS1 c(0.15,0.78)
	            # AHI1 EAF2 c(0.95,0.8)
	            # PREX1 c(0.95,0.8)
	            theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
	            theme(axis.line = element_line(color="black",size=0.5),axis.ticks = element_line(color="black",size=0.5))+
	            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9), legend.title=element_text(family="sans",size=9,face="bold"))
	p2=ggplot(res2,aes(POS,p))+
	            geom_point(col=plotcolor1)+
	            labs(y="-log10(eQTL p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
	            theme_minimal()+
	            ggtitle(paste0("Association plot for ",gene," expression in\n",source," ",ct2))+
	            ylim(0,maxy)+
	            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	p3=ggplot(res3,aes(POS,p))+
	            geom_point(col=plotcolor1)+
	            labs(y="-log10(eQTL p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
	            theme_minimal()+
	            ggtitle(paste0("Association plot for ",gene," expression in\n",source," ",ct3))+
	            ylim(0,maxy)+
	            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	if(nrow(res) > 0){
	ps21=ggplot(ress21,aes(POS,p))+
	            geom_point(col=plotcolor2)+
	            labs(y="-log10(eQTL p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
	            theme_minimal()+
	            ggtitle(paste0("Association plot for ",gene," expression in\n",source2," ",ct1))+
	            ylim(0,max(res1$p,res2$p,res3$p,ress21$p,ress22$p,ress23$p)+0.2)+
	            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	ps22=ggplot(ress22,aes(POS,p))+
	            geom_point(col=plotcolor2)+
	            labs(y="-log10(eQTL p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
	            theme_minimal()+
	            ggtitle(paste0("Association plot for ",gene," expression in\n",source2," ",ct2))+
	            ylim(0,max(res1$p,res2$p,res3$p,ress21$p,ress22$p,ress23$p)+0.2)+
	            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	ps23=ggplot(ress23,aes(POS,p))+
	            geom_point(col=plotcolor2)+
	            labs(y="-log10(eQTL p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
	            theme_minimal()+
	            ggtitle(paste0("Association plot for ",gene," expression in\n",source2," ",ct3))+
	            ylim(0,max(res1$p,res2$p,res3$p,ress21$p,ress22$p,ress23$p)+0.2)+
	            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))

    }

	# Boxplots

	p4 <- ggplot(dat1, aes(x=as.factor(snp), y=gene)) + geom_boxplot(fill=plotcolor1)+ geom_jitter(width=0.2) +
      ylab(paste0(gene," expr. in ",ct1)) + xlab(paste0("N. of ",snp,"*",minallele," alleles")) + 
      theme_minimal() + ggtitle(paste0(gene," expression in\n",source," ",ct1)) + 
      theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	p5 <- ggplot(dat2, aes(x=as.factor(snp), y=gene)) + geom_boxplot(fill=plotcolor1) + geom_jitter(width=0.2) +
      ylab(paste0(gene," expr. in ",ct2)) + xlab(paste0("N. of ",snp,"*",minallele," alleles")) + 
      theme_minimal() + ggtitle(paste0(gene," expression in\n",source," ",ct2)) + 
      theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	p6 <- ggplot(dat3, aes(x=as.factor(snp), y=gene)) + geom_boxplot(fill=plotcolor1) + geom_jitter(width=0.2) +
      ylab(paste0(gene," expr. in ",ct3)) + xlab(paste0("N. of ",snp,"*",minallele," alleles")) + 
      theme_minimal() + ggtitle(paste0(gene," expression in\n",source," ",ct3)) + 
      theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
    if (genefoundinsource2 == T){
    ps24 <- ggplot(dats21, aes(x=as.factor(snp), y=gene)) + geom_boxplot(fill=plotcolor2) + geom_jitter(width=0.2) +
      ylab(paste0(gene," expr. in ",ct1)) + xlab(paste0("N. of ",snp,"*",minallele," alleles")) + 
      theme_minimal() + ggtitle(paste0(gene," expression in\n",source2," ",ct1)) + 
      theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	ps25 <- ggplot(dats22, aes(x=as.factor(snp), y=gene)) + geom_boxplot(fill=plotcolor2) + geom_jitter(width=0.2) +
      ylab(paste0(gene," expr. in ",ct2)) + xlab(paste0("N. of ",snp,"*",minallele," alleles")) + 
      theme_minimal() + ggtitle(paste0(gene," expression in\n",source2," ",ct2)) + 
      theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
	ps26 <- ggplot(dats23, aes(x=as.factor(snp), y=gene)) + geom_boxplot(fill=plotcolor2) + geom_jitter(width=0.2) +
      ylab(paste0(gene," expr. in ",ct3)) + xlab(paste0("N. of ",snp,"*",minallele," alleles")) + 
      theme_minimal() + ggtitle(paste0(gene," expression in\n",source2," ",ct3)) + 
      theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))
    }

    # Forest plots

	all_res = list()

	for (cta in  cell_types){

	  print(cta)
      rescsf = datcsf[datcsf$cell_type == cta & datcsf$GENE == gene,]
      rescsf = data.frame(rescsf)
      if (nrow(rescsf) > 0){
      rescsf = rescsf[rescsf$GENE == gene & grepl(paste0(snp,";"),rescsf$ID) == T,]
      rescsf$CHR = rescsf$X.CHROM
      rescsf = rescsf[rescsf$P != "P" & is.na(rescsf$P) == F & rescsf$P != "",]
      rescsf$BETA = as.numeric(as.character(rescsf$BETA))
      rescsf$SE = as.numeric(as.character(rescsf$SE))
	  res_a = c("All","CSF",cta,rescsf[1,"BETA"],rescsf[1,"SE"],rescsf[1,"A1"])
      names(res_a) = c("Cohort","Source","CellType","BETA","SE","A1")
      print(res_a)
      all_res[[length(all_res)+1]] = res_a}

      respbmc = datpbmc[datpbmc$cell_type == cta & datpbmc$GENE == gene,]
      respbmc = data.frame(respbmc)
      if (nrow(respbmc) > 0){
      respbmc = respbmc[respbmc$GENE == gene & grepl(paste0(snp,";"),respbmc$ID) == T,]
      respbmc$CHR = respbmc$X.CHROM
      respbmc = respbmc[respbmc$P != "P" & is.na(respbmc$P) == F & respbmc$P != "",]
      respbmc$BETA = as.numeric(as.character(respbmc$BETA))
      respbmc$SE = as.numeric(as.character(respbmc$SE))
	  res_a = c("All","PBMC",cta,respbmc[1,"BETA"],respbmc[1,"SE"],respbmc[1,"A1"])
      names(res_a) = c("Cohort","Source","CellType","BETA","SE","A1")
      all_res[[length(all_res)+1]] = res_a}

	}

	all_res = do.call("bind_rows",all_res)
	print(all_res)
	all_res = data.frame(all_res)
	all_res$BETA = as.numeric(as.character(all_res$BETA))
	all_res$SE = as.numeric(as.character(all_res$SE))
	all_res = all_res[is.na(all_res$A1) == F,]
	all_res[all_res$A1 != all_res[all_res$Cohort == cohort & all_res$Source == source & all_res$CellType == ctout1,"A1"],"BETA"] = all_res[all_res$A1 != all_res[all_res$Cohort == cohort & all_res$Source == source & all_res$CellType == ctout1,"A1"],"BETA"]*1
	all_res$LOWER = all_res$BETA - 1.96*all_res$SE
	all_res$UPPER = all_res$BETA + 1.96*all_res$SE

	all_res[all_res$CellType == "Plasma.cells","CellType"] = "ASCs"

	#all_res$CellType = factor(all_res$CellType,ordered=T,levels=c("pDCs","mDCs","CD4.T.cells","CD8.T.cells","Tregs","MAIT.cells","NK.cells","B.cells","ASCs","CD14.Mono","CD16.Mono"))
	all_res$CellType = factor(all_res$CellType,ordered=T,levels=c("CD16.Mono","CD14.Mono","ASCs","B.cells","NK.cells","MAIT.cells","Tregs","CD8.T.cells","CD4.T.cells","mDCs","pDCs"))
	all_res$Cohort = factor(all_res$Cohort,ordered=T,levels=c("All"))
	all_res$Source = factor(all_res$Source,ordered=T,levels=c("CSF","PBMC"))



	p7 = ggplot(all_res,aes(factor(CellType),as.numeric(BETA))) + geom_pointrange(aes(ymin=LOWER,ymax=UPPER,color=Source),position = position_dodge(-0.3))+
	scale_color_manual(values=c("dodgerblue4","firebrick"))+
	coord_flip() + theme_minimal() + 
	ylab("eQTL effect size") + xlab("Cell type") + 
	geom_hline(yintercept=0,color="grey20") +
	theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
	theme(legend.margin=margin(0,5,0,0),legend.box.margin=margin(-10,0,-10,-25))+
	#ggtitle(paste0("Effect sizes for ",snp,"\n in different ","cell types")) +
	theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), strip.text.x = element_text(size=9,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))

    # MS Gwas association plot
	ms_gwas_plot1 = ms_gwas_plot[ms_gwas_plot$BP >= pleft & ms_gwas_plot$BP <= pright,]
	ms_gwas_plot1 = ms_gwas_plot1[ms_gwas_plot1$CHR == res1$X.CHROM[1],]
	ms_gwas_plot1 = ms_gwas_plot1[ms_gwas_plot1$POS >= min(res1$POS) & ms_gwas_plot1$POS <= max(res1$POS),]


    p=ggplot(ms_gwas_plot1,aes(POS,P))+
            geom_point(col="black")+
            labs(y="-log10(p)",x=paste0("Chromosome ",chr," position (Mio BP)"))+
            theme_minimal()+
            ggtitle(paste0("Association plot for\n","MS risk"))+
            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"), legend.text=element_text(family="sans",size=9))

            #ggtitle(paste0(res$gene,"\n",res$cell,"\n",res$source,"\n",res$pheno))
    pms=p

    # Locuscompare plot
    comp1 = res1
    comp1$eQTLP = comp1$p
    comp1$P = NULL
    comp1 = join(comp1,ms_gwas_plot1[,c("BP","P")],by="BP",type="left")
    comp1$MS_P = comp1$P; comp1 = comp1[is.na(comp1$MS_P) == F,]

    p=ggplot(comp1)+
            geom_point(aes(x=MS_P,y=eQTLP,fill=R2d),shape=21,color="black",size=2.2,key_glyph='rect')+
            geom_point(size=3,aes(x=MS_P[R2==1.2],y=eQTLP[R2==1.2]),shape=23,color="black",fill="#9F34F0")+
            scale_fill_manual(values=locuszoomcolors,breaks=c("0.2","0.4","0.6","0.8","1"),labels=c(expression(phantom(x)<=0.2),expression(phantom(x)<=0.4),expression(phantom(x)<=0.6),expression(phantom(x)<=0.7),expression(phantom(x)<="1.0")),drop=FALSE,name=bquote(r^2))+
	        geom_text_repel(data=comp1 %>% head(n=1),mapping=aes(x=MS_P,y=eQTLP,label=SNP),nudge_y=0.1,size=3)+
            labs(y="-log10(eQTL p)",x="-log10(MS GWAS p)")+
            ylim(0,maxy)+
            theme_minimal()+
            theme(legend.position="none")+
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
            theme(axis.line = element_line(color="black",size=0.5),axis.ticks = element_line(color="black",size=0.5))+
            theme(plot.title = element_text(family="sans",size=10,face = "bold",hjust=0.5), text = element_text(size=9,family="sans"), legend.text = element_text(family="sans",size=9), axis.text.x = element_text(family="sans",size=9), axis.text.y = element_text(family="sans",size=9), legend.title = element_text(family="sans",size=9,face = "bold"), axis.title.x = element_text(family="sans",size=9,face = "bold"), axis.title.y = element_text(family="sans",size=9,face = "bold"))

            #ggtitle(paste0(res$gene,"\n",res$cell,"\n",res$source,"\n",res$pheno))
    compp=p

	plots = list()
	plots[[length(plots)+1]] = p1
	plots[[length(plots)+1]] = p2
	plots[[length(plots)+1]] = p3
	plots[[length(plots)+1]] = p4
	plots[[length(plots)+1]] = p5
	plots[[length(plots)+1]] = p6
	if(nrow(res) > 0){
	plots[[length(plots)+1]] = ps21
	plots[[length(plots)+1]] = ps22
	plots[[length(plots)+1]] = ps23
	}else{
	plots[[length(plots)+1]] = NA
	plots[[length(plots)+1]] = NA
	plots[[length(plots)+1]] = NA
	}
    if (genefoundinsource2 == T){
	plots[[length(plots)+1]] = ps24
	plots[[length(plots)+1]] = ps25
	plots[[length(plots)+1]] = ps26
	}else{
	plots[[length(plots)+1]] = NA
	plots[[length(plots)+1]] = NA
	plots[[length(plots)+1]] = NA
	}
	plots[[length(plots)+1]] = p7
	plots[[length(plots)+1]] = pms
	plots[[length(plots)+1]] = compp

	return(plots)
}
