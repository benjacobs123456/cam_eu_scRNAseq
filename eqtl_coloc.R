#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(MASS)
library(SingleCellExperiment)
library(Matrix.utils)
library(reshape2)
library(RNOmni)
library(coloc)

#######################################
# Read in data
#######################################

# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/eqtl/")

files = list.files(pattern="full_results")

files_in = purrr::map(files,function(x){
  read_tsv(x) %>%
    mutate(CHR = X.CHROM) %>%
    filter(!is.na(P))
})

eqtl_res = do.call("bind_rows",files_in)


# read in MS GWAS
ms_chip = read_table("/rds/user/hpcjaco1/hpc-work/discovery_metav3.0.meta")
ms_gwas_hg38 = read_table("ms_gwas_hg38.bed",col_names=F)

# read in hg38 coordinates
ms_gwas_hg38 = ms_chip %>%
  left_join(ms_gwas_hg38 %>% dplyr::rename(SNP = X4),by="SNP") %>%
  filter(!is.na(X3)) %>%
  filter(!is.na(P) & !is.na(OR)) %>%
  mutate(BP = X3) %>%
  dplyr::select(-X1,-X2,-X3) %>%
  mutate(SNP = ifelse(grepl("rs",SNP),SNP,paste0("chr",CHR,":",BP)))

ms_gwas_hg38 = ms_gwas_hg38 %>% mutate(chr_pos = paste0(CHR,":",BP))

# find overlapping snps
eqtl_res = eqtl_res %>% mutate(chr_pos = paste0(CHR,":",POS))
eqtl_res_coloc = eqtl_res %>% filter(chr_pos %in% ms_gwas_hg38$chr_pos)
ms_gwas_coloc = ms_gwas_hg38 %>% filter(chr_pos %in% eqtl_res_coloc$chr_pos)


#######################################
# eQTL analysis big picture
#######################################

nrow(eqtl_res)
nrow(eqtl_res %>% distinct(ID,gene))
nrow(eqtl_res %>% distinct(ID))
nrow(eqtl_res %>% distinct(gene))
nrow(eqtl_res %>% distinct(cell_type,gene))

# add in FDR
eqtl_res = eqtl_res %>%
  mutate(fdr = p.adjust(P,method="fdr")) %>%
    arrange(fdr)


# look at inflation
plots = list()
for(this_cell_type in unique(eqtl_res$cell_type)){
  for(this_source in unique(eqtl_res$source)){
    for(this_pheno in unique(eqtl_res$phenotype)){

      plot_dat = eqtl_res %>%
        filter(source == this_source & cell_type == this_cell_type & phenotype==this_pheno) %>%
        filter(!is.na(P))

      pvals = plot_dat$P
      chisq_stats = qchisq(1-pvals,df=1)
      length(chisq_stats)
      lambda = round(median(chisq_stats)/qchisq(df=1,p=0.5),2)

      observed = sort(-log10(pvals))
      expected = sort(-log10(ppoints(length(pvals))))
      df = data.frame(observed,expected)
      p=ggplot(df,aes(expected,observed))+
        geom_point()+
          theme_minimal()+
          geom_abline(slope=1,intercept=0,color="red")+
          scale_x_continuous(limits=c(0,6))+
          scale_y_continuous(limits=c(0,6))+
          labs(x="Expected -log10(P)",y="Observed -log10(P)")+
          ggtitle(paste0(this_pheno,"\n",this_cell_type,"\n",this_source,"\nGC=",lambda))
      plots[[length(plots)+1]] = p
    }
  }
}
png("qqplots.png",width=16,height=16,res=300,units="in")
do.call("grid.arrange",plots)
dev.off()

# plot sig eGenes for cell type and source
eqtl_plot_dat = eqtl_res %>%
  filter(!is.na(P)) %>%
  filter(fdr<0.1) %>%
  dplyr::group_by(cell_type,source,gene,phenotype) %>%
  dplyr::count() %>%
  ungroup() %>%
  dplyr::count(cell_type,source,phenotype)

p=ggplot(eqtl_plot_dat,
  aes(cell_type,n,fill=phenotype))+
  geom_col(position=position_dodge(),color="black")+
  theme_minimal()+
  coord_flip()+
  facet_wrap(~source)+
  labs(y="Number of eGenes (FDR<10%)",x="Cell type",fill="Phenotype")+
  scale_fill_brewer(palette="Set3")


eqtl_plot_dat2 = eqtl_res %>%
  filter(!is.na(P)) %>%
  filter(fdr<0.1) %>%
  dplyr::group_by(cell_type,source,phenotype) %>%
  dplyr::count()

p2=ggplot(eqtl_plot_dat2,
  aes(cell_type,n,fill=phenotype))+
  geom_col(position=position_dodge(),color="black")+
  theme_minimal()+
  coord_flip()+
  facet_wrap(~source)+
  labs(y="Number of eQTLs (FDR<10%)",x="Cell type",fill="phenotype")+
  scale_fill_brewer(palette="Set3")

png("sig_egenes.png",width=5,height=5,res=300,units="in")
grid.arrange(p,p2,nrow=2)
dev.off()
png("sig_egenes_main.png",width=5,height=3,res=300,units="in")
p
dev.off()

# power
# add in cell numbers
# Read in data
all_combo = readRDS("../datasets/all_combo_with_updated_pheno.rds")

# add full cell id
all_combo@meta.data = all_combo@meta.data %>%
  mutate(full_cell_id = paste0(source,"_",donor.id,"_",rownames(all_combo@meta.data)))
rownames(all_combo@meta.data) = colnames(all_combo)

# filter to just Cam data
all_combo = subset(all_combo,cohort=="Cam")

# get counts
counts = all_combo@meta.data %>% dplyr::count(cell_type,source,phenotype)

# join
plot_dat3 = eqtl_plot_dat %>% left_join(counts,by=c("cell_type","source","phenotype"))

# plot
p = ggplot(plot_dat3,aes(n.y,n.x, shape = source, col = cell_type))+
  geom_point(size=3)+
  scale_x_log10()+scale_y_log10()+labs(x="Number of cells",y="Number of eGenes")+
  theme_minimal()+
  facet_wrap(~phenotype)
png("eqtl_by_cell_no.png",width=8,height=4,res=300,units="in")
p
dev.off()

p2 = ggplot(plot_dat3 %>% filter(phenotype=="MS"),aes(n.y,n.x, shape = source, col = cell_type))+
  geom_point(size=3)+
  scale_x_log10()+scale_y_log10()+labs(x="Number of cells",y="Number of eGenes")+
  theme_minimal()+
  labs(color = "Cell type",shape = "Source")
png("eqtl_by_cell_no_just_ms.png",width=4,height=4,res=300,units="in")
p2
dev.off()

# more big picture numbers
eqtl_res %>% filter(fdr < 0.1) %>% distinct(ID) %>% nrow()
eqtl_res %>% filter(fdr < 0.1) %>% distinct(gene) %>% nrow()

eqtl_res %>% filter(P < 1e-3) %>% distinct(ID) %>% nrow()
eqtl_res %>% filter(P < 1e-3) %>% distinct(gene) %>% nrow()

p_plot = ggplot(eqtl_res,aes(-log10(P)))+geom_histogram(bins=100,color="black")
png("pvals.png",width=4,height=4,res=300,units="in")
p_plot
dev.off()


# coloc
overall_res = list()
cell_types = unique(eqtl_res$cell_type)

for(source_to_test in c("PBMC","CSF")){
  for(cell_type_to_test in cell_types){
    for(phenotype_to_test in c("MS")){

message(source_to_test)
message(cell_type_to_test)

# subset
eqtl_coloc_dat = eqtl_res_coloc %>%
  filter(source == source_to_test & cell_type == cell_type_to_test & phenotype == phenotype_to_test)
ms_gwas_coloc_dat = ms_gwas_coloc %>% filter(chr_pos %in% eqtl_coloc_dat$chr_pos)

# loop through each locus
for(i in c(1:length(unique(eqtl_coloc_dat$gene)))){
  message("doing ", i," of ",length(unique(eqtl_coloc_dat$gene)))
  this_gene = unique(eqtl_coloc_dat$gene)[i]

  coloc_eqtl_gene = eqtl_coloc_dat %>% filter(gene == this_gene) %>% distinct(chr_pos,.keep_all=T)
  coloc_ms_gene = ms_gwas_coloc_dat %>% filter(chr_pos %in% coloc_eqtl_gene$chr_pos) %>% distinct(chr_pos,.keep_all=T)


  # format
  coloc_ms_gene = coloc_ms_gene %>%
    mutate(beta = log(OR)) %>%
    mutate(z = (qnorm(1-P/2))) %>%
    mutate( se = abs(beta / z)) %>%
    mutate( se = ifelse(z == 0,0,se)) %>%
    mutate(varbeta = se ^ 2) %>%
    mutate(varbeta = ifelse(varbeta == 0,min(coloc_ms_gene$varbeta[coloc_ms_gene$varbeta!=0]),varbeta))

  coloc_ms_gene_input = list(
    beta = coloc_ms_gene$beta,
    varbeta = coloc_ms_gene$varbeta,
    snp = coloc_ms_gene$chr_pos,
    position = coloc_ms_gene$BP,
    type="cc")

  coloc_eqtl_gene_input = list(
    beta = coloc_eqtl_gene$BETA,
    varbeta = coloc_eqtl_gene$SE^2,
    snp = coloc_eqtl_gene$chr_pos,
    position = coloc_eqtl_gene$POS,
    type="quant",
    sdY = 1)

    coloc_abf = coloc.abf(coloc_ms_gene_input,coloc_eqtl_gene_input)

    res = coloc_abf$summary %>% t %>% data.frame %>% mutate(gene = this_gene, cell = cell_type_to_test, source = source_to_test, pheno = phenotype_to_test)
    overall_res[[length(overall_res)+1]] = res

    if(coloc_abf$summary["PP.H4.abf"]>0.7){
      p=ggplot(coloc_abf$results,aes(position,SNP.PP.H4))+
        geom_point()+
        labs(y="Posterior probability of shared causal variant",x="Genomic position")+
        theme_minimal()+
        geom_text_repel(data = coloc_abf$results %>% filter(SNP.PP.H4>0.5),mapping = aes(position,SNP.PP.H4,label = snp))+
        ggtitle(paste0(res$gene,"\n",res$cell,"\n",res$source,"\n",res$pheno))
      plot_title = paste0("coloc_",res$gene,"_",res$cell,"_",res$source,"_",res$pheno,".png")
      png(plot_title,res=300,units="in",width=4,height=4)
      print(p)
      dev.off()

    }
}
}}}

overall_res = do.call("bind_rows",overall_res)
write_csv(overall_res,"coloc_results.csv")

overall_res_pbmc = overall_res %>% filter(source=="PBMC")
overall_res_csf = overall_res %>% filter(source=="CSF")
overall_res_pbmc %>% distinct(gene) %>% nrow
overall_res_pbmc %>% filter(PP.H4.abf>0.7) %>% distinct(gene) %>% nrow
overall_res_pbmc %>% filter(PP.H4.abf>PP.H3.abf & PP.H4.abf>PP.H2.abf & PP.H4.abf>PP.H1.abf & PP.H4.abf>PP.H0.abf) %>% distinct(gene) %>% nrow
overall_res_csf %>% distinct(gene) %>% nrow
overall_res_csf %>% filter(PP.H4.abf>0.7) %>% distinct(gene) %>% nrow
overall_res %>% filter(PP.H4.abf>PP.H3.abf & PP.H4.abf>PP.H2.abf & PP.H4.abf>PP.H1.abf & PP.H4.abf>PP.H0.abf) %>% distinct(gene) %>% nrow

# examine in detail
coloc_res = read_csv("coloc_results.csv")
coloc_res$gene %>% unique %>% length
coloc_res %>%
  filter(PP.H4.abf > PP.H3.abf & PP.H4.abf > PP.H2.abf & PP.H4.abf > PP.H1.abf & PP.H4.abf > PP.H0.abf) %>%
  distinct(gene)

# look at overlap with eQTLgen
eqtlgen = read_table("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

# merge with hg38 coordinates
eqtlgen_hg38 = read_table("eqtlgen_hg38.bed",col_names=F)

eqtlgen_hg38_merged = eqtlgen %>%
  left_join(eqtlgen_hg38 %>% dplyr::rename("SNP" = X4),by="SNP") %>%
  filter(!is.na(X3)) %>%
  filter(!is.na(FDR)) %>%
  mutate(BP = X3) %>%
  dplyr::select(-X1,-X2,-X3)


eqtlgen_hg38_merged = eqtlgen_hg38_merged %>% mutate(chr_pos = paste0(SNPChr,":",BP))

# filter eqtl snps to those represented in eqtlgen
eqtls_filtered = eqtl_res %>% mutate(sig_eqtl_gen = ifelse(chr_pos %in% eqtlgen_hg38_merged$chr_pos,"eQTL","not eQTL"))

plot_dat = eqtls_filtered %>%
  filter(phenotype=="MS") %>%
  filter(fdr<0.05) %>%
  dplyr::count(source,cell_type,sig_eqtl_gen) %>%
  group_by(source,cell_type) %>%
  mutate(total = sum(n)) %>%
  mutate(prop = n / total * 100) %>%
  filter(sig_eqtl_gen == "eQTL")

p=ggplot(plot_dat,aes(cell_type,prop,fill=source,label=total))+
  geom_col(color="black",position=position_dodge(width=1))+
  geom_text(position=position_dodge(width=1),hjust=1)+
  scale_fill_brewer(palette="Set3")+
  coord_flip()+
  theme_minimal()+
  labs(y="% of eQTLs (FDR<0.05) which are eQTLs in eQTLgen",x="Cell type")

png("eqtlgen_overlap.png",width=8,height=4,res=300,units="in")
p
dev.off()

# look at 'novel' effects
novel_eqtls = eqtls_filtered %>%
  filter(fdr<0.1 & sig_eqtl_gen == "not eQTL")

write_csv(novel_eqtls,"novel_eqtls.csv")

# disease-specific eqtls
discordant_qtls_res = data.frame()
plots = list()
# correlation between betas
for(cell_type_to_test in unique(eqtl_res$cell_type)){
message(cell_type_to_test)
ms_nind = eqtl_res %>%
    filter(source == "PBMC" & cell_type == cell_type_to_test) %>%
    filter(phenotype %in% c("MS","Noninflammatory")) %>%
    tidyr::pivot_wider(id_cols = c(CHR,POS,ID,REF,ALT,A1,gene),values_from = c(BETA,SE,P,fdr),names_from = phenotype) %>%
    filter(!is.na(BETA_MS) & !is.na(BETA_Noninflammatory)) %>%
  mutate(Z_MS = BETA_MS/SE_MS) %>%
  mutate(Z_Noninflammatory = BETA_Noninflammatory/SE_Noninflammatory) %>%
  mutate(Z_diff = (BETA_MS - BETA_Noninflammatory) / (sqrt(SE_MS^2 + SE_Noninflammatory^2))) %>%
  mutate(P_diff = 2*(1-pnorm(abs(Z_diff)))) %>%
  mutate(P_diff_fdr = p.adjust(P_diff,method="fdr")) %>%
  mutate(cell = cell_type_to_test)

discordant_qtls = ms_nind %>% filter(fdr_MS<0.1 & P_diff_fdr<0.1) %>% mutate(comparison = "MS_NIND")
discordant_qtls_res <<- bind_rows(discordant_qtls_res,discordant_qtls)

ms_nind = ms_nind %>% filter(fdr_MS<0.1)
if(nrow(ms_nind)>2){
  r = cor(ms_nind$BETA_MS,ms_nind$BETA_Noninflammatory)
  message("r = ",r)
  p=ggplot(ms_nind,aes(BETA_MS,BETA_Noninflammatory))+
    geom_point()+
    geom_errorbar(
      mapping = aes(
        ymin = BETA_MS - SE_MS * 1.96,
        ymax = BETA_MS + SE_MS * 1.96,
        x = BETA_Noninflammatory
        ), width=0.1)+
        geom_errorbarh(
          mapping = aes(
            xmin = BETA_Noninflammatory - SE_Noninflammatory * 1.96,
            xmax = BETA_Noninflammatory + SE_Noninflammatory * 1.96,
            y = BETA_MS
            ), height=0.1)+
            theme_minimal()+
            geom_smooth(method="lm")+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.5,2.5))+
            ggtitle(paste0(cell_type_to_test,"\nR = ",round(r,2)))
  plots[[length(plots)+1]] = p
  }
}

png("ms_nind_beta_cor_plot.png",width=8,height=8,res=300,units="in")
do.call("grid.arrange",plots)
dev.off()


# repeat for OIND
# disease-specific eqtls
plots = list()
# correlation between betas
for(cell_type_to_test in unique(eqtl_res$cell_type)){
message(cell_type_to_test)
ms_oind = eqtl_res %>%
    filter(source == "PBMC" & cell_type == cell_type_to_test) %>%
    filter(phenotype %in% c("MS","OIND")) %>%
    tidyr::pivot_wider(id_cols = c(CHR,POS,ID,REF,ALT,A1,gene),values_from = c(BETA,SE,P,fdr),names_from = phenotype) %>%
    filter(!is.na(BETA_MS) & !is.na(BETA_OIND)) %>%
  mutate(Z_MS = BETA_MS/SE_MS) %>%
  mutate(Z_OIND = BETA_OIND/SE_OIND) %>%
  mutate(Z_diff = (BETA_MS - BETA_OIND) / (sqrt(SE_MS^2 + SE_OIND^2))) %>%
  mutate(P_diff = 2*(1-pnorm(abs(Z_diff)))) %>%
  mutate(P_diff_fdr = p.adjust(P_diff,method="fdr")) %>%
  mutate(cell = cell_type_to_test)

discordant_qtls = ms_oind %>% filter(fdr_MS<0.1 & P_diff_fdr<0.1) %>% mutate(comparison = "MS_OIND")
discordant_qtls_res <<- bind_rows(discordant_qtls_res,discordant_qtls)

ms_oind = ms_oind %>% filter(fdr_MS<0.1)
if(nrow(ms_oind)>2){
  r = cor(ms_oind$BETA_MS,ms_oind$BETA_OIND)
  message("r = ",r)
  p=ggplot(ms_oind,aes(BETA_MS,BETA_OIND))+
    geom_point()+
    geom_errorbar(
      mapping = aes(
        ymin = BETA_MS - SE_MS * 1.96,
        ymax = BETA_MS + SE_MS * 1.96,
        x = BETA_OIND
        ), width=0.1)+
        geom_errorbarh(
          mapping = aes(
            xmin = BETA_OIND - SE_OIND * 1.96,
            xmax = BETA_OIND + SE_OIND * 1.96,
            y = BETA_MS
            ), height=0.1)+
            theme_minimal()+
            geom_smooth(method="lm")+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.5,2.5))+
            ggtitle(paste0(cell_type_to_test,"\nR = ",round(r,2)))
  plots[[length(plots)+1]] = p
  }
}

png("ms_oind_beta_cor_plot.png",width=8,height=8,res=300,units="in")
do.call("grid.arrange",plots)
dev.off()

discordant_qtls_res %>%
    filter(sign(BETA_MS) != sign(BETA_OIND))

discordant_qtls_res %>%
    filter(sign(BETA_MS) != sign(BETA_Noninflammatory)) %>%


# print
sig_ms = ms_gwas_hg38 %>% filter(P<5e-8)

eqtl_sig_ms_sig = eqtl_res %>%
  mutate(sig_eqtl_gen = ifelse(chr_pos %in% eqtlgen_hg38_merged$chr_pos,"eQTL","not eQTL")) %>%
  mutate(sig_MS_GWAS = ifelse(chr_pos %in% sig_ms$chr_pos,"MS GWAS hit","Not MS GWAS hit")) %>%
  filter(fdr<0.1 & phenotype=="MS")

write_csv(eqtl_sig_ms_sig,"eqtls_ms_sig.csv")
eqtl_sig_ms_sig_nonmhc = eqtl_sig_ms_sig %>%
  filter(! (CHR==6 & POS >25000000 & POS < 35000000)  )
write_csv(eqtl_sig_ms_sig_nonmhc,"eqtls_ms_sig_nonmhc.csv")

# csf hits
csf_hits = eqtl_sig_ms_sig %>% filter(source=="CSF")
pbmc_hits = eqtl_sig_ms_sig %>% filter(source=="PBMC")

csf_hits %>% distinct(chr_pos) %>% nrow
csf_hits %>% distinct(gene) %>% nrow
write_csv(csf_hits,"csf_eqtls_ms_sig.csv")
csf_hits %>% dplyr::distinct(chr_pos,.keep_all=TRUE) %>% dplyr::count(sig_eqtl_gen) %>% mutate(prop = n / sum(n))
csf_hits %>% filter(fdr<0.05) %>% dplyr::distinct(chr_pos,.keep_all=TRUE) %>% dplyr::count(sig_eqtl_gen) %>% mutate(prop = n / sum(n))
csf_hits %>% dplyr::distinct(chr_pos,.keep_all=TRUE) %>%
  mutate(sig_pbmc = ifelse(chr_pos %in% pbmc_hits$chr_pos,"eQTL","not eQTL")) %>%
  dplyr::count(sig_pbmc) %>%
  mutate(prop = n / sum(n))


# plots
plot_dat_eqtls = eqtl_res %>% left_join(ms_gwas_hg38,by="chr_pos")
do_eqtl_plot = function(snp,gene_name,source_to_test,cell_type_no_space){
cell_type_for_plot = str_replace_all(cell_type_no_space,"_"," ")

# get MS risk allele
this_snp = plot_dat_eqtls %>% filter(gene==gene_name & source == source_to_test & cell_type == cell_type_for_plot & chr_pos == snp) %>% distinct(chr_pos,.keep_all=T)
ea = ifelse(this_snp$OR>1,this_snp$A1.y,this_snp$A2)
nea = ifelse(this_snp$OR>1,this_snp$A2,this_snp$A1.y)

# run PLINK
system(paste0("module load plink;plink2 --pfile /rds/user/hpcjaco1/hpc-work/filtered_scRNAseq_genotypes_plink --snp ",this_snp$ID," --recode AD --out plot_geno"))

plot_genos = read_tsv("plot_geno.raw")

# reorient to ms risk-allele
snp_ea = paste0(this_snp$ID,"_",ea)
snp_nea = paste0(this_snp$ID,"_",nea)

plot_genos$genotypes = if(is.null(plot_genos[[snp_ea]])){
  2 - plot_genos[[snp_nea]]
} else {
  plot_genos[[snp_ea]]
}


pheno_file =  paste0("/rds/user/hpcjaco1/hpc-work/gene_",gene_name,"_source_",source_to_test,"_cell_type_",cell_type_no_space,"_expression.tsv")
pheno = read_tsv(pheno_file)
plot_dat = plot_genos %>% left_join(pheno,by="IID")
n_samples = plot_dat %>% filter(!is.na(expression_z)) %>% nrow


p=ggplot(plot_dat,aes(factor(genotypes),expression_z))+
geom_boxplot()+
geom_jitter(width=0.1)+
theme_minimal()+
labs(y=paste0("Normalised ",gene_name," expression"),x=paste0("No. of ",snp,"-",ea," alleles"))+
ggtitle(paste0(cell_type_no_space,"\n",source_to_test,"\nN=",n_samples,"\nSNP: ",this_snp$ID))

plot_dat = plot_dat_eqtls %>% filter(gene==gene_name & source == source_to_test & cell_type == cell_type_for_plot)

# read in genelist
genelist = read_tsv("/rds/user/hpcjaco1/hpc-work/genes") # downloaded from UCSC
genelist = genelist %>%
  distinct(name2,.keep_all=T) %>%
  mutate(chrom = str_remove(chrom,"chr"))

gene_name_lookup = str_replace(pattern="\\.",replace="-",gene_name)
gene_start = genelist[genelist$name2==gene_name_lookup,][["txStart"]]
gene_end = genelist[genelist$name2==gene_name_lookup,][["txEnd"]]

p1=ggplot(plot_dat,
  aes(POS,-log10(P.x),col=-log10(P.y)))+
  geom_point()+
  theme_minimal()+
  geom_segment(mapping = aes(x=gene_start,xend = gene_end,y=-2,yend=-2),color="orange")+
  annotate("text",x=(gene_start+gene_end)/2,y=-3,label=gene_name,color="orange")+
  scale_x_continuous(limits=c(gene_start-100000,gene_end+100000))+
  scale_y_continuous(limits=c(-3.1,10),breaks = seq(0,10,by=2))+
  geom_hline(yintercept=-log10(1e-5),color="blue",alpha=0.1)+
  scale_colour_gradient(low="purple",high="orange")+
  geom_text_repel(plot_dat %>% filter(P.x<1e-5 & P.y<5e-8),
  mapping=aes(POS,-log10(P.x),label=ID),
  max.overlaps=Inf,
  show.legend=F)+
  labs(color="MS GWAS P (-log10(P))",y = "eQTL P value (-log10(P))",x="Genomic position (hg38)")+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

  plot_dat = plot_dat_eqtls %>%
  filter(gene == gene_name & chr_pos == snp & phenotype=="MS") %>%
    mutate(BETA = ifelse(A1.x !=ea,BETA*-1,BETA))

  p2=ggplot(plot_dat,aes(BETA,cell_type,col=source))+
    geom_point(position=ggstance::position_dodgev(height=0.5))+
    geom_errorbarh(mapping=aes(xmin = BETA-1.96*SE,xmax=BETA+1.96*SE),height=0.1,position=ggstance::position_dodgev(height=0.5))+theme_minimal()+geom_vline(xintercept=0,alpha=0.2)+scale_color_brewer(palette="Set1")+labs(x="eQTL effect",y="Cell type") + ggtitle(paste0(gene_name,"\n",snp))

  outfile = paste0(snp,"_",gene_name,"_",source_to_test,"_",cell_type_no_space,"_eqtl.png")
  png(file = outfile,res=300,units="in",width=6,height=8)
  print(grid.arrange(p,p1,p2))
  dev.off()
}


# run plink from r

do_eqtl_plot(
snp = "1:157709901",
gene_name = "FCRL2",
source_to_test = "PBMC",
cell_type_no_space = "B_cells")


do_eqtl_plot(
snp = "3:121824051",
gene_name = "EAF2",
source_to_test = "PBMC",
cell_type_no_space = "B_cells")

do_eqtl_plot(
snp = "6:31360765",
gene_name = "HLA.B",
source_to_test = "PBMC",
cell_type_no_space = "CD4_T_cells")

do_eqtl_plot(
snp = "6:31360765",
gene_name = "HLA.C",
source_to_test = "PBMC",
cell_type_no_space = "CD4_T_cells")


# pleiotropy
pleio_snps = plot_dat_eqtls %>% filter(fdr<0.1 & phenotype=="MS") %>%
  dplyr::count(chr_pos,gene) %>%
  dplyr::count(chr_pos) %>% filter(n>1)
pleio_snps_dat = plot_dat_eqtls %>% filter(chr_pos %in% pleio_snps$chr_pos & phenotype=="MS")
pleio_snps_dat %>% distinct(chr_pos)
pleio_snps_dat = pleio_snps_dat %>% mutate(in_mhc = ifelse(CHR.x==6 & POS > 25000000 & POS < 35000000,"MHC","Non_MHC"))
pleio_snps_dat %>% distinct(chr_pos,.keep_all=T) %>% dplyr::count(in_mhc)
write_csv(pleio_snps_dat,"pleio_snps.csv")

pleioplot = function(snp_to_test){
this_snp = plot_dat_eqtls %>%
  distinct(chr_pos,.keep_all=T) %>%
  filter(ID == snp_to_test)
ea = ifelse(this_snp$OR>1,this_snp$A1.y,this_snp$A2)
nea = ifelse(this_snp$OR>1,this_snp$A2,this_snp$A1.y)

hla_snp =  plot_dat_eqtls %>%
  filter(ID == snp_to_test) %>%
  filter(phenotype=="MS") %>%
  mutate(BETA = ifelse(A1.x !=ea,BETA*-1,BETA))

p=ggplot(hla_snp,aes(cell_type,BETA,col=gene))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(mapping=aes(x=cell_type,ymin = BETA-1.96*SE, ymax = BETA + 1.96*SE),position=position_dodge(width=0.5),width=0.1) + facet_wrap(~source)+geom_hline(yintercept=0,alpha=0.5)+coord_flip()+labs(y="Beta (eQTL effect size)",x="Cell type")+theme_minimal()+scale_color_brewer(palette="Set1")

outfile = paste0("pleio_plot_",snp,".png")
png(file = outfile,res=300,units="in",width=6,height=8)
print(p)
dev.off()
}

pleioplot("rs3891175")
pleioplot("rs2523578")
