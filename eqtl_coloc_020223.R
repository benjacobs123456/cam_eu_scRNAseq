#######################################
# Load packages
#######################################

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(MASS)
library(reshape2)
library(RNOmni)
library(coloc)

#######################################
# Read in data
#######################################


# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/")

# read in eQTLs 
eqtl_res = readRDS("full_eqtl_results_all.rds")

# read in MS GWAS
ms_chip = read_table("~/discovery_metav3.0.meta")

# liftover to hg38
ms_gwas_hg19 = ms_chip %>%
  mutate(chr = paste0("chr",CHR)) %>%
  mutate(start = BP-1) %>%
  dplyr::select(chr,start,BP,SNP)
write_tsv(ms_gwas_hg19,"~/ms_gwas_hg19",col_names=F)
system("~/liftover/liftOver ~/ms_gwas_hg19 ~/liftover/hg19ToHg38.over.chain.gz ~/ms_gwas_hg38.bed ~/unmapped.bed")
ms_gwas_hg38 = read_table("~/ms_gwas_hg38.bed",col_names=F)

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


# coloc
overall_res = list()
cell_types = unique(eqtl_res$cell_type)

#for(pheno in c("MS","OIND","NIND")){
for(pheno in c("MS")){
  for(source_to_test in c("PBMC","CSF")){
    for(cell_type_to_test in cell_types){

  
  message(source_to_test)
  message(cell_type_to_test)
  message(pheno)
  
  # subset
  eqtl_coloc_dat = eqtl_res_coloc %>%
    filter(source == source_to_test & cell_type == cell_type_to_test &
    phenotype==pheno)
  
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
  
      res = coloc_abf$summary %>%
        t %>%
        data.frame %>%
        mutate(gene = this_gene, cell = cell_type_to_test, source = source_to_test, phenotype = pheno)
      overall_res[[length(overall_res)+1]] = res
  
      if(coloc_abf$summary["PP.H4.abf"]>0.7){
  
        snp_res = coloc_abf$results %>%
          arrange(desc(SNP.PP.H4)) %>%
          data.frame() %>%
          mutate(gene = this_gene,
          cell = cell_type_to_test,
          source = source_to_test,
          phenotype = pheno)
        res_title = paste0("coloc_",res$gene,"_",res$cell,"_",res$source,"_",pheno,".csv")
        write_csv(snp_res,res_title)
  
        p=ggplot(coloc_abf$results,aes(position,SNP.PP.H4))+
          geom_point()+
          labs(y="Posterior probability of shared causal variant",x="Genomic position")+
          theme_minimal()+
          geom_text_repel(data = coloc_abf$results %>% filter(SNP.PP.H4>0.1),mapping = aes(position,SNP.PP.H4,label = snp))+
          ggtitle(paste0(res$gene,"\n",res$cell,"\n",res$source,"\nMS"))
        plot_title = paste0("coloc_",res$gene,"_",res$cell,"_",res$source,"_",pheno,".png")
        png(plot_title,res=300,units="in",width=4,height=4)
        print(p)
        dev.off()

      }
    }
  }
  }}

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
overall_res %>% arrange(desc(PP.H4.abf))

# examine in detail
coloc_res = read_csv("coloc_results.csv")
coloc_res$gene %>% unique %>% length
coloc_res %>%
  filter(PP.H4.abf > PP.H3.abf & PP.H4.abf > PP.H2.abf & PP.H4.abf > PP.H1.abf & PP.H4.abf > PP.H0.abf) %>%
  distinct(gene)

sig_coloc_res = coloc_res %>%
  filter(PP.H4.abf > PP.H3.abf & PP.H4.abf > PP.H2.abf & PP.H4.abf > PP.H1.abf & PP.H4.abf > PP.H0.abf)


# look at overlap with eQTLgen
eqtlgen = read_table("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

# liftOver
eqtlgen_hg19 = eqtlgen %>%
  mutate(chr = paste0("chr",SNPChr)) %>%
  mutate(start = SNPPos-1) %>%
  dplyr::select(chr,start,SNPPos,SNP)
write_tsv(eqtlgen_hg19,"~/eqtlgen_hg19",col_names=F)
system("~/liftover/liftOver ~/eqtlgen_hg19 ~/liftover/hg19ToHg38.over.chain.gz ~/eqtlgen_hg38.bed ~/unmapped.bed")

# merge with hg38 coordinates
eqtlgen_hg38 = read_table("~/eqtlgen_hg38.bed",col_names=F)

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
  filter(fdr<0.05) %>%
  dplyr::count(source,cell_type,phenotype,sig_eqtl_gen) %>%
  group_by(source,cell_type,phenotype) %>%
  mutate(total = sum(n)) %>%
  mutate(prop = n / total * 100) %>%
  filter(sig_eqtl_gen == "eQTL")

p=ggplot(plot_dat,aes(cell_type,prop,fill=source,label=total))+
  geom_col(color="black",position=position_dodge(width=1))+
  geom_text(position=position_dodge(width=1),hjust=1)+
  scale_fill_brewer(palette="Set3")+
  coord_flip()+
  facet_wrap(~phenotype)+
  theme_minimal()+
  labs(y="% of eQTLs (FDR<0.05) which are eQTLs in eQTLgen",x="Cell type")

png("eqtlgen_overlap_wholecohort.png",width=8,height=4,res=300,units="in")
p
dev.off()

# look at 'novel' effects
novel_eqtls = eqtls_filtered %>%
  filter(fdr<0.05 & sig_eqtl_gen == "not eQTL")

eqtls_filtered %>%
  filter(fdr<0.05) %>%
  group_by(source) %>%
  dplyr::count(sig_eqtl_gen) %>%
  mutate(prop = n/sum(n))


write_csv(novel_eqtls,"whole_cohort_novel_eqtls.csv")

# print
sig_ms = ms_gwas_hg38 %>% filter(P<5e-8)

eqtl_sig_ms_sig = eqtl_res %>%
  mutate(sig_eqtl_gen = ifelse(chr_pos %in% eqtlgen_hg38_merged$chr_pos,"eQTL","not eQTL")) %>%
  mutate(sig_MS_GWAS = ifelse(chr_pos %in% sig_ms$chr_pos,"MS GWAS hit","Not MS GWAS hit")) %>%
  filter(fdr<0.1)

write_csv(eqtl_sig_ms_sig,"whole_cohort_eqtls_ms_sig.csv")
eqtl_sig_ms_sig_nonmhc = eqtl_sig_ms_sig %>%
  filter(! (CHR==6 & POS >25000000 & POS < 35000000)  )
write_csv(eqtl_sig_ms_sig_nonmhc,"wholecohort_eqtls_ms_sig_nonmhc.csv")

# csf hits
csf_hits = eqtl_sig_ms_sig %>% filter(source=="CSF")
pbmc_hits = eqtl_sig_ms_sig %>% filter(source=="PBMC")

csf_hits %>% distinct(chr_pos) %>% nrow
csf_hits %>% distinct(gene) %>% nrow
write_csv(csf_hits,"whole_cohort_csf_eqtls_ms_sig.csv")
csf_hits %>% dplyr::distinct(chr_pos,.keep_all=TRUE) %>% dplyr::count(sig_eqtl_gen) %>% mutate(prop = n / sum(n))
csf_hits %>% filter(fdr<0.05) %>% dplyr::distinct(chr_pos,.keep_all=TRUE) %>% dplyr::count(sig_eqtl_gen) %>% mutate(prop = n / sum(n))
csf_hits %>% dplyr::distinct(chr_pos,.keep_all=TRUE) %>%
  mutate(sig_pbmc = ifelse(chr_pos %in% pbmc_hits$chr_pos,"eQTL","not eQTL")) %>%
  dplyr::count(sig_pbmc) %>%
  mutate(prop = n / sum(n))


