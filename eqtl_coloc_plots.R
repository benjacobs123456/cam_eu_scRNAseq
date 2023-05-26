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
library(Matrix.utils)
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

# prepare data
eqtl_res = eqtl_res %>% mutate(chr_pos = paste0(CHR,":",POS))
plot_dat_eqtls = eqtl_res %>% left_join(ms_gwas_hg38,by="chr_pos")

# plots
do_eqtl_plot = function(snp,gene_name,source_to_test,cell_type_no_space,pheno_to_test="MS"){
  cell_type_for_plot = str_replace_all(cell_type_no_space,"_"," ")

  # get MS risk allele
  this_snp = plot_dat_eqtls %>% filter(gene==gene_name & source == source_to_test & cell_type == cell_type_for_plot & chr_pos == snp) %>% distinct(chr_pos,.keep_all=T)
  ea = ifelse(this_snp$OR>1,this_snp$A1.y,this_snp$A2)
  nea = ifelse(this_snp$OR>1,this_snp$A2,this_snp$A1.y)
  message("MS risk allele is ",ea)

  # run PLINK
  system(paste0("module load plink;plink2 --bfile /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc_updated --snp ",this_snp$ID," --recode AD --out plot_geno"))

  plot_genos = read_tsv("plot_geno.raw")

  # reorient to ms risk-allele
  snp_ea = paste0(this_snp$ID,"_",ea)
  snp_nea = paste0(this_snp$ID,"_",nea)

  plot_genos$genotypes = if(is.null(plot_genos[[snp_ea]])){
    2 - plot_genos[[snp_nea]]
  } else {
    plot_genos[[snp_ea]]
  }


  pheno_file =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/gene_",gene_name,"_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_",pheno_to_test,"_expression.tsv")
  pheno = read_tsv(pheno_file)
  plot_dat = plot_genos %>% left_join(pheno,by="IID")
  n_samples = plot_dat %>% filter(!is.na(expression_z)) %>% nrow


  p=ggplot(plot_dat,aes(factor(genotypes),expression_z))+
    geom_boxplot()+
    geom_jitter(width=0.1)+
    theme_minimal()+
    labs(y=paste0("Normalised ",gene_name," expression"),x=paste0("No. of ",snp,"-",ea," alleles"))+
    ggtitle(paste0(cell_type_no_space,"\n",source_to_test,"\nN=",n_samples,"\nSNP: ",this_snp$ID))


  plot_dat = plot_dat_eqtls %>%
    filter(gene==gene_name & source == source_to_test & cell_type == cell_type_for_plot & phenotype==pheno_to_test)

  # read in genelist
  genelist = read_tsv("/rds/user/hpcjaco1/hpc-work/genes.gz") # downloaded from UCSC
  genelist = genelist %>%
    distinct(name2,.keep_all=T) %>%
    mutate(chrom = str_remove(chrom,"chr"))

  gene_name_lookup = str_replace(pattern="\\.",replace="-",gene_name)
  gene_start = genelist[genelist$name2==gene_name_lookup,][["txStart"]]
  gene_end = genelist[genelist$name2==gene_name_lookup,][["txEnd"]]

  max_p = max(-log10(plot_dat$P.x))*1.1
  p1=ggplot(plot_dat,
            aes(POS,-log10(P.x),col=-log10(P.y)))+
    geom_point()+
    theme_minimal()+
    geom_segment(mapping = aes(x=gene_start,xend = gene_end,y=-2,yend=-2),color="orange")+
    annotate("text",x=(gene_start+gene_end)/2,y=-3,label=gene_name,color="orange")+
    scale_x_continuous(limits=c(gene_start-1000000,gene_end+1000000))+
    scale_y_continuous(limits=c(-5.1,max_p),breaks = seq(0,max_p,by=2))+
    geom_hline(yintercept=-log10(1e-5),color="blue",alpha=0.1)+
    scale_colour_gradient(low="purple",high="orange")+
    geom_text_repel(plot_dat %>% filter(P.x<1e-5 & P.y<5e-8),
                    mapping=aes(POS,-log10(P.x),label=ID),
                    max.overlaps=Inf,
                    show.legend=F)+
    labs(color="MS GWAS P (-log10(P))",y = "eQTL P value (-log10(P))",x="Genomic position (hg38)")+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

  plot_dat = plot_dat_eqtls %>%
    filter(gene == gene_name & chr_pos == snp) %>%
    mutate(BETA = ifelse(A1.x !=ea,BETA*-1,BETA))

  p2=ggplot(plot_dat,aes(BETA,cell_type,shape=source,col=phenotype))+
    geom_point(position=ggstance::position_dodgev(height=0.5))+
    geom_errorbarh(mapping=aes(xmin = BETA-1.96*SE,xmax=BETA+1.96*SE),
                   height=0.1,position=ggstance::position_dodgev(height=0.5))+
    theme_minimal()+
    geom_vline(xintercept=0,alpha=0.2)+
    scale_color_brewer(palette="Set1")+
    labs(x="eQTL effect",y="Cell type") +
    ggtitle(paste0(gene_name,"\n",snp))

  outfile = paste0("whole_cohort_",gene_name,"_",source_to_test,"_",cell_type_no_space,"_eqtl.png")
  png(file = outfile,res=300,units="in",width=6,height=8)
  print(grid.arrange(p,p1,p2))
  dev.off()
}


# read in coloc results
sig_coloc = read_csv("coloc_results.csv") %>%
 filter(`PP.H4.abf` >= 0.7)

for(i in c(1:nrow(sig_coloc))){
  message("doing row ",i," of ",nrow(sig_coloc))
  this_gene = sig_coloc$gene[i]
  this_source = sig_coloc$source[i]
  this_cell = sig_coloc$cell[i]
  this_cell_no_space = str_replace_all(this_cell," ","_")
  # read in coloc res
  file = paste0("coloc_",this_gene,"_",this_cell,"_",this_source,"_MS.csv")
  coloc = read_csv(file,col_types = cols(.default = "c"))

  # make eqtl plots

  # run plink from r
  do_eqtl_plot(snp = coloc$snp[1],
               gene_name = this_gene,
               source_to_test = this_source,
               cell_type_no_space = this_cell_no_space)

}

