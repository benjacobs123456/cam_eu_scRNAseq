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

#######################################
# Read in data
#######################################
# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/bcr/")

# Read in data
b_cells = readRDS("../b_cells_post_processing.rds")
rownames(b_cells@meta.data) = colnames(b_cells)

#filter to just celltypist B cell annotations
b_cells = SetIdent(b_cells,value="ann_celltypist_highres")
b_cells = subset(b_cells,cohort=="Cam")

# create file with new IDs for plink
codex = read_table("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/data/genotypes/GT_ID.txt",col_names=FALSE)

# sort out long cambridge IDs
new_ids = lapply(codex$X2, function(x){
  if(grepl("C00",x)){
      new_id = str_remove(pattern="C00TU0",x)
      new_id = str_split(pattern="a",new_id)[[1]][1]
      new_id = paste0("TU",new_id)
      new_id = str_remove(pattern="v1",new_id)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else if(grepl("v",x)){
      new_id = str_remove(pattern="v1",x)
      new_id = str_remove(pattern="v2",new_id)
      return(new_id)
  } else {
      return(x)
  }
})

codex$donor.id = unlist(new_ids)
# read in genelist
genelist = read_tsv("/rds/user/hpcjaco1/hpc-work/genes")
genelist = genelist %>%
filter(name2 %in% rownames(b_cells)) %>%
distinct(name2,.keep_all=T) %>%
mutate(chrom = str_remove(chrom,"chr"))

# initialise results list
expression_res_overall = list()
for(source_to_test in c("CSF","PBMC")){
message(source_to_test)
sc_dat = subset(b_cells,source==source_to_test)
# loop through cell types
results_all_cell_types = list()
for(cell_type_to_test in c("Plasma cells")){
message(cell_type_to_test)
# filter to cell type
rownames(sc_dat@meta.data) = colnames(sc_dat)
sc_dat = subset(sc_dat,subset = ann_celltypist_highres == cell_type_to_test)

# find abundant genes
abundant_genes = rownames(sc_dat)[rowSums(sc_dat@assays$RNA)>1000]
genelist = genelist %>% filter(name2 %in% abundant_genes) %>%
filter(!grepl("^MT",name2)) %>%
filter(!grepl("^RP",name2))

# loop through all genes
overall_results = list()
for(i in c(1:nrow(genelist))){
row = i
message("doing row ",i)
genelist_dat = genelist[row,]
gene_name = genelist_dat$name2
chr = genelist_dat$chrom
start = genelist_dat$txStart
end = genelist_dat$txEnd

dat = sc_dat@assays$RNA[grepl(gene_name,rownames(sc_dat@assays$RNA)),] %>% t %>% data.frame()
gene_name = str_replace(pattern="-",replace=".",string=gene_name)
dat = dat %>% mutate(
cell_name = rownames(dat)) %>%
tidyr::separate(cell_name,sep="_",into=c("source","donor","cellid")) %>%
group_by(donor) %>%
summarise(mean = mean(.data[[gene_name]],na.rm=T)) %>%
mutate(z = (mean - mean(mean) )/ sd(mean) ) %>%
dplyr::select(-mean)

dat = dat %>%
dplyr::rename("donor.id"=donor) %>%
left_join(codex %>% distinct(donor.id,.keep_all=T),by="donor.id") %>%
dplyr::select(X1,z) %>%
dplyr::rename("IID" = X1,"expression_z" = z)

pheno_file =  paste0("/rds/user/hpcjaco1/hpc-work/gene_",gene_name,"_source_",source_to_test,"expression.tsv")
write_tsv(dat,file = pheno_file)

chr_start = start - 1000000
chr_end = end + 1000000

# run plink from r
system(paste0("cd /rds/user/hpcjaco1/hpc-work/;module load plink;plink2 --pfile filtered_scRNAseq_genotypes_plink --pheno ",pheno_file," --glm hide-covar --chr ",chr," --from-bp ",chr_start," --to-bp ",chr_end," --out ",source_to_test,"_",gene_name,"_expression"))

# read in results
expression_res = read_tsv(paste0("/rds/user/hpcjaco1/hpc-work/",source_to_test,"_",gene_name,"_expression.expression_z.glm.linear"),col_types = cols_only(`#CHROM` = col_character(),POS = col_double(),ID=col_character(),REF=col_character(),ALT=col_character(),BETA=col_double(),SE=col_double(),P=col_double()))

# make summary table
snps_tested = expression_res %>% nrow
n_sig_1e3 = expression_res %>% filter(P< 1e-3) %>% nrow
n_sig_1e5 = expression_res %>% filter(P< 1e-5) %>% nrow
sig_snps = expression_res %>% filter(P< 1e-3)
sig_snps = list(sig_snps$ID)

# output complete table with additional info
expression_res = expression_res %>%
mutate(cell_type = cell_type_to_test) %>%
mutate(source = source_to_test) %>%
mutate(gene = gene_name)

expression_res_overall[[length(expression_res_overall)+1]] = expression_res
overall_results[[i]] = list(gene_name,chr,start,end,snps_tested,n_sig_1e3,n_sig_1e5,sig_snps)
}

res = do.call("rbind",overall_results) %>% data.frame()
colnames(res) = c("gene","chr","start","end","snps_tested","n_sig_1e3","n_sig_1e5","sig_snps")
res$cell_type = cell_type_to_test
results_all_cell_types[[length(results_all_cell_types)+1]] = res
}
}

res_df = do.call("rbind",results_all_cell_types) %>% data.frame()
write_tsv(res_df,"eqtls.tsv")
expression_res_overall_df = do.call("bind_rows",expression_res_overall)
write_tsv(expression_res_overall_df,"eqtls_full_results.tsv")

# read in MS GWAS
ms_chip = read_table("/rds/user/hpcjaco1/hpc-work/discovery_metav3.0.meta")
ms_sig_snps = ms_chip %>% filter(P<5e-8)

expression_res_overall_df = expression_res_overall_df %>%
mutate(ID = str_remove(ID,"GSA-")) %>%
mutate(in_ms_chip = ifelse(ID %in% ms_chip$SNP,"present","absent")) %>%
mutate(ms_sig = ifelse(ID %in% ms_sig_snps$SNP,"sig","nonsig")) %>%
left_join(ms_chip %>% dplyr::select(SNP,A1,A2,P,OR) %>%
dplyr::rename("ID" = SNP, "P_MS" = P, "OR_MS" = OR),
by="ID")

chrom_eqtl_plot_dat = expression_res_overall_df %>%
mutate(sig_eqtl = ifelse(P<1e-5,"sig","nonsig")) %>%
dplyr::group_by(`#CHROM`) %>%
dplyr::count(sig_eqtl) %>%
mutate(prop = n/sum(n)) %>%
filter(!is.na(sig_eqtl)) %>%
filter(sig_eqtl=="sig") %>%
mutate(chr = as.numeric(`#CHROM`)) %>%
filter(!is.na(chr)) %>%
filter(chr!=23)

ggplot(chrom_eqtl_plot_dat,
aes(chr,prop))+geom_col()+
scale_x_continuous(breaks=c(1:22))+
labs(y="Proportion of tested SNPs\nwhich are eQTLs at P<1e-5")

qqman::qq(expression_res_overall_df[expression_res_overall_df$`#CHROM`!=6,][[P]])

gene_to_plot="HLA.C"
plot_dat = expression_res_overall_df %>% filter(gene==gene_to_plot)
gene_name_lookup = str_replace(pattern="\\.",replace="-",gene_to_plot)
gene_start = genelist[genelist$name2==gene_name_lookup,][["txStart"]]
gene_end = genelist[genelist$name2==gene_name_lookup,][["txEnd"]]

ggplot(plot_dat,
aes(POS,-log10(P),col=ms_sig))+
geom_point()+
theme_minimal()+
ggtitle(gene_to_plot)+
geom_segment(mapping = aes(x=gene_start,xend = gene_end,y=-2,yend=-2),color="orange")+
annotate("text",x=(gene_start+gene_end)/2,y=-3,label=gene_to_plot,color="orange")+
scale_x_continuous(limits=c(gene_start-100000,gene_end+100000))+
scale_y_continuous(limits=c(-3.1,10),breaks = seq(0,10,by=1))+
geom_hline(yintercept=-log10(5e-8),color="blue",alpha=0.2)+
geom_hline(yintercept=-log10(1e-5),color="blue",alpha=0.2)+
geom_hline(yintercept=-log10(1e-3),color="blue",alpha=0.2)+
scale_color_manual(values = c("grey","red"))+
geom_text_repel(plot_dat %>% filter(P<1e-5 & ms_sig=="sig"),
mapping=aes(POS,-log10(P),label=ID),
max.overlaps=Inf)+
facet_wrap(~source)

snp = "rs2524094"

system(paste0("cd /rds/user/hpcjaco1/hpc-work/;module load plink;plink2 --pfile filtered_scRNAseq_genotypes_plink --snp ",snp," --make-pgen --recode AD --out ",snp,"_genotypes"))

# read genos in
geno_file = paste0("/rds/user/hpcjaco1/hpc-work/",snp,"_genotypes.raw")
genos = read_tsv(geno_file) %>%
left_join(codex %>% dplyr::rename("IID"=X1),by="IID") %>%
dplyr::select(10,7,8)


gene_name = "HLA-C"
csf_pcs = subset(b_cells, source=="CSF" & ann_celltypist_highres=="Plasma cells")
dat = csf_pcs@assays$RNA[grepl(gene_name,rownames(csf_pcs@assays$RNA)),] %>% t %>% data.frame()
gene_name = str_replace(pattern="-",replace=".",string=gene_name)
dat = dat %>% mutate(
cell_name = rownames(dat)) %>%
tidyr::separate(cell_name,sep="_",into=c("source","donor","cellid")) %>%
group_by(donor) %>%
summarise(mean = mean(.data[[gene_name]],na.rm=T)) %>%
mutate(z = (mean - mean(mean) )/ sd(mean) ) %>%
dplyr::select(-mean)

# join with main sc data
dat = dat %>%
dplyr::rename("donor.id"=donor) %>%
filter(donor.id %in% genos$donor.id) %>%
left_join(genos,by="donor.id")


````unix
cd /rds/user/hpcjaco1/hpc-work
module load plink
plink2 --vcf /rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/SAWC_TPS_all_samples_22_10_2020.vcf.gz \
--make-pgen \
--out scRNAseq_genotypes_plink

cd /rds/user/hpcjaco1/hpc-work
module load plink
plink2 --pfile scRNAseq_genotypes_plink \
--maf 0.05 --hwe 1e-5 --geno 0.1 --mind 0.1 \
--out filtered_scRNAseq_genotypes_plink \
--make-pgen

# run association
plink2 --pfile scRNAseq_genotypes_plink --pheno gene_CD86expression.tsv --glm hide-covar --chr 1-4 --maf 0.05 --hwe 1e-5 --geno 0.1 --mind 0.1
````
