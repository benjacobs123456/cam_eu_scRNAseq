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
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/")

# read in eQTLs
eqtl_res = readRDS("full_eqtl_results_all.rds")

#######################################
# eQTL analysis big picture
#######################################
# add in FDR
eqtl_res = eqtl_res %>%
  mutate(fdr = p.adjust(P,method="fdr")) %>%
  arrange(fdr)

nrow(eqtl_res)
nrow(eqtl_res %>% distinct(ID,gene))
nrow(eqtl_res %>% distinct(ID))
nrow(eqtl_res %>% distinct(gene))
nrow(eqtl_res %>% distinct(cell_type,gene))
nrow(eqtl_res %>% filter(fdr<0.05) %>% distinct(ID,gene))
nrow(eqtl_res %>% filter(fdr<0.05) %>% distinct(gene))

# plot sig eGenes for cell type and source
eqtl_plot_dat = eqtl_res %>%
  filter(!is.na(P)) %>%
  filter(fdr<0.05) %>%
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
  labs(y="Number of eGenes (FDR<5%)",x="Cell type")+
  scale_fill_brewer(palette="Set3")

eqtl_plot_dat2 = eqtl_res %>%
  filter(!is.na(P)) %>%
  filter(fdr<0.05) %>%
  dplyr::group_by(cell_type,source,phenotype) %>%
  dplyr::count()

p2=ggplot(eqtl_plot_dat2,
          aes(cell_type,n,fill=phenotype))+
  geom_col(position=position_dodge(),color="black")+
  theme_minimal()+
  coord_flip()+
  facet_wrap(~source)+
  labs(y="Number of eQTLs (FDR<5%)",x="Cell type")+
  scale_fill_brewer(palette="Set3")

png("sig_egenes_whole_cohort.png",width=7,height=5,res=300,units="in")
grid.arrange(p,p2,nrow=2)
dev.off()


#######################################
# QQ plots
#######################################
params = expand.grid(unique(eqtl_res$cell_type),
c("CSF","PBMC"),
c("MS","NIND","OIND")
)

library(doParallel)
registerDoParallel(cl = 4)

# draw qqplots
qq_data = foreach(i=1:nrow(params)) %dopar% {
      message(i)
      # get paramaters
      cell_type_to_test = params$Var1[i]
      phenotype_to_test = params$Var3[i]
      source_to_test = params$Var2[i]

      dat = eqtl_res %>%
        filter(
          phenotype==phenotype_to_test,
          source==source_to_test,
          cell_type == cell_type_to_test)

      obs_pval = -log10(dat$P)
      sorted_obs_pval = sort(obs_pval,decreasing=F)
      theoretical_quantiles = seq_along(sorted_obs_pval)/length(sorted_obs_pval)
      exp_pvals = -log10(pchisq(
        qchisq(theoretical_quantiles,df=1,lower.tail=F),
        df=1,lower.tail=F))
      obs_chisq_stats = qchisq(1-dat$P,df=1)
      lambda = median(obs_chisq_stats)/ qchisq(0.5,df=1)
      df = data.frame(obs_pval,exp_pvals,source_to_test,phenotype_to_test,cell_type_to_test,lambda)
      df
}
qq_data_df = do.call("bind_rows",qq_data)

# downsample for plot
qq_data_df_sampled = qq_data_df %>%
  group_by(cell_type_to_test,phenotype_to_test,source_to_test) %>%
  sample_n(size=50000)

# plot
p=ggplot(qq_data_df_sampled,
  aes(exp_pvals,obs_pval,col=cell_type_to_test,shape=phenotype_to_test))+
  geom_point()+
  facet_wrap(~source_to_test)+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_minimal()+
  labs(x="Expected -log10(P)",y="Observed -log10(P)",color="Cell type",shape="Phenotype")+
  scale_color_brewer(palette="Paired")

# combine plots
png("qqplots.png",width=8,height=5,res=600,units="in")
p
dev.off()

# plot
p2=ggplot(qq_data_df_sampled %>%
 distinct(cell_type_to_test,phenotype_to_test,source_to_test,.keep_all=T),
  aes(phenotype_to_test,lambda,fill=cell_type_to_test))+
  geom_col(color="black",position=position_dodge())+
  facet_wrap(~source_to_test)+
  theme_minimal()+
  labs(x="Phenotyoe",y="Lambda",fill="Cell type")+
  scale_fill_brewer(palette="Paired")+
  scale_y_continuous(limits=c(0,1.5))+
  geom_hline(yintercept=1,alpha=0.5,linetype="dashed",color="red")

# combine plots
png("qqplots_lambdas.png",width=6,height=4,res=600,units="in")
p2
dev.off()

#######################################
# correlations
#######################################
# MS v controls
flat_file_cor = eqtl_res %>%
  dplyr::select(ID,BETA,SE,source,gene,cell_type,phenotype,fdr) %>%
  tidyr::pivot_wider(id_cols = c(ID,cell_type, gene),
                     values_from = c(BETA,SE,fdr),
                     names_from = c(phenotype,source))


# just get eqtls which are significant in any test
sig_flat_file_cor = flat_file_cor %>%
  filter_at(.vars = vars(contains("fdr")),
            any_vars(. < 0.05))

# cor test
cor.test(sig_flat_file_cor$BETA_MS_CSF,sig_flat_file_cor$BETA_MS_PBMC)
cor.test(sig_flat_file_cor$BETA_MS_PBMC,sig_flat_file_cor$BETA_NIND_PBMC)
cor.test(sig_flat_file_cor$BETA_MS_CSF,sig_flat_file_cor$BETA_NIND_CSF)
cor.test(sig_flat_file_cor$BETA_MS_CSF,sig_flat_file_cor$BETA_OIND_CSF)

# make plots


make_comparison_plot = function(comparator_1,comparator_2){
  beta_col1 = paste0("BETA_",comparator_1)
  beta_col2 = paste0("BETA_",comparator_2)
  se_col1 = paste0("SE_",comparator_1)
  se_col2 = paste0("SE_",comparator_2)
  xlab = paste0("Beta (",comparator_1,")")
  ylab = paste0("Beta (",comparator_2,")")

# get cor
  cor = cor.test(sig_flat_file_cor[[beta_col1]],sig_flat_file_cor[[beta_col2]])
  r = round(as.numeric(cor$estimate),2)

  p = ggplot(sig_flat_file_cor,
         aes(.data[[beta_col1]],.data[[beta_col2]],col=cell_type))+
    geom_errorbarh(mapping = aes(xmin = .data[[beta_col1]] - .data[[se_col1]] * 1.96,
                                xmax = .data[[beta_col1]] + .data[[se_col1]] * 1.96,
                                y = .data[[beta_col2]]), height=0.1)+
    geom_errorbar(mapping = aes(ymin = .data[[beta_col2]] - .data[[se_col2]] * 1.96,
                                 ymax = .data[[beta_col2]] + .data[[se_col2]] * 1.96,
                                 x = .data[[beta_col1]]), width=0.1)+
    theme_minimal()+
    scale_color_brewer(palette="Paired")+
    geom_abline(intercept = 0,slope=1,alpha=0.5,color="red")+
    scale_x_continuous(limits = c(-4,4))+
    scale_y_continuous(limits = c(-4,4))+
    labs(x=xlab,y=ylab,color="Cell type")+
    ggtitle(paste0("R = ",r))
}

plot_fx = function(x,filename,plotwidth=4,plotheight=4){
  png(filename = filename,res=600,height=plotheight,width=plotwidth,units="in")
  print(x)
  dev.off()
}

p = make_comparison_plot("MS_PBMC","MS_CSF")
p1 = make_comparison_plot("NIND_PBMC","NIND_CSF")
p2 = make_comparison_plot("OIND_PBMC","OIND_CSF")
plot_fx(grid.arrange(p,p1,p2,nrow=1),"csf_v_pbmc.png",plotwidth=12)


p = make_comparison_plot("MS_PBMC","NIND_PBMC")
p1 = make_comparison_plot("MS_PBMC","OIND_PBMC")
p2 = make_comparison_plot("MS_CSF","NIND_CSF")
p3 = make_comparison_plot("MS_CSF","OIND_CSF")
plot_fx(grid.arrange(p,p1,p2,p3,nrow=2),"ms_v_controls.png",plotwidth=10,plotheight=10)

# repeat for cell types
# MS v controls
flat_file_cor_cells = eqtl_res %>%
  dplyr::select(ID,BETA,SE,source,gene,cell_type,phenotype,fdr) %>%
  tidyr::pivot_wider(id_cols = c(ID, gene),
                     values_from = c(BETA,SE,fdr),
                     names_from = c(cell_type,phenotype,source))


# just get eqtls which are significant in any test
sig_flat_file_cor_cells = flat_file_cor_cells %>%
  filter_at(.vars = vars(contains("fdr")),
            any_vars(. < 0.05))

betas_cells = sig_flat_file_cor_cells %>%
  dplyr::select(contains("BETA") & contains("MS"))

# sort out col names
colnames(betas_cells) = str_remove(colnames(betas_cells),"BETA_")

cormat = cor(betas_cells ,use="complete.obs")
png("corplot_ms.png",res=600,units="in",width=8,height=8)
corrplot::corrplot(cormat,order = "hclust",tl.cex = 1,tl.col="black")
dev.off()


#########################################
# MS-specific
#########################################

# define function
find_specific_eqtls = function(cell_type_to_test,source_to_test){

  # filter dat
  dat_for_specific_eqtls = eqtl_res %>%
    filter(cell_type==cell_type_to_test & source==source_to_test)

  # cast wider & filter to MS sig eqtls
  # filter NAs
  dat_for_specific_eqtls = dat_for_specific_eqtls %>%
    pivot_wider(id_cols = c("X.CHROM","POS","ID","A1","REF","ALT","gene","source"),
                values_from = c("BETA","SE","P","fdr"),
                names_from = phenotype) %>%
    filter(fdr_MS < 0.05)

    # calculate het stat
    het_eqtl = dat_for_specific_eqtls %>%
        mutate(het_z_nind = (BETA_MS - BETA_NIND) / (sqrt ( SE_MS^2 + SE_NIND^2) )  ) %>%
        mutate(het_p_nind = 1 - pnorm(abs(het_z_nind))) %>%
        mutate(het_z_oind = (BETA_MS - BETA_OIND) / (sqrt ( SE_MS^2 + SE_OIND^2) )  ) %>%
        mutate(het_p_oind = 1 - pnorm(abs(het_z_oind))) %>%
        mutate(cell_type = cell_type_to_test)

    # add to overall res list
    overall_res <<- bind_rows(overall_res,het_eqtl)
  }

params = expand.grid(
  unique(eqtl_res$cell_type),
  c("PBMC","CSF")
  )
overall_res = data.frame()

for(i in c(1:nrow(params))){
  message(i)
  find_specific_eqtls(params$Var1[i],params$Var2[i])
}

# add FDR
overall_res = overall_res %>%
mutate(fdr_het_nind = p.adjust(het_p_nind,method="fdr")) %>%
mutate(fdr_het_oind = p.adjust(het_p_oind,method="fdr"))

# plots
p1=ggplot(overall_res,aes(BETA_MS,BETA_NIND,col=-log10(het_p_nind)))+
  geom_point()+
  theme_minimal()+
  scale_color_viridis_c(option = "magma")+
  facet_wrap(~cell_type)+
  geom_abline(slope=1,intercept=0,linetype="dashed")
png("het_plots_nind.png",res=600,width=6,height=6,units="in")
p1
dev.off()

p2=ggplot(overall_res,aes(BETA_MS,BETA_OIND,col=-log10(het_p_oind)))+
  geom_point()+
  theme_minimal()+
  scale_color_viridis_c(option = "magma")+
  facet_wrap(~cell_type)+
  geom_abline(slope=1,intercept=0,linetype="dashed")
png("het_plots_oind.png",res=600,width=6,height=6,units="in")
p2
dev.off()

# het plot
make_disease_specific_plot = function(snp_to_test,cell_type,source_to_test,gene_to_test){
cell_type_no_space = str_replace_all(cell_type," ","_")

# filter
plot_dat = eqtl_res %>%
  filter(ID == snp_to_test & gene == gene_to_test & source==source_to_test & cell_type == cell_type)

# run PLINK
system(paste0("module load plink;plink2 --bfile /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc_updated --snp ",snp_to_test," --recode AD --out plot_geno"))
plot_genos = read_tsv("plot_geno.raw")

# get raw counts
pheno_file1 =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/gene_",
gene_to_test,
"_source_",
source_to_test,
"_cell_type_",
cell_type_no_space,
"_pheno_MS_expression.tsv")
pheno_file2 =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/gene_",
gene_to_test,
"_source_",
source_to_test,
"_cell_type_",
cell_type_no_space,
"_pheno_NIND_expression.tsv")
pheno_file3 =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/gene_",
gene_to_test,
"_source_",
source_to_test,
"_cell_type_",
cell_type_no_space,
"_pheno_OIND_expression.tsv")

pheno1 = read_tsv(pheno_file1) %>% mutate(pheno = "MS")
pheno2 = read_tsv(pheno_file2)%>% mutate(pheno = "NIND")
pheno3 = read_tsv(pheno_file3)%>% mutate(pheno = "OIND")
pheno = bind_rows(pheno1,pheno2,pheno3)


plot_dat = plot_genos %>% left_join(pheno,by="IID")
n_samples = plot_dat %>% filter(!is.na(expression_z)) %>% nrow
plot_dat = plot_dat %>% filter(!is.na(expression_z))

colname = colnames(plot_dat)[grepl(snp_to_test,colnames(plot_dat))][1]

p=ggplot(plot_dat,aes(factor(.data[[colname]]),expression_z,fill=pheno))+
  geom_boxplot(position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5))+
  theme_minimal()+
  scale_fill_brewer(palette="Set1")+
  labs(y=paste0("Normalised ",gene_to_test," expression"))+
  ggtitle(paste0(cell_type_no_space,"\n",source_to_test,"\nN=",n_samples,"\nSNP: ",snp_to_test))
p
}



snp_to_test = "6:31281621:T:C"
cell_type = "Plasma cells"
source_to_test = "CSF"
gene_to_test = "HLA.C"

make_disease_specific_plot()
het_snps_sig = overall_res %>%
  filter( (fdr_het_nind < 0.1 & fdr_het_oind < 0.1))


# find CSF-specific
csf = eqtl_res %>%
  tidyr::pivot_wider(id_cols = c(ID,gene,phenotype),
              values_from = c(BETA,SE,P,fdr),
              names_from = c(cell_type,source))

# find CSF hits in B cells
csf$csf_b_cell_het_p = calculate_het_p(csf$`BETA_B cells_CSF`,csf$`BETA_B cells_PBMC`,csf$`SE_B cells_CSF`,csf$`SE_B cells_PBMC`)
csf$csf_plasma_cell_het_p = calculate_het_p(csf$`BETA_Plasma cells_CSF`,csf$`BETA_Plasma cells_PBMC`,csf$`SE_Plasma cells_CSF`,csf$`SE_Plasma cells_PBMC`)

csf %>%
  filter(`fdr_B cells_CSF`<0.05) %>%
  filter(csf_b_cell_het_p < 1e-5) %>%
  dplyr::select(ID,gene,phenotype,contains("B cells"),csf_b_cell_het_p) %>%
  filter(phenotype=="MS") %>%
  print(n=100) %>%
  arrange(`fdr_B cells_CSF`)

csf %>%
  filter(`fdr_Plasma cells_CSF`<0.05) %>%
  filter(csf_plasma_cell_het_p < 1e-3) %>%
  dplyr::select(ID,gene,phenotype,contains("Plasma cells"),csf_plasma_cell_het_p) %>%
  filter(phenotype=="MS") %>%
  print(n=100) %>%glimpse

