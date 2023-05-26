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
library(SingleCellExperiment)
library(Matrix.utils)
library(reshape2)
library(RNOmni)
library(doParallel)
registerDoParallel(cl = detectCores())

#######################################
# Read in data
#######################################
# get args
args = commandArgs(trailingOnly = TRUE)
message(args)

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/")


cell_types = c("B cells",
"Plasma cells",
"mDCs",
"CD14 Mono",
"CD16 Mono",
"CD4 T cells",
"CD8 T cells",
"Tregs",
"MAIT cells",
"NK cells",
"pDCs"
)


# filter to cell type
cell_type_to_test =  cell_types[as.numeric(args[1])]
cell_type_no_space =  str_replace_all(cell_type_to_test," ","_")

# read in file file
infile = paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/eqtl_sc_data_",cell_type_no_space)

message(paste0("Reading in ",infile))
eqtl_dataset = readRDS(infile)

# read in SNPs
snps_for_eqtl = read_tsv("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc_updated.bim",col_names=F)
message("SNPs for eQTL testing")
message(nrow(snps_for_eqtl))

# fam file for eqtl
fam_for_eqtl = read_table("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc_updated.fam",col_names=F)

# read in genelist
genelist = read_tsv("/rds/user/hpcjaco1/hpc-work/genes.gz") # downloaded from UCSC
genelist = genelist %>%
  filter(name2 %in% rownames(eqtl_dataset)) %>%
  distinct(name2,.keep_all=T) %>%
  mutate(chrom = str_remove(chrom,"chr")) %>%
  filter(chrom != "X" & chrom != "Y")

# initialise results list
expression_res_overall = list()

# eqtl analysis
pheno = args[2]
source_to_test = args[3]

# filter sc data
sc_dat = subset(eqtl_dataset,source==source_to_test & Category == pheno)

message(source_to_test)
message(cell_type_to_test)

rownames(sc_dat@meta.data) = colnames(sc_dat)

# find out how many donors in sc data
message(sc_dat@meta.data %>% distinct(iid) %>% nrow," people in sc data")

# restrict to donors with >=2 cells
donors_to_keep = sc_dat@meta.data %>% dplyr::count(genotyping_id) %>% filter(n>=2)
sc_dat = subset(sc_dat,subset = genotyping_id %in% donors_to_keep$genotyping_id)


# find out how many donors in sc data
message(sc_dat@meta.data %>% distinct(iid) %>% nrow," people in sc data after filtering")
message(sc_dat@meta.data %>% filter(genotyping_id %in% fam_for_eqtl$X2) %>% distinct(iid) %>% nrow," people in sc data & genotype data after filtering")

# convert to sce object
sc_dat.sce = as.SingleCellExperiment(sc_dat)

# get counts
counts = assay(sc_dat.sce,"counts")

# transform to log2(cpm)
libsizes = colSums(counts(sc_dat.sce))
logcounts(sc_dat.sce) = log2(t(t(counts)/libsizes*1e6) + 1)

# get rid of low-expressed genes
sc_dat.sce = sc_dat.sce[rowSums(counts(sc_dat.sce))>1000,]

# get rid of mitochondrial and ribosomal genes
sc_dat.sce = sc_dat.sce[!grepl("^MT",rownames(sc_dat.sce)),]
sc_dat.sce = sc_dat.sce[!grepl("^RP",rownames(sc_dat.sce)),]

# aggregate
groups = colData(sc_dat.sce)[,"genotyping_id"]
logcounts_dat = t(logcounts(sc_dat.sce)) %>%
  data.frame()
logcounts_dat$full_cell_id = rownames(logcounts_dat)

# add donor id
donor_ids = sc_dat@meta.data %>% dplyr::select(genotyping_id)
donor_ids$full_cell_id = rownames(donor_ids)


logcounts_dat = logcounts_dat %>%
  left_join(donor_ids,by="full_cell_id") %>%
  group_by(genotyping_id) %>%
  summarise_at(vars(-full_cell_id),.fun = "mean", na.omit = T)

# get pcs
pcs = prcomp(logcounts_dat %>% dplyr::select(-genotyping_id),
  scale = T,
  center = T)

pcs = pcs$x %>% data.frame()
pcs$donor.id = logcounts_dat$genotyping_id

# write to file
# get genotype ID

pcs = pcs %>%
  dplyr::select(donor.id,PC1,PC2,PC3,PC4) %>%
  dplyr::rename("IID" = donor.id) %>%
  filter(!is.na(IID))

# read in genetic PCs
genetic_pcs = read_table("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/pcs_whole_cohort.eigenvec") %>%
  dplyr::select(IID,PC1,PC2,PC3,PC4)
colnames(genetic_pcs)[-1] = paste0("genetic_",colnames(genetic_pcs)[-1])

pcs = pcs %>%
  left_join(genetic_pcs,by="IID") %>%
  na.omit()

# add age and sex
age_sex = sc_dat@meta.data %>% distinct(genotyping_id,.keep_all=T) %>% dplyr::select(genotyping_id,Sex,Age)
pcs = pcs %>%
  left_join(age_sex %>% dplyr::rename("IID" = genotyping_id), by="IID")

# write pc file
covar_file =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/covars_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_",pheno,".tsv")

# add FIDs
pcs = pcs %>%
  left_join(
    fam_for_eqtl %>%
      dplyr::select(X1,X2) %>%
      dplyr::rename("FID"=X1,"IID"=X2),
      by="IID"
  ) %>%
  dplyr::select(FID,IID,PC1,PC2,PC3,PC4,contains("genetic"),Age,Sex)
write_tsv(pcs,file = covar_file)

# filter genelist
genelist$gene_name = str_replace(pattern="-",replace=".",string=genelist$name2)
genelist = genelist %>% filter(gene_name %in% colnames(logcounts_dat))

# loop through all genes
res = foreach(i=1:nrow(genelist)) %dopar% {
  row = i
  message("doing row ",i, "of ",nrow(genelist))
  genelist_dat = genelist[row,]

  gene_name = genelist_dat$gene_name
  message(gene_name)
  chr = genelist_dat$chrom
  start = genelist_dat$txStart
  end = genelist_dat$txEnd

  # filter data to this gene
  dat = logcounts_dat %>% dplyr::rename("donor" = genotyping_id) %>% dplyr::select(donor,all_of(gene_name))

  # filter out NAs
  dat = na.omit(dat)

  # rank normalise
  dat = dat %>%
    mutate(z = RankNorm(.data[[gene_name]]))

  dat = dat %>%
    dplyr::rename("X1" = donor) %>%
    dplyr::select(X1,z) %>%
    dplyr::rename("IID" = X1,"expression_z" = z) %>%
    filter(!is.na(IID))

  # join with FID
  dat = dat %>%
    left_join(pcs %>% dplyr::select(FID,IID),by="IID") %>%
    dplyr::select(FID,IID,expression_z)

  # write pheno file
  pheno_file =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/gene_",gene_name,"_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_",pheno,"_expression.tsv")
  write_tsv(dat,file = pheno_file)

  # get start and end points for eqtl testing
  chr_start = start - 1000000
  chr_end = end + 1000000

  # check there's >=1 SNP

  snps_for_eqtl_this_gene = snps_for_eqtl %>% filter(X1==chr & X4 > chr_start & X4 < chr_end)
  if(nrow(snps_for_eqtl_this_gene)<1){
    message("Insufficient SNPs. Skipping")
  } else {
  # run plink from r
  system(paste0("cd /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/;module load plink;plink2 --bfile merged_sceqtl_genotypes_qc_updated --covar-variance-standardize Age PC1 PC2 PC3 PC4 genetic_PC1 genetic_PC2 genetic_PC3 genetic_PC4 --pheno ",pheno_file," --glm hide-covar --chr ",chr," --covar ",covar_file," --vif 50 --from-bp ",chr_start," --to-bp ",chr_end," --read-freq /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/allele_frequencies.afreq --out /rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_",pheno,"_expression"))

  # read in results
  if(file.exists(paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_",pheno,"_expression.expression_z.glm.linear"))){
    expression_res = read_tsv(paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_",pheno,"_expression.expression_z.glm.linear"),
    col_types = cols_only(`#CHROM` = col_character(),POS = col_double(),ID=col_character(),A1=col_character(), REF=col_character(),ALT=col_character(),OBS_CT = col_double(),BETA=col_double(),SE=col_double(),P=col_double()))

    # output complete table with additional info
    expression_res = expression_res %>%
    mutate(cell_type = cell_type_to_test) %>%
    mutate(source = source_to_test) %>%
    mutate(gene = gene_name) %>%
    mutate(phenotype = pheno)

    # remove tests with NA P values
    expression_res = expression_res %>% filter(!is.na(P))

    # run conditional analysis if SNPs with P < 1e-5
    # add freq
    freqs = read_table("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/allele_frequencies.afreq")

    expression_res = expression_res %>%
      left_join(freqs %>% filter(ID %in% expression_res$ID) %>%
        dplyr::select(-1,-6),by=c("ID","REF","ALT")) %>%
        filter(!is.na(ALT_FREQS)) %>%
        mutate(a1_af = ifelse(A1 == ALT,ALT_FREQS,1-ALT_FREQS))
    expression_res = expression_res %>%
      mutate(A2 = ifelse(A1 == ALT,REF,ALT))


    if(min(expression_res$P)<1e-5){
      # format to ma format
      expression_res_ma = expression_res %>% dplyr::select(ID,A1,A2,a1_af,BETA,SE,P,OBS_CT)

      # write to file
      ma_file =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/scratch/gene_",gene_name,"_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_",pheno,"_expression_gcta.ma")
      write_tsv(expression_res_ma,file = ma_file)

      system(paste0("module load gcta; gcta64 --bfile /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc_updated --cojo-p 1e-5 --cojo-slct --cojo-file ",ma_file," --out ",ma_file))
      cojo_file = read_table(paste0(ma_file,".jma.cojo"))

      expression_res = expression_res %>%
        mutate(cond_sig = ifelse(ID %in% cojo_file$SNP,"cond_sig","not_cond_sig"))
      # clean up
      system(paste0("rm ",ma_file))
      return(expression_res)
    } else {
      system(paste0("rm ",pheno_file))
      expression_res = expression_res %>%
        mutate(cond_sig = NA)
        return(expression_res)
    }

    }
  }
}

expression_res_overall_df = do.call("bind_rows",res) %>% data.frame()


outfile = paste0(cell_type_no_space,"_",source_to_test,"_",pheno,"_eqtls_full_results.tsv")
write_tsv(expression_res_overall_df,file = outfile)

message("Cleaning up")
system(paste0("rm ./scratch/",cell_type_no_space,"_",source_to_test,"*"))
