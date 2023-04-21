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

#######################################
# Read in data
#######################################
# get args
args = commandArgs(trailingOnly = TRUE)
message(args)

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/")

# filter to cell type
cell_type_to_test =  unique(eqtl_dataset@meta.data$cell_type_crude)[1]
cell_type_no_space =  str_replace_all(cell_type_to_test," ","_")

# read in file file
infile = paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/eqtl_sc_data_",cell_type_no_space)
eqtl_dataset = readRDS(infile)

# read in SNPs
snps_for_eqtl = read_tsv("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/merged_sceqtl_genotypes_qc.bim",col_names=F)
message("SNPs for eQTL testing")
message(nrow(snps_for_eqtl))

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

for(pheno in c("MS","OIND","Noninflammatory")){


  for(source_to_test in c("CSF","PBMC")){

    # filter sc data
    sc_dat = subset(eqtl_dataset,source==source_to_test & phenotype == pheno)

    message(source_to_test)
    message(cell_type_to_test)

    rownames(sc_dat@meta.data) = colnames(sc_dat)

    # restrict to donors with >=2 cells
    donors_to_keep = sc_dat@meta.data %>% dplyr::count(genotyping_id) %>% filter(n>=2)
    sc_dat = subset(sc_dat,subset = genotyping_id %in% donors_to_keep$genotyping_id)

    # convert to sce object
    sc_dat.sce = as.SingleCellExperiment(sc_dat)

    # get counts
    counts = assay(sc_dat.sce,"counts")

    # transform to log2(cpm)
    libsizes = colSums(counts(sc_dat.sce))
    logcounts(sc_dat.sce) = log2(t(t(counts)/libsizes*1e6) + 1)

    # get rid of low-expressed genes
    sc_dat.sce = sc_dat.sce[rowSums(counts(sc_dat.sce))>100,]

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
      dplyr::select(donor.id,PC1:4) %>%
      dplyr::rename("IID" = donor.id) %>%
      filter(!is.na(IID))

    # write pc file
    covar_file =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/covars_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_",pheno,".tsv")
    write_tsv(pcs,file = covar_file)


    # filter genelist
    genelist$gene_name = str_replace(pattern="-",replace=".",string=genelist$name2)
    genelist = genelist %>% filter(gene_name %in% colnames(logcounts_dat))

    # loop through all genes
    overall_results = list()
    for(i in c(1:nrow(genelist))){
      row = i
      message("doing row ",i)
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

      # write pheno file
      pheno_file =  paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/gene_",gene_name,"_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_",pheno,"_expression.tsv")
      write_tsv(dat,file = pheno_file)

      # get start and end points for eqtl testing
      chr_start = start - 100000
      chr_end = end + 100000

      # check there's >=1 SNP

      snps_for_eqtl_this_gene = snps_for_eqtl %>% filter(X1==chr & X4 > chr_start & X4 < chr_end)
      if(nrow(snps_for_eqtl_this_gene)<1){
        message("Insufficient SNPs. Skipping")
        next
      } else {

      # run plink from r
      system(paste0("cd /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/;module load plink;plink2 --bfile merged_sceqtl_genotypes_qc --pheno ",pheno_file," --glm hide-covar --chr ",chr," --covar ",covar_file," --vif 50 --from-bp ",chr_start," --to-bp ",chr_end," --out /rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_",pheno,"_expression"))

      # read in results
      if(file.exists(paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_",pheno,"_expression.expression_z.glm.linear"))){
        expression_res = read_tsv(paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/eqtl/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_",pheno,"_expression.expression_z.glm.linear"),
        col_types = cols_only(`#CHROM` = col_character(),POS = col_double(),ID=col_character(),A1=col_character(), REF=col_character(),ALT=col_character(),BETA=col_double(),SE=col_double(),P=col_double()))

        # output complete table with additional info
        expression_res = expression_res %>%
        mutate(cell_type = cell_type_to_test) %>%
        mutate(source = source_to_test) %>%
        mutate(gene = gene_name) %>%
        mutate(phenotype = pheno)

        expression_res_overall[[length(expression_res_overall)+1]] = expression_res
        }
      }
    }
  }
}

expression_res_overall_df = do.call("bind_rows",expression_res_overall) %>% data.frame()
outfile = paste0(cell_type_no_space,"_eqtls_full_results.tsv")
write_tsv(expression_res_overall_df,file = outfile)
