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

#######################################
# Read in data
#######################################
# get args
args = commandArgs(trailingOnly = TRUE)
message(args)

# set WD
setwd("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/eqtl/")

# Read in data
all_combo = readRDS("../datasets/all_combo_with_updated_pheno.rds")

# add full cell id
all_combo@meta.data = all_combo@meta.data %>%
  mutate(full_cell_id = paste0(source,"_",donor.id,"_",rownames(all_combo@meta.data)))
rownames(all_combo@meta.data) = colnames(all_combo)

# filter to just Cam data
all_combo = subset(all_combo,cohort=="Cam")

# rename cells with full ID
all_combo = RenameCells(all_combo,new.names = all_combo@meta.data$full_cell_id)

# create file with new IDs for plink
codex = read_table("/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/data/genotypes/GT_ID.txt",col_names=FALSE)

# read in SNPs
snps_for_eqtl = read_tsv("/rds/user/hpcjaco1/hpc-work/filtered_scRNAseq_genotypes_plink.pvar",
  skip = 40)

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
genelist = read_tsv("/rds/user/hpcjaco1/hpc-work/genes") # downloaded from UCSC
genelist = genelist %>%
  filter(name2 %in% rownames(all_combo)) %>%
  distinct(name2,.keep_all=T) %>%
  mutate(chrom = str_remove(chrom,"chr")) %>%
  filter(chrom != "X" & chrom != "Y")


# initialise results list
expression_res_overall = list()

# get cell type
cell_type_to_test =  unique(all_combo@meta.data$cell_type_crude)[as.numeric(args[1])]
cell_type_no_space =  str_replace_all(cell_type_to_test," ","_")

for(source_to_test in c("CSF","PBMC")){

  # filter sc data
  sc_dat = subset(all_combo,source==source_to_test)

  message(source_to_test)
  message(cell_type_to_test)

  # filter to cell type
  rownames(sc_dat@meta.data) = colnames(sc_dat)
  sc_dat = subset(sc_dat,subset = cell_type_crude == cell_type_to_test)

  # restrict to donors with >=2 cells
  donors_to_keep = sc_dat@meta.data %>% dplyr::count(donor.id) %>% filter(n>=2)
  sc_dat = subset(sc_dat,subset = donor.id %in% donors_to_keep$donor.id)

  # convert to sce object
  sc_dat.sce = as.SingleCellExperiment(sc_dat)

  # get rid of low-expressed genes
  sc_dat.sce = sc_dat.sce[rowSums(counts(sc_dat.sce))>100,]

  # get rid of mitochondrial and ribosomal genes
  sc_dat.sce = sc_dat.sce[!grepl("^MT",rownames(sc_dat.sce)),]
  sc_dat.sce = sc_dat.sce[!grepl("^RP",rownames(sc_dat.sce)),]

  # get counts
  counts = assay(sc_dat.sce,"counts")

  # transform to log2(cpm)
  libsizes = colSums(counts(sc_dat.sce))
  logcounts(sc_dat.sce) = log2(t(t(counts)/libsizes*1e6) + 1)

  # aggregate
  groups = colData(sc_dat.sce)[,"donor.id"]
  logcounts_dat = t(logcounts(sc_dat.sce)) %>%
    data.frame()
  logcounts_dat$full_cell_id = rownames(logcounts_dat)
  logcounts_dat = logcounts_dat %>%
    tidyr::separate(full_cell_id,sep="_",into = c("source","donor","cell_id"))
  logcounts_dat = logcounts_dat %>%
    group_by(donor) %>%
    summarise_at(vars(-source,-cell_id),.fun = "mean", na.omit = T)

  # get pcs
  pcs = prcomp(logcounts_dat %>% dplyr::select(-donor),
    scale = T,
    center = T)

  pcs = pcs$x %>% data.frame()
  pcs$donor.id = logcounts_dat$donor

  # write to file
  # get genotype ID

  pcs = pcs %>%
    left_join(codex %>% distinct(donor.id,.keep_all=T),by="donor.id") %>%
    dplyr::select(X1,PC1:4) %>%
    dplyr::rename("IID" = X1) %>%
    filter(!is.na(IID))

  # write pc file
  covar_file =  paste0("/rds/user/hpcjaco1/hpc-work/covars_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_all.tsv")
  write_tsv(pcs,file = covar_file)


  # filter genelist
  genelist$gene_name = str_replace(pattern="-",replace=".",string=genelist$name2)
  genelist = genelist %>% filter(gene_name %in% colnames(logcounts_dat)) %>%
    filter(!is.na(as.numeric(chrom)))

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
    dat = logcounts_dat %>% dplyr::select(donor,all_of(gene_name))

    # filter out NAs
    dat = na.omit(dat)

    # rank normalise
    dat = dat %>%
      mutate(z = RankNorm(.data[[gene_name]]))

    dat = dat %>%
      dplyr::rename("donor.id" = donor) %>%
      left_join(codex %>% distinct(donor.id,.keep_all=T),by="donor.id") %>%
      dplyr::select(X1,z) %>%
      dplyr::rename("IID" = X1,"expression_z" = z) %>%
      filter(!is.na(IID))

    # write pheno file
    pheno_file =  paste0("/rds/user/hpcjaco1/hpc-work/gene_",gene_name,"_source_",source_to_test,"_cell_type_",cell_type_no_space,"_pheno_all_expression.tsv")
    write_tsv(dat,file = pheno_file)

    # get start and end points for eqtl testing
    chr_start = start - 100000
    chr_end = end + 100000

    # check there's >=1 SNP

    snps_for_eqtl_this_gene = snps_for_eqtl %>% filter(`#CHROM`==chr & POS > chr_start & POS < chr_end)
    if(nrow(snps_for_eqtl_this_gene)<1){
      message("Insufficient SNPs. Skipping")
      next
    } else {

    # run plink from r
    system(paste0("cd /rds/user/hpcjaco1/hpc-work/;module load plink;plink2 --pfile filtered_scRNAseq_genotypes_plink --pheno ",pheno_file," --glm hide-covar --chr ",chr," --covar ",covar_file," --vif 50 --from-bp ",chr_start," --to-bp ",chr_end," --out ",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_all_expression"))

    # read in results
    if(file.exists(paste0("/rds/user/hpcjaco1/hpc-work/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_all_expression.expression_z.glm.linear"))){
      expression_res = read_tsv(paste0("/rds/user/hpcjaco1/hpc-work/",cell_type_no_space,"_",source_to_test,"_",gene_name,"_pheno_all_expression.expression_z.glm.linear"),
      col_types = cols_only(`#CHROM` = col_character(),POS = col_double(),ID=col_character(),A1=col_character(), REF=col_character(),ALT=col_character(),BETA=col_double(),SE=col_double(),P=col_double()))

      # output complete table with additional info
      expression_res = expression_res %>%
      mutate(cell_type = cell_type_to_test) %>%
      mutate(source = source_to_test) %>%
      mutate(gene = gene_name) %>%
      mutate(phenotype = "all")

      expression_res_overall[[length(expression_res_overall)+1]] = expression_res
      }
    }
  }
}


expression_res_overall_df = do.call("bind_rows",expression_res_overall) %>% data.frame()
outfile = paste0(cell_type_no_space,"_eqtls_full_results_whole_cohort.tsv")
write_tsv(expression_res_overall_df,file = outfile)
