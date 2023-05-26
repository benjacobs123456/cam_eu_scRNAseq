library(tidyverse)
setwd("/home/hpcjaco1/rds/hpc-work/Cambridge_EU_combined/eqtl/")

# read in
eqtl_res = readRDS("full_eqtl_results_all.rds")

# add fdr
eqtl_res = eqtl_res %>%
  mutate(fdr = p.adjust(P,method="fdr"))


# define function
find_specific_eqtls = function(phenotype_to_test,source_to_test){

  # filter dat
  dat_for_specific_eqtls = eqtl_res %>%
    filter(phenotype==phenotype_to_test & source==source_to_test)

  # get cell types
  cell_types = unique(dat_for_specific_eqtls$cell_type)

  # loop through each cell type
  overall_res = data.frame()

  for(i in c(1:length(cell_types))){
    # get cell name
    this_cell_type = cell_types[i]
    message(this_cell_type)

    # filter to sig eqtls for this cell
    eqtls_this_cell = dat_for_specific_eqtls %>%
      filter(cell_type == this_cell_type & fdr <= 0.05)

    # loop through other cell types
    for(j in c(1:length(cell_types))){
      target_cell_type = cell_types[j]
      message(target_cell_type)

      # filter to sig eqtls for this cell
      eqtls_target_cell = ms_pbmc_eqtl %>%
        filter(cell_type == target_cell_type)

      # merge
      merged_eqtls = eqtls_this_cell %>%
        left_join(eqtls_target_cell,by=c("X.CHROM","POS","ID","A1","REF","ALT","gene","source"))

      # calculate het stat
      merged_eqtls = merged_eqtls %>%
          mutate(het_z = (BETA.x - BETA.y) / (sqrt ( SE.x^2 + SE.y^2) )  ) %>%
          mutate(het_p = 1 - pnorm(abs(het_z))) %>%
          dplyr::select(ID,A1,het_p,BETA.x,SE.x,P.x,BETA.y,SE.y,P.y,gene) %>%
          mutate(cell_type1 = this_cell_type,cell_type2 = target_cell_type)

      # add to overall res list
      overall_res <<- bind_rows(overall_res,merged_eqtls)
    }
  }


# spread wider
res_wide = overall_res %>%
filter(cell_type1 != cell_type2) %>%
  pivot_wider(id_cols = c(ID,A1,gene,cell_type1),
    values_from = het_p,
    names_from = cell_type2,
    names_prefix="HET_P_")

# get het eqtls
het_eqtls = res_wide %>%
  filter_at(vars(contains("HET_P")),
  all_vars(. < 0.05 | is.na(.))) %>%
  filter_at(vars(contains("HET_P")),
  any_vars(!is.na(.)))

# add in total n sig
sig_eqtl_per_celltype = dat_for_specific_eqtls %>%
  filter(fdr<0.05) %>%
  dplyr::count(cell_type) %>%
  dplyr::rename("cell_type1"=cell_type)

prop_specific_eqtl = het_eqtls %>%
  dplyr::count(cell_type1) %>%
  left_join(sig_eqtl_per_celltype,by="cell_type1") %>%
  mutate(prop_specific = n.x/n.y)
p1=ggplot(het_eqtls,aes(cell_type1))+
geom_bar(color="black")+
theme_minimal()+
labs(x="Cell type",y="N cell-type-specific eQTLs")+
coord_flip()
p2=ggplot(prop_specific_eqtl,aes(cell_type1,prop_specific*100))+
geom_col(color="black")+
theme_minimal()+
labs(x="Cell type",y="% cell-type-specific eQTLs")+
coord_flip()
png(paste0(phenotype_to_test,"_",source_to_test,"_cell_type_specific.png"),res=600,units="in",width=6,height=4)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

# write to file
write_csv(het_eqtls,paste0(phenotype_to_test,"_",source_to_test,"_het_eqtls.csv"))
}

# run
find_specific_eqtls("MS","CSF")

# het plots
make_het_plot = function(this_snp,this_gene){
plot_dat = eqtl_res %>%
  filter(ID == this_snp & gene == this_gene)
p=ggplot(plot_dat,aes(BETA,cell_type,col=phenotype))+
  geom_point(position=ggstance::position_dodgev(height=0.3))+
  facet_wrap(~source)+
  geom_errorbarh(mapping = aes(xmin = BETA - 1.96*SE, xmax = BETA + 1.96*SE,y=cell_type),height=0.1,position=ggstance::position_dodgev(height=0.3))+
  theme_minimal()+
  ggtitle(paste0("SNP: ",this_snp,"\nGene: ",this_gene))+
  geom_vline(linetype="dashed",color="black",xintercept=0)+
  scale_color_brewer(palette="Set1")+
  labs(y="Cell type")
p
}


make_het_plot("3:111522059:G:A","CD96")
make_het_plot("6:32243441:A:G","HLA.DRB1")
make_het_plot("1:26251494:G:A","CD52")

# coloc hits
files = list.files(pattern="^coloc")
files = files[grepl("_MS.csv",files)]
coloc_filelist = purrr::map(files,function(x){
  read_csv(x,col_types = cols(.default = "c"))
})
coloc_filelist = do.call("bind_rows",coloc_filelist)

coloc_filelist$chrpos = coloc_filelist$snp
coloc_filelist = coloc_filelist %>% filter(source=="PBMC" & phenotype=="MS")
coloc_filelist$cell_type1 = coloc_filelist$cell

het_eqtl_cell_type_specific = het_eqtls %>% tidyr::separate(ID,sep=":",into=c("chr","pos","a1","a2")) %>%
  mutate(chrpos = paste0(chr,":",pos)) %>%
  filter(chrpos %in% coloc_filelist$snp) %>%
  left_join(coloc_filelist %>% dplyr::select(chrpos,SNP.PP.H4,gene,cell_type1),by=c("chrpos","gene","cell_type1")) %>%
  filter(!is.na(as.numeric(SNP.PP.H4)))


ggplot(het_eqtl_cell_type_specific,aes(as.numeric(SNP.PP.H4)))+
geom_histogram()+
theme_minimal()

