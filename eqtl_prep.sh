####################################
# CAMBRIDGE
####################################


module load bcftools-1.9-gcc-5.4.0-b2hdt5n
cd /rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/data/genotypes

# make rename map file
for i in {1..22};
  do
    echo $i chr$i >> mapfile
  done

# rename chrs
bcftools \
  annotate --rename-chrs mapfile SAWC_TPS_all_samples_22_10_2020.vcf.gz -Oz > genos_all_chrs.vcf.gz
bcftools index genos_all_chrs.vcf.gz

for i in {22..1};
  do
    echo $i
    bcftools view genos_all_chrs.vcf.gz -r chr$i -Oz > sc_geno_chr$i\.vcf.gz
  done

cd /rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/imputed_genotype_data/

# imputed with topmed-r2 at https://imputation.biodatacatalyst.nhlbi.nih.gov/start.html#!jobs/job-20221018-083007-316/results
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/725613/eb29179528389c61d5d0bfd071dd0a246f68c6a27b19842096c7a3fc5ebe7e29 | bash

for i in {1..22}
  do
    7z x chr_$i\.zip
  done

# filter by info and MAF
for i in {22..1}
  do
    echo $i
    bcftools index -f chr$i\.dose.vcf.gz
    bcftools filter -i 'INFO/R2>0.9 & INFO/MAF>0.05' chr$i\.dose.vcf.gz -Oz > filtered_chr$i\.vcf.gz
  done

# concatenate
cd /rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/imputed_genotype_data/
rm merge_filelist
for i in {1..22}
  do
    echo filtered_chr$i\.vcf.gz >> merge_filelist
  done

bcftools concat -f merge_filelist -Oz > merged_vcf.vcf.gz

# filter with plink

module load plink
plink2 --vcf merged_vcf.vcf.gz \
--make-pgen \
--out scRNAseq_genotypes_plink

plink2 --pfile scRNAseq_genotypes_plink \
--maf 0.05 --hwe 1e-5 --geno 0.1 --mind 0.1 \
--out /rds/user/hpcjaco1/hpc-work/filtered_scRNAseq_genotypes_plink \
--make-pgen

####################################
# TUM
####################################

# cambridge batch
cd /rds/user/hpcjaco1/hpc-work/tum_genos_round2/
module load plink
plink2 --vcf GSA_2020_formatted_renamed.vcf.gz dosage=GP \
--out cam_tum_plink \
--make-bed

# reset var ids
plink2 --bfile cam_tum_plink \
--set-all-var-ids @:#:\$r:\$a \
--make-bed \
--snps-only \
--out cam_tum_plink_newids  \
--new-id-max-allele-len 577

# liftOver to hg38
awk '{print "chr"$1,$4-1,$4,$2}' cam_tum_plink_newids.bim > cam_tum_plink_newids_hg19.bed
/home/hpcjaco1/liftover/liftOver cam_tum_plink_newids_hg19.bed /home/hpcjaco1/liftover/hg19ToHg38.over.chain.gz cam_tum_plink_newids_hg38.bed unlifted.bed

awk '{print $4,$3}' cam_tum_plink_newids_hg38.bed > updated_snp_map_hg38

plink2 --bfile cam_tum_plink_newids \
--update-map updated_snp_map_hg38 \
--make-bed \
--out cam_tum_plink_hg38

# basic QC
plink2 --bfile cam_tum_plink_hg38 \
--set-all-var-ids @:# \
--mind 0.1 \
--maf 0.05 --hwe 1e-5 --geno 0.1 \
--out /rds/user/hpcjaco1/hpc-work/genotypes/filtered_scRNAseq_genotypes_plink_TUM_batch2 \
--make-pgen


cd /rds/user/hpcjaco1/hpc-work/TUM\ data/data_featherstone/christiane/SC/transfer_CAM/TUM_part2/

# main batch
# convert to PLINK format
awk 'NR>1{print $2,$1,0,$3}' map_file_small.txt > updated_snp_map

module load plink
for i in {1..22};
  do
    plink2 --import-dosage chr$i\_qc.dosage noheader \
    --fam samples.fam \
    --map updated_snp_map \
    --out chr$i \
    --make-bed &
  done

# merge
rm filelist_for_merge
for i in {2..22};
  do
    echo chr$i >> filelist_for_merge
  done

module unload plink
module load plink-1.9-gcc-5.4.0-sm3ojoi

plink --bfile chr1 \
--merge-list filelist_for_merge \
--make-bed \
--out merged_all_chroms

# repeat excluded triallelics
for i in {1..22};
  do
    plink --bfile chr$i \
    --exclude merged_all_chroms-merge.missnp \
    --make-bed \
    --out no_triallelics_chr$i
  done

rm filelist_for_merge
for i in {2..22};
  do
    echo no_triallelics_chr$i >> filelist_for_merge
  done

plink --bfile no_triallelics_chr1 \
--merge-list filelist_for_merge \
--make-bed \
--out merged_all_chroms

# liftOver to hg38
awk '{print "chr"$1,$4-1,$4,$2}' merged_all_chroms.bim > tum_hg19.bed
/home/hpcjaco1/liftover/liftOver tum_hg19.bed /home/hpcjaco1/liftover/hg19ToHg38.over.chain.gz tum_hg38.bed unlifted.bed

awk '{print $4,$3}' tum_hg38.bed > updated_snp_map_hg38

plink --bfile merged_all_chroms \
--update-map updated_snp_map_hg38 \
--make-bed \
--out merged_all_chroms_hg38

# basic QC
module unload plink
module load plink
plink2 --bfile merged_all_chroms_hg38 \
--set-all-var-ids @:# \
--maf 0.05 --hwe 1e-5 --geno 0.1 --mind 0.1 \
--out /rds/user/hpcjaco1/hpc-work/genotypes/filtered_scRNAseq_genotypes_plink_TUM \
--make-pgen

#####################################
# MERGE
#####################################

cd /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes

# first pass merge
module unload plink
module load plink

# covert to plink1
plink2 --pfile /rds/user/hpcjaco1/hpc-work/genotypes/filtered_scRNAseq_genotypes_plink_TUM \
--set-all-var-ids @:#:\$r:\$a \
--snps-only \
--make-bed \
--out plink_TUM

plink2 --pfile /rds/user/hpcjaco1/hpc-work/genotypes/filtered_scRNAseq_genotypes_plink_TUM_batch2 \
--set-all-var-ids @:#:\$r:\$a \
--snps-only \
--make-bed \
--out plink_TUM2

plink2 --pfile /rds/user/hpcjaco1/hpc-work/filtered_scRNAseq_genotypes_plink \
--set-all-var-ids @:#:\$r:\$a \
--make-bed \
--snps-only \
--out plink_CAM

# get allele freqs
plink2 --bfile plink_CAM \
--freq \
--out cam_freqs

plink2 --bfile plink_TUM \
--freq \
--out tum_freqs

plink2 --bfile plink_TUM2 \
--freq \
--out tum_freqs2

````R
library(tidyverse)
setwd("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes")
tum = read_table("plink_TUM.bim",col_names=F) %>%
mutate(chrpos = paste0(X1,":",X4))
tum2 = read_table("plink_TUM2.bim",col_names=F) %>%
mutate(chrpos = paste0(X1,":",X4))
cam = read_table("plink_CAM.bim",col_names=F) %>%
mutate(chrpos = paste0(X1,":",X4))

# combine
all_snps = bind_rows(tum,cam,tum2)

# filter
all_snps = all_snps %>%
  distinct(X1,X4,.keep_all=T) %>% # restrict to unique positions
  filter(chrpos %in% tum$chrpos & chrpos %in% cam$chrpos & chrpos %in% tum2$chrpos) # restrict to snps in all datasets

tum = tum %>% filter(chrpos %in% all_snps$chrpos)
tum2 = tum2 %>% filter(chrpos %in% all_snps$chrpos)
cam = cam %>% filter(chrpos %in% all_snps$chrpos)

# join to look for compatible alleles
combo = tum %>%
  left_join(cam,by=c("chrpos","X1","X4")) %>%
  left_join(tum2,by = c("chrpos","X1","X4") )

# compatible
combo = combo %>%
  filter(X5.x == X5.y & X6.x == X6.y) %>%
  filter(X5.x == X5 & X6.x == X6) %>%
  filter(X5 == X5.y & X6 == X6.y)

# get rid of palindromes
combo = combo %>%
filter(
!(X5.x == "C" & X6.y == "G") &
!(X5.x == "G" & X6.y == "C") &
!(X5.x == "A" & X6.y == "T") &
!(X5.x == "T" & X6.y == "A"))

combo = combo %>%
filter(
!(X5.x == "C" & X6 == "G") &
!(X5.x == "G" & X6 == "C") &
!(X5.x == "A" & X6 == "T") &
!(X5.x == "T" & X6 == "A"))

combo = combo %>%
filter(
!(X5.y == "C" & X6 == "G") &
!(X5.y == "G" & X6 == "C") &
!(X5.y == "A" & X6 == "T") &
!(X5.y == "T" & X6 == "A"))


# read in freqs
cam_freqs = read_table("cam_freqs.afreq") %>% filter(ID %in% combo$X2.x) %>% dplyr::select(ID,ALT_FREQS,ALT) %>% dplyr::rename("cam_FREQ" = ALT_FREQS)
tum_freqs = read_table("tum_freqs.afreq") %>% filter(ID %in% combo$X2.x) %>% dplyr::select(ID,ALT_FREQS,ALT) %>% dplyr::rename("tum_FREQ" = ALT_FREQS)
tum2_freqs = read_table("tum_freqs2.afreq") %>% filter(ID %in% combo$X2.x) %>% dplyr::select(ID,ALT_FREQS,ALT) %>% dplyr::rename("tum2_FREQ" = ALT_FREQS)

# compare freqs
all_freqs = cam_freqs %>%
  left_join(tum_freqs,by=c("ID","ALT")) %>%
  left_join(tum2_freqs,by=c("ID","ALT"))

all_freqs = all_freqs %>%
  mutate(delta_freq_camtum = abs( cam_FREQ - tum_FREQ )) %>%
  mutate(delta_freq_camtum2 = abs( cam_FREQ - tum2_FREQ )) %>%
  mutate(delta_freq_tumtum = abs( tum_FREQ - tum2_FREQ ))


ggplot(all_freqs,aes(cam_FREQ,tum_FREQ))+
geom_point()+
scale_x_continuous(limits=c(0,1))+
scale_y_continuous(limits=c(0,1))+
geom_abline(intercept=0,slope=1)


ggplot(all_freqs,aes(delta_freq))+
geom_histogram()


write_tsv(combo %>% dplyr::select(X2.x),"snps_to_include.tsv",col_names=F)

````

# restrict to these snps
# covert to plink1
plink2 --bfile plink_TUM \
--extract snps_to_include.tsv \
--make-bed \
--out plink_TUM_formerge

plink2 --bfile plink_TUM2 \
--extract snps_to_include.tsv \
--make-bed \
--out plink_TUM2_formerge

plink2 --bfile plink_CAM \
--extract snps_to_include.tsv \
--make-bed \
--out plink_CAM_formerge

# now merge
module unload plink
module load plink-1.9-gcc-5.4.0-sm3ojoi

echo "plink_CAM_formerge" > final_merge_list
echo "plink_TUM2_formerge" >> final_merge_list
plink --bfile plink_TUM_formerge \
--merge-list final_merge_list \
--make-bed \
--out merged_sceqtl_genotypes

# filter by missingness

plink --bfile merged_sceqtl_genotypes \
--geno 0.1 \
--maf 0.01 \
--hwe 1e-5 \
--make-bed \
--out merged_sceqtl_genotypes_qc

cd /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes
plink --bfile merged_sceqtl_genotypes_qc --freq

########################
# QC with 1kg
########################

# convert to plink
for i in {1..22};
  do
    plink --vcf /rds/user/hpcjaco1/hpc-work/reference_1kg/ALL.chr$i\.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    --maf 0.05 \
    --hwe 1e-5 \
    --snps-only \
    --geno 0.1 \
    --mind 0.1 \
    --make-bed \
    --out /rds/user/hpcjaco1/hpc-work/reference_1kg/chr$i &
  done

# attempt to merge
for i in {2..22};
  do
    echo /rds/user/hpcjaco1/hpc-work/reference_1kg/chr$i >> /rds/user/hpcjaco1/hpc-work/reference_1kg/kg_mergelist
  done

cd /rds/user/hpcjaco1/hpc-work/reference_1kg
plink --bfile chr1 \
--merge-list kg_mergelist \
--out merged_kg \
--make-bed

# run again without triallelics
for i in {1..22};
  do
    plink --bfile chr$i \
    --make-bed \
    --exclude merged_kg-merge.missnp \
    --out merge_chr$i &
  done

rm kg_mergelist
for i in {2..22};
  do
    echo merge_chr$i >> kg_mergelist
  done

cd /rds/user/hpcjaco1/hpc-work/reference_1kg
plink --bfile merge_chr1 \
--merge-list kg_mergelist \
--out merged_kg \
--make-bed

# liftover to hg38
# liftOver to hg38
awk '{print "chr"$1,$4-1,$4,$2}' merged_kg.bim > kg_hg19.bed
/home/hpcjaco1/liftover/liftOver kg_hg19.bed /home/hpcjaco1/liftover/hg19ToHg38.over.chain.gz kg_hg38.bed unlifted.bed
awk '{print $4,$3}' kg_hg38.bed > updated_snp_map_hg38

module unload plink
module load plink

cd /rds/user/hpcjaco1/hpc-work/reference_1kg/
plink2 --bfile merged_kg \
--update-map updated_snp_map_hg38 \
--make-bed \
--out kg_merged_all_chroms_hg38

plink2 --bfile kg_merged_all_chroms_hg38 \
--set-all-var-ids @:#:\$r:\$a \
--make-bed \
--out kg_merged_all_chroms_hg38_chrpos

# now filter to SNPs in main dataset
cd /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes

plink --bfile /rds/user/hpcjaco1/hpc-work/reference_1kg/kg_merged_all_chroms_hg38_chrpos \
--extract merged_sceqtl_genotypes_qc.bim \
--out kg_for_merge \
--make-bed

# filter sceqtl genos to these snps only
plink --bfile merged_sceqtl_genotypes_qc \
--extract kg_for_merge.bim \
--out merged_sceqtl_genotypes_for_kg_merge \
--make-bed

# now merge
plink --bfile merged_sceqtl_genotypes_for_kg_merge \
--bmerge kg_for_merge \
--make-bed \
--out sceqtl_kg_merge

# basic filtering
plink --bfile sceqtl_kg_merge \
--geno 0.1 \
--mind 0.1 \
--maf 0.05 \
--hwe 1e-5 \
--out sceqtl_kg_merge_qc \
--make-bed

# prune
plink --bfile sceqtl_kg_merge_qc \
--indep-pairwise 1000 100 0.01 \
--out pruned_snps

plink --bfile sceqtl_kg_merge_qc \
--extract pruned_snps.prune.in \
--out pruned_for_pca \
--make-bed


# PCA
plink --bfile pruned_for_pca \
--pca


# PCA without 1kg



# prune
plink --bfile merged_sceqtl_genotypes_qc \
--indep-pairwise 1000 100 0.01 \
--out pruned_snps_just_tum_cam

plink --bfile merged_sceqtl_genotypes_qc \
--extract pruned_snps_just_tum_cam.prune.in \
--out pruned_for_pca_just_tum_cam \
--make-bed

plink --bfile pruned_for_pca_just_tum_cam \
--pca \
--out pca_just_tum_cam


````R
library(tidyverse)
setwd("/rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes")
df = read_table("plink.eigenvec",col_names=F)
pcs_just_tum_cam = read_table("pca_just_tum_cam.eigenvec",col_names=F)

tum_pheno = read_table("plink_TUM.fam",col_names=F)
tum2_pheno = read_table("plink_TUM2.fam",col_names=F)
cam_pheno = read_table("plink_CAM.fam",col_names=F)

df = df %>%
mutate(cohort = case_when(
  X2 %in% tum_pheno$X2 ~ "TUM",
  X2 %in% tum2_pheno$X2 ~ "TUM2",
  X2 %in% cam_pheno$X2 ~ "CAM",
  !(X2 %in% cam_pheno$X2) & !(X2 %in% tum_pheno$X2) & !(X2 %in% tum2_pheno$X2)~ "1KG"
  ))


kg_pops = read_csv("/rds/user/hpcjaco1/hpc-work/genotypes/20130606_sample_info.csv") %>%
dplyr::select(1,3) %>%
dplyr::rename("X2" = Sample)

# just cohort

pcs_just_tum_cam = pcs_just_tum_cam %>%
mutate(cohort = case_when(
  X2 %in% tum_pheno$X2 ~ "TUM",
  X2 %in% tum2_pheno$X2 ~ "TUM2",
  X2 %in% cam_pheno$X2 ~ "CAM"))
png("pc_plot_just_cohort.png",res=600,width=6,height=4,units="in")
ggplot(pcs_just_tum_cam,aes(X3,X4,col=cohort))+geom_point()+
labs(x="PC1",y="PC2")+theme_minimal()
dev.off()


df = df %>% left_join(kg_pops,by="X2")

png("pc_plot.png",res=600,width=6,height=4,units="in")
ggplot(df,aes(X3,X4,col=Population))+geom_point()+facet_wrap(~cohort)+
labs(x="PC1",y="PC2")+theme_minimal()
dev.off()

superpops = read_tsv("igsr_populations.tsv")
superpops = superpops %>% dplyr::select(`Population code`,`Superpopulation code`)
colnames(superpops) = c("Population","Ancestry")
df = df %>% left_join(superpops,by="Population")
png("pc_plot2.png",res=600,width=6,height=4,units="in")
ggplot(df,aes(X3,X4,col=Ancestry))+geom_point()+facet_wrap(~cohort)+
labs(x="PC1",y="PC2")+theme_minimal()
dev.off()


png("pc_plot3.png",res=600,width=6,height=4,units="in")
ggplot(df %>% filter(is.na(Ancestry) | Ancestry=="EUR" ),aes(X3,X4,col=Population))+geom_point()+facet_wrap(~cohort)+
labs(x="PC1",y="PC2")+theme_minimal()
dev.off()

````


# UPDATE SEX AND FID INFO
# add FIDs
cd /rds/user/hpcjaco1/hpc-work/genotypes/sceqtl_genotypes/
# make update map
awk '{print $1,$2,$2,$2}' merged_sceqtl_genotypes_qc.fam > id_map_for_update

module load plink
plink2 --bfile merged_sceqtl_genotypes_qc \
--update-ids id_map_for_update \
--out tmp_updated_ids \
--make-bed

# update sex info
awk '{print $1,$2,0}' tmp_updated_ids.fam > sex_map_for_update

# set all sex values to missing
plink2 --bfile tmp_updated_ids \
--update-sex sex_map_for_update \
--out merged_sceqtl_genotypes_qc_updated \
--make-bed

# get freqs
plink2 --bfile merged_sceqtl_genotypes_qc_updated \
--freq \
--out allele_frequencies

# get pcs
plink2 --bfile merged_sceqtl_genotypes_qc_updated \
--pca \
--out pcs_whole_cohort
