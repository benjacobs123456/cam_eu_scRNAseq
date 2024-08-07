#! /bin/bash

wd="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
mkdir "${wd}/genotype_data"

# Tools
bcftools="/xx/xxx/bin/bcftools"
tabix="/xx/xxx/bin/htslib-1.9/bin/tabix"
bgzip="/xx/xxx/bin/htslib-1.9/bin/bgzip"
gatk="/xx/xxx/bin/gatk-4.3.0.0/gatk"
plink2="/xx/xxx/bin/plink2"

# CAM TUM
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC --make-bed --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC
$plink2 --bfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC --export vcf --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_gt
$bgzip -c ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_gt.vcf > ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_gt.vcf.gz
$tabix ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_gt.vcf.gz

# 1000G Variant subsetb
awk '{ print $1}' /data_featherstone/shared_data/references/1000Genomes/more/integrated_call_samples_v3.20130502.ALL.panel > ${wd}/genotype_data/QC/1000G_samples.txt
$bcftools view -S ${wd}/genotype_data/QC/1000G_samples.txt --force-samples -O z -o ${wd}/genotype_data/QC/1000g_hg38_subset.vcf.gz ${wd}/genotype_data/PCA/CAMTUM_1000g_hg38.vcf.gz
$tabix ${wd}/genotype_data/QC/1000g_hg38_subset.vcf.gz 

#Merge 
$bcftools isec -n +2 -p ${wd}/genotype_data/QC/shared_variants ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_gt.vcf.gz ${wd}/genotype_data/QC/1000g_hg38_subset.vcf.gz
$bgzip -c ${wd}/genotype_data/QC/shared_variants/0000.vcf > ${wd}/genotype_data/QC/shared_variants/0000.vcf.gz
$bgzip -c ${wd}/genotype_data/QC/shared_variants/0001.vcf > ${wd}/genotype_data/QC/shared_variants/0001.vcf.gz
$tabix ${wd}/genotype_data/QC/shared_variants/0000.vcf.gz
$tabix ${wd}/genotype_data/QC/shared_variants/0001.vcf.gz
$bcftools merge ${wd}/genotype_data/QC/shared_variants/0000.vcf.gz ${wd}/genotype_data/QC/shared_variants/0001.vcf.gz -O z -o ${wd}/genotype_data/QC/CAMTUM_1000g_hg38.vcf.gz
$tabix ${wd}/genotype_data/QC/CAMTUM_1000g_hg38.vcf.gz
awk '{ if ($3 == "EUR") {print $1}}' /data_featherstone/shared_data/references/1000Genomes/more/integrated_call_samples_v3.20130502.ALL.panel > ${wd}/genotype_data/QC/EUR_samples.txt
awk '{ if (!/^#/){ print $1 }}' ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC.psam > ${wd}/genotype_data/QC/merged_samples.txt
cat ${wd}/genotype_data/QC/merged_samples.txt ${wd}/genotype_data/QC/EUR_samples.txt > ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_samples.txt
$bcftools view -S ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_samples.txt --force-samples -O z -o ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38.vcf.gz ${wd}/genotype_data/QC/CAMTUM_1000g_hg38.vcf.gz
cp ${wd}/genotype_data/QC/CAMTUM_1000g_hg38.vcf.gz ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38.vcf.gz


# To plink format
$plink2 --vcf ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38.vcf.gz --make-pgen --out ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38
$plink2 --vcf ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38.vcf.gz --make-pgen --out ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38

# Calculate PCs
$plink2 --pfile ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38 --missing --out ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38_miss
$plink2 --pfile ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38 --exclude range ${wd}/genotype_data/QC/remove_prune.txt --maf 0.05 --hwe 1e-3 --geno 0.05 --indep-pairwise 200 100 0.2 --out ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38_prunelist
$plink2 --pfile ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38 --extract ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38_prunelist.prune.in --pca --out ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38_pca
$plink2 --pfile ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38 --missing --out ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38_miss
$plink2 --pfile ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38 --exclude range ${wd}/genotype_data/QC/remove_prune.txt --maf 0.05 --hwe 1e-3 --geno 0.05 --indep-pairwise 200 100 0.2 --out ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38_prunelist
$plink2 --pfile ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38 --extract ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38_prunelist.prune.in --pca --out ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38_pca

# PCA plots
cohort_info_file="/data_featherstone/shared_data/references/1000Genomes/integrated_call_samples_v3.20130502.ALL.panel"
Rscript ${wd}/Scripts/pca_plots_1000G.R ${wd}/genotype_data/QC/CAMTUM_1000g_EUR_hg38_pca.eigenvec $cohort_info_file EUR ${wd}/genotype_data/QC/PCplot_EUR_afterQC.png
Rscript ${wd}/Scripts/pca_plots_1000G.R ${wd}/genotype_data/QC/CAMTUM_1000g_ALL_hg38_pca.eigenvec $cohort_info_file ALL ${wd}/genotype_data/QC/PCplot_ALL_afterQC.png
