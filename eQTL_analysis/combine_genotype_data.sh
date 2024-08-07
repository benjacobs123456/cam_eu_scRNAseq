wd="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
mkdir "${wd}/genotype_data"

# Tools
bcftools="/xx/xxx/bin/bcftools"
tabix="/xx/xxx/bin/htslib-1.9/bin/tabix"
bgzip="/xx/xxx/bin/htslib-1.9/bin/bgzip"
gatk="/xx/xxx/bin/gatk-4.3.0.0/gatk"
plink2="/xx/xxx/bin/plink2"

#######################################
# TUM GSA 2020/2022 select indidviduals 
#######################################

mkdir "${wd}/genotype_data/TUM"

tum_samples="${wd}/../TUM_genIDs_new.txt"

parallel -q $bcftools view -S $tum_samples --force-samples -e 'INFO<0.7' -e 'MAF<0.001' -O z -o ${wd}/genotype_data/TUM/GSA_{1}.vcf.gz /data_featherstone/shared_data/genetics/GSA_2020-10_2022_Combined/06_VCF/TUM_GSA_2020-10_2022_Combined_sampleqc_varqc_refalign_phased_imputed_postimpqc_chr{1}.vcf.gz ::: {1..22}
parallel -q $tabix ${wd}/genotype_data/TUM/GSA_{1}.vcf.gz ::: {1..22}
ls ${wd}/genotype_data/TUM/GSA*.vcf.gz > ${wd}/genotype_data/TUM/GSA.list
$bcftools concat -f ${wd}/genotype_data/TUM/GSA.list -O z -o ${wd}/genotype_data/TUM/GSA.vcf.gz
rm ${wd}/genotype_data/TUM/GSA.list
$bcftools stats ${wd}/genotype_data/TUM/GSA.vcf.gz > ${wd}/genotype_data/TUM/GSA.stats
$bcftools annotate -x "^FORMAT/GP,INFO" -O z -o ${wd}/genotype_data/TUM/GSA_formatted.vcf.gz ${wd}/genotype_data/TUM/GSA.vcf.gz
$tabix ${wd}/genotype_data/TUM/GSA_formatted.vcf.gz
echo "1 chr1" > ${wd}/genotype_data/TUM/change_chromosome_names.txt
for i in {2..22}; do 
	echo "${i} chr${i}" >> ${wd}/genotype_data/TUM/change_chromosome_names.txt
done
$bcftools annotate --rename-chrs ${wd}/genotype_data/TUM/change_chromosome_names.txt -O z -o ${wd}/genotype_data/TUM/GSA_formatted_renamed.vcf.gz ${wd}/genotype_data/TUM/GSA_formatted.vcf.gz
$tabix ${wd}/genotype_data/TUM/GSA_formatted_renamed.vcf.gz

#######################################
# TUM GSA 2023 select indidviduals 
#######################################

tum_samples="${wd}/../../TUM_GSA2023genIDs.txt"

parallel -q $bcftools view -S $tum_samples --force-samples -e 'INFO<0.7' -e 'MAF<0.001' -O z -o ${wd}/genotype_data/TUM/GSA2023_{1}.vcf.gz /data_featherstone/shared_data/genetics/2023-GSA/06_VCF/TUM_GSA_2023_sampleqc_varqc_refalign_phased_imputed_postimpqc_chr{1}.vcf.gz ::: {1..22}
parallel -q $tabix ${wd}/genotype_data/TUM/GSA2023_{1}.vcf.gz ::: {1..22}
ls ${wd}/genotype_data/TUM/GSA2023*.vcf.gz > ${wd}/genotype_data/TUM/GSA2023.list
$bcftools concat -f ${wd}/genotype_data/TUM/GSA2023.list -O z -o ${wd}/genotype_data/TUM/GSA2023.vcf.gz
rm ${wd}/genotype_data/TUM/GSA2023.list
$bcftools stats ${wd}/genotype_data/TUM/GSA2023.vcf.gz > ${wd}/genotype_data/TUM/GSA2023.stats
$bcftools annotate -x "^FORMAT/GP,INFO" -O z -o ${wd}/genotype_data/TUM/GSA2023_formatted.vcf.gz ${wd}/genotype_data/TUM/GSA2023.vcf.gz
$tabix ${wd}/genotype_data/TUM/GSA2023_formatted.vcf.gz
echo "1 chr1" > ${wd}/genotype_data/TUM/change_chromosome_names.txt
for i in {2..22}; do 
	echo "${i} chr${i}" >> ${wd}/genotype_data/TUM/change_chromosome_names.txt
done
$bcftools annotate --rename-chrs ${wd}/genotype_data/TUM/change_chromosome_names.txt -O z -o ${wd}/genotype_data/TUM/GSA2023_formatted_renamed.vcf.gz ${wd}/genotype_data/TUM/GSA2023_formatted.vcf.gz
$tabix ${wd}/genotype_data/TUM/GSA2023_formatted_renamed.vcf.gz

######################################
# Merge both TUM datasets
######################################

mkdir "${wd}/genotype_data/TUM/shared_variants"

$bcftools isec -n +2 -p ${wd}/genotype_data/TUM/shared_variants ${wd}/genotype_data/TUM/GSA_formatted_renamed.vcf.gz ${wd}/genotype_data/TUM/GSA2023_formatted_renamed.vcf.gz
$bgzip -c ${wd}/genotype_data/TUM/shared_variants/0000.vcf > ${wd}/genotype_data/TUM/shared_variants/0000.vcf.gz
$bgzip -c ${wd}/genotype_data/TUM/shared_variants/0001.vcf > ${wd}/genotype_data/TUM/shared_variants/0001.vcf.gz
$tabix ${wd}/genotype_data/TUM/shared_variants/0000.vcf.gz
$tabix ${wd}/genotype_data/TUM/shared_variants/0001.vcf.gz

$bcftools merge ${wd}/genotype_data/TUM/shared_variants/0000.vcf.gz ${wd}/genotype_data/TUM/shared_variants/0001.vcf.gz -O z -o ${wd}/genotype_data/TUM/TUM_merged.vcf.gz
$tabix ${wd}/genotype_data/TUM/TUM_merged.vcf.gz
$bcftools stats ${wd}/genotype_data/TUM/TUM_merged.vcf.gz > ${wd}/genotype_data/TUM/TUM_merged.stats

######################################
# Liftover to hg83
######################################

# Download hg38.fa.gz from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/ 
# Download hg19ToHg38.over.chain.gz from #http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
# Download resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta from https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
# Download resources-broad-hg38-v0-Homo_sapiens_assembly38.dict from https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
cp ${wd}/genotype_data/TUM/TUM_merged.vcf.gz ${wd}/genotype_data/TUM/TUM.vcf.gz 
cp ${wd}/genotype_data/TUM/TUM_merged.vcf.gz.tbi ${wd}/genotype_data/TUM/TUM.vcf.gz.tbi
gunzip ${wd}/genotype_data/TUM/TUM.vcf.gz 
sed -i 's/Type=Float,Number=G/Number=G,Type=Float/' ${wd}/genotype_data/TUM/TUM.vcf 
$bgzip ${wd}/genotype_data/TUM/TUM.vcf
$tabix ${wd}/genotype_data/TUM/TUM.vcf.gz
$gatk LiftoverVcf -I ${wd}/genotype_data/TUM/TUM.vcf.gz -O ${wd}/genotype_data/TUM/TUM_hg38.vcf.gz -CHAIN ${wd}/genotype_data/TUM/hg19ToHg38.over.chain -REJECT ${wd}/genotype_data/TUM/GSA_2020_hg38_rejected_variants.vcf -REFERENCE_SEQUENCE ${wd}/genotype_data/TUM/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta 
#-RECOVER_SWAPPED_REF_ALT

######################################
# CAM datasets
######################################

cam_samples="${wd}/../CAM_genIDs.txt"

$bcftools view -S $cam_samples --force-samples -e 'R2<0.7' -q 0.001:minor -O z -o ${wd}/genotype_data/CAM/CAM.vcf.gz /data_featherstone/christiane/SC/CAM_TUM/sc_for_tum/merged_vcf.vcf.gz
$bcftools annotate -x "^FORMAT/GP,INFO" -O z -o ${wd}/genotype_data/CAM/CAM_formatted.vcf.gz ${wd}/genotype_data/CAM/CAM.vcf.gz

$tabix ${wd}/genotype_data/CAM/CAM_formatted.vcf.gz

######################################
# Merge TUM and CAM data
######################################

$bcftools isec -n +2 -p ${wd}/genotype_data/shared_variants ${wd}/genotype_data/TUM/TUM_hg38.vcf.gz ${wd}/genotype_data/CAM/CAM_formatted.vcf.gz
$bgzip -c ${wd}/genotype_data/shared_variants/0000.vcf > ${wd}/genotype_data/shared_variants/0000.vcf.gz
$bgzip -c ${wd}/genotype_data/shared_variants/0001.vcf > ${wd}/genotype_data/shared_variants/0001.vcf.gz
$tabix ${wd}/genotype_data/shared_variants/0000.vcf.gz
$tabix ${wd}/genotype_data/shared_variants/0001.vcf.gz

$bcftools merge ${wd}/genotype_data/shared_variants/0000.vcf.gz ${wd}/genotype_data/shared_variants/0001.vcf.gz -O z -o ${wd}/genotype_data/merged.vcf.gz
$tabix ${wd}/genotype_data/merged.vcf.gz
$bcftools stats ${wd}/genotype_data/merged.vcf.gz > ${wd}/genotype_data/merged.stats

######################################
# Get info about strand issues
######################################

gunzip -c ${wd}/genotype_data/CAM/CAM.vcf.gz >  ${wd}/genotype_data/CAM/CAM.vcf
gunzip -c ${wd}/genotype_data/TUM/TUM_hg38.vcf.gz >  ${wd}/genotype_data/TUM/TUM_hg38.vcf

awk '{ if (!/^#/) {print $1"\t"$2"\t"$4"\t"$5}}' ${wd}/genotype_data/CAM/CAM.vcf > ${wd}/genotype_data/shared_variants/CAM.variants
awk '{ if (!/^#/) {print $1"\t"$2"\t"$4"\t"$5}}' ${wd}/genotype_data/TUM/TUM_hg38.vcf > ${wd}/genotype_data/shared_variants/TUM_hg38.variants
Rscript ${wd}/Scripts/check_strands.R ${wd}/genotype_data/shared_variants/CAM.variants ${wd}/genotype_data/shared_variants/TUM_hg38.variants ${wd}/genotype_data/shared_variants/sites.txt ${wd}/genotype_data/shared_variants/strand_issues.txt


######################################
# Convert to plink format
######################################

$plink2 --vcf ${wd}/genotype_data/merged.vcf.gz dosage="GP" --make-pgen --out ${wd}/genotype_data/merged
$plink2 --pfile ${wd}/genotype_data/merged --chr 1-22 --make-pgen --out ${wd}/genotype_data/merged_1to22

######################################
# Get variants with AF that differs from the reference (1000G) MAF
######################################

$plink2 --vcf ${wd}/genotype_data/CAM/CAM_formatted.vcf.gz dosage="GP" --make-pgen --out testCAM
$plink2 --pfile testCAM --set-all-var-ids @:#:$\r:$\a --new-id-max-allele-len 1000 --make-pgen --out testCAM2
$plink2 --pfile testCAM2 --freq --out testCAM2
$plink2 --vcf ${wd}/genotype_data/TUM/TUM_hg38.vcf.gz dosage="GP" --allow-extra-chr --make-pgen --out testTUM
$plink2 --pfile testTUM --allow-extra-chr --set-all-var-ids @:#:$\r:$\a --new-id-max-allele-len 1000 --make-pgen --out testTUM2
$plink2 --pfile testTUM2 --allow-extra-chr --freq --out testTUM2
awk '{ if (!/^#/) {print $3}}' testTUM2.pvar > testsnps.list
$plink2 --vcf /data_featherstone/shared_data/references/1000G_hg38/1000G_hg38.vcf.gz --allow-extra-chr --make-pgen --out test1000G
$plink2 --pfile test1000G --allow-extra-chr --set-all-var-ids @:#:$\r:$\a --new-id-max-allele-len 1000 --make-pgen --out test1000G2
awk '{ if ($3 == "EUR") {print $1}}' /data_featherstone/shared_data/references/1000Genomes/more/integrated_call_samples_v3.20130502.ALL.panel > test1000G_inds.txt
$plink2 --pfile test1000G2 --keep test1000G_inds.txt --make-pgen --out test1000G2
$plink2 --pfile test1000G2 --extract testsnps.list --allow-extra-chr --freq --out test1000G2
$plink2 --pfile ${wd}/genotype_data/merged_1to22 --freq --out ${wd}/genotype_data/merged_1to22
Rscript ${wd}/Scripts/check_MAF.R testCAM2.afreq testTUM2.afreq test1000G2.afreq ${wd}/genotype_data/merged_1to22.afreq ${wd}/genotype_data/merged_keep_MAF.txt ${wd}/genotype_data/merged_MAF_beforeQC.png ${wd}/genotype_data/merged_MAF_afterQC.png
$plink2 --pfile ${wd}/genotype_data/merged_1to22 --extract ${wd}/genotype_data/merged_keep_MAF.txt --make-pgen --out ${wd}/genotype_data/merged_1to22
rm test*

######################################
# Save
######################################

mkdir ${wd}/genotype_data/merged
mv ${wd}/genotype_data/merged*  ${wd}/genotype_data/merged/
mv ${wd}/genotype_data/shared_variants ${wd}/genotype_data/merged/


######################################
# Population plots (with 1000G)
######################################

mkdir mkdir ${wd}/genotype_data/PCA
# Prepare data
ref1000g="/data_featherstone/shared_data/references/1000G_hg38/1000G_hg38.vcf.gz"
#$tabix ${ref1000g}
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22 --make-bed --out ${wd}/genotype_data/merged/merged_1to22
$plink2 --bfile ${wd}/genotype_data/merged/merged_1to22 --export vcf --out ${wd}/genotype_data/merged/merged_1to22_gt
$bgzip -c ${wd}/genotype_data/merged/merged_1to22_gt.vcf > ${wd}/genotype_data/merged/merged_1to22_gt.vcf.gz
$tabix ${wd}/genotype_data/merged/merged_1to22_gt.vcf.gz

#Merge 
$bcftools isec -n +2 -p ${wd}/genotype_data/PCA/shared_variants ${wd}/genotype_data/merged/merged_1to22_gt.vcf.gz $ref1000g
$bgzip -c ${wd}/genotype_data/PCA/shared_variants/0000.vcf > ${wd}/genotype_data/PCA/shared_variants/0000.vcf.gz
$bgzip -c ${wd}/genotype_data/PCA/shared_variants/0001.vcf > ${wd}/genotype_data/PCA/shared_variants/0001.vcf.gz
$tabix ${wd}/genotype_data/PCA/shared_variants/0000.vcf.gz
$tabix ${wd}/genotype_data/PCA/shared_variants/0001.vcf.gz
$bcftools merge ${wd}/genotype_data/PCA/shared_variants/0000.vcf.gz ${wd}/genotype_data/PCA/shared_variants/0001.vcf.gz -O z -o ${wd}/genotype_data/PCA/CAMTUM_1000g_hg38.vcf.gz
$tabix ${wd}/genotype_data/PCA/CAMTUM_1000g_hg38.vcf.gz
awk '{ if ($3 == "EUR") {print $1}}' /data_featherstone/shared_data/references/1000Genomes/more/integrated_call_samples_v3.20130502.ALL.panel > ${wd}/genotype_data/PCA/EUR_samples.txt
awk '{ if (!/^#/){ print $1 }}' ${wd}/genotype_data/merged/merged_1to22.psam > ${wd}/genotype_data/PCA/merged_samples.txt
cat ${wd}/genotype_data/PCA/merged_samples.txt ${wd}/genotype_data/PCA/EUR_samples.txt > ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_samples.txt
$bcftools view -S ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_samples.txt --force-samples -O z -o ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38.vcf.gz ${wd}/genotype_data/PCA/CAMTUM_1000g_hg38.vcf.gz
cp ${wd}/genotype_data/PCA/CAMTUM_1000g_hg38.vcf.gz ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38.vcf.gz


# To plink format
$plink2 --vcf ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38.vcf.gz --make-pgen --out ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38
$plink2 --vcf ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38.vcf.gz --make-pgen --out ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38

# Calculate PCs
$plink2 --pfile ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38 --missing --out ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38_miss
$plink2 --pfile ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38 --exclude range ${wd}/genotype_data/PCA/remove_prune.txt --maf 0.05 --hwe 1e-3 --geno 0.05 --indep-pairwise 200 100 0.2 --out ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38_prunelist
$plink2 --pfile ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38 --extract ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38_prunelist.prune.in --pca --out ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38_pca
$plink2 --pfile ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38 --missing --out ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38_miss
$plink2 --pfile ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38 --exclude range ${wd}/genotype_data/PCA/remove_prune.txt --maf 0.05 --hwe 1e-3 --geno 0.05 --indep-pairwise 200 100 0.2 --out ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38_prunelist
$plink2 --pfile ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38 --extract ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38_prunelist.prune.in --pca --out ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38_pca

# PCA plots
cohort_info_file="/data_featherstone/shared_data/references/1000Genomes/integrated_call_samples_v3.20130502.ALL.panel"
Rscript ${wd}/Scripts/pca_plots_1000G.R ${wd}/genotype_data/PCA/CAMTUM_1000g_EUR_hg38_pca.eigenvec $cohort_info_file EUR ${wd}/genotype_data/PCA/PCplot_EUR_beforeQC.png
Rscript ${wd}/Scripts/pca_plots_1000G.R ${wd}/genotype_data/PCA/CAMTUM_1000g_ALL_hg38_pca.eigenvec $cohort_info_file ALL ${wd}/genotype_data/PCA/PCplot_ALL_beforeQC.png
