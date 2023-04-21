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
