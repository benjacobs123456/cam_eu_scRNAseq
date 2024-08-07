!/bin/bash


wd="/xx/xxx/SC/CAM_TUM/eqtl/Targeted/genotype_data"
plink="/xx/xxx/bin/plink"
plink2="/xx/xxx/bin/plink2"

dataset0="${wd}/merged/merged_1to22"

mkdir ${wd}/QC

dataset="${wd}/QC/merged_1to22"

# Basic SNP QC before doing Sample QC
$plink2 --pfile ${dataset0} --freq --out ${dataset}_qc1
awk '{ if (!/^#/ && $5 < 0.001) {print $2}}' ${dataset}_qc1.afreq > ${dataset}_qc1_freqfilter.txt 
$plink2 --pfile ${dataset0} --geno 0.02 dosage --exclude ${dataset}_qc1_freqfilter.txt --make-pgen --out ${dataset}_qc1

# Missingness
$plink2 --pfile ${dataset}_qc1 --mind 0.02 dosage -make-pgen --out ${dataset}_qc1_miss
comm -23 <(sort ${dataset0}.psam) <(sort ${dataset}_qc1_miss.psam) > "${wd}/QC/remove_missing1.txt"
awk '{print $1}' "${wd}/QC/remove_missing1.txt" > "${wd}/QC/remove_missing.txt"
$plink2 --pfile ${dataset0} --remove "${wd}/QC/remove_missing.txt" --make-pgen --out ${dataset}_noMiss

# QC Heterozygosity
$plink2 --pfile ${dataset}_noMiss --freq --out ${dataset}_noMiss
awk '{ if (!/^#/ && $5 < 0.05) {print $2}}' ${dataset}_noMiss.afreq > ${dataset}_noMiss_freqfilter.txt 
$plink2 --pfile ${dataset}_noMiss --exclude ${dataset}_noMiss_freqfilter.txt --hwe 1e-6 --geno 0.02 dosage --indep-pairwise 200 100 0.2 --make-bed --out ${dataset}_noMiss_prunedHET
$plink2 --pfile ${dataset}_noMiss --extract ${dataset}_noMiss_prunedHET.prune.in --het --out ${dataset}_noMiss_het
Rscript ${wd}/../Scripts/QC_samples_het.R ${dataset}_noMiss_het.het "${wd}/QC/remove_heterozygosity.txt"
$plink2 --pfile ${dataset}_noMiss --remove  "${wd}/QC/remove_heterozygosity.txt" --make-pgen --out ${dataset}_noMiss_noHet

# QC Duplicate and Relatives
$plink2 --pfile ${dataset}_noMiss_noHet  --freq --out ${dataset}_noMiss_noHet 
awk '{ if (!/^#/ && $5 < 0.05) {print $2}}' ${dataset}_noMiss_noHet.afreq > ${dataset}_noMiss_noHet_freqfilter.txt 
$plink2 --pfile ${dataset}_noMiss_noHet --geno 0.02 dosage --hwe 1e-6 --exclude ${dataset}_noMiss_noHet_freqfilter.txt --indep-pairwise 200 100 0.2 --out ${dataset}_noMiss_noHet_prunedREL
$plink2 --pfile ${dataset}_noMiss_noHet --extract ${dataset}_noMiss_noHet_prunedREL.prune.in --king-cutoff 0.125 --out ${dataset}_noMiss_noHet
$plink2 --pfile ${dataset}_noMiss_noHet --remove ${dataset}_noMiss_noHet.king.cutoff.out.id --make-pgen --out ${dataset}_noMiss_noHet_noRel

# QC Population outliers
$plink2 --pfile ${dataset}_noMiss_noHet_noRel --freq --out ${dataset}_noMiss_noHet_noRel
awk '{ if (!/^#/ && $5 < 0.05) {print $2}}' ${dataset}_noMiss_noHet_noRel.afreq > ${dataset}_noMiss_noHet_noRel_freqfilter.txt 
$plink2 --pfile ${dataset}_noMiss_noHet_noRel --exclude range ${wd}/QC/remove_prune.txt --make-pgen --out ${dataset}_noMiss_noHet_noRel1 
$plink2 --pfile ${dataset}_noMiss_noHet_noRel1 --geno 0.02 dosage --hwe 1e-6 --exclude ${dataset}_noMiss_noHet_noRel_freqfilter.txt --indep-pairwise 200 100 0.2 --out ${dataset}_noMiss_noHet_noRel_prunedMDS
$plink2 --pfile ${dataset}_noMiss_noHet_noRel --extract ${dataset}_noMiss_noHet_noRel_prunedMDS.prune.in --pca --out ${dataset}_noMiss_noHet_noRel_pca
Rscript ${wd}/../Scripts/QC_samples_pop.R ${dataset}_noMiss_noHet_noRel_pca.eigenvec ${dataset}_noMiss_noHet_noRel.psam ${wd}/QC/remove_pop ${wd}/QC/MDSplot_beforeQC.png
$plink2 --pfile ${dataset}_noMiss_noHet_noRel --remove  ${wd}/QC/remove_pop1.txt --make-pgen --out ${dataset}_noMiss_noHet_noRel_noOutl
$plink2 --pfile ${dataset}_noMiss_noHet_noRel_noOutl --remove  ${wd}/QC/remove_pop2.txt --make-pgen --out ${dataset}_noMiss_noHet_noRel_noOutl
$plink2 --pfile ${dataset}_noMiss_noHet_noRel_noOutl --remove  ${wd}/QC/remove_pop3.txt --make-pgen --out ${dataset}_noMiss_noHet_noRel_noOutl

# QCd individuals
awk '{print $1}' ${dataset}_noMiss_noHet_noRel_noOutl.psam > "${wd}/QC/QCd_individuals.txt"
$plink2 --pfile ${dataset0} --keep "${wd}/QC/QCd_individuals.txt" --make-pgen --out ${dataset0}_SampleQC

