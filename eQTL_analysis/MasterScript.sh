# Define working directory and releveant files/executables
wd="/xx/xxx/SC/CAM_TUM/eqtl/Targeted"
plink="/xx/xxx/bin/plink"
plink2="/xx/xxx/bin/plink2"
pheno_file="/xx/xxx/SC/Seurat/phenotypes_all.csv"

# Create folders
mkdir ${wd}/CSF/ 
mkdir ${wd}/CSF/MS
mkdir ${wd}/CSF/Control
mkdir ${wd}/CSF/All
mkdir ${wd}/PBMC/ 
mkdir ${wd}/PBMC/MS
mkdir ${wd}/PBMC/Control
mkdir ${wd}/PBMC/All

### Prepare gene expression data
Rscript $wd/Scripts/prepare_expression_data.R ${wd}/CSF/ ${wd}/all_combo_with_updated_pheno.rds CSF 
Rscript $wd/Scripts/prepare_expression_data.R ${wd}/PBMC/ ${wd}/all_combo_with_updated_pheno.rds PBMC 
Rscript $wd/Scripts/cell_numbers.R 

### Calculate PCs
Rscript $wd/Scripts/gex_pca.R ${wd}/CSF/ ${wd}/CSF/all_combo_renormalizedperbatch.RDS CSF
Rscript $wd/Scripts/gex_pca.R ${wd}/CSF/ ${wd}/PBMC/all_combo_renormalizedperbatch.RDS PBMC 

### Prepare information on gene positions 
awk -F"\t" '{ if (!/^#/ && $3 == "gene"){gsub(/chr/,"",$1); split($9,a,";"); $9=a[1]; split($9,b," "); $9=b[2]; $10=a[4]; split($10,c," "); $10=c[2]; gsub(/\"/,"",$9); gsub(/\"/,"",$10); $11=$4-1000000; if ($11 < 0) {$11 = 0}; $12=$4+1000000; print $1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11"\t"$12}}' /xx/xxx/SC/Software/refdata-gex-GRCh38-2020-A/genes/genes.gtf > ${wd}/gene_list1.txt
echo -e "chromosome_name\tstart_position\tend_position\tensembl_gene_id\thgnc_symbol\tStarteQTLWindow\tEndeQTLWindow" > ${wd}/header.txt
cat ${wd}/header.txt ${wd}/gene_list1.txt > ${wd}/gene_list.txt
rm ${wd}/header.txt; rm ${wd}/gene_list1.txt

### Prepare genetic data (cohort selection and merging)
bash $wd/Scripts/combine_genotype_data.sh

### Prepare genetic data - Sample QC
bash $wd/Scripts/Sample_QC.sh

### Prepare genetic data - Genotype QC
awk '{ if (!/^#/ && (length ($4) > 1 || length ($5) > 1)) {print $3}}' ${wd}/genotype_data/QC/merged_1to22_noMiss_noHet_noRel_noOutl.pvar > ${wd}/genotype_data/QC/merged_1to22_noMiss_noHet_noRel_noOutl_indels.txt
$plink2 --pfile  ${wd}/genotype_data/merged/merged_1to22_SampleQC --freq --out  ${wd}/genotype_data/merged/merged_1to22_SampleQC
awk '{ if (!/^#/ && $5 < 0.01) {print $2}}'  ${wd}/genotype_data/merged/merged_1to22_SampleQC.afreq >  ${wd}/genotype_data/merged/merged_1to22_SampleQC_freqfilter.txt 
cat ${wd}/genotype_data/merged/merged_1to22_SampleQC_freqfilter.txt ${wd}/genotype_data/QC/merged_1to22_noMiss_noHet_noRel_noOutl_indels.txt > ${wd}/genotype_data/QC/filterSNPs.txt
rm ${wd}/genotype_data/merged/merged_1to22_SampleQC_freqfilter.txt
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC --exclude ${wd}/genotype_data/QC/filterSNPs.txt --hwe 1e-3 --geno 0.02 dosage --make-pgen --out ${wd}/genotype_data/QC/merged_1to22_SampleQC_VariantQC
cp ${wd}/genotype_data/QC/*VariantQC* ${wd}/genotype_data/merged/

### Get list of ambiguous SNPs
awk '{ if (!/^#/) { if (($4 == "T" && $5 == "A") || ($4 == "A" && $5 == "T") || ($4 == "C" && $5 == "G") || ($4 == "G" && $5 == "C")) {print$3}}}' ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC.pvar > ${wd}/genotype_data/merged/ambiguousSNPs.txt

### Prepare covariates
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC --freq --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC 
awk '{ if (!/^#/ && $5 < 0.05) {print $2}}' ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC.afreq > ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_freqfilter.txt 
rm ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC.afreq
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC --exclude range ${wd}/genotype_data/QC/remove_prune.txt --make-pgen --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC1
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC1 --exclude ${wd}/genotype_data/merged/ambiguousSNPs.txt --make-pgen --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC2
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC2 --geno 0.02 dosage --hwe 1e-3 --exclude ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_freqfilter.txt --indep-pairwise 200 100 0.2 --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_pruned
$plink2 --pfile ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC --extract ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_pruned.prune.in --pca --out ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_pca
rm ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC1*
rm ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_pruned*
#rm ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC*filter*
Rscript $wd/Scripts/prepare_covars.R ${wd}/genotype_data/merged/merged_1to22_SampleQC_VariantQC_pca.eigenvec ${wd}/../all_combo_with_updated_pheno_metadata.csv ${wd}/covariates.txt

### PCA plots
Rscript $wd/Scripts/pca_plots.R ${wd}/covariates.txt ${wd}/pca_plots.png
bash $wd/Scripts/pca_plots_1000G.sh

### Run first eQTL analysis
parallel -j 50 Rscript $wd/Scripts/eQTL_analysis_plink.R {1} {2} {3} {4} ::: "B.cells" "Tregs" "Plasma.cells" "CD4.T.cells" "CD8.T.cells" "NK.cells" "mDCs" "pDCs" "CD14.Mono" "CD16.Mono" "MAIT.cells" ::: {1..22} ::: All MS Control ::: PBMC CSF

# Combine results
for ct in "B.cells" "Tregs" "Plasma.cells" "CD4.T.cells" "CD8.T.cells" "NK.cells" "mDCs" "pDCs" "CD14.Mono" "CD16.Mono" "MAIT.cells"; do
for cohort in "All" "MS" "Control"; do
for source in "CSF" "PBMC"; do
cat ${wd}/${source}/${cohort}/${ct}_chr*results_plink.tsv > ${wd}/${source}/${cohort}/${ct}_results1.tsv
done
done
done

# Combine all results
Rscript $wd/Scripts/combine_results_Targeted.R $wd All
Rscript $wd/Scripts/combine_results_Targeted.R $wd MS
Rscript $wd/Scripts/combine_results_Targeted.R $wd Control
#bash $wd/Scripts/start_combine_results_All.sh

# Get relevant eQTL results
Rscript $wd/Scripts/select_releveant_results.R 

# Find previously described / "novel" eQTLs
Rscript $wd/Scripts/novel_eqtls.R 

# Permutation on "novel" findings
mkdir $wd/Permutation
#parallel -j 65 Rscript $wd/Scripts/Permutation.R ${1} ::: {293..1} 
parallel -j 65 bash $wd/Scripts/start_Permutation.R ${1} ::: {293..1} 
cat ${wd}/Permutation/perm* > ${wd}/permutation_results.txt
Rscript $wd/Scripts/Add_permutation_pvalues.R 

# Find overlap with MS susceptibility
Rscript $wd/Scripts/MSGWAS_overlap.R 

# Find source specific eQTLs
Rscript $wd/Scripts/CSFspecificity.R 

# Find cell type specific eQTLs
Rscript $wd/Scripts/cellspecificity.R 

# Some figures
Rscript $wd/Scripts/single_plots.R rs4810801 PREX1 B.cells CD4.T.cells CD8.T.cells CSF $wd
#Rscript $wd/Scripts/single_plots.R rs61909096 ETS1 CD4.T.cells B.cells CD8.T.cells CSF $wd
Rscript $wd/Scripts/plot_singleGene_CSFvsPBMC.R #for PREX1
Rscript $wd/Scripts/Submission_Figures.R 




