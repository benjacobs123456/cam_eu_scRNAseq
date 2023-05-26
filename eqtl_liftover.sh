# liftover MS GWAS stats
cd /rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/eqtl/

awk 'NR>1{print "chr"$1,$2-1,$2,$3}' /rds/user/hpcjaco1/hpc-work/discovery_metav3.0.meta > ms_gwas_bed_hg19.bed

/home/hpcjaco1/liftover/liftOver ms_gwas_bed_hg19.bed /home/hpcjaco1/liftover/hg19ToHg38.over.chain.gz ms_gwas_hg38.bed unlifted.bed

# liftover eQTLgen
cd /rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/eqtl/
zcat 2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | awk 'NR>1{print "chr"$3,$4-1,$4,$2}'  > eqtlgen_hg19.bed

/home/hpcjaco1/liftover/liftOver eqtlgen_hg19.bed /home/hpcjaco1/liftover/hg19ToHg38.over.chain.gz eqtlgen_hg38.bed unlifted.bed
