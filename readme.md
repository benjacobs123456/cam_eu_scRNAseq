# Single-cell RNA seq analysis of Multiple Sclerosis CSF & blood

- This directory contains scripts used to analyse scRNAseq data from the Cambridge and TUM cohorts.
- All analyses were run on the Cambridge Slurm HPC.
- Unless specified, all scripts were run in R/4.0.3.

# Deconvolution and basic qc
````unix
sbatch gex_deconvolution_qc_step1.sh
````
This script runs `Rscript gex_deconvolution_qc.R` which performs the following QC steps on each batch of GEX data
- Filtering by RNA count & MT%
- Ambient RNA correction with SoupX
- Doublet identification
- Per-batch normalisation with SCTransform

# Integration
````unix
sbatch integration_icelake.sh
````
This script runs `Rscript integration.R` which integrates the post-qc datasets using Harmony.

# UMAP
````unix
sbatch umap_icelake.sh
````
This script runs  `Rscript umap.R` which performs UMAP and Louvain clustering across a range of parameters.

# Cluster biomarkers
````unix
Rscript cluster_biomarkers_icelake.sh
````

This script runs `Rscript find_cluster_biomarkers_and_update_pheno.R` which does the following
- Cleans phenotypes
- Makes some plots exploring clustering  
- Compares annotations of cell types across different methods (SingleR, Azimuth, Celltypist)
- Calculates cluster-specific biomarkers

# CellTypist annotation
````unix
sbatch celltypist.sh
````
- Which runs celltypist in each cluster

# Update cluster IDs
`````unix
sbatch update_clusters.sh
````
- Which runs `Rscript update_cluster_labels.R` to update cluster IDs

# DE & DA
To run DE and DA using the broad clusters:
````unix
sbatch de_icelake.sh
````
- Runs DE tests with `Rscript de_summary_plots.R`
- Summary plots then made with `Rscript  de_da_tests_phenotypes.R`


# GSEA
Then to run GSEA on the broad clusters with those DE results:
````unix
sbatch gsea.sh
````

# Pathway analysis
````unix
sbatch pathway_analysis.sh
````

# Cell-cell communication
````unix
sbatch ccc_per_sample.sh
````
- Which runs LIANA on a per-sample basis
- `Rscript Rscript ccc_overall.R` then combines and explores these results



# Prepare for Immcantation/Dandelion
These scripts prepare the VDJ data for QC with dandelion
````unix
Rscript dandelion_preparation.R TCR 
Rscript dandelion_preparation.R BCR 
Rscript dandelion_preparation2.R TCR 
Rscript dandelion_preparation2.R BCR 
Rscript make_dandelion_metafile.R
````

And then to run dandelion pre-processing:
````unix
sbatch dandelion_bcr.sh
sbatch dandelion_tcr.sh
````

# Filter GEX based on VDJ QC
These scripts then filter the GEX data based on the QC'd VDJ data:
````unix
sbatch dandelion_filtering_tcr.sh
sbatch dandelion_filtering.sh
````

# Analyse VDJ data
````unix
Rscript bcr_analysis.R
Rscript tcr_analysis.R
````

# eQTL analysis
````unix
bash eqtl_prep.sh
sbatch eqtl_split_cell_types_icelake.sh
sbatch eqtl_icelake_all_analyses.sh

# read in data
Rscript read_eqtl_results.R

# Basic plots
Rscript eqtl_basic_plots.R

# run coloc
Rscript eqtl_coloc_020223.R

# Liftover prep for coloc
Rscript eqtl_coloc_plots.R

Rscript eqtl_cell_type_specific.R
````

