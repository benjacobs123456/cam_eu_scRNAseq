# Single-cell RNA seq analysis of Multiple Sclerosis CSF & blood

- This directory contains scripts used to analyse scRNAseq data from the Cambridge and TUM cohorts.
- All analyses were run on the Cambridge Slurm HPC.
- Unless specified, all scripts were run in R/4.0.3.

# Deconvolution and basic qc
````unix
Rscript gex_deconvolution_qc_step1.R
````

This script performs the following QC steps on each batch of GEX data
- Filtering by RNA count & MT%
- Ambient RNA correction with SoupX
- Doublet identification
- Per-batch normalisation with SCTransform

# Integration
````unix
Rscript integration.R
````
This script integrates data across batches using Harmony.

# UMAP
````unix
Rscript umap.R
````
This script runs UMAP and Louvain clustering across a range of paramaters.

# Cluster biomarkers
````unix
Rscript cluster_biomarkers.R
````

This script does the following
- Makes some plots exploring clustering  
- Compares annotations of cell types across different methods (SingleR, Azimuth, Celltypist)
- Calculates cluster-specific biomarkers

# DE & DA
To run DE and DA using the broad clusters:
````unix
Rscript de_da_tests.R
````

# GSEA
Then to run GSEA on the broad clusters with those DE results:
````unix
Rscript gsea.R
````

# Prepare for Immcantation/Dandelion
These scripts prepare the VDJ data for QC with dandelion
````unix
Rscript dandelion_preparation.R TCR
Rscript dandelion_preparation.R BCR
````

And then to run dandelion:
````unix
sbatch dandelion_preprocess_bcr.sh
sbatch dandelion_preprocess_tcr.sh
````

# Filter GEX based on VDJ QC
These scripts then filter the GEX data based on the QC'd VDJ data:
````unix
sbatch dandelion_filtering_tcr.sh
sbatch dandelion_filtering.sh
````

# Analyse VDJ data
````unix
Rscript bcr_analysis_28_06_22.R
Rscript tcr_clonality_exploration.R
````

# Effect of DRB1
````unix 
Rscript eqtl.R
````
