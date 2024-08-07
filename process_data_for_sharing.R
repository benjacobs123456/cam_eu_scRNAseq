library(Seurat)
library(tidyverse)

# read in data
all_combo = readRDS( file = "/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/all_combo_with_updated_pheno.rds")
DefaultAssay(all_combo) = "RNA"
all_combo[["SCT"]] =NULL
all_combo[["umap"]] = NULL

# select relevant columns in metadata
codex = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/pat_codex.csv")
all_combo@meta.data = all_combo@meta.data %>%
  left_join(codex,by="iid") %>%
  dplyr::select(nCount_RNA, nFeature_RNA, batch_id, contains("percent."),source,processing_site, Category,Sex,OCB,fully_anonymous_pseudoid, contains("ann_celltypist"))

rownames(all_combo@meta.data)  = colnames(all_combo)

# save dataset
saveRDS(all_combo,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_gex_data_for_sharing.rds")
write_csv(all_combo@meta.data,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_gex_data_for_sharing_metadata.csv")

# repeat for b cells
# read in data
all_combo = readRDS( file = "/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/BCR/b_cells_post_processing.rds")
DefaultAssay(all_combo) = "RNA"
all_combo[["SCT"]] =NULL
all_combo[["umap"]] = NULL

# select relevant columns in metadata
codex = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/pat_codex.csv")
all_combo@meta.data = all_combo@meta.data %>%
  left_join(codex,by="iid")

all_combo@meta.data = all_combo@meta.data %>%
  dplyr::select(nCount_RNA, batch_id, nFeature_RNA, contains("percent."),source,processing_site, Category,Sex,OCB,fully_anonymous_pseudoid, contains("ann_celltypist"),c(202:264),-cell_id,-sample_id)

rownames(all_combo@meta.data)  = colnames(all_combo)

# save dataset
saveRDS(all_combo,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_bcell_data_for_sharing.rds")
write_csv(all_combo@meta.data,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_bcell_data_for_sharing_metadata.csv")
dat = readRDS("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_bcell_data_for_sharing.rds")

# repeat for t cells
# read in data
all_combo = readRDS( file = "/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/TCR/t_cells_post_processing.rds")
DefaultAssay(all_combo) = "RNA"
all_combo[["SCT"]] =NULL
all_combo[["umap"]] = NULL

# select relevant columns in metadata
codex = read_csv("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/datasets/pat_codex.csv")
all_combo@meta.data = all_combo@meta.data %>%
  left_join(codex,by="iid")

all_combo@meta.data = all_combo@meta.data %>%
  dplyr::select(nCount_RNA, batch_id, nFeature_RNA, contains("percent."),source,processing_site, Category,Sex,OCB,fully_anonymous_pseudoid, contains("ann_celltypist"),c(202:230))

rownames(all_combo@meta.data)  = colnames(all_combo)

# save dataset
saveRDS(all_combo,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_tcell_data_for_sharing.rds")
write_csv(all_combo@meta.data,"/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/sharing_data/anonymised_tcell_data_for_sharing_metadata.csv")
