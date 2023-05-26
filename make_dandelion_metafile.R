library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/dandelion_inputs")


### make meta file

make_metafile = function(vdj){

  # change wd
  wd = paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/dandelion_inputs/",vdj)
  setwd(wd)

  # list all samples
  samples = list.dirs("./",recursive=FALSE, full.names = FALSE)

  # make into df
  meta_file = data.frame(samples)

  # grab sample id
  meta_file = meta_file %>% mutate(new_sample_id = samples)
  meta_file = meta_file %>% cbind(data.frame(str_split_fixed(meta_file$new_sample_id,pattern="_",n=2))) %>%
  dplyr::select(-new_sample_id, - X1) %>%
  dplyr::rename("sample" = samples, "individual" = X2) %>%
  dplyr::select(sample,individual)

  write_csv(meta_file,"meta_file.csv")
}

make_metafile("BCR")
make_metafile("TCR")
