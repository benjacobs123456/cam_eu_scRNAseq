library(tidyverse)
args = commandArgs(trailingOnly=T)
meta = read_csv("meta_file.csv")
this_donor = unique(meta$individual)[as.numeric(args[1])]
message(this_donor)
meta = meta %>% filter(individual==this_donor)
outfile = paste0("meta_file_donor_",args[1])
write_csv(meta,outfile)
