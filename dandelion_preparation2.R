
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

args = commandArgs(trailingOnly = TRUE)

print(args)
if(!args[1] %in% c("BCR","TCR")){
  stop("Argument must be one of BCR or TCR")
}

name_of_vdj = if(args[1]=="BCR"){
  "vdj_b"
} else {
  "vdj_t"
}

# set WD
setwd("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/")

# get sets of files
tum_dir="/rds/user/hpcjaco1/hpc-work/TUM data/data_featherstone/christiane/SC/transfer_CAM/TUM_part1/"
run1_samples=list.files(paste0(tum_dir,"hhz_022021","/")); run1_samples = run1_samples[substr(run1_samples,1,1) == "B"]
run2_samples=list.files(paste0(tum_dir,"hhz_052021","/")); run2_samples = run2_samples[substr(run2_samples,1,1) == "B"]
run3_samples=list.files(paste0(tum_dir,"hhz_082021","/")); run3_samples = run3_samples[substr(run3_samples,1,1) == "B"]
run4_samples=list.files(paste0(tum_dir,"hhz_012022","/")); run4_samples = run4_samples[substr(run4_samples,1,1) == "B"]
run4_samples = run4_samples[grepl("_noforce",run4_samples) == F]
run5_samples=list.files(paste0(tum_dir,"hhz_052022","/")); run5_samples = run5_samples[substr(run5_samples,1,1) == "B"]
run5_samples=run5_samples[run5_samples != "B6P14"]
run6_samples=list.files(paste0(tum_dir,"hhz_082022","/")); run6_samples = run6_samples[substr(run6_samples,1,1) == "B"]
run7_samples=list.files(paste0(tum_dir,"hhz_102022","/")); run7_samples = run7_samples[substr(run7_samples,1,1) == "B"]

outdir = paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/dandelion_inputs/",args[1],"/")

# Get Gex metadata for donor information
metadata = read_csv("./datasets/metadata_for_tum_dandelion.csv")

# split barcode field
metadata = metadata %>%
  tidyr::separate(original_barcode,sep="_",into=c("barcode","extra"))


process_vdj_data = function(run){

  wd = paste0(tum_dir,run)
  setwd(wd)
  folders = list.dirs(recursive=FALSE, full.names = FALSE)
  folders = folders[substr(folders,1,1) == "B"]
  folders = folders[grepl("_noforce",folders) == F]
  folders = folders[folders != "B6P14"]

  files = list()

  # Read in donor info for cells
  # loop through each batch to get individual donor ids

  for (i in 1:length(folders)){

    message("Starting to process batch ",folders[i])

      # read in fasta and filtered annotation files
      contig_file = paste0(folders[i],"/outs/",name_of_vdj,"/filtered_contig_annotations.csv")

      if(file.exists(contig_file)){
        contigs = read_csv(contig_file)
        fasta = read_table(paste0(folders[i],"/outs/",name_of_vdj,"/filtered_contig.fasta"),col_names=FALSE)
        # check that there's >2 B cell contigs in the data
        if(nrow(contigs)<2){
          message("Skipping - insufficient B/T cell contigs")
        } else {

          # filter metadata
          sub_metadata =  metadata %>%
            filter(batch_id == paste0("TUM_",folders[i])) %>%
            dplyr::select(barcode,iid)

          # join contigs and donor decon
          contigs_with_ids = contigs %>%
            left_join(sub_metadata, by="barcode") %>%
            dplyr::rename("donor_id" = iid)

            # filter non-singlets
            contigs_with_ids = contigs_with_ids %>% filter(!donor_id %in% c("unassigned","doublet",NA)) %>% filter(!grepl("donor",donor_id))
            # make list of unique donors for this batch
            unique_donors = unique(contigs_with_ids$donor_id)

          # split contig file into each donor
          for(j in 1:length(unique_donors)){

            # split contig file into donor
            this_donor = contigs_with_ids %>%
              filter(donor_id == unique_donors[j])
            message("Processing ",unique_donors[j])

            # split fasta into a df
            contig_names = fasta[seq(1,nrow(fasta),by=2),]$X1
            contig_sequences = fasta[seq(2,nrow(fasta),by=2),]$X1
            contig_df = data.frame(contig_names,contig_sequences)

            # split barcode names to fasta format
            new_names = list()
            for(k in c(1:length(contig_df$contig_names))){
              int_name = str_split(contig_df$contig_names[k],pattern="_",n=2)[[1]][1]
              full_name = str_remove(int_name,pattern=">")
              new_names[[k]] = full_name
            }

            contig_df$new_names = unlist(new_names)

            # filter to donor
            contig_df = contig_df %>% filter(new_names %in% this_donor$barcode)

            # get back into fasta format
            contig_df = contig_df %>% dplyr::select(-new_names)

            filtered_contig_list = list()
            for(l in 1:nrow(contig_df)){
              filtered_contig_list[[(l*2-1)]] = contig_df[l,1]
              filtered_contig_list[[(l*2)]] = contig_df[l,2]
            }

            new_fasta = data.frame(unlist(filtered_contig_list))

            # go back to old cols
            this_donor = this_donor %>% dplyr::select(colnames(contigs))

            # check that this donor has enough contigs
            if(nrow(this_donor)<2){
              message("Skipping - insufficient B/T cell contigs")
            } else {
              # write to file
              path_to_write = paste0(outdir,"TUM_",unique_donors[j])
              if(file.exists(path_to_write)){
                write_csv(this_donor,paste0(path_to_write,"/filtered_contig_annotations.csv"))
                write.table(new_fasta,paste0(path_to_write,"/filtered_contig.fasta"),col.names=F,quote=F,row.names=F)
              } else {
                system(paste0("mkdir ",path_to_write))
                write_csv(this_donor,paste0(path_to_write,"/filtered_contig_annotations.csv"))
                write.table(new_fasta,paste0(path_to_write,"/filtered_contig.fasta"),col.names=F,quote=F,row.names=F)

              }
            }
          }
        }
      }
  }
}

runs = c("hhz_022021","hhz_052021","hhz_082021","hhz_012022","hhz_052022","hhz_082022","hhz_102022")
mapply(process_vdj_data,runs)
