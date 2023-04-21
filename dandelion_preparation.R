# load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# get command line args (BCR vs TCR)
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

# read in link file
reference = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/EU_ID.txt")
reference$cohort = "EU"
reference$PBMC_PoolSize = as.character(reference$PBMC_PoolSize)
reference_cam = read_tsv("/rds/project/sjs1016/rds-sjs1016-msgen/10X_5prime/GT_ID.txt")
reference_cam$cohort = "Cam"
reference_cam$PBMC_PoolSize = as.character(reference_cam$PBMC_PoolSize)
reference = bind_rows(reference,reference_cam)

# define list of good batches
good_pbmc_batches = unique(reference$PBMC_GEX)
good_csf_batches = unique(reference$CSF_GEX)

# get sets of files
eu_pbmc_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_EU")
eu_pbmc_samples = eu_pbmc_samples[eu_pbmc_samples %in% good_pbmc_batches]
cam_pbmc_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_output/PBMC/")
cam_pbmc_samples = cam_pbmc_samples[cam_pbmc_samples %in% good_pbmc_batches]
eu_csf_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_EU")
eu_csf_samples = eu_csf_samples[eu_csf_samples %in% good_csf_batches]
cam_csf_samples = list.files("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/CRv5_output/CSF/")
cam_csf_samples = cam_csf_samples[cam_csf_samples %in% good_csf_batches]
unpooled_csf_samples = reference[reference$CSF_PoolSize==1,]
unpooled_pbmc_samples = reference[reference$PBMC_PoolSize==1,]

batches_without_decon = list()

# set out directory
outdir = paste0("/rds/user/hpcjaco1/hpc-work/Cambridge_EU_combined/dandelion_inputs/",args[1],"/")

# define processing function
process_vdj_data = function(source, cohort){

# error checking
if(!source %in% c("CSF","PBMC")){
  stop("Source must be one of CSF or PBMC")
}
if(!cohort %in% c("EU","Cam")){
  stop("Source must be one of CSF or PBMC")
}

# set wd
wd = if(cohort == "Cam"){
    message("Looking in Cambridge folders")
    paste0("/rds-d4/project/sjs1016/rds-sjs1016-msgen/10X_5prime/CRv5_output/",source)
  } else if(cohort == "EU") {
    message("Looking in EU folders")
    paste0("/rds-d4/project/sjs1016/rds-sjs1016-msgen/10X_5prime/CRv5_EU/")
  } # set the working directory
setwd(wd)

# get folders
folders = list.dirs(recursive=FALSE, full.names = FALSE)

# filter folders to correct folders given source and cohort
folders = if(cohort == "Cam" & source == "CSF"){
  folders[folders %in% cam_csf_samples]
} else if(cohort == "Cam" & source == "PBMC"){
  folders[folders %in% cam_pbmc_samples]
} else if(cohort == "EU" & source == "CSF"){
  folders[folders %in% eu_csf_samples]
} else if(cohort == "EU" & source == "PBMC"){
  folders[folders %in% eu_pbmc_samples]
}

files = list()

# loop through each batch to get individual donor ids

for (i in 1:length(folders)){
    message("Starting to process batch ",folders[i])
    if(folders[i] %in% unpooled_csf_samples$CSF_GEX | folders[i] %in% unpooled_pbmc_samples$PBMC_GEX){ # if unpooled
      message("Batch ", folders[i]," is unpooled. Processing")

      # read in fasta and filtered annotation files
      contig_file = paste0(folders[i],"/outs/",name_of_vdj,"/filtered_contig_annotations.csv")
      fasta_file = paste0(folders[i],"/outs/",name_of_vdj,"/filtered_contig.fasta")


      if( file.exists(contig_file) & file.exists(fasta_file)){
        contigs = read_csv(contig_file)
        fasta = read_table(fasta_file,col_names=FALSE)

        # check that there's >=2 B/T cell contigs in the data
        if(nrow(contigs)<2){
          message("Skipping - insufficient B/T cell contigs")
        } else {

        donor = if(cohort == "EU"){ # get donor ID
          if(folders[i] %in% unpooled_csf_samples$CSF_GEX){
            unpooled_csf_samples[unpooled_csf_samples$CSF_GEX==folders[i],]$ID
          } else if(folders[i] %in% unpooled_pbmc_samples$PBMC_GEX){
            unpooled_pbmc_samples[unpooled_pbmc_samples$PBMC_GEX==folders[i],]$ID
          }
        } else {
          if(folders[i] %in% unpooled_csf_samples$CSF_GEX){
            unpooled_csf_samples[unpooled_csf_samples$CSF_GEX==folders[i],]$ShortID
          } else if(folders[i] %in% unpooled_pbmc_samples$PBMC_GEX){
            unpooled_pbmc_samples[unpooled_pbmc_samples$PBMC_GEX==folders[i],]$ShortID
          }
        }
        message("Donor: ", donor)

        path_to_write = paste0(outdir,source,"_",donor)
        if(!file.exists(path_to_write)){
          system(paste0("mkdir ",path_to_write))
        }
        system(paste0("cp ",fasta_file," ", path_to_write))
        system(paste0("cp ",contig_file," ", path_to_write))
        message("written files")
        }
      } else {
        message("Those files didn't exist")
        message(source)
        message(folders[i])
      }
    } else { # i.e. if multiplexed

      # read in fasta and filtered annotation files
      contig_file = paste0(folders[i],"/outs/",name_of_vdj,"/filtered_contig_annotations.csv")
      decon_file = paste0("/home/hpcjaco1/rds/rds-sjs1016-msgen/10X_5prime/Donor_Decon/",folders[i],"_donor_ids.tsv")

      if(file.exists(decon_file) & file.exists(contig_file)){
        contigs = read_csv(contig_file)
        decon = read_tsv(decon_file)
        fasta = read_table(paste0(folders[i],"/outs/",name_of_vdj,"/filtered_contig.fasta"),col_names=FALSE)
        # check that there's >2 B cell contigs in the data
        if(nrow(contigs)<2){
          message("Skipping - insufficient B/T cell contigs")
        } else {
          # join contigs and donor decon
          contigs_with_ids = contigs %>% left_join(decon %>% dplyr::rename("barcode"="cell"), by="barcode")
          # filter non-singlets
          contigs_with_ids = contigs_with_ids %>% filter(!donor_id %in% c("unassigned","doublet",NA)) %>% filter(!grepl("donor",donor_id))
          # make list of unique donors for this batch
          unique_donors = unique(contigs_with_ids$donor_id)

          # split contig file into each donor
            for(j in 1:length(unique_donors)){

              # split contig file into donor
              this_donor = contigs_with_ids %>% filter(donor_id == unique_donors[j])
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
                path_to_write = paste0(outdir,source,"_",unique_donors[j])
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
      } else {
          if(!file.exists(decon_file)){
            message("Decon file does not exist for ",folders[i])
            batches_without_decon[[length(batches_without_decon)+1]] = folders[i]
          }
          if(!file.exists(contig_file)){
            message("Contig file does not exist for ",folders[i])
          }
        }
      }
    }
}



sources = c("PBMC","CSF")
cohorts = c("EU","Cam")
mat = expand.grid(sources,cohorts)
mapply(process_vdj_data,mat[,1],mat[,2])
