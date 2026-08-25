#!/bin/Rscript

suppressMessages(library(dplyr))
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)


# read custom lib
{
  bin_path_ <- strsplit(args[1],"=")[[1]][2]
  source(paste0(bin_path_,"/custom_func_lib_.R"))
}

# read args
{
  opts_ <- parse_args(args)
  
  repeats_file_ <- opts_$repeats_file_ %||% stop("repeats_file_ required")
  arrays_file_ <- opts_$arrays_file_ %||% stop("arrays_file_ required")
  min_gap_length_ <- as.numeric(opts_$min_gap_length_) %||% 250
  max_gap_length_ <- as.numeric(opts_$max_gap_length_) %||% 100000
  filter_comb_i_  <- opts_$filter_comb_ %||% NULL
}


# get input data
{
  # get repeats data
  file_i_ <- repeats_file_
  if (file.info(file_i_)$size > 0) {
    repeats_data_ <- data.table::fread(repeats_file_, header = T,sep = ",") %>% data.frame()
    # for data generated with TRASH1
    if(length(colnames(repeats_data_)[grepl("new_class$",colnames(repeats_data_))]) == 0){
      repeats_data_$new_class <- repeats_data_$class
    }
    rm(file_i_)
  }else{
    cat("Missing repeats_data_ file\n")
    quit(status = 1)
  }
  
  # get arrays data
  file_i_ <- arrays_file_
  if (file.info(file_i_)$size > 0) {
    arrays_data_ <- read.table(arrays_file_, header = T,sep = ",")
    rm(file_i_)
  }else{
    cat("Missing arrays_data_ file\n")
    quit(status = 1)
  }
}

root_name_ <- arrays_file_ %>%
  {gsub(".*/","",.)} %>%
  {gsub("_fam_.*_arrays_expanded.csv","",.)} %>%
  {gsub(".fasta","",.)} %>%
  {gsub(".fa","",.)}

print(root_name_)

# parse input files
arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))

# keep relvant arrays
if(!is.null(filter_comb_i_)){
  arrays_data_ <- arrays_data_[which(arrays_data_$filter_comb_ %in% filter_comb_i_),]
}

# get centromeric TR classes
arrays_data_$TRASH_class <- as.character(arrays_data_$TRASH_class)
TRASH_class_ <- arrays_data_$TRASH_class %>% unique() %>% strsplit(.,";") %>% unlist()

# keep non cent family monomers
repeats_data_ <- repeats_data_[-which(repeats_data_$new_class %in% TRASH_class_),]

# intersec non cent family monomers with array boundaries
seqID_ <- arrays_data_$seqID %>% unique()
array_repeats_data_ <- data.frame()
for(i_ in 1:length(seqID_)){
  seqID_i_ <- seqID_[i_]
  print(seqID_i_)
  idx_ <- which(arrays_data_$seqID %in% seqID_i_)
  for(j_ in 1:length(idx_)){
    start_j_ <- arrays_data_$start[idx_[j_]]
    end_j_ <- arrays_data_$end[idx_[j_]]
    arrayID_j_ <- arrays_data_$arrayID[idx_[j_]]
    filter_comb_j_ <- arrays_data_$filter_comb_[idx_[j_]]
    idx_j_ <- which(repeats_data_$seqID %in% seqID_i_ & repeats_data_$start >= start_j_ & repeats_data_$end <= end_j_)
    if(length(idx_j_) > 0){
      repeats_data_j_ <- repeats_data_[idx_j_,] %>%
        select(seqID,start,end,new_class)
      repeats_data_j_$arrayID <- arrayID_j_
      repeats_data_j_$filter_comb_ <- filter_comb_j_
      array_repeats_data_ <- rbind(array_repeats_data_,repeats_data_j_)
    }    
  }
}

# summarise the data
if(nrow(array_repeats_data_) > 0){
  array_repeats_data_ <- array_repeats_data_[with(array_repeats_data_,order(seqID,start)),]
  array_repeats_data_sum_ <- array_repeats_data_ %>%
    mutate(width_ = (end-start)) %>%
    group_by(filter_comb_,seqID,arrayID) %>%
    summarise(total_length_arrayID_ = sum(width_)) %>%
    ungroup() %>%
    group_by(filter_comb_,seqID) %>%
    mutate(total_length_seqID_ = sum(total_length_arrayID_)) %>%
    ungroup() %>%
    group_by(filter_comb_) %>%
    mutate(total_length_ = sum(total_length_seqID_))
  
  # save the results
  write.table(array_repeats_data_sum_,paste0(root_name_,"_non_cent_family_repeats_data_"),quote = F, sep = "\t", col.names = T, row.names = F)
}else{
  cat("No non_cent_family_repeats found\n")
  quit()
}
