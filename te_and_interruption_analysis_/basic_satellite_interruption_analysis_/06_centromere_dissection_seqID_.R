#!/bin/Rscript

suppressMessages(library(dplyr))
options(scipen = 999)
args <- commandArgs(TRUE)

# read custom lib
{
  bin_path_ <- strsplit(args[1],"=")[[1]][2]
  source(paste0(bin_path_,"/custom_func_lib_.R"))
}

# read args
{
  opts_ <- parse_args(args)
  
  arrays_file_ <- opts_$arrays_file_ %||% stop("arrays_file_ required")
  TEanno_isec_gaps_parsed_file_ <- opts_$TEanno_isec_gaps_parsed_file_ %||% stop("TEanno_isec_gaps_parsed_file_ required")
  gap_coordinates_file_ <- opts_$gap_coordinates_file_ %||% stop("gap_coordinates_file_ required")
  non_cent_family_repeats_data_file_ <- opts_$non_cent_family_repeats_data_file_
  karyo_file_ <- opts_$karyo_file_ %||% stop("karyo_file_ required")
  min_gap_length_ <- as.numeric(opts_$min_gap_length_) %||% 250
  max_gap_length_ <- as.numeric(opts_$max_gap_length_) %||% 100000
  filter_comb_i_  <- opts_$filter_comb_ %||% NULL
}

# get input data
{
  # get arrays data
  file_i_ <- arrays_file_
  if (file.info(file_i_)$size > 0) {
    arrays_data_ <- read.table(arrays_file_, header = T,sep = ",")
    rm(file_i_)
  }else{
    cat("Missing arrays_data_ file\n")
    quit(status = 1)
  }
  
  # get TEanno_isec_gaps_parsed_file_ data
  file_i_ <- TEanno_isec_gaps_parsed_file_
  if (file.info(file_i_)$size > 0) {
    TEanno_isec_gaps_parsed_ <- read.table(TEanno_isec_gaps_parsed_file_, header = T,sep = "\t")
    rm(file_i_)
  }else{
    cat("Missing TEanno_isec_gaps_parsed_ file\n")
    quit(status = 1)
  }
  
  # get gap_coordinates_file_ data
  file_i_ <- gap_coordinates_file_
  if (file.info(file_i_)$size > 0) {
    gap_coordinates_ <- read.table(gap_coordinates_file_, header = T,sep = "\t")
    rm(file_i_)
  }else{
    cat("Missing gap_coordinates_ file\n")
    quit(status = 1)
  }
  
  # get non_cent_family_repeats_data_file_ data
  file_i_ <- non_cent_family_repeats_data_file_
  if (file.info(file_i_)$size > 0) {
    non_cent_family_repeats_data_ <- read.table(non_cent_family_repeats_data_file_, header = T,sep = "\t")
    non_cent_family_repeats_data_test_ <- "PASS"
    rm(file_i_)
  }else{
    non_cent_family_repeats_data_test_ <- NA
  }
  
  # get non_cent_family_repeats_data_file_ data
  file_i_ <- karyo_file_
  if (file.info(file_i_)$size > 0) {
    karyo_ <- read.table(karyo_file_, header = F, sep = "\t")
    colnames(karyo_) <- c("seqID_","ChrOnly_seqID_","seq_length_")
    rm(file_i_)
  }else{
    cat("Missing karyo_ file\n")
    quit(status = 1)
  }
  
  TE_class_list_updated_ <- read.table(paste0(bin_path_,"/TE_class_list_updated_"),sep = "\t",header = F)
  colnames(TE_class_list_updated_) <- c("TEanno_cls_","old_TEanno_cls_","TEanno_feature_")
  
}

root_name_ <- TEanno_isec_gaps_parsed_file_ %>%
  {gsub(".*/","",.)} %>%
  {gsub("_TEanno_isec_gaps_parsed_","",.)}
print(root_name_)

arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))

# keep relvant arrays
if(!is.null(filter_comb_i_)){
  arrays_data_ <- arrays_data_[which(arrays_data_$filter_comb_ %in% filter_comb_i_),]
}

TEanno_isec_gaps_parsed_ <- left_join(TEanno_isec_gaps_parsed_,TE_class_list_updated_[,c(1,3)], by = "TEanno_feature_")

# prepare variables
total_length_ <- data.frame()

# calculate centromere space
{
  arrays_space_ <- arrays_data_ %>%
    group_by(seqID) %>%
    summarise(total_length_ = sum(width))
  total_length_df_ <- arrays_space_
  colnames(arrays_space_)[2] <- "total_length_arrays_"
  
  total_length_df_$filter_comb_ <- filter_comb_i_
  total_length_df_$min_gap_length_ <- min_gap_length_
  total_length_df_$max_gap_length_ <- max_gap_length_
  total_length_df_$side_ <- "upper"
  total_length_df_$feature_ <- "centromere_nt_"
  total_length_df_$TEclass_ <- NA
  total_length_df_$step_ <- 1
  
  total_length_ <- rbind(total_length_,total_length_df_)
}

# calculate gap space
{
  gap_space_ <- gap_coordinates_ %>%
    group_by(seqID) %>%
    summarise(total_length_ = sum(gap_length_))
  total_length_df_ <- gap_space_
  colnames(gap_space_)[2] <- "total_length_gaps_"
  
  total_length_df_$filter_comb_ <- filter_comb_i_
  total_length_df_$min_gap_length_ <- min_gap_length_
  total_length_df_$max_gap_length_ <- max_gap_length_
  total_length_df_$side_ <- "upper"
  total_length_df_$feature_ <- "gap_nt_"
  total_length_df_$TEclass_ <- NA
  total_length_df_$step_ <- 2
  
  total_length_ <- rbind(total_length_,total_length_df_)
}

# calculate cent family repeats space
{
  total_length_df_ <- left_join(arrays_space_,gap_space_,by = "seqID") %>%
    group_by(seqID) %>%
    mutate(total_length_ = total_length_arrays_-total_length_gaps_) %>%
    mutate(total_length_ = ifelse(is.na(total_length_) == F, total_length_,total_length_arrays_)) %>%
    select(!c("total_length_gaps_","total_length_arrays_"))
  
  total_length_df_$filter_comb_ <- filter_comb_i_
  total_length_df_$min_gap_length_ <- min_gap_length_
  total_length_df_$max_gap_length_ <- max_gap_length_
  total_length_df_$side_ <- "upper"
  total_length_df_$feature_ <- "TR_nt_"
  total_length_df_$TEclass_ <- NA
  total_length_df_$step_ <- 2
  
  total_length_ <- rbind(total_length_,total_length_df_)
}

# calculate non cent family repeats space
{
  if(!is.na(non_cent_family_repeats_data_test_)){
    
    otherTR_space_ <- non_cent_family_repeats_data_ %>%
      group_by(seqID) %>%
      summarise(total_length_ = sum(total_length_arrayID_))
    total_length_df_ <- otherTR_space_
    colnames(otherTR_space_)[2] <- "total_length_otherTR_"
    
    total_length_df_$filter_comb_ <- filter_comb_i_
    total_length_df_$min_gap_length_ <- min_gap_length_
    total_length_df_$max_gap_length_ <- max_gap_length_
    total_length_df_$side_ <- "upper"
    total_length_df_$feature_ <- "otherTR_nt_"
    total_length_df_$TEclass_ <- NA
    total_length_df_$step_ <- 3
    
    total_length_ <- rbind(total_length_,total_length_df_)
  }else{
    otherTR_space_ <- data.frame("seqID" = unique(total_length_$seqID),
                                   "total_length_" = 0)
    total_length_df_ <- otherTR_space_
    colnames(otherTR_space_)[2] <- "total_length_otherTR_"
    
    total_length_df_$filter_comb_ <- filter_comb_i_
    total_length_df_$min_gap_length_ <- min_gap_length_
    total_length_df_$max_gap_length_ <- max_gap_length_
    total_length_df_$side_ <- "upper"
    total_length_df_$feature_ <- "otherTR_nt_"
    total_length_df_$TEclass_ <- NA
    total_length_df_$step_ <- 3
    
    total_length_ <- rbind(total_length_,total_length_df_)
  }
}

# calculate overall TE space
{
  TEanno_space_ <- TEanno_isec_gaps_parsed_[which(TEanno_isec_gaps_parsed_$TEanno_start_ != -1 & TEanno_isec_gaps_parsed_$gap_length_ >= min_gap_length_ & TEanno_isec_gaps_parsed_$gap_length_ <= max_gap_length_),c(5,6,7,8,9,14)] %>%
    unique() %>%
    group_by(TEanno_seqID_) %>%
    summarise(total_TEanno_length_ = sum(TEanno_length_))
  
  if(nrow(TEanno_space_) > 0){
    colnames(TEanno_space_) <- colnames(TEanno_space_) %>%
      {gsub("TEanno_","",.)} %>%
      {gsub("seqID_","seqID",.)}
    
    total_length_df_ <- TEanno_space_
    colnames(TEanno_space_)[2] <- "total_length_TEanno_"
    
    total_length_df_$filter_comb_ <- filter_comb_i_
    total_length_df_$min_gap_length_ <- min_gap_length_
    total_length_df_$max_gap_length_ <- max_gap_length_
    total_length_df_$side_ <- "upper"
    total_length_df_$feature_ <- "TE_nt_"
    total_length_df_$TEclass_ <- "all"
    total_length_df_$step_ <- 3
    
    total_length_ <- rbind(total_length_,total_length_df_)
  }else{
    
    total_length_df_ <- data.frame("seqID" = unique(total_length_$seqID),
                                   "total_length_" = 0)
    
    TEanno_space_ <- total_length_df_
    colnames(TEanno_space_)[2] <- "total_length_TEanno_"
    
    total_length_df_$filter_comb_ <- filter_comb_i_
    total_length_df_$min_gap_length_ <- min_gap_length_
    total_length_df_$max_gap_length_ <- max_gap_length_
    total_length_df_$side_ <- "upper"
    total_length_df_$feature_ <- "TE_nt_"
    total_length_df_$TEclass_ <- "all"
    total_length_df_$step_ <- 3
    
    total_length_ <- rbind(total_length_,total_length_df_)
  }
}

# calculate unknown space
{
  
  TEanno_and_otherTR_ <- left_join(TEanno_space_,otherTR_space_,by = "seqID") %>%
    group_by(seqID) %>%
    mutate(total_length_TEanno_and_otherTR_ = total_length_TEanno_+total_length_otherTR_) %>%
    mutate(total_length_TEanno_and_otherTR_ = case_when(is.na(total_length_TEanno_and_otherTR_) == F ~ total_length_TEanno_and_otherTR_,
                                     is.na(total_length_otherTR_) == T ~ total_length_TEanno_,
                                     is.na(total_length_TEanno_) == T ~ total_length_otherTR_,
                                     is.na(total_length_TEanno_) == T & is.na(total_length_otherTR_) == T ~ 0,
                                     .default = NA)) %>%
    select(!c("total_length_TEanno_","total_length_otherTR_"))
  
  # quantifying degree of overlap
  {
    sink("TEanno_overlap_within_gaps_")
    left_join(gap_space_,TEanno_and_otherTR_,by = "seqID") %>%
      group_by(seqID) %>%
      mutate(total_length_ = case_when(is.na(total_length_gaps_) == F & is.na(total_length_TEanno_and_otherTR_) == F ~ total_length_gaps_-total_length_TEanno_and_otherTR_,
                                       is.na(total_length_gaps_) == F & is.na(total_length_TEanno_and_otherTR_) == T ~ total_length_gaps_,
                                       .default = NA)) %>%
      filter(total_length_ < 0) %>%
      select(seqID,total_length_) %>% mutate(root_name_ = root_name_)
    sink()
  }
  
  total_length_df_ <- left_join(gap_space_,TEanno_and_otherTR_,by = "seqID") %>%
    group_by(seqID) %>%
    mutate(total_length_ = case_when(is.na(total_length_gaps_) == F & is.na(total_length_TEanno_and_otherTR_) == F ~ total_length_gaps_-total_length_TEanno_and_otherTR_,
                                     is.na(total_length_gaps_) == F & is.na(total_length_TEanno_and_otherTR_) == T ~ total_length_gaps_,
                                     .default = NA)) %>%
    select(!c("total_length_gaps_","total_length_TEanno_and_otherTR_")) %>%
    mutate(total_length_ = ifelse(total_length_ < 0, 0 , total_length_))


  total_length_df_$filter_comb_ <- filter_comb_i_
  total_length_df_$min_gap_length_ <- min_gap_length_
  total_length_df_$max_gap_length_ <- max_gap_length_
  total_length_df_$side_ <- "upper"
  total_length_df_$feature_ <- "unk_nt_"
  total_length_df_$TEclass_ <- NA
  total_length_df_$step_ <- 3
  
  total_length_ <- rbind(total_length_,total_length_df_)
  
}

# calculate TE space per class
{
  total_length_df_ <- TEanno_isec_gaps_parsed_[which(TEanno_isec_gaps_parsed_$TEanno_start_ != -1 & TEanno_isec_gaps_parsed_$gap_length_ >= min_gap_length_ & TEanno_isec_gaps_parsed_$gap_length_ <= max_gap_length_),c(5,6,7,8,9,14,20)] %>% unique() %>%
    group_by(TEanno_cls_,TEanno_seqID_) %>%
    summarise(total_length_ = sum(TEanno_length_))
  
  if(nrow(total_length_df_) > 0){
    colnames(total_length_df_) <- colnames(total_length_df_) %>%
      {gsub("TEanno_seqID_","seqID",.)}
    colnames(total_length_df_) <- colnames(total_length_df_) %>%
      {gsub("TEanno_cls_","TEclass_",.)}
    
    total_length_df_$filter_comb_ <- filter_comb_i_
    total_length_df_$min_gap_length_ <- min_gap_length_
    total_length_df_$max_gap_length_ <- max_gap_length_
    total_length_df_$side_ <- "upper"
    total_length_df_$feature_ <- "TE_nt_"
    total_length_df_$step_ <- 4
    
    total_length_ <- rbind(total_length_,total_length_df_[,c(2,3,4,5,6,7,8,1,9)])
  }
  
}

# calculate TE space per superfamily
{
  total_length_df_ <- TEanno_isec_gaps_parsed_[which(TEanno_isec_gaps_parsed_$TEanno_start_ != -1 & TEanno_isec_gaps_parsed_$gap_length_ >= 250 & TEanno_isec_gaps_parsed_$gap_length_ <= 100000),c(5,6,7,8,9,14)] %>%
    unique() %>%
    group_by(TEanno_feature_,TEanno_seqID_) %>%
    summarise(total_length_ = sum(TEanno_length_))
  
  if(nrow(total_length_df_) > 0){
    colnames(total_length_df_) <- colnames(total_length_df_) %>%
      {gsub("TEanno_seqID_","seqID",.)}
    colnames(total_length_df_) <- colnames(total_length_df_) %>%
      {gsub("TEanno_feature_","TEclass_",.)}
    
    total_length_df_$filter_comb_ <- filter_comb_i_
    total_length_df_$min_gap_length_ <- min_gap_length_
    total_length_df_$max_gap_length_ <- max_gap_length_
    total_length_df_$side_ <- "upper"
    total_length_df_$feature_ <- "TE_nt_"
    total_length_df_$step_ <- 5
    
    total_length_ <- rbind(total_length_,total_length_df_[,c(2,3,4,5,6,7,8,1,9)])
    colnames(total_length_) <- colnames(total_length_) %>%
      {gsub("seqID","seqID_",.)}
  }

}

# calculate centromere and genome wide prop
colnames(total_length_)[1] <- "seqID_"
centromere_dissection_df_ <- left_join(total_length_,karyo_[,c(1,3)],by = "seqID_")

genome_size_ <- karyo_$seq_length_ %>% sum()
centromere_size_ <- centromere_dissection_df_$total_length_[which(centromere_dissection_df_$feature_ == "centromere_nt_")]

centromere_size_df_ <- centromere_dissection_df_[which(centromere_dissection_df_$feature_ == "centromere_nt_"),] %>%
  group_by(seqID_) %>%
  summarise(centromere_size_ = sum(total_length_))
centromere_dissection_df_ <- left_join(centromere_dissection_df_,centromere_size_df_,by = "seqID_")

centromere_dissection_df_ <- centromere_dissection_df_ %>%
  mutate(centromere_prop_ = total_length_/centromere_size_,
         genome_prop_ = total_length_/seq_length_)

# save the results
write.table(centromere_dissection_df_,paste0(root_name_,"_centromere_dissection_data_seqID_"),quote = F, sep = "\t", col.names = T, row.names = F)
