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
total_length_v_ <- vector()
side_ <- "upper"
feature_v_ <- vector()
TEclass_v_ <- vector()
step_v_ <- vector()

# calculate centromere space
{
  feature_ <- "centromere_nt_"
  feature_v_ <- c(feature_v_,feature_)
  
  total_length_ <- arrays_data_$width %>% sum()
  names(total_length_) <- feature_
  total_length_v_ <- c(total_length_v_,total_length_)
  
  TEclass_v_ <- c(TEclass_v_,NA)
  
  step_ <- 1
  step_v_ <- c(step_v_,step_)
}

# calculate gap space
{
  feature_ <- "gap_nt_"
  feature_v_ <- c(feature_v_,feature_)
  
  total_length_ <- gap_coordinates_$gap_length_ %>% sum()
  names(total_length_) <- feature_
  total_length_v_ <- c(total_length_v_,total_length_)
  
  TEclass_v_ <- c(TEclass_v_,NA)
  
  step_ <- 2
  step_v_ <- c(step_v_,step_)
}

# calculate cent family repeats space
{
  feature_ <- "TR_nt_"
  feature_v_ <- c(feature_v_,feature_)
  
  total_length_ <- as.numeric(total_length_v_[which(names(total_length_v_) == "centromere_nt_")]) - as.numeric(total_length_v_[which(names(total_length_v_) == "gap_nt_")])
  names(total_length_) <- paste0(feature_)
  total_length_v_ <- c(total_length_v_,total_length_)
  
  TEclass_v_ <- c(TEclass_v_,NA)
  
  step_ <- 2
  step_v_ <- c(step_v_,step_)
}

# calculate non cent family repeats space
{
  if(!is.na(non_cent_family_repeats_data_test_)){
    feature_ <- "otherTR_nt_"
    feature_v_ <- c(feature_v_,feature_)
    
    total_length_ <- non_cent_family_repeats_data_$total_length_arrayID_ %>% sum()
    names(total_length_) <- paste0(feature_)
    total_length_v_ <- c(total_length_v_,total_length_)
    
    TEclass_ <- NA
    TEclass_v_ <- c(TEclass_v_,TEclass_)
    
    step_ <- 3
    step_v_ <- c(step_v_,step_)
  }else{
    feature_ <- "otherTR_nt_"
    feature_v_ <- c(feature_v_,feature_)
    
    total_length_ <- 0
    names(total_length_) <- paste0(feature_)
    total_length_v_ <- c(total_length_v_,total_length_)
    
    TEclass_ <- NA
    TEclass_v_ <- c(TEclass_v_,TEclass_)
    
    step_ <- 3
    step_v_ <- c(step_v_,step_)
    
  }
}

# calculate overall TE space
{
  feature_ <- "TE_nt_"
  feature_v_ <- c(feature_v_,feature_)
  
  total_length_ <- TEanno_isec_gaps_parsed_[which(TEanno_isec_gaps_parsed_$TEanno_start_ != -1 & TEanno_isec_gaps_parsed_$gap_length_ >= min_gap_length_ & TEanno_isec_gaps_parsed_$gap_length_ <= max_gap_length_),c(5,6,7,8,9,14)] %>%
    unique() %>%
    summarise(total_TEanno_length_ = sum(TEanno_length_))
  
  if(nrow(total_length_) > 0){
     names(total_length_) <- paste0(feature_,"0")
    total_length_v_ <- c(total_length_v_,total_length_)
    
    TEclass_ <- "all"
    TEclass_v_ <- c(TEclass_v_,TEclass_)
    
    step_ <- 3
    step_v_ <- c(step_v_,step_)
  }else{
    total_length_ <- 0
    names(total_length_) <- paste0(feature_,"0")
    total_length_v_ <- c(total_length_v_,total_length_)
    
    TEclass_ <- "all"
    TEclass_v_ <- c(TEclass_v_,TEclass_)
    
    step_ <- 3
    step_v_ <- c(step_v_,step_)
  }
}

# calculate unknown space
{
  feature_ <- "unk_nt_"
  feature_v_ <- c(feature_v_,feature_)
  
  unk_space_ <- as.numeric(total_length_v_[which(names(total_length_v_) == "gap_nt_")]) - (as.numeric(total_length_v_[which(names(total_length_v_) == "TE_nt_0")]) + as.numeric(total_length_v_[which(names(total_length_v_) == "otherTR_nt_")]))
  unk_space_ <- ifelse(unk_space_ < 0, 0, unk_space_)
  
  total_length_ <- unk_space_
  names(total_length_) <- paste0(feature_)
  total_length_v_ <- c(total_length_v_,total_length_)
  
  TEclass_ <- NA
  TEclass_v_ <- c(TEclass_v_,TEclass_)
  
  step_ <- 3
  step_v_ <- c(step_v_,step_)
}

# calculate TE space per class
{
  feature_ <- "TE_nt_"
  df_ <- TEanno_isec_gaps_parsed_[which(TEanno_isec_gaps_parsed_$TEanno_start_ != -1 & TEanno_isec_gaps_parsed_$gap_length_ >= min_gap_length_ & TEanno_isec_gaps_parsed_$gap_length_ <= max_gap_length_),c(5,6,7,8,9,14,20)] %>% unique() %>%
    group_by(TEanno_cls_) %>%
    summarise(total_lenght_ = sum(TEanno_length_))
  
  if(nrow(df_) > 0){
    total_length_ <- df_$total_lenght_
    total_classes_ <- length(df_$TEanno_cls_)
    names(total_length_) <- paste0(feature_,seq(1,total_classes_))
    total_length_v_ <- c(total_length_v_,total_length_)
    
    feature_v_ <- c(feature_v_,replicate(total_classes_,feature_))
    TEclass_v_ <- c(TEclass_v_,df_$TEanno_cls_)
    
    step_ <- 4
    step_v_ <- c(step_v_,replicate(total_classes_,step_))
  }
  
}

# calculate TE space per superfamily
{
  feature_ <- "TE_nt_"
  df_ <- TEanno_isec_gaps_parsed_[which(TEanno_isec_gaps_parsed_$TEanno_start_ != -1 & TEanno_isec_gaps_parsed_$gap_length_ >= 250 & TEanno_isec_gaps_parsed_$gap_length_ <= 100000),c(5,6,7,8,9,14)] %>%
    unique() %>%
    group_by(TEanno_feature_) %>%
    summarise(total_lenght_ = sum(TEanno_length_))
  
  if(nrow(df_) > 0){
    total_length_ <- df_$total_lenght_
    names(total_length_) <- paste0(feature_,seq(total_classes_+1,length(df_$TEanno_feature_)+total_classes_))
    total_length_v_ <- c(total_length_v_,total_length_)
    
    feature_v_ <- c(feature_v_,replicate(length(df_$TEanno_feature_),feature_))
    TEclass_v_ <- c(TEclass_v_,df_$TEanno_feature_)
    
    step_ <- 5
    step_v_ <- c(step_v_,replicate(length(df_$TEanno_feature_),step_))
  }
}

# aggregate the vectors
centromere_dissection_df_ <- data.frame("total_length_" = unlist(total_length_v_),
                                        "filter_comb_" = replicate(length(total_length_v_),filter_comb_i_),
                                        "min_gap_length_" = replicate(length(total_length_v_),min_gap_length_),
                                        "max_gap_length_" = replicate(length(total_length_v_),max_gap_length_),
                                        "side_" = replicate(length(total_length_v_),side_),
                                        "feature_" = feature_v_,
                                        "TEclass_" = TEclass_v_,
                                        "step_" = step_v_)

# calculate centromere and genome wide prop
genome_size_ <- karyo_$seq_length_ %>% sum()
centromere_size_ <- centromere_dissection_df_$total_length_[which(centromere_dissection_df_$feature_ == "centromere_nt_")]

centromere_dissection_df_ <- centromere_dissection_df_ %>%
  mutate(centromere_prop_ = total_length_/centromere_size_,
         genome_prop_ = total_length_/genome_size_)

# save the results
write.table(centromere_dissection_df_,paste0(root_name_,"_centromere_dissection_data_gw_"),quote = F, sep = "\t", col.names = T, row.names = F)
