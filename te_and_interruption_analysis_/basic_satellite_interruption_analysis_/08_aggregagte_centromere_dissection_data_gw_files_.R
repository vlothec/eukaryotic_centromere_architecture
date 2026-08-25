#!/bin/Rscript

suppressMessages(library(dplyr))
options(scipen = 999)
args <- commandArgs(TRUE)

{
  bin_path_ <- "/path/to/bin_"
  path_to_data_ <- "/path/to/data_"
}

source(paste0(bin_path_,"/custom_func_lib_.R"))

# read args
# get input data
TE_class_list_updated_ <- read.table(paste0(bin_path_,"/TE_class_list_updated_"),sep = "\t",header = F)
colnames(TE_class_list_updated_) <- c("TEanno_cls_","old_TEanno_cls_","TEanno_feature_")
TE_class_list_updated_red_ <- data.frame("TEclass_" = unique(TE_class_list_updated_$TEanno_cls_))

DToL_info_ <- read.table(paste0(bin_path_,"/DToL_info_.txt"),sep = "\t",header = F)
colnames(DToL_info_) <- c("subset_","fasta_name_","phylotaxa_", "sub_phylotaxa_i_", "sub_phylotaxa_ii_", "genome_size_")
DToL_info_$dtolID_ <- DToL_info_$fasta_name_ %>%
  {gsub(".fasta","",.)} %>%
  {gsub(".fa","",.)} %>%
  {gsub("rosCan_S27_v1","rosCan",.)}

centromere_dissection_data_gw_file_list_ <- list.files(path = path_to_data_, pattern = "_centromere_dissection_data_gw_$",recursive = T)
centromere_dissection_data_gw_all_ <- data.frame()
for(i_ in 1:length(centromere_dissection_data_gw_file_list_)){
  centromere_dissection_data_gw_file_ <- centromere_dissection_data_gw_file_list_[i_]
  file_i_ <- paste0(path_to_data_,"/",centromere_dissection_data_gw_file_)
  if (file.info(file_i_)$size > 0) {
    dtolID_ <- centromere_dissection_data_gw_file_ %>%
      {gsub("_gap_landscape_/.*","",.)}
    print(dtolID_)
    centromere_dissection_data_gw_ <- read.table(paste0(path_to_data_,"/",centromere_dissection_data_gw_file_),header = T,sep = "\t")
    rm(file_i_)
    
    # get dtolID_
    centromere_dissection_data_gw_ <- full_join(centromere_dissection_data_gw_,TE_class_list_updated_red_,by ="TEclass_")
    centromere_dissection_data_gw_$dtolID_ <- dtolID_
    # get phylotaxa
    centromere_dissection_data_gw_ <- left_join(centromere_dissection_data_gw_,DToL_info_[,c(7,3)],by = "dtolID_")
    centromere_dissection_data_gw_$total_length_[which(is.na(centromere_dissection_data_gw_$total_length_ == T))] <- 0
    centromere_dissection_data_gw_$step_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- 4
    
    centromere_dissection_data_gw_$filter_comb_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- unique(centromere_dissection_data_gw_$filter_comb_[which(centromere_dissection_data_gw_$total_length_ != 0)])
    centromere_dissection_data_gw_$min_gap_length_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- unique(centromere_dissection_data_gw_$min_gap_length_[which(centromere_dissection_data_gw_$total_length_ != 0)])
    centromere_dissection_data_gw_$max_gap_length_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- unique(centromere_dissection_data_gw_$max_gap_length_[which(centromere_dissection_data_gw_$total_length_ != 0)])
    centromere_dissection_data_gw_$side_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- unique(centromere_dissection_data_gw_$side_[which(centromere_dissection_data_gw_$total_length_ != 0)])
    centromere_dissection_data_gw_$feature_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- "TE_nt_"
    centromere_dissection_data_gw_$centromere_prop_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- 0
    centromere_dissection_data_gw_$genome_prop_[which(centromere_dissection_data_gw_$total_length_ == 0)] <- 0
    
    centromere_dissection_data_gw_$phyloclade_ <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
                                                                               x == "Monocots" ~ "plant",
                                                                               x == "Aves" ~ "chordate",
                                                                               x == "Mammalia" ~ "chordate",
                                                                               x == "Actinopterygii" ~ "chordate",
                                                                               x == "Hymenoptera" ~ "invertebrate",
                                                                               x == "Annelida" ~ "invertebrate",
                                                                               x == "Coleoptera" ~ "invertebrate",
                                                                               x == "Diptera" ~ "invertebrate",
                                                                               x == "Hemiptera" ~ "invertebrate",
                                                                               x == "Cnidaria" ~ "invertebrate",
                                                                               x == "Lepidoptera" ~ "invertebrate",
                                                                               x == "Neuroptera" ~ "invertebrate",
                                                                               x == "Tunicata" ~ "chordate",
                                                                               x == "Reptilia" ~ "chordate",
                                                                               x == "Bryozoa" ~ "invertebrate",
                                                                               x == "Mollusca" ~ "invertebrate",
                                                                               x == "Plecoptera" ~ "invertebrate",
                                                                               x == "Echinodermata" ~ "invertebrate",
                                                                               x == "Blattodea" ~ "invertebrate",
                                                                               x == "Ephemeroptera" ~ "invertebrate",
                                                                               x == "Porifera" ~ "invertebrate",
                                                                               x == "Arthropoda" ~ "invertebrate",
                                                                               x == "Nemertea" ~ "invertebrate",
                                                                               .default = NA_character_),centromere_dissection_data_gw_$phylotaxa_)
    
    # calculate pure TE proportion
    centromere_dissection_data_gw_ <- centromere_dissection_data_gw_ %>%
      group_by(step_) %>%
      mutate(pure_TE_prop_ = total_length_/sum(total_length_))
    # aggregate
    centromere_dissection_data_gw_all_ <- rbind(centromere_dissection_data_gw_all_,centromere_dissection_data_gw_)
  }else{
    cat("Missing centromere_dissection_data_gw_ file\n")
    quit(status = 1)
  }
}

write.table(centromere_dissection_data_gw_all_,"centromere_dissection_data_gw_mono_aggregated_", col.names = T, row.names = F, sep = "\t", quote = F)
