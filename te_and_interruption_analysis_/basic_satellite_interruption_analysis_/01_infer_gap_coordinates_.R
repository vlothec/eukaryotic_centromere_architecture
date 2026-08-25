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
  {gsub("_arrays_expanded.csv","",.)} %>%
  {gsub("_fam_.*","",.)} %>%  
  {gsub(".fasta","",.)} %>%
  {gsub(".fa","",.)}
print(root_name_)

arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))

# keep relvant arrays
if(!is.null(filter_comb_i_)){
  arrays_data_ <- arrays_data_[which(arrays_data_$filter_comb_ %in% filter_comb_i_),]
}

# get centromeric TR classes
arrays_data_$TRASH_class <- as.character(arrays_data_$TRASH_class)
TRASH_class_ <- arrays_data_$TRASH_class %>% unique() %>% strsplit(.,";") %>% unlist()


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
    idx_j_ <- which(repeats_data_$seqID %in% seqID_i_ & repeats_data_$start >= start_j_ & repeats_data_$end <= end_j_ & repeats_data_$new_class %in% TRASH_class_)
    if(length(idx_j_) > 0){
      repeats_data_j_ <- repeats_data_[idx_j_,] %>%
        select(seqID,start,end,new_class)
      repeats_data_j_$arrayID <- arrayID_j_
      repeats_data_j_$filter_comb_ <- filter_comb_j_
      array_repeats_data_ <- rbind(array_repeats_data_,repeats_data_j_)
    }
  }
}

array_repeats_data_ <- array_repeats_data_[with(array_repeats_data_,order(seqID,start)),]
# infer gap coordinates
gap_data_ <- array_repeats_data_ %>%
  group_by(seqID,arrayID) %>%
  mutate(start_ = end+1,
         end_ = lead(start)-1)
gap_data_ <- gap_data_[,c(1,5,6,7,8)]

gap_data_ <- left_join(gap_data_,arrays_data_[,c(16,2,3)],by ="arrayID")
colnames(gap_data_)[c(6,7)] <- c("array_start_","array_end_")
gap_data_$end_[which(is.na(gap_data_$end_) == T)] <- gap_data_$array_end_[which(is.na(gap_data_$end_) == T)]
gap_data_ <- gap_data_[,-c(6,7)]

gap_data_ <- gap_data_ %>%
  mutate(gap_length_ = (end_-start_)+1)

# gap summary
# instead of keeping gaps meeting min and max sizes
# we flag them accordingly
gap_data_ <- gap_data_[which(gap_data_$gap_length_ > 0),] %>%
  mutate(gap_filter_ = case_when(gap_length_ >= min_gap_length_ & gap_length_ <= max_gap_length_ ~ "target",
                             gap_length_ < min_gap_length_ ~ "shorter",
                             gap_length_ > max_gap_length_ ~ "longer",
                             .default = NA))

# summarise gap_data_
gap_data_sum_ <- gap_data_ %>%
  unique() %>%
  group_by(seqID,arrayID,filter_comb_,gap_filter_) %>%
  summarise(n_ = n(),
            gaps_mean_size_ = mean(gap_length_),
            sd_ = sd(gap_length_))

write.table(gap_data_,paste0(root_name_,"_gap_coordinates_all_.tsv"),quote = F, col.names = T, row.names = F, sep = "\t")
write.table(gap_data_sum_,paste0(root_name_,"_gap_coordinates_sum_.tsv"),quote = F, col.names = T, row.names = F, sep = "\t")

# keep gaps meeting min and max positions
gap_data_ <- gap_data_[which(gap_data_$gap_length_ >= min_gap_length_ & gap_data_$gap_length_ <= max_gap_length_),-c(ncol(gap_data_))]
write.table(gap_data_,paste0(root_name_,"_gap_coordinates_.tsv"),quote = F, col.names = T, row.names = F, sep = "\t")


