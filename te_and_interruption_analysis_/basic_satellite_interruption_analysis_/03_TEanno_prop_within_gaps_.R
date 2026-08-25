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
  
  TEanno_isec_gaps_file_ <- opts_$TEanno_isec_gaps_file_ %||% stop("TEanno_isec_gaps_file_ required")
}


# get input data
{
  # get repeats data
  file_i_ <- TEanno_isec_gaps_file_
  if (file.info(file_i_)$size > 0) {
    TEanno_isec_gaps_ <- read.delim(TEanno_isec_gaps_file_, header=FALSE, comment.char = "")
    colnames(TEanno_isec_gaps_) <- c("gap_seqID_","gap_start_","gap_end_","gap_length_id_","TEanno_seqID_","TEanno_start_","TEanno_end_","TEanno_feature_","TEanno_score_","TEanno_strand_","overlap_","filter_comb_")
    rm(file_i_)
  }else{
    cat("Missing TEanno_isec_gaps_ file\n")
    quit(status = 1)
  }
}

root_name_ <- TEanno_isec_gaps_file_ %>%
  {gsub("_TEanno_isec_gaps_","",.)}


TEanno_isec_gaps_$gap_length_ <- TEanno_isec_gaps_$gap_length_id_ %>%
  {gsub(".*_","",.)} %>%
  {as.numeric(.)}

# additional vars
TEanno_isec_gaps_ <- TEanno_isec_gaps_ %>%
  mutate(TEanno_length_ = (TEanno_end_-TEanno_start_),
         prop_gap_covered_raw_ = overlap_/gap_length_,
         prop_TEanno_covered_raw_ = overlap_/((TEanno_end_-TEanno_start_))) %>%
  group_by(gap_length_id_) %>%
  mutate(prop_gap_covered_ = sum(overlap_)/unique(gap_length_),
         prop_TEanno_covered_ = sum(overlap_)/sum(TEanno_length_),
         prop_gap_unannotated_ = ifelse((gap_length_ - sum(overlap_)) > 0, (gap_length_ - sum(overlap_))/gap_length_, 0))

write.table(TEanno_isec_gaps_, paste0(root_name_,"_TEanno_isec_gaps_parsed_"), row.names = F, col.names = T, sep = "\t", quote=F)
