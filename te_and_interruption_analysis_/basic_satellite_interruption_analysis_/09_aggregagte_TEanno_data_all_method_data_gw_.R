
# panel C

bin_path_ <- "/path/to/bin_"
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

path_to_data_ <- "/path/to/data_"
TEanno_file_list_ <- list.files(path = path_to_data_, pattern = "_edta_filtered.csv.reassigned$",recursive = T)
arrays_file_list_ <- list.files(path = path_to_data_, pattern = "_arrays_expanded.csv$",recursive = T)

TEanno_data_all_ <- data.frame()
for(i_ in 1:length(arrays_file_list_)){
  arrays_file_ <- arrays_file_list_[i_]
  
  file_i_ <- paste0(path_to_data_,"/",arrays_file_)
  if (file.info(file_i_)$size > 0) {
    dtolID_ <- arrays_file_ %>%
      {gsub("(.*)_arrays_expanded.csv$","\\1",.)} %>%
      {gsub(".fasta","",.)} %>%
      {gsub(".fa","",.)} %>%
      {gsub("rosCan_S27_v1","rosCan",.)}
    
    print(dtolID_)
    rm(file_i_)
    
    arrays_data_ <- read.table(paste0(path_to_data_,"/",arrays_file_), header = T,sep = ",")
    if(ncol(arrays_data_) > 17){
      arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))
      filter_comb_i_ <- "TRUE_FALSE"
      # keep relevant arrays
      if(!is.null(filter_comb_i_)){
        arrays_data_ <- arrays_data_[which(arrays_data_$filter_comb_ %in% filter_comb_i_),]
      }
      
      TEanno_file_ <- TEanno_file_list_[grepl(dtolID_,TEanno_file_list_)]
      if(length(TEanno_file_) > 0){
        
        TEanno_data_ <- data.table::fread(paste0(path_to_data_,"/",TEanno_file_), header = F,sep = ",") %>% data.frame()
        header_check_ <- TEanno_data_[1,]
        if(length(grep("Method",header_check_)) > 0){
          TEanno_data_ <- TEanno_data_[-c(1),]
        }
        colnames(TEanno_data_) <- c("row","V1","V2","V3","V4","V5","V6","V7","V8","ID","Name","Classification","Sequence_ontology","Identity","Method","TSD","TIR","motif","tsd","oldV3","overlapping_bp","width","overlapping_percentage","reassignment_data_")
        colnames(TEanno_data_)[c(2,4,5,6)] <- c("seqID","TEanno_feature_","start_","end_")
        TEanno_data_$start_ <- as.numeric(TEanno_data_$start_)
        TEanno_data_$end_ <- as.numeric(TEanno_data_$end_)
        TEanno_data_$Identity <- as.numeric(TEanno_data_$Identity)

        if(length(which(TEanno_data_$TEanno_feature_ != "long_terminal_repeat" & TEanno_data_$TEanno_feature_ != "target_site_duplication")) > 0){
          TEanno_data_ <- TEanno_data_[which(TEanno_data_$TEanno_feature_ != "long_terminal_repeat" & TEanno_data_$TEanno_feature_ != "target_site_duplication"),]
        }
        
        TEanno_data_$centromere_status_ <- F
        TEanno_data_$arrayID_ <- NA_character_
        for(k_ in 1:nrow(arrays_data_)){
          print(k_)
          array_seqID_ <- arrays_data_$seqID[k_]
          array_start_ <- arrays_data_$start[k_]
          array_end_ <- arrays_data_$end[k_]
          array_arrayID_ <- arrays_data_$arrayID[k_]
          relevant_loci_ <- which(TEanno_data_$seqID %in% array_seqID_ & TEanno_data_$start_ >= array_start_ & TEanno_data_$end_ <= array_end_)
          TEanno_data_$centromere_status_[relevant_loci_] <- T
          TEanno_data_$arrayID_[relevant_loci_] <- array_arrayID_
        }
        
        TEanno_data_method <- TEanno_data_ %>%
          mutate(TEanno_length_ = (end_-start_)+1) %>%
          group_by(centromere_status_,Method) %>%
          summarise(n_ = n(),
                    total_TEanno_length_ = sum(TEanno_length_)) %>%
          group_by(centromere_status_) %>%
          mutate(prop_centromere_status_n_ = n_/sum(n_),
                 prop_centromere_status_length_ = total_TEanno_length_/sum(total_TEanno_length_))
        
        TEanno_data_method$centromere_status_ <- TEanno_data_method$centromere_status_ %>%
          {gsub(FALSE,"arm",.)} %>%
          {gsub(TRUE,"cent",.)}
        
        TEanno_data_method$dtolID_ <- dtolID_
        TEanno_data_method <- left_join(TEanno_data_method,DToL_info_[,c(7,3)],by = "dtolID_")
        
        TEanno_data_all_ <- rbind(TEanno_data_all_,TEanno_data_method)
        
      }
    }
  }else{
    cat("Missing TEanno_data_ file\n")
    quit(status = 1)
  }
}

TEanno_data_all_$phyloclade_ <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
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
                                                             .default = NA_character_),TEanno_data_all_$phylotaxa_)

write.table(TEanno_data_all_,"TEanno_data_all_method_data_gw_aggregated_", col.names = T, row.names = F, sep = "\t", quote = F)
