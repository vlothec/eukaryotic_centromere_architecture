# panel D
library(dplyr)
options(scipen=999,width=999)

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
intact_file_list_ <- list.files(path = path_to_data_, pattern = ".intact.gff3$",recursive = T)
arrays_file_list_ <- list.files(path = path_to_data_, pattern = "_arrays_expanded.csv$",recursive = T)
karyotype_file_list_ <- list.files(path = path_to_data_, pattern = "_karyotype_ext_$",recursive = T)

phyloclade_v_ <- c("Dicots" = "plant","Monocots" = "plant","Aves" = "chordate","Mammalia" = "chordate","Actinopterygii" = "chordate","Hymenoptera" = "invertebrate","Annelida" = "invertebrate","Coleoptera" = "invertebrate","Diptera" = "invertebrate","Hemiptera" = "invertebrate","Cnidaria" = "invertebrate","Lepidoptera" = "invertebrate","Neuroptera" = "invertebrate","Tunicata" = "chordate","Reptilia" = "chordate","Bryozoa" = "invertebrate","Mollusca" = "invertebrate","Plecoptera" = "invertebrate","Echinodermata" = "invertebrate","Blattodea" = "invertebrate","Ephemeroptera" = "invertebrate","Porifera" = "invertebrate","Arthropoda" = "invertebrate","Nemertea" = "invertebrate")


intact_data_all_ <- data.frame()
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

      # load karyo data
      karyotype_file_ <- karyotype_file_list_[grepl(gsub("\\..*","",dtolID_),karyotype_file_list_)]
      if(length(karyotype_file_) > 0){
        karyotype_data_ <- read.table(paste0(path_to_data_,"/",karyotype_file_),sep = "\t",header = F)
        colnames(karyotype_data_) <- c("ENA_","ChrOnly","length_")
      }
      
      # load intact data
      intact_file_ <- intact_file_list_[grepl(dtolID_,intact_file_list_)]
      if(length(intact_file_) > 1){
        intact_file_ <- intact_file_[!grepl("ChrOnly",intact_file_)]
      }
      if(length(intact_file_) > 0){
        intact_data_ <- data.table::fread(paste0(path_to_data_,"/",intact_file_), header = F,sep = "\t") %>% data.frame()
        colnames(intact_data_) <- c("seqID","method_anno_","TEanno_feature_","start_","end_","score_","strand_","phase_","attributes_")
        
        seqID_type_ <- which(intact_data_$seqID %in% karyotype_data_$ENA_) %>% length()
        if(seqID_type_ == 0){
          colnames(karyotype_data_)[2] <- "seqID"
          intact_data_ <- left_join(intact_data_,karyotype_data_[,c(1,2)])
          intact_data_$seqID <- intact_data_$ENA_
          intact_data_ <- intact_data_[,-ncol(intact_data_)]
        }
        
        intact_data_$centromere_status_ <- F
        intact_data_$arrayID_ <- NA_character_
        for(k_ in 1:nrow(arrays_data_)){
          print(k_)
          array_seqID_ <- arrays_data_$seqID[k_]
          array_start_ <- arrays_data_$start[k_]
          array_end_ <- arrays_data_$end[k_]
          array_arrayID_ <- arrays_data_$arrayID[k_]
          relevant_loci_ <- which(intact_data_$seqID %in% array_seqID_ & intact_data_$start_ >= array_start_ & intact_data_$end_ <= array_end_)
          if(length(relevant_loci_) > 0){
            intact_data_$centromere_status_[relevant_loci_] <- T
            intact_data_$arrayID_[relevant_loci_] <- array_arrayID_
          }
        }
        
        LTR_piden_data_ <- intact_data_ %>%
          {.[grepl("LTR_",.$TEanno_feature_),]} %>%
          {.[!grepl("non_LTR_",.$TEanno_feature_),]}
        
        if(nrow(LTR_piden_data_) > 0){
          
          LTR_piden_data_$ltr_identity <- LTR_piden_data_$attributes_ %>% gsub(".*dentity=([^;]*);.*","\\1",.) %>% as.numeric()
          
          LTR_piden_data_$centromere_status_ <- ifelse(LTR_piden_data_$centromere_status_ == TRUE, "cent","arm")
          LTR_piden_data_$centromere_status_  <- LTR_piden_data_$centromere_status_  %>%
            {factor(.,levels = unique(LTR_piden_data_$centromere_status_ ))}
          
          LTR_piden_data_$dtolID_ <- dtolID_
          LTR_piden_data_ <- left_join(LTR_piden_data_,DToL_info_[,c(7,3)],by = "dtolID_")
          
          phylotaxa_ <- unique(LTR_piden_data_$phylotaxa_)
          LTR_piden_data_$phyloclade_ <- phyloclade_v_[which(names(phyloclade_v_) == phylotaxa_)]
        }else{
          
          LTR_piden_data_ <- intact_data_[1,]
          
          LTR_piden_data_$ltr_identity <- NA
          LTR_piden_data_$TEanno_feature_ <- NA
          LTR_piden_data_$start_ <- NA
          LTR_piden_data_$end_ <- NA
          LTR_piden_data_$score_ <- NA
          LTR_piden_data_$strand_ <- NA
          LTR_piden_data_$attributes_ <- NA
          LTR_piden_data_$arrayID_ <- NA
          LTR_piden_data_$centromere_status_ <- NA
          LTR_piden_data_$dtolID_ <- dtolID_
          LTR_piden_data_ <- left_join(LTR_piden_data_,DToL_info_[,c(7,3)],by = "dtolID_")
          
          phylotaxa_ <- unique(LTR_piden_data_$phylotaxa_)
          LTR_piden_data_$phyloclade_ <- phyloclade_v_[which(names(phyloclade_v_) == phylotaxa_)]
        }
        intact_data_all_ <- rbind(intact_data_all_,LTR_piden_data_)
      }
    }
  }else{
    cat("Missing intact_data_ file\n")
    # quit(status = 1)
  }
}

nrow(intact_data_all_)
write.table(intact_data_all_,"LTR_piden_data_gw_aggregated_", col.names = T, row.names = F, sep = "\t", quote = F)

