# drGeuRiva1.hap1.1
{
  rm(list = ls())
  
  # get centromere repeat
  repeats_file_ <- "drGeuRiva1.hap1.1_repeats_filtered.csv"
  file_i_ <- repeats_file_
  if (file.info(file_i_)$size > 0) {
    repeats_data_ <- data.table::fread(repeats_file_, header = T,sep = ",") %>% data.frame()
    # for data generated with TRASH1
    if(length(colnames(repeats_data_)[grepl("new_class$",colnames(repeats_data_))]) == 0){
      repeats_data_$new_class <- repeats_data_$class
    }
    rm(file_i_)
  }
  # get arrays data
  arrays_file_ <- "drGeuRiva1.hap1.1_arrays_expanded.csv"
  file_i_ <- arrays_file_
  if (file.info(file_i_)$size > 0) {
    arrays_data_ <- read.table(arrays_file_, header = T,sep = ",")
    rm(file_i_)
  }
  # get TE data
  TEanno_file_ <- "drGeuRiva1.hap1.1_cent_feature_coord.txt"
  file_i_ <- TEanno_file_
  if (file.info(file_i_)$size > 0) {
    TEanno_data_ <- read.table(TEanno_file_, header = T,sep = "\t")
    rm(file_i_)
  }
  
  
  root_name_ <- arrays_file_ %>%
    {gsub(".*/","",.)} %>%
    {gsub("_arrays_expanded.csv","",.)} %>%
    {gsub("_fam_.*","",.)} %>%  
    {gsub(".fasta","",.)} %>%
    {gsub(".fa","",.)}
  print(root_name_)
  
  arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))
  
  filter_comb_i_ <- "TRUE_FALSE"
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
  array_repeats_data_$atr_13 <- (array_repeats_data_$end-array_repeats_data_$start)+1
  colnames(array_repeats_data_)[1] <- "chr"
  
  binwidth_ <- 10000
  array_repeats_data_sum_ <- array_repeats_data_ %>%
    mutate(diststart_ = findInterval(start,seq(1,max(start),length.out = max(start)/binwidth_),all.inside = T)) %>%
    group_by(chr,diststart_) %>%
    summarise(n_ = n(),
              total_length_ = sum(atr_13),
              prop_bin_ = sum(atr_13)/binwidth_)
  
  
  max_y_ <- max(TEanno_data_$atr_13,na.rm = T)
  pdf("drGeuRiva1.hap1.1_cent_feature_coord_plot_.pdf",height = 6,width = 9,onefile = T)
  p_ <- ggplot() +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade != "Athila"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13), color = "grey85",linewidth = .45) +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade == "Athila" & TEanno_data_$strand == "?"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13, color = strand),linewidth = .45) +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade == "Athila" & TEanno_data_$strand != "?"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13, color = strand),linewidth = .45) +
    geom_tile(data = array_repeats_data_sum_,aes(x = (diststart_*binwidth_)/1000000, y = -1500, fill = prop_bin_), color = NA,height = 1500) +
    scale_color_manual(values = c("+" = "red", "-" = "blue", "?" = "grey75")) +
    scale_fill_gradient(low = "lightgreen", high = "darkgreen", limits = range(array_repeats_data_sum_$prop_bin_)) +
    facet_wrap(vars(chr),ncol = 5, nrow = 5,scales = "free_x") +
    coord_cartesian(ylim = c(NA,max_y_))+
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, colour = "black", linewidth = 1),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = NA),
          axis.title = element_blank())
  print(p_)
  dev.off()
  
}

# drGeuUrba1.1
{
  rm(list = ls())
  
  # get centromere repeat
  repeats_file_ <- "drGeuUrba1.1_repeats_filtered.csv"
  repeats_data_ <- read.table(repeats_file_,sep = ",",header = T)
  cent_fam_ <- "159_6"
  repeats_data_ <- repeats_data_[which(repeats_data_$new_class %in% cent_fam_),]
  
  cent_file_ <- "drGeuUrba1_cent.coord"
  cent_data_ <- read.table(cent_file_,sep = "\t")
  colnames(cent_data_) <- c("chr","start","end","feature")
  colnames(cent_data_)[c(2,3)] <- c("start.cent","end.cent")
  
  cent_data_$length.cent <- mapply(function(x,y) (y-x)+1, cent_data_$start.cent, cent_data_$end.cent)
  
  repeats_data_$cent_status_ <- "OUT"
  for(i_ in 1:nrow(cent_data_)){
    seqID_ <- cent_data_$chr[i_]
    print(seqID_)
    start_i_ <- cent_data_$start.cent[i_]
    end_i_ <- cent_data_$end.cent[i_]
    idx_ <- which(repeats_data_$seqID == seqID_ & repeats_data_$start >= start_i_ & repeats_data_$end <= end_i_)
    if(length(idx_) > 0){
      repeats_data_$cent_status_[idx_] <- "IN"
    }
  }
  
  array_repeats_data_ <- repeats_data_[which(repeats_data_$cent_status_ == "IN"),]
  array_repeats_data_ <- array_repeats_data_[with(array_repeats_data_,order(seqID,start)),]
  array_repeats_data_$atr_13 <- (array_repeats_data_$end-array_repeats_data_$start)+1
  colnames(array_repeats_data_)[1] <- "chr"
  
  binwidth_ <- 10000
  array_repeats_data_sum_ <- array_repeats_data_ %>%
    mutate(diststart_ = findInterval(start,seq(1,max(start),length.out = max(start)/binwidth_),all.inside = T)) %>%
    group_by(chr,diststart_) %>%
    summarise(n_ = n(),
              total_length_ = sum(atr_13),
              prop_bin_ = sum(atr_13)/binwidth_)
  
  # get arrays data
  # get TE data
  TEanno_file_ <- "drGeuUrba1.1_cent_feature_coord_june2025_extended_rectified.txt"
  file_i_ <- TEanno_file_
  if (file.info(file_i_)$size > 0) {
    TEanno_data_ <- read.table(TEanno_file_, header = T,sep = "\t")
    rm(file_i_)
  }
  
  chr_list_ <- unique(TEanno_data_$chr)
  chr_list_ <- chr_list_[order(as.numeric(gsub("SUPER_","",chr_list_)))]
  
  TEanno_data_$chr <- TEanno_data_$chr %>%
    {factor(.,levels = chr_list_)}
  
  array_repeats_data_sum_$chr <- array_repeats_data_sum_$chr %>%
    {factor(.,levels = chr_list_)}
  
  TEanno_data_$atr_13 <- (TEanno_data_$end-TEanno_data_$start)+1
  max_y_ <- max(TEanno_data_$atr_13,na.rm = T)
  max_y_ <- 18107
  
  range(array_repeats_data_sum_$prop_bin_)
  
  pdf("drGeuUrba1.1_cent_feature_coord_plot_.pdf",height = 6,width = 9,onefile = T)
  p_ <- ggplot() +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade != "Athila"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13), color = "grey75",linewidth = .45) +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade == "Athila" & TEanno_data_$strand == "?"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13, color = strand)) +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade == "Athila" & TEanno_data_$strand != "?"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13, color = strand)) +
    geom_tile(data = array_repeats_data_sum_,aes(x = (diststart_*binwidth_)/1000000, y = -1500, fill = prop_bin_), color = NA,height = 1500) +
    scale_color_manual(values = c("+" = "red", "-" = "blue", "?" = "grey75")) +
    scale_fill_gradient(low = "lightgreen", high = "darkgreen", limits = c(0,1)) +
    facet_wrap(vars(chr),ncol = 5, nrow = 5,scales = "free_x") +
    coord_cartesian(ylim = c(NA,max_y_))+
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, colour = "black", linewidth = 1),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = NA),
          axis.title = element_blank())
  print(p_)
  dev.off()
  
  
  pdf("drGeuUrba1.1_cent_feature_coord_plot_.pdf",height = 6,width = 9,onefile = T)
  p_ <- ggplot() +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade != "Athila"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13), color = "grey75",linewidth = .45) +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade == "Athila" & TEanno_data_$strand == "?"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13, color = strand)) +
    geom_segment(data = TEanno_data_[which(TEanno_data_$Clade == "Athila" & TEanno_data_$strand != "?"),],aes(x = start/1000000, xend = start/1000000, y = 0, yend = atr_13, color = strand)) +
    geom_tile(data = TEanno_data_[which(TEanno_data_$Clade != "Athila"),],aes(x = start/1000000, y = -1500, fill = strand), color = NA, height = 1500) +
    scale_color_manual(values = c("+" = "red", "-" = "blue", "?" = "grey75")) +
    scale_fill_manual(values = c("+" = "transparent", "-" = "transparent", "?" = "transparent")) +
    facet_wrap(vars(chr),ncol = 5, nrow = 5,scales = "free_x") +
    coord_cartesian(ylim = c(NA,max_y_))+
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, colour = "black", linewidth = 1),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = NA),
          axis.title = element_blank())
  print(p_)
  dev.off()
}

