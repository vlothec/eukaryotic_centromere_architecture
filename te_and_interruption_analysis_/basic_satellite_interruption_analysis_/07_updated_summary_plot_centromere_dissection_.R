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
  TEanno_file_ <- opts_$TEanno_file_ %||% stop("TEanno_file_ required")
  intact_file_ <- opts_$intact_file_ %||% stop("intact_file_ required")
  karyo_file_ <- opts_$karyo_file_ %||% stop("karyo_file_ required")
  gap_coordinates_file_ <- opts_$gap_coordinates_file_ %||% stop("gap_coordinates_file_ required")
  gap_coordinates_all_file_ <- opts_$gap_coordinates_all_file_ %||% stop("gap_coordinates_all_file_ required")
  centromere_dissection_data_gw_file_ <- opts_$centromere_dissection_data_gw_file_ %||% stop("centromere_dissection_data_gw_file_ required")
  centromere_dissection_data_seqID_file_ <- opts_$centromere_dissection_data_seqID_file_ %||% stop("centromere_dissection_data_seqID_file_ required")
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
  
  # get TEanno_file_ data
  file_i_ <- TEanno_file_
  if (file.info(file_i_)$size > 0) {
    TEanno_data_ <- data.table::fread(TEanno_file_, header = F,sep = ",") %>% data.frame()
    header_check_ <- TEanno_data_[1,]
    if(length(grep("Method",header_check_)) > 0){
      TEanno_data_ <- TEanno_data_[-c(1),]
    }
    colnames(TEanno_data_) <- c("row","V1","V2","V3","V4","V5","V6","V7","V8","ID","Name","Classification","Sequence_ontology","Identity","Method","TSD","TIR","motif","tsd","oldV3","overlapping_bp","width","overlapping_percentage","reassignment_data_")
    colnames(TEanno_data_)[c(2,4,5,6)] <- c("seqID","TEanno_feature_","start_","end_")
    TEanno_data_$start_ <- as.numeric(TEanno_data_$start_)
    TEanno_data_$end_ <- as.numeric(TEanno_data_$end_)
    TEanno_data_$Identity <- as.numeric(TEanno_data_$Identity)
    rm(file_i_)
  }else{
    cat("Missing TEanno_data_ file\n")
    quit(status = 1)
  }
  
  # get intact_file_ data
  file_i_ <- intact_file_
  if (file.info(file_i_)$size > 0) {
    intact_data_ <- read.table(intact_file_, header = F, sep = "\t")
    colnames(intact_data_) <- c("seqID","method_anno_","TEanno_feature_","start_","end_","score_","strand_","phase_","attributes_")
    rm(file_i_)
  }else{
    cat("Missing intact_data_ file\n")
    quit(status = 1)
  }
  
  # get karyo_file_ data
  file_i_ <- karyo_file_
  if (file.info(file_i_)$size > 0) {
    karyo_ <- read.table(karyo_file_, header = F, sep = "\t")
    colnames(karyo_) <- c("ENA_","seqID_","seq_length_")
    rm(file_i_)
  }else{
    cat("Missing karyo_ file\n")
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
  
  # get gap_coordinates_all_file_ data
  file_i_ <- gap_coordinates_all_file_
  if (file.info(file_i_)$size > 0) {
    gap_coordinates_all_ <- read.table(gap_coordinates_all_file_, header = T,sep = "\t")
    rm(file_i_)
  }else{
    cat("Missing gap_coordinates_all_ file\n")
    quit(status = 1)
  }
  
  # get centromere_dissection_data_gw_file_ data
  file_i_ <- centromere_dissection_data_gw_file_
  if (file.info(file_i_)$size > 0) {
    centromere_dissection_data_gw_ <- read.table(centromere_dissection_data_gw_file_,header = T,sep = "\t")
    
    idx_ <- which(centromere_dissection_data_gw_$TEclass_ == "repeat_fragment" & centromere_dissection_data_gw_$total_length_ > 0)
    if(length(idx_)>0){
      centromere_dissection_data_gw_$TEclass_[idx_] <- "repeat_region"
    }
    
    centromere_dissection_data_gw_$TEclass_[which(centromere_dissection_data_gw_$TEclass_ == "Class_I_LTR")] <- "Class_I_LTR"
    centromere_dissection_data_gw_$TEclass_[which(centromere_dissection_data_gw_$TEclass_ == "Class_II_DNA_element")] <- "Class_II_DNA"
    centromere_dissection_data_gw_$TEclass_[which(centromere_dissection_data_gw_$TEclass_ == "Class_I_non_LTR_TPRT")] <- "Class_I_TPRT"
    centromere_dissection_data_gw_$TEclass_[which(centromere_dissection_data_gw_$TEclass_ == "Class_I_non_LTR_other")] <- "Class_I_Other"
    
    rm(file_i_)
  }else{
    cat("Missing centromere_dissection_data_gw_ file\n")
    quit(status = 1)
  }
  
  # get centromere_dissection_data_seqID_file_ data
  file_i_ <- centromere_dissection_data_seqID_file_
  if (file.info(file_i_)$size > 0) {
    centromere_dissection_data_seqID_ <- read.table(centromere_dissection_data_seqID_file_,header = T,sep = "\t")
    rm(file_i_)
  }else{
    cat("Missing centromere_dissection_data_seqID_ file\n")
    quit(status = 1)
  }
  
  TE_class_list_updated_ <- read.table(paste0(bin_path_,"/TE_class_list_updated_"),sep = "\t",header = F)
  colnames(TE_class_list_updated_) <- c("TEanno_cls_","old_TEanno_cls_","TEanno_feature_")
  
}

dtol_info_ <- read.table(paste0(bin_path_,"/DToL_info_.txt"),header = F,sep = "\t")
colnames(dtol_info_) <- c("batch_","assembly_","clade_","genus_","species_","genome_size_")
dtol_info_ <- dtol_info_ %>%
  mutate(binomial_name_ = paste0(gsub("([A-Za-z]).*","\\1",genus_),". ",species_),
         long_binomial_name_ = paste0(genus_," ",species_),
         dtolID_ = gsub(".fa.*","",assembly_))

root_name_ <- centromere_dissection_data_gw_file_ %>%
  {gsub(".*/","",.)} %>%
  {gsub("_centromere_dissection_data_gw_","",.)}
print(root_name_)

long_binomial_name_ <- dtol_info_$long_binomial_name_[which(dtol_info_$dtolID_ %in% root_name_)]

arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))

# keep relevant arrays
if(!is.null(filter_comb_i_)){
  arrays_data_ <- arrays_data_[which(arrays_data_$filter_comb_ %in% filter_comb_i_),]
}

arrays_data_$TRASH_class <- as.character(arrays_data_$TRASH_class)

# annotate centromere boundaries split files
TEanno_data_$centromere_status_ <- F
TEanno_data_$arrayID_ <- NA_character_
for(i_ in 1:nrow(arrays_data_)){
  print(i_)
  array_seqID_ <- arrays_data_$seqID[i_]
  array_start_ <- arrays_data_$start[i_]
  array_end_ <- arrays_data_$end[i_]
  array_arrayID_ <- arrays_data_$arrayID[i_]
  rm(relevant_loci_)
  relevant_loci_ <- which(TEanno_data_$seqID %in% array_seqID_ & TEanno_data_$start_ >= array_start_ & TEanno_data_$end_ <= array_end_)
  TEanno_data_$centromere_status_[relevant_loci_] <- T
  TEanno_data_$arrayID_[relevant_loci_] <- array_arrayID_
}

idx_ <- which(gap_coordinates_all_$gap_filter_ != "target")
if(length(idx_) > 0){
  for(i_ in 1:length(idx_)){
  # print(i_)
  gap_seqID_ <- gap_coordinates_all_$seqID[idx_[i_]]
  gap_start_ <- gap_coordinates_all_$start_[idx_[i_]]
  gap_end_ <- gap_coordinates_all_$end_[idx_[i_]]
  rm(relevant_loci_)
  relevant_loci_ <- which(TEanno_data_$seqID %in% gap_seqID_ & TEanno_data_$start_ >= gap_start_ & TEanno_data_$end_ <= gap_end_)
  TEanno_data_$centromere_status_[relevant_loci_] <- F
  }
}

# annotate centromere boundaries intact files
intact_data_$centromere_status_ <- F
intact_data_$arrayID_ <- NA_character_

seqID_list_ <- intact_data_$seqID %>% unique()
seqID_test_ <- length(which(karyo_$ENA_ %in% seqID_list_))

for(i_ in 1:nrow(arrays_data_)){
    if(seqID_test_ > 0){
      print(i_)
      array_seqID_ <- arrays_data_$seqID[i_]
      array_start_ <- arrays_data_$start[i_]
      array_end_ <- arrays_data_$end[i_]
      array_arrayID_ <- arrays_data_$arrayID[i_]
      rm(relevant_loci_)
      relevant_loci_ <- which(intact_data_$seqID %in% array_seqID_ & intact_data_$start_ >= array_start_ & intact_data_$end_ <= array_end_)
      intact_data_$centromere_status_[relevant_loci_] <- T
      intact_data_$arrayID_[relevant_loci_] <- array_arrayID_
    }else{
      print(i_)
      array_ENA_seqID_ <- arrays_data_$seqID[i_]
      array_seqID_ <- karyo_$seqID_[which(karyo_$ENA_ %in% array_ENA_seqID_)]
      array_start_ <- arrays_data_$start[i_]
      array_end_ <- arrays_data_$end[i_]
      array_arrayID_ <- arrays_data_$arrayID[i_]
      rm(relevant_loci_)
      relevant_loci_ <- which(intact_data_$seqID %in% array_seqID_ & intact_data_$start_ >= array_start_ & intact_data_$end_ <= array_end_)
      intact_data_$centromere_status_[relevant_loci_] <- T
      intact_data_$arrayID_[relevant_loci_] <- array_arrayID_
    }
}

idx_ <- which(gap_coordinates_all_$gap_filter_ != "target")
if(length(idx_) > 0){
  for(i_ in 1:length(idx_)){
    if(seqID_test_ > 0){
      # print(i_)
      gap_seqID_ <- gap_coordinates_all_$seqID[idx_[i_]]
      gap_start_ <- gap_coordinates_all_$start_[idx_[i_]]
      gap_end_ <- gap_coordinates_all_$end_[idx_[i_]]
      rm(relevant_loci_)
      relevant_loci_ <- which(intact_data_$seqID %in% gap_seqID_ & intact_data_$start_ >= gap_start_ & intact_data_$end_ <= gap_end_)
      intact_data_$centromere_status_[relevant_loci_] <- F
    }else{
      # print(i_)
      gap_ENA_seqID_ <- gap_coordinates_all_$seqID[idx_[i_]]
      gap_seqID_ <- karyo_$seqID_[which(karyo_$ENA_ %in% gap_ENA_seqID_)]
      gap_start_ <- gap_coordinates_all_$start_[idx_[i_]]
      gap_end_ <- gap_coordinates_all_$end_[idx_[i_]]
      rm(relevant_loci_)
      relevant_loci_ <- which(intact_data_$seqID %in% gap_seqID_ & intact_data_$start_ >= gap_start_ & intact_data_$end_ <= gap_end_)
      intact_data_$centromere_status_[relevant_loci_] <- F
    }
  }
}


# centromere space overview
{
  centromere_dissection_data_seqID_ <- centromere_dissection_data_seqID_ %>%
    rowwise() %>%
    mutate(var_ = paste0(step_,"_",feature_))
  
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "1_centromere_nt_")] <- "1_cent"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "2_TR_nt_")] <- "2_centTR"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "2_gap_nt_")] <- "3_gaps"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "3_TE_nt_")] <- "4_TEanno"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "3_otherTR_nt_")] <- "5_otherTR"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "3_unk_nt_")] <- "6_unk"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "4_TE_nt_")] <- "7_TEanno"
  centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$var_ == "5_TE_nt_")] <- "8_TEanno"
}


# define grid structure
mat_ <- rbind(c(1, 2, 2, 2, 3),
              c(4, 4, 6, 6, 6),
              c(5, 5, 6, 6, 6),
              c(7, 7, 7, 7, 7),
              c(8, 8, 8, 8, 8))
seqID_v_ <- arrays_data_$seqID %>% unique()
seqID_panel_ <- seq(9,length(seqID_v_)+9,1)
for(k_ in 1:length(seqID_panel_)){
  mat_k_ <- c(seqID_panel_[k_],seqID_panel_[k_],seqID_panel_[k_],seqID_panel_[k_],seqID_panel_[k_])
  mat_ <- rbind(mat_,mat_k_)
}
last_null_plot_ <- seqID_panel_[length(seqID_panel_)]+1
mat_ <- rbind(mat_,
              c(last_null_plot_, last_null_plot_, last_null_plot_, last_null_plot_, last_null_plot_))

total_width_cm_ <- 10*2.54

total_fixed_height_ <- 8
total_var_height_ <- ((length(seqID_v_)+1)*32)/42
total_height_ <-  total_fixed_height_ + total_var_height_



total_height_cm_ <- total_height_*2.54
fixed_panel_height_cm_ <- 8*2.54

var_panel_height_cm_ <- (total_height_cm_-fixed_panel_height_cm_)
var_panel_height_cm_ <- var_panel_height_cm_/(length(seqID_v_)+1)


pdf(paste0("new_layout_",root_name_,"_centromere_dissection_overview_.pdf"),width = 10, height = total_height_,onefile = T)

# abs
layout(mat = mat_,
       heights = c(lcm(fixed_panel_height_cm_*0.425),
                   lcm(fixed_panel_height_cm_*0.210),
                   lcm(fixed_panel_height_cm_*0.235),
                   lcm(fixed_panel_height_cm_*0.050),
                   lcm(var_panel_height_cm_),
                   sapply(replicate(length(seqID_v_),var_panel_height_cm_),lcm),
                   lcm(fixed_panel_height_cm_*0.075)),
       widths = c(lcm(total_width_cm_*0.3),
                  lcm(total_width_cm_*0.1),
                  lcm(total_width_cm_*0.1),
                  lcm(total_width_cm_*0.2),
                  lcm(total_width_cm_*0.2))
       )


# common var all plots
col_v_ <- c("TRASH" = "#D83067","all_TEs" = "#487CE3", "unknown" = "#DDDDDD")
cex_main_ <- .75

# empty panel
# par(mar = c(0,0,0,0))
# plot.new()

# centromere main numbers
{
  n_karyo_ <- nrow(karyo_)
  genome_size_ <- karyo_$seq_length_ %>% sum()
  cent_fam_ <- arrays_data_$TRASH_class %>% unique()
  n_cent_ <- arrays_data_$seqID %>% unique() %>% length()
  
  total_n_gaps_ <- nrow(gap_coordinates_all_)
  surveyed_n_gaps_ <- nrow(gap_coordinates_)
  omited_n_gaps_ <- nrow(gap_coordinates_all_[which(gap_coordinates_all_$gap_filter_ != "target"),])
  
  total_gap_space_ <- sum(gap_coordinates_all_$gap_length_)
  surveyed_gap_space_ <- sum(gap_coordinates_$gap_length_)
  omitted_gap_space_ <- sum(gap_coordinates_all_$gap_length_[which(gap_coordinates_all_$gap_filter_ != "target")])
  
  
  plot(c(0:6.5),
       c(0:6.5),
       type = "n", axes = FALSE, xlab = "", ylab = "",
       main = long_binomial_name_,
       adj = 0,
       font = 4)
  
  # text(0, 6, paste0(long_binomial_name_),
  #      adj = 0,font = 2)
  text(0, 6, paste0(root_name_),
       adj = 0,font = 1)
  text(0, 5.5, paste0(genome_size_, " bp"),
       adj = 0)
  
  # text(0, 4.5, paste0(cent_fam_, " satellite family"),
  #      adj = 0)
  text(0, 4.5, paste0(cent_fam_),
       adj = 0)
  
  longer_label_width_ <- paste0("satellite-based chromosomes                 ") %>% strwidth()
  tab_space_ <- 0
  text(0, 4, paste0("satellite-based chromosomes"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 4, paste0(n_cent_,"/",n_karyo_),
       adj = 1)
  
  # gaps
  text(0, 3, paste0("num of interruptions"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 3, paste0(total_n_gaps_),
       adj = 1)
  
  text(0, 2.5, paste0("surveyed"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 2.5, paste0(surveyed_n_gaps_),
       adj = 1)
  
  text(0, 2, paste0("omitted"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 2, paste0(omited_n_gaps_),
       adj = 1)
  
  # space
  text(0, 1, paste0("interruption space"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 1, paste0(total_gap_space_),
       adj = 1)
  
  text(0, 0.5, paste0("surveyed"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 0.5, paste0(surveyed_gap_space_),
       adj = 1)
  
  text(0, 0, paste0("omitted"),
       adj = 0)
  text(0+longer_label_width_+tab_space_, 0, paste0(omitted_gap_space_),
       adj = 1)
  # box()
  
  # label
  text(grconvertX(0.005, from = "nfc", to = "user"),
       grconvertY(0.985, from = "nfc", to = "user"), 
       "a", xpd = NA, cex = 1.4, font = 2)
}

# centromere space overview: dots
par(mar = c(6, 4, 3, 4)-0.25) # bottom, left, top, right
{
  pch_ <- 21
  max_y_ <- max(centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$feature_ != "centromere_nt_")])
  range_x_ <- range(centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")])
  if(n_digits_(range_x_[1]) == n_digits_(range_x_[2])){
    max_x_ <- range_x_[2]
    min_x_ <- 1*10^(n_digits_(range_x_[1])-1)
  }else{
    max_x_ <- range_x_[2]
    min_x_ <- range_x_[1]
  }
  plot(centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")],
       replicate(length(which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")),max_y_),
       type = "n",
       ann = F,
       ylim = c(10,max_y_),
       xlim = c(min_x_,max_x_),
       log = "xy",
       xaxt = "n",
       yaxt = "n")
  abline(v = centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")],
         ann = F,
         col = "grey95",
         lty = 1,
         xaxt = "n",
         yaxt = "n")
  points(centromere_dissection_data_seqID_$centromere_size_[which(centromere_dissection_data_seqID_$var_ == "2_centTR")],
         centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$var_ == "2_centTR")],
         pch = pch_,
         col = col_v_[1],
         ann = F,
         xaxt = "n",
         yaxt = "n")
  points(centromere_dissection_data_seqID_$centromere_size_[which(centromere_dissection_data_seqID_$var_ == "5_otherTR")],
         centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$var_ == "5_otherTR")],
         pch = pch_,
         col = "#34c975",
         ann = F,
         xaxt = "n",
         yaxt = "n")
  points(centromere_dissection_data_seqID_$centromere_size_[which(centromere_dissection_data_seqID_$var_ == "4_TEanno")],
         centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$var_ == "4_TEanno")],
         pch = pch_,
         col = col_v_[2],
         ann = F,
         xaxt = "n",
         yaxt = "n")
  points(centromere_dissection_data_seqID_$centromere_size_[which(centromere_dissection_data_seqID_$var_ == "6_unk")],
         centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$var_ == "6_unk")],
         pch = pch_,
         col = "grey30",
         ann = F,
         xaxt = "n",
         yaxt = "n")
  boxplot(centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")],
          at = 20,
          ann = F,
          horizontal = T,
          frame.plot = F,
          boxwex = 0.75,
          pch = 16,
          lty = 1,
          add = T,
          xaxt = "n",
          yaxt = "n")
  
  box()
  exp_ <- n_digits_(max_y_)
  
  # Add custom axis
  axis(1, at = 10^(0:exp_), labels = 10^(0:exp_))  # major ticks
  axis(2, at = 10^(0:exp_), labels = 10^(0:exp_))  # major ticks
  
  # minor ticks
  range_x_ <- c(min_x_,max_x_)
  minor_ticks_x_ <- unlist(lapply(0:6, function(e) (1:9) * 10^e))
  minor_ticks_x_ <- minor_ticks_x_[which(minor_ticks_x_ >= range_x_[1] & minor_ticks_x_ <= range_x_[2])]
  axis(1, at = minor_ticks_x_, labels = FALSE, tcl = -0.2, col = "black")
  
  range_y_ <- range(c(10, max_y_))
  minor_ticks_y_ <- unlist(lapply(0:6, function(e) (1:9) * 10^e))
  minor_ticks_y_ <- minor_ticks_y_[which(minor_ticks_y_ >= range_y_[1] & minor_ticks_y_ <= range_y_[2])]
  axis(2, at = minor_ticks_y_, labels = FALSE, tcl = -0.2, col = "black")
  
  # axis title
  # title(xlab = "Centromere Length - log(bp)",
        # ylab = "Feature Length - log(bp)")
  mtext("Centromere Length - log(bp)", side = 1, line = 3, cex = 0.75)
  mtext("Feature Length - log(bp)", side = 2, line = 3, cex = 0.75)
  
  # label
  text(grconvertX(0.005, from = "nfc", to = "user"),
       grconvertY(0.985, from = "nfc", to = "user"), 
       "b", xpd = NA, cex = 1.4, font = 2)
  
}

# centromere space overview: boxplot
# par(mar = c(6,6,6,6)) # bottom, left, top, right
par(mar = c(6, 4, 3, 4)-0.25) # bottom, left, top, right
{
  x_labels_ <- centromere_dissection_data_seqID_$var_[which(centromere_dissection_data_seqID_$step_ <= 3 & centromere_dissection_data_seqID_$step_ > 1 & centromere_dissection_data_seqID_$var_ != "3_gaps")]
  
  x_labels_[which(x_labels_ == "2_centTR")] <- "Satellites"
  x_labels_[which(x_labels_ == "4_TEanno")] <- "Transposons"
  x_labels_[which(x_labels_ == "5_otherTR")] <- "Other Satellites"
  x_labels_[which(x_labels_ == "6_unk")] <- "Unannotated"
  
  x_labels_ <- x_labels_ %>%
    {factor(.,levels = c("Satellites","Transposons","Other Satellites","Unannotated"))}
  
  boxplot(centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$step_ <= 3 & centromere_dissection_data_seqID_$step_ > 1 & centromere_dissection_data_seqID_$var_ != "3_gaps")] ~ x_labels_,
          horizontal = F,
          ylim = c(10,max_y_),
          log = "y",
          las=2,
          ann = F,
          col = c(col_v_[1],col_v_[2],"#34c975","grey30"),
          frame.plot = F,
          boxwex = 0.5,
          pch = 16,
          lty = 1,
          yaxt = "n")
  
  # Add custom axis
  axis(2, at = 10^(0:6), labels = 10^(0:6))  # major ticks
  
  # minor ticks
  range_y_ <- range(c(10, max_y_))
  minor_ticks_y_ <- unlist(lapply(0:6, function(e) (1:9) * 10^e))
  minor_ticks_y_ <- minor_ticks_y_[which(minor_ticks_y_ >= range_y_[1] & minor_ticks_y_ <= range_y_[2])]
  axis(2, at = minor_ticks_y_, labels = FALSE, tcl = -0.2, col = "black")
  
  box()
  
  # axis title
  # title(ylab = "Feature Length - log(bp)")
  mtext("Feature Length - log(bp)", side = 2, line = 3, cex = 0.75)
  
  # label
  # label
  text(grconvertX(0.005, from = "nfc", to = "user"),
       grconvertY(0.985, from = "nfc", to = "user"), 
       "c", xpd = NA, cex = 1.4, font = 2)
  
}

# struc/homology 
{
  par(mar = c(6, 4, 3, 4)-0.25) # bottom, left, top, right
  line_val_ <- 1.35
  
  # method
  TEanno_data_method <- TEanno_data_ %>%
    mutate(TEanno_length_ = (end_-start_)+1) %>%
    group_by(centromere_status_,Method) %>%
    summarise(n_ = n(),
              total_TEanno_length_ = sum(TEanno_length_)) %>%
    group_by(centromere_status_) %>%
    mutate(prop_centromere_status_n_ = n_/sum(n_),
           prop_centromere_status_length_ = total_TEanno_length_/sum(total_TEanno_length_))
  
  all_ <- expand.grid(centromere_status_ = c(T,F),
                      Method = c("homology","structural"))
  TEanno_data_method <- full_join(TEanno_data_method,all_)
  TEanno_data_method$n_[is.na(TEanno_data_method$n_ == T)] <- 0
  TEanno_data_method$total_TEanno_length_[is.na(TEanno_data_method$total_TEanno_length_ == T)] <- 0
  TEanno_data_method$prop_centromere_status_n_[is.na(TEanno_data_method$prop_centromere_status_n_ == T)] <- 0
  TEanno_data_method$prop_centromere_status_length_[is.na(TEanno_data_method$prop_centromere_status_length_ == T)] <- 0
  
  TEanno_data_mat_ <- TEanno_data_method$prop_centromere_status_length_[which(TEanno_data_method$Method == "homology")]
  col_names_ <- TEanno_data_method$centromere_status_[which(TEanno_data_method$Method == "homology")] %>%
    {gsub(FALSE,"arm",.)} %>%
    {gsub(TRUE,"cent",.)}
  names(TEanno_data_mat_) <- col_names_
  TEanno_data_mat_ <- rbind(TEanno_data_mat_,
                            TEanno_data_method$prop_centromere_status_length_[which(TEanno_data_method$Method != "homology")])
  rownames(TEanno_data_mat_) <- c("homology","structural")
  
  offset_struct_ <- 0.035
  bp_ <- barplot(TEanno_data_mat_,
                 names.arg = colnames(TEanno_data_mat_),
                 col = c("#B5D6E4","#000A83"),
                 border = NA,
                 xlim = c(0,1+(offset_struct_*10)),
                 beside = F,
                 horiz = T,
                 xaxt = "n",
                 yaxt = "n",
                 main = " ",
                 width = 0.5,
                 space = 1,
                 las = 1)
  
  n_ <- round(TEanno_data_method$total_TEanno_length_[which(TEanno_data_method$Method == "homology")]/1000000,digits = 2)
  text(
    y = bp_+0.05,
    # x = 1 + offset_struct_,
    x = 0.10,
    adj = 0,
    labels = paste0(n_, " Mb fragment"),
    cex = 0.65
  )
  n_ <- round(TEanno_data_method$total_TEanno_length_[which(TEanno_data_method$Method != "homology")]/1000000,digits = 2)
  text(
    y = bp_+0.05,
    # x = 1 + offset_struct_,
    x = 0.10,
    adj = 0,
    labels = paste0("\n\n", n_, " Mb intact"),
    cex = 0.65
  )
  
  axis(
    side = 1,
    at = seq(0,1,.25),
    labels = seq(0,1,.25),
    line = line_val_
  )
  
  axis(
    side = 2,
    at = bp_,
    labels = colnames(TEanno_data_mat_),
    lwd = 0
  )
  
  legend(
    grconvertX(0.0875, from = "nfc", to = "user"),
    grconvertY(0.9800, from = "nfc", to = "user"), 
    legend = c("fragment","intact"),
    col = c("#B5D6E4","#000A83"),
    pch = 15,
    bty = "n",
    pt.cex = 1.25, xpd = NA,
    horiz = T
  )
  
  # axis title
  # title(xlab = "transposon type (%)")
  mtext("Transposon Type (%)", at = c(0.5),side = 1, line = 4, cex = 0.75)
  
  # label
  # label
  text(grconvertX(0.005, from = "nfc", to = "user"),
       grconvertY(0.985, from = "nfc", to = "user"), 
       "d", xpd = NA, cex = 1.4, font = 2)
  
}

# pident
{
  par(mar = c(2, 8, 3, 12)-.25)
  line_val_ <- .5
  LTR_piden_data_ <- intact_data_ %>%
    {.[grepl("LTR_",.$TEanno_feature_),]} %>%
    {.[!grepl("non_LTR_",.$TEanno_feature_),]}
  
  LTR_piden_data_$ltr_identity <- LTR_piden_data_$attributes_ %>% gsub(".*ltr_identity=([^;]*);.*","\\1",.) %>% as.numeric()
  
  LTR_piden_data_$centromere_status_ <- ifelse(LTR_piden_data_$centromere_status_ == TRUE, "cent","arm")
  LTR_piden_data_$centromere_status_  <- LTR_piden_data_$centromere_status_  %>%
    {factor(.,levels = unique(LTR_piden_data_$centromere_status_ ))}
  
  if(length(LTR_piden_data_$ltr_identity[which(is.na(LTR_piden_data_$ltr_identity) == F)]) > 0){
    min_pident_val_ <- min(LTR_piden_data_$ltr_identity,na.rm = T)
    min_pident_val_ <- 0.85
    offset_ <- 0.0025
    bp_ <- boxplot(LTR_piden_data_$ltr_identity ~ LTR_piden_data_$centromere_status_,
                   horizontal = F,
                   ylab = paste0(" "),
                   xlab = "",
                   col = c("snow","snow"),
                   ylim = c(min_pident_val_,max(LTR_piden_data_$ltr_identity,na.rm = T)+(offset_*5)),
                   xaxt = "n",
                   yaxt = "n",
                   main = " ",
                   frame.plot = F,
                   boxwex = 0.5,
                   pch = 16,
                   lty = 1)
    
    n_ <- LTR_piden_data_$centromere_status_ %>% table
    text(
      x = 1:length(n_),
      # x = bp_$stats[5,] + offset_,
      y = min_pident_val_,
      labels = paste0("n: ", n_)
    )
    
    min_floor_ <- (floor(min(LTR_piden_data_$ltr_identity,na.rm = T)*100))/100
    min_floor_ <- (floor(min_pident_val_*100))/100
    axis(
      side = 2,
      at = seq(min_floor_, 1, by = (1-min_floor_)/4),
      labels = paste0(round(seq(min_floor_,1, by = (1-min_floor_)/4),digits = 4)*100,"%"),
      line = line_val_,
      las = 1
    )
    
    axis(
      1,
      at = unique(LTR_piden_data_$centromere_status_),
      labels = unique(LTR_piden_data_$centromere_status_),
      lwd = 0
    )
    # axis title
    # title(ylab = "LTR Identity")
    mtext("LTR Identity",side = 2, line = 6, cex = 0.75)
    
    # label
    # label
    text(grconvertX(0.005, from = "nfc", to = "user"),
         grconvertY(0.985, from = "nfc", to = "user"), 
         "e", xpd = NA, cex = 1.4, font = 2)
    
  }else{
    plot.new()
  }
  
}

# TEclasses
{
  par(mar = c(2, 4, 3, 4)-0.25) # bottom, left, top, right
  # par(mar = c(0, 0, 0, 0))
  line_val_ <- 1
  
  TEanno_cls_ <- TE_class_list_updated_$TEanno_cls_ %>% unique()
  TEanno_cls_[which(TEanno_cls_ == "Class_I_LTR")] <- "Class_I_LTR"
  TEanno_cls_[which(TEanno_cls_ == "Class_II_DNA_element")] <- "Class_II_DNA"
  TEanno_cls_[which(TEanno_cls_ == "Class_I_non_LTR_TPRT")] <- "Class_I_TPRT"
  TEanno_cls_[which(TEanno_cls_ == "Class_I_non_LTR_other")] <- "Class_I_Other"
  
  TEanno_cls_ <- TEanno_cls_[!grepl("repeat_fragment",TEanno_cls_)]
  
  TEanno_cls_ <- TEanno_cls_ %>%
    {gsub("_"," ",.)}
  
  centromere_dissection_data_gw_$TEclass_[which(centromere_dissection_data_gw_$step_ == 4)] <- centromere_dissection_data_gw_$TEclass_[which(centromere_dissection_data_gw_$step_ == 4)] %>%
    {gsub("_"," ",.)}
  
  # names_arg_ <- c("centromere"," ","TR","gaps","  ","unknown","otherTR","TEanno","   ",TEanno_cls_)
  names_arg_ <- c("centromere"," ","Satellites","Interruptions","  ","Unannotated","Other Satellites","Transposons","   ",TEanno_cls_)
  names_arg_ <- rev(names_arg_)
  
  # cent size
  x_cent_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "centromere_nt_")]
  x_cent_ <- ifelse(length(x_cent_) == 0, 0, x_cent_)
  # cent accounted by gaps
  x_gap_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "gap_nt_")]
  x_gap_ <- ifelse(length(x_gap_) == 0, 0, x_gap_)
  # cent accounted by TR
  x_TR_ <- x_cent_-x_gap_
  x_TR_ <- ifelse(length(x_TR_) == 0, 0, x_TR_)
  # cent accounted by TEanno
  x_all_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$TEclass_ == "all")]
  x_all_ <- ifelse(length(x_all_) == 0, 0, x_all_)
  # cent non accounted
  x_unknown_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "unk_nt_")]
  x_unknown_ <- ifelse(length(x_unknown_) == 0, 0, x_unknown_)
  
  # cent accounted by TEanno
  x_otherTR_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "otherTR_nt_")]
  x_otherTR_ <- ifelse(length(x_otherTR_) == 0, 0, x_otherTR_)
  
  max_x_ <- max(centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "centromere_nt_")],na.rm = T)
  min_x_ <- 0
  
  # all categories y axis
  bp_ <- barplot(replicate(length(names_arg_),x_cent_),
                 names.arg = names_arg_,
                 horiz = T,
                 beside = F,
                 main = "",
                 width = 0.65,
                 space = 1,
                 col = c(replicate(length(TEanno_cls_),col_v_[3]),NA,col_v_[3],col_v_[3],col_v_[3],NA,col_v_[3],col_v_[3],NA,col_v_[3]),
                 border = NA,
                 xlim = c(min_x_,max_x_),
                 las = 1,
                 axisnames = T)
  
  # centromere
  barplot(replicate(length(names_arg_),x_cent_),
          names.arg = names_arg_,
          horiz = T,
          beside = F,
          main = "",
          width = 0.65,
          space = 1,
          border = NA,
          col = c(replicate(length(names_arg_)-1,NA),"black"),
          xlim = c(min_x_,max_x_),
          las = 1,
          axisnames = T,
          xaxt = "n",
          yaxt = "n",
          add = T)
  # TR
  barplot(replicate(length(names_arg_),x_TR_),
          names.arg = names_arg_,
          horiz = T,
          beside = F,
          main = "",
          width = 0.65,
          space = 1,
          col = c(replicate(length(names_arg_)-3,NA),col_v_[1],NA,NA),
          border = NA,
          xlim = c(min_x_,max_x_),
          las = 1,
          axisnames = T,
          xaxt = "n",
          yaxt = "n",
          add = T)
  # gaps
  barplot(replicate(length(names_arg_),x_gap_),
          names.arg = names_arg_,
          horiz = T,
          beside = F,
          main = "",
          width = 0.65,
          space = 1,
          col = c(replicate(length(names_arg_)-4,NA),"grey65",NA,NA,NA),
          border = NA,
          xlim = c(min_x_,max_x_),
          las = 1,
          axisnames = T,
          xaxt = "n",
          yaxt = "n",
          add = T)
  
  # unknown
  barplot(replicate(length(names_arg_),x_unknown_),
          names.arg = names_arg_,
          horiz = T,
          beside = F,
          main = "",
          width = 0.65,
          space = 1,
          col = c(replicate(length(names_arg_)-6,NA),"grey65",NA,NA,NA,NA,NA),
          border = NA,
          xlim = c(min_x_,max_x_),
          las = 1,
          axisnames = T,
          xaxt = "n",
          yaxt = "n",
          add = T)
  
  # other TR
  barplot(replicate(length(names_arg_),x_otherTR_),
          names.arg = names_arg_,
          horiz = T,
          beside = F,
          main = "",
          width = 0.65,
          space = 1,
          col = c(replicate(length(names_arg_)-7,NA),"#34c975",NA,NA,NA,NA,NA,NA),
          border = NA,
          xlim = c(min_x_,max_x_),
          las = 1,
          axisnames = T,
          xaxt = "n",
          yaxt = "n",
          add = T)
  
  # TEanno
  barplot(replicate(length(names_arg_),x_all_),
          names.arg = names_arg_,
          horiz = T,
          beside = F,
          main = "",
          width = 0.65,
          space = 1,
          col = c(replicate(length(names_arg_)-8,NA),col_v_[2],NA,NA,NA,NA,NA,NA,NA),
          border = NA,
          xlim = c(min_x_,max_x_),
          las = 1,
          axisnames = T,
          xaxt = "n",
          yaxt = "n",
          add = T)
  
  
  for(j_ in 1:length(TEanno_cls_)){
    TEanno_cls_j_ <- TEanno_cls_[j_]
    x_TEanno_cls_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$TEclass_ %in% TEanno_cls_j_ & centromere_dissection_data_gw_$step_ == 4)]
    if(length(x_TEanno_cls_) > 0){
      fct_i_ <- 9+j_
      fct_j_ <- 8+j_
      TEanno_cls_col_v_ <- c(replicate(length(names_arg_)-fct_i_,NA),col_v_[2],replicate(fct_j_,NA)) %>% unlist()
      barplot(replicate(length(names_arg_),x_TEanno_cls_),
              names.arg = names_arg_,
              horiz = T,
              beside = F,
              main = "",
              width = 0.65,
              space = 1,
              col = TEanno_cls_col_v_,
              border = NA,
              xlim = c(min_x_,max_x_),
              las = 1,
              axisnames = T,
              xaxt = "n",
              yaxt = "n",
              add = T)
    }
  }
  
  # axis title
  # title(xlab = "Centromere Length (bp)")
  mtext("Centromere Length (bp)", side = 1, line = 3.5, cex = 0.75)
  
  # label
  text(grconvertX(0.000, from = "nfc", to = "user"),
       grconvertY(0.985, from = "nfc", to = "user"), 
       "f", xpd = NA, cex = 1.4, font = 2)
  
  
}

# empty panel
par(mar = c(0,0,0,0))
plot.new()

# centromere dissection genome wide
{
  par(mar = c(2, 4, 3, 4)-0.25)
  line_val_ <- 1
  
  # cent size
  x_cent_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "centromere_nt_")]
  x_cent_ <- ifelse(length(x_cent_) == 0, 0, x_cent_)
  # cent accounted by gaps
  x_gap_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "gap_nt_")]
  x_gap_ <- ifelse(length(x_gap_) == 0, 0, x_gap_)
  # cent accounted by TR
  x_TR_ <- x_cent_-x_gap_
  x_TR_ <- ifelse(length(x_TR_) == 0, 0, x_TR_)
  # cent accounted by TEanno
  x_all_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$TEclass_ == "all")]
  x_all_ <- ifelse(length(x_all_) == 0, 0, x_all_)
  # we add up centromere TR space to plot the value stacked-wise
  x_TEanno_ <- ifelse(x_all_ == 0, 0, x_all_+x_TR_)
  
  # cent accounted by TEanno
  x_otherTR_ <- centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_gw_$feature_ == "otherTR_nt_")]
  x_otherTR_ <- ifelse(length(x_otherTR_) == 0, 0, x_TEanno_+x_otherTR_)
  
  max_x_ <- max(centromere_dissection_data_gw_$total_length_[which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")],na.rm = T)
  # x left limit is defined so as to allocate the karypotype
  min_x_ <- (max_x_*.25)*(-1)
  
  # cent size
  bp_ <- barplot(x_cent_,
                 names.arg = centromere_dissection_data_gw_$filter_comb_[which(centromere_dissection_data_gw_$feature_ == "centromere_nt_")],
                 horiz = T,
                 beside = F,
                 xaxt = "n",
                 yaxt = "n",
                 main = paste0("all ",root_name_),
                 las = 1,
                 width = 0.65,
                 space = 1,
                 col = col_v_[3],
                 border = NA,
                 xlim = c(min_x_,max_x_),
                 cex.main = cex_main_,
                 adj = 0)
  
  # plot cent accounted by other TR
  barplot(x_otherTR_,
          names.arg = centromere_dissection_data_gw_$filter_comb_[which(centromere_dissection_data_gw_$TEclass_ == "all")],
          horiz = T,
          add = T,
          beside = F,
          xaxt = "n",
          yaxt = "n",
          main = paste0("all ",root_name_),
          las = 1,
          width = 0.65,
          space = 1,
          col = "#34c975",
          border = NA,
          cex.main = cex_main_,
          adj = 0)
  
  # plot cent accounted by TEanno
  barplot(x_TEanno_,
          names.arg = centromere_dissection_data_gw_$filter_comb_[which(centromere_dissection_data_gw_$TEclass_ == "all")],
          horiz = T,
          add = T,
          beside = F,
          xaxt = "n",
          yaxt = "n",
          main = paste0("all ",root_name_),
          las = 1,
          width = 0.65,
          space = 1,
          col = col_v_[2],
          border = NA,
          cex.main = cex_main_,
          adj = 0)
  
  # plot cent accounted by TR
  barplot(x_TR_,
          names.arg = centromere_dissection_data_gw_$filter_comb_[which(centromere_dissection_data_gw_$feature_ == "gap_nt_")],
          horiz = T,
          add = T,
          beside = F,
          xaxt = "n",
          yaxt = "n",
          main = paste0("all ",root_name_),
          las = 1,
          width = 0.65,
          space = 1,
          col = col_v_[1],
          border = NA,
          cex.main = cex_main_,
          adj = 0)
  
  text(
    y = bp_,
    x = 0-(max_x_*0.0275),
    adj = 1,
    labels = paste0(round((centromere_dissection_data_gw_$genome_prop_[which(centromere_dissection_data_gw_$feature_ == "centromere_nt_")])*100,digits = 1),"%")
  )
  
  axis(side = 1,
       at = pretty(c(0,max_x_)),
       line = line_val_)
  
  # label
  text(grconvertX(0.005, from = "nfc", to = "user"),
       grconvertY(0.985, from = "nfc", to = "user"), 
       "g", xpd = NA, cex = 1.4, font = 2)
  
}

# centromere dissection by seqID 
{
  par(mar = c(2, 4, 3, 4)-0.25)
  line_val_ <- 1
  
  # seqID_v_ <- unique(karyo_$seqID_)
  max_x_ <- max(centromere_dissection_data_seqID_$total_length_[which(centromere_dissection_data_seqID_$feature_ == "centromere_nt_")],na.rm = T)
  # x left limit is defined so as to allocate the karypotype
  min_x_ <- (max_x_*.25)*(-1)
  
  for(i_ in 1:length(seqID_v_)){
    seqID_ <- seqID_v_[i_]
    centromere_dissection_data_seqID_i_ <- centromere_dissection_data_seqID_[which(centromere_dissection_data_seqID_$seqID_ %in% seqID_),]
    
    # cent size
    x_cent_ <- centromere_dissection_data_seqID_i_$total_length_[which(centromere_dissection_data_seqID_i_$feature_ == "centromere_nt_")]
    x_cent_ <- ifelse(length(x_cent_) == 0, 0, x_cent_)
    # cent accounted by gaps
    x_gap_ <- centromere_dissection_data_seqID_i_$total_length_[which(centromere_dissection_data_seqID_i_$feature_ == "gap_nt_")]
    x_gap_ <- ifelse(length(x_gap_) == 0, 0, x_gap_)
    # cent accounted by TR
    x_TR_ <- x_cent_-x_gap_
    x_TR_ <- ifelse(length(x_TR_) == 0, 0, x_TR_)
    # cent accounted by TEanno
    x_all_ <- centromere_dissection_data_seqID_i_$total_length_[which(centromere_dissection_data_seqID_i_$TEclass_ == "all")]
    x_all_ <- ifelse(length(x_all_) == 0, 0, x_all_)
    # we add up centromere TR space to plot the value stacked-wise
    x_TEanno_ <- ifelse(x_all_ == 0, 0, x_all_+x_TR_)
    # cent accounted by TEanno
    x_otherTR_ <- centromere_dissection_data_seqID_i_$total_length_[which(centromere_dissection_data_seqID_i_$feature_ == "otherTR_nt_")]
    x_otherTR_ <- ifelse(length(x_otherTR_) == 0, 0, x_TEanno_+x_otherTR_)
    
    # plot cent length
    if(x_cent_ > 0){
      bp_ <- barplot(x_cent_,
                     names.arg = seqID_,
                     horiz = T,
                     beside = F,
                     xaxt = "n",
                     yaxt = "n",
                     main = paste0("\n",seqID_),
                     las = 1,
                     width = 0.65,
                     space = 1,
                     col = col_v_[3],
                     border = NA,
                     xlim = c(min_x_,max_x_),
                     cex.main = cex_main_,
                     adj = 0)
    }
    
    # plot cent accounted by otherTR
    if(x_otherTR_ > 0){
      bp_ <- barplot(x_otherTR_,
                     names.arg = seqID_,
                     horiz = T,
                     add = T,
                     beside = F,
                     xaxt = "n",
                     yaxt = "n",
                     main = paste0("\n",seqID_),
                     las = 1,
                     width = 0.65,
                     space = 1,
                     col = "#34c975",
                     border = NA,
                     xlim = c(min_x_,max_x_),
                     cex.main = cex_main_,
                     adj = 0)
    }
    
    # plot cent accounted by TEanno
    if(x_TEanno_ > 0){
      bp_ <- barplot(x_TEanno_,
                     names.arg = seqID_,
                     horiz = T,
                     add = T,
                     beside = F,
                     xaxt = "n",
                     yaxt = "n",
                     main = paste0("\n",seqID_),
                     las = 1,
                     width = 0.65,
                     space = 1,
                     col = col_v_[2],
                     border = NA,
                     xlim = c(min_x_,max_x_),
                     cex.main = cex_main_,
                     adj = 0)
    }
    
    # plot cent accounted by TR
    if(x_TR_ > 0){
      bp_ <- barplot(x_TR_,
                     names.arg = seqID_,
                     horiz = T,
                     add = T,
                     beside = F,
                     xaxt = "n",
                     yaxt = "n",
                     main = paste0("\n",seqID_),
                     las = 1,
                     width = 0.65,
                     space = 1,
                     col = col_v_[1],
                     border = NA,
                     xlim = c(min_x_,max_x_),
                     cex.main = cex_main_,
                     adj = 0)
    }
    
    axis(side = 1,
         at = pretty(c(0,max_x_)),
         line = line_val_)
    
    # calculate relative positions for start/end chromosome and arrays
    karyo_i_ <- karyo_[which(karyo_$ENA_ %in% seqID_),]
    max_karyo_length_ <- max(karyo_$seq_length_,na.rm = T)
    
    min_karyo_x_ <- (max_x_*.225)*(-1)
    max_karyo_x_ <- (max_x_*.025)*(-1)
    
    rel_length_ <- max_karyo_x_-min_karyo_x_
    karyo_end_ <- ((karyo_i_$seq_length_*rel_length_)/max_karyo_length_)
    karyo_end_ <- karyo_end_+min_karyo_x_
    
    rect(min_karyo_x_,
         bp_-.125,
         karyo_end_,
         bp_+.125,
         col = "grey75",
         border = NA)
    
    if(length(which(arrays_data_$seqID %in% seqID_)) > 0){
      arrays_data_i_ <- arrays_data_[which(arrays_data_$seqID %in% seqID_),]
      
      array_start_ <- ((arrays_data_i_$start*rel_length_)/max_karyo_length_)
      array_start_ <- array_start_+min_karyo_x_
      array_end_ <- ((arrays_data_i_$end*rel_length_)/max_karyo_length_)
      array_end_ <- array_end_+min_karyo_x_
      
      rect(array_start_,
           bp_-.125,
           array_end_,
           bp_+.125,
           col = "grey35",
           border = NA)
    }
    
    text(
      y = bp_,
      x = 0-(max_x_*0.0275),
      adj = 1,
      labels = paste0(round((centromere_dissection_data_seqID_i_$genome_prop_[which(centromere_dissection_data_seqID_i_$feature_ == "centromere_nt_")])*100,digits = 1),"%")
    )
    
  }
}

# empty panel
par(mar = c(0,0,0,0))
plot.new()

dev.off()
