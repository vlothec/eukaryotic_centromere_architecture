library(tidyverse)
rm(list = ls())

pattern_ <- "_TEanno_isec_gaps_parsed_"
file_list_ <- list.files(pattern = pattern_)
species_list_ <- file_list_ %>%
  {gsub(pattern_,"",.)} %>%
  unique

target_species_ <- "drHedHeli1|ddEupPepu3|daMisOron1|daLinVulg1|icHarAxyr1|ibEctPall1|icCocSept1"
species_list_ <- species_list_[grepl(target_species_,species_list_)]

TE_class_list_updated_ <- read.table("TE_class_list_updated_",sep = "\t",header = F)
colnames(TE_class_list_updated_) <- c("TEanno_cls_","old_TEanno_cls_","TEanno_feature_")

col_pal_ <- c("Class_I_LTR" = "#03a9fc", "Class_I_non_LTR_other" = "#fc6b03", "Class_I_non_LTR_TPRT" = "#fc3f65", "Class_II_DNA_element" = "#fcc11e", "TE_unclass" = "black", "repeat_region" = "purple", "rRNA_gene" = "cyan")


# c("LTR_retrotransposon" = "darkblue", "Class_I_Retrotransposon" = "orange", "nonLTR_retrotransposon" = "brown", "Class_II_DNA_Transposon" = "#F8D36B", "repeat_region" = "red", "rRNA_gene" = "cyan", "TE_unclass" = "black")

for(i_ in 1:length(species_list_)){
  species_ <- species_list_[i_]
  print(paste0(i_,": ",species_))
  file_i_ <- file_list_[grepl(species_,file_list_)]
  if(file.info(file_i_)$size > 0){
    data_ <- read.table(file_i_, header = T, sep = "\t")
    
    gaps_n_ <- length(unique(data_$gap_length_id_))
    data_ <- left_join(data_,TE_class_list_updated_[,c(1,3)],by = "TEanno_feature_")
    
    data_ <- data_ %>%
      mutate(gap_rel_start_ = gap_start_ - gap_start_,
             gap_rel_end_ = gap_end_ - gap_start_,
             TEanno_rel_start_ = TEanno_start_ - gap_start_,
             TEanno_rel_end_ = TEanno_end_ - gap_start_) %>%
      mutate(TEanno_rel_start_ = ifelse(TEanno_rel_start_ < 0, 0, TEanno_rel_start_),
             TEanno_rel_end_ = ifelse(TEanno_rel_end_ < 0, 0, TEanno_rel_end_),
             TEanno_rel_end_ = ifelse(TEanno_rel_end_ > gap_length_, gap_length_, TEanno_rel_end_),
             label_ = paste0(gap_length_id_," - ",sprintf("%02d",as.numeric(as.factor(gap_seqID_))))) %>%
      group_by(gap_length_id_)
    
    
    data_ <- data_[with(data_,order(gap_length_)),]
    data_$gap_length_id_ <- data_$gap_length_id_ %>%
      {factor(.,levels = unique(data_$gap_length_id_))}
    data_$label_ <- data_$label_ %>%
      {factor(.,levels = unique(data_$label_))}
    
    
    karyotype_file_ <- list.files(pattern = "_karyotype_ext_$")
    karyotype_file_ <- karyotype_file_[grepl(gsub("\\..*","",species_),karyotype_file_)]
    if(length(karyotype_file_) > 0){
      karyotype_data_ <- read.table(paste0(karyotype_file_),sep = "\t",header = F)
      colnames(karyotype_data_) <- c("ENA_","ChrOnly","length_")
    }
    
    path_to_intact_ <- ""
    intact_file_ <- list.files(path = path_to_intact_, pattern = ".intact.gff3$")
    intact_file_ <- intact_file_[grepl(species_,intact_file_)]
    
    intact_data_ <- read.table(paste0(path_to_intact_,intact_file_), sep = "\t", header = T)
    colnames(intact_data_) <- c("seqID","method_anno_","TEanno_feature_","start_","end_","score_","strand_","phase_","attributes_")
    
    intact_data_ <- intact_data_[!grepl("long_terminal_repeat|target_site_duplication|repeat_region",intact_data_$TEanno_feature_),]
    # intact_data_ <- intact_data_[grepl("_LTR_retrotransposon|^LTR_retrotransposon",intact_data_$TEanno_feature_),]
    
    # minimum intact length == 2000 bp
    intact_data_ <- intact_data_ %>%
      mutate(length_ = abs(end_-start_)+1)
    intact_data_ <- intact_data_[which(intact_data_$length_ >= 2000),]
    
    seqID_type_ <- which(intact_data_$seqID %in% karyotype_data_$ENA_) %>% length()
    if(seqID_type_ == 0){
      colnames(karyotype_data_)[2] <- "seqID"
      intact_data_ <- left_join(intact_data_,karyotype_data_[,c(1,2)])
      intact_data_$seqID <- intact_data_$ENA_
      intact_data_ <- intact_data_[,-ncol(intact_data_)]
    }
    
    data_ <- data_[which(data_$TEanno_start_ > 0),]
    gaps_n_ <- length(unique(data_$gap_length_id_))
    
    data_$intact_ <- FALSE
    idx_ <- which(data_$TEanno_start_ > 0)
    # for(h_ in 1:length(idx_)){
    for(h_ in 1:nrow(data_)){
        
      if(data_$TEanno_start_[h_] > 0){
        idx_h_ <- h_
        
        TEanno_seqID_ <- data_$gap_seqID_[idx_h_]
        TEanno_start_ <- data_$gap_start_[idx_h_]+1
        TEanno_end_ <- data_$gap_end_[idx_h_]
        gap_length_ <- data_$gap_length_[idx_h_]
        
        jdx_ <- which(intact_data_$seqID == TEanno_seqID_ & intact_data_$start_ >= TEanno_start_-10 & intact_data_$end_ <= TEanno_end_+10 )
        
        if(length(jdx_) > 0){
          data_$intact_[idx_h_] <- TRUE
        }
      }
    }
    
    table(data_$intact_)
    length(unique(data_$gap_length_id_[which(data_$intact_ == T)]))
    
    solo_file_ <- list.files(pattern = "_wga_soloLTR.gff3$")
    solo_file_ <- solo_file_[grepl(species_,solo_file_)]
    
    if(nrow(file.info(solo_file_)) > 0){
      solo_data_ <- read.table(paste0(solo_file_), sep = "\t", header = T)
      colnames(solo_data_) <- c("seqID","method_anno_","solo_feature_","start_","end_","score_","strand_","phase_","attributes_")
      
      seqID_type_ <- which(solo_data_$seqID %in% karyotype_data_$ENA_) %>% length()
      if(seqID_type_ == 0){
        colnames(karyotype_data_)[2] <- "seqID"
        solo_data_ <- left_join(solo_data_,karyotype_data_[,c(1,2)])
        solo_data_$seqID <- solo_data_$ENA_
        solo_data_ <- solo_data_[,-ncol(solo_data_)]
      }
      
      solo_data_ <- solo_data_ %>%
        mutate(length_ = abs(end_-start_)+1)
      
      binwidth_ <- 1
      max_val_ <- max(solo_data_$length_,na.rm = T)
      myBreaks_ <- seq(0,max_val_,length.out = round(max_val_/binwidth_))
      hist(solo_data_$length_,
           breaks = myBreaks_,
           col = "blue",
           border = NA,
           xlim = c(0,1000))
      
      #  minimum soloLTR length == 150 bp
      solo_data_ <- solo_data_[which(solo_data_$length_ >= 150),]
      
      data_$solo_ <- FALSE
      # idx_ <- which(data_$TEanno_start_ > 0)
      for(h_ in 1:nrow(data_)){
        
        idx_h_ <- h_
        
        if(data_$TEanno_start_[h_] > 0){
          TEanno_seqID_ <- data_$gap_seqID_[idx_h_]
          TEanno_start_ <- data_$gap_start_[idx_h_]+1
          TEanno_end_ <- data_$gap_end_[idx_h_]
          gap_length_ <- data_$gap_length_[idx_h_]
          
          jdx_ <- which(solo_data_$seqID == TEanno_seqID_ & solo_data_$start_ >= TEanno_start_-10 & solo_data_$end_ <= TEanno_end_+10 & gap_length_-solo_data_$length_ <= 100 & gap_length_-solo_data_$length_ >= -100)
          
          if(length(jdx_) > 0){
            data_$solo_[idx_h_] <- TRUE
          }
        }
      }
      
      table(data_$solo_)
      length(unique(data_$gap_length_id_[which(data_$solo_ == T)]))
    }
    
    # axis limits
    if (length(unique(data_$gap_length_id_)) > 636) {
      ylim <- c(0, length(unique(data_$gap_length_id_)))
    } else {
      ylim <- c(-(636-length(unique(data_$gap_length_id_))), length(unique(data_$gap_length_id_)))
    }
    data_$gap_length_id_num_ <- data_$gap_length_id_ %>%
      as.numeric()
    
    fct_height_ <- (gaps_n_*3)/300
    pdf(paste0(species_,"_mosaic_gaps_for_manuscript_.pdf"),onefile = T, width = 2.15,height = fct_height_)
    plot_ <- ggplot(data_) +
      geom_segment(aes(x = 0, xend = gap_rel_end_, y = gap_length_id_, yend = gap_length_id_), color = "grey90", linewidth = .35) +
      geom_segment(data = data_[which(data_$TEanno_start_ > 0),], aes(x = TEanno_rel_start_, xend = TEanno_rel_end_, y = gap_length_id_, yend = gap_length_id_, color = TEanno_cls_), linewidth = .35) +
      geom_text(data = data_[which(data_$intact_ == T),],aes(x = -1000, y = gap_length_id_,label = "-"), color = "black", size = 1, vjust = 0.34) +
      geom_text(data = data_[which(data_$solo_ == T),],aes(x = -1000, y = gap_length_id_,label = "-"), color = "grey80", size = 1, vjust = 0.34) +
      scale_color_manual(values = col_pal_) +
      labs(title = paste0(species_,"\nn: ",gaps_n_)) +
      # coord_cartesian(ylim = ylim) +
      coord_cartesian(xlim = c(-1000,25000)) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none",
            panel.grid = element_blank(),
            axis.line.x = element_line(linewidth = .65),
            axis.ticks.x = element_line(linewidth = .65),
            text = element_text(size = 3.5))
    print(plot_)
    dev.off()
  }
}


