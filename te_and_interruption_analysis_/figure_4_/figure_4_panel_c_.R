# G. rivale
{
  library(tidyverse)
  options(scipen=999)
  
  # proportions
  
  intact.prop <- 0.156741
  frag.prop <- 0.189503
  prot.prop <-  0.0007585 # PROT found in GAPS
  non_athila_LTRRT.prop <- 0.0112407+0.00300261
  unclassified_LTRRT.prop <- 0.000112797+0.00339904
  nonLTR.prop <- 0.000394139+0.00991287
  # drGeuRiva1.hap1.1_centromere_dissection_data_gw_
  satellite.prop <- 0.657470731502205
  unnanotated.prop <- 0.0112906551548804
  
  overlap_degree_ <- abs(1-(satellite.prop+(intact.prop+frag.prop+prot.prop+non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop)+unnanotated.prop))
  
  recalc_prop_ <- rbind(intact.prop,frag.prop,prot.prop,non_athila_LTRRT.prop,unclassified_LTRRT.prop,nonLTR.prop) %>%
   {./0.3743062} %>%
   {.*overlap_degree_}

  intact.prop <- (intact.prop - recalc_prop_[1,] )
  frag.prop <- (frag.prop - recalc_prop_[2,] )
  non_athila_LTRRT.prop <- ( non_athila_LTRRT.prop - recalc_prop_[4,] ) 
  unclassified_LTRRT.prop <- ( unclassified_LTRRT.prop - recalc_prop_[5,] ) 
  nonLTR.prop <- ( nonLTR.prop - recalc_prop_[6,] ) 
  
  definetly_unnanotated_ <- 1-(satellite.prop+(intact.prop+frag.prop+prot.prop+non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop)+unnanotated.prop)
  
  cent_feature_prop <- data.frame("chr" = "all_CENT",
                                  "ft.comb.red.ii" = c("other TEs","GETHILA","sat"),
                                  "total.cent.feature.space" = c((non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop),(intact.prop+frag.prop+prot.prop),satellite.prop))
  
  text.size <- 12
  col_v_ <- c("TRASH" = "#D83067","all_TEs" = "#487CE3", "unknown" = "#DDDDDD")
  
  pdf("drGeuRiva1.hap1.1_cent_feature_prop_barplot_all_anno.pdf", width = 9, height = 1.15, onefile = T)
  ggplot(cent_feature_prop) +
  
    # bars
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = 1),aes(y = chr, x = total.cent.feature.space), fill = "grey92", color = "black", width = 0.2) +
    geom_col(aes(y = chr, x = total.cent.feature.space, fill = fct_inorder(ft.comb.red.ii)), color = "black", width = 0.2,show.legend=TRUE) +
    
    # sat
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.prop)),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +
    
    # ATHILA intact
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.prop+intact.prop)),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +
    
    # ATHILA intact
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.prop+intact.prop+frag.prop)),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +
    
    # ATHILA frag
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.prop+intact.prop+frag.prop)),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +
    
    # ATHILA frag
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.prop+intact.prop+frag.prop+non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop+prot.prop)),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +
    scale_x_continuous(expand = expansion(0), labels = c(0,25,50,75,100)) +
    scale_y_discrete(expand = expansion(add = c(.3, .3))) +
    scale_fill_manual(values = rev(c("GETHILA" = "#487CE3",
                                     "other TEs" = "grey50",
                                     "unknown" = "grey92",
                                     sat = "#D83067")),drop = FALSE) +
    coord_cartesian(expand = T,clip = "off") +
    theme_minimal() + 
    guides(
      fill = guide_legend(override.aes = list(linewidth = .5,byrow = TRUE))
    ) +
    theme(axis.title = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_blank(),
          legend.position = "right",
          legend.key = element_rect(),
          legend.key.size = unit(.85,"line"),
          legend.text = element_text(),
          legend.title = element_blank(),
          legend.spacing.y = unit(5,"line"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97",colour = NA),
          text = element_text(size = text.size)) 
  
  dev.off()
  
  # now the same but using real space
  centromere_space_ <- 56756671
  satellite.space_ <- satellite.prop*centromere_space_
  unnanotated.space_ <- unnanotated.prop*centromere_space_
  
  intact.space_ <- intact.prop*centromere_space_
  frag.space_ <- frag.prop*centromere_space_
  non_athila_LTRRT.space_ <- non_athila_LTRRT.prop*centromere_space_
  unclassified_LTRRT.space_ <- unclassified_LTRRT.prop*centromere_space_
  nonLTR.space_ <- nonLTR.prop*centromere_space_
  prot.space_ <- 0
  
  text.size <- 12
  pdf("drGeuRiva1.hap1.1_cent_feature_prop_barplot_all_anno_total_space_.pdf", width = 9, height = 1.15, onefile = T)
  ggplot() +
  
    # bars
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = centromere_space_),aes(y = chr, x = total.cent.feature.space), fill = "grey92", color = "black", width = 0.2) +
    
    # other TEs
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.space_+intact.space_+frag.space_+non_athila_LTRRT.space_+unclassified_LTRRT.space_+nonLTR.space_+prot.space_),feature_ = "other_TEs"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    # sat
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.space_+intact.space_+frag.space_),feature_ = "sat"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    # ATHILA intact+frag
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (intact.space_+frag.space_),feature_ = "GETHILA"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    # ATHILA intact
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (intact.space_),feature_ = "GETHILA"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    scale_x_continuous(expand = expansion(0), labels = ~ paste0(./1000000,"Mb")) +
    scale_y_discrete(expand = expansion(add = c(.3, .3))) +
    scale_fill_manual(values = rev(c("GETHILA" = "#487CE3",
                                     "other_TEs" = "grey50",
                                     "unknown" = "grey92",
                                     sat = "#D83067")),drop = FALSE) +
    coord_cartesian(expand = T,clip = "off",xlim = c(0,60000000)) +
    theme_minimal() + 
    guides(
      fill = guide_legend(override.aes = list(linewidth = .5,byrow = TRUE))
    ) +
    theme(axis.title = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.y = element_blank(),
          legend.position = "right",
          legend.key = element_rect(),
          legend.key.size = unit(.85,"line"),
          legend.text = element_text(),
          legend.title = element_blank(),
          legend.spacing.y = unit(5,"line"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97",colour = NA),
          text = element_text(size = text.size)) 
  
  dev.off()

}

# G. urbanum
{
  # the following snippet was run in the server to get the satellite total space
  
  rm(list = ls())
  # proportions
  
  intact.prop <- 0.290919
  frag.prop <- 0.584015
  prot.prop <- 0.02194767 # PROT found in GAPS
  non_athila_LTRRT.prop <- 0.0358482
  unclassified_LTRRT.prop <- 0.0163805
  nonLTR.prop <- 0.0123491
  satellite.prop <- 0.0087367
  unnanotated.prop <- 0
  
  cent_feature_prop <- data.frame("chr" = "all_CENT",
                                  "ft.comb.red.ii" = c("other TEs","GETHILA"),
                                  "total.cent.feature.space" = c(0.0646,(intact.prop+frag.prop+prot.prop)))
  
  text.size <- 12
  
  pdf("drGeuUrba1.1_cent_feature_prop_barplot_all_anno.pdf", width = 9, height = 1.15, onefile = T)
  ggplot(cent_feature_prop) +
    # grid
    # ATHILA intact
    geom_vline(xintercept = intact.prop, linetype = "solid",color = "white",linewidth = 1.25) +
    
    # all ATHILA
    geom_vline(xintercept = (intact.prop+frag.prop+prot.prop), linetype = "solid",color = "white",linewidth = 1.25) +
    
    # any TE
    geom_vline(xintercept = (intact.prop+frag.prop+non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop+prot.prop), linetype = "solid",color = "white",linewidth = 1.25) +
    
    # bars
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = 1),aes(y = chr, x = total.cent.feature.space), fill = "grey92", color = "black", width = 0.2) +
    geom_col(aes(y = chr, x = total.cent.feature.space, fill = fct_inorder(ft.comb.red.ii)), color = "black", width = 0.2,show.legend=TRUE) +
    
    # ATHILA intact
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = intact.prop),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +
    
    scale_x_continuous(expand = expansion(0), labels = c(0,25,50,75,100)) +
    scale_y_discrete(expand = expansion(add = c(.3, .3))) +
    scale_fill_manual(values = rev(c("GETHILA" = "#487CE3",
                                     "other TEs" = "grey50",
                                     "unknown" = "grey92",
                                     sat = "#D83067")),drop = FALSE) +
    coord_cartesian(expand = T,clip = "off") +
    theme_minimal() + 
    guides(
      fill = guide_legend(override.aes = list(linewidth = .5,byrow = TRUE))
    ) +
    theme(axis.title = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_blank(),
          legend.position = "right",
          legend.key = element_rect(),
          legend.key.size = unit(.85,"line"),
          legend.text = element_text(),
          legend.title = element_blank(),
          legend.spacing.y = unit(5,"line"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97",colour = NA),
          text = element_text(size = text.size)) 
  
  dev.off()
  
  
  # now the same but using real space
  centromere_space_ <- 18000000
  satellite.space_ <- 157260
  unnanotated.space_ <- unnanotated.prop*centromere_space_
  
  intact.space_ <- intact.prop*centromere_space_
  frag.space_ <- frag.prop*centromere_space_
  non_athila_LTRRT.space_ <- non_athila_LTRRT.prop*centromere_space_
  unclassified_LTRRT.space_ <- unclassified_LTRRT.prop*centromere_space_
  nonLTR.space_ <- nonLTR.prop*centromere_space_
  prot.space_ <- prot.prop*centromere_space_
  
  text.size <- 12
  pdf("drGeuUrba1.1_cent_feature_prop_barplot_all_anno_total_space_.pdf", width = 9, height = 1.15, onefile = T)
  ggplot() +
    
    # bars
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = centromere_space_),aes(y = chr, x = total.cent.feature.space), fill = "grey92", color = "black", width = 0.2) +
    
    # other TEs
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.space_+intact.space_+frag.space_+non_athila_LTRRT.space_+unclassified_LTRRT.space_+nonLTR.space_+prot.space_),feature_ = "other_TEs"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    # sat
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (satellite.space_+intact.space_+frag.space_),feature_ = "sat"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    # ATHILA intact+frag
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (intact.space_+frag.space_),feature_ = "GETHILA"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    # ATHILA intact
    geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = (intact.space_),feature_ = "GETHILA"),aes(y = chr, x = total.cent.feature.space, fill = feature_), color = "black", width = 0.2) +
    
    scale_x_continuous(expand = expansion(0), labels = ~ paste0(./1000000,"Mb")) +
    scale_y_discrete(expand = expansion(add = c(.3, .3))) +
    scale_fill_manual(values = rev(c("GETHILA" = "#487CE3",
                                     "other_TEs" = "grey50",
                                     "unknown" = "grey92",
                                     "sat" = "#D83067")),drop = FALSE) +
    coord_cartesian(expand = T,clip = "off",xlim = c(0,60000000)) +
    theme_minimal() + 
    guides(
      fill = guide_legend(override.aes = list(linewidth = .5,byrow = TRUE))
    ) +
    theme(axis.title = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.y = element_blank(),
          legend.position = "right",
          legend.key = element_rect(),
          legend.key.size = unit(.85,"line"),
          legend.text = element_text(),
          legend.title = element_blank(),
          legend.spacing.y = unit(5,"line"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97",colour = NA),
          text = element_text(size = text.size)) 
  
  dev.off()
  
}
