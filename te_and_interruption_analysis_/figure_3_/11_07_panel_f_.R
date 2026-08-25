# panel E2

TE_class_list_updated_ <- read.table("TE_class_list_updated_",sep = "\t",header = F)
colnames(TE_class_list_updated_) <- c("TEanno_cls_","old_TEanno_cls_","TEclass_")

centromere_dissection_data_gw_all_ <- read.table("centromere_dissection_data_gw_mono_aggregated_", sep = "\t", header = T)
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$feature_ == "TE_nt_" & centromere_dissection_data_gw_all_$step_ == 3),]

dtolID_order_ <- readLines("dtolID_order_")

centromere_dissection_data_gw_all_$phyloclade_ <- centromere_dissection_data_gw_all_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}

centromere_dissection_data_gw_all_$dtolID_ <- centromere_dissection_data_gw_all_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(dtolID_,phyloclade_)),]

# dummy_
{
  dummy_ <- data.frame(distinct(centromere_dissection_data_gw_all_[,c(11,12,13)]), feature_ = "TE_nt_")
  dummy_ <- left_join(data.frame("dtolID_" = dtolID_order_),dummy_)
  dummy_$phylotaxa_[which(dummy_$dtolID_ == "idBibMarc1.1")] <- "Diptera"
  dummy_$phyloclade_[which(dummy_$dtolID_ == "idBibMarc1.1")] <- "invertebrate"
  dummy_$phylotaxa_[which(dummy_$dtolID_ == "fToxJac2.1")] <- "Actinopterygii"
  dummy_$phyloclade_[which(dummy_$dtolID_ == "fToxJac2.1")] <- "chordate"
  dummy_$feature_[which(is.na(dummy_$feature_) == T)] <- "TE_nt_"

  dummy_$dtolID_ <- dummy_$dtolID_ %>%
    {factor(.,levels = dtolID_order_)}
}

text.size <- 6.5
bar.width <- 0.425
linewidth_ <- 0


centromere_dissection_data_gw_all_$total_length_[which(centromere_dissection_data_gw_all_$dtolID_ %in% c("laApoDist1.hap1.1","lsIriFoet1.hap1.1"))] <- 20000000
pdf(paste0("draft_plots_/figure_3_panel_f_.pdf"),onefile = T, width = 1.5,height = 12.275)
ggplot(dummy_) +
  geom_point(data = dummy_, aes(y = dtolID_, x = feature_), color = "transparent",stroke = linewidth_) +
  # geom_col(aes(x = max(centromere_dissection_data_gw_all_$total_length_)/1000000, y = dtolID_), fill = "grey95", color = NA, width = .5) +
  # geom_col(aes(x = total_length_/1000000, y = dtolID_, fill = feature_), color = NA, width = .5) +
  geom_point(data = centromere_dissection_data_gw_all_, aes(y = dtolID_, x = feature_, color = feature_, size = total_length_),stroke = linewidth_,alpha=.5) +
  scale_color_manual(values = c("TE_nt_" = "black")) + 
  scale_size_continuous(range = c(.5,6.05),guide = NULL) +
  facet_grid(phyloclade_ ~ feature_, scales = "free", space = "free") +
  scale_y_discrete(expand = expansion(add = .85)) +
  # guides(color = "none") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        axis.title = element_blank(),
        axis.line = element_line(colour = NA, linewidth = linewidth_),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        panel.background = element_rect(fill = NA, color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey98"),
        panel.ontop = F,
        # panel.spacing.x = unit(8,"pt"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        plot.title = element_text(hjust = 0),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(size = text.size))
dev.off()

