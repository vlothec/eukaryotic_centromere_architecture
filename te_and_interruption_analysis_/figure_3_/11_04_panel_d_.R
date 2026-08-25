# panel D
rm(list = ls())
LTR_piden_data_gw_aggregated_ <- data.table::fread("LTR_piden_data_gw_aggregated_", sep = "\t", header = T) %>% data.frame()

dtolID_order_ <- readLines("dtolID_order_")
setdiff(unique(LTR_piden_data_gw_aggregated_$dtolID_), dtolID_order_)
setdiff(dtolID_order_,unique(LTR_piden_data_gw_aggregated_$dtolID_))

LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_[which(LTR_piden_data_gw_aggregated_$dtolID_ != "ihAelAcum1.1"),]

LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_[which(LTR_piden_data_gw_aggregated_$dtolID_ %in% dtolID_order_),]

LTR_piden_data_gw_aggregated_$dtolID_ <- LTR_piden_data_gw_aggregated_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_[with(LTR_piden_data_gw_aggregated_,order(dtolID_)),]

LTR_piden_data_gw_aggregated_$phyloclade_ <- LTR_piden_data_gw_aggregated_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}

# LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_[-which(is.na(LTR_piden_data_gw_aggregated_$ltr_identity) == T  & is.na(LTR_piden_data_gw_aggregated_$start_) == F),]
LTR_piden_data_gw_aggregated_$centromere_status_[which(is.na(LTR_piden_data_gw_aggregated_$centromere_status_) == TRUE)] <- "arm"

# dummy_
{
  dummy_ <- data.frame(distinct(LTR_piden_data_gw_aggregated_[,c(13,15,14)]), ltr_identity = 0.99)
  dummy_ <- left_join(data.frame("dtolID_" = dtolID_order_),dummy_)
  dummy_$ltr_identity[which(is.na(dummy_$ltr_identity) == T)] <- 0.99
  dummy_$phyloclade_[which(is.na(dummy_$phyloclade_) == T)] <- "invertebrate"
  dummy_$phylotaxa_[which(is.na(dummy_$phylotaxa_) == T)] <- "Coleoptera"
  
  dummy_$dtolID_ <- dummy_$dtolID_ %>%
    {factor(.,levels = dtolID_order_)}
}

LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_[grepl("structural",LTR_piden_data_gw_aggregated_$attributes_),]

LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_ %>%
  group_by(dtolID_,centromere_status_) %>%
  slice_sample(n = 1000)

text.size <- 6.5
bar.width <- 0.425

pdf(paste0("draft_plots_/figure_3_panel_d_.pdf"),onefile = T, width = 1.8,height = 12.15)
png(paste0("draft_plots_/figure_3_panel_d_.png"),width = 1.8, height = 12.15, units = "in",res = 900)
ggplot(dummy_) +
  geom_jitter(aes(x = ltr_identity, y = dtolID_), color = "transparent", size = .05) +
  # geom_jitter(data = LTR_piden_data_gw_aggregated_[which(LTR_piden_data_gw_aggregated_$centromere_status_ == "arm"),], aes(x = ltr_identity, y = dtolID_, color = centromere_status_),width = 0, size = .05) +
  # geom_jitter(data = LTR_piden_data_gw_aggregated_, aes(x = ltr_identity, y = dtolID_, color = centromere_status_), alpha = .5, position = position_jitterdodge(dodge.width = .675), stroke = 0, size = .25) +
  geom_jitter(data = LTR_piden_data_gw_aggregated_[which(LTR_piden_data_gw_aggregated_$centromere_status_ == "arm"),], aes(x = ltr_identity, y = dtolID_, color = centromere_status_), alpha = .5, stroke = 0, size = .45) +
  geom_jitter(data = LTR_piden_data_gw_aggregated_[which(LTR_piden_data_gw_aggregated_$centromere_status_ != "arm"),], aes(x = ltr_identity, y = dtolID_, color = centromere_status_), alpha = .5, stroke = 0, size = .45) +
  scale_color_manual(values = c("arm" = "#B5D6E4", "cent" = "#000A83")) +
  scale_y_discrete(expand = expansion(add = .85)) +
  facet_grid(phyloclade_ ~ ., scales = "free_y", space = "free_y",drop = F) +
  guides() +
  coord_cartesian(expand = F,clip = "off", xlim = c(0.825,1)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey95"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = text.size))
dev.off()

