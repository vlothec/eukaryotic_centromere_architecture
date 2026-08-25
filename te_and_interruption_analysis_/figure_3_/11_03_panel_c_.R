# panel c

library(tidyverse)
rm(list = ls())

TEanno_data_all_method_data_gw_aggregated_ <- read.table("TEanno_data_all_method_data_gw_aggregated_", sep = "\t", header = T)

dtolID_order_ <- readLines("dtolID_order_")
setdiff(unique(TEanno_data_all_method_data_gw_aggregated_$dtolID_), dtolID_order_)
setdiff(dtolID_order_,unique(TEanno_data_all_method_data_gw_aggregated_$dtolID_))

TEanno_data_all_method_data_gw_aggregated_ <- TEanno_data_all_method_data_gw_aggregated_[which(TEanno_data_all_method_data_gw_aggregated_$dtolID_ != "ihAelAcum1.1"),]

TEanno_data_all_method_data_gw_aggregated_$dtolID_ <- TEanno_data_all_method_data_gw_aggregated_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
TEanno_data_all_method_data_gw_aggregated_ <- TEanno_data_all_method_data_gw_aggregated_[with(TEanno_data_all_method_data_gw_aggregated_,order(dtolID_)),]

TEanno_data_all_method_data_gw_aggregated_$phyloclade_ <- TEanno_data_all_method_data_gw_aggregated_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}
TEanno_data_all_method_data_gw_aggregated_ <- TEanno_data_all_method_data_gw_aggregated_[complete.cases(TEanno_data_all_method_data_gw_aggregated_),]

# dummy 
dummy_ <- data.frame(dtolID_ = c("fToxJac2.1","idBibMarc1.1"), phyloclade_ = c("chordate","invertebrate"), centromere_status_ = c("cent","cent"))
dummy_$phyloclade_ <- dummy_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}

text.size <- 6.5
bar.width <- 0.425

pdf(paste0("draft_plots_/figure_3_panel_c_.pdf"),onefile = T, width = 2,height = 12.15)
ggplot(TEanno_data_all_method_data_gw_aggregated_) +
  geom_col(aes(y = dtolID_, x = prop_centromere_status_length_,fill = Method), color = NA, width = .5) +
  geom_col(data = dummy_,aes(y = dtolID_, x = 1),fill = "grey90", color = NA, width = .5) +
  scale_fill_manual(values = c("homology" = "#B5D6E4", "structural" = "#000A83")) +
  scale_y_discrete(expand = expansion(add = .85)) +
  facet_grid(phyloclade_ ~ centromere_status_, scales = "free_y", space = "free_y") +
  guides(color = "none") +
  coord_cartesian(expand = F,clip = "on", xlim = c(0,1)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.background = element_rect(fill = "grey97", color = NA),
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





