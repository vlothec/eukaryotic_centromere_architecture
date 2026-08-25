# panel B
centromere_dissection_data_gw_all_ <- read.table("centromere_dissection_data_gw_mono_aggregated_", sep = "\t", header = T)
dtolID_v_ <- centromere_dissection_data_gw_all_$dtolID_ %>% unique()

centromere_dissection_data_gw_all_$phyloclade_ <- centromere_dissection_data_gw_all_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}

centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(centromere_prop_)),]
dtolID_order_ <- centromere_dissection_data_gw_all_$dtolID_[which(centromere_dissection_data_gw_all_$feature_ == "TR_nt_")] %>% unique()

centromere_dissection_data_gw_all_$dtolID_ <- centromere_dissection_data_gw_all_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(dtolID_,phyloclade_)),]

text.size <- 6.5
bar.width <- 0.425

centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$feature_ == "centromere_nt_"),]

pdf(paste0("draft_plots_/figure_3_panel_b_.pdf"),onefile = T, width = 1.45,height = 12.15)
ggplot(centromere_dissection_data_gw_all_) +
  geom_col(aes(x = max(centromere_dissection_data_gw_all_$total_length_)/1000000, y = dtolID_), fill = "grey95", color = NA, width = .5) +
  geom_col(aes(x = total_length_/1000000, y = dtolID_, fill = feature_), color = NA, width = .5) +
  scale_fill_manual(values = c("centromere_nt_" = "black")) + 
  facet_grid(phyloclade_ ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(expand = expansion(add = .85)) +
  guides(color = "none") +
  coord_cartesian(expand = F,clip = "on") +
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
        strip.text.y = element_blank(),
        text = element_text(size = text.size))
dev.off()

