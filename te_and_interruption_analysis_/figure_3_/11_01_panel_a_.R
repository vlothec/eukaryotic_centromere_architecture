
# panel A
library(tidyverse)

centromere_dissection_data_gw_all_ <- read.table("centromere_dissection_data_gw_mono_aggregated_", sep = "\t", header = T)
dtolID_v_ <- centromere_dissection_data_gw_all_$dtolID_ %>% unique()

dtol_info_ <- read.table("DToL_info_.txt",header = F,sep = "\t")
colnames(dtol_info_) <- c("batch_","assembly_","clade_","genus_","species_","genome_size_")
dtol_info_ <- dtol_info_ %>%
  mutate(binomial_name_ = paste0(gsub("([A-Za-z]).*","\\1",genus_),". ",species_),
         dtolID_ = gsub(".fa.*","",assembly_))
centromere_dissection_data_gw_all_ <- left_join(centromere_dissection_data_gw_all_,dtol_info_[,c(7,8)])
binomial_name_ <- dtol_info_$binomial_name_
names(binomial_name_) <- dtol_info_$dtolID_

centromere_dissection_data_gw_all_$phyloclade_ <- centromere_dissection_data_gw_all_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(phyloclade_)),]

centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(centromere_prop_)),]
dtolID_order_ <- centromere_dissection_data_gw_all_$dtolID_[which(centromere_dissection_data_gw_all_$feature_ == "TR_nt_")] %>% unique()
writeLines(dtolID_order_,"dtolID_order_")

centromere_dissection_data_gw_all_$dtolID_ <- centromere_dissection_data_gw_all_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(dtolID_)),]

text.size <- 6.5
bar.width <- 0.425
col_pal_ <- c("plant" = "#A8BA70","chordate" = "#C6517D", "invertebrate" = "#3C3C7A")

pdf(paste0("draft_plots_/figure_3_panel_a_.pdf"),onefile = T, width = 1.85,height = 12.15)
ggplot(centromere_dissection_data_gw_all_) +
  geom_col(data = centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$step_ == 1),], aes(x = centromere_prop_, y = dtolID_), color = NA, fill = "#DDDDDD", width = .5) +
  geom_col(data = centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$step_ >= 2 & centromere_dissection_data_gw_all_$step_ < 4 & centromere_dissection_data_gw_all_$feature_ != "unk_nt_" & centromere_dissection_data_gw_all_$feature_ != "gap_nt_"),], aes(x = centromere_prop_, y = dtolID_, fill = feature_), color = NA, width = .5) +
  facet_grid(phyloclade_ ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("TR_nt_" = "#D83067", "TE_nt_" = "#487CE3", "unk_nt_" = "#DDDDDD", "otherTR_nt_" = "#34c975", "centromere_nt_" = "#DDDDDD")) +
  # scale_x_continuous(labels = ~ ifelse(. < 0, " ",.), expand = expansion(0,0)) +
  scale_y_discrete(expand = expansion(add = .85), labels = binomial_name_) +
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
        strip.text.y = element_blank(),
        text = element_text(size = text.size))
dev.off()
