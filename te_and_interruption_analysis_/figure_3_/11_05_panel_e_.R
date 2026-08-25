# panel E

TE_class_list_updated_ <- read.table("TE_class_list_updated_",sep = "\t",header = F)
colnames(TE_class_list_updated_) <- c("TEanno_cls_","old_TEanno_cls_","TEclass_")

centromere_dissection_data_gw_all_ <- read.table("centromere_dissection_data_gw_mono_aggregated_", sep = "\t", header = T)

# gap space
{
  gap_space_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$feature_ == "gap_nt_"),]
  gap_space_$TEclass_ <- "gap_space_"
  gap_space_$TEanno_cls_ <- "gap_space_"
  max_gap_space_ <- max(gap_space_$total_length_, na.rm = T)
}

centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$step_ == 5 & centromere_dissection_data_gw_all_$TEclass_ != "all"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(is.na(centromere_dissection_data_gw_all_$TEclass_) == F),]


dtolID_order_ <- readLines("dtolID_order_")
setdiff(unique(centromere_dissection_data_gw_all_$dtolID_), dtolID_order_)
setdiff(dtolID_order_,unique(centromere_dissection_data_gw_all_$dtolID_))

centromere_dissection_data_gw_all_$dtolID_ <- centromere_dissection_data_gw_all_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(dtolID_)),]

centromere_dissection_data_gw_all_$phyloclade_ <- centromere_dissection_data_gw_all_$phyloclade_ %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(phyloclade_,dtolID_)),]

dtolID_v_ <- centromere_dissection_data_gw_all_$dtolID_ %>% unique()


centromere_dissection_data_gw_all_ <- left_join(centromere_dissection_data_gw_all_,TE_class_list_updated_[,c(1,3)], by = "TEclass_")


centromere_dissection_data_gw_all_$TEclass_[which(centromere_dissection_data_gw_all_$TEclass_ == "repeat_region")] <- "TE_unclass"
centromere_dissection_data_gw_all_$TEclass_[which(centromere_dissection_data_gw_all_$TEclass_ == "repeat_fragment")] <- "TE_unclass"

centromere_dissection_data_gw_all_$TEanno_cls_[which(centromere_dissection_data_gw_all_$TEanno_cls_ == "repeat_region")] <- "TE_unclass"
centromere_dissection_data_gw_all_$TEanno_cls_[which(centromere_dissection_data_gw_all_$TEanno_cls_ == "repeat_fragment")] <- "TE_unclass"

centromere_dissection_data_gw_all_$TEclass_raw_ <- centromere_dissection_data_gw_all_$TEclass_
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_ %>%
  rowwise() %>%
  mutate(TEclass_ = case_when(length(grep("_TIR_",TEclass_raw_)) > 0 & step_ == 5 ~ "TIR_DDE",
                              length(grep("LINE",TEclass_raw_)) > 0 & step_ == 5 ~ "LINE",
                              length(grep("SINE",TEclass_raw_)) > 0 & step_ == 5 ~ "SINE",
                              length(grep("olinton",TEclass_raw_)) > 0 & step_ == 5 ~ "Polinton",
                              length(grep("Crypton_Tyrosine_Recombinase",TEclass_raw_)) > 0 & step_ == 5 ~ "Crypton",
                              length(grep("Penelope_retrotransposon",TEclass_raw_)) > 0 & step_ == 5 ~ "Penelope",
                              length(grep("DIRS_YR_retrotransposon",TEclass_raw_)) > 0 & step_ == 5 ~ "YR",
                              length(grep("Tyrosine_Recombinase_Elements",TEclass_raw_)) > 0 & step_ == 5 ~ "YR",
                              length(grep("DIRS_YR_retrotransposon",TEclass_raw_)) > 0 & step_ == 5 ~ "YR",
                              length(grep("Bel_Pao_LTR_retrotransposon",TEclass_raw_)) > 0 & step_ == 5 ~ "BEL/Pao",
                              length(grep("non_LTR_retrotransposon",TEclass_raw_)) > 0 & step_ == 5 ~ "other_nonLTR",
                              .default = TEclass_raw_))

# following classes won't be shown in the final plot
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "rRNA_gene"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "DNA_transposon"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "pararetrovirus"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "other_nonLTR"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "Caulimoviridae"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "LTR_retrotransposon"),]
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$TEclass_ != "TE_unclass"),]

table(centromere_dissection_data_gw_all_$TEclass_) %>% data.frame()

centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_ %>%
  group_by(TEclass_, step_, dtolID_) %>%
  mutate(pure_TE_prop_raw_ = pure_TE_prop_) %>%
  mutate(pure_TE_prop_ = sum(pure_TE_prop_raw_)) %>%
  group_by(TEclass_, step_) %>%
  mutate(sum_total_length_ = sum(total_length_))
centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$pure_TE_prop_ > 0.01),]


centromere_dissection_data_gw_all_[,c(15,7,18)] %>%
  distinct() %>%
  arrange(TEanno_cls_) %>%
  group_by(TEanno_cls_) %>%
  mutate(total_sum_total_length_ = sum(sum_total_length_)) %>%
  arrange(total_sum_total_length_,sum_total_length_) 


centromere_dissection_data_gw_all_$TEanno_cls_ <- centromere_dissection_data_gw_all_$TEanno_cls_ %>%
  {factor(.,levels = c("Class_I_LTR","Class_II_DNA_element","Class_I_non_LTR_TPRT","Class_I_non_LTR_other"))}
centromere_dissection_data_gw_all_$TEclass_ <- centromere_dissection_data_gw_all_$TEclass_ %>%
  {factor(.,levels = c("Gypsy_LTR_retrotransposon","Copia_LTR_retrotransposon","Retrovirus","BEL/Pao","TIR_DDE","helitron","Polinton","Crypton","LINE","SINE","Penelope","YR"))}

centromere_dissection_data_gw_all_ <- centromere_dissection_data_gw_all_[with(centromere_dissection_data_gw_all_,order(dtolID_,TEanno_cls_,TEclass_,sum_total_length_)),]

# dummy
{
dummy_ <- data.frame(distinct(centromere_dissection_data_gw_all_[,c(11,12,13)]), TEanno_cls_ = "Class_I_LTR", TEclass_ = "Gypsy_LTR_retrotransposon")
dummy_ <- left_join(data.frame("dtolID_" = dtolID_order_),dummy_)
dummy_$phylotaxa_[which(dummy_$dtolID_ == "idBibMarc1.1")] <- "Diptera"
dummy_$phyloclade_[which(dummy_$dtolID_ == "idBibMarc1.1")] <- "invertebrate"
dummy_$phylotaxa_[which(dummy_$dtolID_ == "fToxJac2.1")] <- "Actinopterygii"
dummy_$phyloclade_[which(dummy_$dtolID_ == "fToxJac2.1")] <- "chordate"
dummy_$TEanno_cls_[which(is.na(dummy_$TEanno_cls_) == T)] <- "Class_I_LTR"
dummy_$TEclass_[which(is.na(dummy_$TEclass_) == T)] <- "Gypsy_LTR_retrotransposon"

dummy_$dtolID_ <- dummy_$dtolID_ %>%
  {factor(.,levels = dtolID_order_)}
}

col_pal_ <- c("plant" = "#A8BA70","chordate" = "#C6517D", "invertebrate" = "#3C3C7A")
col_pal_ <- c("Class_I_LTR" = "#03a9fc", "Class_I_non_LTR_other" = "#fc6b03", "Class_I_non_LTR_TPRT" = "#fc3f65", "Class_II_DNA_element" = "#fcc11e", "TE_unclass" = "grey70")

text.size <- 6.5
linewidth_ <- 0

pdf(paste0("draft_plots_/figure_3_panel_e_.pdf"),onefile = T, width = 3.45,height = 12.275)
ggplot(dummy_) +
  geom_point(data = dummy_, aes(y = dtolID_, x = TEclass_), color = "transparent",  stroke = linewidth_) +
  geom_point(data = centromere_dissection_data_gw_all_[which(centromere_dissection_data_gw_all_$pure_TE_prop_ > 0),], aes(y = dtolID_, x = TEclass_, color = TEanno_cls_, size = pure_TE_prop_, alpha = pure_TE_prop_),stroke = linewidth_) +
  # geom_point(data = gap_space_[which(gap_space_$dtolID_ %in% dtolID_order_),], aes(y = dtolID_, x = TEclass_, size = total_length_/max_gap_space_, alpha = centromere_prop_/max_gap_space_), color = "black",stroke = linewidth_) +
  scale_color_manual(values = col_pal_) +
  scale_size_continuous(range = c(.5,6.05),guide = NULL) +
  scale_alpha_continuous(range = c(.75,.45),guide = NULL) +
  scale_y_discrete(expand = expansion(add = .85)) +
  theme_minimal() +
  facet_grid(phyloclade_ ~ fct_inorder(TEanno_cls_), scales = "free", space = "free") +
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

