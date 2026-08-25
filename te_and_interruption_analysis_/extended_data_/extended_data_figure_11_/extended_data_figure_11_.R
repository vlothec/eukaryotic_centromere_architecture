#!/bin/Rscript

# one tree
options(scipen = 999)

library(cowplot); library(ggnewscale); library(ggpubr); library(ggtree); library(ggtreeExtra); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

rm(list = ls())

tree_file_ <- "all_plant.Gypsy_Athila.retree2.maxiterate0.mafft.fasttree"
tree <- read.tree(file = paste0(tree_file_))
tree$tip.label[grepl("_sce",tree$tip.label)]
tree <- treeio::root(tree,outgroup = tree$tip.label[grepl("_sce",tree$tip.label)])

tree$tip.label <- {tree$tip.label %>% gsub("\\&","_",.)}

data_for_tree_ <- "all_plant.Gypsy_Athila.retree2.maxiterate0.mafft.fasttree_data_for_tree_.tsv"
q <- read.delim(data_for_tree_,sep = "\t",header = T) %>%
  data.frame()

q$ID <- q$rownames.tag
q <- q %>%
  mutate(plant_type_ = case_when(gsub("([a-zA-Z]).*","\\1",dtolID_) == "d" ~ "dicot",
                                 gsub("([a-zA-Z]).*","\\1",dtolID_) == "r" ~ "dicot",
                                 gsub("([a-zA-Z]).*","\\1",dtolID_) == "l" ~ "monocot",
                                 gsub("([a-zA-Z]).*","\\1",dtolID_) == "c" ~ "bryophyte",
                                 .default = NA_character_))

# tree aesthetics
layout.i <- "fan"
offset.v <- -0.05
p <- ggtree(tree, layout = layout.i,size=.075, open.angle=3.5) +
  geom_nodelab(aes(label = node))

pdf(paste0(tree_file_,"_single_plot_node_numbers_.pdf"),height = 30, width = 30)
print(p)
dev.off()


if(tree_file_ == tree_file_){
  relevant.nodes <- c(48348,57269,57139,57270,57137,54042,57271,59890,60036,53165,52540)
  d <- p$data
  d <- d[!d$isTip,] %>%
    arrange(x)
  d$label <- as.numeric(d$label)
  d$label <- round(d$label * 100)
  d <- d[which(d$node %in% relevant.nodes),]
  d <- d[which(d$label >= 90),]
}else{
  d <- p$data
  d <- d[!d$isTip,] %>%
    arrange(x)
  d$label <- as.numeric(d$label)
  d$label <- round(d$label * 100)
  # following line will pull out the 20 closest node labels to root
  d <- d[which(d$node %in% d$node[c(1:20)]),]
}


# define colors
# aesthetics
{
  # define colors
  plant_type_color_ <- c("monocot" = "#2EC4B6",
                         "dicot" = "#D4AF37",
                         "bryophyte" = "#FF6B6B")
  
  plant_type_color_ <- c("monocot" = "black",
                         "dicot" = "grey45",
                         "bryophyte" = "grey90")
  
  dark_col_v_ <- c("Satellite" = "#EA346F",
                   "Holocentric" = "#447B3B",
                   "Transposon" = "#4E87F7",
                   "Satellite_Transposon" = "#7943E3",
                   "Monocentric_unspecified" = "darkgreen",
                   "Unknown" = "black",
                   "ref" = "black",
                   "non_target_" = "white",
                   "target_" = "grey92")
  
  light_col_v_ <- c("Satellite" = "#fad9e4",
                    "Holocentric" = "#d3f7cd",
                    "Transposon" = "#d4e1fa",
                    "Satellite_Transposon" = "#d2befa",
                    "Monocentric_unspecified" = "lightgreen",
                    "Unknown" = "grey98",
                    "ref" = "grey98",
                    "OUT" = "grey98")
  
}

rownames(q) <- tree$tip.label

sat_ <- paste(q$seqID_header_[which(q$centromere_status_ == "IN" & q$class_ == "Satellite")],collapse = "|")
transp_ <- paste(q$seqID_header_[which(q$centromere_status_ == "IN" & q$class_ == "Transposon")],collapse = "|")
sat_transp_ <- paste(q$seqID_header_[which(q$centromere_status_ == "IN" & q$class_ == "Satellite_Transposon")],collapse = "|")
holo_ <- paste(q$seqID_header_[which(q$centromere_status_ == "IN" & q$class_ == "Holocentric")],collapse = "|")

offset_ <- 0.050
main_ <- "angiosperm ATHILA"
p <- ggtree(tree, layout = layout.i,size=.075, open.angle=3.5) +
  geom_tiplab(align = TRUE, aes(subset=grepl(sat_,label,fixed=FALSE)==TRUE,label = ""),color = light_col_v_[grep("^Satellite$",names(light_col_v_))], linetype = "solid", linesize = 0.05,offset = offset_) +
  geom_tiplab(align = TRUE, aes(subset=grepl(transp_,label,fixed=FALSE)==TRUE,label = ""),color = light_col_v_[grep("^Transposon$",names(light_col_v_))], linetype = "solid", linesize = 0.05,offset = offset_) +
  geom_tiplab(align = TRUE, aes(subset=grepl(sat_transp_,label,fixed=FALSE)==TRUE,label = ""),color = light_col_v_[grep("Satellite_Transposon",names(light_col_v_))], linetype = "solid", linesize = 0.05,offset = offset_) +
  geom_tiplab(align = TRUE, aes(subset=grepl(holo_,label,fixed=FALSE)==TRUE,label = ""),color = light_col_v_[grep("Holocentric",names(light_col_v_))], linetype = "solid", linesize = 0.05,offset = offset_) +
  # geom_tiplab(align = TRUE, aes(subset=grepl(paste(q$seqID_header_[which(q$centromere_status_ == "IN" & q$class_ == "Unknown")],collapse = "|"),label,fixed=FALSE)==TRUE,label = ""),color = light_col_v_[grep("Unknown",names(light_col_v_))], linetype = "solid", linesize = 0.05,offset = .2) +
  geom_label(data=d, aes(label=label), color = "black",size.unit = "pt",size = 4.5,label.padding = unit(0.125, "lines"), label.r = unit(0.00,"lines"),label.size = 0,fill = "transparent", alpha = .75) +
  labs(title = paste0(main_, "\nn: ",nrow(q))) +
  theme(text = element_text(size = 3.75))

p.1 <- p +
  geom_fruit(data = q,
             geom = "geom_tile",
             offset = 0,
             mapping = aes(y = ID, x = hmm_comb_, fill = class_),
             pwidth=0.05,
             color = "transparent",
             axis.params=list(),
             grid.params=list(vline = F,
                              color = "transparent")) +
  scale_fill_manual(values = dark_col_v_, guide = "none") +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color()


{
  dark_col_v_ <- c("Satellite" = "#EA346F",
                   "Holocentric" = "#447B3B",
                   "Transposon" = "#4E87F7",
                   "Satellite_Transposon" = "#7943E3",
                   "Monocentric_unspecified" = "darkgreen",
                   "Unknown" = "black",
                   "ref" = "black",
                   "non_target_" = "grey95",
                   "target_" = "black")
  dark_col_v_v2_ <- c("Satellite" = "#EA346F",
                   "Holocentric" = "#447B3B",
                   "Transposon" = "#4E87F7",
                   "Satellite_Transposon" = "#7943E3",
                   "Monocentric_unspecified" = "darkgreen",
                   "Unknown" = "black",
                   "ref" = "black",
                   "non_target_" = "transparent",
                   "target_" = "black")
}

dtolID_list_ <- q$dtolID_ %>% unique()
# independent invasion + drGeuUrba1
dtolID_list_ <- dtolID_list_[grepl("laApoDist1|drGeuRiva1|drGeuUrba1|rosCan_S27_v1|ddBarVulg1|lpLuzSylv1|lpJunEffu1|lpSchTriq1|lpLuzPall1|ddPopNigr1|drPotAnse1|lpCarDepa1|dmThaMinu1|ddAraThal",dtolID_list_)]

offset_ <- 0.025
p.3 <- p.1
for(i_ in 1:length(dtolID_list_)){
  dtolID_ <- dtolID_list_[i_]
  print(dtolID_)
  print(Sys.time())
  q$traget_species_ <- ifelse(q$dtolID_ == dtolID_ & q$centromere_status_ == "IN", q$class_, "non_target_")

  p.4 <- p.3 +
    geom_fruit(data = q[which(q$traget_species_ != "non_target_"),],
               geom = "geom_tile",
               offset = offset_,
               mapping = aes(y = ID, x = hmm_comb_, fill = traget_species_, color = traget_species_),
               pwidth=0.095,
               size = .100,
               axis.params=list(),
               grid.params=list(vline = F,
                                color = "grey95")) +
    scale_fill_manual(values = dark_col_v_, guide = "none") +
    scale_color_manual(values = dark_col_v_v2_, guide = "none") + 
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color()
  
  p.3 <- p.4
}

pdf(paste0(tree_file_,"_single_plot_rings_.pdf"),height = 6.5, width = 6.5)
print(p.3)
dev.off()

png(paste0(tree_file_,"_single_plot_rings_.png"),height = 6.5, width = 6.5,units = "in",res = 1200)
print(p.3)
dev.off()
