library(cowplot); library(ggnewscale); library(ggpubr); library(ggtree); library(ggtreeExtra); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

rm(list = ls())

tree_file_ <- "drGeuUrba1_and_drGeuRiva1_plant_Gypsy_Athila_.globalpair.max100.mafft.fasttree"

tree_ <- read.tree(file = paste0(tree_file_))
tree_$tip.label[grepl("_sce",tree_$tip.label)]
tree_ <- treeio::root(tree_,outgroup = tree_$tip.label[grepl("_sce",tree_$tip.label)])

tree_$tip.label <- {tree_$tip.label %>% gsub("\\&","_",.)}

data_att_list_ <- strsplit(tree_$tip.label,"/")
max_length_ <- max(lengths(data_att_list_))
data_att_list_tmp_ <- lapply(data_att_list_,function(x) {
  length(x) <- max_length_
  x
  })

data_att_ <- do.call(rbind,data_att_list_tmp_) %>%
  as.data.frame()
colnames(data_att_) <- c("seqID_header_","cls_","hmm_comb_","root_name_","def_one_","def_two_")
# if ATHILA root

data_att_$feature_ <- mapply(function(x) gsub(".*#","",x),data_att_$seqID_header_)
data_att_$seqID_ <- mapply(function(x) gsub("(.*)_[0-9]+__.*","\\1",x),data_att_$seqID_header_)
data_att_$start_ <- mapply(function(x) gsub(".*_([0-9]+)__.*#.*","\\1",x),data_att_$seqID_header_) %>% as.numeric()
data_att_$end_ <- mapply(function(x) gsub(".*_[0-9]+__([0-9]+)-.*#.*","\\1",x),data_att_$seqID_header_) %>% as.numeric()
data_att_$length_ <- mapply(function(x,y) (y-x)+1,data_att_$start_,data_att_$end_) %>% as.numeric()
data_att_$method_ <- mapply(function(x) gsub(".*_[0-9]+__[0-9]+-(.*)#.*","\\1",x),data_att_$seqID_header_)

rownames(data_att_) <- tree_$tip.label

drGeuUrba1_cent <- read.delim("drGeuUrba1_cent.coord", header=FALSE)
colnames(drGeuUrba1_cent) <- c("seqID_","start_","end_")
data_att_$cent_status_ <- "OUT"
for(i_ in 1:nrow(drGeuUrba1_cent)){
  root_name_i_ <- "drGeuUrba1"
  seqID_ <- drGeuUrba1_cent$seqID_[i_]
  start_i_ <- drGeuUrba1_cent$start_[i_]
  end_i_ <- drGeuUrba1_cent$end_[i_]
  idx_ <- which(data_att_$seqID_ == seqID_ & data_att_$start_ >= start_i_ & data_att_$end_ <= end_i_ & data_att_$root_name_ == root_name_i_)
  if(length(idx_) > 0){
    data_att_$cent_status_[idx_] <- "IN"
  }
}

drGeuRiva1_cent <- read.delim("drGeuRiva1_cent.coord", header=FALSE)
colnames(drGeuRiva1_cent) <- c("seqID_","start_","end_")
for(i_ in 1:nrow(drGeuRiva1_cent)){
  root_name_i_ <- "drGeuRiva1"
  seqID_ <- drGeuRiva1_cent$seqID_[i_]
  start_i_ <- drGeuRiva1_cent$start_[i_]
  end_i_ <- drGeuRiva1_cent$end_[i_]
  idx_ <- which(data_att_$seqID_ == seqID_ & data_att_$start_ >= start_i_ & data_att_$end_ <= end_i_ & data_att_$root_name_ == root_name_i_)
  if(length(idx_) > 0){
    data_att_$cent_status_[idx_] <- "IN"
  }
}


data_att_ <- data_att_ %>%
  mutate(drGeuRiva1_cent_ = case_when(root_name_ == "drGeuRiva1" & cent_status_ == "IN" ~ "IN",
                                      root_name_ == "drGeuRiva1" & cent_status_ == "OUT" ~ "OUT",
                                      .default = "non_relevant_"),
         drGeuUrba1_cent_ = case_when(root_name_ == "drGeuUrba1" & cent_status_ == "IN" ~ "IN",
                                      root_name_ == "drGeuUrba1" & cent_status_ == "OUT" ~ "OUT",
                                      .default = "non_relevant_"))

# tree aesthetics
layout_i_ <- "circular"
layout_i_ <- "rectangular"
offset_ <- .17

# cols to be plot as heatmap
feature_list_ <- colnames(data_att_)[c(14,15)]
color_map_list_ <- list("drGeuRiva1_cent_" = c("IN" = "#ff4343", "OUT" = "grey85", "non_relevant_" = "transparent"),
                        "drGeuUrba1_cent_" = c("IN" = "#015b87", "OUT" = "grey85", "non_relevant_" = "transparent"))

# plot basic tree 
p_ <- ggtree(tree_, layout = layout_i_,size=.75, open.angle=7.5)

relevant.nodes <- c(3584,5449,5451,5448,5452,5448,5424)
d <- p_$data
d <- d[!d$isTip,] %>%
  arrange(x)
d$label <- as.numeric(d$label)
d$label <- round(d$label * 100)
d <- d[which(d$node %in% relevant.nodes),]
d <- d[which(d$label >= 90),]

p_ <- p_ +
  geom_label(data=d, aes(label=label), color = "black",size.unit = "pt",size = 3.25,label.padding = unit(0.125, "lines"), label.r = unit(0.00,"lines"),label.size = 0,fill = "transparent", alpha = .75)

if (is.null(feature_list_)){
  return (p_)
}

p2_ <- p_
if(length(feature_list_)>1){
  for(feature in 1:length(feature_list_)){
    
    feature_i_ <- feature_list_[feature]
    color_map_i_ <- color_map_list_[[which(grepl(feature_i_,names(color_map_list_)) == T)]]
    
    data_att_i_ <- data_att_[,which(str_detect(colnames(data_att_),feature_i_))] %>% as.data.frame()
    rownames(data_att_i_) <- rownames(data_att_)
    
    feature_i_type_ <- class(data_att_i_[,1])
    n_levels_i_ <- length(unique(data_att_i_[,1]))
    
    p2_ <- p2_ +
      new_scale_fill() +
      new_scale_color()
    
    p2_ <- gheatmap(p2_, data_att_i_, offset=feature*offset_, width=.095,
                   colnames_angle=0, colnames_offset_y = .05,color = NULL)
    
    if (feature_i_type_ == 'numeric'){
      p2_ <- p2_ +
        scale_color_brewer(type="div", palette=cmap)
    } else {
      p2_ <- p2_ +
        scale_fill_manual(values=color_map_i_, name=feature_i_,na.value = "transparent")
    }
    
  }
}

p_feature_ <- p2_
p_feature_ <- p_feature_ +
  theme(legend.position = "none")

# Extended Data Figure 11
{
  data_att_ %>% str()
  
  data_att_$rownames_tag_ <- rownames(data_att_)
  
  data_att_$rownames_tag_ <- data_att_$rownames_tag_ %>%
    {factor(.,levels = unique(data_att_$rownames_tag_))}
  data_att_ <- data_att_[with(data_att_,rev(order(rownames_tag_))),]
    
  tree.order.chr <- get_taxa_name(p_) 
  data_att_$rownames_tag_ <- data_att_$rownames_tag_ %>% factor(., levels = tree.order.chr)
  data_att_ <- data_att_[order(data_att_$rownames_tag_),]
  
  seqID_order_ <- unique(data_att_$seqID_)[order(as.numeric(gsub("[^0-9]+","",unique(data_att_$seqID_))))]
  seqID_order_ <- c(as.character(seqID_order_[grep("Chr",seqID_order_)]),as.character(seqID_order_[grep("SUPER",seqID_order_)]),"M34549_1_sce")
  
  data_att_$seqID_ <- data_att_$seqID_ %>%
    {factor(.,levels = seqID_order_)}
  
  approx_karyo_ <- data_att_ %>%
    group_by(seqID_) %>%
    summarise(min_start_ = 0,
              max_end_ = max(end_,na.rm = T))
  
  data_att_ <- data_att_[with(data_att_,order(rownames_tag_,seqID_)),]
  
  p_chr_feature_ <- ggplot(data_att_[which(data_att_$root_name_ != "unk"),]) +
    geom_rect(data = approx_karyo_[which(approx_karyo_$seqID_ != "M34549_1_sce"),], aes(xmin = min_start_, xmax = max_end_, ymin = -Inf, ymax = Inf), color = NA, fill = "grey96") +
    geom_point(data = data_att_[which(data_att_$root_name_ != "unk"),],aes(x = start_, y = rownames_tag_, color = root_name_, alpha = cent_status_), stroke = 0, size = 3.75) +
    scale_color_manual(values = c("drGeuRiva1" = "#ff4343", "drGeuUrba1" = "#015b87")) +
    scale_alpha_manual(values = c("IN" = 1, "OUT" = .245)) +
    scale_y_discrete(limits = rev) +
    # facet_grid(. ~ seqID_) +
    facet_grid(. ~ seqID_, scales = "free_x",space = "free_x") +
    coord_cartesian(clip = "on") +
    theme(panel.background = element_rect(fill = "grey92", color = NA),
          panel.spacing.x = unit(5.05,"pt"),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          panel.ontop = F)
  }


pdf(file = paste0("drGeuUrba1_and_drGeuRiva1_plant_Gypsy_Athila_.globalpair.max100.mafft.fasttree_plot_",layout_i_,"_plus_chr_.pdf"),width = 36,height = 12)
plot_grid(p_feature_,p_chr_feature_,ncol = 2,nrow = 1,rel_widths = c(.125,1),align = "hv",axis = "tb")
dev.off()
