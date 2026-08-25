library(dplyr)
library(data.table)
options(scipen = 999)

rm(list = ls())

# intact+frag
path_to_anno_data_ <- ""
list_anno_file_ <- list.files(path = path_to_anno_data_, pattern = ".cls.tsv$")

path_to_array_data_ <- ""
list_array_file_ <- list.files(path = path_to_array_data_, pattern = "_arrays_expanded.csv$")
list_karyotype_file_ <- list.files(path = path_to_array_data_, pattern = "_karyotype_ext_$",recursive = T)

list_dtolID_ <- list_array_file_ %>%
  {gsub("_arrays_expanded.csv","",.)}

anno_data_sum_all_ <- data.frame()
counter_ <- 0
for(i_ in 1:length(list_dtolID_)){
  root_name_ <- list_dtolID_[i_]
  print(root_name_)
  
  anno_data_ <- fread(paste0(path_to_anno_data_,list_anno_file_[grepl(gsub("\\.hap.*|\\.fa.*","",root_name_),list_anno_file_)]), header = T, comment.char = "",sep = "\t") %>%
    data.frame()
  
  anno_data_$feature <- mapply(function(x) gsub(".*#","",x),anno_data_$X.TE)
  anno_data_$seqID <- mapply(function(x) gsub("(.*):.*","\\1",x),anno_data_$X.TE)
  anno_data_$start_ <- mapply(function(x) gsub(".*:([0-9]+)\\.\\..*#.*","\\1",x),anno_data_$X.TE) %>% as.numeric()
  anno_data_$end_ <- mapply(function(x) gsub(".*:.*\\.\\.([0-9]+)-.*#.*","\\1",x),anno_data_$X.TE) %>% as.numeric()
  anno_data_$length_ <- mapply(function(x,y) (y-x)+1,anno_data_$start_,anno_data_$end_) %>% as.numeric()
  anno_data_$method_ <- mapply(function(x) gsub(".*:.*\\.\\.[0-9]+-(.*)#.*","\\1",x),anno_data_$X.TE)
  anno_data_$n_hmm_ <- mapply(function(x) length(unlist(gregexpr("\\|",x))),anno_data_$Domains)
  anno_data_$n_hmm_[which(anno_data_$Domains == "none")] <- 0
  
  # clean the data little bit
  # keep only sites with EDTA == TEsorter
  test.consistency <- mapply(function(x,y) ifelse(length(grep(x,y)) > 0,T,F),anno_data_$Order,anno_data_$feature)
  anno_data_ <- anno_data_[test.consistency,]
  nrow(anno_data_)
  
  # load arrays data
  arrays_data_ <- read.table(paste0(path_to_array_data_,list_array_file_[grepl(root_name_,list_array_file_)]), header = T, sep = ",")
  if(ncol(arrays_data_) > 17){
    arrays_data_$filter_comb_ <- apply(arrays_data_[,c(17,18)],1,function (x) paste(x[1],x[2],sep = "_"))
    filter_comb_i_ <- "TRUE_FALSE"
    # keep relevant arrays
    if(!is.null(filter_comb_i_)){
      arrays_data_ <- arrays_data_[which(arrays_data_$filter_comb_ %in% filter_comb_i_),]
    }
  }else{
    print("Holo!")
    counter_ <- counter_+1
    next
  }
  
  # load karyo data
  karyotype_file_ <- list_karyotype_file_[grepl(gsub("\\..*","",root_name_),list_karyotype_file_)]
  if(length(karyotype_file_) > 0){
    karyotype_data_ <- read.table(paste0(path_to_array_data_,"/",karyotype_file_),sep = "\t",header = F)
    colnames(karyotype_data_) <- c("ENA_","ChrOnly","length_")
    karyotype_data_$new_ENA_ <- karyotype_data_$ENA_ %>%
      {gsub("\\|","_",.)}
  }
    
  seqID_type_ <- which(anno_data_$seqID %in% karyotype_data_$ENA_) %>% length()
  if(seqID_type_ == 0){
    if(length(grep("ENA_",anno_data_$seqID)) > 0){
      colnames(karyotype_data_)[4] <- "seqID"
      anno_data_ <- left_join(anno_data_,karyotype_data_[,c(1,4)])
      anno_data_$seqID <- anno_data_$ENA_
      anno_data_ <- anno_data_[,-ncol(anno_data_)]
    }else{
      colnames(karyotype_data_)[2] <- "seqID"
      anno_data_ <- left_join(anno_data_,karyotype_data_[,c(1,2)])
      anno_data_$seqID <- anno_data_$ENA_
      anno_data_ <- anno_data_[,-ncol(anno_data_)]
    }
  }
  
  anno_data_$centromere_status_ <- "OUT"
  anno_data_$arrayID_ <- NA_character_
  for(k_ in 1:nrow(arrays_data_)){
    print(k_)
    array_seqID_ <- arrays_data_$seqID[k_]
    array_start_ <- arrays_data_$start[k_]
    array_end_ <- arrays_data_$end[k_]
    array_arrayID_ <- arrays_data_$arrayID[k_]
    relevant_loci_ <- which(anno_data_$seqID %in% array_seqID_ & anno_data_$start_ >= array_start_ & anno_data_$end_ <= array_end_)
    anno_data_$centromere_status_[relevant_loci_] <- "IN"
    anno_data_$arrayID_[relevant_loci_] <- array_arrayID_
  }
  
  anno_data_sum_ <- anno_data_ %>%
    group_by(Superfamily,Clade,centromere_status_) %>%
    summarise(n_ = n(),
              space_ = sum(length_),
              dtolID_ = root_name_)
  
  anno_data_sum_all_ <- rbind(anno_data_sum_all_,anno_data_sum_)
}


write.table(anno_data_sum_all_, "all.edta_filtered_reassigned.rexdb-plant.cls.tsv.tally_", sep = "\t", col.names = T, row.names = F, quote = F)

# intact+frag
# Extended Data Figure 9
{
  library(tidyverse)
  library(gridExtra)
  rm(list = ls())

  data.i <- read.table("all.edta_filtered_reassigned.rexdb-plant.cls.tsv.tally_",header = T)
  colnames(data.i) <- c("Superfamily","Clade","cent.occ","n.count","nt.sum","species")
  
  data.i %>% str
  # remove mixture elements
  data.i <- data.i[-which(data.i$Clade == "mixture"),]
  data.i <- data.i[-which(data.i$Clade == "unknown"),]
  data.i <- data.i[-which(data.i$Clade == "chromo-outgroup"),]
  data.i <- data.i[-which(data.i$Clade == "chromo-unclass"),]
  data.i <- data.i[-which(data.i$Clade == "non-chromo-outgroup"),]
    
  data.i <- data.i %>%
    group_by(species,cent.occ) %>%
    mutate(nt.sum.prop = nt.sum/sum(nt.sum),
           n.count.prop = n.count/sum(n.count),
           com.label = paste0(gsub("\\..*","",species)," - (n: ",sum(n.count),")"))
  
  superfamily_and_lineage_list_ <- read.table("superfamily_and_lineage_list_", sep = " ", header = F)
  colnames(superfamily_and_lineage_list_) <- c("Superfamily","Clade")
  superfamily_and_lineage_list_$Clade[which(superfamily_and_lineage_list_$Clade == "unknown")] <- paste0(superfamily_and_lineage_list_$Clade[which(superfamily_and_lineage_list_$Clade == "unknown")],"_",superfamily_and_lineage_list_$Superfamily[which(superfamily_and_lineage_list_$Clade == "unknown")])
  
  data.i$Clade[which(data.i$Clade == "unknown")] <- paste0(data.i$Clade[which(data.i$Clade == "unknown")],"_",data.i$Superfamily[which(data.i$Clade == "unknown")])
  
  data.i$Clade <- data.i$Clade %>%
    {factor(.,levels = unique(data.i$Clade))}
  data.i <- data.i[which(data.i$Superfamily %in% c("Copia","Gypsy")),]
  data.i$Superfamily <- data.i$Superfamily %>%
    {factor(.,levels = c("Copia","Gypsy"))}
  
  species.list <- unique(data.i$species)
  data.i <- data.i[with(data.i,order(Clade,Superfamily)),]
  
  text.size <- 12
  p.1 <- list()
  p.2 <- list()
  for(sp in 1:length(species.list)){
    sp.i <- species.list[sp]
    linewidth <- 0.5
    
    if(sp == 1){
      font.col <- "black"
    }else{
      font.col <- "transparent"
    }
    
    if(length(which(data.i$species == sp.i & data.i$cent.occ == "IN")) > 0){
      p.1[[length(p.1) + 1]] <- ggplot(data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$cent.occ == "IN")]))) +
        geom_col(aes(x = species, y = 1), fill = "grey95",color = "transparent", width=4.5) +
        geom_col(data = data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$species == sp.i & data.i$cent.occ == "IN")])),aes(x = species, y = 1), fill = "transparent",color = "black", width=4.5) +
        geom_col(data = data.i[which(data.i$species == sp.i & data.i$cent.occ == "IN"),], aes(x = species, y = nt.sum.prop,fill = Superfamily),color = "black", width=4.5) +
        geom_text(data = data.i[which(data.i$species == sp.i & data.i$cent.occ == "IN"),], aes(x = species, y = 0, label = paste0("n: ",n.count,"\n",round(nt.sum/1000,digits = 1)," kb")), size.unit = "pt", size = 5.15) +
        scale_color_manual(values = c("transparent","black")) +
        scale_fill_manual(values = c("Copia" = "#D83067","Gypsy" = "#487CE3")) +
        coord_polar("y", start=0) +
        facet_grid(com.label~Clade,labeller = labeller(.rows = ~ gsub(" - ","\n",.))) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.line.y = element_line(colour = NA, linewidth = linewidth),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = NA, colour = NA, linewidth = linewidth*2),
              panel.ontop = F,
              plot.margin = margin(t=0,r=0,b=0,l=0),
              plot.title = element_text(hjust = 0),
              strip.text.y = element_text(angle = 0,hjust = 0,size = 9),
              strip.text.x = element_text(angle = 0,hjust = .5,vjust = .5,size = 9,color = font.col),
              text = element_text(size = text.size))
      
      p.2[[length(p.2) + 1]] <- ggplot(data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$cent.occ == "IN")]))) +
        geom_col(aes(x = species, y = 1), fill = "grey95",color = "transparent", width=4.5) +
        geom_col(data = data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$species == sp.i & data.i$cent.occ == "IN")])),aes(x = species, y = 1), fill = "transparent",color = "black", width=4.5) +
        geom_col(data = data.i[which(data.i$species == sp.i & data.i$cent.occ == "IN"),], aes(x = species, y = n.count.prop,fill = Superfamily),color = "black", width=4.5) +
        geom_text(data = data.i[which(data.i$species == sp.i & data.i$cent.occ == "IN"),], aes(x = species, y = 0, label = paste0("n: ",n.count,"\n",round(nt.sum/1000,digits = 1)," kb")), size.unit = "pt", size = 5.15) +
        scale_color_manual(values = c("transparent","black")) +
        scale_fill_manual(values = c("Copia" = "#D83067","Gypsy" = "#487CE3")) +
        coord_polar("y", start=0) +
        facet_grid(com.label~Clade,labeller = labeller(.rows = ~ gsub(" - ","\n",.))) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.line.y = element_line(colour = NA, linewidth = linewidth),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = NA, colour = NA, linewidth = linewidth*2),
              panel.ontop = F,
              plot.margin = margin(t=0,r=0,b=0,l=0),
              plot.title = element_text(hjust = 0),
              strip.text.y = element_text(angle = 0,hjust = 0,size = 9),
              strip.text.x = element_text(angle = 0,hjust = .5,vjust = .5,size = 9,color = font.col),
              text = element_text(size = text.size))
    }else{
      print(sp.i)
    }
  }
  
  
  pdf("intact_plus_frag_pie_test_no_unk_.pdf",height = length(p.2),width = length(unique(data.i$Clade[which(data.i$cent.occ == "IN")])),onefile = T)
  print(do.call("grid.arrange", c(p.2, ncol=1, nrow=length(p.1))))
  dev.off()
  
}


# Figure 3 Panel d

library(tidyverse)
library(gridExtra)
rm(list = ls())

data.i <- read.table("all.edta_filtered_reassigned.rexdb-plant.cls.tsv.tally_17072026_",header = T)
colnames(data.i) <- c("Superfamily","Clade","cent.occ","n.count","nt.sum","species")

data.i %>% str
# remove mixture elements
data.i <- data.i[-which(data.i$Clade == "mixture"),]
data.i <- data.i[-which(data.i$Clade == "unknown"),]
data.i <- data.i[-which(data.i$Clade == "chromo-outgroup"),]
data.i <- data.i[-which(data.i$Clade == "chromo-unclass"),]
data.i <- data.i[-which(data.i$Clade == "non-chromo-outgroup"),]

data.i <- data.i %>%
  group_by(species,cent.occ) %>%
  mutate(nt.sum.prop = nt.sum/sum(nt.sum),
         n.count.prop = n.count/sum(n.count),
         com.label = paste0(gsub("\\..*","",species)," - (n: ",sum(n.count),")"))


superfamily_and_lineage_list_ <- read.table("superfamily_and_lineage_list_", sep = " ", header = F)
colnames(superfamily_and_lineage_list_) <- c("Superfamily","Clade")
superfamily_and_lineage_list_$Clade[which(superfamily_and_lineage_list_$Clade == "unknown")] <- paste0(superfamily_and_lineage_list_$Clade[which(superfamily_and_lineage_list_$Clade == "unknown")],"_",superfamily_and_lineage_list_$Superfamily[which(superfamily_and_lineage_list_$Clade == "unknown")])

data.i$Clade[which(data.i$Clade == "unknown")] <- paste0(data.i$Clade[which(data.i$Clade == "unknown")],"_",data.i$Superfamily[which(data.i$Clade == "unknown")])

target_lineages_ <- c("Ale","Angela","Bianca","SIRE","TAR","Tork","Galadriel","CRM","Tekay","Reina","Athila","Ogre","Retand")

data.i <- data.i[which(data.i$Clade %in% target_lineages_),]
data.i$Clade <- data.i$Clade %>%
  {factor(.,levels = target_lineages_)}
data.i <- data.i[which(data.i$Superfamily %in% c("Copia","Gypsy")),]
data.i$Superfamily <- data.i$Superfamily %>%
  {factor(.,levels = c("Copia","Gypsy"))}

species.list <- unique(data.i$species)
data.i <- data.i[with(data.i,order(Clade,Superfamily)),]

data.i <- data.i %>%
  group_by(Clade, cent.occ) %>%
  mutate(species_aff_ = n())

pdf("figure_3_panel_d_.pdf",height = 3.5,width = 3.5,onefile = T)
ggplot(data.i[which(data.i$cent.occ == "IN"),]) +
  geom_boxplot(aes(x = Clade, y = nt.sum.prop),fill = "white",outliers = F,width = .625,linewidth = .265) +
  geom_jitter(aes(x = Clade, y = nt.sum.prop, color = Superfamily, size = nt.sum/1000000),stroke = 0, alpha = .5, height = 0,width = .15) +
  geom_text(data = distinct(data.i[which(data.i$cent.occ == "IN"),c(1,2,3,10)]), aes(x = Clade, y = 1, label = paste0(" ",species_aff_," ")), size.unit = "pt",size = 4.85) +
  scale_color_manual(values = c("Copia" = "#D83067","Gypsy" = "#487CE3")) +
  guides(colour = "none") +
  coord_cartesian(clip = "off") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        panel.background = element_rect(fill = "grey96",color = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "white"),
        strip.background = element_blank(),
        text = element_text(size = 7.5))
dev.off()  
