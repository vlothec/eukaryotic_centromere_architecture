# panel c stats

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


taxa.groups <- TEanno_data_all_method_data_gw_aggregated_ %>%
  select(dtolID_,phyloclade_) %>%
  distinct()

TEanno_data_all_method_data_gw_aggregated_red_ <- complete(TEanno_data_all_method_data_gw_aggregated_[which(TEanno_data_all_method_data_gw_aggregated_$Method == "structural"),c(7,1,6)],dtolID_,centromere_status_,fill = list(value = 0))
TEanno_data_all_method_data_gw_aggregated_red_$prop_centromere_status_length_[which(is.na(TEanno_data_all_method_data_gw_aggregated_red_$prop_centromere_status_length_) == TRUE)] <- 0
TEanno_data_all_method_data_gw_aggregated_red_ <- left_join(TEanno_data_all_method_data_gw_aggregated_red_,taxa.groups,by = "dtolID_")



# following lines are intended to assess variance differs
# and thus the correct test to apply is Welch Two Sample t-test
{
  TEanno_data_all_method_data_gw_aggregated_red_ %>%
    group_by(centromere_status_) %>%
    summarise(avg.prop.nt.sum = mean(prop_centromere_status_length_),
              median.prop.nt.sum = median(prop_centromere_status_length_),
              st.prop.nt.sum = sd(prop_centromere_status_length_))
  
  boxplot(prop_centromere_status_length_ ~ centromere_status_,
          data = TEanno_data_all_method_data_gw_aggregated_red_)
  
  var.test(prop_centromere_status_length_ ~ centromere_status_,
           data = TEanno_data_all_method_data_gw_aggregated_red_)
}

sink("stats_figure_3_panel_C_all.txt", split = TRUE)

print("all")

t.i <- t.test(TEanno_data_all_method_data_gw_aggregated_red_$prop_centromere_status_length_ ~ TEanno_data_all_method_data_gw_aggregated_red_$centromere_status_,
              paired = FALSE)
print(t.i)

rm(t.i)
for(tax in 1:3){
  tax.i <- unique(TEanno_data_all_method_data_gw_aggregated_red_$phyloclade_)[tax]
  TEanno_data_all_method_data_gw_aggregated_red_subset_ <- TEanno_data_all_method_data_gw_aggregated_red_[which(TEanno_data_all_method_data_gw_aggregated_red_$phyloclade_ == tax.i),] 
  print(tax.i)
  t.i <- t.test(TEanno_data_all_method_data_gw_aggregated_red_subset_$prop_centromere_status_length_ ~ TEanno_data_all_method_data_gw_aggregated_red_subset_$centromere_status_,
                paired = FALSE)
  print(t.i)
  rm(t.i)
}

sink()
