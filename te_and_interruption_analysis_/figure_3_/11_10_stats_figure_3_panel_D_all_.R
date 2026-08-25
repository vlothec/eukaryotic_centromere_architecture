# panel D stats

library(tidyverse)
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

LTR_piden_data_gw_aggregated_$centromere_status_[which(is.na(LTR_piden_data_gw_aggregated_$centromere_status_) == TRUE)] <- "arm"
LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_[grepl("structural",LTR_piden_data_gw_aggregated_$attributes_),]


counts.f <- function(x){
  hist(x, breaks=seq(0,1,0.001), plot=F)$counts %>%
    {paste0(.,collapse = ";")}
}
mids.f <- function(x){
  hist(x, breaks=seq(0,1,0.001), plot=F)$mids %>%
    {paste0(.,collapse = ";")}
}

tally.data_LTR_piden_data_gw_aggregated_ <- LTR_piden_data_gw_aggregated_ %>%
  group_by(dtolID_,centromere_status_) %>%
  summarise(obs.values = paste0(ltr_identity, collapse = ";"),
            counts.hist = counts.f(ltr_identity),
            mids.hist = mids.f(ltr_identity),
            lower.whisker.bp = boxplot.stats(ltr_identity)$stats[1],
            lower.hinge.bp = boxplot.stats(ltr_identity)$stats[2],
            median.bp = boxplot.stats(ltr_identity)$stats[3],
            upper.hinge.bp = boxplot.stats(ltr_identity)$stats[4],
            upper.whisker.bp = boxplot.stats(ltr_identity)$stats[5],
            n.obs = boxplot.stats(ltr_identity)$n)
tally.data_LTR_piden_data_gw_aggregated_$gs.i <- NA
tally.data_LTR_piden_data_gw_aggregated_ %>% ungroup() %>% str()

write.table(tally.data_LTR_piden_data_gw_aggregated_,paste0("tally_data_LTR_piden_data_gw_aggregated_"),row.names = F, col.names = T, sep = "\t", quote = F)


# extract obs.value
tally.data_LTR_piden_data_gw_aggregated_obs.values <- tally.data_LTR_piden_data_gw_aggregated_ %>%
  select(dtolID_,centromere_status_,obs.values) %>%
  mutate(obs.values.list = lapply(strsplit(obs.values, ";", TRUE), as.numeric))

tally.data_LTR_piden_data_gw_aggregated_obs.values <- do.call('rbind', do.call('Map', c(data.frame, tally.data_LTR_piden_data_gw_aggregated_obs.values)))
tally.data_LTR_piden_data_gw_aggregated_obs.values <- tally.data_LTR_piden_data_gw_aggregated_obs.values[,-c(3)]

taxa.groups <- LTR_piden_data_gw_aggregated_[,c(13,15)] %>%
  distinct()

tally.data_LTR_piden_data_gw_aggregated_obs.values_red_ <- tally.data_LTR_piden_data_gw_aggregated_obs.values %>%
  group_by(dtolID_,centromere_status_) %>%
  summarise(mean.obs.values = mean(obs.values.list,na.rm = T))

tally.data_LTR_piden_data_gw_aggregated_obs.values_red_ <- left_join(tally.data_LTR_piden_data_gw_aggregated_obs.values_red_,taxa.groups,by = "dtolID_")

# hist
{
  myBreaks <- seq(0,1,0.001)
  data_ <- tally.data_LTR_piden_data_gw_aggregated_obs.values
  data_ <- tally.data_LTR_piden_data_gw_aggregated_obs.values_red_
  hist(data_[which(data_$centromere_status_ == "arm"),c(3)],
       breaks = myBreaks,
       col = "lightskyblue",
       border = NA,
       xlim = c(0.80,1))
  hist(data_[which(data_$centromere_status_ != "arm"),c(3)],
       breaks = myBreaks,
       col = "blue",
       border = NA,
       xlim = c(0.80,1),
       add = T)
}

# following lines are intended to assess that variance differs
{
  tally.data_LTR_piden_data_gw_aggregated_obs.values_red_ %>%
    group_by(centromere_status_) %>%
    summarise(avg.mean.obs.values = mean(mean.obs.values),
              median.mean.obs.values = median(mean.obs.values),
              st.mean.obs.values = sd(mean.obs.values))
  
  tally.data_LTR_piden_data_gw_aggregated_obs.values_red_ %>%
    group_by(centromere_status_,phyloclade_) %>%
    summarise(avg.mean.obs.values = mean(mean.obs.values),
              median.mean.obs.values = median(mean.obs.values),
              st.mean.obs.values = sd(mean.obs.values))
  
  boxplot(mean.obs.values ~ centromere_status_,
          data = tally.data_LTR_piden_data_gw_aggregated_obs.values_red_)
  
  var.test(mean.obs.values ~ centromere_status_,
           data = tally.data_LTR_piden_data_gw_aggregated_obs.values_red_)
  
}

sink("stats_figure_3_panel_D_all.txt", split = TRUE)

print("all")

t.i <- t.test(tally.data_LTR_piden_data_gw_aggregated_obs.values_red_$mean.obs.values ~ tally.data_LTR_piden_data_gw_aggregated_obs.values_red_$centromere_status_,
              paired = FALSE)
print(t.i)
rm(t.i)

for(tax in 1:3){
  tax.i <- unique(tally.data_LTR_piden_data_gw_aggregated_obs.values_red_$phyloclade_)[tax]
  tally.data_LTR_piden_data_gw_aggregated_obs.values_red_subset_ <- tally.data_LTR_piden_data_gw_aggregated_obs.values_red_[which(tally.data_LTR_piden_data_gw_aggregated_obs.values_red_$phyloclade_ == tax.i),] 
  print(tax.i)
  t.i <- t.test(tally.data_LTR_piden_data_gw_aggregated_obs.values_red_subset_$mean.obs.values ~ tally.data_LTR_piden_data_gw_aggregated_obs.values_red_subset_$centromere_status_,
                paired = FALSE)
  print(t.i)
  rm(t.i)
}

sink()
