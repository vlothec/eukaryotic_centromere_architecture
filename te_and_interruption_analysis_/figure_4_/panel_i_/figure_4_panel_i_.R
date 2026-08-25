# R


# input data as provided in the rep: drGeuUrba1.1_stitched_strata_*.anno
list.files.v <- list.files(path = "data/",pattern = ".anno")

data.i <- data.frame()
for(i in list.files.v){
  data.i.i <- read.delim(paste0("data/",i), header=FALSE)
  colnames(data.i.i) <- c("coord","chr","method","feature","start","end","score","strand","phase","attributes","te.id","piden")
  data.i.i$subset <- i
  data.i <- rbind(data.i,data.i.i)
}

data.i <- data.i[!grepl("unexp",data.i$te.id),]

table(data.i$subset)
data.i$subset[which(data.i$subset == "drGeuUrba1.1_stitched_strata_1.anno")] <- "strata1"
data.i$subset[which(data.i$subset == "drGeuUrba1.1_stitched_strata_2.anno")] <- "strata2"
data.i$subset[which(data.i$subset == "drGeuUrba1.1_stitched_strata_3.anno")] <- "strata3"

for(k in 1:7){
  col.i <- mapply(function (id,k) data.frame(str_split(id,"/",simplify = T))[,c(k)],data.i$te.id,k)
  data.i$col.i <- col.i
  colnames(data.i)[ncol(data.i)] <- letters[k]
}

fam <- mapply(function(d) ifelse(length(grep("Athila",d)) > 0, "ATHILA", "nonATHILA"),data.i$d)
data.i$fam <- fam

comb.cat <- mapply(function (cent.occ,fam) paste(cent.occ,fam,sep = "/"),data.i$c,data.i$fam)
table(comb.cat)
data.i$comb.cat <- comb.cat

data.i$comb.cat.1 <- mapply(function(x) ifelse(length(grep("OUT/",x)) > 0, "OUT/all", x), data.i$comb.cat)
table(data.i$comb.cat.1)

data.i$subset <- data.i$subset %>%
  {factor(.,levels = unique(data.i$subset)[c(3,1,2)])}

table(data.i$comb.cat.1,data.i$f)

# if only exp sites
data.i <- data.i[which(data.i$f == "exp"),]

# to plot only ATHILA/IN data
comb.cat.list <- unique(comb.cat)
comb.cat.list <- comb.cat.list[c(3)]
subset.v <- levels(data.i$subset)

rm(list = ls(pattern = "sub_plot_"))
rm(list = ls(pattern = "plot_"))
for(i in 1:length(comb.cat.list)){
  comb.cat.i <- comb.cat.list[i]
  
  for(k in 1:length(subset.v)){
    subset.i <- subset.v[k]
    
    stat.d <- data.i[which(data.i$comb.cat == comb.cat.i & data.i$subset == subset.i),] %>%
      group_by(subset) %>%
      summarise(n.counts = n(),
                median.prop = median(piden),
                comb.cat.i = comb.cat.i)
    
    p.1 <- ggplot(data.i[which(data.i$comb.cat == comb.cat.i & data.i$subset == subset.i),]) +
      geom_jitter(aes(x = piden, y = subset), color ="grey90",width = 0) +
      geom_boxplot(aes(x = piden, y = subset), fill = "transparent",outliers = F) +
      geom_text(data = stat.d,x=0.9515, y=stat.d$subset, label = paste0("n: ", stat.d$n.counts,"\nmedian: ",stat.d$median.prop),hjust = 0, size.unit = "pt", size = 6, vjust = .5) +
      coord_cartesian(xlim = c(0.95,1),
                      clip = "on",
                      expand = T) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.ticks.x = element_line(),
            axis.text.y = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = NA,colour = "black",linewidth = 0.5*2),
            panel.grid = element_blank(),
            panel.ontop = T,
            plot.title = element_text(hjust = .5))
    
    assign(paste0("sub_plot_",k),p.1)
  }
  
  p.2 <- cowplot::plot_grid(sub_plot_1,sub_plot_2,sub_plot_3,ncol = 1,align = "hv",axis = "l")
  assign(paste0("plot_",i),p.2)
}

p.all <- cowplot::plot_grid(plot_1,ncol = 1,align = "hv",axis = "l")

pdf("figure_4_panel_C.pdf",height = 2.75, width = (7.75*length(comb.cat.list))/4)
print(p.all)
dev.off()

