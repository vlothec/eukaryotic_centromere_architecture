data <- read.csv("/home/pwlodzimierz/ToL/Metadata/Supplementary_Table_2_updated_v2_namefinal.csv")

# data <- read.csv(r"(C:\Users\Piotr Włodzimierz\Desktop\ToL\temp/Supplementary_Table_2_updated_v2_namefinal.csv)")

data$col <- ""

data$col[data$class == "invert"] <- "#3f37c9"
data$col[data$class == "chordate"] <- "#f72585"
data$col[data$class == "plant"] <- "#8ac926"

data$mean.gc <- 100 * data$mean.gc

fit <- lm(x ~ y, data = data.frame(x = data$Chromosome..GC, y = data$mean.gc))

setwd("/home/pwlodzimierz/ToL/upload_files/73_genome_vs_Statelite_GC")

fit_summary <-summary(fit)

pdf(file = "chr_gc_vs_sat_gc.pdf")
plot(data$Chromosome..GC, data$mean.gc, col = data$col, pch = 16, 
     main = paste0("R2 of x~y trendline (blue line) = ", round(fit_summary$r.squared, 2)))
lines(x=c(0,100), y = c(0,100), col = "red", lwd = 2)
abline(fit, col = "blue", lwd = 2)
dev.off()


pdf(file = "chr_gc_vs_sat_gc_extra_padding.pdf")
plot(data$Chromosome..GC, data$mean.gc, col = data$col, pch = 16, xlim = c(18,55), ylim = c(0,75),
     main = paste0("R2 of x~y trendline (blue line) = ", round(fit_summary$r.squared, 2)))
lines(x=c(0,100), y = c(0,100), col = "red", lwd = 2)
abline(fit, col = "blue", lwd = 2)

dev.off()

library(ggplot2)
ggplot(data, aes(Chromosome..GC, mean.gc)) +
  geom_point() +
  geom_smooth(method = "lm", se = T)
