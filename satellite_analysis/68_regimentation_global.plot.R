library(ggplot2)


{
  bins_signal <- function(values, bins) {
    
    bins_width <- length(values) / bins
    bins_borders <- c(seq(1, length(values), bins_width), length(values) + 1)
    bins_borders <- round(bins_borders)
    
    unlist(sapply(1:bins, function(X) mean(values[bins_borders[X] : (bins_borders[X+1] - 1)])))
    
  }
  
  
}




setwd("/home/pwlodzimierz/ToL/upload_files/67_regimentation/data")

species <- list.dirs(path = ".", full.names = F, recursive = F)


chromosome_data_per_repeat <- data.frame(genome = vector(mode = "character"),
                                         chromosome = vector(mode = "character"),
                                         repeat_name = vector(mode = "character"),
                                         mean_HOR_score = vector(mode = "numeric"),
                                         total_regimented_windows = vector(mode = "numeric"),
                                         total_non_regimented_windows = vector(mode = "numeric"),
                                         whole_chromosome_regimentation_call = vector(mode = "numeric"),
                                         whole_chromosome_regimentation_p_value = vector(mode = "numeric"))


min_regimentation_bins_for_averaging <- 20
bins <- 20
species_HOR_signal <- list()
species_regimentation_averaged_signal <- list()

for(i in seq_along(species)) {
  cat(i, "\n")
  
  setwd(paste0("/home/pwlodzimierz/ToL/upload_files/67_regimentation/data/", species[i]))
  
  sats_names <- list.dirs(path = ".", full.names = F, recursive = F)
  
  chromosomes_HOR_signal <- list()
  chromosomes_regimentation_averaged_signal <- list()
  chr_signals_ID <- 1
  
  for(j in seq_along(sats_names)) {
    
    setwd(paste0("/home/pwlodzimierz/ToL/upload_files/67_regimentation/data/", species[i], "/", sats_names[j]))
    
    
    data_files <- list.files(path = ".", full.names = T)
    
    
    for(k in seq_along(data_files)) {
      load(data_files[k]) # data "a"
      
      # > str(a)
      # List of 8
      # $           :'data.frame':     1 obs. of  6 variables:
      # ..$ estimated_period: num 11
      # ..$ peak_power      : num 1.94e-07
      # ..$ spectral_entropy: num 0.848
      # ..$ peak_prominence : num 1.77
      # ..$ p_value         : num 0.00201
      # ..$ call            : chr "heterogeneous"
      # $           :'data.frame':     30 obs. of  7 variables:
      # ..$ window_start    : num [1:30] 1 56.2 111.5 166.8 222 ...
      # ..$ window_end      : num [1:30] 227 283 338 393 448 ...
      # ..$ estimated_period: num [1:30] 11.2 11.8 10.7 25 25 ...
      # ..$ peak_prominence : num [1:30] 1.582 1.434 0.837 4.09 4.237 ...
      # ..$ p_value         : num [1:30] 0.261 0.999 0.966 0.993 0.998 ...
      # ..$ call            : chr [1:30] "heterogeneous" "heterogeneous" "heterogeneous" "heterogeneous" ...
      # ..$ window_mid_bp   : int [1:30] 121781809 121791137 121800659 121809955 121819237 121828278 121837449 121846784 121861761 121870683 ...
      # $           :'data.frame':     2 obs. of  4 variables:
      # ..$ start           : num [1:2] 1 1192
      # ..$ end             : num [1:2] 1191 1811
      # ..$ call            : chr [1:2] "heterogeneous" "insufficient_data"
      # ..$ estimated_period: logi [1:2] NA NA
      # $           :'data.frame':     1 obs. of  3 variables:
      # ..$ period  : num 0
      # ..$ repeats : num 1811
      # ..$ fraction: num 100
      # $ chromosome: chr "CM000663.2"
      # $ genome    : chr "GCA_000001405"
      # $ sattelite : chr "170_2"
      # $ repeats   :'data.frame':     1811 obs. of  2 variables:
      # ..$ start    : int [1:1811] 103860766 103860797 103860967 103861138 103861310 103861480 103861651 103861820 103861991 103862157 ...
      # ..$ HOR_score: num [1:1811] 0 0 0 0 0 0 0 0 0 0 ...
      
      repeats <- a[["repeats"]]
      regimentation <- a[[2]]
      
      new_row <- data.frame(genome = a[["genome"]],
                            chromosome = a[["chromosome"]],
                            repeat_name = a[["sattelite"]],
                            mean_HOR_score = mean(repeats$HOR_score, na.rm = T),
                            total_regimented_windows = sum(regimentation$call == "periodic"),
                            total_non_regimented_windows = sum(regimentation$call != "periodic"),
                            whole_chromosome_regimentation_call = a[[1]]$call,
                            whole_chromosome_regimentation_p_value = a[[1]]$p_value)
      
      
      chromosome_data_per_repeat <- rbind(chromosome_data_per_repeat, new_row)
      
      if(nrow(regimentation) < min_regimentation_bins_for_averaging) next
      
      binned_HOR_score <- bins_signal(values = repeats$HOR_score, bins = bins)
      binned_regimentation <- bins_signal(values = ifelse(regimentation$call == "periodic", 1, 0), bins = bins)
      
      chromosomes_HOR_signal[[as.character(chr_signals_ID)]] <- binned_HOR_score
      chromosomes_regimentation_averaged_signal[[as.character(chr_signals_ID)]] <- binned_regimentation
      
      rev <- sample(c(0,1), 1)
      
      if(rev == 1) {
        chromosomes_HOR_signal[[as.character(chr_signals_ID)]] <- rev(binned_HOR_score)
        chromosomes_regimentation_averaged_signal[[as.character(chr_signals_ID)]] <- rev(binned_regimentation)
      }
      
      chr_signals_ID <- chr_signals_ID + 1
      
    }
    
  }
  
  if(length(chromosomes_HOR_signal) == 0) next
  rev <- sample(c(0,1), 1)
  if(rev == 0) {
    species_HOR_signal[[species[i]]] <- colMeans(do.call(rbind, chromosomes_HOR_signal), na.rm = T)
    species_regimentation_averaged_signal[[species[i]]] <- colMeans(do.call(rbind, chromosomes_regimentation_averaged_signal), na.rm = T)
    
  } else {
    species_HOR_signal[[species[i]]] <- rev(colMeans(do.call(rbind, chromosomes_HOR_signal), na.rm = T))
    species_regimentation_averaged_signal[[species[i]]] <- rev(colMeans(do.call(rbind, chromosomes_regimentation_averaged_signal), na.rm = T))
    
  }
  
  
  
}

setwd("/home/pwlodzimierz/ToL/upload_files/68_regimentation_summaries")

write.csv(x = chromosome_data_per_repeat, file = "chromosome_data_per_repeat.csv")






### per genome periodicity


# cen_families <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")

species_HOR_score <- vector(mode = "numeric", length = length(species))
species_regimentation_score <- vector(mode = "numeric", length = length(species))

chromosome_data_per_repeat$regimentation_score <- 100 * chromosome_data_per_repeat$total_regimented_windows / 
  (chromosome_data_per_repeat$total_regimented_windows + chromosome_data_per_repeat$total_non_regimented_windows)

for(i in seq_along(species)) {
  if(!(species[i] %in% chromosome_data_per_repeat$genome)) {
    print(species[i])
    next
  }
  
  cat(i, "\n")
  vals <- chromosome_data_per_repeat$mean_HOR_score[chromosome_data_per_repeat$genome == species[i]]
  if(length(vals) == 0) species_HOR_score[i] <- NA else species_HOR_score[i] <- mean(vals, na.rm = T)
  
  cat(i, "\n")
  vals <- chromosome_data_per_repeat$regimentation_score[chromosome_data_per_repeat$genome == species[i]]
  if(length(vals) == 0) species_regimentation_score[i] <- NA else species_regimentation_score[i] <- mean(vals, na.rm = T)
  
  
  
  
  
}


df <- data.frame(
  HOR = species_HOR_score,
  regimentation = species_regimentation_score
)

# Fit the model to extract R²
fit <- lm(regimentation ~ HOR, data = df)
r2 <- summary(fit)$r.squared
r2_label <- paste0("R^2 == ", round(r2, 3))

p <- ggplot(df, aes(x = HOR, y = regimentation)) +
  geom_point(shape = 16) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  annotate("text",
           x = 5, y = 98,
           label = r2_label,
           parse = TRUE,
           hjust = 0, vjust = 1,
           size = 5) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "HOR score", y = "Regimentation score") +
  theme_bw()

ggsave("species HOR vs regimentation score.pdf", plot = p, width = 5, height = 5)


### make df per repeat family
cen_families <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")

all_families <- TRUE
if(all_families) {
  cen_families$HOR_score <- 0
  cen_families$regimentation_score <- 0
  
  
  for(i in 1 : nrow(cen_families)) {
    
    sat_fasta_short <- strsplit(cen_families$fasta[i], split = "[.]")[[1]][1]
    sub_df <- chromosome_data_per_repeat[grep(sat_fasta_short, chromosome_data_per_repeat$genome),]
    sub_df <- sub_df[grep(cen_families$TRASH_new_class[i], sub_df$repeat_name),]
    
    cat(i, "\n")
    vals <- sub_df$mean_HOR_score
    if(length(vals) == 0) cen_families$HOR_score[i] <- NA else cen_families$HOR_score[i] <- mean(vals, na.rm = T)
    
    vals <- sub_df$regimentation_score
    if(length(vals) == 0) cen_families$regimentation_score[i] <- NA else cen_families$regimentation_score[i] <- mean(vals, na.rm = T)
    
    
  }
  
} else {
  
  cen_families$HOR_score <- NA
  cen_families$regimentation_score <- NA
}

### make a double plot


cen_families$fasta <- unlist(lapply(cen_families$fasta, function(X) strsplit(X, split = "[.]")[[1]][1]))

# cen_families <- rbind(cen_families, data.frame(X = 285,
#                                                fasta ="GCA_000001405",
#                                                TRASH_new_class = "170_2" ,
#                                                count_total = 0 ,
#                                                is_holocentric = FALSE ))

# --- Assign background colours per row -------------------------
cen_families$bg_col <- "#99999990"  # default
cen_families$bg_col[  1:154] <- "#3f37c990"
cen_families$bg_col[155:181] <- "#f7258590"
cen_families$bg_col[182:265] <- "#8ac92690"
cen_families$bg_col[263] <- "#99999990"

cen_families$order <- 0  # default
cen_families$order[  1:154] <- 2
cen_families$order[155:181] <- 3
cen_families$order[182:265] <- 4
cen_families$order[263] <- 1



if(!all_families) {
  remove_rows <- NULL
  for(i in 1 : nrow(cen_families)) {
    which_id <- grep(cen_families$fasta[i], species)
    if(length(which_id) != 1) {
      print(cen_families$fasta[i])
      remove_rows <- c(remove_rows, i) 
      next
    }
    cen_families$HOR_score[i] <- species_HOR_score[which_id]
    cen_families$regimentation_score[i] <- species_regimentation_score[which_id]
    if(i == nrow(cen_families)) break
    if(cen_families$fasta[i] == cen_families$fasta[i + 1]) remove_rows <- c(remove_rows, i)
  }
  
  
  cen_families <- cen_families[-remove_rows, ]
  
}

cen_families <- cen_families[order(cen_families$HOR_score, decreasing = FALSE),]
cen_families <- cen_families[order(cen_families$order, decreasing = FALSE),]

write.csv(cen_families, file = "cen_families.csv")

# modify the data frame for plotting, remove holocentrics and set to 0 the 5kb repeat
# that does not have any HORs
cen_families <- cen_families[cen_families$is_holocentric == FALSE,]
cen_families$HOR_score[is.na(cen_families$HOR_score)] <- 0
cen_families$regimentation_score[is.na(cen_families$regimentation_score)] <- 0

cen_families <- cen_families[order(cen_families$HOR_score, decreasing = FALSE),]
cen_families <- cen_families[order(cen_families$order, decreasing = FALSE),]

# cen_families <- cen_families[cen_families$HOR_score != 0,]
# --- Dynamic dimensions ----------------------------------------
n           <- nrow(cen_families)
plot_width  <- max(20, n * 0.25)
plot_height <- 17

x_idx <- seq_len(n)

# --- Data list for the 3 panels --------------------------------
panel_data <- list(
  cen_families$HOR_score,
  cen_families$regimentation_score
)

panel_ylab <- c(
  "HOR score (%)",
  "Regimentation score (%)"
)

# panel_col <- c("#4C8BE0", "#E07B4C", "#4CB87E")

# --- Label rows ------------------------------------------------
label_rows <- list(
  as.character(cen_families$fasta),
  as.character(cen_families$TRASH_new_class),
  as.character(cen_families$common_canonical_HOR_dist)
)

cex_lab    <- 1.6
row_lines  <- c(1.0, 14.0)
bottom_mar <- 30

# --- Equal plot-area heights via layout ------------------------
fig_lines  <- plot_height * 6
extra_frac <- bottom_mar / fig_lines
h_weights  <- c(1, 1 + extra_frac * 3)

# --- Helper: draw one panel ------------------------------------
draw_panel <- function(i, ymin = 0, ymax = 100) {
  is_bottom <- (i == 2)
  par(
    mar  = c(
      if (is_bottom) bottom_mar else 0,
      5,
      if (i == 1) 1 else 0,
      1
    ),
    xaxs = "i",
    yaxs = "i"
  )
  
  plot(NA,
       xlim    = c(0.5, n + 0.5), ylim = c(ymin, ymax),
       xlab    = "", ylab = panel_ylab[i],
       xaxt    = "n", yaxt = "s", las = 1, main = "",
       cex.lab = 3)
  
  # Grid lines
  abline(h = seq(20, 80, by = 20), col = "grey80", lty = "dotted")
  
  # Bars
  bar_w <- 0.45
  rect(x_idx - bar_w, 0, x_idx + bar_w, panel_data[[i]],
       col = cen_families$bg_col, border = NA)
  
  # box()
  
  if (is_bottom) {
    axis(1, at = x_idx, labels = FALSE, tick = TRUE)
    for (r in seq_len(length(label_rows))) {
      axis(1,
           at       = x_idx,
           labels   = label_rows[[r]],
           tick     = FALSE,
           las      = 2,
           cex.axis = cex_lab,
           line     = row_lines[r],
           col = cen_families$bg_col)
    }
  } else {
    axis(1, at = x_idx, labels = FALSE, tick = TRUE)
  }
}

# --- Shared plot call ------------------------------------------
do_plot <- function() {
  layout(matrix(1:2, nrow = 2, ncol = 1), heights = h_weights[1:2])
  for (i in 1:2) {
    
    plot_ymin <- 0
    plot_ymax <- 100
    
    draw_panel(i, plot_ymin, plot_ymax)
  }
}


# --- PDF -------------------------------------------------------
pdf("/home/pwlodzimierz/ToL/upload_files/68_regimentation_summaries/cen_families_2_panel.pdf",
    width = plot_width, height = plot_height)
do_plot()
dev.off()
message("Saved PDF  (", round(plot_width, 1), " x ", plot_height, " in)")





# --- Function to compute mean and 95% CI at each position ---
summarise_positions <- function(list_of_vectors) {
  mat <- do.call(rbind, list_of_vectors)
  mat[mat == 0] <- NA  # treat zeros as missing (drop this line if zeros are real values)
  
  means <- colMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 2, sd, na.rm = TRUE)
  ns    <- apply(mat, 2, function(x) sum(!is.na(x)))
  se    <- sds / sqrt(ns)
  
  # 95% CI using t-distribution (more accurate than assuming normal/z for small n)
  t_crit <- qt(0.975, df = ns - 1)
  ci_lower <- means - t_crit * se
  ci_upper <- means + t_crit * se
  
  data.frame(position = seq_along(means), mean = means,
             lower = ci_lower, upper = ci_upper)
}

# --- Compute summaries ---
hor_summary <- summarise_positions(species_HOR_signal)
regimentation_summary <- summarise_positions(species_regimentation_averaged_signal)

# --- Plotting function: colored line + shaded 95% CI ribbon ---
plot_with_ci <- function(df, line_col, ribbon_col, title, ylab) {
  plot(df$position, df$mean, type = "n",
       ylim = range(df$lower, df$upper, na.rm = TRUE),
       xlab = "Position along array", ylab = ylab,
       main = title, xaxt = "n")
  axis(1, at = df$position)
  
  # Shaded CI ribbon
  polygon(c(df$position, rev(df$position)),
          c(df$lower, rev(df$upper)),
          col = ribbon_col, border = NA)
  
  # Mean line + points
  lines(df$position, df$mean, col = line_col, lwd = 2)
  points(df$position, df$mean, col = line_col, pch = 19)
}

# --- Generate PDF ---
pdf("HOR score and regimentation along arrays.pdf")
par(mfrow = c(2, 1))

plot_with_ci(hor_summary,
             line_col = "steelblue", ribbon_col = adjustcolor("steelblue", alpha.f = 0.25),
             title = "Global HOR signal along array",
             ylab = "Mean HOR signal")

plot_with_ci(regimentation_summary,
             line_col = "darkorange", ribbon_col = adjustcolor("darkorange", alpha.f = 0.25),
             title = "Global regimentation signal along array",
             ylab = "Mean regimentation signal")

dev.off()




