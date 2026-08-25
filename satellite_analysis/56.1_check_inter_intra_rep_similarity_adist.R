
min_reps_per_chromosome <- 100
bins <- 30
sample_per_bin <- 10
sample_per_chromosome <- 100



repeat_files <- list.files(path = r"(/home/pwlodzimierz/ToL/upload_files/2026_filtered_TRASH)", pattern = "_repeats_filtered.csv", full.names = T)

satellite_metadata <- read.csv(file = r"(/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv)")

chromosome_metadata_full <- read.csv(file = r"(/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv)")



taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)#
print(i)

for(i in i) {
  
  print(i)
  
  assembly_l <- strsplit(basename(repeat_files[i]), split = "_")[[1]][1]
  assembly <- strsplit(assembly_l, split = ".fa")[[1]][1]
  
  satellite_metadata <- satellite_metadata[grep(assembly, satellite_metadata$fasta, ignore.case = T),]
  
  
  
  if(nrow(satellite_metadata) == 0) {
    print(paste0("Did not find cen satellite for ", assembly))
    next
  }
  
  
  chromosome_metadata <- chromosome_metadata_full[grep(assembly, chromosome_metadata_full$assembly.name),]
  
  chromosomes <- chromosome_metadata$chromosome.name[chromosome_metadata$is.chr == 1]
  
  cat(assembly, "chromosomes no:", length(chromosomes))
  
  repeats <- read.csv(file = repeat_files[i])
  
  centromeric_classes <- satellite_metadata$TRASH_new_class
  
  
  for(cen_class in centromeric_classes) {
    
    if(file.exists(paste0("/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/single_species_data/similarity_matrix_", assembly_l, "_c_", cen_class, ".rds"))) {
      next
    }
    if(file.exists(paste0("/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/single_species_data/similarity_matrix_", assembly, "_c_", cen_class, ".rds"))) {
      next
    }
    cat("\n", cen_class, "")
    repeats_cen_class <- repeats[repeats$new_class %in% cen_class,]
    
    if(mean(repeats_cen_class$width) > 500) {
      
      min_reps_per_chromosome <- round(min_reps_per_chromosome / 2)
      sample_per_bin <- round(sample_per_bin / 2)
      sample_per_chromosome <- round(sample_per_chromosome / 2)
    }
    
    
    similarity_matrix <- matrix(nrow = length(chromosomes), ncol = length(chromosomes))
    spacial_similarity_genome <- data.frame()
    
    for(j in seq_along(chromosomes)) {
      cat("\n", j, ": ")
      reps_1 <- repeats_cen_class[repeats_cen_class$seqID == chromosomes[j],]
      if(nrow(reps_1) < min_reps_per_chromosome) next
      
      reps_1$bin <- sort(rep(1:bins, length.out = nrow(reps_1)))
      
      seq_1_sample <- sample(reps_1$sequence, size = sample_per_chromosome, replace = T)
      
      
      for(k in j : length(chromosomes)) {
        cat(k, "")
        
        
        reps_2 <- repeats_cen_class[repeats_cen_class$seqID == chromosomes[k],]
        if(nrow(reps_2) < min_reps_per_chromosome) next
        
        reps_2$bin <- sort(rep(1:bins, length.out = nrow(reps_2)))
        
        seq_2_sample <- sample(reps_2$sequence, size = sample_per_chromosome, replace = T)
        
        similarity_matrix[j,k] <- 100 * (1 - mean(adist(seq_1_sample, seq_2_sample)) / mean(nchar(c(seq_1_sample, seq_2_sample))))
        
        spacial_chr_vector <- rep(0, bins)
        for(l in 1 : bins) {
          seq_bin1_sample <- sample(reps_1$sequence[reps_1$bin == l], size = sample_per_bin, replace = T)
          seq_bin2_sample <- sample(reps_2$sequence[reps_2$bin == l], size = sample_per_bin, replace = T)
          spacial_chr_vector[l] <- 100 * (1 - mean(adist(seq_bin1_sample, seq_bin2_sample)) / mean(nchar(c(seq_bin1_sample, seq_bin2_sample))))
        }
        spacial_similarity_genome <- rbind(spacial_similarity_genome, spacial_chr_vector)
      }
      
    }
    save(similarity_matrix, file = paste0("/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/single_species_data/similarity_matrix_", assembly, "_c_", cen_class, ".rds"))
    save(spacial_similarity_genome, file = paste0("/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/single_species_data/spacial_similarity_genome_", assembly, "_c_", cen_class, ".rds"))
    
  }
  remove(repeats, chromosomes, assembly)
  invisible(gc())
  cat(" done \n")
  
  
  
  
}

summarising <- FALSE

if(summarising) {
  
  min_reps_per_chromosome <- 100
  bins <- 30
  sample_per_bin <- 10
  sample_per_chromosome <- 100
  
  
  
  repeat_files <- list.files(path = r"(/home/pwlodzimierz/ToL/upload_files/2026_filtered_TRASH)", pattern = "_repeats_filtered.csv", full.names = T)
  
  satellite_metadata <- read.csv(file = r"(/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv)")
  
  chromosome_metadata_full <- read.csv(file = r"(/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv)")
  
  
  
  output_matrices <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/single_species_data/", pattern = "similarity_matrix", full.names = T)
  output_spacial <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/single_species_data/", pattern = "spacial_similarity_genome", full.names = T)
  
  
  
  
  similarity_matrices <- list()
  spacial_similarity <- list()
  
  for(i in output_matrices) {
    assembly <- basename(i)
    assemblyt <- strsplit(assembly, split = "similarity_matrix_")[[1]][2]
    assemblyt <- strsplit(assemblyt, split = ".rds")[[1]][1]
    assembly <- strsplit(assemblyt, split = "_c_")[[1]][1]
    cen_class <- strsplit(assemblyt, split = "_c_")[[1]][2]
    if(sum(grepl(assembly,satellite_metadata$fasta))) {
      load(i)
      similarity_matrices[[paste0(assembly, " ", cen_class)]] <- similarity_matrix
    } else {
      print(assembly)
    }
  }
  for(i in output_spacial) {
    assembly <- basename(i)
    assemblyt <- strsplit(assembly, split = "spacial_similarity_genome_")[[1]][2]
    assemblyt <- strsplit(assemblyt, split = ".rds")[[1]][1]
    assembly <- strsplit(assemblyt, split = "_c_")[[1]][1]
    cen_class <- strsplit(assemblyt, split = "_c_")[[1]][2]
    if(sum(grepl(assembly,satellite_metadata$fasta))) {
      load(i)
      spacial_similarity[[paste0(assembly, " ", cen_class)]] <- colMeans(spacial_similarity_genome)
    } else {
      print(assembly)
    }
  }
  
  # check which repeat families are not done
  found <- names(similarity_matrices)
  all <- unlist(sapply(1:nrow(satellite_metadata), function(X) paste(satellite_metadata$fasta[X], satellite_metadata$TRASH_new_class[X], sep = " ")))
  all <- as.character(unlist(sapply(all, function(X) paste(strsplit(X, split = ".fa")[[1]][1], strsplit(X, split = ".fa")[[1]][2], sep = ""))))
  all[!(all %in% found)]
  
  
  draw_track <- function(values, y_pos, col_rgb, max_val = 0, track_h = 0.9) {
    if(max_val == 0) max_val <- max(values, na.rm = TRUE)
    if(max_val <= 0) return(invisible(NULL))
    alphas <- values / max_val
    nonzero <- 1:bins
    for(k in nonzero) {
      col_with_alpha <- rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alphas[k])
      rect(xleft  = k-1,
           xright = k,
           ybottom = y_pos - track_h/2,
           ytop    = y_pos + track_h/2,
           col = col_with_alpha, border = NA)
    }
  }
  
  
  spacial_similarity_df <- as.data.frame(do.call(rbind, spacial_similarity))
  spacial_similarity_df <- spacial_similarity_df[!is.na(spacial_similarity_df[,1]),]
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/spacial_adist_all.pdf)", width = 10, height = 10)
  plot(NULL, NULL, xlim = c(0,bins), ylim = c(-50, nrow(spacial_similarity_df)),
       xlab = "centromeric bins",
       ylab = "genome index, at -30 mean normalised plotted coverage")
  for(i in 1 : nrow(spacial_similarity_df)) {
    draw_track(spacial_similarity_df[i,] - min(spacial_similarity_df[i,]), i, c(0, 0, 0))
  }
  draw_track(colMeans(spacial_similarity_df) - min(colMeans(spacial_similarity_df)), -30, c(0, 0, 0), track_h = 10)
  dev.off()
  
  
  
  ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/spacial_adist_all_lineplot.pdf)", width = 10, height = 10)
  plot(NULL, NULL, xlim = c(0,bins), ylim = c(min(spacial_similarity_df), max(spacial_similarity_df)),
       xlab = "centromeric bins",
       ylab = "similarity between bins")
  for(i in 1 : nrow(spacial_similarity_df)) {
    vals <- spacial_similarity_df[i,]
    vals <- ma(c(vals[1],vals[1],vals,vals[length(vals)],vals[length(vals)]))[3:(length(vals)+2)]
    lines(1:bins, vals)
  }
  
  lines(1:bins, colMeans(spacial_similarity_df), lwd = 2, col = "blue")
  dev.off()
  
  
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/spacial_adist_all_scatterplot.pdf)", width = 10, height = 10)
  plot(NULL, NULL, xlim = c(1,ceiling(bins/2)), ylim = c(-10, 25),
       xlab = "edge to middle of the chromosome in bins",
       ylab = "similarity between bins")
  for(i in 1 : nrow(spacial_similarity_df)) {
    vals <- as.numeric(((spacial_similarity_df[i,] + rev(spacial_similarity_df[i,])) / 2)[1:ceiling(bins/2)])
    vals <- ma(c(vals[1],vals[1],vals,vals[length(vals)],vals[length(vals)]))[3:(length(vals)+2)]
    vals <- vals - vals[1]
    lines(1:ceiling(bins/2), vals)
  }
  
  vals <- ((colMeans(spacial_similarity_df) + rev(colMeans(spacial_similarity_df))) / 2)[1:ceiling(bins/2)]
  vals <- ma(c(vals[1],vals[1],vals,vals[length(vals)],vals[length(vals)]))[3:(length(vals)+2)]
  vals <- vals - vals[1]
  lines(1:ceiling(bins/2), vals, lwd = 2, col = "blue")
  dev.off()
  
  
  
  
  
  
  
  
  
  inter_intra_values <- data.frame()
  
  for(i in 1 : length(similarity_matrices)) {
    inter_vals <- similarity_matrices[[i]][upper.tri(similarity_matrices[[i]])]
    intra_vals <- diag(similarity_matrices[[i]])
    inter_intra_values <- rbind(inter_intra_values, data.frame(mean(inter_vals, na.rm = T), mean(intra_vals, na.rm = T)))
    rownames(inter_intra_values)[i] <- names(similarity_matrices)[i]
    
  }
  
  colnames(inter_intra_values) <- c("inter_sim_perc", "intra_sim_perc")
  
  phyla <- read.csv("/home/pwlodzimierz/ToL/Metadata/fasta_phyla.csv")
  
  inter_intra_values$clade <- ""
  for(i in 1 : nrow(inter_intra_values)) {
    assembly <- rownames(inter_intra_values)[i]
    assembly <- strsplit(assembly, split = " ")[[1]][1]
    assembly <- strsplit(assembly, split = "[.]")[[1]][1]
    inter_intra_values$clade[i] <- phyla$phyla[grep(assembly, phyla$fasta)]
  }
  
  
  
  
  
  fast_min_variance <- function(vec) {
    s <- sort(vec)
    gaps <- diff(s)
    split_at <- which.max(gaps)
    g1 <- s[1:split_at]
    g2 <- s[(split_at+1):length(s)]
    mean(c(ifelse(length(g1) < 2, 0, var(g1)),
           ifelse(length(g2) < 2, 0, var(g2))))
  }
  
  inter_intra_values$variance <- NA
  inter_intra_values$sd <- NA
  inter_intra_values$variance_normalised_by_mean_value <- NA
  inter_intra_values$mean_of_two_minimised_variances <- NA
  inter_intra_values$mean_of_two_minimised_variances_norm_by_mean_val <- NA
  inter_intra_values$mean_of_two_var_norm_normalised_by_var_norm <- NA
  inter_intra_values$diff_of_two_min_var_and_var <- NA
  for(i in 1 : length(similarity_matrices)) {
    similarity_matrices[[i]][similarity_matrices[[i]] < 0] = 0
    inter_vals <- similarity_matrices[[i]][upper.tri(similarity_matrices[[i]])]
    intra_vals <- diag(similarity_matrices[[i]])
    inter_vals <- inter_vals[!is.na(inter_vals)]
    if(length(inter_vals) < 2) next
    
    
    inter_intra_values$variance[i] <- var(inter_vals[inter_vals != 0], na.rm = T)
    inter_intra_values$sd[i] <- sd(inter_vals[inter_vals != 0], na.rm = T)
    inter_intra_values$variance_normalised_by_mean_value[i] <- inter_intra_values$variance[i] / inter_intra_values$inter_sim_perc[i]
    inter_intra_values$mean_of_two_minimised_variances[i] <- fast_min_variance(inter_vals)
    inter_intra_values$mean_of_two_minimised_variances_norm_by_mean_val[i] <- inter_intra_values$mean_of_two_minimised_variances[i] / inter_intra_values$inter_sim_perc[i]
    inter_intra_values$mean_of_two_var_norm_normalised_by_var_norm[i] <- inter_intra_values$mean_of_two_minimised_variances_norm_by_mean_val[i] / inter_intra_values$variance_normalised_by_mean_value[i]
    inter_intra_values$diff_of_two_min_var_and_var[i] <- inter_intra_values$variance[i] / inter_intra_values$mean_of_two_minimised_variances[i]
    
  }
  inter_intra_values$id <- 1 : nrow(inter_intra_values)
  
  inter_intra_values_t <- inter_intra_values
  inter_intra_values_t <- inter_intra_values_t[!is.na(inter_intra_values_t$variance),]
  
  
  
  inter_intra_values$col <- "black"
  inter_intra_values$col[inter_intra_values$clade == "chordate"] <- "#f72585cc"
  inter_intra_values$col[inter_intra_values$clade == "planta"] <- "#8ac926cc"
  inter_intra_values$col[inter_intra_values$clade == "invertebrate"] <- "#3f37c9cc"
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/inter_intra_similarity.pdf)", width = 10, height = 10)
  
  plot(inter_intra_values$inter_sim_perc, inter_intra_values$intra_sim_perc, col = inter_intra_values$col, pch = 20, cex = 3,
       xlab = "Mean interchromosomal similarity", ylab = "Mean intrachromosomal similarity", 
       xlim = c(20,100), ylim = c(20,100))
  

  lines(x = 0:900, y = 0:900, lty = 2, lwd = 2)
  dev.off()
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/sd_vs_var.pdf)", width = 10, height = 10)
  plot(inter_intra_values$variance, inter_intra_values$sd)
  dev.off()
  
  # scatter col by variance
  inter_intra_values$col <- colorRampPalette(c("#0099eebb", "#ee9900bb"))(length(inter_intra_values$sd))[rank(inter_intra_values$sd)]
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/inter_intra_similarity_col_by_sd.pdf)", width = 10, height = 10)
  
  plot(x = c(-100,100), y = c(-100,100), col = "red", type = "l", lwd = 2, lty = 1,
       xlab = "Mean interchromosomal similarity", ylab = "Mean intrachromosomal similarity", 
       xlim = c(0,100), ylim = c(0,100))
  points(inter_intra_values$inter_sim_perc, inter_intra_values$intra_sim_perc, col = inter_intra_values$col, pch = 20, cex = 3)
  
  {
    # Legend: colour ramp for variance
    legend_cols <- colorRampPalette(c("#0099eebb", "#ee9900bb"))(100)
    legend_x <- grconvertX(0.95, from = "npc", to = "user")
    legend_y <- grconvertY(0.05, from = "npc", to = "user")
    bar_width  <- diff(grconvertX(c(0, 0.03), "npc", "user"))
    bar_height <- diff(grconvertY(c(0, 0.25), "npc", "user"))
    n <- length(legend_cols)
    step <- bar_height / n
    
    for (i in seq_len(n)) {
      rect(xleft   = legend_x - bar_width,
           xright  = legend_x,
           ybottom = legend_y + (i - 1) * step,
           ytop    = legend_y + i * step,
           col     = legend_cols[i], border = NA)
    }
    
    rect(legend_x - bar_width, legend_y,
         legend_x, legend_y + bar_height,
         col = NA, border = "black", lwd = 0.8)
    
    var_range <- range(inter_intra_values$sd[!is.na(inter_intra_values$sd)])
    text(legend_x - bar_width / 2, legend_y - diff(grconvertY(c(0, 0.015), "npc", "user")),
         labels = formatC(var_range[1], format = "g", digits = 3),
         adj = c(0.5, 1), cex = 0.75)
    text(legend_x - bar_width / 2, legend_y + bar_height + diff(grconvertY(c(0, 0.015), "npc", "user")),
         labels = formatC(var_range[2], format = "g", digits = 3),
         adj = c(0.5, 0), cex = 0.75)
    text(2+legend_x - bar_width / 2, legend_y + bar_height / 2,
         labels = "SD, % identity", srt = 90, adj = c(0.5, 1.8), cex = 0.8)
    
  }
  dev.off()
  
  cor.test(inter_intra_values$sd, (inter_intra_values$inter_sim_perc - inter_intra_values$intra_sim_perc))
  
  
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/inter_intra_similarity_mixed_cols.pdf)", width = 10, height = 10)
  
  inter_intra_values$col <- colorRampPalette(c("black", "white"))(length(inter_intra_values$variance))[rank(inter_intra_values$variance)]
  plot(inter_intra_values$inter_sim_perc, inter_intra_values$intra_sim_perc, col = inter_intra_values$col, pch = 19, cex = 2,
       xlab = "Mean interchromosomal similarity", ylab = "Mean intrachromosomal similarity", 
       xlim = c(20,100), ylim = c(20,100))
  
  inter_intra_values$col[inter_intra_values$clade == "chordate"] <- "#f72585cc"
  inter_intra_values$col[inter_intra_values$clade == "planta"] <- "#8ac926cc"
  inter_intra_values$col[inter_intra_values$clade == "invertebrate"] <- "#3f37c9cc"
  points(inter_intra_values$inter_sim_perc, inter_intra_values$intra_sim_perc, col = inter_intra_values$col, pch = 1, cex = 2)
  points(inter_intra_values$inter_sim_perc, inter_intra_values$intra_sim_perc, col = inter_intra_values$col, pch = 1.1, cex = 2)
  
  
  lines(x = 0:900, y = 0:900, lty = 2, lwd = 2)
  dev.off()
  
  
  
  library(plot.matrix)
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/matrices_sorted_by_variance.pdf)", 
      width = 10, height = 10, onefile = T)
  inter_intra_values_t <- inter_intra_values_t[order(inter_intra_values_t$variance),]
  for(i in 1 : nrow(inter_intra_values_t)) {
    m_clamped <- pmin(pmax(similarity_matrices[[inter_intra_values_t$id[i]]], 20), 100)
    rownames(m_clamped) <- paste0("Chr ", 1:nrow(m_clamped))
    colnames(m_clamped) <- paste0("Chr ", 1:ncol(m_clamped))
    pastel_pal <- colorRampPalette(c("#ffffcc", "#a1dab4", "#41b6c4", "#c994c7"))
    
    keep <- which(!is.na(diag(m_clamped)))
    m_filtered <- m_clamped[keep, keep]
    
    par(mar = c(5.1, 5.1, 4.1, 5.1))
    plot(m_filtered,
         col    = pastel_pal,
         breaks = seq(20, 100, length.out = 11),  # 10 colour bins between 20–100
         na.col = "grey80",
         key    = list(side = 4, las = 1),
         xlab   = "", ylab   = "",
         main   = rownames(inter_intra_values_t)[i])
    
  }
  dev.off()
  
  ### internal similarity vs variance
  
  pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/log_variance_vs_intra.pdf)", 
      width = 10, height = 10, onefile = T)
  plot(log(inter_intra_values_t$variance), log(abs(inter_intra_values_t$intra_sim_perc)),
       xlab = "log variance in interchromosomal similarity",
       ylab = "log mean intrachromosomal sim",
       main = "r2 = 0.347, p=10e-14")
  dev.off()
  
  library(ggplot2)
  
  x <- log(inter_intra_values_t$variance)
  y <- log(abs(inter_intra_values_t$intra_sim_perc))
  
  # Pearson vs Spearman
  pearson  <- cor.test(x, y, method = "pearson")
  spearman <- cor.test(x, y, method = "spearman")
  
  cat("Pearson  r =", round(pearson$estimate, 3),
      " | r2 =", round(pearson$estimate^2, 3),
      " | p =", format(pearson$p.value, scientific = TRUE), "\n")
  cat("Spearman rho =", round(spearman$estimate, 3),
      " | p =", format(spearman$p.value, scientific = TRUE), "\n")
  
  # Identify leverage of upper-right outliers
  df <- data.frame(x = x, y = y)
  x_thresh <- quantile(x, 0.90)
  y_thresh <- quantile(y, 0.90)
  
  pearson_no_outliers <- cor.test(df$x[df$x < x_thresh & df$y < y_thresh],
                                  df$y[df$x < x_thresh & df$y < y_thresh],
                                  method = "pearson")
  
  cat("\nPearson r (bottom 90% of both axes) =",
      round(pearson_no_outliers$estimate, 3),
      " | r2 =", round(pearson_no_outliers$estimate^2, 3),
      " | p =", format(pearson_no_outliers$p.value, scientific = TRUE), "\n")
  
  # Plot with outlier labels
  df$outlier <- T
  df$label <- ifelse(df$outlier, rownames(df), NA)
  
  ggplot(df, aes(x, y)) +
    geom_point(aes(color = outlier), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "#D85A30", linewidth = 0.8) +
    geom_smooth(data = subset(df, !outlier), method = "lm", se = FALSE,
                color = "#1D9E75", linewidth = 0.8, linetype = "dashed") +
    geom_text(aes(label = label), size = 2.8, vjust = -0.8, na.rm = TRUE) +
    scale_color_manual(values = c("FALSE" = "#534AB7", "TRUE" = "#D85A30"),
                       labels = c("core", "top-10% outlier")) +
    labs(x = "variance in interchromosomal similarity",
         y = "mean intra similarity",
         color = NULL,
         caption = paste0("Pearson r² = ", round(pearson$estimate^2, 3),
                          " | Spearman rho = ", round(spearman$estimate, 3))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  ggsave(filename = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/variance_vs_intra_pearson_vs_spearman.pdf)")
  
  {
    
    
    x <- log(inter_intra_values_t$sd)
    y <- log(abs(inter_intra_values_t$intra_sim_perc - inter_intra_values_t$inter_sim_perc))
    
    # Pearson vs Spearman
    pearson  <- cor.test(x, y, method = "pearson")
    spearman <- cor.test(x, y, method = "spearman")
    
    cat("Pearson  r =", round(pearson$estimate, 3),
        " | r2 =", round(pearson$estimate^2, 3),
        " | p =", format(pearson$p.value, scientific = TRUE), "\n")
    cat("Spearman rho =", round(spearman$estimate, 3),
        " | p =", format(spearman$p.value, scientific = TRUE), "\n")
    
    
    ### log delta vs sd
    pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/log_delta_vs_sd.pdf)", 
        width = 10, height = 10, onefile = T)
    plot(log(inter_intra_values_t$sd), log(abs(inter_intra_values_t$intra_sim_perc - inter_intra_values_t$inter_sim_perc)),
         xlab = "log variance in interchromosomal similarity",
         ylab = "log mean delta sim",
         main = paste0("Spearman rho =", round(spearman$estimate, 3),
                       " | p =", format(spearman$p.value, scientific = TRUE)))
    dev.off()
    
    
    
  }
  
  {
    
    
    x <- (inter_intra_values_t$sd)
    y <- (abs(inter_intra_values_t$intra_sim_perc - inter_intra_values_t$inter_sim_perc))
    
    # Pearson vs Spearman
    pearson  <- cor.test(x, y, method = "pearson")
    spearman <- cor.test(x, y, method = "spearman")
    
    cat("Pearson  r =", round(pearson$estimate, 3),
        " | r2 =", round(pearson$estimate^2, 3),
        " | p =", format(pearson$p.value, scientific = TRUE), "\n")
    cat("Spearman rho =", round(spearman$estimate, 3),
        " | p =", format(spearman$p.value, scientific = TRUE), "\n")
    ###  delta vs sd
    pdf(file = r"(/home/pwlodzimierz/ToL/upload_files/56.1_inter_intra_similarity/delta_vs_sd.pdf)", 
        width = 10, height = 10, onefile = T)
    plot((inter_intra_values_t$sd), (abs(inter_intra_values_t$intra_sim_perc - inter_intra_values_t$inter_sim_perc)),
         xlab = "variance in interchromosomal similarity",
         ylab = "mean delta sim",
         main = paste0("Spearman rho =", round(spearman$estimate, 3),
                       " | p =", format(spearman$p.value, scientific = TRUE)))
    dev.off()
    
    
    
    
  }
}

