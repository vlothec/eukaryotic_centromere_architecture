## 0. Prepare the data and choose which species to sue
## 1. Make histograms of distance between repeats sizes to find out the best array size definition
## 2. Calculate internal similarity measures within arrays, between arrays same chr, between arrays 
##    different chr: boxplot for each species
## 3. Scatter plots of repeat content vs chromosome size and array number vs chromosome size
# 4. HOR scores averaged within arrays
## 5. Strand switching statistics

array_pairs_to_sample_for_similarity <- 1500
repeats_to_sample_per_array <- 10
bin_size_simple_plot <- 50000


library(GenomicRanges)
library(seqinr)
library(ggplot2)

chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

species_to_analyse <- c("lpCarDepa1.1", "lpSchLacu1", "lpLuzSylv1.1", "iiLimLuna2.1", "iiLimMarm1.1",
                        "iiLimRhom1.1", 
                        "lpEleQuin1.1",
                        "lpCypFusc1.hap1.1", "lpSchTriq1.1", "lpSchTabe1.hap1.1", "lpEriAngu1.1", "lpSciSylv1.hap1.1",
                        "lpLuzPall1.1")


repeats_to_analyse <- list("81_2", c("183_4", "58_10"), c("124_1", "174_2"), "353_2", '166_3',
                           "161_5",
                           "80_1",
                           "325_1", c("183_2",	"58_1"), c("183_3",	"58_1"), "31_5", "31_4",
                           c("125_3",	"104_2"))

# species_to_analyse <- c("lpCarDepa1.1", "lpSchLacu1", "lpLuzSylv1.1", "iiLimLuna2.1", "iiLimMarm1.1",
#                         "iiLimRhom1.1", "ioIscEleg1.1", "ihAelAcum1.1", "ihIcePurc2.1", "iuLoeVari1.hap1.1",
#                         "iuPsoGibb1.1", "igLabMino1.1", "iiLimAuri1.1", "ilPieNapi4.1", "lpEleQuin1.1",
#                         "lpCypFusc1.hap1.1", "lpSchTriq1.1", "lpSchTabe1.hap1.1", "lpEriAngu1.1", "lpSciSylv1.hap1.1",
#                         "lpLuzPall1.1")
# 
# 
# repeats_to_analyse <- list("81_2", c("183_4", "58_10"), c("124_1", "174_2"), "353_2", '166_3',
#                            "161_5", "315_60", c("154_72",	"155_1"), "287_11", c("140_5", "549_9",	"146_8"),
#                            "166_1", c("161_1",	"160_3"), "169_1", c("667_2",	"79_1"), "80_1",
#                            "325_1", c("183_2",	"58_1"), c("183_3",	"58_1"), "31_5", "31_4",
#                            c("125_3",	"104_2"))

# species_to_analyse <- c("lpSchLacu1", "lpSchLacu1", "lpLuzSylv1.1", "lpLuzSylv1.1", 
#                         "ihAelAcum1.1", "ihAelAcum1.1", "iuLoeVari1.hap1.1", "iuLoeVari1.hap1.1", "iuLoeVari1.hap1.1",
#                         "igLabMino1.1", "igLabMino1.1", "ilPieNapi4.1", "ilPieNapi4.1", 
#                         "lpSchTriq1.1", "lpSchTriq1.1", "lpSchTabe1.hap1.1", "lpSchTabe1.hap1.1", 
#                         "lpLuzPall1.1", "lpLuzPall1.1")
# 
# 
# repeats_to_analyse <- list("183_4", "58_10", "124_1", "174_2", 
#                            "154_72",	"155_1", "140_5", "549_9",	"146_8", 
#                            "161_1",	"160_3", "667_2",	"79_1",
#                            "183_2",	"58_1", "183_3",	"58_1", 
#                            "125_3",	"104_2")


min_distances <- rep(50000,30)

if(F) {
  for(i in seq_along(species_to_analyse)) {
    cat(i, "\n")
    setwd(list.files("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", pattern = species_to_analyse[i], full.names = T))
    
    repeats <- read.csv(list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE, recursive = FALSE))
    
    repeats <- repeats[repeats$new_class %in% repeats_to_analyse[[i]], ]
    
    repeats <- repeats[order(repeats$start),]
    repeats <- repeats[order(repeats$seqID),]
    
    dist_to_next <- repeats$start[2:nrow(repeats)] - repeats$end[1 : (nrow(repeats) - 1)]
    dist_to_next <- dist_to_next[dist_to_next > 5000]
    dist_to_next <- dist_to_next[dist_to_next < 500000]
    
    pdf(paste0(species_to_analyse[i], "_distances_between_repeats_histogram.pdf"), 
        width = 12)
    hist(dist_to_next, breaks = seq(0,500000, by = 5000), xlim = c(0,500000))
    abline(v = min_distances[i], col = "red")
    dev.off()
    
    pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_histograms_holocentric_distances_between_repeats/", species_to_analyse[i], "_distances_between_repeats_histogram.pdf"), 
        width = 12)
    hist(dist_to_next, breaks = seq(0,500000, by = 5000), xlim = c(0,500000))
    abline(v = min_distances[i], col = "red")
    dev.off()
    
  }
}



colors <- c(
  "#E63946", "#FF9B9B",  # red pair
  "#2196F3", "#90CAF9",  # blue pair
  "#2E7D32", "#A5D6A7",  # green pair
  "#F4A000", "#FFD980"   # amber pair
)

for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  setwd(list.files("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", pattern = species_to_analyse[i], full.names = T))
  
  chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv")
  chr_sizes <- chr_sizes[grep(species_to_analyse[i], chr_sizes$assembly.name), ]
  chr_sizes <- chr_sizes[chr_sizes$is.chr == 1, ]
  
  chr_sizes <- chr_sizes[,-1]
  chr_sizes <- unique(chr_sizes)
  
  
  repeats_all <- read.csv(list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE, recursive = FALSE))
  repeats_all <- repeats_all[repeats_all$new_class %in% repeats_to_analyse[[i]], ]
  
  ### simple plot
  if (F){ 
    pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/33_simple_repeat_plot/", species_to_analyse[i], ".pdf"), width = 8, height = (nrow(chr_sizes)/2 + 1)) # 1 inch for top and bottom margins and 1 per plot
    
    line_height <- par("cin")[2]  
    outer_top_in    <- 0.3  # inches
    outer_bottom_in <- 0.7  # inches
    par(
      mfrow = c(nrow(chr_sizes), 1),      # chr_sizes rows, 1 column
      mar  = c(0, 4, 0, 2),  
      oma = c(outer_bottom_in / line_height,   # bottom
              0,
              outer_top_in    / line_height,   # top
              0)
    )
    
    for(z in 1 : nrow(chr_sizes)) {
      
      plot(x = NA, y = NA, xlim = c(0,chr_sizes$size[z]), ylim = c(0,100), 
           xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "h",
           col = "#59a46c")
      
      bins <- tileGenome(setNames(chr_sizes$size[z], chr_sizes$chromosome.name[z]),
                         tilewidth = bin_size_simple_plot, cut.last.tile.in.chrom = T)
      repeats_chr <- repeats_all[repeats_all$seqID == chr_sizes$chromosome.name[z],]
      
      for(zz in 1 : length(repeats_to_analyse[[i]])) {
        repeats_class <- repeats_chr[repeats_chr$new_class %in% repeats_to_analyse[[i]][zz],]
        
        for(s in 1:2) {
          repeats_s <- repeats_class[repeats_class$strand == c("+","-")[s],]
          if(nrow(repeats_s) == 0) next
          
          rep_gr <- GRanges(repeats_s$seqID, IRanges(repeats_s$start, repeats_s$end))
          
          hits      <- findOverlaps(bins, rep_gr)
          overlap_w <- width(pintersect(bins[queryHits(hits)], rep_gr[subjectHits(hits)]))
          
          agg <- tapply(overlap_w, queryHits(hits), sum)
          bins$coverage <- 0
          bins$coverage[as.integer(names(agg))] <- agg
          bins$coverage <- 100 * bins$coverage / width(bins)
          
          
          points(x = mid(bins)[bins$coverage != 0], y = bins$coverage[bins$coverage != 0], 
                 type = "h",
                 col = colors[zz+s-1])
        }
        
      }
      
      axis(side = 1, at = seq(0,chr_sizes$size[z], by = 1000000), labels = F, tick = T)
    }
    labels <- NULL
    for(z in seq_along(repeats_to_analyse[[i]])) {
      labels <- c(labels, paste0(repeats_to_analyse[[i]][z], "(+)"), paste0(repeats_to_analyse[[i]][z], "(-)"))
    }
    
    par(mfg = c(1, 1))
    
    
    usr <- par("usr")
    x <- usr[1]
    y <- usr[4]  # top of plot area
    
    for (j in seq_along(labels)) {
      text(x, y, labels[j], col = colors[j], cex = 1,
           adj = c(0, -0.5),   # adj[2] pushes above the top border
           xpd = NA)            # allows drawing outside plot region
      # advance x by width of previous label
      x <- x + strwidth(labels[j], cex = 1) * 1.15  # 1.15 adds small gap
    }
    
    dev.off()
    
  }
  
  ### array definitions
  for(z in 1 : length(repeats_to_analyse[[i]])) {
    
    repeats <- repeats_all[repeats_all$new_class %in% repeats_to_analyse[[i]][z], ]
    
    repeats <- repeats[order(repeats$start),]
    repeats <- repeats[order(repeats$seqID),]
    
    repeats <- repeats[repeats$seqID %in% chr_sizes$chromosome.name, ]
    
    dist_to_next <- repeats$start[2:nrow(repeats)] - repeats$end[1 : (nrow(repeats) - 1)]
    
    array_ends <- which(dist_to_next > min_distances[i])
    
    array_ends <- c(array_ends, which(dist_to_next < 0))
    
    array_ends <- c(array_ends, length(dist_to_next))
    
    repeats$end_array <- FALSE
    repeats$end_array[array_ends] <- TRUE
    repeats$end_array[nrow(repeats)] <- TRUE
    
    repeats$start_array <- FALSE
    repeats$start_array[array_ends + 1] <- TRUE
    repeats$start_array[1] <- TRUE
    
    repeats$width <- repeats$end - repeats$start + 1
    
    for(j in unique(repeats$seqID)) {
      repeats$start_array[repeats$seqID == j][1] = TRUE
      repeats$end_array[repeats$seqID == j][nrow(repeats[repeats$seqID == j,])] = TRUE
    }
    
    array_starts <- which(repeats$start_array)
    array_ends <- which(repeats$end_array)
    
    repeats$array_ID <- NA
    
    arrays <- NULL
    repeats$holo_array_ID <- 0
    for(j in seq_along(array_starts)) {
      array_repeats <- repeats[array_starts[j] : array_ends[j], ]
      repeats$holo_array_ID[array_starts[j] : array_ends[j]] = j
      arrays <- rbind(arrays, data.frame(holo_array_ID = j,
                                         chromosome = names(sort(table(array_repeats$seqID), decreasing = TRUE))[1], 
                                         start = min(array_repeats$start), 
                                         end = max(array_repeats$end), 
                                         rep_number = nrow(array_repeats), 
                                         length = max(array_repeats$end) - min(array_repeats$start), 
                                         repeat_total_length = sum(array_repeats$width), 
                                         plust_strand_perc = length(which(array_repeats$strand == "+")) / nrow(array_repeats),
                                         chr_size = chr_sizes$size[chr_sizes$chromosome.name == names(sort(table(array_repeats$seqID), decreasing = TRUE))[1]]))
    }
    
    
    arrays <- arrays[arrays$rep_number >= 2,]
    
    write.csv(arrays, file = paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"), row.names = FALSE)
    
    write.csv(arrays, file = paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"), row.names = FALSE)
    
    write.csv(repeats, file = paste0("/home/pwlodzimierz/ToL/upload_files/33_repeats_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_repeats.csv"), row.names = FALSE)
    
    print("### repeat content vs chromosome size and array number vs chromosome size ###")
    
    if(F) {
      ### repeat content vs chromosome size and array number vs chromosome size ###
      {
      chr_sizes$repeat_array_number <- NA
      chr_sizes$repeat_bp <- NA
      
      for(j in 1 : nrow(chr_sizes)) {
        chr_arrays <- arrays[arrays$chromosome == chr_sizes$chromosome.name[j], ]
        chr_sizes$repeat_array_number[j] = nrow(chr_arrays)
        chr_sizes$repeat_bp[j] = sum(chr_arrays$repeat_total_length)
      }
      
      pdf(paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_array_no_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_array_number, chr_sizes$size,
           xlab = "Repeat array number", ylab = "Chromosome size",
           main = "Repeat array number vs Chromosome size")
      
      # Add trendline
      model1 <- lm(size ~ repeat_array_number, data = chr_sizes)
      abline(model1, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result1 <- cor.test(chr_sizes$repeat_array_number, chr_sizes$size)
      r1 <- round(cor_result1$estimate, 2)
      p1 <- cor_result1$p.value
      sig1 <- if (p1 < 0.001) "***" else if (p1 < 0.01) "**" else if (p1 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r1, ", p = ", signif(p1, 2), " (", sig1, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      pdf(paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_repeat_bp_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_bp, chr_sizes$size,
           xlab = "Repeat base pairs", ylab = "Chromosome size",
           main = "Repeat base pairs vs Chromosome size")
      
      # Add trendline
      model2 <- lm(size ~ repeat_bp, data = chr_sizes)
      abline(model2, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result2 <- cor.test(chr_sizes$repeat_bp, chr_sizes$size)
      r2 <- round(cor_result2$estimate, 2)
      p2 <- cor_result2$p.value
      sig2 <- if (p2 < 0.001) "***" else if (p2 < 0.01) "**" else if (p2 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r2, ", p = ", signif(p2, 2), " (", sig2, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      
      
      
      pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_scatters_holocentrics_against_chr_size/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_array_no_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_array_number, chr_sizes$size,
           xlab = "Repeat array number", ylab = "Chromosome size",
           main = "Repeat array number vs Chromosome size")
      
      # Add trendline
      model1 <- lm(size ~ repeat_array_number, data = chr_sizes)
      abline(model1, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result1 <- cor.test(chr_sizes$repeat_array_number, chr_sizes$size)
      r1 <- round(cor_result1$estimate, 2)
      p1 <- cor_result1$p.value
      sig1 <- if (p1 < 0.001) "***" else if (p1 < 0.01) "**" else if (p1 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r1, ", p = ", signif(p1, 2), " (", sig1, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_scatters_holocentrics_against_chr_size/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_repeat_bp_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_bp, chr_sizes$size,
           xlab = "Repeat base pairs", ylab = "Chromosome size",
           main = "Repeat base pairs vs Chromosome size")
      
      # Add trendline
      model2 <- lm(size ~ repeat_bp, data = chr_sizes)
      abline(model2, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result2 <- cor.test(chr_sizes$repeat_bp, chr_sizes$size)
      r2 <- round(cor_result2$estimate, 2)
      p2 <- cor_result2$p.value
      sig2 <- if (p2 < 0.001) "***" else if (p2 < 0.01) "**" else if (p2 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r2, ", p = ", signif(p2, 2), " (", sig2, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      
      
      
      
      
      
      
      
      
      
    } ### repeat content vs chromosome size and array number vs chromosome size
    
    }
    
    print("### strand switching ###")
    
    if(F) {
      arrays$array_strands <- "mixed"
      
      arrays$array_strands[arrays$plust_strand_perc >= 0.6] <- "+"
      arrays$array_strands[arrays$plust_strand_perc <= 0.4] <- "-"
      
      strand_vec <- arrays$array_strands
      strand_clean <- strand_vec[strand_vec %in% c("+", "-")]
      
      transitions <- paste0(head(strand_clean, -1), "->", tail(strand_clean, -1))
      
      # Tabulate transitions
      transition_table <- table(transitions)
      
      # Reformat into matrix
      transition_matrix <- matrix(0, nrow = 2, ncol = 2,
                                  dimnames = list(from = c("+", "-"), to = c("+", "-")))
      for (tr in names(transition_table)) {
        parts <- strsplit(tr, "->")[[1]]
        transition_matrix[parts[1], parts[2]] <- transition_table[tr]
      }
      chisq.test(transition_matrix)
      
      capture.output(chisq.test(transition_matrix), file = paste0(species_to_analyse[i], ".txt"))
      # 
      capture.output(chisq.test(transition_matrix), file = paste0("/home/pwlodzimierz/ToL/upload_files/33_chisq_test_results/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z],  ".txt"))
      
      arrays$array_strands <- "mixed"
      
      arrays$array_strands[arrays$plust_strand_perc >= 0.6] <- "+"
      arrays$array_strands[arrays$plust_strand_perc <= 0.4] <- "-"
      
      strand_vec <- arrays$array_strands
      strand_clean <- strand_vec[strand_vec %in% c("+", "-")]
      
      strand_binary <- ifelse(strand_clean == "+", 1, 0)
      
      # diff calculates difference between two elements of the vector, so it can be 0 if there's no
      # strand switch, -1 or 1 when there's 0-1 or 1-0 switches. Sum counts how many non-zero 
      # values are in the vector, so how many strand transistions happened
      count_switches <- function(vec) {
        sum(diff(vec) != 0)
      }
      
      # Count number of switches
      observed_switches <- count_switches(strand_binary)
      
      # permutate and get switches
      n_simulations <- 10000
      simulated_switches <- replicate(n_simulations, {
        random_vec <- sample(strand_binary)  # keep the 0/1 ratio
        count_switches(random_vec)
      })
      
      
      
      pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/33_strand_swithc_histograms/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], ".pdf"))
      hist(simulated_switches, breaks = 50, col = "lightgray",
           main = "",
           xlab = "Number of Switches")
      abline(v = observed_switches, col = "red", lwd = 2)
      
      mtext(text = paste0(species_to_analyse[i], " ", repeats_to_analyse[[i]][z]), side = 3, line = 1, adj = 0, cex = 1.5)
      mtext(text = paste0("n = ", length(strand_binary)), side = 3, line = 1, adj = 1, cex = 1.5)
      
      # Calculate critical values for p = 0.05 (two-tailed)
      critical_values <- quantile(simulated_switches, probs = c(0.025, 0.975))
      abline(v = critical_values, col = "blue", lty = 2, lwd = 1.5)  # Add dashed lines
      
      p_value <- mean(abs(simulated_switches - mean(simulated_switches)) >= abs(observed_switches - mean(simulated_switches)))
      
      cat("Observed switches:", observed_switches, "\n")
      cat("Mean simulated switches:", mean(simulated_switches), "\n")
      cat("P-value:", p_value, "\n")
      
      text(x = (par("usr")[1] + (par("usr")[2] - par("usr")[1])*0.1), y = (par("usr")[3] + (par("usr")[4] - par("usr")[3])*0.9), labels = paste0("P-value: "))
      text(x = (par("usr")[1] + (par("usr")[2] - par("usr")[1])*0.1), y = (par("usr")[3] + (par("usr")[4] - par("usr")[3])*0.86), labels = paste0(p_value))
      dev.off()
      
      
      
      
      
    } 
  }
  
  
}


# array width and spacing analysis ...

array_widths <- list()
array_spacings <- list()


for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  arrays <- data.frame()
  for(z in 1 : length(repeats_to_analyse[[i]])) {
    arrayst = read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"))
    arrayst$rep <- repeats_to_analyse[[i]][z]
    arrayst <- unique(arrayst)
    arrays = rbind(arrays, arrayst)
    
  }
  
  arrays$width <- arrays$end - arrays$start + 1
  arrays <- arrays[order(arrays$start, decreasing = FALSE),]
  arrays <- arrays[order(arrays$chromosome, decreasing = FALSE),]
  arrays$spacing <- c(arrays$start[2:nrow(arrays)] - arrays$end[1 : (nrow(arrays) - 1)], 0)
  arrays$spacing[arrays$spacing < min_distances[1]] <- 0
  # arrays$spacing <- log10(arrays$spacing)
  # arrays$width <- log10(arrays$width)
  # arrays$spacing[arrays$spacing == -Inf] <- NA
  arrays$spacing[arrays$spacing == 0] <- NA
  
  array_widths[[species_to_analyse[i]]] <- arrays$width[!is.na(arrays$width)]
  array_spacings[[species_to_analyse[i]]] <- arrays$spacing[!is.na(arrays$spacing)]
  
  
  
  
  
  
  
}

if(F) { # calculate to change order of plots 
  median_width <- NULL
  median_spacing <- NULL
  for(i in seq_along(array_widths)) {
    median_width <- c(median_width, median(array_widths[[i]]))
    median_spacing <- c(median_spacing, median(array_spacings[[i]]))
  }
  cor.test(median_width, median_spacing)
  
  reorder <- order(median_width)
  species_to_analyse <- species_to_analyse[reorder]
  repeats_to_analyse <- repeats_to_analyse[reorder]
  
}



# ... and similarity within vs between arrays

interarray_similarity <- list()
intraarray_similarity <- list()

for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  setwd(list.files("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", pattern = species_to_analyse[i], full.names = T))
  
  chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv")
  chr_sizes <- chr_sizes[grep(species_to_analyse[i], chr_sizes$assembly.name), ]
  chr_sizes <- chr_sizes[chr_sizes$is.chr == 1, ]
  
  chr_sizes <- chr_sizes[,-1]
  chr_sizes <- unique(chr_sizes)
  
  intraarray_similarity_species <- NULL
  interarray_similarity_species <- NULL
 
  for(z in 1 : length(repeats_to_analyse[[i]])) {
    arrays = read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"))
    repeats <- read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/33_repeats_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_repeats.csv"))
    
    arrays$width <- arrays$end - arrays$start + 1
    arrays <- arrays[order(arrays$start, decreasing = FALSE),]
    arrays <- arrays[order(arrays$chromosome, decreasing = FALSE),]
    
    
    
    ### similarity within and between arrays ###
    print("### similarity within arrays ###")
    
    sample_arrays <- sample(repeats$holo_array_ID, size = array_pairs_to_sample_for_similarity, replace = TRUE)
    
    
    similarity_score = rep(0, array_pairs_to_sample_for_similarity)
    
    for(j in 1 : array_pairs_to_sample_for_similarity) {
      # cat("array", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
      repeats_A <- repeats[repeats$holo_array_ID == sample_arrays[j], ]
      repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
      dist_matrix <- adist(repeats_A_sample$sequence)
      similarity_score[j] <- 100-100 * mean(dist_matrix[upper.tri(dist_matrix)]) /
        mean(repeats_A_sample$width)
    }
    
    intraarray_similarity_species <- c(intraarray_similarity_species, similarity_score)
    
    
    print("### similarity between arrays ###")
    
    sample_arrays_A <- sample(repeats$holo_array_ID, size = array_pairs_to_sample_for_similarity, replace = TRUE)
    sample_arrays_B <- sample(repeats$holo_array_ID, size = array_pairs_to_sample_for_similarity, replace = TRUE)
    
    # resample only the entries that collided with sample_arrays_A
    mismatch <- sample_arrays_B == sample_arrays_A
    while (any(mismatch)) {
      sample_arrays_B[mismatch] <- sample(repeats$holo_array_ID, 
                                          size = sum(mismatch), 
                                          replace = TRUE)
      mismatch <- sample_arrays_B == sample_arrays_A
    }
    
    similarity_score = rep(0, array_pairs_to_sample_for_similarity)
    
    
    for(j in 1 : array_pairs_to_sample_for_similarity) {
      # cat("chr", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
      repeats_A <- repeats[repeats$holo_array_ID == sample_arrays_A[j], ]
      repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
      repeats_B <- repeats[repeats$holo_array_ID == sample_arrays_B[j], ]
      repeats_B_sample <- repeats_B[sample(1 : nrow(repeats_B), replace = TRUE, size = repeats_to_sample_per_array),]
      similarity_score[j] <- 100-100*mean(adist(repeats_A_sample$sequence, repeats_B_sample$sequence)) /
        mean(c(repeats_A_sample$width, repeats_B_sample$width))
    }
    
    interarray_similarity_species <- c(interarray_similarity_species, similarity_score)
    
    mean(interarray_similarity_species)
    mean(intraarray_similarity_species)
  }
  interarray_similarity[[species_to_analyse[i]]] <- sample(interarray_similarity_species, array_pairs_to_sample_for_similarity, replace = F)
  intraarray_similarity[[species_to_analyse[i]]] <- sample(intraarray_similarity_species, array_pairs_to_sample_for_similarity, replace = F)
  
  
  
  
  
}

for(i in 1 : 13) {
  cat(
    mean(intraarray_similarity[[i]]), 
    mean(interarray_similarity[[i]]), mean(intraarray_similarity[[i]]) - mean(interarray_similarity[[i]]),
    "\n")
}



### plot array width and spacing violin plots 


widths_df <- do.call(rbind, lapply(names(array_widths), function(sp) {
  data.frame(width = array_widths[[sp]], species = sp, stringsAsFactors = FALSE)
}))
widths_df$species <- factor(widths_df$species, levels = names(array_widths))

spacings_df <- do.call(rbind, lapply(names(array_spacings), function(sp) {
  data.frame(spacing = array_spacings[[sp]], species = sp, stringsAsFactors = FALSE)
}))
spacings_df$species <- factor(spacings_df$species, levels = names(array_spacings))

theme_arrays <- function() {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, face = "italic"),
      plot.title         = element_text(size = 14, face = "bold", margin = margin(b = 6)),
      plot.subtitle      = element_text(size = 11, color = "grey40", margin = margin(b = 10)),
      legend.position    = "none",
      plot.margin        = margin(12, 16, 12, 12)
    )
}


p_widths <- ggplot(widths_df, aes(x = species, y = width, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA, width = 0.9) +
  geom_boxplot(width = 0.08, outlier.shape = NA,
               fill = "white", color = "grey30", linewidth = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_si("bp")),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title    = "Array width distribution by species",
    subtitle = "Violin = full distribution · Box = IQR + median  |  log10 y-axis",
    x        = NULL,
    y        = "Array width (bp)"
  ) +
  theme_arrays()

p_spacings <- ggplot(spacings_df, aes(x = species, y = spacing, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA, width = 0.9) +
  geom_boxplot(width = 0.08, outlier.shape = NA,
               fill = "white", color = "grey30", linewidth = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_si("bp")),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title    = "Inter-array spacing distribution by species",
    subtitle = "Violin = full distribution · Box = IQR + median  |  log10 y-axis",
    x        = NULL,
    y        = "Spacing between arrays (bp)"
  ) +
  theme_arrays()


setwd("/home/pwlodzimierz/ToL/upload_files/33_figure_plots")
ggsave("array_widths.pdf",   p_widths,   width = 10, height = 5, device = cairo_pdf)
ggsave("array_spacings.pdf", p_spacings, width = 10, height = 5, device = cairo_pdf)

intraarray_similarity_df <- do.call(rbind, lapply(names(intraarray_similarity), function(sp) {
  data.frame(intraarray_similarity = intraarray_similarity[[sp]], species = sp, stringsAsFactors = FALSE)
}))
intraarray_similarity_df$species <- factor(intraarray_similarity_df$species, levels = names(intraarray_similarity))

interarray_similarity_df <- do.call(rbind, lapply(names(interarray_similarity), function(sp) {
  data.frame(interarray_similarity = interarray_similarity[[sp]], species = sp, stringsAsFactors = FALSE)
}))
interarray_similarity_df$species <- factor(interarray_similarity_df$species, levels = names(interarray_similarity))

p_intraarray_similarity <- ggplot(intraarray_similarity_df, aes(x = species, y = intraarray_similarity, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA, width = 0.9) +
  geom_boxplot(width = 0.08, outlier.shape = NA,
               fill = "white", color = "grey30", linewidth = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title    = "Intra-array sequence similarity by species",
    subtitle = "Violin = full distribution \u00b7 Box = IQR + median",
    x        = NULL,
    y        = "Intra-array similarity"
  ) +
  theme_arrays()

p_interarray_similarity <- ggplot(interarray_similarity_df, aes(x = species, y = interarray_similarity, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA, width = 0.9) +
  geom_boxplot(width = 0.08, outlier.shape = NA,
               fill = "white", color = "grey30", linewidth = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title    = "Inter-array sequence similarity by species",
    subtitle = "Violin = full distribution \u00b7 Box = IQR + median",
    x        = NULL,
    y        = "Inter-array similarity"
  ) +
  theme_arrays()

ggsave("intraarray_similarity.pdf", p_intraarray_similarity, width = 10, height = 5, device = cairo_pdf)
ggsave("interarray_similarity.pdf", p_interarray_similarity, width = 10, height = 5, device = cairo_pdf)




# check whether different

species_list <- intersect(names(intraarray_similarity), names(interarray_similarity))

sim_test_results <- do.call(rbind, lapply(species_list, function(sp) {
  intra <- intraarray_similarity[[sp]]
  inter <- interarray_similarity[[sp]]
  
  if (length(intra) < 2 || length(inter) < 2) {
    return(data.frame(
      species = sp,
      n_intra = length(intra),
      n_inter = length(inter),
      median_intra = if (length(intra)) median(intra, na.rm = TRUE) else NA,
      median_inter = if (length(inter)) median(inter, na.rm = TRUE) else NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  test <- wilcox.test(intra, inter, exact = FALSE)
  
  data.frame(
    species = sp,
    n_intra = length(intra),
    n_inter = length(inter),
    median_intra = median(intra, na.rm = TRUE),
    median_inter = median(inter, na.rm = TRUE),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))

# correct for multiple testing across species
sim_test_results$p_adj <- p.adjust(sim_test_results$p_value, method = "BH")
sim_test_results$significant <- sim_test_results$p_adj < 0.05
sim_test_results <- sim_test_results[order(sim_test_results$p_adj), ]

print(sim_test_results)




### repeat similarity within arrays, chromosomes and genomes
## plot as a triple violin with significance marks















### metaplots

# Function to calculate density for a region using GenomicRanges coming from the Rhynchospora paper
# calc_density <- function(features_starts, features_ends, bin_starts, bin_ends) {
#   if(length(features_starts) == 0) return(rep(0, length(bin_starts)))
#   
#   # Create GRanges objects for bins and features
#   bins_gr <- GRanges(seqnames = "chr1", 
#                      ranges = IRanges(start = bin_starts, end = bin_ends))
#   
#   features_gr <- GRanges(seqnames = "chr1",
#                          ranges = IRanges(start = features_starts, end = features_ends))
#   
#   # Find overlaps and calculate coverage
#   density <- sapply(1:length(bins_gr), function(i) {
#     # Get overlapping features for this bin
#     overlaps <- findOverlaps(bins_gr[i], features_gr)
#     
#     if(length(overlaps) == 0) return(0)
#     
#     # Calculate intersection of bin with each overlapping feature
#     feature_indices <- subjectHits(overlaps)
#     bin_range <- ranges(bins_gr[i])
#     
#     total_overlap <- sum(sapply(feature_indices, function(idx) {
#       feature_range <- ranges(features_gr[idx])
#       intersect_range <- intersect(bin_range, feature_range)
#       width(intersect_range)
#     }))
#     
#     bin_length <- width(bins_gr[i])
#     return(100 * total_overlap / bin_length)  # Return percentage
#   })
#   
#   return(density)
# }
calc_density <- function(features_gr, bin_starts, bin_ends) {
  n_bins <- length(bin_starts)
  if (length(features_gr) == 0) return(rep(0, n_bins))
  
  bins_gr <- GRanges(seqnames = "chr1",
                     ranges = IRanges(start = bin_starts, end = bin_ends))
  
  overlaps <- findOverlaps(bins_gr, features_gr)
  
  if (length(overlaps) == 0) return(rep(0, n_bins))
  
  bin_hits     <- queryHits(overlaps)
  feature_hits <- subjectHits(overlaps)
  
  clipped <- pintersect(bins_gr[bin_hits], features_gr[feature_hits])
  overlap_widths <- width(clipped)
  
  # Sum overlapping bp per bin, then express as % of bin width
  total_overlap <- tapply(overlap_widths, bin_hits, sum)
  
  density <- rep(0, n_bins)
  density[as.integer(names(total_overlap))] <-
    100 * total_overlap / width(bins_gr)[as.integer(names(total_overlap))]
  
  return(density)
}

library(GenomicRanges)


genomes_means_vectors <- list()

# the plot contains of 100 bins in the flanks, each of them 250 bp, so going on to 25 kbp, 100 even array bins and 100 downstream bins
bins_flanks <- 100
bins_flank_size <- 250 # bp
min_array_size <- 100 # bp
for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  setwd(list.files("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", pattern = species_to_analyse[i], full.names = T))
  
  
  arrays <- data.frame()
  for(z in 1 : length(repeats_to_analyse[[i]])) {
    arrayst = read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"))
    arrayst$rep <- repeats_to_analyse[[i]][z]
    arrayst <- unique(arrayst)
    arrays = rbind(arrays, arrayst)

  }
  
  arrays$width <- arrays$end - arrays$start + 1
  arrays <- arrays[order(arrays$start, decreasing = FALSE),]
  arrays <- arrays[order(arrays$chromosome, decreasing = FALSE),]
  arrays$spacing <- c(arrays$start[2:nrow(arrays)] - arrays$end[1 : (nrow(arrays) - 1)], 0)
  arrays$spacing[arrays$spacing < min_distances[1]] <- 0
  arrays$spacing[arrays$spacing == 0] <- NA
  
  arrays <- arrays[!is.na(arrays$spacing),]


  repeats <- read.csv(list.files(pattern = "_repeats_filtered.csv", full.names = T))
  edta <- read.csv(list.files(pattern = "_edta_filtered.csv", full.names = T)[1])
  helixer <- read.csv(list.files(pattern = "_helixer_filtered.csv", full.names = T))
  
  helixer <- helixer[helixer$V3 == "CDS",]
  
  setwd("/home/pwlodzimierz/ToL/upload_files/33_figure_plots/single_metaplots")
  
  
  
  edta$V4 <- as.numeric(edta$V4)
  edta$V5 <- as.numeric(edta$V5)
  helixer$V4 <- as.numeric(helixer$V4)
  helixer$V5 <- as.numeric(helixer$V5)
  
  edta <- edta[!is.na(edta$V4),]
  edta <- edta[!is.na(edta$V5),]
  helixer <- helixer[!is.na(helixer$V4),]
  helixer <- helixer[!is.na(helixer$V4),]
  
  
  setwd(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", strsplit(species_to_analyse[i], split = "[.]")[[1]][1]))
  reps_files <- list.files(path = ".", pattern = "HOR_scored_repeats_with", full.names = T)

  features_genome_vectors <- list()
  for(j in 1 : length(unique(arrays$chromosome))) {
    chromosome <- unique(arrays$chromosome)[j]
    arrays_chr <- which(arrays$chromosome == chromosome)
    
    repeats_chr <- repeats[repeats$seqID == chromosome,]
    edta_chr <- edta[edta$V1 == chromosome,]
    helixer_chr <- helixer[helixer$V1 == chromosome,]
    
    edta_gr    <- GRanges("chr1", IRanges(edta_chr$V4,       edta_chr$V5))
    helixer_gr <- GRanges("chr1", IRanges(helixer_chr$V4,    helixer_chr$V5))
    repeats_gr <- GRanges("chr1", IRanges(repeats_chr$start, repeats_chr$end))
    
    reps_with_hors <- reps_files[grep(paste0(chromosome, ".csv"), reps_files, fixed = TRUE)]
    
    if(length(reps_with_hors) != 0) {
      reps_with_hors_dat <- data.frame()
      for(k in 1 : length(reps_with_hors)) {
        if(file.size(reps_with_hors[k]) != 0) reps_with_hors_dat <- rbind(reps_with_hors_dat, read.csv(file = reps_with_hors[k]))
      }
    }
    
    
    for(k in arrays_chr) {
      cat(i, "/", length(species_to_analyse), j, "/", length(unique(arrays$chromosome)), k, "/", nrow(arrays), "\n")
      # calculate bin positions
      upstream_bin_starts <- rev(arrays$start[k] - (1:bins_flanks) * 250)
      upstream_bin_ends   <- upstream_bin_starts + 249
      
      array_length <- arrays$end[k] - arrays$start[k] + 1
      if (array_length < min_array_size) next
      array_bin_size <- array_length / 100
      
      array_bin_starts      <- arrays$start[k] + floor((0:99)  * array_bin_size)
      array_bin_ends        <- arrays$start[k] + floor((1:100) * array_bin_size) - 1
      array_bin_ends[100]   <- arrays$end[k]
      
      downstream_bin_starts <- arrays$end[k] + (1:bins_flanks) * 250
      downstream_bin_ends   <- downstream_bin_starts + 249
      
      # Key change 5: concatenate all bins → single calc_density call per feature type
      all_starts <- c(upstream_bin_starts, array_bin_starts, downstream_bin_starts)
      all_ends   <- c(upstream_bin_ends,   array_bin_ends,   downstream_bin_ends)
      
      te_density     <- calc_density(edta_gr,    all_starts, all_ends)
      gene_density   <- calc_density(helixer_gr, all_starts, all_ends)
      repeat_density <- calc_density(repeats_gr, all_starts, all_ends)
      
      # HOR_score
      
      HOR_score_bins <- rep(NA, 300)
      if(nrow(reps_with_hors_dat) != 0) {
        bins_gr    <- GRanges("chr1", IRanges(all_starts, all_ends))
        starts_gr <- GRanges("chr1", IRanges(reps_with_hors_dat$start, width = 1))
        hits      <- findOverlaps(starts_gr, bins_gr, select = "first")
        
        reps_with_hors_dat$bin <- hits
        
        mean_by_bin <- tapply(reps_with_hors_dat$HOR_score, reps_with_hors_dat$bin, mean, na.rm = TRUE)
        
        HOR_score_bins[as.integer(names(mean_by_bin))] <- mean_by_bin
        HOR_score_bins[all_ends <= 1] = NA
        
        
      }
      
      
      
      # randomize order to reduce noise 
      if (sample(0:1, 1) == 1) reverse = TRUE else reverse = FALSE
      features_genome_vectors[[as.character(k)]]$te_density <- if (reverse) te_density else rev(te_density)
      features_genome_vectors[[as.character(k)]]$gene_density <- if (reverse) gene_density else rev(gene_density)
      features_genome_vectors[[as.character(k)]]$repeat_density <- if (reverse) repeat_density else rev(repeat_density)
      features_genome_vectors[[as.character(k)]]$HOR_score <- if (reverse) HOR_score_bins else rev(HOR_score_bins)
    }
    
    
    
  }
  
  te_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$te_density)), na.rm = TRUE)
  gene_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$gene_density)), na.rm = TRUE)
  repeat_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$repeat_density)), na.rm = TRUE)
  hor_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$HOR_score)), na.rm = TRUE)
  
  genomes_means_vectors[[species_to_analyse[i]]]$te <- te_per_bin_mean
  genomes_means_vectors[[species_to_analyse[i]]]$gene <- gene_per_bin_mean
  genomes_means_vectors[[species_to_analyse[i]]]$reps <- repeat_per_bin_mean
  genomes_means_vectors[[species_to_analyse[i]]]$hor <- hor_per_bin_mean
  
  
  mean_arrays_width <- mean(arrays$width)
  # mean_arrays_width <- mean_arrays_width/10
  if(mean_arrays_width > 150000) mean_arrays_width <- 150000
  
  pdf(file = paste0("single_plot_", species_to_analyse[i], ".pdf")) 
  plot(NULL,NULL, 
       xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size), 
       ylim = c(0,100),
       main = species_to_analyse[i])
  
  plot_bin_starts <- c(0:99 * bins_flank_size + 1, round(100 * bins_flank_size + 0:99 * (mean_arrays_width/100)), round(100 * bins_flank_size + 99 * (mean_arrays_width/100)) + 0:99 * bins_flank_size + 1)
  plot_bin_ends <- c(plot_bin_starts[2:length(plot_bin_starts)] - 1, plot_bin_starts[length(plot_bin_starts)] + 249)
  plot_bin_mids <- unlist(lapply(1 : 300, function(X) round(mean(c(plot_bin_ends[X], plot_bin_starts[X])))))
  
  lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
  lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
  lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
  abline(v = plot_bin_starts[101])
  abline(v = plot_bin_ends[200])
  dev.off()
}

save(genomes_means_vectors, 
     file = "/home/pwlodzimierz/ToL/upload_files/33_figure_plots/genomes_means_vectors.rds")
load("/home/pwlodzimierz/ToL/upload_files/33_figure_plots/genomes_means_vectors.rds")



mean_arrays_width = 25000


### metaplot
n_genomes <- length(genomes_means_vectors)
n_bins <- 300

# Robustly extract a field as a numeric vector of length n_bins,
# substituting NA-vectors for any NULL/missing/wrong-length entries.
# This preserves genome order -- a bad entry just becomes an all-NA row.
safe_extract <- function(lst, field, n_bins) {
  do.call(rbind, lapply(lst, function(x) {
    val <- x[[field]]
    if (is.null(val) || length(val) != n_bins) {
      rep(NA_real_, n_bins)
    } else {
      as.numeric(val)  # coerces NaN/NA uniformly, keeps them as missing
    }
  }))
}

te_matrix <- safe_extract(genomes_means_vectors, "te", n_bins)
gene_matrix <- safe_extract(genomes_means_vectors, "gene", n_bins)
repeat_matrix <- safe_extract(genomes_means_vectors, "reps", n_bins)
HOR_matrix <- safe_extract(genomes_means_vectors, "hor", n_bins)

# Helper: mean/SE per bin, ignoring NA/NaN, using the per-bin valid count
compute_mean_ci <- function(mat, ci_mult = 1.96) {
  bin_mean <- colMeans(mat, na.rm = TRUE)
  bin_sd   <- apply(mat, 2, sd, na.rm = TRUE)
  n_valid  <- apply(mat, 2, function(x) sum(!is.na(x)))  # is.na() also catches NaN
  bin_se   <- bin_sd / sqrt(n_valid)
  bin_ci   <- ci_mult * bin_se
  list(mean = bin_mean, ci = bin_ci, n = n_valid)
}

te_stats <- compute_mean_ci(te_matrix)
gene_stats <- compute_mean_ci(gene_matrix)
repeat_stats <- compute_mean_ci(repeat_matrix)
HOR_stats <- compute_mean_ci(HOR_matrix)

te_per_bin_mean <- te_stats$mean
gene_per_bin_mean <- gene_stats$mean
repeat_per_bin_mean <- repeat_stats$mean
HOR_per_bin_mean <- HOR_stats$mean
HOR_per_bin_mean[1:100] <- 0
HOR_per_bin_mean[201:300] <- 0

te_ci <- te_stats$ci
gene_ci <- gene_stats$ci
repeat_ci <- repeat_stats$ci
HOR_ci <- HOR_stats$ci
HOR_ci[1:100] <- 0
HOR_ci[201:300] <- 0

pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/33_figure_plots/metaplot_holocentromeres_HOR_scaling_5.pdf"))
plot(NULL, NULL,
     xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size),
     ylim = c(0, 100),
     main = "all holocentrics")

plot_bin_starts <- c(0:99 * bins_flank_size + 1, round(100 * bins_flank_size + 0:99 * (mean_arrays_width/100)), round(100 * bins_flank_size + 99 * (mean_arrays_width/100)) + 0:99 * bins_flank_size + 1)
plot_bin_ends <- c(plot_bin_starts[2:length(plot_bin_starts)] - 1, plot_bin_starts[length(plot_bin_starts)] + bins_flank_size - 1)
plot_bin_mids <- unlist(lapply(1 : 300, function(X) round(mean(c(plot_bin_ends[X], plot_bin_starts[X])))))

draw_ci_ribbon <- function(x, mean_vals, ci_vals, col) {
  valid <- !is.na(mean_vals) & !is.na(ci_vals)
  if (sum(valid) < 2) return(invisible(NULL))
  xv <- x[valid]
  upper <- (mean_vals + ci_vals)[valid]
  lower <- (mean_vals - ci_vals)[valid]
  polygon(x = c(xv, rev(xv)),
          y = c(upper, rev(lower)),
          col = adjustcolor(col, alpha.f = 0.2),
          border = NA)
}

HOR_scaling = 5

draw_ci_ribbon(plot_bin_mids, te_per_bin_mean, te_ci, "blue")
draw_ci_ribbon(plot_bin_mids, gene_per_bin_mean, gene_ci, "green")
draw_ci_ribbon(plot_bin_mids, repeat_per_bin_mean, repeat_ci, "red")
draw_ci_ribbon(plot_bin_mids, HOR_per_bin_mean * HOR_scaling, HOR_ci * HOR_scaling, "orange")

lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
lines(x = plot_bin_mids, y = HOR_per_bin_mean * HOR_scaling, col = "orange")

abline(v = plot_bin_starts[101])
abline(v = plot_bin_ends[200])
dev.off()



pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/33_figure_plots/metaplot_holocentromeres_rescaled.pdf"))
plot(NULL, NULL,
     xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size),
     ylim = c(0, 100),
     main = "all holocentrics")

plot_bin_starts <- c(0:99 * bins_flank_size + 1, round(100 * bins_flank_size + 0:99 * (mean_arrays_width/100)), round(100 * bins_flank_size + 99 * (mean_arrays_width/100)) + 0:99 * bins_flank_size + 1)
plot_bin_ends <- c(plot_bin_starts[2:length(plot_bin_starts)] - 1, plot_bin_starts[length(plot_bin_starts)] + bins_flank_size - 1)
plot_bin_mids <- unlist(lapply(1 : 300, function(X) round(mean(c(plot_bin_ends[X], plot_bin_starts[X])))))

draw_ci_ribbon <- function(x, mean_vals, ci_vals, col) {
  valid <- !is.na(mean_vals) & !is.na(ci_vals)
  if (sum(valid) < 2) return(invisible(NULL))
  xv <- x[valid]
  upper <- (mean_vals + ci_vals)[valid]
  lower <- (mean_vals - ci_vals)[valid]
  polygon(x = c(xv, rev(xv)),
          y = c(upper, rev(lower)),
          col = adjustcolor(col, alpha.f = 0.2),
          border = NA)
}

draw_ci_ribbon(plot_bin_mids, te_per_bin_mean, te_ci, "blue")
draw_ci_ribbon(plot_bin_mids, gene_per_bin_mean, gene_ci, "green")
draw_ci_ribbon(plot_bin_mids, repeat_per_bin_mean, repeat_ci, "red")
draw_ci_ribbon(plot_bin_mids, HOR_per_bin_mean, HOR_ci, "orange")

lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
lines(x = plot_bin_mids, y = HOR_per_bin_mean, col = "orange")

abline(v = plot_bin_starts[101])
abline(v = plot_bin_ends[200])
dev.off()



if(F) {
  ### metaplot
  te_per_bin_mean <- vector(mode = "numeric", length = 300)
  gene_per_bin_mean <- vector(mode = "numeric", length = 300)
  repeat_per_bin_mean <- vector(mode = "numeric", length = 300)
  for(i in 1 : length(genomes_means_vectors)) {
    te_per_bin_mean <- te_per_bin_mean + genomes_means_vectors[[j]]$te
    gene_per_bin_mean <- gene_per_bin_mean + genomes_means_vectors[[j]]$gene
    repeat_per_bin_mean <- repeat_per_bin_mean + genomes_means_vectors[[j]]$reps
  }
  te_per_bin_mean <- te_per_bin_mean / 21
  gene_per_bin_mean <- gene_per_bin_mean / 21
  repeat_per_bin_mean <- repeat_per_bin_mean / 21
  
  
  mean_arrays_width = 15000
  
  pdf(file = paste0("metaplot_holocentric_arrays.pdf")) 
  plot(NULL,NULL, 
       xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size), 
       ylim = c(0,100),
       main = "all holocentrics")
  
  plot_bin_starts <- c(0:99 * bins_flank_size + 1, round(100 * bins_flank_size + 0:99 * (mean_arrays_width/100)), round(100 * bins_flank_size + 99 * (mean_arrays_width/100)) + 0:99 * bins_flank_size + 1)
  plot_bin_ends <- c(plot_bin_starts[2:length(plot_bin_starts)] - 1, plot_bin_starts[length(plot_bin_starts)] + 249)
  plot_bin_mids <- unlist(lapply(1 : 300, function(X) round(mean(c(plot_bin_ends[X], plot_bin_starts[X])))))
  
  lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
  lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
  lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
  abline(v = plot_bin_starts[101])
  abline(v = plot_bin_ends[200])
  dev.off()

}
























