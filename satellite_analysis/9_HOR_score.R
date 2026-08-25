setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out") 

plot_regim = FALSE

suppressWarnings(suppressPackageStartupMessages(library(data.table)))

suppressWarnings(suppressPackageStartupMessages(library(IRanges)))

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}

cen_families <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")

names(cen_families) <- c("ID", "assembly", "repeat_name", "repeat_count", "is_holocentric")

summary_df <- cbind(cen_families, data.frame(repeat_mean_width = rep("", nrow(cen_families)),
                                           repeat_total_bp = rep("", nrow(cen_families)),
                                           chr_no_analysed = rep("", nrow(cen_families)),
                                           chr_names_analysed = rep("", nrow(cen_families)),
                                           mean_repetitiveness = rep("", nrow(cen_families))))

hor_dirs <- paste0("./", unlist(lapply(summary_df$assembly, function(X) strsplit(X, split = "[.]")[[1]][1])))

rep_dir_a <- "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/"


max_hor_distance <- 60
bin_size <- 100 # repeats per bin
threshold_canonical_HOR <- bin_size * 0.5 # half of the repeats must have a signal

repeat_analysis <- FALSE



taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 265
print(i)

# for every genome
# for(i in seq_along(hor_dirs)) 
{
  
  # if(grepl("GRCh38", summary_df$assembly[i])) {
  #   repeats <- read.csv(file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/GCA_000001405.29_GRCh38.p14_genomic.fa/GCA_000001405.29_GRCh38.p14_genomic.fa_repeats_with_seq.csv")
  #   
  # } else {
  #   repeats <- read.csv(file = paste0(rep_dir_a, summary_df$assembly[i], ".fa/", summary_df$assembly[i], "_repeats_filtered.csv"))
  #   
  # }
  
  cat(i, summary_df$assembly[i], "\n")
  
  hor_files <- list.files(path = hor_dirs[i])
  hor_files <- hor_files[grep(summary_df$repeat_name[i], hor_files)]
  if(length(hor_files) == 0) {
    cat("No HOR files found\n")
    stop()
  }
  
  all_hors_files <- hor_files[grep("HORs_", hor_files)]
  all_hors_files <- all_hors_files[grep("csv", all_hors_files)]
  repeats_with_hors_files <- hor_files[grep("repeats_with_hors_", hor_files)]
  repeats_with_hors_files <- repeats_with_hors_files[!grepl("HOR_scored_repeats", repeats_with_hors_files)]
  if(length(repeats_with_hors_files) == 0) {
    cat("No repeats with HOR files found\n")
    stop()
  }
  if(length(all_hors_files) == 0) {
    cat("No HOR files found\n")
    stop()
  }
  repeat_class <- summary_df$repeat_name[i]
  
  rep_out <- "/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/"
  dir_out <- paste0(rep_out, hor_dirs[i])
  if(!dir.exists(dir_out)) dir.create(dir_out, showWarnings = F)
  
  for(j in seq_along(repeats_with_hors_files)) {
    
    if(!file.exists(paste0(dir_out, "/HOR_scored_", repeats_with_hors_files[j]))) {
      
      chromosome_name <- strsplit(repeats_with_hors_files[j], split = paste0(repeat_class, "_"))[[1]][2]
      chromosome_name <- strsplit(chromosome_name, split = ".csv")[[1]][1]
      cat(i, "/", length(hor_dirs), summary_df$assembly[i], j, "/", length(repeats_with_hors_files), 
          chromosome_name, repeat_class, "")
      
      
      # if(file.size(paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))/1048576 > 5000) {
      #   cat(" file bigger than 5000 MB, will process later \n")
      #   next
      # }
      
      repeats_chr <- read.csv(file = paste0(hor_dirs[i], "/", repeats_with_hors_files[j]))
      
      rep_out <- "/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/plots/"
      
      
      if(!file.exists(paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))) {
        cat("\n")
        next
      }
      hors <- read.csv(file = paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))
      
      
      if(F) { # legacy code for comparison
        
        # expand both block columns
        long_data1 <- data.table(
          element = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_A[X] : hors$end_A[X], hors$block.size.in.units[X])  )),
          partner = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_B[X] : hors$end_B[X], each = hors$block.size.in.units[X])   ))
        )
        
        long_data2 <- data.table(
          element = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_B[X] : hors$end_B[X], hors$block.size.in.units[X])  )),
          partner = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_A[X] : hors$end_A[X], each = hors$block.size.in.units[X])   ))
        )
        
        long_data1 <- unique(long_data1)
        long_data2 <- unique(long_data2)
        
        
        # combine
        long_data <- rbind(long_data1, long_data2)
        long_data <- unique(long_data)
        
        remove(long_data1, long_data2)
        
        # remove self-interactions
        long_data <- long_data[element != partner]
        
        # caluclate unique interactors
        unique_interactions <- long_data[, .(num_interactors = uniqueN(partner)), by = element]
        
        # get results for each repeat
        result_vector <- integer(nrow(repeats_chr))
        result_vector[unique_interactions$element] <- unique_interactions$num_interactors
        
        # normalise to get HOR score
        repeats_chr$HOR_score <- 100 * result_vector / nrow(repeats_chr)
        
        
        
        repeats_chr$hors_formed_count = 0
      }
      # Work only with the four columns we need; convert to integer to save RAM
      hors_dt <- as.data.table(hors)[, .(
        start_A = as.integer(start_A),
        end_A   = as.integer(end_A),
        start_B = as.integer(start_B),
        end_B   = as.integer(end_B),
        bs      = as.integer(block.size.in.units)
      )]
      
      # Build the pair table by expanding only the *ranges*, not the block size
      # (we only care about unique element-partner pairs, so multiplicity is irrelevant)
      make_pairs <- function(sA, eA, sB, eB) {
        # one row per HOR
        data.table(
          element = Map(seq.int, sA, eA),
          partner = Map(seq.int, sB, eB)
        )[, .(
          element = unlist(element),
          partner = unlist(partner)
        )]
      }
      
      # Direction A→B
      pairs1 <- make_pairs(hors_dt$start_A, hors_dt$end_A,
                           hors_dt$start_B, hors_dt$end_B)
      # Direction B→A
      pairs2 <- make_pairs(hors_dt$start_B, hors_dt$end_B,
                           hors_dt$start_A, hors_dt$end_A)
      
      long_data <- unique(rbind(pairs1, pairs2))
      rm(pairs1, pairs2, hors_dt); gc()
      
      # Drop self-pairs
      long_data <- long_data[element != partner]
      
      # Count unique partners per element
      unique_interactions <- long_data[, .(num_interactors = uniqueN(partner)), by = element]
      rm(long_data); gc()
      
      # Write into the result vector
      result_vector <- integer(nrow(repeats_chr))
      result_vector[unique_interactions$element] <- unique_interactions$num_interactors
      rm(unique_interactions)
      
      repeats_chr$HOR_score <- 100 * result_vector / nrow(repeats_chr)
      rm(result_vector); gc()
      
      # --- 2. hors_formed_count (coverage) – already efficient --------------------
      ir <- IRanges(
        start = c(hors$start_A, hors$start_B),
        end   = c(hors$end_A,   hors$end_B)
      )
      n <- nrow(repeats_chr)
      cov <- coverage(ir, width = n)
      repeats_chr$hors_formed_count <- as.integer(cov)
      rm(ir, cov); gc()
      
      repeats_chr$hors_formed_tot_rep_normalised <-
        100 * repeats_chr$hors_formed_count / nrow(repeats_chr)
      
      # write
      write.csv(repeats_chr,
                file = paste0(dir_out, "/HOR_scored_", repeats_with_hors_files[j]),
                row.names = FALSE)
      
      
     
      
      
      
      
      
    }
    
    if(plot_regim) {
      
      
      if(!repeat_analysis) {
        if(file.exists(file.path(rep_out,
                                 paste0("new_HORs_lines_simple_", summary_df$assembly[i], "-", repeat_class, "-", chromosome_name, "-repeat_number_", nrow(repeats_chr), ".png", collapse = "")))) {
          if(file.size(file.path(rep_out,
                                 paste0("new_HORs_lines_simple_", summary_df$assembly[i], "-", repeat_class, "-", chromosome_name, "-repeat_number_", nrow(repeats_chr), ".png", collapse = ""))) > 0) {
            cat("\n")
            next
          }
        }
        
      }
      cat("plotting...\n")
      
      
      png(filename = file.path(rep_out, 
                               paste0("new_HORs_lines_simple_", summary_df$assembly[i], "-", repeat_class, "-", chromosome_name, "-repeat_number_", nrow(repeats_chr), ".png", collapse = "")), 
          width = 3000, height = 5000, pointsize = 40)
      
      par(mar = c(4, 4, 4, 4), oma = c(1, 1, 1, 1))
      
      par(fig=c(0,1,0.4,1))
      ax.len <- nrow(repeats_chr)
      
      hors = hors[order(hors$SNV_per_kbp, decreasing = TRUE), ]
      ### FILTERING
      # hors <- hors[hors$SNV_per_kbp < 30,]
      
      
      unit.name <- "repeat ID"
      plot(x = NULL, y = NULL,
           xlab = paste0(chromosome_name, ", ", unit.name),
           ylab = paste0(chromosome_name, ", ", unit.name),
           xlim = c(0, ax.len),
           ylim = c(0, ax.len),
           pch = 19, cex = 0.1,
           main = paste0("HORs ", summary_df$assembly[i], " ", repeat_class, " ", chromosome_name, " ", nrow(repeats_chr), " repeats"))
      
      SNV_per_kbp_in_red = 20
      
      colours_SNV <- colorRampPalette(c("green", "yellow", "red"))(length(hors$SNV_per_kbp)) [findInterval(hors$SNV_per_kbp, seq(0, SNV_per_kbp_in_red, length.out = length(hors$SNV_per_kbp)))]
      while (TRUE) {
        if (ax.len >= 100000) {lwd_plot <- 2; break}
        if (ax.len >= 50000) {lwd_plot <- 3; break}
        if (ax.len >= 25000) {lwd_plot <- 4; break}
        if (ax.len >= 12500) {lwd_plot <- 5; break}
        lwd_plot <- 6
        break
      }
      for (k in seq_len(nrow(hors))) {
        lines(x = c(hors$start_A[k], hors$end_A[k]),
              y = c(hors$start_B[k], hors$end_B[k]),
              pch = 19, lwd = lwd_plot, col = colours_SNV[k])
      }
      
      
      ### check the most common HOR distance (between 3 and max_hor_distance)
      
      par(fig=c(0,1,0.2,0.4), new = TRUE)
      bins <- round(nrow(repeats_chr) / bin_size)
      if(bins < 1) bins <- 1
      
      hors_distance <- hors[(hors$start_B - hors$start_A) <= max_hor_distance,]
      
      hors_distance$dist <- hors_distance$start_B - hors_distance$start_A
      hors_distance <- hors_distance[order(hors_distance$start_A, decreasing = FALSE), ]
      
      hors_distances_counts <- rep(list(rep(0,bins)), max_hor_distance)
      hors_distances_bin_repeats <- ceiling((0:bins)*(nrow(repeats_chr)/bins))[-1] - ceiling((0:bins)*(nrow(repeats_chr)/bins))[-(bins+1)]
      
      # for(k in seq_len(nrow(hors_distance))) {
      #   distance <- hors_distance$start_B[k] - hors_distance$start_A[k]
      #   for(l in 0 : (hors_distance$block.size.in.units[k]-1)) {
      #     which_bin <- ceiling((hors_distance$start_A[k] + l) / (nrow(repeats_chr) / bins))
      #     hors_distances_counts[[distance]][which_bin] = hors_distances_counts[[distance]][which_bin] + 1
      #     
      #     which_bin <- ceiling((hors_distance$start_B[k] + l) / (nrow(repeats_chr) / bins))
      #     hors_distances_counts[[distance]][which_bin] = hors_distances_counts[[distance]][which_bin] + 1
      #     
      #   }
      #   
      # }
      ### vectorised implementation of the above:
      
      distances <- hors_distance$start_B - hors_distance$start_A
      block_sizes <- hors_distance$block.size.in.units
      
      row_indices <- rep(seq_len(nrow(hors_distance)), times = block_sizes)
      offsets     <- sequence(block_sizes) - 1
      
      which_bin_A <- ceiling((hors_distance$start_A[row_indices] + offsets) / bin_size)
      which_bin_B <- ceiling((hors_distance$start_B[row_indices] + offsets) / bin_size)
      
      dist_values <- distances[row_indices]
      
      for (d in unique(dist_values)) {
        mask <- dist_values == d
        
        bin_A_counts <- tabulate(which_bin_A[mask], nbins = bins)
        bin_B_counts <- tabulate(which_bin_B[mask], nbins = bins)
        
        hors_distances_counts[[d]] <- hors_distances_counts[[d]] + bin_A_counts + bin_B_counts
      }
      
      # each vector in hors_distances_counts is per bin value of HORs at specific distance.
      # since there are 100 repeats per bin, each bin has 100 possible HOR seeds to start
      # a new HOR. with max distance of 60, a theoretical score of a window that constitutes
      # of perfect 60-mer set is 100, but a perfect 30-mer will have a score of 200, as each
      # repeat/seed gets scored once on it's "own" HOR, and once on the "previous" one
      # This means the score can go as high as 6000, if it's a stretch of identical repeats.
      # In reality, the score above 100 should be a very strong indicator of a canonical
      # HOR. Each potential canonical HOR signal can be "strengthened" by adding the 
      # signal from its multiplication. For example, a 2-mer signal at distance 2, can be
      # strengthened by adding the signal from distance 4. Although the signal from
      # distance 4 is already accounted in distance 2, since it is that "previous"
      # signal. So adding these values would not make sense. What would make sense, is
      # deducting the value of the first spotted canonical HOR from its multiplications,
      # so that if a window has more than one canonical HOR, it can be spotted, without 
      # counting the multiplications over and over again. 
      # Since the HOR identification is already normalising for repeat size, there can be
      # a global threshold value above which the canonical HOR is identified, with comparable
      # results between species.
      
      
      matrix_canonical_hors <- matrix(nrow = bins, ncol = max_hor_distance, data = F)
      temp_hors_distances_counts <- hors_distances_counts
      for(d in 1 : max_hor_distance) {
        for(bin_id in 1 : bins) {
          if(temp_hors_distances_counts[[d]][bin_id] > threshold_canonical_HOR) {
            matrix_canonical_hors[bin_id, d] <- TRUE
            d_mult <- (1:max_hor_distance) * d
            d_mult <- d_mult[d_mult <= max_hor_distance]
            d_mult <- d_mult[-1]
            if(length(d_mult) == 0) next
            for(d2 in d_mult) {
              temp_hors_distances_counts[[d2]][bin_id] <- temp_hors_distances_counts[[d2]][bin_id] - temp_hors_distances_counts[[d]][bin_id]
              if(temp_hors_distances_counts[[d2]][bin_id] < 0) temp_hors_distances_counts[[d2]][bin_id] <- 0
            }
            
          }
        }
      }
      vector_canonical_HORs <- vector(mode = "character", length = bins)
      for(k in seq_along(vector_canonical_HORs)) {
        vector_canonical_HORs[k] <- paste0(which(matrix_canonical_hors[k,]), collapse = ";")
      }
      table_canonical_HORs <- table(vector_canonical_HORs)
      table_canonical_HORs <- sort(table_canonical_HORs, decreasing = TRUE)
      out_dir_mat <- "/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/distance_data/matrices/matrix_canonical_hors"
      out_dir_vect <- "/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/distance_data/short_vectors/vector_canonical_hors"
      out_dir_mat <- paste(out_dir_mat, summary_df$assembly[i], summary_df$repeat_name[i], chromosome_name, sep = ";")
      out_dir_vect <- paste(out_dir_vect, summary_df$assembly[i], summary_df$repeat_name[i], chromosome_name, sep = ";")
      write.csv(x = matrix_canonical_hors, file = paste0(out_dir_mat, ".csv"))
      write.csv(x = table_canonical_HORs, file = paste0(out_dir_vect, ".csv"))
      
      plot(x = NULL, y = NULL, xlim = c(0,bins+1), ylim = c(0,max_hor_distance), 
           xlab = "repeat bin", ylab = "HOR distance frequency")
      for(k in 1 : max_hor_distance) {
        points(1:bins, rep(k,bins), cex = hors_distances_counts[[k]]/hors_distances_bin_repeats/3, pch = 15, col = "grey")
      }
      for(k in 1 : max_hor_distance) {
        points(1:bins, rep(k,bins), cex = temp_hors_distances_counts[[k]]/hors_distances_bin_repeats/3, pch = 15)
      }
      
      for(k in 1 : max_hor_distance) {
        temp_hors_distances_counts[[k]][!matrix_canonical_hors[,k]] = 0
        points(1:bins, rep(k,bins), cex = temp_hors_distances_counts[[k]]/hors_distances_bin_repeats/3, pch = 15, col = "red")
      }
      abline(h = 1:(max_hor_distance+1) - 0.5, v = 1:(bins+1)- 0.5, col = "gray")
      abline(h = seq(-0.5,max_hor_distance+0.5,by=5), col = "#0088ee", cex = 2)
      abline(h = seq(0.5,max_hor_distance+0.5,by=5), col = "#0088ee", cex = 2)
      abline(h = seq(-0.5,max_hor_distance+0.5,by=10), col = "#0000ee", cex = 4)
      abline(h = seq(0.5,max_hor_distance+0.5,by=10), col = "#0000ee", cex = 4)
      
      
      
      # add HOR score plot
      par(fig=c(0,1,0.0,0.2), new = TRUE)
      
      plot(y = repeats_chr$HOR_score, 
           x = seq_along(repeats_chr$HOR_score), type = "p", 
           pch = 16, col = "grey", ylim = c(0,100),
           xlab = "repeat ID", ylab = "HOR normalised count, %")
      points(y = ma(x = repeats_chr$HOR_score, n = round(nrow(repeats_chr) / bins)), 
             x = seq_along(repeats_chr$HOR_score), type = "l", ylim = c(0,100))
      
      dev.off()
      cat("\n")
      
    }
    
    
    
  }
  
  
  
  
}











