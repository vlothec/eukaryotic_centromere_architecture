

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
library(plot.matrix)
library(ape)

repeats_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/2026_filtered_TRASH", full.names = T)

satellites <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")
metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv")
phyla <- read.csv("/home/pwlodzimierz/ToL/Metadata/fasta_phyla.csv")

metadata <- metadata[metadata$is.chr == 1,]

species <- unique(satellites$fasta)


taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 166
print(i)

plots_no <- 50
long_rep_plots_no <- 5
rep_per_plot <- 50


{
  estimate_base_distance <- function(distances, min_count = 2) {
    # optional: drop very rare values as likely noise
    tab <- table(distances)
    keep <- as.numeric(names(tab)[tab >= min_count])
    d <- distances[distances %in% keep]
    if (length(d) == 0) d <- distances  # fallback if filtering removes everything
    
    # start from the smallest plausible value
    b <- min(d)
    
    repeat {
      k <- round(d / b)
      k[k == 0] <- 1  # guard against div-by-zero
      b_new <- sum(k * d) / sum(k^2)   # least-squares fit of d ~ k*b through origin
      if (abs(b_new - b) < 1e-6) break
      b <- b_new
    }
    b
  }
  
  
  find_best_period <- function(div_matrix) {
    div_vals <- as.numeric(div_matrix)
    
    lowest_div_threshold <- sort(div_vals)[round(length(div_vals) / 20)] 
    idx <- which(div_matrix <= lowest_div_threshold, arr.ind = TRUE)
    
    distances <- as.numeric(sapply(1 : nrow(idx), function(X) idx[X,1] - idx[X,2]))
    distances <- distances[distances > 0]
    
    tab_dist <- table(distances)
    tab_dist <- sort(tab_dist, decreasing = T)
    
    if(length(tab_dist) == 0) return(rep_per_plot)
    if(as.numeric(names(tab_dist)[1]) == 1) return(rep_per_plot)
    estimate_base_distance(distances)
  }
  
  
  
  
  align_and_plot_tree <- function(repeats_seq, sample_id = "", k = NULL, outdir = "trees") {
    # dir.create(outdir, showWarnings = FALSE)
    
    seqs <- as.character(repeats_seq)
    names(seqs) <- seq_along(seqs)  # labels = indices
    
    fasta_in  <- tempfile(fileext = ".fasta")
    fasta_out <- tempfile(fileext = ".fasta")
    writeLines(paste0(">", names(seqs), "\n", seqs), fasta_in)
    
    if (nzchar(Sys.which("mafft"))) {
      system2("mafft", args = c("--auto", "--quiet", fasta_in), stdout = fasta_out)
    } else {
      library(Biostrings); library(msa)
      dna <- readDNAStringSet(fasta_in)
      aln <- msaClustalOmega(dna)
      writeXStringSet(as(aln, "DNAStringSet"), filepath = fasta_out)
    }
    
    aln  <- read.FASTA(fasta_out)
    d    <- dist.dna(aln, model = "raw", pairwise.deletion = TRUE)
    tree <- njs(d)
    
    # --- sanitize distance matrix for hclust ---
    bad <- !is.finite(d)
    if (any(bad)) {
      warning(sprintf("[%s] %d non-finite distance(s) in dist.dna output; imputing with max observed distance",
                      sample_id, sum(bad)))
      d[bad] <- if (any(is.finite(d))) max(d[is.finite(d)]) else 0
    }
    
    hc <- tryCatch(
      hclust(d, method = "average"),
      error = function(e) {
        warning(sprintf("[%s] hclust failed even after sanitizing: %s", sample_id, conditionMessage(e)))
        NULL
      }
    )
    
    if (is.null(k)) k <- 10
    
    if (!is.null(hc)) {
      clusters <- cutree(hc, k = k)
    } else {
      clusters <- setNames(rep(1L, length(seqs)), names(seqs))
    }
    names(clusters) <- names(seqs)
    
    pal <- rainbow(max(clusters))
    
    label_cols <- pal[clusters[as.character(seq_along(seqs))]]
    names(label_cols) <- NULL
    
    tip_cols <- pal[clusters[tree$tip.label]]
    tip_cols <- adjustcolor(tip_cols, alpha.f = 0.6)
    
    # png(file.path(outdir, paste0(sample_id, "_tree.png")), width = 900, height = 900, res = 150)
    plot(tree, type = "unrooted", cex = 1.6, main = sample_id, show.tip.label = FALSE)
    tiplabels(pch = 19, col = tip_cols, cex = 2.6)
    add.scale.bar(length = 0.05)  # fixed value, same across all trees
    # dev.off()
    
    unlink(c(fasta_in, fasta_out))
    list(tree = tree, clusters = clusters, hc = hc, label_cols = label_cols)
  }
  
  
}



for(i in i) {
  species_name <- strsplit(species[i], split = "[.]")[[1]][1]
  species_sats <- satellites[satellites$fasta == species[i],]
  
  if(species_sats$is_holocentric[1]) {
    cat(i, species_name, "holocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
    next
  }
  cat(i, species_name, "monocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
  
  if(i == 166) {
    repeats <- read.csv(file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_output/GCA_000001405.29_GRCh38.p14_genomic.fa/GCA_000001405.29_GRCh38.p14_genomic.fa_repeats_with_seq.csv")
  } else {
    repeats <- read.csv(file = repeats_files[grep(species_name, repeats_files)])
  }
  
  chromosomes_metadata <- metadata[grep(species_name, metadata$assembly.name),]
  chrs <- chromosomes_metadata$chromosome.name
  if(i == 166) chrs <-  c("CM000663.2", "CM000664.2", "CM000665.2", "CM000666.2", "CM000667.2",
                          "CM000668.2", "CM000669.2", "CM000670.2", "CM000671.2", "CM000672.2",
                          "CM000673.2", "CM000674.2", "CM000675.2", "CM000676.2", "CM000677.2",
                          "CM000678.2", "CM000679.2", "CM000680.2", "CM000681.2", "CM000682.2",
                          "CM000683.2", "CM000684.2", "CM000685.2", "CM000686.2")
  
  repeats <- repeats[repeats$seqID %in% chrs,]
  
  for(j in 1 : nrow(species_sats)) {
    
    
    repeats_family <- repeats[repeats$new_class == species_sats$TRASH_new_class[j],]
    
    if(nrow(repeats_family) == 0) stop("No repeats in the family?")
    
    chromosomes <- table(repeats_family$seqID)
    
    chromosomes <- chromosomes[chromosomes >= 10]
    
    sample_chrs <- sample(1:length(chromosomes), plots_no, replace = TRUE)
    if(i == 166) sample_chrs <- sample(1:length(chromosomes), plots_no, replace = TRUE)
    
    mean_rep_width <- mean(nchar(repeats_family$sequence))
    if(mean_rep_width > 500) sample_chrs <- sample(1:length(chromosomes), long_rep_plots_no, replace = TRUE)
    
    
    pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/66_HOR_matrix_plots/", species_name, "_", species_sats$TRASH_new_class[j], ".pdf"), 
        width = 10, height = 5, onefile = T)
    par(mfrow = c(1,2), oma = c(0,0,0,0))
    
    
    for(k in sample_chrs) {
      try({
        cat(k, "\n")
        
        rep_per_plot_t <- rep_per_plot
        if(chromosomes[k] < rep_per_plot) {
          rep_per_plot_t <- chromosomes[k]
        }
        
        first_rep_ID <- sample(1 : (chromosomes[k] - rep_per_plot_t), 1)
        
        repeats_seq <- repeats_family$sequence[repeats_family$seqID == names(chromosomes[k])][first_rep_ID : (first_rep_ID + rep_per_plot_t - 1)]
        
        mean_rep_width <- mean(nchar(repeats_seq))
        
        adist_matrix <- adist(repeats_seq)
        
        div_matrix <- 100 * adist_matrix / mean_rep_width
        div_matrix[div_matrix > 20] = 20
        
        best_period <- find_best_period(div_matrix)
        
        if(best_period >= rep_per_plot_t) best_period <- NULL
        
        ### PLOTTING: explicit 3-panel layout 
        layout(matrix(c(1,2,3), nrow = 1), widths = c(4, 4, 0.6))
        
        ### 1. tree
        tree_data <- align_and_plot_tree(repeats_seq, k = best_period)
        label_cols <- tree_data$label_cols
        # label_cols <- adjustcolor(label_cols, alpha.f = 0.8)
        
        ### 2. matrix (key disabled so it doesn't override the layout)
        par(mar = c(1,1,4,1))
        plot(div_matrix,
             col = topo.colors,
             breaks = 0:20,
             border = NA,
             key = NULL,
             axis.col = list(side = 1, labels = FALSE),
             axis.row = list(side = 2, labels = FALSE),
             main = paste0(species_sats$TRASH_new_class[j], " ", names(chromosomes[k]), " ",
                           first_rep_ID, ":", (first_rep_ID + rep_per_plot_t - 1)))
        
        usr <- par("usr")
        x_offset <- diff(usr[1:2]) * 0.01
        y_offset <- diff(usr[3:4]) * 0.01
        
        points(x = 1:ncol(div_matrix), y = rep(usr[3] - y_offset, ncol(div_matrix)),
               pch = 16, col = label_cols, xpd = NA, cex = 2)
        points(x = rep(usr[1] - x_offset, nrow(div_matrix)), y = 1:nrow(div_matrix),
               pch = 16, col = label_cols, xpd = NA, cex = 2)
        ### 3. manual color key, since plot.matrix's own key is disabled
        par(mar = c(3,2,2,2))
        key_breaks <- 0:20
        key_cols <- topo.colors(length(key_breaks) - 1)
        plot(0, 0, type = "n", xlim = c(0,1), ylim = range(key_breaks),
             axes = FALSE, xlab = "", ylab = "")
        rasterImage(as.raster(matrix(rev(key_cols), ncol = 1)),
                    0, min(key_breaks), 1, max(key_breaks))
        axis(4, at = key_breaks, labels = as.integer(key_breaks), las = 1)
      })
    }
    
    dev.off()
    
  }
  
  
}


