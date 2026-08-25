

library(GenomicRanges)
library(seqinr)
library(ggplot2)

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


repeats_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/2026_filtered_TRASH", full.names = T)
edta_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/2026_EDTA_filtered", full.names = T)
helixer_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/2026_helixer_filtered", full.names = T)
cen_coords_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/55_cen_coords", pattern = "genome_metadata", full.names = T)

satellites <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")
metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv")
phyla <- read.csv("/home/pwlodzimierz/ToL/Metadata/fasta_phyla.csv")

bins_flanks <- 100
bins_flank_size <- 25000 # bp


species <- unique(satellites$fasta)

genomes_means_vectors <- list()
mean_centromeres_sizes <- NULL
for(i in 1 : length(species)) {
  # try({
    species_name <- strsplit(species[i], split = "[.]")[[1]][1]
    species_sats <- satellites[satellites$fasta == species[i],]
    
    if(species_sats$is_holocentric[1]) {
      cat(i, species_name, "holocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
      next
    }
    cat(i, species_name, "monocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
    
    
    repeats <- read.csv(file = repeats_files[grep(species_name, repeats_files)])
    edta <- read.csv(file = edta_files[grep(species_name, edta_files)])
    helixer <- read.csv(file = helixer_files[grep(species_name, helixer_files)])
    if(species_name == "rosCan_S27_v1") {
      cen_coords <- read.csv(file = "/home/pwlodzimierz/ToL/upload_files/55_cen_coords/rosCan_genome_metadata.csv")
    } else {
      cen_coords <- read.csv(file = cen_coords_files[grep(species_name, cen_coords_files)])
    }
    
    
    # helixer <- helixer[helixer$V3 == "CDS",] #TODO consider using gene vs CDS
    helixer <- helixer[helixer$V3 == "CDS",]
    
    edta$V4 <- as.numeric(edta$V4)
    edta$V5 <- as.numeric(edta$V5)
    helixer$V4 <- as.numeric(helixer$V4)
    helixer$V5 <- as.numeric(helixer$V5)
    
    edta <- edta[!is.na(edta$V4),]
    edta <- edta[!is.na(edta$V5),]
    helixer <- helixer[!is.na(helixer$V4),]
    helixer <- helixer[!is.na(helixer$V4),]
    
    features_genome_vectors <- list()
    cat(nrow(cen_coords))
    for(j in 1 : nrow(cen_coords)) {
      cat("", j)
      chromosome <- cen_coords$chromosome[j]
      chromosome_size <- cen_coords$size[j]
      
      repeats_chr <- repeats[repeats$seqID == chromosome,]
      edta_chr <- edta[edta$V1 == chromosome,]
      helixer_chr <- helixer[helixer$V1 == chromosome,]
      
      if(nrow(repeats_chr) == 0) next
      if(nrow(edta_chr) != 0) {
        edta_gr    <- GRanges("chr1", IRanges(edta_chr$V4,       edta_chr$V5))
      } else {
        edta_gr    <- GRanges()
      }
      if(nrow(helixer_chr) != 0) {
        helixer_gr <- GRanges("chr1", IRanges(helixer_chr$V4,    helixer_chr$V5))
      } else {
        edta_gr    <- GRanges()
      }
      repeats_gr <- GRanges("chr1", IRanges(repeats_chr$start, repeats_chr$end))
      
      
      
      # calculate bin positions
      upstream_bin_starts <- rev(cen_coords$centromere_start[j] - (1:bins_flanks) * bins_flank_size)
      upstream_bin_ends   <- upstream_bin_starts + bins_flank_size - 1
      
      upstream_bin_starts[upstream_bin_starts < 1] <- 0
      upstream_bin_ends[upstream_bin_ends <= 1]    <- 1
      upstream_bin_ends[upstream_bin_ends <= 1]    <- 1
      
      centromere_length <- cen_coords$centromere_end[j] - cen_coords$centromere_start[j] + 1
      array_bin_size <- centromere_length / 100
      
      array_bin_starts      <- cen_coords$centromere_start[j] + floor((0:99)  * array_bin_size)
      array_bin_ends        <- cen_coords$centromere_start[j] + floor((1:100) * array_bin_size) - 1
      array_bin_ends[100]   <- cen_coords$centromere_end[j]
      
      downstream_bin_starts <- cen_coords$centromere_end[j] + (1:bins_flanks) * bins_flank_size
      downstream_bin_ends   <- downstream_bin_starts + bins_flank_size - 1
      
      downstream_bin_starts[downstream_bin_starts >= chromosome_size] <- chromosome_size - 1
      downstream_bin_ends[downstream_bin_ends >= chromosome_size] <- chromosome_size
      
      # Key change 5: concatenate all bins → single calc_density call per feature type
      all_starts <- c(upstream_bin_starts, array_bin_starts, downstream_bin_starts)
      all_ends   <- c(upstream_bin_ends,   array_bin_ends,   downstream_bin_ends)
      
      te_density     <- calc_density(edta_gr,    all_starts, all_ends)
      gene_density   <- calc_density(helixer_gr, all_starts, all_ends)
      repeat_density <- calc_density(repeats_gr, all_starts, all_ends)
      
      te_density[all_ends <= 1] = NA
      gene_density[all_ends <= 1] = NA
      repeat_density[all_ends <= 1] = NA
      
      te_density[all_starts >= chromosome_size] = NA
      gene_density[all_starts >= chromosome_size] = NA
      repeat_density[all_starts >= chromosome_size] = NA
      
      # # TODO: HOR_score
      # HOR_score[1 : bins_flanks]
      # HOR_score[(bins_flanks + 1) : (bins_flanks + 100)]
      # HOR_score[(bins_flanks + 101) : (bins_flanks + 200)]
      
      
      # randomize order to reduce noise 
      if (sample(0:1, 1) == 1) reverse = TRUE else reverse = FALSE
      features_genome_vectors[[as.character(j)]]$te_density <- if (reverse) te_density else rev(te_density)
      features_genome_vectors[[as.character(j)]]$gene_density <- if (reverse) gene_density else rev(gene_density)
      features_genome_vectors[[as.character(j)]]$repeat_density <- if (reverse) repeat_density else rev(repeat_density)
      # features_genome_vectors[[as.character(k)]]$HOR_score <- if (reverse) HOR_score else rev(HOR_score)
      
    }
    cat("\n")
    
    te_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$te_density)), na.rm = TRUE)
    gene_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$gene_density)), na.rm = TRUE)
    repeat_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$repeat_density)), na.rm = TRUE)
    
    genomes_means_vectors[[species_name]]$te <- te_per_bin_mean
    genomes_means_vectors[[species_name]]$gene <- gene_per_bin_mean
    genomes_means_vectors[[species_name]]$reps <- repeat_per_bin_mean
    
    cen_coords$cen_size <- cen_coords$centromere_end - cen_coords$centromere_start + 1
    mean_arrays_width <- mean(cen_coords$cen_size)
    # mean_arrays_width <- mean_arrays_width/10
    if(mean_arrays_width > 5000000) mean_arrays_width <- 5000000
    
    setwd("/home/pwlodzimierz/ToL/upload_files/33.1_metaplots_metacentrics")
    pdf(file = paste0("single_plot_", species_name, ".pdf")) 
    plot(NULL,NULL, 
         xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size), 
         ylim = c(0,100),
         main = species_name)
    
    plot_bin_starts <- c(0: (bins_flanks-1) * bins_flank_size + 1, 
                         round(bins_flanks * bins_flank_size + 0:(bins_flanks-1) * (mean_arrays_width/bins_flanks)), 
                         round(bins_flanks * bins_flank_size + (bins_flanks-1) * (mean_arrays_width/bins_flanks)) + 0:(bins_flanks-1) * bins_flank_size + 1)
    plot_bin_ends <- c(plot_bin_starts[2:length(plot_bin_starts)] - 1, 
                       plot_bin_starts[length(plot_bin_starts)] + bins_flank_size - 1)
    plot_bin_mids <- unlist(lapply(1 : (3*bins_flanks), 
                                   function(X) round(mean(c(plot_bin_ends[X], plot_bin_starts[X])))))
    
    lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
    lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
    lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
    abline(v = plot_bin_starts[bins_flanks + 1])
    abline(v = plot_bin_ends[2 * bins_flanks])
    dev.off()
    
  # })
  
  mean_centromeres_sizes <- c(mean_centromeres_sizes, mean(cen_coords$cen_size))
}

save(genomes_means_vectors, 
     file = "/home/pwlodzimierz/ToL/upload_files/33.1_metaplots_metacentrics/genomes_means_vectors.rds")
load("/home/pwlodzimierz/ToL/upload_files/33.1_metaplots_metacentrics/genomes_means_vectors.rds")

for(i in 1 : length(species)) {
 
  species_name <- strsplit(species[i], split = "[.]")[[1]][1]
  species_sats <- satellites[satellites$fasta == species[i],]
  
  if(species_sats$is_holocentric[1]) {
    cat(i, species_name, "holocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
    next
  }
  cat(i, species_name, "monocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
  
  if(!dir.exists(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", species_name))) next
  
  setwd(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", species_name))
  reps_files <- list.files(path = ".", pattern = "HOR_scored_repeats_with", full.names = T)
  
  
  if(species_name == "rosCan_S27_v1") {
    cen_coords <- read.csv(file = "/home/pwlodzimierz/ToL/upload_files/55_cen_coords/rosCan_genome_metadata.csv")
  } else {
    cen_coords <- read.csv(file = cen_coords_files[grep(species_name, cen_coords_files)])
  }
  
  features_genome_vectors <- list()
  cat(nrow(cen_coords))
  for(j in 1 : nrow(cen_coords)) {
    cat("", j)
    chromosome <- cen_coords$chromosome[j]
    chromosome_size <- cen_coords$size[j]
    
    repeats_files <- reps_files[grep(paste0(chromosome, ".csv"), reps_files, fixed = T)]
    if(length(repeats_files) == 0) next
    
    repeats_chr <- data.frame()
    for(k in 1 : length(repeats_files)) {
      if(file.size(repeats_files[k]) < 10) next
      repeats_chr <- rbind(repeats_chr, read.csv(file = repeats_files[k]))
    }
    
    if(nrow(repeats_chr) == 0) next
    
    upstream_bin_starts <- rev(cen_coords$centromere_start[j] - (1:bins_flanks) * bins_flank_size)
    upstream_bin_ends   <- upstream_bin_starts + bins_flank_size - 1
    
    upstream_bin_starts[upstream_bin_starts < 1] <- 0
    upstream_bin_ends[upstream_bin_ends <= 1]    <- 1
    upstream_bin_ends[upstream_bin_ends <= 1]    <- 1
    
    centromere_length <- cen_coords$centromere_end[j] - cen_coords$centromere_start[j] + 1
    array_bin_size <- centromere_length / 100
    
    array_bin_starts      <- cen_coords$centromere_start[j] + floor((0:99)  * array_bin_size)
    array_bin_ends        <- cen_coords$centromere_start[j] + floor((1:100) * array_bin_size) - 1
    array_bin_ends[100]   <- cen_coords$centromere_end[j]
    
    downstream_bin_starts <- cen_coords$centromere_end[j] + (1:bins_flanks) * bins_flank_size
    downstream_bin_ends   <- downstream_bin_starts + bins_flank_size - 1
    
    downstream_bin_starts[downstream_bin_starts >= chromosome_size] <- chromosome_size - 1
    downstream_bin_ends[downstream_bin_ends >= chromosome_size] <- chromosome_size
    
    all_starts <- c(upstream_bin_starts, array_bin_starts, downstream_bin_starts)
    all_ends   <- c(upstream_bin_ends,   array_bin_ends,   downstream_bin_ends)
    
    
    bins_gr    <- GRanges("chr1", IRanges(all_starts, all_ends))
    starts_gr <- GRanges("chr1", IRanges(repeats_chr$start, width = 1))
    hits      <- findOverlaps(starts_gr, bins_gr, select = "first")
    
    repeats_chr$bin <- hits
    
    HOR_score_bins <- rep(NA, 300)
    mean_by_bin <- tapply(repeats_chr$HOR_score, repeats_chr$bin, mean, na.rm = TRUE)
    
    HOR_score_bins[as.integer(names(mean_by_bin))] <- mean_by_bin
    
    
    
    
    HOR_score_bins [all_ends <= 1] = NA
    HOR_score_bins [all_starts >= chromosome_size] = NA
    
    
    
    # randomize order to reduce noise 
    if (sample(0:1, 1) == 1) reverse = TRUE else reverse = FALSE
    features_genome_vectors[[as.character(j)]]$HOR_score_bins  <- if (reverse) HOR_score_bins  else rev(HOR_score_bins )
    
    
  }
  cat("\n")
  
  hor_per_bin_mean <- colMeans(do.call(rbind, lapply(features_genome_vectors, function(x) x$HOR_score)), na.rm = TRUE)
  
  genomes_means_vectors[[species_name]]$hor <- hor_per_bin_mean


}

save(genomes_means_vectors, 
     file = "/home/pwlodzimierz/ToL/upload_files/33.1_metaplots_metacentrics/genomes_means_vectors.rds")
load("/home/pwlodzimierz/ToL/upload_files/33.1_metaplots_metacentrics/genomes_means_vectors.rds")



# ### metaplot
# te_per_bin_mean <- vector(mode = "numeric", length = 300)
# gene_per_bin_mean <- vector(mode = "numeric", length = 300)
# repeat_per_bin_mean <- vector(mode = "numeric", length = 300)
# HOR_per_bin_mean <- vector(mode = "numeric", length = 300)
# for(i in 1 : length(genomes_means_vectors)) {
#   te_per_bin_mean <- te_per_bin_mean + genomes_means_vectors[[i]]$te
#   gene_per_bin_mean <- gene_per_bin_mean + genomes_means_vectors[[i]]$gene
#   repeat_per_bin_mean <- repeat_per_bin_mean + genomes_means_vectors[[i]]$reps
#   HOR_per_bin_mean <- HOR_per_bin_mean + genomes_means_vectors[[i]]$hor
# }
# te_per_bin_mean <- te_per_bin_mean / length(genomes_means_vectors)
# gene_per_bin_mean <- gene_per_bin_mean / length(genomes_means_vectors)
# repeat_per_bin_mean <- repeat_per_bin_mean / length(genomes_means_vectors)
# HOR_per_bin_mean <- HOR_per_bin_mean / length(genomes_means_vectors)
# 
# 
# 
# pdf(file = paste0("metaplot_monocentromeres.pdf")) 
# plot(NULL,NULL, 
#      xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size), 
#      ylim = c(0,100),
#      main = "all monocentrics")
# 
# plot_bin_starts <- c(0:99 * bins_flank_size + 1, round(100 * bins_flank_size + 0:99 * (mean_arrays_width/100)), round(100 * bins_flank_size + 99 * (mean_arrays_width/100)) + 0:99 * bins_flank_size + 1)
# plot_bin_ends <- c(plot_bin_starts[2:length(plot_bin_starts)] - 1, plot_bin_starts[length(plot_bin_starts)] + bins_flank_size - 1)
# plot_bin_mids <- unlist(lapply(1 : 300, function(X) round(mean(c(plot_bin_ends[X], plot_bin_starts[X])))))
# 
# lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
# lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
# lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
# lines(x = plot_bin_mids, y = HOR_per_bin_mean, col = "orange")
# abline(v = plot_bin_starts[101])
# abline(v = plot_bin_ends[200])
# dev.off()


mean_arrays_width = 2500000


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

HOR_scaling = 5

pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/33.1_metaplots_metacentrics/metaplot_monocentromeres_2.5_HOR_scaling_5.pdf"))
plot(NULL, NULL,
     xlim = c(0, bins_flanks*bins_flank_size + mean_arrays_width + bins_flanks*bins_flank_size),
     ylim = c(0, 100),
     main = "all monocentrics")

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
draw_ci_ribbon(plot_bin_mids, HOR_per_bin_mean * HOR_scaling, HOR_ci * HOR_scaling, "orange")

lines(x = plot_bin_mids, y = te_per_bin_mean, col = "blue")
lines(x = plot_bin_mids, y = gene_per_bin_mean, col = "green")
lines(x = plot_bin_mids, y = repeat_per_bin_mean, col = "red")
lines(x = plot_bin_mids, y = HOR_per_bin_mean * HOR_scaling, col = "orange")

abline(v = plot_bin_starts[101])
abline(v = plot_bin_ends[200])
dev.off()
