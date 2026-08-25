
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))



### functions

{
  
  
  plot_hist <- function(dist_data, plot_name = "histogram", min_dist = 2, max_dist = 1000) {
    par(mar = c(3,3,1,3))
    
    dist_data <- dist_data[dist_data >= min_dist]
    dist_data <- dist_data[dist_data <= max_dist]
    
    hist(dist_data, breaks = min_dist : max_dist, xlim = c(min_dist, max_dist), main = plot_name)
    mtext("Distance between HOR blocks", side = 1, line = 2, cex = 0.5)
    mtext("Frequency", side = 2, line = 2)
    
    
  }
  
  
  # -----------------------------------------------------------------------------
  # 1. Native-resolution counts, edge-corrected for the N-d combinatorial ceiling
  #    Raw count at distance d is capped by (N - d) possible pairs -- dividing
  #    by this removes the purely combinatorial decay so what's left reflects
  #    biological structure, not array-length artifacts.
  # -----------------------------------------------------------------------------
  build_count_signal <- function(dist_vec, rep_no, max_d = NULL) {
    max_d <- if (is.null(max_d)) rep_no - 1 else max_d
    counts <- tabulate(dist_vec, nbins = max_d)          # exact bin = 1 repeat
    d <- seq_len(max_d)
    possible_pairs <- pmax(rep_no - d, 1)                # avoid divide-by-zero
    data.frame(d = d, rate = counts / possible_pairs)
  }
  
  # -----------------------------------------------------------------------------
  # 2. Periodogram on the edge-corrected signal
  #    detrend removes any residual smooth trend (e.g. gradual similarity
  #    decay along the array) so remaining power reflects oscillation only
  # -----------------------------------------------------------------------------
  run_periodogram <- function(sig_df) {
    spec <- spec.pgram(sig_df$rate, taper = 0.1, detrend = TRUE,
                       demean = TRUE, plot = FALSE)
    list(freq = spec$freq, power = spec$spec)   # freq already in cycles/repeat
  }
  
  # -----------------------------------------------------------------------------
  # 3. Spectral entropy (0 = concentrated/periodic, 1 = spread/heterogeneous)
  # -----------------------------------------------------------------------------
  spectral_entropy <- function(power) {
    p <- power[power > 0] / sum(power)
    -sum(p * log(p)) / log(length(p))
  }
  
  # -----------------------------------------------------------------------------
  # 4. Peak prominence vs local background, restricted to a biologically
  #    plausible period range [min_period, max_period]. Frequencies outside
  #    this range (e.g. the very-low-frequency artifact peaks seen in
  #    heterogeneous arrays) are excluded from the search entirely.
  #
  #    The fundamental frequency is chosen via a Harmonic Product Spectrum
  #    (HPS) rather than a simple global max. A comb-like periodic signal has
  #    real power not just at the fundamental but at all its harmonics
  #    (2x, 3x, 4x...), and with limited data noise can make a harmonic bin
  #    momentarily outrank the true fundamental. HPS instead scores each
  #    candidate fundamental by its support across its *whole* harmonic
  #    series -- the true fundamental's series is consistently strong, while
  #    a false (harmonic) candidate only borrows a subset of it -- which
  #    reliably resolves this "octave ambiguity".
  # -----------------------------------------------------------------------------
  harmonic_product_freq <- function(freq, power, candidates, n_harmonics = 5) {
    max_f <- max(freq)
    h <- n_harmonics
    valid <- candidates[candidates * h <= max_f]
    while (length(valid) == 0 && h > 1) {        # relax only if nothing qualifies
      h <- h - 1
      valid <- candidates[candidates * h <= max_f]
    }
    if (length(valid) == 0) valid <- candidates  # last resort: single-bin comparison
    
    scores <- vapply(valid, function(f0) {
      harm <- f0 * seq_len(h)                    # every candidate scored on the SAME h terms
      amp <- approx(freq, power, xout = harm, rule = 2)$y
      mean(log(amp + 1e-300))                    # mean, not sum -- avoids penalizing/rewarding by term count
    }, numeric(1))
    valid[which.max(scores)]
  }
  
  peak_prominence <- function(freq, power, min_period = 2, max_period = 100,
                              n_harmonics = 5) {
    period <- 1 / freq
    in_range <- period >= min_period & period <= max_period
    f <- freq[in_range]; p <- power[in_range]
    if (length(f) < 2) {
      return(list(
        peak_freq = NA_real_,
        peak_power = NA_real_,
        prominence = NA_real_,
        in_range_power = numeric(0)
      ))
    }
    f0 <- harmonic_product_freq(freq, power, f, n_harmonics)
    i <- which.min(abs(f - f0))
    excl <- max(1, i - 3):min(length(p), i + 3)
    list(peak_freq = f[i], peak_power = p[i],
         prominence = p[i] / median(p[-excl]), in_range_power = p)
  }
  
  # -----------------------------------------------------------------------------
  # 4b. Fisher's g-test: exact p-value for the strongest periodogram ordinate
  #    under the null of no periodicity. Unlike a fixed prominence_cutoff, this
  #    threshold is self-calibrating -- it accounts for how many frequencies
  #    (m) were tested, so it's valid whether the window has 5 cycles or 500.
  #    g = max(power) / sum(power); reference: Fisher (1929) exact test.
  # -----------------------------------------------------------------------------
  # fisher_g_pvalue <- function(power) {
  #   m <- length(power)
  #   g <- max(power) / sum(power)
  #   K <- floor(1 / g)
  #   k <- seq_len(min(K, m))
  #   terms <- (-1)^(k - 1) * choose(m, k) * pmax(1 - k * g, 0)^(m - 1)
  #   min(max(sum(terms), 0), 1)
  # }
  
  fisher_g_pvalue <- function(power) {
    
    power <- power[is.finite(power)]
    if (length(power) < 2)
      return(NA_real_)
    
    s <- sum(power)
    if (s <= 0)
      return(NA_real_)
    
    m <- length(power)
    g <- max(power) / s
    if (!is.finite(g) || g <= 0)
      return(NA_real_)
    
    K <- floor(1 / g)
    if (K < 1)
      return(1)
    
    k <- seq_len(min(K, m))
    terms <- (-1)^(k - 1) *
      choose(m, k) *
      pmax(1 - k * g, 0)^(m - 1)
    min(max(sum(terms), 0), 1)
  }
  
  # -----------------------------------------------------------------------------
  # 5. Master function
  #    Note: p_value (fisher_g_pvalue) is computed on the raw in-range power
  #    vector, i.e. it still tests the global max ordinate -- this is an
  #    omnibus "is anything significant here" test and is unaffected by which
  #    frequency HPS ultimately reports as the fundamental. estimated_period
  #    is what HPS corrects; p_value/call answer a different question
  #    ("periodic at all?") and remain valid regardless.
  # -----------------------------------------------------------------------------
  analyze_periodicity <- function(dist_vec, rep_no, min_period = 2, max_period = 100) {
    sig  <- build_count_signal(dist_vec, rep_no)
    spec <- run_periodogram(sig)
    prom <- peak_prominence(spec$freq, spec$power, min_period, max_period)
    data.frame(
      estimated_period = 1 / prom$peak_freq,   # in repeat units
      peak_power       = prom$peak_power,
      spectral_entropy = spectral_entropy(spec$power),
      peak_prominence  = prom$prominence,
      p_value          = fisher_g_pvalue(prom$in_range_power)
    )
  }
  
  # -----------------------------------------------------------------------------
  # 6. Classifier: is there a significant peak within the plausible period
  #    range? p_value from Fisher's g-test replaces the old fixed
  #    prominence_cutoff -- self-calibrated to window/sample size, so the same
  #    alpha (e.g. 0.01) is valid across very different window sizes.
  # -----------------------------------------------------------------------------
  classify_periodicity <- function(dist_vec, rep_no, min_period = 2, max_period = 100,
                                   alpha = 0.001) {
    res <- analyze_periodicity(dist_vec, rep_no, min_period, max_period)
    res$call <- "heterogeneous"
    res$call[!is.na(res$p_value) & res$p_value <= alpha] <- "periodic"
    res
  }
  
  # -----------------------------------------------------------------------------
  # 7. Sliding-window localization along the array (by repeat ID, using
  #    start_A/start_B -- no genomic coordinates needed).
  #    Each window is treated as its own mini-array: local N = window_size,
  #    so the same edge correction/periodogram logic applies unchanged, just
  #    restricted to pairs whose start_A falls inside the window.
  # -----------------------------------------------------------------------------
  window_track <- function(hors, rep_no, window_size = 60, step = 30,
                           min_period = 2, max_period = 100,
                           alpha = 0.001, min_pairs = 20) {
    optimal_window_no <- round(rep_no / window_size)
    optimal_window_no <- max(optimal_window_no, 1)
    window_size <- rep_no / optimal_window_no
    if(window_size >= rep_no) {
      starts <- 1 
    } else {
      starts <- seq(1, rep_no - window_size, by = step)
    }
    
    if(starts[length(starts)] < rep_no - window_size) starts <- c(starts, rep_no - window_size)
    out <- lapply(starts, function(w) {
      sub <- hors[hors$start_A >= w & hors$start_A < w + window_size, ]
      if (nrow(sub) < min_pairs) {
        return(data.frame(window_start = w, window_end = w + window_size,
                          estimated_period = NA, peak_prominence = NA,
                          p_value = NA, call = "insufficient_data"))
      }
      res <- classify_periodicity(sub$dist, window_size, min_period,
                                  max_period, alpha)
      cbind(window_start = w, window_end = w + window_size,
            res[, c("estimated_period", "peak_prominence", "p_value", "call")])
    })
    do.call(rbind, out)
  }
  
  # -----------------------------------------------------------------------------
  # 8. Collapse the window track into contiguous regions.
  #
  #    Step 1 -- smooth flickers: a short run of "heterogeneous" calls flanked
  #    on both sides by "periodic" windows reporting the SAME period is treated
  #    as a local significance dip (e.g. weaker local signal, not a real
  #    biological gap) and reclassified as periodic. Flanks must agree on
  #    period -- a run flanked by two DIFFERENT periods is left alone, since
  #    that's exactly the signature of a genuine boundary between two
  #    differently-periodic regions, not noise.
  #
  #    Step 2 -- merge: adjacent windows are grouped into contiguous regions
  #    whenever they agree on call (and, for periodic windows, on period
  #    rounded to the nearest integer) -- this is what lets two directly-
  #    adjacent regions with different periods, with no gap, still be split.
  # -----------------------------------------------------------------------------
  smooth_flickers <- function(track, max_flicker = 3) {
    calls <- track$call; periods <- track$estimated_period
    r <- rle(calls == "heterogeneous")
    idx <- 1
    for (i in seq_along(r$lengths)) {
      len <- r$lengths[i]; rng <- idx:(idx + len - 1)
      if (r$values[i] && len <= max_flicker && idx > 1 && idx + len <= nrow(track)) {
        before <- calls[idx - 1]; after <- calls[idx + len]
        bp <- periods[idx - 1]; ap <- periods[idx + len]
        if (before == "periodic" && after == "periodic" &&
            !is.na(bp) && !is.na(ap) && bp == ap) {
          calls[rng] <- "periodic"; periods[rng] <- bp
        }
      }
      idx <- idx + len
    }
    track$call <- calls; track$estimated_period <- periods
    track
  }
  
  summarize_regions <- function(track, max_flicker = 3) {
    track <- smooth_flickers(track, max_flicker)
    key <- ifelse(track$call == "periodic", round(track$estimated_period), NA)
    grp <- cumsum(track$call != c(track$call[1], head(track$call, -1)) |
                    (!is.na(key) & key != c(key[1], head(key, -1))))
    agg <- split(track, grp)
    agg <- do.call(rbind, lapply(agg, function(r) data.frame(
      start = min(r$window_start), end = max(r$window_end),
      call = r$call[1],
      estimated_period = if (r$call[1] == "periodic") median(r$estimated_period) else NA
    )))
    for(i in seq_len(nrow(agg))) {
      if(i == nrow(agg)) break
      if(agg$end[i] > agg$start[i + 1]) {
        mid <- floor(mean(c(agg$end[i], agg$start[i + 1])))
        
        agg$end[i] <- mid
        agg$start[i + 1] <- mid + 1
        
      }
    }
    agg
  }
  
  summarise_chromosome <- function(regions) {
    regions$estimated_period[is.na(regions$estimated_period)] = 0
    regions$estimated_period <- round(regions$estimated_period)
    regions$width <- regions$end - regions$start + 1
    periods <- table(regions$estimated_period)
    if(length(periods) == 0) return(data.frame(period = 0, repeats = sum(regions$width), fraction = 100))
    
    periods <- as.data.frame(periods)
    names(periods) <- c("period", "repeats")
    periods$period <- as.numeric(as.character(periods$period))
    periods$repeats <- 0
    
    for(i in seq_len(nrow(periods))) {
      periods$repeats[i] <- sum(regions$width[regions$estimated_period == periods$period[i]])
    }
    
    sum_repeats <- sum(regions$width)
    periods$fraction <- round(100 * periods$repeats / sum_repeats, 2)
    
    periods
  }
  
  # -----------------------------------------------------------------------------
  # Example usage:
  # result <- classify_periodicity(hors$dist, rep_no, max_period = 100)
  # print(result)
  #
  # track   <- window_track(hors, rep_no, window_size = 60, step = 30)
  # regions <- summarize_regions(track)
  # print(regions)   # start/end (repeat IDs), call, estimated_period per region
  #
  # Notes:
  #  - window_size should comfortably span several periods (e.g. >= 5x the
  #    expected period) for a stable local estimate; step controls boundary
  #    resolution (smaller step = finer localization, more compute).
  #  - Region boundaries are only as precise as window_size/step; for exact
  #    breakpoints, follow up with a finer-resolution scan near a detected
  #    transition.
  
  periodicity <- function(hors, rep_no) {
    result <- classify_periodicity(hors$dist, rep_no, max_period = 100)
    if (is.na(result$estimated_period)) {
      window_size <- 1000
    } else if (result$estimated_period < 50) {
      window_size <- round(result$estimated_period * 20)
    } else {
      window_size <- 1000
    }
    
    track   <- window_track(hors, rep_no, window_size = window_size, step = window_size / 4)
    regions <- summarize_regions(track)
    summary <- summarise_chromosome(regions)
    print(result)
    print(summary)
    list(result,track,regions,summary)
  }
  
  
  regimentation_plotter <- function(repeats, plotting_data, plot_name) {
    par(mar = c(3,3,1,3))
    plot(x = repeats$start, y = repeats$HOR_score, ylim = c(0,100), type = "h", main = plot_name,
         ylab = "HOR score")
    axis(4, at = seq(0, 100, 10))   
    mtext("Coordinates (bp)", side = 1, line = 2, cex = 0.5)
    mtext("HOR score", side = 2, line = 2)
    mtext("HOR period", side = 4, line = 2)
    
    gr_rep <- GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = repeats$start, 
                       end = repeats$end)
    )
    gr_rep <- reduce(gr_rep)
    for(i in seq_along(gr_rep)) {
      colline = "#0000ee"
      
      lines(x = c(start(gr_rep)[i], end(gr_rep)[i]), 
            y = rep(-2, 2), 
            col = colline,
            lwd = 4,
            lend = 1)
    }
    
    regions <- plotting_data[[3]]
    regions$estimated_period[is.na(regions$estimated_period)] = 0
    
    gr_regions <- GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = repeats$start[regions$start], 
                       end = repeats$end[regions$end - 1]),
      estimated_period = regions$estimated_period,
      call = regions$call
    )
    
    hits <- findOverlaps(gr_regions, gr_rep)
    
    gr_regions_split <- pintersect(gr_regions[queryHits(hits)], gr_rep[subjectHits(hits)])
    mcols(gr_regions_split) <- mcols(gr_regions)[queryHits(hits), , drop = FALSE]
    
    for(i in seq_along(gr_regions_split)) {
      colline = "#eeee4499"
      if(mcols(gr_regions_split)$call[i] == "periodic") colline = "#eeee44dd"
      
      lines(x = c(start(gr_regions_split)[i], end(gr_regions_split)[i]), 
            y = rep(mcols(gr_regions_split)$estimated_period[i], 2), 
            col = colline,
            lwd = 4,
            lend = 1)
    }
    
    gr_regions_split
    
  }
  
  
  
}



repeats_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/2026_filtered_TRASH", full.names = T)

satellites <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")
metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv")
phyla <- read.csv("/home/pwlodzimierz/ToL/Metadata/fasta_phyla.csv")

metadata <- metadata[metadata$is.chr == 1,]

species <- unique(satellites$fasta)


taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 166
print(i)


species_name <- strsplit(species[i], split = "[.]")[[1]][1]
species_sats <- satellites[satellites$fasta == species[i],]

if(species_sats$is_holocentric[1]) {
  cat(i, species_name, "holocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")
  next
}
cat(i, species_name, "monocentric \n")#, file = "/home/pwlodzimierz/ToL/git_ToL/out.txt")


chromosomes_metadata <- metadata[grep(species_name, metadata$assembly.name, fixed = T),]
chromosomes <- chromosomes_metadata$chromosome.name
if(i == 166) chromosomes <-  c("CM000663.2", "CM000664.2", "CM000665.2", "CM000666.2", "CM000667.2",
                               "CM000668.2", "CM000669.2", "CM000670.2", "CM000671.2", "CM000672.2",
                               "CM000673.2", "CM000674.2", "CM000675.2", "CM000676.2", "CM000677.2",
                               "CM000678.2", "CM000679.2", "CM000680.2", "CM000681.2", "CM000682.2",
                               "CM000683.2", "CM000684.2", "CM000685.2", "CM000686.2")

dir.create(path = paste0("/home/pwlodzimierz/ToL/upload_files/67_regimentation/data/", species_name), showWarnings = F)
for(j in 1 : nrow(species_sats)) {
  
  setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out/", species_name))
  hor_files <- list.files(path = ".", pattern = "HORs")
  hor_files <- hor_files[grep(".csv", hor_files, fixed = T)]
  hor_files <- hor_files[grep(species_sats$TRASH_new_class[j], hor_files, fixed = T)]
  
  if(length(hor_files) == 0) next
  
  setwd(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", species_name))
  repeat_files <- list.files(path = ".")
  repeat_files <- repeat_files[grep(species_sats$TRASH_new_class[j], repeat_files, fixed = T)]
  
  if(length(repeat_files) == 0) next
  
  
  
  # count chromosomes that can be plotted
  chr_with_data <- 0
  for(k in seq_along(chromosomes)) {
    hor_file <- hor_files[grep(paste0(chromosomes[k], ".csv"), hor_files, fixed = T)]
    repeat_file <- repeat_files[grep(paste0(chromosomes[k], ".csv"), repeat_files, fixed = T)]
    if(length(hor_file) != 1) next
    if(length(repeat_file) != 1) next
    if(file.size(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out/", species_name, "/", hor_file)) < 10) next
    if(file.size(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", species_name, "/", repeat_file)) < 10) next
    chr_with_data <- chr_with_data + 1
    
  }
  
  setwd("/home/pwlodzimierz/ToL/upload_files/67_regimentation/species")
  # 
  # if(file.exists(paste0("regimentation_per_chr_", species_name, "_", species_sats$TRASH_new_class[j], ".pdf"))) {
  #   print("File already done, skipping")
  #   next
  # }
  
  pdf(file = paste0("regimentation_per_chr_", species_name, "_", species_sats$TRASH_new_class[j], ".pdf"), width = 12, height = 3 * chr_with_data)
  par(mfrow = c(2*chr_with_data, 1))
  
  dir.create(path = paste0("/home/pwlodzimierz/ToL/upload_files/67_regimentation/data/", species_name, "/", species_sats$TRASH_new_class[j]), showWarnings = F)
  
  
  
  for(k in seq_along(chromosomes)) {
    try({
      cat(species_name, species_sats$TRASH_new_class[j], k, "/", chr_with_data, " chr with no data:", length(chromosomes) - chr_with_data, "\n")
      
      hor_file <- hor_files[grep(paste0(chromosomes[k], ".csv"), hor_files, fixed = T)]
      repeat_file <- repeat_files[grep(paste0(chromosomes[k], ".csv"), repeat_files, fixed = T)]
      
      if(length(hor_file) != 1) next
      if(length(repeat_file) != 1) next
      
      if(file.size(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out/", species_name, "/", hor_file)) < 10) next
      if(file.size(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", species_name, "/", repeat_file)) < 10) next
      
      hors <- read.csv(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out/", species_name, "/", hor_file))
      repeats <- read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs/", species_name, "/", repeat_file))
      
      if(nrow(hors) < 10) next
      if(nrow(repeats) < 10) next
      
      rep_no <- nrow(repeats)
      hors$dist <- hors$start_B - hors$start_A
      
      a <- periodicity(hors, rep_no)
      
      a[["chromosome"]] <- chromosomes[k]
      a[["genome"]] <- species_name
      a[["sattelite"]] <- species_sats$TRASH_new_class[j]
      a[["repeats"]] <- repeats[,c("start", "HOR_score")]
      
      a[[2]]$window_mid_bp <- repeats$start[round(a[[2]]$window_start + (a[[2]]$window_end - a[[2]]$window_start)/2)]
      
      save(a, file = paste0("/home/pwlodzimierz/ToL/upload_files/67_regimentation/data/", 
                            species_name, "/", species_sats$TRASH_new_class[j], "/", chromosomes[k],
                            ".rds"))
      
      plot_hist(hors$dist, paste0(chromosomes[k], ": ", a[[1]]$call, ", periodic fraction: ", 100 - a[[4]]$fraction[a[[4]]$period == 0]))
      split_regions <- regimentation_plotter(repeats, a, "")
     
    })
    
    
  }
  
  dev.off()
  
  
}





