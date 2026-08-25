


#############

# # 
# # Satellite-species
# # a. Centromeres include centromeric repeats (check histogram of array lengths to decide a minimum size in a 
# 1-25 kb range) and any gaps of up to 100 kb (50 kb for holocentrics), single region is called a centromeric array 
# and more than one centromeric array can be annotated per chromosome, together they from a centromeric region. Gaps 
# longer than 100 kbp that divide two arrays might be included in the centromere region manually, if they are 
# occupied by a high density centrophilic (as defined within the remaining gaps) elements, so far no gaps like this 
# were found.
# # b. Pericentromeres include non-centromeric regions that are characterised as abundant in LTRs, compared to 
# the rest of the chromosome, as calculated using window-based TE density scoring divided by Otsu method. In some 
# cases when LTR distribution is mostly flat and automatic division struggles to make a clear choice, pericentromeres 
# are annotated as regions adjacent to centromeric arrays, with the size being equal to 10% of the chromosome size 
# (5% on each side)
# # c. Arms: anything that's left
# # 
# # TE-species
# #   a. Centromeres include distinguished and uniform regions of a specific TE superfamily enrichment, annotated 
# manually only for the most prominent cases (Geum, Malus etc). Remaining species have their centromere unannotated. 
# Centromeric gaps are regions between centromeric TEs.
# #   b. Pericentromeres include non-centromeric regions that are characterised as abundant in LTRs, compared to the 
# rest of the chromosome, as calculated using window-based TE density scoring divided by Otsu method. In some cases 
# when LTR distribution is mostly flat and automatic division struggles to make a clear choice, pericentromeres are 
# annotated as regions adjacent to the chromosome-wide LTR peak (most likely centromere locus), with the size being 
# equal to 10% of the chromosome size (5% on each side)
# #   c. Arms: anything that's left
# # 
# # Holocentric species
# # a. Centromeres include centromeric repeats and any gaps of up to 50 kb, single region is called a centromeric 
# island and all islands together from a centromeric region.
# # b. Pericentromeres are annotated as regions adjacent to centromeric islands, with the size being equal to the 
# size of the centromeric array (50% on ech side)
# # c. Arms: anything that's left
# # 

gap_TE_coverage_count <- 1

# find CENTROMERIC arrays, keep those where gaps are not more than:

max_gap_distance <- 100000 # base pairs for monocentrics and
max_gap_dist_holo <- 50000 # base pairs for holocentric


# within arrays, consider a gap of at least
min_internal_gap_distance <- 250 # base pairs


te_score_calculation_window = 100000 # the window is as close as possible to this value, 
# while ensuring each window is the same size




cen_TE_grep_pattern <- "(?i)^(?!.*non[-]?LTR).*LTR"

plotting_window = 100000
average_plots <- TRUE


##############


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}


if(FALSE) { # DEV SETTINGS, SKIP TO "helixer_file" line 138 
  i = 24
  setwd("C:/Users/Piotr Włodzimierz/Desktop/ToL/temp_data/ddAraThal4.1.fa")
  assembly_name <- "ddAraThal4.1"
  fasta_name <- "ddAraThal4.1.fa"
  source("../aux_fun.R")
  
  genome_metadata <- read.csv("../../chr.no.and.sizes.full.Ian.csv")
  genome_metadata <- genome_metadata[genome_metadata$assembly.name == fasta_name, ]
  chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr == 1]
  chromosomes_lengths <- genome_metadata$size[genome_metadata$is.chr == 1]
  non_chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr != 1]
  non_chr_lengths <- genome_metadata$size[genome_metadata$is.chr != 1]
  
  # satellite_organisation_metadata <- read.csv(file = "../../satellite_organisation.csv")
  satellite_metadata <- read.csv("../../curated_satellites_metadata_on_chromosomes_only_may.csv")
  
  fasta <- read.fasta(file = "C:/Users/Piotr Włodzimierz/Desktop/ToL/assemblies/ddAraThal4.1.fa")
}

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")


data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2025", recursive = FALSE, full.names = TRUE))
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]


# satellite_organisation_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/satellite_organisation.csv")

satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
fasta_name = strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2]
fasta_name_short = strsplit(fasta_name, split = ".fa")[[1]][1]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])
fasta <- read.fasta(file = assembly_files[assembly_file])

setwd(data_directories[i])
print(getwd())

if(file.exists(paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/", fasta_name, "_genome_metadata.csv"))) {
  print("done before")
  quit()
}
  


genome_metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")
genome_metadata <- genome_metadata[genome_metadata$assembly.name == fasta_name, ]
chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr == 1]
chromosomes_lengths <- genome_metadata$size[genome_metadata$is.chr == 1]
non_chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr != 1]
non_chr_lengths <- genome_metadata$size[genome_metadata$is.chr != 1]


helixer_file = list.files(pattern = "helixer_filtered.gff", full.names = TRUE)
if(length(helixer_file) != 1) {warning(paste0(i, "No genes!")); setwd(".."); quit(save = "no", status = 1)}

repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}

# array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
# if(length(array_file) != 1) {warning(paste0(i, " no arrays!")); setwd(".."); quit(save = "no", status = 1)}

classes_file = list.files(pattern = "_classes_merged_filtered", full.names = TRUE)
if(length(classes_file) != 1) {warning(paste0(i, " no classes!")); setwd(".."); quit(save = "no", status = 1)}


edta_file = list.files(pattern = paste0(assembly_name, "_edta_filtered.csv"), full.names = TRUE)
edta_file <- edta_file[!grepl("reassigned", edta_file)]
if(length(edta_file) != 1) {warning(paste0(i, " no edta!")); setwd(".."); quit(save = "no", status = 1)}


print("load annotations")
repeats = read.csv(file = repeat_file)
# arrays = read.csv(file = array_file)
classes = read.csv(file = classes_file)
classes$num_ID <- 1 : nrow(classes)
edta = read.csv(file = edta_file, header = FALSE)
genes <- read.table(file = helixer_file, header = FALSE, sep = "\t", skip = 4)
genes <- genes[genes$V3 == "CDS", ]
print("annotations loaded")

names(edta)
names(genes)
print(nrow(genes))

edta$V5 <- as.numeric(edta$V5)
edta$V6 <- as.numeric(edta$V6)



empty_array <- data.frame(seqID = vector(mode = "character", length = 0), 
                          start = vector(mode = "numeric"), 
                          end = vector(mode = "numeric"), 
                          strand = vector(mode = "character"), 
                          strand_percentage = vector(mode = "numeric"), 
                          width = vector(mode = "numeric"),
                          satellite_name = vector(mode = "character"), 
                          TRASH_class = vector(mode = "character"), 
                          repeats_number = vector(mode = "numeric"), 
                          mean_repeat_width = vector(mode = "numeric"), 
                          mean_repeat_width_SD = vector(mode = "numeric"),
                          mean_repeat_pairwise_similarity = vector(mode = "numeric"), 
                          mean_repeat_pairwise_similarity_SD = vector(mode = "numeric"),
                          gaps_number = vector(mode = "numeric"), 
                          gaps_mean_size = vector(mode = "numeric"), 
                          relative_position_along_chromosome = vector(mode = "numeric"),
                          relative_size_against_chromosome_size = vector(mode = "numeric"),
                          arrayID = vector(mode = "character"))


repeats$is_on_chromosome <- FALSE
repeats$is_on_chromosome[repeats$seqID %in% chromosomes] <- TRUE


{
  genome_metadata$analysed_centromeric_Satellite_name <- ""
  genome_metadata$analysed_centromeric_Trash_name <- ""
  genome_metadata$Architecture <- ""
  genome_metadata$centromere_start <- 0
  genome_metadata$centromere_end <- 0
  genome_metadata$pericentromere_start <- 0
  genome_metadata$pericentromere_end <- 0
  genome_metadata$centromere_arrays_starts <- ""
  genome_metadata$centromere_arrays_ends <- ""
}

cat(fasta_name, "prepped for analysis \n")
which_sats <- grep(fasta_name_short, satellite_metadata$fasta)

if(length(which_sats) != 0) {
  
  cat("satellites found \n")
  repeats$is_centromeric <- FALSE
  
  genome_satellite_metadata <- satellite_metadata[which_sats,]
  
  centromeric_classes <- NULL
  
  for(k in seq_len(nrow(genome_satellite_metadata))) {
    centromeric_classes <- c(centromeric_classes, genome_satellite_metadata$TRASH_new_class[k])
  }
  
  
  genome_metadata$analysed_centromeric_Satellite_name <- paste(genome_satellite_metadata$TRASH_new_class, collapse = ";")
  genome_metadata$analysed_centromeric_Trash_name <- paste(centromeric_classes, collapse = ";")
  
  repeats$is_centromeric[repeats$new_class %in% centromeric_classes] <- TRUE
  
  cat(sum(repeats$is_centromeric), " centromeric repeats \n")
  
  arrays_genome <- NULL
  
  repeats$arrayID = 0
  
  if(!genome_satellite_metadata$is_holocentric[1]) {
    ########## Treat as satellite monocentric
    
    ### ARRAYS mono done
    for(j in seq_len(nrow(genome_metadata))) {
      # for(j in 1) {
      cat("MONOCENTRIC identifying arrays on chromosome", genome_metadata$chromosome.name[j],
          genome_metadata$assembly.name[j], "\n")
      repeats_chromosome <- repeats[repeats$seqID == genome_metadata$chromosome.name[j], ]
      
      repeats_chromosome <- repeats_chromosome[repeats_chromosome$is_centromeric,]
      
      if(nrow(repeats_chromosome) < 2) next
      
      
      repeats_chromosome <- repeats_chromosome[order(repeats_chromosome$start, decreasing = FALSE), ]
      
      chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
      
      repeats_chromosome$dist_to_next = 0
      repeats_chromosome$dist_to_next[1:(nrow(repeats_chromosome) - 1)] = 
        repeats_chromosome$start[2:nrow(repeats_chromosome)] - 
        repeats_chromosome$end[1:(nrow(repeats_chromosome) - 1)] - 1
      repeats_chromosome$dist_to_next[nrow(repeats_chromosome)] = 999999999
      
      big_gaps <- repeats_chromosome[repeats_chromosome$dist_to_next > max_gap_distance,]
      
      arrays <-  empty_array
      
      
      arrays[nrow(arrays) + 1,] = list(genome_metadata$chromosome.name[j],
                                       repeats_chromosome$start[1], 
                                       big_gaps$end[1], 
                                       "",0,0,"","",0,0,0,0,0,0,0,0,0,0)
      
      if(nrow(big_gaps) > 1) {
        for(k in 1 : (nrow(big_gaps) - 1)) {
          arrays[nrow(arrays) + 1,] = list(genome_metadata$chromosome.name[j],
                                           big_gaps$end[k] + big_gaps$dist_to_next[k] + 1, 
                                           big_gaps$end[k + 1], 
                                           "",0,0,"","",0,0,0,0,0,0,0,0,0,0)
        }
      } 
      
      # remove arrays consisting of one repeat (small bug in this version of TRASH which made some repeats extend over the array edge, which meant when arrays were removed, 
      # in the 1_filter_TRASH.R, those repeats stayed as they were technically outside of arrays)
      arrays_to_remove <- NULL
      if(nrow(arrays) > 1) {
        for(k in 1 : nrow(arrays)) {
          repeats_in_array <- length(which(repeats_chromosome$start > arrays$start[k] & repeats_chromosome$end < arrays$end[k]))
          if(repeats_in_array < 2) {
            arrays_to_remove <- c(arrays_to_remove, k)
          }
        }
      }
      if(length(arrays_to_remove) != 0) {
        arrays <- arrays[-arrays_to_remove, ]
      }
      
      # if more than 1 array, check if gaps contain TEs and merge if needed
      if(nrow(arrays) > 1) {
        arrays$to_merge_with_next <- FALSE
        for(k in 1 : (nrow(arrays) - 1)) {
          distance <- arrays$start[k+1] - arrays$end[k] - 1
          
          if(distance > (max_gap_distance / gap_TE_coverage_count)) next # too much, no matter what
          
          TEs_in_proximity <- chromosome_edta[chromosome_edta$V5 > (arrays$end[k] - 50000) & 
                                                chromosome_edta$V5 < (arrays$start[k+1] + 50000), ]
          if(nrow(TEs_in_proximity) == 0) next
          
          TE_coordinates_in_proximity <- unique(unlist(lapply(X = seq_len(nrow(TEs_in_proximity)), 
                                                              function(X) TEs_in_proximity$V5[X] : TEs_in_proximity$V6[X] )))
          TE_coordinates_in_gap <- sum(TE_coordinates_in_proximity >= arrays$end[k] & 
                                         TE_coordinates_in_proximity < arrays$start[k+1]) 
          
          distance <- distance - TE_coordinates_in_gap + (TE_coordinates_in_gap * gap_TE_coverage_count )
          
          if(distance <= max_gap_distance) {
            # merging!
            arrays$to_merge_with_next[k] = TRUE
            cat("found one to merge!\n")
          }
        }
        if(sum(arrays$to_merge_with_next) != 0) {
          for(k in nrow(arrays) : 2) {
            if(arrays$to_merge_with_next[k - 1]) {
              arrays$end[k-1] = arrays$end[k]
              arrays <- arrays[-k,]
            }
          }
        }
      }
      arrays$to_merge_with_next <- NULL
      
      if(nrow(arrays) == 0) next
      for(k in seq_len(nrow(arrays))) {
        cat("calculating characteristics of array", k, "/", nrow(arrays), 
            genome_metadata$assembly.name[j], "\n")
        
        array_repeats <- repeats_chromosome[repeats_chromosome$start >= arrays$start[k] & 
                                              repeats_chromosome$start < arrays$end[k], ]
        
        
        
        arrays$strand_percentage[k] <- 100 * sum(array_repeats$strand == "+") / nrow(array_repeats)
        if(arrays$strand_percentage[k] >= 50) {
          arrays$strand[k] = "+"
        } else {
          arrays$strand[k] = "-"
          arrays$strand_percentage[k] <- 100 - arrays$strand_percentage[k]
        }
        
        arrays$width[k] <- arrays$end[k] - arrays$start[k]
        arrays$satellite_name[k] <- genome_metadata$analysed_centromeric_Satellite_name[1]
        arrays$TRASH_class[k] <- genome_metadata$analysed_centromeric_Trash_name[1]
        arrays$repeats_number[k] <- nrow(array_repeats)
        
        arrays$mean_repeat_width[k] <- mean(array_repeats$width)
        arrays$mean_repeat_width_SD[k] <- sd(array_repeats$width)
        arrays$relative_position_along_chromosome[k] <- 100 * mean(c(arrays$start[k], arrays$end[k])) / genome_metadata$size[genome_metadata$chromosome.name == arrays$seqID[k]]
        arrays$relative_size_against_chromosome_size[k] <- 100 * arrays$width[k] / genome_metadata$size[genome_metadata$chromosome.name == arrays$seqID[k]]
        arrays$arrayID[k] <- paste0("array_", arrays$seqID[k], "_", arrays$satellite_name[k], "_", k)
        
        repeats$arrayID[repeats$seqID == genome_metadata$chromosome.name[j] & 
                          repeats$start >= arrays$start[k] & 
                          repeats$start < arrays$end[k]] <- arrays$arrayID[k]
        
        
        array_repeats$dist_to_next[nrow(array_repeats)] = 0
        arrays$gaps_number[k] = sum(array_repeats$dist_to_next > min_internal_gap_distance)
        arrays$gaps_mean_size[k] = mean(array_repeats$dist_to_next[array_repeats$dist_to_next > min_internal_gap_distance])
        
      }
      
      arrays_genome <- rbind(arrays_genome, arrays)
    } ### ARRAYS mono
    
    
    ### PERICENTROMERES mono  done
    repeats_expanded <- NULL
    arrays_expanded <- NULL
    for(j in seq_len(nrow(genome_metadata))) {
      # for(j in 1) {
      
      if(genome_metadata$is.chr[j] != 1) {
        
        cat("not a chromosome, no centromere will be identified\n")
        next
      }
      arrays <- arrays_genome[arrays_genome$seqID == genome_metadata$chromosome.name[j],]
      if(nrow(arrays) == 0) next
      
      
      repeats_chromosome <- repeats[repeats$seqID == genome_metadata$chromosome.name[j], ]
      repeats_chromosome <- repeats_chromosome[repeats_chromosome$is_centromeric,]
      
      other_repeats <- repeats[repeats$seqID == genome_metadata$chromosome.name[j], ]
      other_repeats <- other_repeats[!other_repeats$is_centromeric,]
      
      if(nrow(repeats_chromosome) > 0) {
        repeats_chromosome <- repeats_chromosome[order(repeats_chromosome$start, decreasing = FALSE), ]
        repeats_chromosome$in_centromere <- FALSE
      }
      if(nrow(repeats_chromosome) < 2) next
      
      
      
      cat("MONOCENTRIC identifying pericentromeres on chromosome", genome_metadata$chromosome.name[j],
          genome_metadata$assembly.name[j], "\n")
      
      chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
      
      # organisation = satellite_organisation_metadata$chromosome.acro...meta[satellite_organisation_metadata$assembly.name == fasta_name & 
      #                                                                         satellite_organisation_metadata$chromosome.name ==  genome_metadata$chromosome.name[j]]
      # 
      # 
      # if(length(organisation) == 1) genome_metadata$Architecture[j] <- organisation
      
      
      window_starts <- round(seq(1, genome_metadata$size[j], length.out = round(genome_metadata$size[j] / te_score_calculation_window)))
      window_ends <- window_starts - 1
      window_starts = window_starts[-length(window_starts)]
      window_ends = window_ends[-1]
      window_widths <- window_ends - window_starts
      
      
      if(length(window_starts) == 0) {
        window_starts = 1
        window_ends = genome_metadata$size[j]
        window_widths = genome_metadata$size[j] - 1
      }
      chromosome_edta$V5 <- as.numeric(chromosome_edta$V5)
      chromosome_edta$V6 <- as.numeric(chromosome_edta$V6)
      
      
      chromosome_centromeric_edta <- chromosome_edta[grep(pattern = cen_TE_grep_pattern, x = chromosome_edta$V4, perl = T),]
      
      
      cat("edta all", nrow(chromosome_edta), "and centromeric TEs with a regex:  ", cen_TE_grep_pattern, "  at count", nrow(chromosome_centromeric_edta), "\n")
      
      
      gr1 <- IRanges(window_starts, window_ends)
      gr2 <- IRanges(chromosome_centromeric_edta$V5, chromosome_centromeric_edta$V6)
      gr3 <- IRanges(repeats_chromosome$start[repeats_chromosome$is_centromeric], repeats_chromosome$end[repeats_chromosome$is_centromeric])
      gr4 <- IRanges(other_repeats$start, other_repeats$end)
      
      genes_chromosome <- genes[genes$V1 == genome_metadata$chromosome.name[j], ]
      gr5 <- IRanges(genes_chromosome$V4, genes_chromosome$V5)
      
      window_overlapping_bp <- NULL
      for(k in seq_along(gr1)) {
        window_overlapping_bp <- c(window_overlapping_bp, 
                                   (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr2))$width) + 0) * 2.5 # EDTA lineage multiplier
                                   + (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width) + 0) * 2   # CEN repeats multiplier
                                   + (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr4))$width) + 0) * 1  # All other repeats multiplier
                                   - (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr5))$width) + 0) * 0.5  # GENES deduction with multiplier
        )
      }
      window_overlapping_bp[window_overlapping_bp < 0] = 0
      window_overlapping_perc <- window_overlapping_bp / window_widths
      window_overlapping_perc[window_overlapping_perc > 1] = 1
      
      if(length(window_overlapping_perc) > 1) {
        window_overlapping_perc <-  ma(c(window_overlapping_perc[1],
                                         window_overlapping_perc[1],
                                         window_overlapping_perc,
                                         window_overlapping_perc[length(window_overlapping_perc)],
                                         window_overlapping_perc[length(window_overlapping_perc)]))[3:(length(window_overlapping_perc)+2)]
      }
      if(length(window_overlapping_perc) > 1) {
        window_overlapping_perc <-  ma(c(window_overlapping_perc[1],
                                         window_overlapping_perc[1],
                                         window_overlapping_perc,
                                         window_overlapping_perc[length(window_overlapping_perc)],
                                         window_overlapping_perc[length(window_overlapping_perc)]))[3:(length(window_overlapping_perc)+2)]
      }
      if(length(window_overlapping_perc) > 1) {
        window_overlapping_perc <-  ma(c(window_overlapping_perc[1],
                                         window_overlapping_perc[1],
                                         window_overlapping_perc,
                                         window_overlapping_perc[length(window_overlapping_perc)],
                                         window_overlapping_perc[length(window_overlapping_perc)]))[3:(length(window_overlapping_perc)+2)]
      }
      
      
      # otsu threshold implementation from EBImage::otsu
      range = c(0, 1)
      levels = 256
      breaks = seq(range[1], range[2], length.out = levels + 1)
      h = hist.default(window_overlapping_perc, breaks = breaks, plot = FALSE)
      counts = as.double(h$counts)
      mids = as.double(h$mids)
      len = length(counts)
      w1 = cumsum(counts)
      w2 = w1[len] + counts - w1
      cm = counts * mids
      m1 = cumsum(cm)
      m2 = m1[len] + cm - m1
      var = w1 * w2 * (m2/w2 - m1/w1)^2
      maxi = which(var == max(var, na.rm = TRUE))
      otsu_mid = (mids[maxi[1]] + mids[maxi[length(maxi)]])/2
      
      if(is.na(otsu_mid)) {
        otsu_mid <- 0.8
      } else if(!is.numeric(otsu_mid)) {
        otsu_mid <- 0.8
      } else if(otsu_mid <= 0) {
        otsu_mid <- 0.8
      }
      
      
      centromeric_repeat_weight <- round(median(c(chromosome_centromeric_edta$V5, chromosome_centromeric_edta$V6)))
      if(nrow(repeats_chromosome) > 0) {
        centromeric_repeat_weight <- round(median(repeats_chromosome$start))
      } 
      
      
      windows_df <- data.frame(start = window_starts, end = window_ends,
                               overlapping_bp = window_overlapping_bp,
                               overlapping_perc = window_overlapping_perc)
      windows_df$above_otsu <- FALSE
      windows_df$above_otsu[windows_df$overlapping_perc >= otsu_mid] = TRUE
      windows_df$centromeric_repeat_weight <- "no"
      windows_df$centromeric_repeat_weight[which(windows_df$start <= centromeric_repeat_weight & 
                                                   windows_df$end >= centromeric_repeat_weight)] = "yes"
      
      windows_df$above_otsu[windows_df$centromeric_repeat_weight == "yes"] = TRUE
      
      if(nrow(windows_df) > 2) {
        for(z in 1 : (nrow(windows_df) - 2)){
          if(windows_df$above_otsu[z] & !windows_df$above_otsu[z+1] & windows_df$above_otsu[z+2]) {
            windows_df$above_otsu[z+1] <- TRUE
          }
        }
      }
      
      
      
      windows_df$is_pericentromere <- FALSE
      windows_df$is_pericentromere[windows_df$centromeric_repeat_weight == "yes"] = TRUE
      k0 = which(windows_df$is_pericentromere)
      k = k0
      while(k <= nrow(windows_df) & windows_df$above_otsu[k]) {
        windows_df$is_pericentromere[k] <- TRUE
        k = k + 1
      }
      k = k0
      while(k > 1 & windows_df$above_otsu[k] ) {
        windows_df$is_pericentromere[k] <- TRUE
        k = k - 1
      }
      
      genome_metadata$pericentromere_start[j] <- min(windows_df$start[windows_df$is_pericentromere])
      genome_metadata$pericentromere_end[j] <- max(windows_df$end[windows_df$is_pericentromere])
      
      
      pericentromere_gr <- IRanges(genome_metadata$pericentromere_start[j], genome_metadata$pericentromere_end[j])
      centromeric_arrays <- NULL
      noncentromeric_arrays <- NULL
      
      arrays$in_centromere <- FALSE
      arrays$edge_filter <- FALSE
      for(l in seq_len(nrow(arrays))) {
        if((length(intersect(pericentromere_gr, IRanges(arrays$start[l], arrays$end[l]))) > 0) |
           centromeric_repeat_weight %in% (arrays$start[l] : arrays$end[l])) {
          arrays$in_centromere[l] <- TRUE
          centromeric_arrays <- rbind(centromeric_arrays, arrays[l,])
        } else {
          noncentromeric_arrays <- rbind(noncentromeric_arrays, arrays[l,])
        }
      }
      
      if(length(centromeric_arrays) == 0) {
        centromeric_arrays <- arrays
        noncentromeric_arrays <- NULL
      }
      
      
      if(nrow(repeats_chromosome) > 0) {
        repeats_chromosome$in_centromere[repeats_chromosome$arrayID %in% centromeric_arrays$arrayID] = TRUE
      }
      
      ### filter out centromeric arrays that are at the same time:
      # 1. smaller than 10% of all centromeric repeat length on this chromosome
      # 2. found at the edge of the repeats IDs (25% on each side of the all repeats count)
      
      if(nrow(arrays) != 0) {
        total_cen_bp <- sum(repeats_chromosome$width)
        which_smaller_than_10 <- which((arrays$repeats_number * arrays$mean_repeat_width) < (total_cen_bp/10))
        which_in_25_edges <- NULL
        if(length(which_smaller_than_10) != 0) {
          which_25th_repeat <-  round(nrow(repeats_chromosome) / 4)
          which_75th_repeat <-  round(nrow(repeats_chromosome) / 4) * 3
          
          left_edge_coordinates <- (repeats_chromosome$start[1]:repeats_chromosome$end[which_25th_repeat])
          right_edge_coordinates <- (repeats_chromosome$start[which_75th_repeat]:repeats_chromosome$end[nrow(repeats_chromosome)])
          
          which_in_25_edges <- which( ((arrays$start %in% left_edge_coordinates) & (arrays$end %in% left_edge_coordinates)) | 
                                        ((arrays$start %in% right_edge_coordinates) & (arrays$end %in% right_edge_coordinates)) )
        }
        if(length(intersect(which_in_25_edges, which_smaller_than_10)) != 0) {
          arrays$edge_filter[intersect(which_in_25_edges, which_smaller_than_10)] <- TRUE
        }
      }
      
      
      
      genome_metadata$centromere_start[j] <- min(centromeric_arrays$start)
      genome_metadata$centromere_end[j] <- max(centromeric_arrays$end)
      
      genome_metadata$centromere_arrays_starts[j] <- paste0(centromeric_arrays$start, collapse = ";")
      genome_metadata$centromere_arrays_ends[j] <- paste0(centromeric_arrays$end, collapse = ";")
      
      
      if(nrow(repeats_chromosome) > 0) {
        repeats_expanded <- rbind(repeats_expanded, repeats_chromosome)
      }
      
      arrays_expanded <- rbind(arrays_expanded, arrays)
      
    }  ### PERICENTROMERES mono
    
    
    ### PLOT done
    {
      pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/",
                        "cen_locations_plot_mono_", fasta_name, ".pdf"), 
          width = 10, height = 6, onefile = TRUE)
      
      which.chromosomes.in.table <- which(genome_metadata$is.chr == 1)
      for(j in which.chromosomes.in.table) {
        # for(j in 1) {
        cat("MONOCENTRIC plot", j, assembly_name, "\n")
        
        
        arrays <- arrays_expanded[arrays_expanded$seqID == genome_metadata$chromosome.name[j],]
        chromosome_repeats <- repeats[repeats$is_centromeric & repeats$seqID == genome_metadata$chromosome.name[j],]
        chromosome_repeats_full <- repeats[repeats$seqID == genome_metadata$chromosome.name[j],]
        chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
        chromosome_genes <- genes[genes$V1 == genome_metadata$chromosome.name[j],]
        
        cat("genes no ", nrow(genes), " chr genes ", nrow(chromosome_genes), "\n")
        
        pericentromere_gr <- IRanges(genome_metadata$pericentromere_start[j], genome_metadata$pericentromere_end[j])
        
        centromeric_arrays <- arrays[!arrays$edge_filter & arrays$in_centromere,]
        centromeric_arrays_on_edge <- arrays[arrays$edge_filter & arrays$in_centromere,]
        noncentromeric_arrays <- arrays[!arrays$in_centromere,]
        
        
        window_starts <- round(seq(1, genome_metadata$size[j], length.out = round(genome_metadata$size[j] / plotting_window)))
        window_ends <- window_starts - 1
        window_starts = window_starts[-length(window_starts)]
        window_ends = window_ends[-1]
        window_widths <- window_ends - window_starts
        
        
        if(length(window_starts) == 0) {
          window_starts = 1
          window_ends = genome_metadata$size[j]
          window_widths = genome_metadata$size[j] - 1
        }
        window_middles <- window_starts + window_widths/2
        
        
        gr1 <- IRanges(window_starts, window_ends)
        gr2 <- IRanges(chromosome_edta$V5, chromosome_edta$V6)
        gr3 <- IRanges(chromosome_repeats$start, chromosome_repeats$end)
        gr4 <- IRanges(chromosome_repeats_full$start, chromosome_repeats_full$end)
        gr5 <- IRanges(chromosome_genes$V4, chromosome_genes$V5)
        
        window_genes <- NULL
        window_tes <- NULL
        window_repeats <- NULL
        window_cen_repeats <- NULL
        for(k in seq_along(gr1)) {
          window_genes <- c(window_genes, 
                            (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr5))$width)))
          window_tes <- c(window_tes, 
                          (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr2))$width)))
          window_repeats <- c(window_repeats, 
                              (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr4))$width))
                              - (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width)))
          window_cen_repeats <- c(window_cen_repeats, 
                                  (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width)))
        }
        window_genes <- window_genes / window_widths
        window_genes[window_genes > 1] = 1
        window_genes <- window_genes * 100
        window_genes[window_genes < 0] = 0
        
        window_tes <- window_tes / window_widths
        window_tes[window_tes > 1] = 1
        window_tes <- window_tes * 100
        window_tes[window_tes < 0] = 0
        
        window_repeats <- window_repeats / window_widths
        window_repeats[window_repeats > 1] = 1
        window_repeats <- window_repeats * 100
        window_repeats[window_repeats < 0] = 0
        
        window_cen_repeats <- window_cen_repeats / window_widths
        window_cen_repeats[window_cen_repeats > 1] = 1
        window_cen_repeats <- window_cen_repeats * 100
        window_cen_repeats[window_cen_repeats < 0] = 0
        
        if(average_plots &  (length(window_genes) > 1)) {
          for(doit in 1 : 3) {
            window_genes <- ma(c(window_genes[1],
                                 window_genes[1],
                                 window_genes,
                                 window_genes[length(window_genes)],
                                 window_genes[length(window_genes)]))[3:(length(window_genes)+2)]
            
            window_tes <-  ma(c(window_tes[1],
                                window_tes[1],
                                window_tes,
                                window_tes[length(window_tes)],
                                window_tes[length(window_tes)]))[3:(length(window_tes)+2)]
            window_repeats <-  ma(c(window_repeats[1],
                                    window_repeats[1],
                                    window_repeats,
                                    window_repeats[length(window_repeats)],
                                    window_repeats[length(window_repeats)]))[3:(length(window_repeats)+2)]
            window_cen_repeats <-  ma(c(window_cen_repeats[1],
                                        window_cen_repeats[1],
                                        window_cen_repeats,
                                        window_cen_repeats[length(window_cen_repeats)],
                                        window_cen_repeats[length(window_cen_repeats)]))[3:(length(window_cen_repeats)+2)]
          }
          
        }
        window_genes[window_genes <= 0] = 0
        window_tes[window_tes <= 0] = 0
        window_repeats[window_repeats <= 0] = 0
        window_cen_repeats[window_cen_repeats <= 0] = 0
        
        plot(NA,NA,xlim = c(0,genome_metadata$size[j]), ylim = c(0,110), 
             xlab = "Coordinates, bp", ylab = "% window occupancy",
             main = paste("Monocentric", assembly_name, genome_metadata$chromosome.name[j], sep = " "))
        mtext(text = c("genes", "TEs",  paste0("centromeric ", assembly_name, " repeats"),  "other repeats"), 
              adj = c(.1,.33,.66,.9),
              col = c("#009900", "#000099", "#eeee00", "#990000"),
              side = 3, line = 0.2)
        # genes
        lines(y = window_genes, x = window_middles, col = "#009900", lwd = 3)
        polygon(y = c(0,window_genes, 0), border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#00990050")
        # TEs
        lines(y = window_tes, x = window_middles, col = "#000099", lwd = 3)
        polygon(y = c(0,window_tes, 0),  border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#00009950")
        # cen repeats
        lines(y = window_cen_repeats, x = window_middles, col = "#eeee00", lwd = 3)
        polygon(y = c(0,window_cen_repeats, 0),  border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#eeee0050")
        # other repeats
        lines(y = window_repeats, x = window_middles, col = "#990000", lwd = 3)
        polygon(y = c(0,window_repeats, 0),  border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#99000050")
        # pericentromere
        abline(h = c(105,108), lty = 5)
        polygon(x = c(genome_metadata$pericentromere_start[j], genome_metadata$pericentromere_start[j], genome_metadata$pericentromere_end[j], genome_metadata$pericentromere_end[j]),
                y = c(105,108,108,105), col = "#999999", border = F)
        axis(side = 2, at = 106.5, labels = "peric", las = 2, tck = 0)
        # centromere
        abline(h = c(110,113), lty = 5)
        if(length(noncentromeric_arrays) != 0) {
          for(k in seq_len(nrow(noncentromeric_arrays))) {
            polygon(x = c(noncentromeric_arrays$start[k], noncentromeric_arrays$start[k], noncentromeric_arrays$end[k], noncentromeric_arrays$end[k]),
                    y = c(110,113,113,110), col = "#000099", border = F)
          }
        }
        if(length(centromeric_arrays_on_edge) != 0) {
          for(k in seq_len(nrow(centromeric_arrays_on_edge))) {
            polygon(x = c(centromeric_arrays_on_edge$start[k], centromeric_arrays_on_edge$start[k], centromeric_arrays_on_edge$end[k], centromeric_arrays_on_edge$end[k]),
                    y = c(110,113,113,110), col = "#999999", border = F)
          }
        }
        if(length(centromeric_arrays) != 0) {
          for(k in seq_len(nrow(centromeric_arrays))) {
            polygon(x = c(centromeric_arrays$start[k], centromeric_arrays$start[k], centromeric_arrays$end[k], centromeric_arrays$end[k]),
                    y = c(110,113,113,110), col = "#000000", border = F)
          }
        }
        axis(side = 2, at = 111.5, labels = "cen", las = 2, tck = 0)
        
        
      }
      
      dev.off()
      
    } ### PLOT mono
    
    cat("plot done\n")
    
  } else {
    ########## Treat as satellite holocentric
    
    ### ARRAYS holo done
    for(j in seq_len(nrow(genome_metadata))) {
      # for(j in 1) {
      cat("HOLOCENTRIC identifying arrays onchromosome", genome_metadata$chromosome.name[j],
          genome_metadata$assembly.name[j], "\n")
      repeats_chromosome <- repeats[repeats$seqID == genome_metadata$chromosome.name[j], ]
      
      repeats_chromosome <- repeats_chromosome[repeats_chromosome$is_centromeric,]
      
      if(nrow(repeats_chromosome) < 2) next
      
      if(nrow(repeats_chromosome) > 0) {
        repeats_chromosome <- repeats_chromosome[order(repeats_chromosome$start, decreasing = FALSE), ]
        
        repeats_chromosome$in_centromere <- FALSE
      }
      
      
      chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
      
      
      repeats_chromosome$dist_to_next = 0
      repeats_chromosome$dist_to_next[1:(nrow(repeats_chromosome) - 1)] = 
        repeats_chromosome$start[2:nrow(repeats_chromosome)] - 
        repeats_chromosome$end[1:(nrow(repeats_chromosome) - 1)] - 1
      repeats_chromosome$dist_to_next[nrow(repeats_chromosome)] = 999999999
      
      big_gaps <- repeats_chromosome[repeats_chromosome$dist_to_next > max_gap_dist_holo,]
      
      arrays <-  empty_array
      
      
      arrays[nrow(arrays) + 1,] = list(genome_metadata$chromosome.name[j],
                                       repeats_chromosome$start[1], 
                                       big_gaps$end[1], 
                                       "",0,0,"","",0,0,0,0,0,0,0,0,0,0)
      
      if(nrow(big_gaps) > 1) {
        for(k in 1 : (nrow(big_gaps) - 1)) {
          arrays[nrow(arrays) + 1,] = list(genome_metadata$chromosome.name[j],
                                           big_gaps$end[k] + big_gaps$dist_to_next[k] + 1, 
                                           big_gaps$end[k + 1], 
                                           "",0,0,"","",0,0,0,0,0,0,0,0,0,0)
        }
      }
      
      # remove arrays consisting of one repeat (small bug in this version of TRASH which made some repeats extend over the array edge, which meant when arrays were removed, 
      # in the 1_filter_TRASH.R, those repeats stayed as they were technically outside of arrays)
      arrays_to_remove <- NULL
      if(nrow(arrays) > 1) {
        for(k in 1 : nrow(arrays)) {
          repeats_in_array <- length(which(repeats_chromosome$start > arrays$start[k] & repeats_chromosome$end < arrays$end[k]))
          if(repeats_in_array < 2) {
            arrays_to_remove <- c(arrays_to_remove, k)
          }
        }
      }
      if(length(arrays_to_remove) != 0) {
        arrays <- arrays[-arrays_to_remove, ]
      }
      
      
      
      # if more than 1 array, check if gaps contain TEs and merge if needed
      arrays$to_merge_with_next <- FALSE
      if(nrow(arrays) > 1) {
        for(k in 1 : (nrow(arrays) - 1)) {
          distance <- arrays$start[k+1] - arrays$end[k] - 1
          
          if(distance > (max_gap_dist_holo / gap_TE_coverage_count)) next # too much, no matter what
          
          TEs_in_proximity <- chromosome_edta[chromosome_edta$V5 > (arrays$end[k] - 50000) & 
                                                chromosome_edta$V5 < (arrays$start[k+1] + 50000), ]
          if(nrow(TEs_in_proximity) == 0) next
          
          TE_coordinates_in_proximity <- unique(unlist(lapply(X = seq_len(nrow(TEs_in_proximity)), 
                                                              function(X) TEs_in_proximity$V5[X] : TEs_in_proximity$V6[X] )))
          TE_coordinates_in_gap <- sum(TE_coordinates_in_proximity >= arrays$end[k] & 
                                         TE_coordinates_in_proximity < arrays$start[k+1]) 
          
          distance <- distance - TE_coordinates_in_gap + (TE_coordinates_in_gap * gap_TE_coverage_count )
          
          if(distance <= max_gap_dist_holo) {
            # merging!
            arrays$to_merge_with_next[k] = TRUE
            cat("found one to merge!\n")
          }
        }
        if(sum(arrays$to_merge_with_next) != 0) {
          for(k in nrow(arrays) : 2) {
            if(arrays$to_merge_with_next[k - 1]) {
              arrays$end[k-1] = arrays$end[k]
              arrays <- arrays[-k,]
            }
          }
        }
      }
      arrays$to_merge_with_next <- NULL
      
      
      for(k in seq_len(nrow(arrays))) {
        cat("calculating characteristics of array", k, "/", nrow(arrays), 
            genome_metadata$assembly.name[j], "\n")
        
        array_repeats <- repeats_chromosome[repeats_chromosome$start >= arrays$start[k] & 
                                              repeats_chromosome$start < arrays$end[k], ]
        
        
        
        arrays$strand_percentage[k] <- 100 * sum(array_repeats$strand == "+") / nrow(array_repeats)
        if(arrays$strand_percentage[k] >= 50) {
          arrays$strand[k] = "+"
        } else {
          arrays$strand[k] = "-"
          arrays$strand_percentage[k] <- 100 - arrays$strand_percentage[k]
        }
        
        arrays$width[k] <- arrays$end[k] - arrays$start[k]
        arrays$satellite_name[k] <- genome_metadata$analysed_centromeric_Satellite_name[1]
        arrays$TRASH_class[k] <- genome_metadata$analysed_centromeric_Trash_name[1]
        arrays$repeats_number[k] <- nrow(array_repeats)
        
        arrays$mean_repeat_width[k] <- mean(array_repeats$width)
        arrays$mean_repeat_width_SD[k] <- sd(array_repeats$width)
        arrays$relative_position_along_chromosome[k] <- 100 * mean(c(arrays$start[k], arrays$end[k])) / genome_metadata$size[genome_metadata$chromosome.name == arrays$seqID[k]]
        arrays$relative_size_against_chromosome_size[k] <- 100 * arrays$width[k] / genome_metadata$size[genome_metadata$chromosome.name == arrays$seqID[k]]
        arrays$arrayID[k] <- paste0("array_", arrays$seqID[k], "_", arrays$satellite_name[k], "_", k)
        
        repeats$arrayID[repeats$seqID == genome_metadata$chromosome.name[j] & 
                          repeats$start >= arrays$start[k] & 
                          repeats$start < arrays$end[k]] <- arrays$arrayID[k]
        
       
        
        
        array_repeats$dist_to_next[nrow(array_repeats)] = 0
        arrays$gaps_number[k] = sum(array_repeats$dist_to_next > min_internal_gap_distance)
        arrays$gaps_mean_size[k] = mean(array_repeats$dist_to_next[array_repeats$dist_to_next > min_internal_gap_distance])
        
      }
      
      arrays_genome <- rbind(arrays_genome, arrays)
    } ### ARRAYS holo
    
    
    
    ### PERICENTROMERES holo done
    repeats_expanded <- NULL
    arrays_expanded <- NULL
    for(j in seq_len(nrow(genome_metadata))) {
      # for(j in 1) {
      if(genome_metadata$is.chr[j] != 1) next
      cat("HOLOCENTRIC identifying pericentromeres on chromosome", genome_metadata$chromosome.name[j],
          genome_metadata$assembly.name[j], "\n")
      repeats_chromosome <- repeats[repeats$seqID == genome_metadata$chromosome.name[j], ]
      repeats_chromosome <- repeats_chromosome[repeats_chromosome$is_centromeric,]
      if(nrow(repeats_chromosome) < 2) next
      
      repeats_chromosome <- repeats_chromosome[order(repeats_chromosome$start, decreasing = FALSE), ]
      
      repeats_chromosome$in_centromere <- FALSE
      
      
      chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
      
      arrays <- arrays_genome[arrays_genome$seqID == genome_metadata$chromosome.name[j], ]
      
      if(nrow(arrays) == 0) next
      
      
      
      # holocentric centromere starts and ends
      genome_metadata$centromere_arrays_starts[j] = ""
      genome_metadata$centromere_arrays_ends[j] = ""
      
      cen_starts <- NULL
      cen_ends <- NULL
      for(k in seq_len(nrow(arrays))) {
        cen_starts <- c(cen_starts, (arrays$start[k]))
        cen_ends <- c(cen_ends, (arrays$end[k]))
      }
      
      genome_metadata$centromere_arrays_starts[j] <- paste(cen_starts, collapse = ";")
      genome_metadata$centromere_arrays_ends[j] <- paste(cen_ends, collapse = ";")
      
      # holocentric pericentromere starts and ends
      genome_metadata$pericentromere_start[j] = ""
      genome_metadata$pericentromere_end[j] = ""
      
      peri_starts <- NULL
      peri_ends <- NULL
      for(k in seq_len(nrow(arrays))) {
        peri_starts <- c(peri_starts, (arrays$start[k] - round(0.5*arrays$width[k])))
        peri_ends <- c(peri_ends, (arrays$end[k] + round(0.5*arrays$width[k])))
      }
      for(k in seq_along(peri_starts)) {
        if((k + 1) <= length(peri_starts)) {
          if(peri_ends[k] > peri_starts[k+1]) {
            middle <- mean(c(arrays$end[k], arrays$start[k+1]))
            peri_ends[k] <- middle
            peri_starts[k+1] <- middle+1
          }
        }
      }
      if(peri_starts[1] < 1) peri_starts[1] = 1
      if(peri_ends[length(peri_ends)] > genome_metadata$size[j]) peri_ends[length(peri_ends)] = genome_metadata$size[j]
      
      
      genome_metadata$pericentromere_start[j] <- paste(peri_starts, collapse = ";")
      genome_metadata$pericentromere_end[j] <- paste(peri_ends, collapse = ";")
      
      ## ##
      
      pericentromere_gr <- IRanges(peri_starts, peri_ends)
      centromeric_arrays <- arrays
      arrays$in_centromere <- TRUE
      noncentromeric_arrays <- NULL
      
      
      
      genome_metadata$centromere_start[j] <- min(centromeric_arrays$start)
      genome_metadata$centromere_end[j] <- max(centromeric_arrays$end)
      
      
      if(nrow(repeats_chromosome) > 0) {
        repeats_chromosome$in_centromere[repeats_chromosome$is_centromeric] = TRUE
        repeats_expanded <- rbind(repeats_expanded, repeats_chromosome)
      }
      
      arrays_expanded <- rbind(arrays_expanded, arrays)
      
      
      
    } ### PERICENTROMERES holo 
    
    
    ### plot holo
    
    {
      pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/",
                        "cen_locations_plot_holo_", fasta_name, ".pdf"), 
          width = 10, height = 6, onefile = TRUE)
      
      which.chromosomes.in.table <- which(genome_metadata$is.chr == 1)
      for(j in which.chromosomes.in.table) {
        # for(j in 1) {
        cat("HOLOCENTRIC plot", j, assembly_name, "\n")
        
        arrays <- arrays_genome[arrays_genome$seqID == genome_metadata$chromosome.name[j],]
        chromosome_repeats <- repeats[repeats$is_centromeric & repeats$seqID == genome_metadata$chromosome.name[j],]
        chromosome_repeats_full <- repeats[repeats$seqID == genome_metadata$chromosome.name[j],]
        chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
        chromosome_genes <- genes[genes$V1 == genome_metadata$chromosome.name[j],]
        
        cat("genes no ", nrow(genes), " chr genes ", nrow(chromosome_genes), "\n")
        
        pericentromere_gr <- IRanges(as.numeric(strsplit(genome_metadata$pericentromere_start[j], split = ";")[[1]]), as.numeric(strsplit(genome_metadata$pericentromere_end[j], split = ";")[[1]]))
        centromeric_arrays <- arrays
        
        
        window_starts <- round(seq(1, genome_metadata$size[j], length.out = round(genome_metadata$size[j] / plotting_window)))
        window_ends <- window_starts - 1
        window_starts = window_starts[-length(window_starts)]
        window_ends = window_ends[-1]
        window_widths <- window_ends - window_starts
        
        
        if(length(window_starts) == 0) {
          window_starts = 1
          window_ends = genome_metadata$size[j]
          window_widths = genome_metadata$size[j] - 1
        }
        window_middles <- window_starts + window_widths/2
        
        
        gr1 <- IRanges(window_starts, window_ends)
        gr2 <- IRanges(chromosome_edta$V5, chromosome_edta$V6)
        gr3 <- IRanges(chromosome_repeats$start, chromosome_repeats$end)
        gr4 <- IRanges(chromosome_repeats_full$start, chromosome_repeats_full$end)
        gr5 <- IRanges(chromosome_genes$V4, chromosome_genes$V5)
        
        window_genes <- NULL
        window_tes <- NULL
        window_repeats <- NULL
        window_cen_repeats <- NULL
        for(k in seq_along(gr1)) {
          window_genes <- c(window_genes, 
                            (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr5))$width)))
          window_tes <- c(window_tes, 
                          (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr2))$width)))
          window_repeats <- c(window_repeats, 
                              (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr4))$width))
                              - (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width)))
          window_cen_repeats <- c(window_cen_repeats, 
                                  (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width)))
        }
        window_genes <- window_genes / window_widths
        window_genes[window_genes > 1] = 1
        window_genes <- window_genes * 100
        window_genes[window_genes < 0] = 0
        
        window_tes <- window_tes / window_widths
        window_tes[window_tes > 1] = 1
        window_tes <- window_tes * 100
        window_tes[window_tes < 0] = 0
        
        window_repeats <- window_repeats / window_widths
        window_repeats[window_repeats > 1] = 1
        window_repeats <- window_repeats * 100
        window_repeats[window_repeats < 0] = 0
        
        window_cen_repeats <- window_cen_repeats / window_widths
        window_cen_repeats[window_cen_repeats > 1] = 1
        window_cen_repeats <- window_cen_repeats * 100
        window_cen_repeats[window_cen_repeats < 0] = 0
        
        
        if(average_plots &  (length(window_genes) > 1)) {
          for(doit in 1 : 3) {
            window_genes <- ma(c(window_genes[1],
                                 window_genes[1],
                                 window_genes,
                                 window_genes[length(window_genes)],
                                 window_genes[length(window_genes)]))[3:(length(window_genes)+2)]
            
            window_tes <-  ma(c(window_tes[1],
                                window_tes[1],
                                window_tes,
                                window_tes[length(window_tes)],
                                window_tes[length(window_tes)]))[3:(length(window_tes)+2)]
            window_repeats <-  ma(c(window_repeats[1],
                                    window_repeats[1],
                                    window_repeats,
                                    window_repeats[length(window_repeats)],
                                    window_repeats[length(window_repeats)]))[3:(length(window_repeats)+2)]
            window_cen_repeats <-  ma(c(window_cen_repeats[1],
                                        window_cen_repeats[1],
                                        window_cen_repeats,
                                        window_cen_repeats[length(window_cen_repeats)],
                                        window_cen_repeats[length(window_cen_repeats)]))[3:(length(window_cen_repeats)+2)]
          }
          
        }
        window_genes[window_genes <= 0] = 0
        window_tes[window_tes <= 0] = 0
        window_repeats[window_repeats <= 0] = 0
        window_cen_repeats[window_cen_repeats <= 0] = 0
        
        plot(NA,NA,xlim = c(0,genome_metadata$size[j]), ylim = c(0,110), 
             xlab = "Coordinates, bp", ylab = "% window occupancy",
             main = paste("Holocentric", assembly_name, genome_metadata$chromosome.name[j], sep = " "))
        mtext(text = c("genes", "TEs",  paste0("centromeric ", assembly_name, " repeats"),  "other repeats"), 
              adj = c(.1,.33,.66,.9),
              col = c("#009900", "#000099", "#eeee00", "#990000"),
              side = 3, line = 0.2)
        # genes
        lines(y = window_genes, x = window_middles, col = "#009900", lwd = 3)
        polygon(y = c(0,window_genes, 0), border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#00990050")
        # TEs
        lines(y = window_tes, x = window_middles, col = "#000099", lwd = 3)
        polygon(y = c(0,window_tes, 0),  border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#00009950")
        # cen repeats
        lines(y = window_cen_repeats, x = window_middles, col = "#eeee00", lwd = 3)
        polygon(y = c(0,window_cen_repeats, 0),  border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#eeee0050")
        # other repeats
        lines(y = window_repeats, x = window_middles, col = "#990000", lwd = 3)
        polygon(y = c(0,window_repeats, 0),  border = F,
                x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#99000050")
        # pericentromere
        abline(h = c(105,108), lty = 5)
        if(is.character(genome_metadata$pericentromere_start[j]) & genome_metadata$pericentromere_start[j] != "") {
          starts <- as.numeric(strsplit(genome_metadata$pericentromere_start[j], split = ";")[[1]])
          ends <- as.numeric(strsplit(genome_metadata$pericentromere_end[j], split = ";")[[1]])
          for(k in seq_along(starts)) {
            polygon(x = c(starts[k], starts[k], ends[k], ends[k]),
                    y = c(105,108,108,105), col = "#660044", border = F)
          }
        }
        axis(side = 2, at = 106.5, labels = "peric", las = 2, tck = 0)
        # centromere
        abline(h = c(110,113), lty = 5)
        if(is.character(genome_metadata$centromere_arrays_starts[j]) & genome_metadata$centromere_arrays_starts[j] != "") {
          starts <- as.numeric(strsplit(genome_metadata$centromere_arrays_starts[j], split = ";")[[1]])
          ends <- as.numeric(strsplit(genome_metadata$centromere_arrays_ends[j], split = ";")[[1]])
          for(k in seq_along(starts)) {
            polygon(x = c(starts[k], starts[k], ends[k], ends[k]),
                    y = c(110,113,113,110), col = "#bb0088", border = F)
          }
        }
        
        
        axis(side = 2, at = 111.5, labels = "cen", las = 2, tck = 0)
        
        
      }
      
      dev.off()
      
    } ### plot holo
    
    cat("plot done\n")
    
    
    
  }
  
  write.csv(x = arrays_expanded , file = paste0("./", fasta_name, "_centromeric_arrays.csv"), row.names = FALSE )
  
  write.csv(x = repeats_expanded , file = paste0("./", fasta_name, "_repeats_filtered_array_assigned.csv"), row.names = FALSE )
  
  write.csv(x = genome_metadata , file = paste0("./",  fasta_name, "_genome_metadata.csv"), row.names = FALSE )
  
  write.csv(x = arrays_expanded , file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/", fasta_name, "_centromeric_arrays.csv"), row.names = FALSE )
  
  write.csv(x = repeats_expanded , file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/", fasta_name, "_repeats_filtered_array_assigned.csv"), row.names = FALSE )
  
  write.csv(x = genome_metadata , file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/",  fasta_name, "_genome_metadata.csv"), row.names = FALSE )
  
  
  
  
  
} else {
  ########## Treat as TE based
  
  ### NO ARRAYS
  
  
  
  
  ### PERICENTROMERES TE done
  repeats_expanded <- NULL
  for(j in seq_len(nrow(genome_metadata))) {
    # for(j in 1) {
    if(genome_metadata$is.chr[j] != 1) next
    cat("TE ONLY identifying pericentromeres on chromosome", genome_metadata$chromosome.name[j],
        genome_metadata$assembly.name[j], "\n")
    repeats_chromosome <- repeats[repeats$seqID == genome_metadata$chromosome.name[j], ]
    if(nrow(repeats_chromosome) != 0) {
      repeats_chromosome <- repeats_chromosome[order(repeats_chromosome$start, decreasing = FALSE), ]
    }
    
    
    
    chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
    
    
    
    
    cat("identifying pericentromeres on chromosome", genome_metadata$chromosome.name[j],
        genome_metadata$assembly.name[j], "\n")
    
    # organisation = satellite_organisation_metadata$chromosome.acro...meta[satellite_organisation_metadata$assembly.name == fasta_name & 
    #                                                                         satellite_organisation_metadata$chromosome.name ==  genome_metadata$chromosome.name[j]]
    # 
    # if(length(organisation) == 1) genome_metadata$Architecture[j] <- organisation
    
    
    window_starts <- round(seq(1, genome_metadata$size[j], length.out = round(genome_metadata$size[j] / te_score_calculation_window)))
    window_ends <- window_starts - 1
    window_starts = window_starts[-length(window_starts)]
    window_ends = window_ends[-1]
    window_widths <- window_ends - window_starts
    
    if(length(window_starts) == 0) {
      window_starts = 1
      window_ends = genome_metadata$size[j]
      window_widths = genome_metadata$size[j] - 1
    }
    
    chromosome_centromeric_edta <- chromosome_edta[grep(pattern = cen_TE_grep_pattern, x = chromosome_edta$V4, perl = T),]
    
    cat("edta all", nrow(chromosome_edta), "and centromeric TEs with a regex:  ", cen_TE_grep_pattern, "  at count", nrow(chromosome_centromeric_edta), "\n")
    
    
    gr1 <- IRanges(window_starts, window_ends)
    gr2 <- IRanges(chromosome_centromeric_edta$V5, chromosome_centromeric_edta$V6)
    gr3 <- IRanges(repeats_chromosome$start, repeats_chromosome$end)
    
    genes_chromosome <- genes[genes$V1 == genome_metadata$chromosome.name[j], ]
    gr5 <- IRanges(genes_chromosome$V4, genes_chromosome$V5)
    
    window_overlapping_bp <- NULL
    for(k in seq_along(gr1)) {
      window_overlapping_bp <- c(window_overlapping_bp, 
                                 (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr2))$width) + 0) * 2.5 # EDTA lineage multiplier
                                 + (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width) + 0) * 1   # all repeats multiplier
                                 - (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr5))$width) + 0) * 0.5  # GENES deduction with multiplier
      )
    }
    window_overlapping_bp[window_overlapping_bp < 0] = 0
    window_overlapping_perc <- window_overlapping_bp / window_widths
    window_overlapping_perc[window_overlapping_perc > 1] = 1
    
    if(length(window_overlapping_perc) > 1) {
      window_overlapping_perc <-  ma(c(window_overlapping_perc[1],
                                       window_overlapping_perc[1],
                                       window_overlapping_perc,
                                       window_overlapping_perc[length(window_overlapping_perc)],
                                       window_overlapping_perc[length(window_overlapping_perc)]))[3:(length(window_overlapping_perc)+2)]
    }
    if(length(window_overlapping_perc) > 1) {
      window_overlapping_perc <-  ma(c(window_overlapping_perc[1],
                                       window_overlapping_perc[1],
                                       window_overlapping_perc,
                                       window_overlapping_perc[length(window_overlapping_perc)],
                                       window_overlapping_perc[length(window_overlapping_perc)]))[3:(length(window_overlapping_perc)+2)]
    }
    if(length(window_overlapping_perc) > 1) {
      window_overlapping_perc <-  ma(c(window_overlapping_perc[1],
                                       window_overlapping_perc[1],
                                       window_overlapping_perc,
                                       window_overlapping_perc[length(window_overlapping_perc)],
                                       window_overlapping_perc[length(window_overlapping_perc)]))[3:(length(window_overlapping_perc)+2)]
    }
    
    # otsu threshold implementation from EBImage::otsu
    range = c(0, 1)
    levels = 256
    breaks = seq(range[1], range[2], length.out = levels + 1)
    h = hist.default(window_overlapping_perc, breaks = breaks, plot = FALSE)
    counts = as.double(h$counts)
    mids = as.double(h$mids)
    len = length(counts)
    w1 = cumsum(counts)
    w2 = w1[len] + counts - w1
    cm = counts * mids
    m1 = cumsum(cm)
    m2 = m1[len] + cm - m1
    var = w1 * w2 * (m2/w2 - m1/w1)^2
    maxi = which(var == max(var, na.rm = TRUE))
    otsu_mid = (mids[maxi[1]] + mids[maxi[length(maxi)]])/2
    
    
    
    ### pericentromere middle point calculation
    
    window_starts2 <- round(seq(1, genome_metadata$size[j], length.out = 50))
    window_ends2 <- window_starts2 - 1
    window_starts2 = window_starts2[-length(window_starts2)]
    window_ends2 = window_ends2[-1]
    window_widths2 <- window_ends2 - window_starts2
    gr1 <- IRanges(window_starts2, window_ends2)
    
    window_overlapping_bp2 <- NULL
    for(k in seq_along(gr1)) {
      window_overlapping_bp2 <- c(window_overlapping_bp2, 
                                  (window_widths2[k] - sum(as.data.frame(setdiff(gr1[k], gr2))$width) + 0) * 2 # EDTA lineage multiplier
                                  + (window_widths2[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width) + 0) * 0.75   # all repeats multiplier
      )
    }
    window_overlapping_perc2 <- window_overlapping_bp2 / window_widths2
    window_overlapping_perc2[window_overlapping_perc2 > 1] = 1
    
    if(length(window_overlapping_perc2) > 1) {
      window_overlapping_perc2 <-  ma(c(window_overlapping_perc2[1],
                                        window_overlapping_perc2[1],
                                        window_overlapping_perc2,
                                        window_overlapping_perc2[length(window_overlapping_perc2)],
                                        window_overlapping_perc2[length(window_overlapping_perc2)]))[3:(length(window_overlapping_perc2)+2)]
    }
    if(length(window_overlapping_perc2) > 1) {
      window_overlapping_perc2 <-  ma(c(window_overlapping_perc2[1],
                                        window_overlapping_perc2[1],
                                        window_overlapping_perc2,
                                        window_overlapping_perc2[length(window_overlapping_perc2)],
                                        window_overlapping_perc2[length(window_overlapping_perc2)]))[3:(length(window_overlapping_perc2)+2)]
    }
    
    
    which_top_window <- which.max(window_overlapping_perc2)
    centromeric_repeat_weight <- round(window_starts2[which_top_window] + window_widths2[which_top_window]/2)
    
    
    
    windows_df <- data.frame(start = window_starts, end = window_ends,
                             overlapping_bp = window_overlapping_bp,
                             overlapping_perc = window_overlapping_perc)
    windows_df$above_otsu <- FALSE
    windows_df$above_otsu[windows_df$overlapping_perc >= otsu_mid] = TRUE
    windows_df$centromeric_repeat_weight <- "no"
    windows_df$centromeric_repeat_weight[which((centromeric_repeat_weight >= windows_df$start) & 
                                                 (centromeric_repeat_weight <= windows_df$end))] = "yes"
    
    windows_df$above_otsu[windows_df$centromeric_repeat_weight == "yes"] = TRUE
    
    if(nrow(windows_df) > 2) {
      for(z in 1 : (nrow(windows_df) - 2)){
        if(windows_df$above_otsu[z] & !windows_df$above_otsu[z+1] & windows_df$above_otsu[z+2]) {
          windows_df$above_otsu[z+1] <- TRUE
        }
      }
    }
    
    
    
    if(!windows_df$above_otsu[windows_df$centromeric_repeat_weight == "yes"]) {
      cat("pericentromeres couldn't be identified", genome_metadata$assembly.name[j], genome_metadata$chromosome.name[j], j, "/", nrow(genome_metadata), "\n")
    } else {
      windows_df$is_pericentromere <- FALSE
      windows_df$is_pericentromere[windows_df$centromeric_repeat_weight == "yes"] = TRUE
      k0 = which(windows_df$is_pericentromere)
      k = k0
      while(k <= nrow(windows_df) & windows_df$above_otsu[k]) {
        windows_df$is_pericentromere[k] <- TRUE
        k = k + 1
      }
      k = k0
      while(k > 1 & windows_df$above_otsu[k] ) {
        windows_df$is_pericentromere[k] <- TRUE
        k = k - 1
      }
      
      genome_metadata$pericentromere_start[j] <- min(windows_df$start[windows_df$is_pericentromere])
      genome_metadata$pericentromere_end[j] <- max(windows_df$end[windows_df$is_pericentromere])
    }
    
    
    if(nrow(repeats_chromosome) > 0) {
      repeats_expanded <- rbind(repeats_expanded, repeats_chromosome)
    }
  } ### PERICENTROMERES TE 
  
  
  
  ### PLOT TE done
  
  {
    pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/",
                      "cen_locations_plot_other_", fasta_name, ".pdf"), 
        width = 10, height = 6, onefile = TRUE)
    
    which.chromosomes.in.table <- which(genome_metadata$is.chr == 1)
    for(j in which.chromosomes.in.table) {
      # for(j in 1) {
      cat("TE ONLY plot", j, assembly_name, "\n")
      
      chromosome_repeats_full <- repeats[repeats$seqID == genome_metadata$chromosome.name[j],]
      chromosome_edta <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
      chromosome_genes <- genes[genes$V1 == genome_metadata$chromosome.name[j],]
      
      cat("genes no ", nrow(genes), " chr genes ", nrow(chromosome_genes), "\n")
      
      pericentromere_gr <- IRanges(as.numeric(genome_metadata$pericentromere_start[j]), as.numeric(genome_metadata$pericentromere_end[j]))
      
      window_starts <- round(seq(1, genome_metadata$size[j], length.out = round(genome_metadata$size[j] / plotting_window)))
      window_ends <- window_starts - 1
      window_starts = window_starts[-length(window_starts)]
      window_ends = window_ends[-1]
      window_widths <- window_ends - window_starts
      
      
      if(length(window_starts) == 0) {
        window_starts = 1
        window_ends = genome_metadata$size[j]
        window_widths = genome_metadata$size[j] - 1
      }
      window_middles <- window_starts + window_widths/2
      
      
      gr1 <- IRanges(window_starts, window_ends)
      gr2 <- IRanges(chromosome_edta$V5, chromosome_edta$V6)
      # gr3 <- IRanges(chromosome_repeats$start, chromosome_repeats$end)
      gr4 <- IRanges(chromosome_repeats_full$start, chromosome_repeats_full$end)
      gr5 <- IRanges(chromosome_genes$V4, chromosome_genes$V5)
      
      window_genes <- NULL
      window_tes <- NULL
      window_repeats <- NULL
      window_cen_repeats <- NULL
      for(k in seq_along(gr1)) {
        window_genes <- c(window_genes, 
                          (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr5))$width)))
        window_tes <- c(window_tes, 
                        (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr2))$width)))
        window_repeats <- c(window_repeats, 
                            (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr4))$width))
                            # - (window_widths[k] - sum(as.data.frame(setdiff(gr1[k], gr3))$width))
        )
      }
      window_genes <- window_genes / window_widths
      window_genes[window_genes > 1] = 1
      window_genes <- window_genes * 100
      window_genes[window_genes < 0] = 0
      
      window_tes <- window_tes / window_widths
      window_tes[window_tes > 1] = 1
      window_tes <- window_tes * 100
      window_tes[window_tes < 0] = 0
      
      window_repeats <- window_repeats / window_widths
      window_repeats[window_repeats > 1] = 1
      window_repeats <- window_repeats * 100
      window_repeats[window_repeats < 0] = 0

      
      if(average_plots &  (length(window_genes) > 1)) {
        for(doit in 1 : 3) {
          window_genes <- ma(c(window_genes[1],
                               window_genes[1],
                               window_genes,
                               window_genes[length(window_genes)],
                               window_genes[length(window_genes)]))[3:(length(window_genes)+2)]
          
          window_tes <-  ma(c(window_tes[1],
                              window_tes[1],
                              window_tes,
                              window_tes[length(window_tes)],
                              window_tes[length(window_tes)]))[3:(length(window_tes)+2)]
          window_repeats <-  ma(c(window_repeats[1],
                                  window_repeats[1],
                                  window_repeats,
                                  window_repeats[length(window_repeats)],
                                  window_repeats[length(window_repeats)]))[3:(length(window_repeats)+2)]
        }
        
      }
      window_genes[window_genes <= 0] = 0
      window_tes[window_tes <= 0] = 0
      window_repeats[window_repeats <= 0] = 0
      # window_cen_repeats[window_cen_repeats <= 0] = 0
      
      plot(NA,NA,xlim = c(0,genome_metadata$size[j]), ylim = c(0,110), 
           xlab = "Coordinates, bp", ylab = "% window occupancy",
           main = paste("TE based", assembly_name, genome_metadata$chromosome.name[j], sep = " "))
      mtext(text = c("genes", "TEs", "other repeats"), 
            adj = c(.1,.33,.66,.9),
            col = c("#009900", "#000099", "#990000"),
            side = 3, line = 0.2)
      # genes
      lines(y = window_genes, x = window_middles, col = "#009900", lwd = 3)
      polygon(y = c(0,window_genes, 0), border = F,
              x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#00990050")
      # TEs
      lines(y = window_tes, x = window_middles, col = "#000099", lwd = 3)
      polygon(y = c(0,window_tes, 0),  border = F,
              x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#00009950")
      
      # other repeats
      lines(y = window_repeats, x = window_middles, col = "#990000", lwd = 3)
      polygon(y = c(0,window_repeats, 0),  border = F,
              x = c(window_middles[1],window_middles, window_middles[length(window_middles)]), col="#99000050")
      # pericentromere
      abline(h = c(105,108), lty = 5)
      polygon(x = c(genome_metadata$pericentromere_start[j], genome_metadata$pericentromere_start[j], genome_metadata$pericentromere_end[j], genome_metadata$pericentromere_end[j]),
              y = c(105,108,108,105), col = "#228822", border = F)
      axis(side = 2, at = 106.5, labels = "peric", las = 2, tck = 0)
      # centromere
      abline(h = c(110,113), lty = 5)
      
      axis(side = 2, at = 111.5, labels = "cen", las = 2, tck = 0)
      
      
    }
    
    dev.off()
    
  } ### PLOT TE
  
  cat("plot done\n")
  
  
  write.csv(x = repeats_expanded , file = paste0("./", fasta_name, "_repeats_filtered_array_assigned.csv"), row.names = FALSE )
  
  write.csv(x = genome_metadata , file = paste0("./", fasta_name, "_genome_metadata.csv"), row.names = FALSE )
  
  write.csv(x = repeats_expanded , file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/", fasta_name, "_repeats_filtered_array_assigned.csv"), row.names = FALSE )
  
  write.csv(x = genome_metadata , file = paste0("/home/pwlodzimierz/ToL/upload_files/55_cen_coords/", fasta_name, "_genome_metadata.csv"), row.names = FALSE )
  
  
  
  
  
}











