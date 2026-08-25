

# the idea is to get the value of each TE family density within centromeric arrays, around them and outside of them,
# then the values within and around are normalised by the value outside, which is the "baseline" for the enrichment
# Both can be then represented on a plot where the middle is value 0, left is -10 and right is +10, as the values 
# will be on the log scale (so from -10 bilion to +10 bilion relative enrichment). 
# To get a single value, each region will be counted together from the whole genome (chromosomes only), alternative
# would be to get the results from each chromosome and average them out, but that might be problematic when
# there are chromosomes that are very short or don't have any cen repeats attached. The "around" and "outside"
# transition should also be well picked. Flat 50 kbp sounds good, should be enough to capture the proximity of 
# arrays without being too big bridging array interspace in the holocentric species. If the gap between 2 arrays
# is less than 50 kb, it's internal, if it's more, then it is divided between around and outside. 
# Last issue are the non-family repeats in these regions. They would affect the calculations a lot, and
# we know they exist especially around the centromeres, or in the politypic species, so the non-cen repeat arrays
# should be excluded from the counted space. 





taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i) # i = 435

max_gap_size <- 50000




edta_classes <- list(
  # class I (retrotransposons)
  ## LTR retrotransposons
  c("Gypsy_LTR_retrotransposon"),
  c("Copia_LTR_retrotransposon"),
  c("Bel_Pao_LTR_retrotransposon"),
  c("TRIM_LTR_retrotransposon"),
  c("Caulimoviridae"),
  c("Retrovirus", "LTR_retrotransposon", "long_terminal_repeat"),
  ## Non-LTR retrotransposons
  c("LINE_element"),
  c("SINE_element"),
  c("Penelope_retrotransposon"),
  c("DIRS_YR_retrotransposon"),
  c("non_LTR_retrotransposon"),
  # class II (DNA transposons)
  ## TIRs
  c("Kolobok_TIR_transposon" , "Ginger_TIR_transposon", "Academ_TIR_transposon", "Novosib_TIR_transposon", "Sola_TIR_transposon", "Merlin_TIR_transposon", "IS3EU_TIR_transposon", "PiggyBac_TIR_transposon", "hAT_TIR_transposon", "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon", "Dada_TIR_transposon", "CACTA_TIR_transposon", "Zisupton_TIR_transposon", "PIF_Harbinger_TIR_transposon"),
  ## other class II
  c("DNA_transposon"),
  c("helitron"),
  c("MITE"),
  c("Maverick_Polinton", "polinton"),
  # other, recombinase element based
  c("Tyrosine_Recombinase_Elements", "Crypton_Tyrosine_Recombinase"),
  # others
  c("TE", "TE_unclass"),
  # likely not TEs, remove for plotting?
  c("repeat_region", "SUPER", "Sequence_Ontology", "rRNA_gene", "target_site_duplication", "chr"))
edta_classes_colours <-  c(
  "#E31A1C",  # red
  "#D55E00",  # reddish-orange
  "#F5793A",  # bright orange
  "#FF7F00",  # orange
  "#FDBF6F",  # peach
  "#F0E442",  # yellow
  "#6A3D9A",  # dark purple
  "#A95AA1",  # purple
  "#CC79A7",  # pink
  "#DDA0DD",  # light pinkish purple (plum)
  "#CAB2D6",  # lavender
  "#1F78B4",  # blue
  "#B2DF8A",  # light green
  "#33A02C",  # green
  "#009E73",  # teal green
  "#56B4E9",  # light blue
  "#7F7F7F",  # grey
  "#000000",   # black
  "#000000"   # black
)


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))


suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(Biostrings))
suppressMessages(library(GenomicRanges))



setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")


data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2025", recursive = FALSE, full.names = TRUE))
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
setwd(data_directories[i])
print(getwd())

assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])

genome_metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full_2025.csv")
genome_metadata <- genome_metadata[grep(assembly_name, genome_metadata$assembly.name),]
if(nrow(genome_metadata) == 0) stop(paste0("No genome metadata in ", assembly_name))
genome_metadata <- genome_metadata[genome_metadata$is.chr == 1,]

satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/cen_satellite_families_march_2026.csv")
satellite_metadata <- satellite_metadata[grep(assembly_name, satellite_metadata$fasta),]
if(nrow(satellite_metadata) == 0) {
  (print(paste0("No satellites in ", assembly_name)))
  next
}



no_heli = FALSE
helixer_file = list.files(pattern = "helixer_filtered.csv", full.names = TRUE)
if(length(helixer_file) != 1) {
  print(paste0(i, "No genes!"))
  no_heli = TRUE
}

repeats_all_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeats_all_file) != 1) {print(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}

array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
if(length(array_file) != 1) {print(paste0(i, " no arrays!")); setwd(".."); quit(save = "no", status = 1)}

no_edta <- FALSE
if(!no_edta) {
  edta_file = list.files(pattern = paste0(assembly_name, "_edta_filtered.csv"), full.names = TRUE)
  if(length(edta_file) > 1) {
    edta_file <- edta_file[grep("reassigned", edta_file)]
  }
  if(length(edta_file) != 1) {print(paste0(i, " no edta!")); no_edta = TRUE}
  
}

repeats_all <- read.csv(repeats_all_file)
arrays <- read.csv(array_file)
helixer <- read.csv(helixer_file)
helixer <- helixer[helixer$V3 == "gene",]
edta <- read.csv(edta_file, header = F)

edta$V5 <- as.numeric(edta$V5)
edta$V6 <- as.numeric(edta$V6)
edta <- edta[!is.na(edta$V5),]
edta <- edta[!is.na(edta$V6),]

helixer$V4 <- as.numeric(helixer$V4)
helixer$V5 <- as.numeric(helixer$V5)
helixer <- helixer[!is.na(helixer$V4),]
helixer <- helixer[!is.na(helixer$V5),]

repeats <- repeats_all[repeats_all$new_class %in% satellite_metadata$TRASH_new_clas,]

for(repeat_ID in 1 : nrow(satellite_metadata)) {
  repeats_cen_sat <- repeats[repeats$new_class == satellite_metadata$TRASH_new_class[repeat_ID], ]
  og_cen_sat_names <- unique(repeats_cen_sat$class)
  arrays_cen_sat <- arrays[arrays$class %in% og_cen_sat_names, ]
  
  within_bp_total <- 0
  within_TE_bp_total <- 0
  within_TE_bp_C1_LTR <- 0
  within_TE_bp_C1_nonLTR <- 0
  within_TE_bp_C2_TIR <- 0
  within_TE_bp_C2_other <- 0
  within_TE_bp_other <- 0
  within_TE_bp_Gypsy <- 0
  within_TE_bp_Copia <- 0
  within_TE_bp_Line <- 0
  within_TE_bp_Sine <- 0
  within_gene_bp_total <- 0
  
  around_bp_total <- 0
  around_TE_bp_total <- 0
  around_TE_bp_C1_LTR <- 0
  around_TE_bp_C1_nonLTR <- 0
  around_TE_bp_C2_TIR <- 0
  around_TE_bp_C2_other <- 0
  around_TE_bp_other <- 0
  around_TE_bp_Gypsy <- 0
  around_TE_bp_Copia <- 0
  around_TE_bp_Line <- 0
  around_TE_bp_Sine <- 0
  around_gene_bp_total <- 0
  
  outside_bp_total <- 0
  outside_TE_bp_total <- 0
  outside_TE_bp_C1_LTR <- 0
  outside_TE_bp_C1_nonLTR <- 0
  outside_TE_bp_C2_TIR <- 0
  outside_TE_bp_C2_other <- 0
  outside_TE_bp_other <- 0
  outside_TE_bp_Gypsy <- 0
  outside_TE_bp_Copia <- 0
  outside_TE_bp_Line <- 0
  outside_TE_bp_Sine <- 0
  outside_gene_bp_total <- 0
  
  
  
  for(chr_ID in 1 : nrow(genome_metadata)) {
    cat(satellite_metadata$TRASH_new_class[repeat_ID], genome_metadata$chromosome.name[chr_ID], chr_ID, "/", nrow(genome_metadata), "\n")
    
    arrays_cen_sat_chr <- arrays_cen_sat[arrays_cen_sat$seqID == genome_metadata$chromosome.name[chr_ID],]
    repeats_cen_sat_chr <- arrays_cen_sat_chr[arrays_cen_sat_chr$seqID == genome_metadata$chromosome.name[chr_ID],]
    
    repeats_all_chr <- repeats_all[repeats_all$seqID == genome_metadata$chromosome.name[chr_ID],]
    edta_all_chr <- edta[edta$V2 == genome_metadata$chromosome.name[chr_ID],]
    helixer_all_chr <- helixer[helixer$V1 == genome_metadata$chromosome.name[chr_ID],]
    
    gr_repeats <- GRanges(genome_metadata$chromosome.name[chr_ID], 
                          IRanges(start = repeats_all_chr$start,
                                  end   = repeats_all_chr$end))
    
    gr_edta <- GRanges(genome_metadata$chromosome.name[chr_ID], 
                       IRanges(start = edta_all_chr$V5,
                               end   = edta_all_chr$V6))
    
    gr_helixer <- GRanges(genome_metadata$chromosome.name[chr_ID], 
                          IRanges(start = helixer_all_chr$V4,
                                  end   = helixer_all_chr$V5))
    
    if(nrow(arrays_cen_sat_chr) == 0) {
      gr_within <- GRanges()
      gr_around <- GRanges()
      gr_outside <- GRanges(genome_metadata$chromosome.name[chr_ID], IRanges(1, genome_metadata$size[chr_ID]))
      gr_outside <- setdiff(gr_outside, gr_repeats)
    } else {
      
      # genomic ranges within arrays:
      arrays_cen_sat_chr$dist_to_next <- 999999999
      if(nrow(arrays_cen_sat_chr) > 1) arrays_cen_sat_chr$dist_to_next[-nrow(arrays_cen_sat_chr)] <- arrays_cen_sat_chr$start[-1] - arrays_cen_sat_chr$end[-nrow(arrays_cen_sat_chr)]
      arrays_cen_sat_chr$internal_gap_after <- FALSE
      arrays_cen_sat_chr$internal_gap_after[arrays_cen_sat_chr$dist_to_next < max_gap_size] <- TRUE
      if(sum(arrays_cen_sat_chr$internal_gap_after) > 0) {
        gr_within <- GRanges(genome_metadata$chromosome.name[chr_ID], 
                             IRanges(arrays_cen_sat_chr$end[arrays_cen_sat_chr$internal_gap_after] + 1, 
                                     arrays_cen_sat_chr$start[which(arrays_cen_sat_chr$internal_gap_after) + 1] - 1))
        gr_within <- setdiff(gr_within, gr_repeats)
      } else {
        gr_within <- GRanges()
      }
      
      
      # genomic ranges around arrays:
      around_df <- data.frame()
      for(j in seq_len(nrow(arrays_cen_sat_chr))) {
        
        if(j == 1) { # for the first one
          if(arrays_cen_sat_chr$start[1] >= max_gap_size) {
            around_df <- rbind(around_df, data.frame(start = (arrays_cen_sat_chr$start[1] - max_gap_size - 1), end = arrays_cen_sat_chr$start[1] - 1))
            around_df$start[around_df$start < 1] <- 1
          }
        } 
        
        if (i == nrow(arrays_cen_sat_chr)) { # for the last one
          if((genome_metadata$size[chr_ID] - arrays_cen_sat_chr$end[nrow(arrays_cen_sat_chr)]) >= max_gap_size) {
            around_df <- rbind(around_df, data.frame(start = arrays_cen_sat_chr$end[nrow(arrays_cen_sat_chr)] + 1, end = arrays_cen_sat_chr$end[nrow(arrays_cen_sat_chr)] + max_gap_size + 1))
            around_df$end[around_df$end >  genome_metadata$size[chr_ID]] <- genome_metadata$size[chr_ID]
          }
        } 
        
        if(nrow(arrays_cen_sat_chr) > 1) {
          if(j == nrow(arrays_cen_sat_chr)) next # done already
          
          if(arrays_cen_sat_chr$dist_to_next[j] > max_gap_size) { # gap between this and next is big enough to fit one or two "around" regions
            if(arrays_cen_sat_chr$dist_to_next[j] > (2*max_gap_size)) { # two regions can be fit
              around_df <- rbind(around_df, data.frame(start = arrays_cen_sat_chr$end[j]+1, end = arrays_cen_sat_chr$end[j] + max_gap_size + 1))
              around_df <- rbind(around_df, data.frame(start = arrays_cen_sat_chr$start[j + 1] - max_gap_size - 1, end = arrays_cen_sat_chr$start[j + 1] - 1))
            } else { # only one region fits
              around_df <- rbind(around_df, data.frame(start = arrays_cen_sat_chr$end[j]+1, end = arrays_cen_sat_chr$start[j + 1] - 1))
            }
          }
          
        }
      }
      
      if(nrow(around_df) != 0) {
        gr_around <- GRanges(genome_metadata$chromosome.name[chr_ID], 
                             IRanges(around_df$start, 
                                     around_df$end))
        gr_around <- setdiff(gr_around, gr_repeats)
      } else {
        gr_around <- GRanges()
      }
      
      
      # genomic ranges outside arrays:
      
      gr_all_repeats_around_within <- reduce(c(gr_within, gr_around, gr_repeats))
      
      gr_chromosome <- GRanges(genome_metadata$chromosome.name[chr_ID],
                               IRanges(start = 1, 
                                       end = genome_metadata$size[chr_ID]))
      
      gr_outside <- setdiff(gr_chromosome, gr_all_repeats_around_within) # all ranges that are not repeats, gr_within or gr_around
      
    }
    
    get_te_subset_ranges <- function(te_family_names, range_default) {
      edta_subset <- edta_all_chr[grep(paste0(te_family_names, collapse = "|"), edta_all_chr$V4), ]
      if(nrow(edta_subset) == 0) return(0)
      sum(width(intersect(GRanges(genome_metadata$chromosome.name[chr_ID], 
                                  IRanges(start = edta_subset$V5,
                                          end   = edta_subset$V6)), range_default)))
    }
    
    within_bp_total <- within_bp_total + sum(width(gr_within))
    within_TE_bp_total <- within_TE_bp_total + sum(width(intersect(gr_edta, gr_within)))
    within_TE_bp_C1_LTR <- within_TE_bp_C1_LTR + get_te_subset_ranges(unlist(edta_classes[1:6]), gr_within)
    within_TE_bp_C1_nonLTR <- within_TE_bp_C1_nonLTR + get_te_subset_ranges(unlist(edta_classes[7:11]), gr_within)
    within_TE_bp_C2_TIR <- within_TE_bp_C2_TIR + get_te_subset_ranges(unlist(edta_classes[12]), gr_within)
    within_TE_bp_C2_other <- within_TE_bp_C2_other + get_te_subset_ranges(unlist(edta_classes[13:16]), gr_within)
    within_TE_bp_other <- within_TE_bp_other + get_te_subset_ranges(unlist(edta_classes[17:18]), gr_within)
    within_TE_bp_Gypsy <- within_TE_bp_Gypsy + get_te_subset_ranges(unlist(edta_classes[1]), gr_within)
    within_TE_bp_Copia <- within_TE_bp_Copia + get_te_subset_ranges(unlist(edta_classes[2]), gr_within)
    within_TE_bp_Line <- within_TE_bp_Line + get_te_subset_ranges(unlist(edta_classes[7]), gr_within)
    within_TE_bp_Sine <- within_TE_bp_Sine + get_te_subset_ranges(unlist(edta_classes[8]), gr_within)
    within_gene_bp_total <- within_gene_bp_total + sum(width(intersect(gr_helixer, gr_within)))
    
    
    around_bp_total <- around_bp_total + sum(width(gr_around))
    around_TE_bp_total <- around_TE_bp_total + sum(width(intersect(gr_edta, gr_around)))
    around_TE_bp_C1_LTR <- around_TE_bp_C1_LTR + get_te_subset_ranges(unlist(edta_classes[1:6]), gr_around)
    around_TE_bp_C1_nonLTR <- around_TE_bp_C1_nonLTR + get_te_subset_ranges(unlist(edta_classes[7:11]), gr_around)
    around_TE_bp_C2_TIR <- around_TE_bp_C2_TIR + get_te_subset_ranges(unlist(edta_classes[12]), gr_around)
    around_TE_bp_C2_other <- around_TE_bp_C2_other + get_te_subset_ranges(unlist(edta_classes[13:16]), gr_around)
    around_TE_bp_other <- around_TE_bp_other + get_te_subset_ranges(unlist(edta_classes[17:18]), gr_around)
    around_TE_bp_Gypsy <- around_TE_bp_Gypsy + get_te_subset_ranges(unlist(edta_classes[1]), gr_around)
    around_TE_bp_Copia <- around_TE_bp_Copia + get_te_subset_ranges(unlist(edta_classes[2]), gr_around)
    around_TE_bp_Line <- around_TE_bp_Line + get_te_subset_ranges(unlist(edta_classes[7]), gr_around)
    around_TE_bp_Sine <- around_TE_bp_Sine + get_te_subset_ranges(unlist(edta_classes[8]), gr_around)
    around_gene_bp_total <- around_gene_bp_total + sum(width(intersect(gr_helixer, gr_around)))
    
    
    outside_bp_total <- outside_bp_total + sum(width(gr_outside))
    outside_TE_bp_total <- outside_TE_bp_total + sum(width(intersect(gr_edta, gr_outside)))
    outside_TE_bp_C1_LTR <- outside_TE_bp_C1_LTR + get_te_subset_ranges(unlist(edta_classes[1:6]), gr_outside)
    outside_TE_bp_C1_nonLTR <- outside_TE_bp_C1_nonLTR + get_te_subset_ranges(unlist(edta_classes[7:11]), gr_outside)
    outside_TE_bp_C2_TIR <- outside_TE_bp_C2_TIR + get_te_subset_ranges(unlist(edta_classes[12]), gr_outside)
    outside_TE_bp_C2_other <- outside_TE_bp_C2_other + get_te_subset_ranges(unlist(edta_classes[13:16]), gr_outside)
    outside_TE_bp_other <- outside_TE_bp_other + get_te_subset_ranges(unlist(edta_classes[17:18]), gr_outside)
    outside_TE_bp_Gypsy <- outside_TE_bp_Gypsy + get_te_subset_ranges(unlist(edta_classes[1]), gr_outside)
    outside_TE_bp_Copia <- outside_TE_bp_Copia + get_te_subset_ranges(unlist(edta_classes[2]), gr_outside)
    outside_TE_bp_Line <- outside_TE_bp_Line + get_te_subset_ranges(unlist(edta_classes[7]), gr_outside)
    outside_TE_bp_Sine <- outside_TE_bp_Sine + get_te_subset_ranges(unlist(edta_classes[8]), gr_outside)
    outside_gene_bp_total <- outside_gene_bp_total + sum(width(intersect(gr_helixer, gr_outside)))
    
    
    
  }
  
  enrichment_df <- data.frame(count_of = c("all_analysed_regions", "all_TE", "TE_C1_LTR", "TE_C1_nonLTR", 
                                           "TE_C2_TIR", "TE_C2_nonTIR", "TE_other", "TE_Gypsy", "TE_Copa", 
                                           "TE_Line", "TE_Sine", "all_gene"),
                              bp_count_within_cen_array = c(within_bp_total, within_TE_bp_total, within_TE_bp_C1_LTR, within_TE_bp_C1_nonLTR, 
                                                            within_TE_bp_C2_TIR, within_TE_bp_C2_other, within_TE_bp_other, within_TE_bp_Gypsy, within_TE_bp_Copia,
                                                            within_TE_bp_Line, within_TE_bp_Sine, within_gene_bp_total),
                              bp_count_around_cen_array = c(around_bp_total, around_TE_bp_total, around_TE_bp_C1_LTR, around_TE_bp_C1_nonLTR, 
                                                            around_TE_bp_C2_TIR, around_TE_bp_C2_other, around_TE_bp_other, around_TE_bp_Gypsy, around_TE_bp_Copia,
                                                            around_TE_bp_Line, around_TE_bp_Sine, around_gene_bp_total),
                              bp_count_outside_cen_array = c(outside_bp_total, outside_TE_bp_total, outside_TE_bp_C1_LTR, outside_TE_bp_C1_nonLTR, 
                                                            outside_TE_bp_C2_TIR, outside_TE_bp_C2_other, outside_TE_bp_other, outside_TE_bp_Gypsy, outside_TE_bp_Copia,
                                                            outside_TE_bp_Line, outside_TE_bp_Sine, outside_gene_bp_total))
  
  enrichment_df$fraction_within_cen_array <- 100 * enrichment_df$bp_count_within_cen_array / enrichment_df$bp_count_within_cen_array[1]
  enrichment_df$fraction_around_cen_array <- 100 * enrichment_df$bp_count_around_cen_array / enrichment_df$bp_count_around_cen_array[1]
  enrichment_df$fraction_outside_cen_array <- 100 * enrichment_df$bp_count_outside_cen_array / enrichment_df$bp_count_outside_cen_array[1]
  
  enrichment_df$normalised_fraction_within_cen_array <- enrichment_df$fraction_within_cen_array / (enrichment_df$fraction_outside_cen_array + 0.000000001)
  enrichment_df$normalised_fraction_around_cen_array <- enrichment_df$fraction_around_cen_array / (enrichment_df$fraction_outside_cen_array + 0.000000001)
  
  enrichment_df$normalised_log10_fraction_within_cen_array <- log10(enrichment_df$normalised_fraction_within_cen_array)
  enrichment_df$normalised_log10_fraction_around_cen_array <- log10(enrichment_df$normalised_fraction_around_cen_array)
  
  write.csv(x = enrichment_df, file = paste0("/home/pwlodzimierz/ToL/upload_files/58_TE_gene_repeat_enrichment/data_files/", assembly_name, ";", satellite_metadata$TRASH_new_class[repeat_ID], ".csv"))
  
  cat("data frame saved:",  paste0("/home/pwlodzimierz/ToL/upload_files/58_TE_gene_repeat_enrichment/data_files/", assembly_name, ";", satellite_metadata$TRASH_new_class[repeat_ID], ".csv\n"))
}

































