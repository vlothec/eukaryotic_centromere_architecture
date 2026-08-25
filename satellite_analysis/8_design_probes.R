#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))

# here is an updated list of the species available for wet-lab experiments 
# (we have them in cultivation, and/or collected their seeds):
#
# - Ailanthus altissima
# - Ballota nigra
# - Chamaenerion angustifolium
# - Filipendula ulmaria
# - Geum urbanum
# - Linaria vulgaris
# - Solanum dulcamara
# 
# There are species from you updated list(s) that can easily be added:
#   
# - Schoenoplectus lacustris (grows in a pool located 50 m from my lab ;-))
# - Malus domestica
# - Hedera helix
# - Quercus robur
# - Alnus glutinosa
# - and maybe some other

avail_species <- c("drAilAlti1.1.fa","drHedHeli1.1.fa", "drChaAngu1.1.fa", "lpSchLacu1.1.fa","daBalNigr1.1.fa", 
                   "dhQueRobu3.1.fa", "daSolDulc1.1.fa", 
                   "drMalDome5.1.fa", "dhAlnGlut1.1.fa",
                   "drGeuUrba1.1.fa", "daLinVulg1.1.fa", "drFilUlma1.1.fa")

repeat_classes <- list(c("42_16"),c("160_1"),c("215_3"),c("183_4"),c("86_2"),
                       c("146_1"),c("183_1"),
                       c("345_5"),c(),
                       c("169_36"),c("266_11"),c("177_2"))

setwd("/home/pwlodzimierz/ToL/fish_probes")

metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.csv")

genomes <- unique(metadata$assembly.name)

avail_species %in% genomes

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

# for(i in seq_along(avail_species)) {
{
  cat(avail_species[i], "\n")
  repeats <- read.csv(file = paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_output/", 
                                    avail_species[i], "/", 
                                    strsplit(avail_species[i], split = ".fa")[[1]][1], 
                                    "_repeats_filtered.csv"))
  for(j in seq_along(repeat_classes[[i]])) {
    repeats_class <- repeats[repeats$new_class == repeat_classes[[i]][j],]
    
    # # extend the sequence
    # for(k in 1 : (nrow(repeats_class) - 1)) {
    #   if(repeats_class$start[k + 1] - repeats_class$end[k] < 10) {
    #     repeats_class$sequence[k] <- paste0(repeats_class$sequence[k], repeats_class$sequence[k+1])
    #   }
    # }
    cat("writing fasta", "\n")
    
    seq_to_align <- 1 : nrow(repeats_class)
    
    if(length(seq_to_align) > 5000) {
      seq_to_align <- sample(seq_to_align, 5000)
    }
    
    write.fasta(sequences = as.list(repeats_class$sequence), 
                names = 1:nrow(repeats_class), file.out = paste0("/home/pwlodzimierz/ToL/fish_probes/temp_fasta_to_align", i, ".fasta"))
    
    mafft_bat <- "mafft"
    
    cat("aligning", "\n")
    system2(command = mafft_bat , args = paste0(" --inputorder --kimura 1 --large --retree 1 ", 
                                                paste0("/home/pwlodzimierz/ToL/fish_probes/temp_fasta_to_align", i, ".fasta"),
                                                " > ", 
                                                paste0("/home/pwlodzimierz/ToL/fish_probes/temp_fasta_aligned", i,".fasta")))
    cat("reading alignment and generating frequency table", "\n")
    
    alignment <- read.fasta(file = paste0("/home/pwlodzimierz/ToL/fish_probes/temp_fasta_aligned", i, ".fasta"), 
                            seqonly = FALSE, as.string = FALSE)
    
    freq_table <- data.frame(a = rep(0,length(alignment[[1]])), 
                             c = rep(0,length(alignment[[1]])), 
                             g = rep(0,length(alignment[[1]])), 
                             t = rep(0,length(alignment[[1]])), 
                             none = rep(0,length(alignment[[1]])))
    for(l in seq_along(alignment)) {
      freq_table$a[alignment[[l]] %in% "a"] <- freq_table$a[alignment[[l]] %in% "a"] + 1
      freq_table$c[alignment[[l]] %in% "c"] <- freq_table$c[alignment[[l]] %in% "c"] + 1
      freq_table$t[alignment[[l]] %in% "t"] <- freq_table$t[alignment[[l]] %in% "t"] + 1
      freq_table$g[alignment[[l]] %in% "g"] <- freq_table$g[alignment[[l]] %in% "g"] + 1
      freq_table$none[alignment[[l]] %in% "-"] <- freq_table$none[alignment[[l]] %in% "-"] + 1
      
    }
    freq_table <- freq_table / length(alignment) * 100
    
    freq_table <- freq_table[freq_table$none < 35,]
    
    freq_table_og <- freq_table
    
    consensus_options <- NULL
    consensus_scores <- NULL
    consensus_sizes <- NULL
    
    min_size <- min(35, round(1*mean(repeats_class$width)))
    max_size <- min(135, round(1*mean(repeats_class$width)))
    # min_size <- min(35, round(2*mean(repeats_class$width)))
    # max_size <- min(135, round(2*mean(repeats_class$width)))
    
    cat("generating candidate probes", "\n")
    for(probe_size in min_size : max_size) {
      cat(i, j, "probe size", probe_size, "\n")
      
      consensus_temp <- NULL
      max_score_temp <- 0
      which_max_score_temp <- 0
      for(l in 1 : nrow(freq_table)) {
        max_n <- (l + probe_size - 1)
        if(max_n > nrow(freq_table)) max_n <- nrow(freq_table)
        score_temp = 0
        for(m in l : max_n) {
          score_temp <- score_temp + max(freq_table[m,])
        }
        if(score_temp > max_score_temp) {
          max_score_temp <- score_temp
          which_max_score_temp <- l
        }
      }
      max_n <- (which_max_score_temp + probe_size - 1)
      if(max_n > nrow(freq_table)) max_n <- nrow(freq_table)
      for(m in which_max_score_temp : max_n) {
        consensus_temp <- c(consensus_temp, c("a","c","g","t","-")[which.max(freq_table[m,])])
      }
      consensus_options <- c(consensus_options, paste(consensus_temp, collapse = ""))
      consensus_scores <- c(consensus_scores, max_score_temp)
      consensus_sizes <- c(consensus_sizes, probe_size)
      
      
      
      second_tops <- NULL
      second_tops_freq <- NULL
      for(l in 1 : nrow(freq_table)) {
        second_tops <- c(second_tops, c("a","c","g","t","-")[order(as.matrix(freq_table[l,(1:4)]))[3]])
        second_tops_freq <- c(second_tops_freq, freq_table[l,(1:4)][order(as.matrix(freq_table[l,(1:4)]))[3]])
      }
      second_tops_freq <- as.numeric(second_tops_freq)
      
      for(alt_cons in 1 : nrow(freq_table)) {
        freq_table <- freq_table_og
        if(second_tops_freq[alt_cons] < 40) next
        freq_table[alt_cons,][which.max(freq_table[alt_cons,])] = 0
        
        consensus_temp <- NULL
        max_score_temp <- 0
        which_max_score_temp <- 0
        for(l in 1 : nrow(freq_table)) {
          max_n <- (l + probe_size - 1)
          if(max_n > nrow(freq_table)) max_n <- nrow(freq_table)
          score_temp = 0
          for(m in l : max_n) {
            score_temp <- score_temp + max(freq_table[m,])
          }
          if(score_temp > max_score_temp) {
            max_score_temp <- score_temp
            which_max_score_temp <- l
          }
        }
        max_n <- (which_max_score_temp + probe_size - 1)
        if(max_n > nrow(freq_table)) max_n <- nrow(freq_table)
        for(m in which_max_score_temp : max_n) {
          consensus_temp <- c(consensus_temp, c("a","c","g","t","-")[which.max(freq_table[m,])])
        }
        consensus_options <- c(consensus_options, paste(consensus_temp, collapse = ""))
        consensus_scores <- c(consensus_scores, max_score_temp)
        consensus_sizes <- c(consensus_sizes, probe_size)
      }
      
    }
    ### deduplicate just in case
    cat("deduplicating", "\n")
    
    consensus_scores <- consensus_scores[order(consensus_options)]
    consensus_sizes <- consensus_sizes[order(consensus_options)]
    consensus_options <- consensus_options[order(consensus_options)]
    
    consensus_scoresb = consensus_scores
    consensus_sizesb=consensus_sizes
    consensus_optionsb=consensus_options
    consensus_scores=consensus_scoresb
    consensus_sizes=consensus_sizesb
    consensus_options=consensus_optionsb
    
    deduplicate_ids <- NULL
    for(k in seq_along(consensus_options)[-1]) {
      if(consensus_options[k-1] == consensus_options[k]) {
        deduplicate_ids <- c(deduplicate_ids, k)
      }
    }
    if(length(deduplicate_ids) > 0) {
      consensus_scores <- consensus_scores[-deduplicate_ids]
      consensus_sizes <- consensus_sizes[-deduplicate_ids]
      consensus_options <- consensus_options[-deduplicate_ids]
    }
    
    ### get size normalised scores, push towards longer probes
    size_normalised_consensus_scores <- consensus_scores / (consensus_sizes)
    
    consensus_sizes <- consensus_sizes[order(size_normalised_consensus_scores, decreasing = TRUE)]
    consensus_options <- consensus_options[order(size_normalised_consensus_scores, decreasing = TRUE)]
    consensus_scores <- consensus_scores[order(size_normalised_consensus_scores, decreasing = TRUE)]
    size_normalised_consensus_scores <- size_normalised_consensus_scores[order(size_normalised_consensus_scores, decreasing = TRUE)]
    
    consensus_scores <- consensus_scores / consensus_sizes
    
    # use no more than 100 top hits
    if(length(consensus_sizes) > 100) {
      consensus_sizes <- consensus_sizes[1:100]
      consensus_options <- consensus_options[1:100]
      consensus_scores <- consensus_scores[1:100]
      size_normalised_consensus_scores <- size_normalised_consensus_scores[1:100]
    }
    
    
    ### score the % of repeats each probe hits with 0, 1, 2, 3, 5, 10 allowed mismatches
    ### score the % of repeats per genome and per chromosome it's present on
    ### plot the hits and scores
    cat("plotting", "\n")
    
    chromosomes <- unique(repeats_class$seqID)
    
    counts <- table(repeats_class$seqID)
    
    pdf(file = paste0(strsplit(avail_species[i], split = ".fa")[[1]][1], "_", repeat_classes[[i]][j], ".pdf"), 
        width = 10, height = 4 * length(chromosomes))
    
    par(mfrow = c(length(chromosomes),2), oma = c(1,1,1,1), mar = c(4,4,4,1))
    nf <- layout(
      matrix(c(1:(2*length(chromosomes))), ncol=2, byrow=TRUE), 
      widths=c(4,1)
    )
    
    cat(counts, "\n")
    
    for(k in seq_along(chromosomes)) {
      cat(i, "/", length(avail_species), 
          j, "/", length(repeat_classes[[i]]), 
          "plotting", k, "/", length(chromosomes), "\n")
      chromosome_repeats <- repeats_class[repeats_class$seqID == chromosomes[k], ]
      chromosome_repeats$start <- chromosome_repeats$start / 1000000
      plot(x = chromosome_repeats$start, y = rep(-5, length(chromosome_repeats$start)), 
           xlim = c(min(chromosome_repeats$start), max(chromosome_repeats$start)), 
           ylim = c(-5,length(consensus_options)), col = "blue", pch = 18,
           xlab = "Coordinates, Mbp", ylab = "repeat markers and probe hits")
      
      abline(h = seq(0.5,100.5, by = 10), lwd = 0.5, col = "#000000")
      
      hitsall <- adist(x = consensus_options, chromosome_repeats$sequence, costs = list(ins = 1, del = 1, sub = 1))
      
      
      for(l in seq_along(consensus_options)) {
        hits <- hitsall[l,]
        
        points(x = chromosome_repeats$start[which(hits <= 25)], y = rep(l, length(which(hits <= 25))), pch = 16, col = "#bbbbbb", cex = 0.5)
        points(x = chromosome_repeats$start[which(hits <= 10)], y = rep(l, length(which(hits <= 10))), pch = 16, col = "#999999", cex = 0.5)
        points(x = chromosome_repeats$start[which(hits <= 5)], y = rep(l, length(which(hits <= 5))), pch = 16, col = "#666666", cex = 0.5)
        points(x = chromosome_repeats$start[which(hits <= 2)], y = rep(l, length(which(hits <= 2))), pch = 16, col = "#444444", cex = 0.5)
        points(x = chromosome_repeats$start[which(hits <= 1)], y = rep(l, length(which(hits <= 1))), pch = 16, col = "#222222", cex = 0.5)
        points(x = chromosome_repeats$start[which(hits == 0)], y = rep(l, length(which(hits == 0))), pch = 16, col = "#000000", cex = 0.5)
      }
      
      
      
      plot(y=NULL,x=NULL,xlab="",ylab="", xlim = c(0,10), 
           ylim = c(-5,length(consensus_options)), 
           xaxt = "n", yaxt = "n", bty="n")
      text(1,-5, paste0(nrow(chromosome_repeats), " total repeats of the class on this chromosome"), cex = 0.5, pos = 4)
      for(l in seq_along(consensus_options)) {
        hits <- hitsall[l,]
        
        text(1,l, paste("len: ", consensus_sizes[l], " ", 
                        sum(hits <= 25), "<=25 ", 
                        sum(hits <= 10), "<=10 ", 
                        sum(hits <= 5), "<=5 ", 
                        sum(hits <= 2), "<=2 ", 
                        sum(hits <= 1), "<=1 ", 
                        sum(hits == 0), "==0", sep = ""), cex = 0.35, pos = 4)
      }
      
    }
    dev.off()
    
    remove(alignment)
  } # j: each repeat class of the genome
  remove(repeats)
  remove(repeats_class)
} # i: each genome












