bins = 30


ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}


####################################
# taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# i = as.numeric(taskid)# 1 to 193
# print(i)
setwd("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs") 
hor_dirs <- list.dirs(path = ".", recursive = FALSE)



for(i in 1 : length(hor_dirs)) {
  cat(i, "\n")
  
  
  
  assembly <- strsplit(hor_dirs[i], split = "/")[[1]][2]
  
  
  setwd("/home/pwlodzimierz/ToL/upload_files/9_HOR_periods/repeat_files_with_HORs") 
  setwd(hor_dirs[i])
  
  repeats_files <- list.files(path = ".", pattern = "HOR_scored_repeats_with_hors")
  
  total_HOR_score <- rep(0, bins)
  total_repeats_per_bin <- rep(0, bins)
  
  for(rep_file in repeats_files) {
    if(file.size(rep_file) == 0) next
    
    chr_repeats <- read.csv(file = rep_file)
    
    starts <- round(seq(0 , nrow(chr_repeats), length.out = bins + 1))
    ends <- starts[-1]
    starts <- starts[-length(starts)] + 1
    starts[starts > nrow(chr_repeats)] <- nrow(chr_repeats)
    ends[ends > nrow(chr_repeats)] <- nrow(chr_repeats)
    
    chr_repeats$HOR_score[is.na(chr_repeats$HOR_score)] = 0
    
    n <- nrow(chr_repeats)
    if (n == 0) next
    b <- ceiling(seq_len(n) * bins / n)
    total_HOR_score <- total_HOR_score +
      as.vector(tapply(chr_repeats$HOR_score, factor(b, 1:bins), sum, default = 0))
    total_repeats_per_bin <- total_repeats_per_bin + tabulate(b, bins)
    if(sum(is.na(total_HOR_score)) > 0) stop()
    
    
  }
  if(sum(total_repeats_per_bin) == 0) {
    normalised_HOR_score <- rep(0, bins)
  } else {
    normalised_HOR_score <- total_HOR_score / (total_repeats_per_bin + 0)
    
    normalised_HOR_score[is.na(normalised_HOR_score)] <- 0
    
    normalised_HOR_score <- ma(c(normalised_HOR_score[1], normalised_HOR_score[1], normalised_HOR_score, normalised_HOR_score[length(normalised_HOR_score)], normalised_HOR_score[length(normalised_HOR_score)]))
    normalised_HOR_score <- normalised_HOR_score[!is.na(normalised_HOR_score)]
    
  }
  
  
  # pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/15_HOR_score_along/hor_score_along_", assembly, ".pdf"), 
  #     width = 7, height = 4)
  # 
  # plot(1 : bins, normalised_HOR_score, type = "l", lwd = 2, 
  #      ylab = "HOR score, %", xlab = paste0(bins, " bins"), main = paste0(" hor score along ", assembly, " chromosomes"))
  # 
  # dev.off()
  # 
  # 
  # pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/15_HOR_score_along/plots_same_ylim/hor_score_along_", assembly, ".pdf"), 
  #     width = 7, height = 4)
  # 
  # plot(1 : bins, normalised_HOR_score, type = "l", lwd = 2, 
  #      ylab = "HOR score, %", xlab = paste0(bins, " bins"), main = paste0(" hor score along ", assembly, " chromosomes"),
  #      ylim = c(0,8))
  # 
  # dev.off()
  
  
  write.csv(x = normalised_HOR_score, file = paste0("/home/pwlodzimierz/ToL/upload_files/15_HOR_score_along/scores_csvs/scores_csvs_", assembly, ".csv"), row.names = FALSE)
  
  
  
}






















































































































