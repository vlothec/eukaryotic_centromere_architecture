



library(ggplot2)
library(patchwork)

data_dir <- "/home/pwlodzimierz/ToL/upload_files/58_TE_gene_repeat_enrichment/data_files/" 
files <- sort(list.files(data_dir, pattern = "\\.csv$", full.names = TRUE))

clamp  <- function(x, lo = -3, hi = 3) pmax(pmin(x, hi), lo)
is_oob <- function(x, lo = -3, hi = 3) x < lo | x > hi

rows_per_page <- 15
n_files  <- length(files)
n_panels <- NULL

# build one row of panels per file
row_panels <- vector("list", n_files)

for (fi in seq_along(files)) {
  cat(fi, "\n")
  df <- read.csv(files[fi])
  df <- df[-1,]
  row_label <- tools::file_path_sans_ext(basename(files[fi]))
  n_panels  <- nrow(df)
  
  panels <- vector("list", n_panels)
  for (i in 1:n_panels) {
    xw <- df$normalised_log10_fraction_within_cen_array[i]
    xa <- df$normalised_log10_fraction_around_cen_array[i]
    
    pts <- data.frame(
      x   = clamp(c(xw, xa)),
      grp = c("within", "around"),
      oob = is_oob(c(xw, xa))
    )
    pts$colour <- with(pts, ifelse(grp == "within",
                                   ifelse(oob, "#0066ff", "#6baed6"),
                                   ifelse(oob, "#ff2200", "#fc8d59")))
    
    panels[[i]] <- ggplot(pts, aes(x = x, colour = colour)) +
      geom_segment(x = -3, xend = 3, y = 0, yend = 0, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
      geom_segment(x = 0,  xend = 0, y = -1, yend = 1, color = "black", linewidth = 0.6, inherit.aes = FALSE) +
      geom_segment(aes(xend = x, y = -0.5, yend = 0.5), linewidth = 1.2, alpha = 0.6) +
      scale_color_identity() +
      scale_x_continuous(limits = c(-3, 3), breaks = -3:3,
                         labels = c("0.001","0.01","0.1","1","10","100","1000")) +
      scale_y_continuous(limits = c(-1, 1)) +
      labs(title = NULL, x = NULL, y = if (i == 1) row_label else NULL) +
      theme_minimal(base_size = 7) +
      theme(
        axis.line.x  = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid   = element_blank(),
        plot.title   = element_text(size = 12, hjust = 0.5),
        axis.title.y = element_text(size = 12, angle = 0, vjust = 0.5,
                                    margin = margin(r = 10)),
        aspect.ratio = 0.15,
        plot.margin  = unit(c(1, 1, 1, 8), "mm")
      )
  }
  row_panels[[fi]] <- panels
}

# add column titles to first row of first page, tick labels to last row of each page
add_titles   <- function(panels, df) {
  for (i in seq_along(panels))
    panels[[i]] <- panels[[i]] +
      labs(title = df$count_of[i]) +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
  panels
}
add_ticklabs <- function(panels) {
  for (i in seq_along(panels))
    panels[[i]] <- panels[[i]] +
      theme(axis.text.x  = element_text(size = 12, angle = 0, hjust = 0.5),
            axis.ticks.x = element_line(color = "black"))
  panels
}

# read first file's count_of for titles
titles_df <- read.csv(files[1])

pdf("cen_enrichment_plot.pdf", width = n_panels * 2, height = rows_per_page * 0.5)

page_starts <- seq(1, n_files, by = rows_per_page)
for (ps in page_starts) {
  cat(ps, "\n")
  pe        <- min(ps + rows_per_page - 1, n_files)
  page_rows <- row_panels[ps:pe]
  n_page    <- length(page_rows)
  
  # titles on first row only if this is the first page
  if (ps == 1) page_rows[[1]] <- add_titles(page_rows[[1]], titles_df)
  
  # tick labels on last row of every page
  page_rows[[n_page]] <- add_ticklabs(page_rows[[n_page]])
  
  all_panels <- unlist(page_rows, recursive = FALSE)
  page_plot  <- wrap_plots(all_panels, nrow = n_page, ncol = n_panels)
  print(page_plot)
}

dev.off()
message("Saved: cen_enrichment_plot.pdf")