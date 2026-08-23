#!/usr/bin/env Rscript
# plot_entropy_325sp.R
# Publication-ready Shannon entropy figures for the 325-sp dataset.
#
# Panel C (positional CenpA entropy by group):
#   Reads cenpa_entropy_split_phylo.tsv from the 327-sp entropia-split folder.
#   The group-level TSV is not affected by removing 2 species.
#
# Panel D companion (CenpA vs H3 split entropy):
#   Reads split_entropy_bnni_gap085.tsv from 04_cenpa_phylogeny/split_entropy/.
#   This was computed on the 325-sp bnni tree.
#
# Outputs → 03_entropy/figures/ and top-level figures/

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
ENT_DIR <- file.path(BASE, "PhylogeneticProfiling/20_cenpa_analysis/20.1-march/entropia-split")
SP_DIR  <- file.path(PUB_DIR, "04_cenpa_phylogeny/split_entropy")
FIG_DIR <- file.path(PUB_DIR, "03_entropy/figures")
TOP_FIG <- file.path(PUB_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TOP_FIG, showWarnings = FALSE, recursive = TRUE)

groups_pal <- c(
  Fungi         = "#6A1B9A",
  Invertebrates = "#E65100",
  Vertebrates   = "#1565C0",
  Viridiplantae = "#2E7D32"
)

# ── Panel C: positional CenpA entropy by taxonomic group ─────────────────────
ent <- read.delim(file.path(ENT_DIR, "cenpa_entropy_split_phylo.tsv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("CenpA group entropy: %d positions, groups: %s\n",
            nrow(ent), paste(setdiff(colnames(ent), "pos"), collapse = ", ")))

ent_long <- ent %>%
  pivot_longer(-pos, names_to = "group", values_to = "entropy") %>%
  mutate(group = factor(group, levels = names(groups_pal)))

# Raw line plot
p_c_raw <- ggplot(ent_long, aes(x = pos, y = entropy, colour = group)) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  scale_colour_manual(values = groups_pal, name = "Taxonomic group") +
  scale_x_continuous(name = "Alignment position (trimmed, gap <= 85%)") +
  scale_y_continuous(name = "Normalised Shannon entropy",
                     limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "CENP-A/CENH3 positional entropy by taxonomic group (325-sp)",
    subtitle = sprintf("428 CenpA sequences; %d trimmed alignment positions", max(ent$pos))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

# Smoothed (10-position rolling average)
roll_mean <- function(x, k = 10) as.numeric(stats::filter(x, rep(1/k, k), sides = 2))

ent_smooth <- ent_long %>%
  group_by(group) %>%
  arrange(pos) %>%
  mutate(entropy_sm = roll_mean(entropy, k = 10)) %>%
  ungroup()

p_c_smooth <- ggplot(ent_smooth, aes(x = pos, y = entropy_sm, colour = group)) +
  geom_line(linewidth = 0.65, alpha = 0.90, na.rm = TRUE) +
  scale_colour_manual(values = groups_pal, name = "Taxonomic group") +
  scale_x_continuous(name = "Alignment position (trimmed, gap <= 85%)") +
  scale_y_continuous(name = "Normalised Shannon entropy (10-pos rolling avg)",
                     limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Smoothed positional entropy (window = 10 positions)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

# Combined Panel C
p_c <- p_c_raw / p_c_smooth +
  plot_annotation(
    title = "Figure 1C — CENP-A/CENH3 positional Shannon entropy (325-sp)",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(FIG_DIR, "panel_C_entropy_raw_325sp.pdf"),    p_c_raw,    width = 14, height = 4.5)
ggsave(file.path(FIG_DIR, "panel_C_entropy_raw_325sp.png"),    p_c_raw,    width = 14, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "panel_C_entropy_smooth_325sp.pdf"), p_c_smooth, width = 14, height = 4.5)
ggsave(file.path(FIG_DIR, "panel_C_entropy_smooth_325sp.png"), p_c_smooth, width = 14, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "panel_C_entropy_combined_325sp.pdf"), p_c, width = 14, height = 8)
ggsave(file.path(FIG_DIR, "panel_C_entropy_combined_325sp.png"), p_c, width = 14, height = 8, dpi = 300)
cat("Saved: panel_C_entropy (raw + smoothed + combined)\n")

# Copy to top-level figures
file.copy(file.path(FIG_DIR, "panel_C_entropy_raw_325sp.pdf"),
          file.path(TOP_FIG, "panel_C_entropy_raw_325sp.pdf"), overwrite = TRUE)
file.copy(file.path(FIG_DIR, "panel_C_entropy_raw_325sp.png"),
          file.path(TOP_FIG, "panel_C_entropy_raw_325sp.png"), overwrite = TRUE)

# ── Panel D: CenpA vs H3 split entropy with secondary structure + masking ──────
# Uses per-group masking (NaN where group's own gap > 85%) — same as Python pergroup plot
se_pg <- read.delim(file.path(SP_DIR, "split_entropy_bnni_pergroup_gap085.tsv"),
                    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "nan")
cat(sprintf("CenpA vs H3 split entropy (per-group masked): %d positions\n", nrow(se_pg)))
cat(sprintf("  Masked positions — CENPA: %d  H3: %d\n",
            sum(is.na(se_pg$CENPA)), sum(is.na(se_pg$H3))))

# Helix positions (projected to gap<=85% trimmed coordinates)
helices <- read.delim(file.path(SP_DIR, "helix_positions_gap085.tsv"),
                      stringsAsFactors = FALSE)
cat("Alpha helices:\n")
print(helices)

histone_pal <- c(CENPA = "#C62828", H3 = "#1565C0")

# Masked regions: convert NaN positions to rectangles
mask_cenpa <- se_pg$pos[is.na(se_pg$CENPA)]
mask_h3    <- se_pg$pos[is.na(se_pg$H3)]

# Top 10% CenpA positions (above 90th percentile, ignoring NA)
cenpa_thresh <- quantile(se_pg$CENPA, 0.90, na.rm = TRUE)
top_cenpa    <- se_pg[!is.na(se_pg$CENPA) & se_pg$CENPA >= cenpa_thresh, ]

p_d_split <- ggplot(se_pg, aes(x = pos)) +
  # Masked regions (light shading)
  {if (length(mask_h3) > 0)
    annotate("rect",
             xmin = mask_h3 - 0.5, xmax = mask_h3 + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "#1565C0", alpha = 0.07)
  } +
  {if (length(mask_cenpa) > 0)
    annotate("rect",
             xmin = mask_cenpa - 0.5, xmax = mask_cenpa + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "#C62828", alpha = 0.07)
  } +
  # Entropy lines
  geom_line(aes(y = H3,    colour = "H3"),    linewidth = 0.7, alpha = 0.85, na.rm = TRUE) +
  geom_line(aes(y = CENPA, colour = "CENPA"), linewidth = 0.8, alpha = 0.90, na.rm = TRUE) +
  # Top 10% CenpA dots
  geom_point(data = top_cenpa, aes(y = CENPA),
             colour = "#7B0000", size = 1.2, alpha = 0.75, shape = 19) +
  # Alpha-helix bands as data-driven rects (two strips at top of y-axis)
  geom_rect(
    data = transform(helices[helices$histone == "CENPA", ],
                     ymin = 0.935, ymax = 0.990),
    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
    fill = "#C62828", alpha = 0.45, colour = NA, inherit.aes = FALSE
  ) +
  geom_rect(
    data = transform(helices[helices$histone == "H3", ],
                     ymin = 0.875, ymax = 0.930),
    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
    fill = "#1565C0", alpha = 0.45, colour = NA, inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = histone_pal,
    labels = c(CENPA = "CENP-A (n=418)", H3 = "H3-like (n=901)"),
    name   = NULL
  ) +
  scale_x_continuous(name = "Alignment position (per-group mask, gap <= 85%)") +
  scale_y_continuous(name = "Normalised Shannon entropy",
                     limits = c(0, 1.01), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  labs(
    title    = "CENP-A vs H3 split entropy — bnni tree (325-sp)",
    subtitle = paste0(
      "Per-group masking: NaN where group gap > 85% (CENPA: ", length(mask_cenpa),
      " masked pos; H3: ", length(mask_h3), " masked pos).  ",
      "Top dark dots = top 10% CENP-A entropy positions."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(FIG_DIR, "panel_D_split_entropy_325sp.pdf"), p_d_split, width = 13, height = 4.5)
ggsave(file.path(FIG_DIR, "panel_D_split_entropy_325sp.png"), p_d_split, width = 13, height = 4.5, dpi = 300)
file.copy(file.path(FIG_DIR, "panel_D_split_entropy_325sp.pdf"),
          file.path(TOP_FIG, "panel_D_split_entropy_325sp.pdf"), overwrite = TRUE)
file.copy(file.path(FIG_DIR, "panel_D_split_entropy_325sp.png"),
          file.path(TOP_FIG, "panel_D_split_entropy_325sp.png"), overwrite = TRUE)
cat("Saved: panel_D_split_entropy_325sp\n")

# ── Combined C + D (single publication figure) ────────────────────────────────
p_combined <- p_c_raw / p_d_split +
  plot_annotation(
    title    = "Figure 1 C+D — CENP-A/CENH3 entropy analysis (325-sp)",
    subtitle = "C: per-group positional entropy.  D: CENP-A vs H3 split entropy.",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(FIG_DIR, "panel_CD_entropy_combined_325sp.pdf"), p_combined, width = 14, height = 9)
ggsave(file.path(FIG_DIR, "panel_CD_entropy_combined_325sp.png"), p_combined, width = 14, height = 9, dpi = 300)
file.copy(file.path(FIG_DIR, "panel_CD_entropy_combined_325sp.pdf"),
          file.path(TOP_FIG, "panel_CD_entropy_combined_325sp.pdf"), overwrite = TRUE)
file.copy(file.path(FIG_DIR, "panel_CD_entropy_combined_325sp.png"),
          file.path(TOP_FIG, "panel_CD_entropy_combined_325sp.png"), overwrite = TRUE)
cat("Saved: panel_CD_entropy_combined_325sp\n")

cat("\nAll entropy figures in:", FIG_DIR, "\n")
