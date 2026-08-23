#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/04_cenpa_phylogeny"
SP_DIR  <- file.path(BASE, "split_entropy")
FIG_DIR <- file.path(BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

C_UW    <- "#cfd8dc"
C_WT    <- "#1565c0"
C_SIG   <- "#c62828"
C_H3    <- "#90caf9"
C_CENPA <- "#f48fb1"

# ── Load data ─────────────────────────────────────────────────────────────────
uw_raw <- readLines(file.path(SP_DIR, "groupsim/groupsim_gap085.txt"))
uw_raw <- uw_raw[!grepl("^#", uw_raw) & nchar(trimws(uw_raw)) > 0]
uw <- do.call(rbind, lapply(uw_raw, function(l) {
  p <- strsplit(l, "\t")[[1]]
  if (length(p) < 2) return(NULL)
  data.frame(pos = as.integer(p[1]) + 1L,
             score_uw = suppressWarnings(as.numeric(p[2])))
})) %>% filter(!is.na(pos))

wt <- read.delim(file.path(SP_DIR,
                            "groupsim_weighted/groupsim_cenpa_h3_clade_gap085.tsv"),
                 stringsAsFactors = FALSE)
names(wt)[names(wt)=="groupsim_clade"] <- "score_wt"
names(wt)[names(wt)=="z_clade"]        <- "z_score"

df <- merge(uw, wt[,c("pos","score_wt","z_score","sig_clade")],
            by = "pos", all = TRUE)
df <- df[order(df$pos),]
cat("Positions:", nrow(df), "| sig:", sum(df$sig_clade==1, na.rm=TRUE), "\n")

helix <- read.delim(file.path(SP_DIR, "helix_positions_gap085.tsv"),
                    stringsAsFactors = FALSE)

annot <- data.frame(
  pos   = c(79, 96, 101, 114, 119, 131),
  label = c("Gln/diverse","Phe/Trp","Val/Leu",
            "Gly/diverse","Thr/Ala","Ile/Leu"),
  stringsAsFactors = FALSE
)
annot <- merge(annot, df[,c("pos","score_wt")], by="pos", all.x=TRUE)

# same limits for ALL panels so patchwork aligns correctly
x_min <- min(df$pos, na.rm=TRUE) - 0.5
x_max <- max(df$pos, na.rm=TRUE) + 0.5

# ── Panel 1: helix track (CENP-A only) ───────────────────────────────────────
helix_cenpa <- helix[helix$histone == "CENPA", ]

p_helix <- ggplot(helix_cenpa) +
  geom_rect(aes(xmin=start-0.5, xmax=end+0.5, ymin=0.1, ymax=0.9),
            fill=C_CENPA, alpha=0.85, colour=NA) +
  annotate("text", x=x_min+1, y=0.5, label="CENP-A helices",
           hjust=0, size=3.2, colour="#880e4f") +
  scale_x_continuous(limits=c(x_min,x_max), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_void(base_size=9) +
  theme(plot.margin=margin(2,2,0,2))

# ── Panel 2: GroupSim bars ────────────────────────────────────────────────────
df_uw <- df[!is.na(df$score_uw),]
df_wt <- df[!is.na(df$score_wt),]
df_wt$fill_col <- ifelse(df_wt$pos %in% annot$pos, "sig", "wt")
z2_thresh <- min(df$score_wt[df$sig_clade==1], na.rm=TRUE)

df_uw$fill_col <- "uw"

p_main <- ggplot() +
  geom_col(data=df_uw, aes(x=pos-0.22, y=score_uw, fill=fill_col),
           width=0.40, alpha=1.0) +
  geom_col(data=df_wt, aes(x=pos+0.22, y=score_wt, fill=fill_col),
           width=0.40, alpha=0.92) +
  geom_segment(data=annot,
               aes(x=pos+0.22, xend=pos+0.22,
                   y=pmin(score_wt,1.0), yend=1.12),
               colour=C_SIG, linewidth=0.5) +
  geom_label(data=annot,
             aes(x=pos+0.22, y=1.14 + rep(c(0,0.09),3), label=label),
             size=2.5, colour=C_SIG, fill="white",
             label.size=0.35, label.padding=unit(0.12,"cm"),
             hjust=0.5, vjust=0) +
  scale_fill_manual(
    values = c(uw="grey",  sig=C_SIG, wt=C_WT),
    labels = c(uw="Unweighted", sig="Significant (z >= 2)", wt="Clade-weighted"),
    name   = NULL,
    breaks = c("uw","wt","sig")          # legend order
  ) +
  scale_x_continuous(limits=c(x_min,x_max), expand=c(0,0)) +
  scale_y_continuous(name="GroupSim score", limits=c(0,1.40),
                     breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(x=NULL) +
  theme_bw(base_size=11) +
  theme(axis.text.x         = element_blank(),
        axis.ticks.x        = element_blank(),
        panel.grid.minor    = element_blank(),
        panel.grid.major.x  = element_blank(),
        legend.position     = "right",
        legend.text         = element_text(size=10),
        legend.key.size     = unit(0.5,"cm"),
        legend.spacing.y    = unit(0.3,"cm"),
        plot.margin         = margin(0,2,0,2))

# ── Panel 3: z-score strip ────────────────────────────────────────────────────
df_z <- df[!is.na(df$z_score),]

p_zscore <- ggplot(df_z, aes(x=pos, y=0.5, fill=z_score)) +
  geom_tile(width=1, height=1) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#c62828",
                       midpoint=0, limits=c(-4,4), name="z-score",
                       guide=guide_colorbar(barwidth=1.0, barheight=5,
                                            title.position="top",
                                            title.theme=element_text(size=10),
                                            label.theme=element_text(size=9))) +
  scale_x_continuous(limits=c(x_min, x_max), expand=c(0,0),
                     name="Alignment position (trimmed, gap <= 85%)") +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_void(base_size=9) +
  theme(axis.title.x=element_text(size=10, margin=margin(t=4)),
        axis.text.x=element_text(size=9),
        axis.ticks.x=element_line(linewidth=0.3),
        legend.position="right",
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        panel.background=element_rect(fill="#cccccc", colour=NA),
        plot.margin=margin(0,2,2,2))

# ── Stack panels with patchwork ───────────────────────────────────────────────
p_out <- p_helix / p_main / p_zscore +
  plot_layout(guides="collect", heights=c(0.7, 6, 0.5)) +
  plot_annotation(
    title    = "CENP-A/CENH3 (n=422) vs H3-like (n=897)  |  clade-weighted GroupSim",
    subtitle = "Grey = unweighted  |  Blue = clade-weighted  |  Red = z >= 2  |  H3>CENPA label format",
    theme    = theme(plot.title    = element_text(face="bold", size=11),
                     plot.subtitle = element_text(size=8, colour="grey40"),
                     plot.background = element_rect(fill="white", colour=NA))
  )

ggsave(file.path(FIG_DIR,"groupsim_cenpa_vs_h3_gap085_pub.pdf"),
       p_out, width=13, height=6.5, bg="white")
ggsave(file.path(FIG_DIR,"groupsim_cenpa_vs_h3_gap085_pub.png"),
       p_out, width=13, height=6.5, dpi=300, bg="white")
cat("Saved: groupsim_cenpa_vs_h3_gap085_pub.{pdf,png}\n")
