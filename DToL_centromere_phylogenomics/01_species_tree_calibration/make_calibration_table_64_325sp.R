#!/usr/bin/env Rscript
# make_calibration_table_64_325sp.R
# Comprehensive supplementary calibration table for ALL 64 constraints actually
# used to date the 325-sp chronogram (supersedes the 51-row fossil-only table).
# Sources: over_calib.tsv (constraint bounds) + calibration_nodes_timetree.tsv
#          (real TimeTree median/range from the API + provenance) + calib64 (labels).
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(gridExtra); library(grid)})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR <- file.path(BASE, "DToL_phylogenomics_publication_325genomes/01_species_tree")
QC_DIR  <- file.path(PUB_DIR, "outputs/calibration_qc")
OUT_DIR <- file.path(PUB_DIR, "outputs")
FIG_DIR <- file.path(PUB_DIR, "figures")

oc  <- read.delim(file.path(PUB_DIR, "over_calib.tsv"), stringsAsFactors = FALSE)
ntt <- read.delim(file.path(QC_DIR, "calibration_nodes_timetree.tsv"), stringsAsFactors = FALSE)
c64 <- read.delim(file.path(QC_DIR, "calib64_constraints_for_panelA.tsv"), stringsAsFactors = FALSE)

clean_label <- function(x) {
  x <- gsub("Metazoa_node16","Avemetatarsalia/Ornithodira node",x)
  x <- gsub("Metazoa_node17","Aculeata/Tipuloidea node",x)
  x <- gsub("Vir_node10","Alismatales/Poales node",x)
  x <- gsub("Vir_node11","Poaceae/Asteraceae node",x)
  x <- gsub("Vir_node12","Fabales/Caryophyllales node",x)
  x <- gsub("Vir_node16","Asteraceae node",x)
  x <- gsub("Vir_node17","Moss node",x)
  x <- gsub("Vir_node1","Potamogetonaceae/Droseraceae",x)
  x <- gsub("Marchantiophyta","Liverwort node",x)
  x <- gsub("Vir_node5","Cyperaceae clade",x)
  x <- gsub("Vir_node8","Caryophyllales/Poales node",x)
  gsub("_"," ",x)
}
num <- function(x){x<-as.character(x); suppressWarnings(as.numeric(ifelse(x %in% c("None","NA",""),NA,x)))}

clade_order <- c("Root","Opisthokonta","Metazoa","Viridiplantae","Fungi")
d <- ntt %>%
  mutate(age_min = num(oc$age_min[match(label, oc$label)]),
         age_max = num(oc$age_max[match(label, oc$label)]),
         disp    = c64$disp_label[match(label, c64$label)],
         clade_f = factor(clade, levels = clade_order),
         tt_median = num(tt_median), tt_adjusted = num(tt_adjusted),
         tt_ci_low = num(tt_ci_low), tt_ci_high = num(tt_ci_high),
         our_age = num(our_age)) %>%
  arrange(clade_f, age_min)

src_short <- function(s) dplyr::case_when(
  grepl("^TimeTree", s)     ~ "TimeTree",
  TRUE ~ s)

supp <- transmute(d,
  Clade            = as.character(clade_f),
  Node             = clean_label(disp),
  Taxon_A          = taxonA,
  Taxon_B          = taxonB,
  Constraint_min_Mya = round(age_min, 2),
  Constraint_max_Mya = round(age_max, 2),
  Calibrated_age_Mya = round(our_age, 2),
  TimeTree_median_Mya = ifelse(is.na(tt_median), NA, round(tt_median, 2)),
  TimeTree_adjusted_Mya = ifelse(is.na(tt_adjusted), NA, round(tt_adjusted, 2)),
  TimeTree_range_low  = ifelse(is.na(tt_ci_low),  NA, round(tt_ci_low, 2)),
  TimeTree_range_high = ifelse(is.na(tt_ci_high), NA, round(tt_ci_high, 2)),
  Source           = src_short(source))

out_tsv <- file.path(OUT_DIR, "calibration_table_325sp_supp_62.tsv")
write.table(supp, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
cat("Saved:", out_tsv, "  (", nrow(supp), "nodes )\n")
print(table(supp$Clade)); cat("TimeTree:", sum(supp$Source=="TimeTree"),
    " | Literature:", sum(grepl("^Lit", supp$Source)), "\n")

# ── rendered table figure (PDF + PNG) ─────────────────────────────────────────
disp <- supp
for (cc in c("Constraint_min_Mya","Constraint_max_Mya","Calibrated_age_Mya",
             "TimeTree_median_Mya","TimeTree_adjusted_Mya","TimeTree_range_low","TimeTree_range_high"))
  disp[[cc]] <- ifelse(is.na(disp[[cc]]), "—", sprintf("%.1f", disp[[cc]]))
disp$TimeTree_range <- ifelse(disp$TimeTree_range_low=="—", "—",
                              paste0(disp$TimeTree_range_low, "–", disp$TimeTree_range_high))
disp2 <- disp[, c("Clade","Node","Taxon_A","Taxon_B","Constraint_min_Mya","Constraint_max_Mya",
                  "Calibrated_age_Mya","TimeTree_median_Mya","TimeTree_adjusted_Mya","TimeTree_range","Source")]
names(disp2) <- c("Clade","Node","Taxon A","Taxon B","Min\n(Ma)","Max\n(Ma)",
                  "Calib.\nage","TT\nmedian","TT\nadjusted","TT range","Source")

clade_bg <- c(Root="#ECEFF1",Opisthokonta="#CFD8DC",Metazoa="#E3F2FD",
              Viridiplantae="#E8F5E9",Fungi="#F3E5F5")
th <- ttheme_minimal(
  base_size = 6.5,
  core = list(fg_params = list(hjust = 0, x = 0.02,
               fontface = ifelse(grepl("^Lit", disp2$Source), 3, 1)),
              bg_params = list(fill = clade_bg[disp2$Clade], col = "grey85", lwd = 0.4)),
  colhead = list(fg_params = list(fontface = 2, fontsize = 7),
                 bg_params = list(fill = "#37474F", col = NA),
                 fg_params2 = NULL))
th$colhead$fg_params$col <- "white"
g <- tableGrob(disp2, rows = NULL, theme = th)
title <- textGrob("Supplementary Table - Calibration constraints used for the 325-species DToL chronogram (all 62 nodes)",
                  gp = gpar(fontface = "bold", fontsize = 10), x = 0.01, hjust = 0)
sub <- textGrob(paste0("Min/Max = calibration constraint bounds imposed on chronos.  Calib. age = calibrated node age (chronos correlated).  ",
                       "TT median / TT range = TimeTree pairwise median and range (API).  All 62 constraints TimeTree-derived."),
                gp = gpar(fontsize = 6.5, col = "grey35"), x = 0.01, hjust = 0)
tab <- arrangeGrob(g, top = sub)
full <- arrangeGrob(title, tab, heights = unit.c(unit(1.4,"lines"), unit(1,"null")))

for (ext in c("pdf","png")) {
  out <- file.path(FIG_DIR, paste0("calibration_table_325sp_supp_62.", ext))
  ggsave(out, full, width = 11, height = 16, dpi = 300,
         bg = "white", limitsize = FALSE)
  cat("Saved:", out, "\n")
}
