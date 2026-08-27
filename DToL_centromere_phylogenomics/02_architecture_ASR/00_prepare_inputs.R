#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(readxl)
  library(dplyr)
  library(stringr)
})

# Run this script from ASR_March2026_327species/
root <- normalizePath(".")

# Source files
excel_path <- file.path(dirname(root), "DTOL_327_master_March.xlsx")
itol_path <- "/home/jg2070/Desktop/iTOL/SPECIESTREE_May2025_193species/branch_symbol_anno.tsv"

# Trees (use calibrated chronograms if present)
trees <- list(
  full = file.path(dirname(root), "../full_327species_calibrated.nex"),
  metazoa = file.path(dirname(root), "../metazoa_calibrated.nex"),
  viridiplantae = file.path(dirname(root), "../viridiplantae_calibrated.nex")
)

fallback_trees <- list(
  full = file.path(dirname(root), "../grouped_buscos_phylooutput/04_concatenated/327species_255buscos_run_notpartitioned.contree"),
  metazoa = file.path(dirname(root), "../grouped_buscos_Metazoa_phylooutput/04_concatenated/concatenated_alignment_Metazoa_runnotpartitioned.contree"),
  viridiplantae = file.path(dirname(root), "../grouped_buscos_Viridiplantae_phylooutput/04_concatenated/Viridiplantae_viridiplantaebusco_run_notpartitioned.contree")
)

out_root <- file.path(root, "inputs")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

strip_ext <- function(x) str_replace(x, "\\.(fa|fasta)$", "")

# Known label fixes
label_fix_map <- c(
  "daTanVulg1.hap2.1.fasta" = "daTanVulg1.hap1.1.fa",
  "rosCan_S27_v1.fasta.F2B.ChrOnly.fa" = "rosCan_S27_v1"
)

# Classification -> architecture mapping
class_to_arch <- tibble::tibble(
  classification = c(
    "Unknown",
    "Transposon",
    "Satellite",
    "Holocentric",
    "Satellite/transposon",
    "Monocentric sequence unknown"
  ),
  architecture = c(
    "Unknown",
    "Transposon",
    "Satellite",
    "Holocentric",
    "Satellite/transposon",
    "Monocentric sequence unknown"
  )
)

# iTOL style mapping
itol <- read.delim(itol_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
arch_map_raw <- itol %>%
  select(architecture, symbol, size, colour, border, width) %>%
  distinct()

# If an architecture is missing in iTOL, provide sensible defaults
arch_defaults <- tibble::tibble(
  architecture = c("Holocentric", "Monocentric sequence unknown", "Satellite/transposon"),
  symbol = c("triangle left", "rectangle", "circle"),
  size = c(2, 2, 2),
  colour = c("#fb5609", "#808080", "#8338ec"),
  border = c("yes", "yes", "yes"),
  width = c(1, 1, 1)
)

arch_map <- bind_rows(arch_map_raw, arch_defaults) %>%
  distinct(architecture, .keep_all = TRUE)

# Metadata
meta <- read_excel(excel_path, sheet = "Sheet1") %>%
  select(fasta, genus, species, classification) %>%
  mutate(fasta_fixed = ifelse(fasta %in% names(label_fix_map),
                              unname(label_fix_map[fasta]),
                              fasta),
         fasta_base = strip_ext(fasta_fixed)) %>%
  left_join(class_to_arch, by = "classification") %>%
  left_join(arch_map, by = "architecture")

# Build inputs for each tree
report <- list()

read_any_tree <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (grepl("\\.nex$", path, ignore.case = TRUE)) {
    cat("Reading NEXUS:", path, "\n")
    return(read.nexus(path))
  }
  cat("Reading tree:", path, "\n")
  read.tree(path)
}

for (nm in names(trees)) {
  tree_path <- trees[[nm]]
  tree <- read_any_tree(tree_path)
  if (is.null(tree)) {
    tree_path <- fallback_trees[[nm]]
    tree <- read_any_tree(tree_path)
  }

  # Fix and normalize tip labels
  tip_fixed <- ifelse(tree$tip.label %in% names(label_fix_map),
                      unname(label_fix_map[tree$tip.label]),
                      tree$tip.label)
  tip_base <- strip_ext(tip_fixed)

  tree$tip.label <- tip_base

  # Subset metadata to this tree
  meta_sub <- meta %>%
    filter(fasta_base %in% tip_base) %>%
    distinct(fasta_base, .keep_all = TRUE)

  # Any missing tips
  missing <- setdiff(tip_base, meta_sub$fasta_base)

  # Create anno file (iTOL-like columns)
  anno <- meta_sub %>%
    mutate(
      Species = fasta_base,
      Species.1 = fasta_base,
      symbol = ifelse(is.na(symbol), "rectangle", symbol),
      size = ifelse(is.na(size), 2, size),
      colour = ifelse(is.na(colour), "#808080", colour),
      border = ifelse(is.na(border), "yes", border),
      width = ifelse(is.na(width), 1, width),
      architecture = ifelse(is.na(architecture), "Unknown", architecture)
    ) %>%
    select(Species, Species.1, symbol, size, colour, border, width, architecture)

  out_dir <- file.path(out_root, nm)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  write.tree(tree, file = file.path(out_dir, "tree_renamed.nw"))
  write.table(anno, file = file.path(out_dir, "branch_symbol_anno.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  report[[nm]] <- list(
    tips = length(tip_base),
    matched = nrow(meta_sub),
    missing = missing
  )
}

# Write report
report_path <- file.path(out_root, "prepare_inputs_report.txt")
con <- file(report_path, "wt")
for (nm in names(report)) {
  cat("Dataset:", nm, "\n", file = con)
  cat("Tips:", report[[nm]]$tips, "\n", file = con)
  cat("Matched:", report[[nm]]$matched, "\n", file = con)
  cat("Missing tips:", length(report[[nm]]$missing), "\n", file = con)
  if (length(report[[nm]]$missing) > 0) {
    cat("Missing labels:\n", file = con)
    cat(paste(report[[nm]]$missing, collapse = ", "), "\n", file = con)
  }
  cat("\n", file = con)
}
close(con)

cat("Inputs written to:", out_root, "\n")
cat("Report:", report_path, "\n")
