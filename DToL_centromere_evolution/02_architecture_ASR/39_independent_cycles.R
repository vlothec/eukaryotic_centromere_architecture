#!/usr/bin/env Rscript
# 39_independent_cycles.R
#
# Count INDEPENDENT evolutionary cycles (not just tip counts).
# Strategy: for each reversal tip, find the deepest ancestral node where
# the mid-state (Trans for STS; Sat for TST) first appears on the root→tip
# path. Multiple tips sharing the same "entry node" = same cycle event.
# N unique entry nodes = N independent cycles.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(readxl); library(stringr); library(patchwork)
})
set.seed(42)

root      <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
cache_dir <- file.path(root, "outputs", "cyclical_ARD_irrevH", "cache")
out_dir   <- file.path(root, "outputs", "reversals")

state_levels <- c("H","Mixed","Sat","Trans")
state_colors <- c(H="#f46d43", Mixed="#8e44ad", Sat="#ff2d87", Trans="#4f84e8")
datasets       <- c("metazoa","viridiplantae","full","metazoa_treepl","viridiplantae_treepl","full_treepl")
dataset_labels <- c(full="Full tree", metazoa="Metazoa", viridiplantae="Viridiplantae", full_treepl="Full tree (treePL)", metazoa_treepl="Metazoa (treePL)", viridiplantae_treepl="Viridiplantae (treePL)")

# ── Annotations ──────────────────────────────────────────────────────────────
anno <- read.csv(file.path(root, "inputs", "full", "branch_symbol_anno.tsv"),
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
anno$tree_label <- gsub(" ", "", anno$Species)
map_arch <- c("Holocentric"="H","Satellite/transposon"="Mixed",
              "Satellite"="Sat","Transposon"="Trans")
label_fix_map <- c("daTanVulg1.hap2.1.fasta"="daTanVulg1.hap1.1.fa",
                   "rosCan_S27_v1.fasta.F2B.ChrOnly.fa"="rosCan_S27_v1")
meta <- suppressMessages(
  read_excel(file.path(dirname(root),"DTOL_327_master_March.xlsx"), sheet="Sheet1")
) %>%
  select(fasta, genus, species, classification, taxa1) %>%
  mutate(fasta_fixed = ifelse(fasta %in% names(label_fix_map),
                              unname(label_fix_map[fasta]), fasta),
         fasta_base  = str_replace(fasta_fixed, "\\.(fa|fasta)$", ""),
         sp_name     = paste(genus, species)) %>%
  distinct(fasta_base, .keep_all=TRUE)
sp_name_map <- setNames(meta$sp_name, meta$fasta_base)
taxa1_map   <- setNames(meta$taxa1,   meta$fasta_base)

# ── Helpers ───────────────────────────────────────────────────────────────────
load_dataset <- function(ds) {
  tf <- file.path(root, "inputs", ds, "tree_renamed.nw")
  if (!file.exists(tf)) return(NULL)
  tr <- read.tree(tf); tr <- midpoint.root(tr); tr <- multi2di(tr, random=FALSE)
  ch <- map_arch[anno$architecture[match(tr$tip.label, anno$tree_label)]]
  names(ch) <- tr$tip.label
  ch <- factor(ch, levels=state_levels); names(ch) <- tr$tip.label
  kp <- !is.na(ch); tr <- drop.tip(tr, tr$tip.label[!kp])
  ch <- factor(as.character(ch[kp]), levels=state_levels); names(ch) <- tr$tip.label
  list(tree=tr, char=ch, n=length(ch))
}

build_parent <- function(tree) {
  n <- length(tree$tip.label) + tree$Nnode
  p <- setNames(rep(NA_integer_, n), seq_len(n))
  for (i in seq_len(nrow(tree$edge))) p[tree$edge[i,2]] <- tree$edge[i,1]
  p
}

# ── Core function: find entry nodes for reversal tips ─────────────────────────
find_entry_nodes <- function(tree, all_states, tip_state, mid_state) {
  # For each tip with tip_state, check if root→tip path contains
  # ...tip_state...mid_state...tip_state (a reversal).
  # If yes, record the DEEPEST (most-ancestral) node where mid_state first
  # appears on that path. That is the "entry node" for the cycle event.
  # Tips sharing the same entry node = same independent event.

  n_tips  <- length(tree$tip.label)
  parent  <- build_parent(tree)

  entry_node_map <- list()   # tip_label -> entry_node_index

  for (ti in seq_len(n_tips)) {
    if (as.character(all_states[ti]) != tip_state) next

    # Walk root → tip
    path <- integer(0); curr <- ti
    while (!is.na(curr)) { path <- c(curr, path); curr <- parent[curr] }
    # path: root ... tip (left to right)
    st <- as.character(all_states[path])

    # Compress to runs
    rle_res <- rle(st); rv <- rle_res$values
    if (length(rv) < 3) next

    # Check for pattern: tip_state → mid_state → tip_state
    found <- any(rv[seq_len(length(rv)-2)] == tip_state &
                 rv[seq_len(length(rv)-2)+1] == mid_state &
                 rv[seq_len(length(rv)-2)+2] == tip_state)
    if (!found) next

    # Find the FIRST index in the full path where mid_state appears
    first_mid_pos <- which(st == mid_state)[1]
    if (is.na(first_mid_pos)) next

    entry_node <- path[first_mid_pos]
    entry_node_map[[tree$tip.label[ti]]] <- entry_node
  }
  entry_node_map
}

# ─────────────────────────────────────────────────────────────────────────────
# Main analysis
# ─────────────────────────────────────────────────────────────────────────────
all_results <- list()

for (ds in datasets) {
  cat("\n====", dataset_labels[ds], "====\n")
  dat <- load_dataset(ds); if(is.null(dat)) next
  cf  <- file.path(cache_dir, paste0("simmap_ARD_irrevH_", ds, ".rds"))
  if (!file.exists(cf)) { cat("  No simmap cache\n"); next }
  sm_list <- readRDS(cf)
  sm_sum  <- summary(sm_list)
  ace_mat <- sm_sum$ace
  n_tips  <- length(dat$tree$tip.label)

  # Infer ace dimensions — handle both tip+node and node-only formats
  if (nrow(ace_mat) == n_tips + dat$tree$Nnode) {
    ace_internal <- ace_mat[(n_tips+1):nrow(ace_mat), , drop=FALSE]
  } else {
    ace_internal <- ace_mat
  }

  # MAP states: tips = observed, internal = MAP from ace
  tip_states  <- as.character(dat$char)
  node_states <- apply(ace_internal, 1, function(r) colnames(ace_internal)[which.max(r)])
  node_probs  <- apply(ace_internal, 1, max)
  all_states  <- c(tip_states, node_states)
  all_states  <- factor(all_states, levels=state_levels)

  # ── STS reversals (Sat→Trans→Sat) ─────────────────────────────────────────
  sts_map <- find_entry_nodes(dat$tree, all_states, "Sat", "Trans")
  sts_n_tips   <- length(sts_map)
  sts_n_indep  <- length(unique(unlist(sts_map)))

  # ── TST reversals (Trans→Sat→Trans) ───────────────────────────────────────
  tst_map <- find_entry_nodes(dat$tree, all_states, "Trans", "Sat")
  tst_n_tips  <- length(tst_map)
  tst_n_indep <- length(unique(unlist(tst_map)))

  # ── Per-clade breakdown: which taxa are in each independent cycle? ──────────
  # Group STS reversal tips by entry node
  if (sts_n_indep > 0) {
    sts_groups <- tapply(names(sts_map), unlist(sts_map), function(tips) {
      taxa <- coalesce(taxa1_map[tips], tips)
      sp   <- coalesce(sp_name_map[tips], tips)
      paste0(unique(taxa), " (", length(tips), " tips)")
    })
    cat(sprintf("  STS: %d tips → %d independent cycles\n", sts_n_tips, sts_n_indep))
    for (i in seq_along(sts_groups))
      cat(sprintf("    Cycle %d: %s\n", i, sts_groups[[i]]))
  } else {
    cat("  STS: 0 reversals\n")
  }

  if (tst_n_indep > 0) {
    tst_groups <- tapply(names(tst_map), unlist(tst_map), function(tips) {
      taxa <- coalesce(taxa1_map[tips], tips)
      sp   <- coalesce(sp_name_map[tips], tips)
      paste0(unique(taxa), " (", length(tips), " tips)")
    })
    cat(sprintf("  TST: %d tips → %d independent cycles\n", tst_n_tips, tst_n_indep))
    for (i in seq_along(tst_groups))
      cat(sprintf("    Cycle %d: %s\n", i, tst_groups[[i]]))
  } else {
    cat("  TST: 0 reversals\n")
  }

  all_results[[ds]] <- list(
    dataset   = dataset_labels[ds],
    n_tips    = n_tips,
    sts_tips  = sts_n_tips,  sts_indep = sts_n_indep,
    tst_tips  = tst_n_tips,  tst_indep = tst_n_indep,
    sts_groups = if(sts_n_indep>0) sts_groups else NULL,
    tst_groups = if(tst_n_indep>0) tst_groups else NULL
  )
}

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n\n========= INDEPENDENT CYCLE SUMMARY =========\n")
summ <- bind_rows(lapply(all_results, function(r) {
  data.frame(
    dataset       = r$dataset,
    n_tips        = r$n_tips,
    STS_tips      = r$sts_tips,
    STS_indep     = r$sts_indep,
    TST_tips      = r$tst_tips,
    TST_indep     = r$tst_indep,
    total_indep   = r$sts_indep + r$tst_indep,
    stringsAsFactors = FALSE
  )
}))
print(as.data.frame(summ))
write.table(summ, file.path(out_dir,"independent_cycles_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ── Figure ────────────────────────────────────────────────────────────────────
plot_df <- summ %>%
  filter(dataset != "Full tree") %>%
  pivot_longer(c(STS_tips, STS_indep, TST_tips, TST_indep),
               names_to="metric", values_to="value") %>%
  mutate(
    rev_type = ifelse(grepl("STS", metric), "Sat→Trans→Sat", "Trans→Sat→Trans"),
    count_type = ifelse(grepl("tips", metric), "Tip count\n(overcounts shared ancestry)",
                        "Independent cycles\n(unique entry nodes)")
  ) %>%
  filter(value > 0)

p <- ggplot(plot_df, aes(x=count_type, y=value, fill=rev_type)) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(position="dodge", width=0.6, alpha=0.9) +
  geom_text(aes(label=value), position=position_dodge(0.6), vjust=-0.4, size=4) +
  scale_fill_manual(values=c("Sat→Trans→Sat"="#d73027","Trans→Sat→Trans"="#313695"),
                    name="Reversal type") +
  labs(title="Independent evolutionary cycles vs raw tip counts",
       subtitle="An 'independent cycle' = unique entry node into the mid-state\n(tips sharing an entry node = one cycle event)",
       x=NULL, y="Count") +
  theme_bw(base_size=13) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank(),
        legend.position="bottom", axis.text.x=element_text(size=10))

ggsave(file.path(out_dir,"independent_cycles.pdf"), p, width=10, height=6)
ggsave(file.path(out_dir,"independent_cycles.png"), p, width=10, height=6, dpi=300)

cat("\nOutputs in:", out_dir, "\n")
