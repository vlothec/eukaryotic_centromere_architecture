#!/usr/bin/env Rscript
# 43_cycles_tree_mk_parsimony.R
#
# Plot cycles on tree using TWO ancestral reconstruction methods:
#   1. Mk marginal ASR  -- ace() with ARD_irrevH Q (proper likelihood)
#   2. Maximum parsimony -- ancestral.pars() (Fitch parsimony, model-free)
#
# Cycle entry nodes are identified under each method independently.
# Tips coloured by observed state; internal nodes by MAP state from each method.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(ggtree); library(readxl); library(stringr)
  library(patchwork)
})
set.seed(42)

root    <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
out_dir <- file.path(root, "outputs", "cycles_mk_parsimony")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

state_levels <- c("H","Mixed","Sat","Trans")
state_colors <- c(H="#f46d43", Mixed="#8e44ad", Sat="#ff2d87", Trans="#4f84e8")
state_labels <- c(H="Holocentric", Mixed="Sat/Transposon", Sat="Satellite", Trans="Transposon")
datasets       <- c("metazoa_chronos_correlated","viridiplantae_chronos_correlated")
dataset_labels <- c(metazoa_chronos_correlated="Metazoa (chronos correlated)", viridiplantae_chronos_correlated="Viridiplantae (chronos correlated)")

# ── Annotations ──────────────────────────────────────────────────────────────
anno <- read.csv(file.path(root,"inputs","full","branch_symbol_anno.tsv"),
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
anno$tree_label <- gsub(" ","",anno$Species)
map_arch <- c("Holocentric"="H","Satellite/transposon"="Mixed",
              "Satellite"="Sat","Transposon"="Trans")
label_fix_map <- c("daTanVulg1.hap2.1.fasta"="daTanVulg1.hap1.1.fa",
                   "rosCan_S27_v1.fasta.F2B.ChrOnly.fa"="rosCan_S27_v1")
meta <- suppressMessages(
  read_excel(file.path(dirname(root),"DTOL_327_master_March.xlsx"), sheet="Sheet1")
) %>%
  select(fasta, genus, species, taxa1) %>%
  mutate(fasta_fixed = ifelse(fasta %in% names(label_fix_map),
                              unname(label_fix_map[fasta]), fasta),
         fasta_base  = str_replace(fasta_fixed,"\\.(fa|fasta)$",""),
         sp_name     = paste(genus, species)) %>%
  distinct(fasta_base, .keep_all=TRUE)
sp_name_map <- setNames(meta$sp_name, meta$fasta_base)
taxa1_map   <- setNames(meta$taxa1,   meta$fasta_base)

# ── Helpers ───────────────────────────────────────────────────────────────────
load_dataset <- function(ds) {
  tr <- read.tree(file.path(root,"inputs",ds,"tree_renamed.nw"))
  tr <- midpoint.root(tr); tr <- multi2di(tr, random=FALSE)
  ch <- map_arch[anno$architecture[match(tr$tip.label,anno$tree_label)]]
  names(ch) <- tr$tip.label
  ch <- factor(ch, levels=state_levels); kp <- !is.na(ch)
  tr <- drop.tip(tr, tr$tip.label[!kp])
  ch <- factor(as.character(ch[kp]), levels=state_levels); names(ch) <- tr$tip.label
  list(tree=tr, char=ch, n=length(ch))
}

build_parent <- function(tree) {
  n <- length(tree$tip.label) + tree$Nnode
  p <- setNames(rep(NA_integer_,n), seq_len(n))
  for (i in seq_len(nrow(tree$edge))) p[tree$edge[i,2]] <- tree$edge[i,1]
  p
}

find_entry_nodes <- function(tree, all_states, tip_state, mid_state) {
  n_tips <- length(tree$tip.label); parent <- build_parent(tree)
  entry <- list()
  for (ti in seq_len(n_tips)) {
    if (as.character(all_states[ti]) != tip_state) next
    path <- integer(0); curr <- ti
    while (!is.na(curr)) { path <- c(curr,path); curr <- parent[curr] }
    st <- as.character(all_states[path]); rv <- rle(st)$values
    if (length(rv)<3) next
    found <- any(rv[seq_len(length(rv)-2)]==tip_state &
                 rv[seq_len(length(rv)-2)+1]==mid_state &
                 rv[seq_len(length(rv)-2)+2]==tip_state)
    if (!found) next
    fm <- which(st==mid_state)[1]
    if (!is.na(fm) && path[fm] != path[1]) entry[[tree$tip.label[ti]]] <- path[fm]
  }
  entry
}

# ── Design matrix: ARD_irrevH ─────────────────────────────────────────────────
dm_cyc <- matrix(0L,4,4,dimnames=list(state_levels,state_levels))
k <- 1L
for (i in state_levels) for (j in state_levels)
  if (i!=j) { dm_cyc[i,j] <- k; k <- k+1L }
dm_cyc["H","Mixed"] <- 0L; dm_cyc["H","Sat"] <- 0L; dm_cyc["H","Trans"] <- 0L

# ── Make tree plot (shared helper) ────────────────────────────────────────────
make_cycle_tree <- function(tree, tip_ann, node_ann_xy, layout="rectangular",
                            title_str="", subtitle_str="") {

  tr_obj <- if (layout=="fan")
    ggtree(tree, layout="fan", open.angle=15, linewidth=0.25, color="grey50")
  else
    ggtree(tree, linewidth=0.25, color="grey50")

  p <- tr_obj %<+% tip_ann +
    # Tips
    geom_tippoint(aes(color=obs_state), size=2, alpha=0.9) +
    # Internal node MAP state (small, only confident ones, excluding root)
    geom_point(data = node_ann_xy %>% filter(map_prob > 0.75, !is_root),
               aes(x=x, y=y, color=map_state),
               size=1.2, alpha=0.5, inherit.aes=FALSE) +
    # Root node: large diamond coloured by MAP state + probability label
    geom_point(data = node_ann_xy %>% filter(is_root),
               aes(x=x, y=y, color=map_state),
               shape=18, size=9, inherit.aes=FALSE) +
    geom_label(data = node_ann_xy %>% filter(is_root),
               aes(x=x, y=y, label=root_label),
               size=2.4, vjust=2.2, fill="white", label.size=0.3,
               color="black", inherit.aes=FALSE) +
    # Cycle entry: numbered circles
    geom_point(data = node_ann_xy %>% filter(!is.na(cycle_num)),
               aes(x=x, y=y, fill=cycle_type),
               shape=21, size=7, color="white", stroke=1.5,
               inherit.aes=FALSE) +
    geom_text(data = node_ann_xy %>% filter(!is.na(cycle_num)),
              aes(x=x, y=y, label=cycle_label),
              size=3, fontface="bold", color="white", inherit.aes=FALSE) +
    scale_color_manual(values=state_colors, name="State",
                       labels=state_labels, breaks=names(state_labels)) +
    scale_fill_manual(values=c("Sat->Trans->Sat"="#d73027",
                                "Trans->Sat->Trans"="#313695"),
                      name="Cycle", na.value=NA) +
    labs(title=title_str, subtitle=subtitle_str) +
    theme_tree2() +
    theme(plot.title=element_text(face="bold",size=12),
          plot.subtitle=element_text(size=9,color="grey40"),
          legend.position="right")

  if (layout=="fan")
    p <- p + geom_tiplab(aes(label=sp_name), size=1.5, offset=0.5, hjust=0)
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN LOOP
# ─────────────────────────────────────────────────────────────────────────────
for (ds in datasets) {
  cat("\n====", dataset_labels[ds], "====\n")
  dat <- load_dataset(ds)
  n_tips <- length(dat$tree$tip.label)

  # Numeric character vector for ape (ace needs named numeric or factor)
  char_num <- as.integer(dat$char)
  names(char_num) <- dat$tree$tip.label

  # ── 1. Mk marginal ASR via phytools rerootingMethod ──────────────────────
  cat("  Fitting ARD_irrevH for Mk ASR...\n")
  fit_mk <- fitMk(dat$tree, dat$char, model=dm_cyc,
                  states=state_levels, control=list(maxit=3000))
  Q_mk <- matrix(0,4,4,dimnames=list(state_levels,state_levels))
  for (i in 1:4) for (j in 1:4)
    if (dm_cyc[i,j]>0) Q_mk[i,j] <- fit_mk$rates[dm_cyc[i,j]]
  diag(Q_mk) <- -rowSums(Q_mk)

  # phytools ancr: marginal ASR for non-symmetric Q matrices (ARD)
  # ancr() takes a fitted fitMk object and returns marginal node probabilities
  cat("  Running phytools ancr for marginal Mk ASR...\n")
  ancr_res <- ancr(fit_mk)
  # ancr_res$ace: internal nodes x states marginal probability matrix
  mk_probs <- ancr_res$ace
  if (!all(colnames(mk_probs) == state_levels))
    mk_probs <- mk_probs[, state_levels, drop=FALSE]
  mk_map     <- apply(mk_probs, 1, function(r) state_levels[which.max(r)])
  mk_maxprob <- apply(mk_probs, 1, max)

  cat("  Mk MAP states at internal nodes (top 10):\n")
  print(head(data.frame(node=n_tips+seq_len(dat$tree$Nnode),
                        state=mk_map, prob=round(mk_maxprob,2)), 10))

  # ── 2. Maximum parsimony (Fitch algorithm) ────────────────────────────────
  cat("  Running Fitch parsimony reconstruction...\n")

  fitch_asr <- function(tree, char) {
    n_tips  <- length(tree$tip.label)
    n_total <- n_tips + tree$Nnode
    parent_of <- setNames(rep(NA_integer_, n_total), seq_len(n_total))
    for (i in seq_len(nrow(tree$edge)))
      parent_of[tree$edge[i,2]] <- tree$edge[i,1]
    children <- vector("list", n_total)
    for (i in seq_len(nrow(tree$edge)))
      children[[tree$edge[i,1]]] <- c(children[[tree$edge[i,1]]], tree$edge[i,2])

    # Tip states
    state_sets <- vector("list", n_total)
    for (ti in seq_len(n_tips))
      state_sets[[ti]] <- as.character(char[ti])

    # Downpass (post-order: visit children before parents)
    post_order <- rev(unique(tree$edge[,1]))
    for (nd in post_order) {
      cs <- lapply(children[[nd]], function(ch) state_sets[[ch]])
      cs <- cs[!sapply(cs, is.null)]
      if (!length(cs)) next
      inter <- Reduce(intersect, cs)
      state_sets[[nd]] <- if (length(inter) > 0) inter else Reduce(union, cs)
    }

    # Uppass (pre-order: parents before children)
    root <- n_tips + 1L
    tip_freq <- sort(table(as.character(char)), decreasing=TRUE)
    root_states <- state_sets[[root]]
    # Pick most frequent among root's possible states
    best <- names(tip_freq)[names(tip_freq) %in% root_states][1]
    state_sets[[root]] <- if (!is.na(best)) best else root_states[1]

    pre_order <- c(root, tree$edge[,2][tree$edge[,2] > n_tips])
    for (nd in pre_order) {
      if (nd == root) next
      p <- parent_of[nd]
      if (is.na(p)) next
      inter <- intersect(state_sets[[p]], state_sets[[nd]])
      state_sets[[nd]] <- if (length(inter) > 0) inter[1] else state_sets[[nd]][1]
    }

    # Return internal node MAP states (single string per node)
    sapply((n_tips+1):n_total, function(i)
      if (!is.null(state_sets[[i]])) state_sets[[i]][1] else NA_character_)
  }

  pars_map     <- fitch_asr(dat$tree, dat$char)
  pars_maxprob <- rep(1.0, length(pars_map))   # parsimony gives a single MAP state

  cat("  Parsimony MAP states at internal nodes (top 10):\n")
  print(head(data.frame(node=n_tips+seq_len(dat$tree$Nnode),
                        state=pars_map), 10))

  # ── Build all_states for cycle detection ─────────────────────────────────
  tip_states_vec <- as.character(dat$char)
  mk_all    <- factor(c(tip_states_vec, mk_map),    levels=state_levels)
  pars_all  <- factor(c(tip_states_vec, pars_map),  levels=state_levels)

  # ── Find cycles under each method ────────────────────────────────────────
  for (method_name in c("Mk","Parsimony")) {
    all_states_m <- if (method_name=="Mk") mk_all else pars_all
    node_map_m   <- if (method_name=="Mk") mk_map  else pars_map
    node_prob_m  <- if (method_name=="Mk") mk_maxprob else pars_maxprob

    sts_map_m <- find_entry_nodes(dat$tree, all_states_m, "Sat",   "Trans")
    tst_map_m <- find_entry_nodes(dat$tree, all_states_m, "Trans", "Sat")
    sts_entry <- unique(unlist(sts_map_m))
    tst_entry <- unique(unlist(tst_map_m))

    cat(sprintf("\n  [%s] STS: %d tips, %d independent cycles\n",
                method_name, length(sts_map_m), length(sts_entry)))
    cat(sprintf("  [%s] TST: %d tips, %d independent cycles\n",
                method_name, length(tst_map_m), length(tst_entry)))

    # Cycle key
    for (i in seq_along(sts_entry)) {
      en <- sts_entry[i]
      tips <- names(sts_map_m)[unlist(sts_map_m)==en]
      cat(sprintf("    STS Cycle %d: %s [%d tips]\n", i,
                  paste(unique(coalesce(taxa1_map[tips],tips)),collapse=", "),
                  length(tips)))
    }
    for (i in seq_along(tst_entry)) {
      en <- tst_entry[i]
      tips <- names(tst_map_m)[unlist(tst_map_m)==en]
      cat(sprintf("    TST Cycle %d: %s [%d tips]\n", i+length(sts_entry),
                  paste(unique(coalesce(taxa1_map[tips],tips)),collapse=", "),
                  length(tips)))
    }

    # ── Build node_ann ──────────────────────────────────────────────────────
    sts_cn <- setNames(seq_along(sts_entry), as.character(sts_entry))
    tst_cn <- setNames(seq_along(tst_entry)+length(sts_entry),
                       as.character(tst_entry))

    root_node <- n_tips + 1L
    node_ann <- data.frame(
      node      = (n_tips+1):(n_tips+dat$tree$Nnode),
      map_state = node_map_m,
      map_prob  = node_prob_m,
      stringsAsFactors=FALSE
    ) %>% mutate(
      is_root    = node == root_node,
      root_label = ifelse(is_root,
                          sprintf("Root: %s\np=%.2f", map_state, map_prob),
                          NA_character_),
      is_sts = node %in% sts_entry,
      is_tst = node %in% tst_entry,
      cycle_num = case_when(
        is_sts ~ sts_cn[as.character(node)],
        is_tst ~ tst_cn[as.character(node)],
        TRUE   ~ NA_integer_),
      cycle_label = ifelse(!is.na(cycle_num), as.character(cycle_num), NA_character_),
      cycle_type  = case_when(
        is_sts ~ "Sat->Trans->Sat",
        is_tst ~ "Trans->Sat->Trans",
        TRUE   ~ NA_character_)
    )

    tip_ann <- data.frame(
      label    = dat$tree$tip.label,
      obs_state= as.character(dat$char),
      sp_name  = coalesce(sp_name_map[dat$tree$tip.label], dat$tree$tip.label),
      stringsAsFactors=FALSE)

    # Get layout coordinates
    tr_rect <- ggtree(dat$tree, linewidth=0.25)
    layout_rect <- tr_rect$data
    node_ann_rect <- node_ann %>%
      left_join(layout_rect %>% select(node, x, y), by="node")

    tr_fan <- ggtree(dat$tree, layout="fan", open.angle=15, linewidth=0.25)
    layout_fan <- tr_fan$data
    node_ann_fan <- node_ann %>%
      left_join(layout_fan %>% select(node, x, y), by="node")

    # ── Rectangular ──────────────────────────────────────────────────────────
    n_cycles_total <- length(sts_entry) + length(tst_entry)
    p_rect <- make_cycle_tree(
      dat$tree, tip_ann, node_ann_rect, layout="rectangular",
      title_str = paste0(dataset_labels[ds]," — ", method_name,
                         " ASR: ", n_cycles_total," independent cycles"),
      subtitle_str = paste0("Red circles = Sat->Trans->Sat (",length(sts_entry),
                            ") | Blue = Trans->Sat->Trans (",length(tst_entry),
                            ")\nSmall dots = internal node MAP state (prob>0.75)")
    )
    h_rect <- max(8, n_tips * 0.065 + 2)
    ggsave(file.path(out_dir, paste0(ds,"_",tolower(method_name),"_rectangular.pdf")),
           p_rect, width=14, height=h_rect)
    ggsave(file.path(out_dir, paste0(ds,"_",tolower(method_name),"_rectangular.png")),
           p_rect, width=14, height=h_rect, dpi=200)

    # ── Fan ───────────────────────────────────────────────────────────────────
    p_fan <- make_cycle_tree(
      dat$tree, tip_ann, node_ann_fan, layout="fan",
      title_str = paste0(dataset_labels[ds]," — ",method_name," ASR (fan)"))
    ggsave(file.path(out_dir, paste0(ds,"_",tolower(method_name),"_fan.pdf")),
           p_fan, width=16, height=16)
    ggsave(file.path(out_dir, paste0(ds,"_",tolower(method_name),"_fan.png")),
           p_fan, width=16, height=16, dpi=200)

    cat(sprintf("  Saved %s plots\n", method_name))
  }
}

cat("\nOutputs in:", out_dir, "\n")
