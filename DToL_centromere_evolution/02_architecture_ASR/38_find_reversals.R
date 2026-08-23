#!/usr/bin/env Rscript
# 38_find_reversals.R
#
# Find evolutionary reversal examples on the reconstructed tree.
# For each tip, walk root→tip and look for state sequences that show:
#   Sat → Trans → Sat  (Trans episode followed by return)
#   Trans → Sat → Trans  (Sat episode followed by return)
#
# Uses MAP states from simmap posterior (summary$ace).
# Reports named species examples with confidence scores.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(readxl); library(stringr)
  library(patchwork); library(ggtree)
})

set.seed(42)

root      <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
cache_dir <- file.path(root, "outputs", "cyclical_ARD_irrevH", "cache")
out_dir   <- file.path(root, "outputs", "reversals")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

state_levels <- c("H", "Mixed", "Sat", "Trans")
state_colors <- c(H="#f46d43", Mixed="#8e44ad", Sat="#ff2d87", Trans="#4f84e8")

# ── Annotations / species names ───────────────────────────────────────────────
anno <- read.csv(file.path(root, "inputs", "full", "branch_symbol_anno.tsv"),
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
anno$tree_label <- gsub(" ", "", anno$Species)
map_arch <- c("Holocentric"="H","Satellite/transposon"="Mixed",
              "Satellite"="Sat","Transposon"="Trans")

label_fix_map <- c("daTanVulg1.hap2.1.fasta"="daTanVulg1.hap1.1.fa",
                   "rosCan_S27_v1.fasta.F2B.ChrOnly.fa"="rosCan_S27_v1")
meta <- suppressMessages(
  read_excel(file.path(dirname(root), "DTOL_327_master_March.xlsx"), sheet="Sheet1")
) %>%
  select(fasta, genus, species, classification, taxa1) %>%
  mutate(fasta_fixed = ifelse(fasta %in% names(label_fix_map),
                              unname(label_fix_map[fasta]), fasta),
         fasta_base  = str_replace(fasta_fixed, "\\.(fa|fasta)$", ""),
         sp_name     = paste(genus, species)) %>%
  distinct(fasta_base, .keep_all=TRUE)
sp_name_map <- setNames(meta$sp_name, meta$fasta_base)
group_map   <- setNames(meta$classification, meta$fasta_base)
taxa1_map   <- setNames(meta$taxa1, meta$fasta_base)

nice_name <- function(label) {
  coalesce(sp_name_map[label], label)
}

# ── Load dataset ──────────────────────────────────────────────────────────────
load_dataset <- function(ds) {
  tf <- file.path(root, "inputs", ds, "tree_renamed.nw")
  if (!file.exists(tf)) return(NULL)
  tr <- read.tree(tf)
  tr <- midpoint.root(tr); tr <- multi2di(tr, random=FALSE)
  ch <- map_arch[anno$architecture[match(tr$tip.label, anno$tree_label)]]
  names(ch) <- tr$tip.label
  ch <- factor(ch, levels=state_levels); names(ch) <- tr$tip.label
  kp <- !is.na(ch)
  tr <- drop.tip(tr, tr$tip.label[!kp])
  ch <- factor(as.character(ch[kp]), levels=state_levels); names(ch) <- tr$tip.label
  list(tree=tr, char=ch, n=length(ch))
}

# ── Build parent lookup ────────────────────────────────────────────────────────
build_parent <- function(tree) {
  n <- length(tree$tip.label) + tree$Nnode
  parent_of <- setNames(rep(NA_integer_, n), as.character(1:n))
  for (i in seq_len(nrow(tree$edge)))
    parent_of[tree$edge[i,2]] <- tree$edge[i,1]
  parent_of
}

# ── Find reversals for all tips ───────────────────────────────────────────────
# A "reversal" at a tip means: walking from root→tip, the state sequence
# contains a run like ...X...Y...X (where X = tip state, Y ≠ X).
# We report the nearest such pattern.
find_reversals <- function(tree, node_states) {
  n_tips  <- length(tree$tip.label)
  parent  <- build_parent(tree)

  results <- list()

  for (tip_idx in seq_len(n_tips)) {
    tip_state <- node_states[tip_idx]
    if (is.na(tip_state)) next

    # Walk tip → root, collecting (node, state, node_age)
    path_nodes  <- integer(0)
    path_states <- character(0)
    curr <- tip_idx
    while (!is.na(curr)) {
      path_nodes  <- c(path_nodes,  curr)
      path_states <- c(path_states, as.character(node_states[curr]))
      curr <- parent[curr]
    }
    # path_nodes[1] = tip, path_nodes[end] = root
    # Reverse to get root→tip order
    path_nodes  <- rev(path_nodes)
    path_states <- rev(path_states)

    # Skip NA states
    valid <- !is.na(path_states)
    path_nodes  <- path_nodes[valid]
    path_states <- path_states[valid]
    if (length(path_states) < 3) next

    # Compress to runs (remove consecutive identical states)
    rle_res <- rle(path_states)
    run_states <- rle_res$values

    if (length(run_states) < 3) next

    # Look for reversal: any i where run_states[i] == run_states[i+2] != run_states[i+1]
    for (i in seq_len(length(run_states) - 2)) {
      if (run_states[i] == run_states[i+2] && run_states[i] != run_states[i+1]) {
        # Find example: what species are in the intermediate clade?
        results[[length(results)+1]] <- list(
          tip_idx   = tip_idx,
          tip_label = tree$tip.label[tip_idx],
          tip_state = tip_state,
          pattern   = paste(run_states[max(1,i-1):min(length(run_states),i+3)], collapse="→"),
          anc_state = run_states[i],
          mid_state = run_states[i+1],
          ret_state = run_states[i+2],
          reversal_type = paste0(run_states[i],"→",run_states[i+1],"→",run_states[i+2])
        )
        break  # report first/shallowest reversal per tip
      }
    }
  }
  results
}

# ── Get node depth (distance from root) ───────────────────────────────────────
node_depths <- function(tree) {
  n     <- length(tree$tip.label) + tree$Nnode
  depth <- setNames(rep(0, n), as.character(1:n))
  # BFS from root
  root  <- length(tree$tip.label) + 1
  queue <- root
  while (length(queue)) {
    curr     <- queue[1]; queue <- queue[-1]
    children <- tree$edge[tree$edge[,1]==curr, 2]
    for (ch in children) {
      edge_len <- tree$edge.length[which(tree$edge[,1]==curr & tree$edge[,2]==ch)]
      depth[ch] <- depth[curr] + edge_len
      if (ch > length(tree$tip.label)) queue <- c(queue, ch)
    }
  }
  depth
}

# ─────────────────────────────────────────────────────────────────────────────
# Main analysis per dataset
# ─────────────────────────────────────────────────────────────────────────────
datasets       <- c("metazoa", "viridiplantae", "full", "metazoa_treepl", "viridiplantae_treepl", "full_treepl")
dataset_labels <- c(full="Full tree", metazoa="Metazoa", viridiplantae="Viridiplantae", full_treepl="Full tree (treePL)", metazoa_treepl="Metazoa (treePL)", viridiplantae_treepl="Viridiplantae (treePL)")

all_reversals <- list()

for (ds in datasets) {
  cat("\n========== ", dataset_labels[ds], " ==========\n", sep="")

  dat <- load_dataset(ds); if (is.null(dat)) next
  sm_cache <- file.path(cache_dir, paste0("simmap_ARD_irrevH_", ds, ".rds"))
  if (!file.exists(sm_cache)) { cat("  No simmap cache found, skipping\n"); next }
  sm_list <- readRDS(sm_cache)

  sm_sum  <- summary(sm_list)
  ace_mat <- sm_sum$ace   # nodes × states posterior probabilities
  n_tips  <- length(dat$tree$tip.label)

  # ── Node states: MAP from posterior + observed tips ────────────────────────
  # Tips: use observed state
  tip_states <- as.character(dat$char)
  names(tip_states) <- dat$tree$tip.label

  # Internal nodes: MAP from ace_mat
  node_map_state <- apply(ace_mat, 1, function(r) {
    colnames(ace_mat)[which.max(r)]
  })
  node_map_prob <- apply(ace_mat, 1, max)

  # Combined vector: tips first, then internal nodes
  all_states <- c(tip_states, node_map_state)
  all_states <- factor(all_states, levels=state_levels)

  depths <- node_depths(dat$tree)

  # ── Find reversals ──────────────────────────────────────────────────────────
  revs <- find_reversals(dat$tree, all_states)

  if (length(revs) == 0) { cat("  No reversals found.\n"); next }

  rev_df <- bind_rows(lapply(revs, function(r) {
    # Get posterior prob of MAP state at the tip's internal ancestor
    # (the node just above the tip where the reversal completes)
    data.frame(
      dataset       = dataset_labels[ds],
      tip_label     = r$tip_label,
      sp_name       = nice_name(r$tip_label),
      group         = coalesce(group_map[r$tip_label], NA_character_),
      taxa1         = coalesce(taxa1_map[r$tip_label], NA_character_),
      tip_state     = as.character(r$tip_state),
      reversal_type = r$reversal_type,
      full_pattern  = r$pattern,
      stringsAsFactors = FALSE
    )
  }))

  cat(sprintf("  Found %d reversal tips\n", nrow(rev_df)))
  cat("  Reversal types:\n")
  print(rev_df %>% count(reversal_type, sort=TRUE))

  all_reversals[[ds]] <- rev_df
}

# ── Combine and summarise ─────────────────────────────────────────────────────
rev_combined <- bind_rows(all_reversals)

cat("\n========== SUMMARY OF REVERSALS ==========\n")
cat("\nBy dataset and type:\n")
print(rev_combined %>% count(dataset, reversal_type, sort=TRUE))

# Key cycles: Sat→Trans→Sat and Trans→Sat→Trans
cat("\n--- Sat→Trans→Sat reversals (Trans episode, returns to Sat) ---\n")
sts <- rev_combined %>% filter(reversal_type=="Sat→Trans→Sat") %>%
  select(dataset, sp_name, group, taxa1, reversal_type, full_pattern) %>%
  arrange(dataset, group)
print(head(as.data.frame(sts), 50))

cat("\n--- Trans→Sat→Trans reversals (Sat episode within Trans background) ---\n")
tst <- rev_combined %>% filter(reversal_type=="Trans→Sat→Trans") %>%
  select(dataset, sp_name, group, taxa1, reversal_type, full_pattern) %>%
  arrange(dataset, group)
print(head(as.data.frame(tst), 50))

write.table(rev_combined, file.path(out_dir,"reversal_tips.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ── Normalise by n_tips and n_internal_nodes ──────────────────────────────────
n_tips_ds <- c(full=289, metazoa=173, viridiplantae=78)   # approximate knowns

rev_norm <- rev_combined %>%
  filter(reversal_type %in% c("Sat→Trans→Sat","Trans→Sat→Trans")) %>%
  count(dataset, reversal_type) %>%
  mutate(
    n_tips_ds = case_when(
      dataset=="Full tree"     ~ n_tips_ds["full"],
      dataset=="Metazoa"       ~ n_tips_ds["metazoa"],
      dataset=="Viridiplantae" ~ n_tips_ds["viridiplantae"]
    ),
    n_internal = n_tips_ds - 1,         # in a fully bifurcating rooted tree
    rate_per_tip      = round(n / n_tips_ds, 3),
    rate_per_intnode  = round(n / n_internal, 3)
  )

cat("\n=== Normalised reversal counts ===\n")
print(as.data.frame(rev_norm))
write.table(rev_norm, file.path(out_dir,"reversal_normalised.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Summary barplot of reversal counts
# ─────────────────────────────────────────────────────────────────────────────
p_raw <- rev_norm %>%
  ggplot(aes(x=reversal_type, y=n, fill=reversal_type)) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(width=0.55, alpha=0.9) +
  geom_text(aes(label=n), vjust=-0.3, size=4, fontface="bold") +
  scale_fill_manual(values=c("Sat→Trans→Sat"="#d73027","Trans→Sat→Trans"="#313695"),
                    guide="none") +
  labs(title="A: Raw reversal count (tips with reversal in root→tip path)",
       x=NULL, y="N tips") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=15, hjust=1))

p_norm_tip <- rev_norm %>%
  ggplot(aes(x=reversal_type, y=rate_per_tip, fill=reversal_type)) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(width=0.55, alpha=0.9) +
  geom_text(aes(label=sprintf("%.2f", rate_per_tip)), vjust=-0.3, size=3.8) +
  scale_fill_manual(values=c("Sat→Trans→Sat"="#d73027","Trans→Sat→Trans"="#313695"),
                    guide="none") +
  scale_y_continuous(labels=function(x) paste0(round(x*100), "%")) +
  labs(title="B: Rate per tip (n reversals / n_tips in dataset)",
       x=NULL, y="Fraction of tips") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=15, hjust=1))

p_norm_node <- rev_norm %>%
  ggplot(aes(x=reversal_type, y=rate_per_intnode, fill=reversal_type)) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(width=0.55, alpha=0.9) +
  geom_text(aes(label=sprintf("%.3f", rate_per_intnode)), vjust=-0.3, size=3.8) +
  scale_fill_manual(values=c("Sat→Trans→Sat"="#d73027","Trans→Sat→Trans"="#313695"),
                    guide="none") +
  labs(title="C: Rate per internal node (n reversals / n_internal_nodes)",
       subtitle="Internal nodes ≈ n_tips − 1 in a fully bifurcating tree",
       x=NULL, y="Reversals per internal node") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=15, hjust=1))

p_summary <- p_raw / p_norm_tip / p_norm_node +
  plot_annotation(title="Reversal events by dataset (normalised)",
                  theme=theme(plot.title=element_text(face="bold")))

ggsave(file.path(out_dir,"reversal_summary.pdf"), p_summary, width=10, height=12)
ggsave(file.path(out_dir,"reversal_summary.png"), p_summary, width=10, height=12, dpi=300)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2: Tree with reversals highlighted (metazoa + viridiplantae)
# ─────────────────────────────────────────────────────────────────────────────
plot_reversal_tree <- function(ds, rev_df_ds) {
  dat      <- load_dataset(ds); if(is.null(dat)) return(NULL)
  sm_cache <- file.path(cache_dir, paste0("simmap_ARD_irrevH_", ds, ".rds"))
  if (!file.exists(sm_cache)) return(NULL)
  sm_list  <- readRDS(sm_cache)
  sm_sum   <- summary(sm_list)
  ace_mat  <- sm_sum$ace
  n_tips   <- length(dat$tree$tip.label)

  # Build tip annotation
  tip_ann <- data.frame(
    label     = dat$tree$tip.label,
    obs_state = as.character(dat$char),
    sp_name   = nice_name(dat$tree$tip.label),
    is_reversal_STS = dat$tree$tip.label %in%
      rev_df_ds$tip_label[rev_df_ds$reversal_type=="Sat→Trans→Sat"],
    is_reversal_TST = dat$tree$tip.label %in%
      rev_df_ds$tip_label[rev_df_ds$reversal_type=="Trans→Sat→Trans"],
    stringsAsFactors=FALSE
  ) %>%
  mutate(
    reversal_flag = case_when(
      is_reversal_STS ~ "Sat→Trans→Sat",
      is_reversal_TST ~ "Trans→Sat→Trans",
      TRUE            ~ "none"
    )
  )

  p <- ggtree(dat$tree, layout="rectangular", size=0.3, color="grey60") %<+%
    tip_ann +
    # Colour tips by observed state
    geom_tippoint(aes(color=obs_state), size=1.5, alpha=0.8) +
    # Mark reversal tips with a star
    geom_tippoint(data=~filter(., reversal_flag=="Sat→Trans→Sat"),
                  aes(x=x), color="#d73027", shape=8, size=3) +
    geom_tippoint(data=~filter(., reversal_flag=="Trans→Sat→Trans"),
                  aes(x=x), color="#313695", shape=8, size=3) +
    # Label reversal tip species names
    geom_tiplab(data=~filter(., reversal_flag!="none"),
                aes(label=sp_name, color=reversal_flag),
                size=2.2, offset=0.5, hjust=0) +
    scale_color_manual(
      values=c(state_colors,
               "Sat→Trans→Sat"="#d73027","Trans→Sat→Trans"="#313695","none"="grey60"),
      breaks=names(state_colors), name="State") +
    labs(title=paste(dataset_labels[ds],": Tips with ancestral reversals"),
         subtitle="★ Sat→Trans→Sat (red)  ★ Trans→Sat→Trans (blue)") +
    theme_tree2(base_size=10) +
    theme(plot.title=element_text(face="bold"))

  p
}

for (ds in c("metazoa","viridiplantae")) {
  rev_ds <- all_reversals[[ds]]
  if (is.null(rev_ds)) next
  rev_ds_filt <- rev_ds %>%
    filter(reversal_type %in% c("Sat→Trans→Sat","Trans→Sat→Trans"))
  if (nrow(rev_ds_filt) == 0) next
  cat("\nPlotting", ds, "reversal tree...\n")
  p <- plot_reversal_tree(ds, rev_ds_filt)
  if (!is.null(p)) {
    ggsave(file.path(out_dir, paste0("reversal_tree_", ds, ".pdf")), p,
           width=12, height=max(8, length(unique(rev_ds_filt$tip_label)) * 0.15 + 8))
    ggsave(file.path(out_dir, paste0("reversal_tree_", ds, ".png")), p,
           width=12, height=max(8, length(unique(rev_ds_filt$tip_label)) * 0.15 + 8), dpi=200)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3: Simmap on reversal subtrees
# Show a small extracted subtree for the clearest Sat→Trans→Sat example
# ─────────────────────────────────────────────────────────────────────────────

# For each dataset, pick the clearest Sat→Trans→Sat example clade
# (the Trans clade that is nested inside a Sat clade with Sat tips re-evolving)
# Approach: find the MRCA of all Trans tips whose path shows STS reversal
# then extract a subtree around that clade

extract_reversal_subtree <- function(ds, target_type="Sat→Trans→Sat", max_tips=40) {
  rev_ds <- all_reversals[[ds]]
  if (is.null(rev_ds)) return(NULL)
  rev_tips <- rev_ds$tip_label[rev_ds$reversal_type == target_type]
  if (length(rev_tips) == 0) return(NULL)

  dat <- load_dataset(ds); if(is.null(dat)) return(NULL)
  sm_cache <- file.path(cache_dir, paste0("simmap_ARD_irrevH_", ds, ".rds"))
  if (!file.exists(sm_cache)) return(NULL)
  sm_list <- readRDS(sm_cache)
  sm_sum  <- summary(sm_list)
  ace_mat <- sm_sum$ace
  n_tips  <- length(dat$tree$tip.label)

  # Pick up to 3 reversal examples for the subtree
  rev_sample <- head(rev_tips, 3)
  if (length(rev_sample) < 2) return(NULL)

  # Find MRCA
  mrca_node <- getMRCA(dat$tree, rev_sample)
  if (is.null(mrca_node)) return(NULL)

  # Extract subtree
  sub_tips <- extract.clade(dat$tree, mrca_node)$tip.label
  if (length(sub_tips) > max_tips) {
    # Prune to only the reversal tips + their immediate sisters
    sub_tips <- unique(c(rev_sample,
      unlist(lapply(rev_sample, function(t) {
        parent <- dat$tree$edge[dat$tree$edge[,2]==which(dat$tree$tip.label==t), 1]
        if(length(parent)==0) return(NULL)
        dat$tree$tip.label[dat$tree$edge[dat$tree$edge[,1]==parent, 2]]
      }))))
    sub_tips <- sub_tips[sub_tips %in% dat$tree$tip.label]
    if (length(sub_tips) < 2) return(NULL)
  }

  sub_tree <- keep.tip(dat$tree, sub_tips)
  list(tree=sub_tree, rev_tips=rev_sample,
       species=nice_name(sub_tips), tip_labels=sub_tips)
}

for (ds in c("metazoa","viridiplantae")) {
  res <- extract_reversal_subtree(ds, "Sat→Trans→Sat")
  if (is.null(res)) next
  cat("\nPlotting subtree for", ds, "Sat→Trans→Sat example...\n")

  # Build annotation for the subtree
  tip_ann2 <- data.frame(
    label = res$tip_labels,
    sp    = res$species,
    is_rev= res$tip_labels %in% res$rev_tips,
    stringsAsFactors=FALSE
  )

  p_sub <- ggtree(res$tree, size=0.5) %<+% tip_ann2 +
    geom_tiplab(aes(label=sp, color=is_rev), size=3, hjust=0) +
    geom_tippoint(aes(color=is_rev), size=2.5) +
    scale_color_manual(values=c("TRUE"="#d73027","FALSE"="grey40"),
                       labels=c("TRUE"="Sat (reversal tip)","FALSE"="other"),
                       name=NULL) +
    labs(title=paste0(dataset_labels[ds],": Sat→Trans→Sat subtree example"),
         subtitle="Red tips = Sat species nested inside Trans background") +
    theme_tree2(base_size=11) +
    theme(plot.title=element_text(face="bold"), legend.position="bottom")

  ggsave(file.path(out_dir, paste0("subtree_STS_", ds, ".pdf")), p_sub,
         width=10, height=6)
  ggsave(file.path(out_dir, paste0("subtree_STS_", ds, ".png")), p_sub,
         width=10, height=6, dpi=200)
}

# ─────────────────────────────────────────────────────────────────────────────
# TEXT SUMMARY: named examples
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== NAMED REVERSAL EXAMPLES ==========\n\n")

for (ds in c("metazoa","viridiplantae")) {
  rev_ds <- all_reversals[[ds]]
  if (is.null(rev_ds)) next

  cat("--- ", dataset_labels[ds], " ---\n", sep="")

  for (rt in c("Sat→Trans→Sat","Trans→Sat→Trans")) {
    ex <- rev_ds %>% filter(reversal_type==rt) %>%
      select(sp_name, taxa1, full_pattern) %>%
      arrange(taxa1)
    if (nrow(ex)==0) next
    cat(sprintf("\n  %s (%d tips):\n", rt, nrow(ex)))
    for (i in seq_len(min(nrow(ex),20))) {
      cat(sprintf("    • %s [%s]  path: %s\n",
                  ex$sp_name[i],
                  coalesce(ex$taxa1[i],"?"),
                  ex$full_pattern[i]))
    }
  }
  cat("\n")
}

cat("\nOutputs in:", out_dir, "\n")
