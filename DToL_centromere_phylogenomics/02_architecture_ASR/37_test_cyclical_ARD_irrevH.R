#!/usr/bin/env Rscript
# 37_test_cyclical_ARD_irrevH.R
#
# Alpha-omega cyclical model: ARD_irrevH (H irreversible, direct Sat↔Trans ALLOWED)
# Mixed is its own state — NOT assumed to be an obligate intermediate.
#
# TESTS:
# T1. ARD_irrevH vs ARD_irrevH_terminalTrans (Trans→Sat=0, Trans→Mixed=0)
# T2. Rate asymmetry: direct Sat→Trans vs Trans→Sat
# T3. Quasi-stationary distribution (non-H states, closed submatrix)
# T4. Root state posterior from 100 stochastic maps
# T5. Posterior transition counts (mean ± SD per simmap history)
# T6. Independent Trans origins/returns per simmap history

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(readxl); library(stringr)
  library(patchwork); library(scales)
})

set.seed(42)

root      <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
out_dir   <- file.path(root, "outputs", "cyclical_ARD_irrevH")
cache_dir <- file.path(out_dir, "cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive=TRUE)

state_levels <- c("H", "Mixed", "Sat", "Trans")
state_colors <- c(H="#f46d43", Mixed="#8e44ad", Sat="#ff2d87", Trans="#4f84e8")
datasets       <- c("full", "metazoa", "viridiplantae", "full_treepl", "metazoa_treepl", "viridiplantae_treepl")
dataset_labels <- c(full="Full tree", metazoa="Metazoa", viridiplantae="Viridiplantae", full_treepl="Full tree (treePL)", metazoa_treepl="Metazoa (treePL)", viridiplantae_treepl="Viridiplantae (treePL)")
aicc_calc <- function(aic, k, n) aic + 2*k*(k+1)/(n-k-1)

# ── Annotations ──────────────────────────────────────────────────────────────
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
  select(fasta, genus, species) %>%
  mutate(fasta_fixed = ifelse(fasta %in% names(label_fix_map),
                              unname(label_fix_map[fasta]), fasta),
         fasta_base  = str_replace(fasta_fixed, "\\.(fa|fasta)$", ""),
         sp_name     = paste(genus, species)) %>%
  distinct(fasta_base, .keep_all=TRUE)

# ── Design matrices ──────────────────────────────────────────────────────────
make_ard_dm <- function() {
  dm <- matrix(0L, 4, 4, dimnames=list(state_levels, state_levels))
  k <- 1L
  for (i in state_levels) for (j in state_levels)
    if (i != j) { dm[i,j] <- k; k <- k+1L }
  dm
}
# ARD_irrevH: H→X = 0 for all X; direct Sat↔Trans freely estimated
dm_cyc <- make_ard_dm()
dm_cyc["H","Mixed"] <- 0L; dm_cyc["H","Sat"] <- 0L; dm_cyc["H","Trans"] <- 0L

# ARD_irrevH_terminalTrans: Trans is also a sink (Trans→Sat=0, Trans→Mixed=0)
dm_term <- dm_cyc
dm_term["Trans","Sat"]   <- 0L
dm_term["Trans","Mixed"] <- 0L

# ── Load dataset ─────────────────────────────────────────────────────────────
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
  ch <- factor(as.character(ch[kp]), levels=state_levels)
  names(ch) <- tr$tip.label
  list(tree=tr, char=ch, n=length(ch))
}

# ── Reconstruct Q from fitMk + design matrix ─────────────────────────────────
reconstruct_Q <- function(fit, dm) {
  Q <- matrix(0, 4, 4, dimnames=list(state_levels, state_levels))
  for (i in 1:4) for (j in 1:4)
    if (dm[i,j] > 0) Q[i,j] <- fit$rates[dm[i,j]]
  diag(Q) <- -rowSums(Q)
  Q
}

# ── Count transitions in one simmap ──────────────────────────────────────────
count_transitions <- function(sm) {
  tv <- character(0)
  for (b in seq_along(sm$maps)) {
    st <- names(sm$maps[[b]])
    if (length(st) > 1) tv <- c(tv, paste0(head(st,-1),"->",tail(st,-1)))
  }
  if (!length(tv)) return(table(character(0)))
  table(tv)
}

# ── Quasi-stationary distribution (closed non-H submatrix) ───────────────────
quasi_stationary <- function(Q) {
  nm <- c("Mixed","Sat","Trans")
  Qs <- Q[nm, nm]; diag(Qs) <- 0; diag(Qs) <- -rowSums(Qs)
  ev  <- eigen(t(Qs))
  idx <- which.min(abs(Re(ev$values)))
  pi  <- abs(Re(ev$vectors[,idx]))
  setNames(pi / sum(pi), nm)
}

# ═══════════════════════════════════════════════════════════════════════════════
# T1: Fit models
# ═══════════════════════════════════════════════════════════════════════════════
cat("=== T1: ARD_irrevH vs ARD_irrevH_terminalTrans ===\n")
Q_list <- list(); t1_res <- list()

for (ds in datasets) {
  cat("\n[", ds, "]\n")
  dat <- load_dataset(ds); if (is.null(dat)) next
  fc <- fitMk(dat$tree, dat$char, model=dm_cyc,  states=state_levels, control=list(maxit=3000))
  ft <- fitMk(dat$tree, dat$char, model=dm_term, states=state_levels, control=list(maxit=3000))
  kc <- length(unique(dm_cyc[dm_cyc   > 0]))
  kt <- length(unique(dm_term[dm_term  > 0]))
  ac <- aicc_calc(2*kc - 2*fc$logLik, kc, dat$n)
  at <- aicc_calc(2*kt - 2*ft$logLik, kt, dat$n)
  lr <- 2*(fc$logLik - ft$logLik)
  pv <- pchisq(lr, df=kc-kt, lower.tail=FALSE)
  cat(sprintf("  ARD_irrevH:    logL=%.3f  AICc=%.2f  k=%d\n", fc$logLik, ac, kc))
  cat(sprintf("  TerminalTrans: logL=%.3f  AICc=%.2f  k=%d\n", ft$logLik, at, kt))
  cat(sprintf("  LRT Δ=%.2f  p=%.4f → %s\n", lr, pv,
              ifelse(pv<0.05,"CYCLICAL SUPPORTED","not significant")))
  Q_list[[ds]] <- reconstruct_Q(fc, dm_cyc)
  t1_res[[ds]] <- list(ds=ds, logL_c=fc$logLik, logL_t=ft$logLik,
                       kc=kc, kt=kt, ac=ac, at=at, dAICc=at-ac, lr=lr, p=pv)
}

cat("\n=== Fitted Q matrices ===\n")
for (ds in datasets) { cat("\n",dataset_labels[ds],":\n",sep=""); print(round(Q_list[[ds]],4)) }

t1_df <- bind_rows(lapply(t1_res, function(r)
  data.frame(dataset=dataset_labels[r$ds], logL_cyc=round(r$logL_c,3),
             logL_term=round(r$logL_t,3), k_cyc=r$kc, k_term=r$kt,
             AICc_cyc=round(r$ac,2), AICc_term=round(r$at,2),
             dAICc=round(r$dAICc,2), LRT=round(r$lr,3), LRT_p=round(r$p,4),
             cyclical_supported=r$p<0.05)))
print(t1_df)
write.table(t1_df, file.path(out_dir,"T1_model_comparison.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# T2: Rate asymmetry
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== T2: Rate asymmetry ===\n")
t2_df <- bind_rows(lapply(datasets, function(ds) {
  Q <- Q_list[[ds]]; if(is.null(Q)) return(NULL)
  data.frame(dataset=dataset_labels[ds],
             Sat_to_Trans=Q["Sat","Trans"], Trans_to_Sat=Q["Trans","Sat"],
             Sat_to_Mixed=Q["Sat","Mixed"], Mixed_to_Sat=Q["Mixed","Sat"],
             Trans_to_Mixed=Q["Trans","Mixed"], Mixed_to_Trans=Q["Mixed","Trans"],
             direct_ratio=ifelse(Q["Trans","Sat"]>1e-15,Q["Sat","Trans"]/Q["Trans","Sat"],NA_real_))
}))
print(t2_df)
write.table(t2_df, file.path(out_dir,"T2_rate_asymmetry.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# T3: Quasi-stationary distribution
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== T3: Quasi-stationary distribution ===\n")
t3_df <- bind_rows(lapply(datasets, function(ds) {
  Q <- Q_list[[ds]]; if(is.null(Q)) return(NULL)
  pi <- quasi_stationary(Q)
  data.frame(dataset=dataset_labels[ds], state=names(pi), proportion=pi)
}))
print(t3_df)
write.table(t3_df, file.path(out_dir,"T3_stationary_dist.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Stochastic maps
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Running stochastic maps (100 per dataset) ===\n")
sm_all <- list()
for (ds in datasets) {
  cf <- file.path(cache_dir, paste0("simmap_ARD_irrevH_", ds, ".rds"))
  if (file.exists(cf)) {
    cat("  Loading cached:", ds, "\n"); sm_all[[ds]] <- readRDS(cf); next
  }
  dat <- load_dataset(ds); Q_ds <- Q_list[[ds]]
  if (is.null(dat) || is.null(Q_ds)) next
  cat("  Running 100 maps for", ds, "...\n")
  sm <- make.simmap(dat$tree, dat$char, model="fixed", Q=Q_ds, nsim=100, message=FALSE)
  saveRDS(sm, cf); sm_all[[ds]] <- sm
  cat("  Done.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# T4: Root state posterior
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== T4: Root state posterior ===\n")
t4_df <- bind_rows(lapply(datasets, function(ds) {
  sm <- sm_all[[ds]]; if(is.null(sm)) return(NULL)
  sm_sum <- summary(sm)
  root_row <- sm_sum$ace[length(sm[[1]]$tip.label)+1, , drop=FALSE]
  data.frame(dataset=dataset_labels[ds],
             state=colnames(root_row),
             posterior=as.numeric(root_row[1,]))
}))
print(t4_df)
write.table(t4_df, file.path(out_dir,"T4_root_state.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# T5: Posterior transition counts
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== T5: Posterior transition counts ===\n")
t5_df <- bind_rows(lapply(datasets, function(ds) {
  sm <- sm_all[[ds]]; if(is.null(sm)) return(NULL)
  tabs <- lapply(sm, count_transitions)
  all_keys <- sort(unique(unlist(lapply(tabs, names))))
  if (!length(all_keys)) return(NULL)
  mat <- do.call(rbind, lapply(tabs, function(t) {
    v <- setNames(rep(0, length(all_keys)), all_keys); v[names(t)] <- as.numeric(t); v
  }))
  data.frame(dataset=dataset_labels[ds], transition=all_keys,
             mean=round(colMeans(mat),2), sd=round(apply(mat,2,sd),2),
             q025=round(apply(mat,2,quantile,0.025),1),
             q975=round(apply(mat,2,quantile,0.975),1))
}))
print(t5_df)
write.table(t5_df, file.path(out_dir,"T5_transition_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# T6: Trans origins and returns per history
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== T6: Trans origins and returns ===\n")
t6_df <- bind_rows(lapply(datasets, function(ds) {
  sm <- sm_all[[ds]]; if(is.null(sm)) return(NULL)
  origins <- sapply(sm, function(s) {
    tv <- character(0)
    for (b in seq_along(s$maps)) { st <- names(s$maps[[b]])
      if (length(st)>1) tv <- c(tv, paste0(head(st,-1),"->",tail(st,-1))) }
    sum(tv %in% c("Sat->Trans","Mixed->Trans"))
  })
  returns <- sapply(sm, function(s) {
    tv <- character(0)
    for (b in seq_along(s$maps)) { st <- names(s$maps[[b]])
      if (length(st)>1) tv <- c(tv, paste0(head(st,-1),"->",tail(st,-1))) }
    sum(tv %in% c("Trans->Sat","Trans->Mixed"))
  })
  data.frame(dataset=dataset_labels[ds],
             mean_origins=round(mean(origins),2), sd_origins=round(sd(origins),2),
             mean_returns=round(mean(returns),2), sd_returns=round(sd(returns),2),
             frac_with_return=round(mean(returns>0),3))
}))
print(t6_df)
write.table(t6_df, file.path(out_dir,"T6_trans_origins.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

# T1: AICc comparison
p_t1 <- t1_df %>%
  pivot_longer(c(AICc_cyc, AICc_term), names_to="model", values_to="AICc") %>%
  mutate(
    model_lab = ifelse(model=="AICc_cyc","ARD_irrevH\n(Trans→Sat allowed)",
                       "TerminalTrans\n(Trans is sink)"),
    sig_lab = case_when(LRT_p<0.001~"*** p<0.001", LRT_p<0.01~"** p<0.01",
                        LRT_p<0.05~"* p<0.05", TRUE~sprintf("p=%.3f",LRT_p))
  ) %>%
  ggplot(aes(x=model_lab, y=AICc, fill=model_lab)) +
  facet_wrap(~dataset, nrow=1, scales="free_y") +
  geom_col(width=0.55, alpha=0.9) +
  geom_text(data=~filter(., model=="AICc_cyc"),
            aes(x=1.5, y=-Inf, label=sig_lab), inherit.aes=FALSE,
            vjust=-0.5, size=3.5, color="#333333", fontface="italic") +
  scale_fill_manual(values=c("ARD_irrevH\n(Trans→Sat allowed)"="#2171b5",
                              "TerminalTrans\n(Trans is sink)"="#cb181d"), guide="none") +
  labs(title="T1: Cyclical (ARD_irrevH) vs Terminal-Trans — AICc comparison",
       subtitle="Lower = better. Cyclical allows Trans→Sat; Terminal does not.",
       x=NULL, y="AICc") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"T1_model_comparison.pdf"), p_t1, width=10, height=5)
ggsave(file.path(out_dir,"T1_model_comparison.png"), p_t1, width=10, height=5, dpi=300)

# T2: Rate asymmetry
t2_long <- t2_df %>%
  select(dataset, Sat_to_Trans, Trans_to_Sat, Sat_to_Mixed,
         Mixed_to_Sat, Trans_to_Mixed, Mixed_to_Trans) %>%
  pivot_longer(-dataset, names_to="rate_name", values_to="rate") %>%
  filter(rate > 1e-15) %>%
  mutate(
    trans_lab = recode(rate_name,
      Sat_to_Trans="Sat→Trans", Trans_to_Sat="Trans→Sat",
      Sat_to_Mixed="Sat→Mixed", Mixed_to_Sat="Mixed→Sat",
      Trans_to_Mixed="Trans→Mixed", Mixed_to_Trans="Mixed→Trans"),
    is_direct = rate_name %in% c("Sat_to_Trans","Trans_to_Sat"),
    direction = ifelse(rate_name %in% c("Sat_to_Trans","Sat_to_Mixed","Mixed_to_Trans"),
                       "Forward (→Trans)", "Return (→Sat)"),
    fill_col  = case_when(
      rate_name=="Sat_to_Trans"   ~ "#d73027",
      rate_name=="Trans_to_Sat"   ~ "#313695",
      direction=="Forward (→Trans)" ~ "#fc8d59",
      TRUE                           ~ "#74add1"
    )
  )

p_t2 <- ggplot(t2_long, aes(x=trans_lab, y=rate, fill=I(fill_col))) +
  facet_wrap(~dataset, nrow=1, scales="free_x") +
  geom_col(width=0.65, alpha=0.9) +
  geom_text(aes(label=sprintf("%.4f", rate)), vjust=-0.3, size=2.8) +
  scale_y_log10(labels=label_number(accuracy=0.0001)) +
  labs(title="T2: Fitted transition rates — ARD_irrevH",
       subtitle="Dark bars = direct Sat↔Trans (the cycle). Orange = forward, blue = return.",
       x=NULL, y="Rate (per MY, log scale)") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"T2_rate_asymmetry.pdf"), p_t2, width=13, height=6)
ggsave(file.path(out_dir,"T2_rate_asymmetry.png"), p_t2, width=13, height=6, dpi=300)

# T3: Quasi-stationary distribution
p_t3 <- t3_df %>%
  mutate(state=factor(state, levels=c("Sat","Mixed","Trans"))) %>%
  ggplot(aes(x=state, y=proportion, fill=state)) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(width=0.6, alpha=0.9) +
  geom_text(aes(label=sprintf("%.1f%%", proportion*100)), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=state_colors[c("Sat","Mixed","Trans")], guide="none") +
  scale_y_continuous(labels=percent_format(), expand=expansion(mult=c(0,0.15))) +
  labs(title="T3: Quasi-stationary distribution (non-H states)",
       subtitle="Expected relative time in each state at stationarity",
       x=NULL, y="Proportion") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"T3_stationary_dist.pdf"), p_t3, width=10, height=5)
ggsave(file.path(out_dir,"T3_stationary_dist.png"), p_t3, width=10, height=5, dpi=300)

# T4: Root state posterior
p_t4 <- t4_df %>%
  mutate(state=factor(state, levels=state_levels)) %>%
  filter(!is.na(state)) %>%
  ggplot(aes(x=state, y=posterior, fill=state)) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(width=0.6, alpha=0.9) +
  geom_text(aes(label=sprintf("%.2f", posterior)), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=state_colors, guide="none") +
  scale_y_continuous(limits=c(0,1.15)) +
  labs(title="T4: Root state posterior probability (100 stochastic maps)",
       x=NULL, y="Posterior probability") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"T4_root_state.pdf"), p_t4, width=10, height=5)
ggsave(file.path(out_dir,"T4_root_state.png"), p_t4, width=10, height=5, dpi=300)

# T5: Transition counts — focus on key transitions
key_tr <- c("Sat->Trans","Trans->Sat","Sat->Mixed","Mixed->Sat",
            "Trans->Mixed","Mixed->Trans")
t5_plot <- t5_df %>%
  filter(transition %in% key_tr) %>%
  mutate(
    trans_lab = recode(transition,
      "Sat->Trans"="Sat→Trans","Trans->Sat"="Trans→Sat",
      "Sat->Mixed"="Sat→Mixed","Mixed->Sat"="Mixed→Sat",
      "Trans->Mixed"="Trans→Mixed","Mixed->Trans"="Mixed→Trans"),
    type = case_when(
      transition %in% c("Sat->Trans","Trans->Sat") ~ "Direct cycle (Sat↔Trans)",
      TRUE ~ "Via Mixed"
    ),
    fill_col = case_when(
      transition=="Sat->Trans"   ~ "#d73027",
      transition=="Trans->Sat"   ~ "#313695",
      transition %in% c("Sat->Mixed","Mixed->Trans") ~ "#fc8d59",
      TRUE ~ "#74add1"
    )
  )

p_t5 <- ggplot(t5_plot, aes(x=trans_lab, y=mean, fill=I(fill_col))) +
  facet_wrap(~dataset, nrow=1, scales="free_x") +
  geom_col(width=0.65, alpha=0.9) +
  geom_errorbar(aes(ymin=pmax(0,mean-sd), ymax=mean+sd), width=0.25) +
  geom_text(aes(label=sprintf("%.1f",mean)), vjust=-0.7, size=3) +
  labs(title="T5: Posterior transition counts (mean ± SD, 100 simmap histories)",
       subtitle="Dark = direct Sat↔Trans cycle; light = transitions involving Mixed",
       x=NULL, y="Mean N transitions per history") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"T5_transition_counts.pdf"), p_t5, width=13, height=6)
ggsave(file.path(out_dir,"T5_transition_counts.png"), p_t5, width=13, height=6, dpi=300)

# T6: Origins vs returns
p_t6 <- t6_df %>%
  pivot_longer(c(mean_origins,mean_returns), names_to="type", values_to="mean") %>%
  mutate(sd_val=ifelse(type=="mean_origins",sd_origins,sd_returns),
         type_lab=recode(type,mean_origins="→Trans\n(entries)",mean_returns="Trans→\n(exits)"),
         fill_col=ifelse(type=="mean_origins","#d73027","#313695")) %>%
  ggplot(aes(x=type_lab, y=mean, fill=I(fill_col))) +
  facet_wrap(~dataset, nrow=1) +
  geom_col(width=0.55, alpha=0.9) +
  geom_errorbar(aes(ymin=pmax(0,mean-sd_val), ymax=mean+sd_val), width=0.25) +
  geom_text(aes(label=sprintf("%.1f",mean)), vjust=-0.5, size=3.5) +
  geom_text(data=t6_df, aes(x=1.5, y=0, label=sprintf("%.0f%% maps\nhave return",
                                                         frac_with_return*100)),
            inherit.aes=FALSE, vjust=-0.3, size=3, color="#555555", fontface="italic") +
  labs(title="T6: Trans state entries and exits per simmap history",
       subtitle="Entries = Sat→Trans + Mixed→Trans;  Exits = Trans→Sat + Trans→Mixed",
       x=NULL, y="Mean N events per history") +
  theme_bw(base_size=12) +
  theme(strip.text=element_text(face="bold"), panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"T6_trans_origins.pdf"), p_t6, width=10, height=5)
ggsave(file.path(out_dir,"T6_trans_origins.png"), p_t6, width=10, height=5, dpi=300)

# Combined
p_combined <- (p_t1 / p_t2 / (p_t3 | p_t4) / (p_t5 | p_t6)) +
  plot_annotation(title="Alpha-omega cyclical model tests (ARD_irrevH, direct Sat↔Trans)",
                  theme=theme(plot.title=element_text(face="bold", size=14)))

ggsave(file.path(out_dir,"combined_cyclical_tests.pdf"), p_combined, width=16, height=22)
ggsave(file.path(out_dir,"combined_cyclical_tests.png"), p_combined, width=16, height=22, dpi=300)

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n========= SUMMARY =========\n")
cat("\nT1 — Cyclical vs Terminal-Trans:\n")
for (r in t1_res) cat(sprintf("  %-14s: dAICc=%+.2f  LRT p=%.4f → %s\n",
    dataset_labels[r$ds], r$dAICc, r$p, ifelse(r$p<0.05,"CYCLICAL SUPPORTED","not significant")))

cat("\nT2 — Direct Sat↔Trans rates:\n")
for (i in seq_len(nrow(t2_df)))
  cat(sprintf("  %-14s: Sat→Trans=%.5f  Trans→Sat=%.5f  ratio=%.2f\n",
      t2_df$dataset[i], t2_df$Sat_to_Trans[i], t2_df$Trans_to_Sat[i],
      t2_df$direct_ratio[i]))

cat("\nT5 — Mean N(Sat→Trans) and N(Trans→Sat):\n")
t5_cyc <- t5_df %>% filter(transition %in% c("Sat->Trans","Trans->Sat"))
for (i in seq_len(nrow(t5_cyc)))
  cat(sprintf("  %-14s  %s: mean=%.2f [%.1f–%.1f]\n",
      t5_cyc$dataset[i], t5_cyc$transition[i], t5_cyc$mean[i],
      t5_cyc$q025[i], t5_cyc$q975[i]))

cat("\nT6 — Trans entries/exits:\n")
print(t6_df %>% select(dataset, mean_origins, mean_returns, frac_with_return))

cat("\nOutputs in:", out_dir, "\n")
