#!/usr/bin/env Rscript
# 40_all_models_table.R
# Combined AICc table: ER, SYM, ARD, ARD_irrevH, TerminalTrans
# across full / metazoa / viridiplantae datasets.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(readxl); library(stringr)
})

root     <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
out_dir  <- file.path(root, "outputs", "all_models")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

state_levels <- c("H","Mixed","Sat","Trans")
datasets       <- c("full_chronos_correlated","metazoa_chronos_correlated","viridiplantae_chronos_correlated")
dataset_labels <- c(full_chronos_correlated="Full tree (chronos correlated)", metazoa_chronos_correlated="Metazoa (chronos correlated)", viridiplantae_chronos_correlated="Viridiplantae (chronos correlated)")

aicc_calc <- function(logL, k, n) {
  aic <- 2*k - 2*logL
  aic + 2*k*(k+1)/(n - k - 1)
}

# ── Annotations ──────────────────────────────────────────────────────────────
anno <- read.csv(file.path(root,"inputs","full","branch_symbol_anno.tsv"),
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
anno$tree_label <- gsub(" ","",anno$Species)
map_arch <- c("Holocentric"="H","Satellite/transposon"="Mixed",
              "Satellite"="Sat","Transposon"="Trans")
label_fix_map <- c("daTanVulg1.hap2.1.fasta"="daTanVulg1.hap1.1.fa",
                   "rosCan_S27_v1.fasta.F2B.ChrOnly.fa"="rosCan_S27_v1")

load_dataset <- function(ds) {
  tf <- file.path(root,"inputs",ds,"tree_renamed.nw")
  if (!file.exists(tf)) return(NULL)
  tr <- read.tree(tf); tr <- midpoint.root(tr); tr <- multi2di(tr, random=FALSE)
  ch <- map_arch[anno$architecture[match(tr$tip.label, anno$tree_label)]]
  names(ch) <- tr$tip.label
  ch <- factor(ch, levels=state_levels); names(ch) <- tr$tip.label
  kp <- !is.na(ch); tr <- drop.tip(tr, tr$tip.label[!kp])
  ch <- factor(as.character(ch[kp]), levels=state_levels); names(ch) <- tr$tip.label
  list(tree=tr, char=ch, n=length(ch))
}

# ── Design matrices ───────────────────────────────────────────────────────────
# ER: all rates equal
dm_ER <- matrix(1L, 4, 4, dimnames=list(state_levels,state_levels))
diag(dm_ER) <- 0L

# SYM: i→j = j→i
dm_SYM <- matrix(0L, 4, 4, dimnames=list(state_levels,state_levels))
k <- 1L
for (i in 1:3) for (j in (i+1):4) {
  dm_SYM[i,j] <- k; dm_SYM[j,i] <- k; k <- k+1L
}

# ARD: all different
dm_ARD <- matrix(0L, 4, 4, dimnames=list(state_levels,state_levels))
k <- 1L
for (i in state_levels) for (j in state_levels)
  if (i!=j) { dm_ARD[i,j] <- k; k <- k+1L }

# ARD_irrevH: H is irreversible
dm_irrevH <- dm_ARD
dm_irrevH["H","Mixed"] <- 0L; dm_irrevH["H","Sat"] <- 0L; dm_irrevH["H","Trans"] <- 0L

# TerminalTrans: ARD_irrevH + Trans→Sat=0, Trans→Mixed=0
dm_term <- dm_irrevH
dm_term["Trans","Sat"] <- 0L; dm_term["Trans","Mixed"] <- 0L

models <- list(
  ER          = list(dm=dm_ER,    desc="Equal rates",                           k=1),
  SYM         = list(dm=dm_SYM,   desc="Symmetric rates",                       k=6),
  ARD         = list(dm=dm_ARD,   desc="All rates different",                   k=12),
  ARD_irrevH  = list(dm=dm_irrevH,desc="ARD + H irreversible",                  k=9),
  TerminalTrans=list(dm=dm_term,  desc="ARD_irrevH + Trans irreversible (sink)", k=7)
)

# ── Fit all models ────────────────────────────────────────────────────────────
results <- list()

for (ds in datasets) {
  cat("\n===", dataset_labels[ds], "===\n")
  dat <- load_dataset(ds); if(is.null(dat)) next

  for (mn in names(models)) {
    cat(sprintf("  Fitting %-15s ...", mn))
    m   <- models[[mn]]
    fit <- fitMk(dat$tree, dat$char, model=m$dm,
                 states=state_levels, control=list(maxit=3000))
    k_eff <- length(unique(m$dm[m$dm > 0]))   # effective free params
    ac    <- aicc_calc(fit$logLik, k_eff, dat$n)
    cat(sprintf(" logL=%.3f  AICc=%.2f  k=%d\n", fit$logLik, ac, k_eff))
    results[[paste(ds, mn, sep="__")]] <- data.frame(
      dataset = dataset_labels[ds],
      model   = mn,
      desc    = m$desc,
      k       = k_eff,
      logL    = round(fit$logLik, 3),
      AICc    = round(ac, 2),
      stringsAsFactors = FALSE
    )
  }
}

# ── Assemble and compute dAICc + weights ─────────────────────────────────────
tbl <- bind_rows(results) %>%
  group_by(dataset) %>%
  mutate(
    dAICc  = round(AICc - min(AICc), 2),
    weight = {
      dA <- AICc - min(AICc)
      w  <- exp(-0.5 * dA)
      round(w / sum(w), 3)
    },
    best   = dAICc == 0
  ) %>%
  ungroup() %>%
  arrange(dataset, AICc)

cat("\n\n========= COMBINED MODEL TABLE =========\n")
print(as.data.frame(tbl %>% select(dataset, model, desc, k, logL, AICc, dAICc, weight)))

write.table(tbl, file.path(out_dir,"all_models_aicc.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ── Figure: AICc heatmap / dot-plot ──────────────────────────────────────────
model_order <- c("ER","SYM","ARD","ARD_irrevH","TerminalTrans")
ds_order    <- c("Full tree","Metazoa","Viridiplantae")

p <- tbl %>%
  mutate(
    model   = factor(model,   levels=model_order),
    dataset = factor(dataset, levels=ds_order),
    label   = ifelse(best, paste0("★\n",sprintf("%.1f",AICc)),
                           sprintf("%.1f\n(+%.1f)", AICc, dAICc))
  ) %>%
  ggplot(aes(x=dataset, y=model, fill=dAICc)) +
  geom_tile(color="white", linewidth=0.8) +
  geom_text(aes(label=label, color=best), size=3.2, lineheight=0.9) +
  scale_fill_gradient(low="#2166ac", high="#f7f7f7",
                      name="ΔAICc\n(from best)", limits=c(0,NA)) +
  scale_color_manual(values=c("TRUE"="black","FALSE"="grey40"), guide="none") +
  scale_y_discrete(limits=rev(model_order),
                   labels=c(ER="ER", SYM="SYM", ARD="ARD",
                             ARD_irrevH="ARD_irrevH\n(H irreversible)",
                             TerminalTrans="TerminalTrans\n(Trans+H irreversible)")) +
  labs(title="AICc comparison across all models and datasets",
       subtitle="★ = best model per dataset.  Cell value = AICc (+ΔAICc from best)",
       x=NULL, y=NULL) +
  theme_bw(base_size=13) +
  theme(axis.text.x=element_text(face="bold"),
        axis.text.y=element_text(hjust=1),
        panel.grid=element_blank(),
        legend.position="right")

ggsave(file.path(out_dir,"all_models_aicc_heatmap.pdf"), p, width=9, height=5)
ggsave(file.path(out_dir,"all_models_aicc_heatmap.png"), p, width=9, height=5, dpi=300)

# ── AICc weight dot plot ──────────────────────────────────────────────────────
p2 <- tbl %>%
  mutate(model=factor(model, levels=rev(model_order)),
         dataset=factor(dataset, levels=ds_order)) %>%
  ggplot(aes(x=weight, y=model, color=dataset, shape=dataset)) +
  geom_vline(xintercept=0, color="grey80") +
  geom_point(size=4, alpha=0.9) +
  geom_text(aes(label=sprintf("%.3f",weight)), hjust=-0.25, size=3, show.legend=FALSE) +
  scale_color_brewer(palette="Set1", name=NULL) +
  scale_shape_manual(values=c(16,17,15), name=NULL) +
  scale_x_continuous(limits=c(0,1.1), labels=function(x) sprintf("%.2f",x)) +
  labs(title="AICc model weights",
       x="Akaike weight", y=NULL) +
  theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")

ggsave(file.path(out_dir,"all_models_weights.pdf"), p2, width=8, height=5)
ggsave(file.path(out_dir,"all_models_weights.png"), p2, width=8, height=5, dpi=300)

cat("\nOutputs in:", out_dir, "\n")
