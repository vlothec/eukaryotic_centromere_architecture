#!/usr/bin/env Rscript
# 40b_all_models_3state.R
# 3-state ASR model comparison (H / Sat / Trans; Mixed + Unknown pruned to NA)
# on the chronos-correlated chronogram, full 8-model set:
#   ER, SYM, ARD, ARD_irrevH, ARD_irrevH_noDirectST, ARD_irrevH_symST,
#   ARD_irrevH_noSatToTrans, ARD_irrevH_noTransToSat
# Reads this repo's own 02_asr/inputs/; writes outputs/all_models_3state/.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(ggplot2)
})

# repo 02_asr root = dir containing this script's parent
args <- commandArgs(trailingOnly = FALSE)
sf   <- sub("^--file=", "", args[grep("^--file=", args)])
asr_root <- if (length(sf)) normalizePath(file.path(dirname(sf), "..")) else getwd()
in_dir   <- file.path(asr_root, "inputs")
out_dir  <- file.path(asr_root, "outputs", "all_models_3state")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

states3   <- c("H","Sat","Trans")
map_arch3 <- c(Holocentric = "H", Satellite = "Sat", Transposon = "Trans")  # rest -> NA
datasets  <- c(full_chronos_correlated = "Full tree",
               metazoa_chronos_correlated = "Metazoa",
               viridiplantae_chronos_correlated = "Viridiplantae")

aicc_calc <- function(logL, k, n) 2*k - 2*logL + 2*k*(k+1)/(n - k - 1)

load3 <- function(ds) {
  ann <- read.delim(file.path(in_dir, ds, "branch_symbol_anno.tsv"), stringsAsFactors = FALSE)
  colnames(ann)[c(1,8)] <- c("tip","architecture")
  tr  <- read.tree(file.path(in_dir, ds, "tree_renamed.nw"))
  # haplotype-label fix: tree labels daTanVulg as hap2 but the annotation is hap1
  # (Tanacetum vulgare = Satellite); without this it drops to NA and n is 276 not 277.
  tr$tip.label <- sub("daTanVulg1.hap2.1", "daTanVulg1.hap1.1", tr$tip.label, fixed = TRUE)
  tr  <- midpoint.root(tr); tr <- multi2di(tr, random = FALSE)
  ch  <- factor(map_arch3[ann$architecture[match(tr$tip.label, ann$tip)]], levels = states3)
  names(ch) <- tr$tip.label
  kp  <- !is.na(ch); tr <- drop.tip(tr, tr$tip.label[!kp])
  ch  <- factor(as.character(ch[kp]), levels = states3); names(ch) <- tr$tip.label
  list(tree = tr, char = ch, n = length(ch))
}

# ── design matrices (unique index per free rate) ──────────────────────────────
mk_ard <- function() { dm <- matrix(0L,3,3,dimnames=list(states3,states3)); k<-1L
  for(i in states3) for(j in states3) if(i!=j){dm[i,j]<-k; k<-k+1L}; dm }
mk_irrevH <- function(){ dm<-mk_ard(); dm["H","Sat"]<-0L; dm["H","Trans"]<-0L; dm }
mk_irrevH_noDirectST <- function(){ dm<-mk_irrevH(); dm["Sat","Trans"]<-0L; dm["Trans","Sat"]<-0L; dm }
mk_irrevH_noSatToTrans <- function(){ dm<-mk_irrevH(); dm["Sat","Trans"]<-0L; dm }
mk_irrevH_noTransToSat <- function(){ dm<-mk_irrevH(); dm["Trans","Sat"]<-0L; dm }

models <- list(
  ER                      = "ER",
  SYM                     = "SYM",
  ARD                     = "ARD",
  ARD_irrevH              = mk_irrevH(),
  ARD_irrevH_noDirectST   = mk_irrevH_noDirectST(),
  ARD_irrevH_noSatToTrans = mk_irrevH_noSatToTrans(),
  ARD_irrevH_noTransToSat = mk_irrevH_noTransToSat()
)

# ── fit ───────────────────────────────────────────────────────────────────────
res <- list()
for (ds in names(datasets)) {
  cat("\n===", datasets[ds], "===\n")
  dat <- load3(ds)
  cat(sprintf("  n=%d  (H=%d Sat=%d Trans=%d)\n", dat$n,
              sum(dat$char=="H"), sum(dat$char=="Sat"), sum(dat$char=="Trans")))
  for (mn in names(models)) {
    fit <- tryCatch(fitMk(dat$tree, dat$char, model = models[[mn]],
                          states = states3, control = list(maxit = 3000)),
                    error = function(e) NULL)
    if (is.null(fit)) { cat(sprintf("  %-24s FAILED\n", mn)); next }
    logL <- as.numeric(fit$logLik)
    idx  <- fit$index.matrix; k <- length(unique(idx[!is.na(idx) & idx > 0]))
    ac   <- aicc_calc(logL, k, dat$n)
    cat(sprintf("  %-24s k=%d  logL=%.2f  AICc=%.2f\n", mn, k, logL, ac))
    res[[paste(ds,mn)]] <- data.frame(dataset = datasets[ds], model = mn,
      k = k, logL = round(logL,3), AIC = round(2*k-2*logL,2), AICc = round(ac,2),
      stringsAsFactors = FALSE)
  }
}

tbl <- bind_rows(res) %>%
  group_by(dataset) %>%
  mutate(dAICc  = round(AICc - min(AICc), 2),
         w_AICc = round(exp(-0.5*(AICc-min(AICc))) / sum(exp(-0.5*(AICc-min(AICc)))), 3),
         dAIC   = round(AIC - min(AIC), 2),
         w_AIC  = round(exp(-0.5*(AIC-min(AIC)))  / sum(exp(-0.5*(AIC-min(AIC)))),  3),
         best   = dAICc == 0) %>%
  ungroup() %>% arrange(dataset, AICc)

write.table(tbl, file.path(out_dir, "all_models_3state_aicc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n========= 3-STATE MODEL TABLE =========\n")
print(as.data.frame(tbl %>% select(dataset, model, k, logL, AICc, dAICc, w_AICc)))

# ── weight dot plot (AIC + AICc facets, red one-directional labels) ───────────
# order models by total AICc weight across datasets: best-supported at the BOTTOM
# (ARD_irrevH lowest). model_order is top->bottom = ascending total weight.
model_order <- tbl %>% group_by(model) %>% summarise(tw = sum(w_AICc), .groups = "drop") %>%
  arrange(tw) %>% pull(model)
ds_pal   <- c(`Full tree`="#E41A1C", Metazoa="#377EB8", Viridiplantae="#2E7D32")
ds_shape <- c(`Full tree`=16, Metazoa=17, Viridiplantae=15)

plt <- tbl %>%
  tidyr::pivot_longer(c(w_AIC, w_AICc), names_to="criterion", values_to="w") %>%
  mutate(criterion = ifelse(criterion=="w_AIC","AIC","AICc"),
         model     = factor(model, levels=model_order),
         dataset   = factor(dataset, levels=names(ds_pal))) %>%
  ggplot(aes(w, model, colour=dataset, shape=dataset)) +
  geom_vline(xintercept=0, colour="grey85") +
  geom_point(size=3.4, alpha=0.92) +
  scale_colour_manual(values=ds_pal, name=NULL) +
  scale_shape_manual(values=ds_shape, name=NULL) +
  scale_x_continuous(limits=c(0,1.05), breaks=c(0,.25,.5,.75,1)) +
  scale_y_discrete(limits=rev(model_order)) +
  facet_wrap(~criterion, ncol=2) +
  labs(title="3-state ASR model weights — chronos-correlated (325-sp)",
       subtitle="One-directional models (noSatToTrans/noTransToSat) carry ~0 weight: bidirectional Sat<->Trans cycling",
       x="Akaike weight", y=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        legend.position="bottom", strip.text=element_text(face="bold"),
        plot.subtitle=element_text(size=9, colour="grey40"))

ggsave(file.path(out_dir,"all_models_3state_weights.pdf"), plt, width=10, height=5)
ggsave(file.path(out_dir,"all_models_3state_weights.png"), plt, width=10, height=5, dpi=300)

# ── Rate-flow diagrams (best model per dataset) ───────────────────────────────
suppressPackageStartupMessages({ library(igraph); library(ggraph); library(patchwork) })

node_pos <- data.frame(state=states3, x=c(0,-1,1), y=c(1,0,0),
  label=c("Holo-\ncentric","Satellite","Transposon"),
  color=c("#2d7d32","#E53935","#FB8C00"), stringsAsFactors=FALSE)

fit_Q <- function(ds_id, model_spec) {
  dat <- load3(ds_id)
  fit <- fitMk(dat$tree, dat$char, model=model_spec, states=states3, control=list(maxit=3000))
  as.Qmatrix(fit)
}
flow_one <- function(Q, title_str) {
  ed <- expand.grid(from=states3, to=states3, stringsAsFactors=FALSE) %>%
    filter(from!=to) %>% mutate(rate=mapply(function(i,j) Q[i,j], from, to)) %>%
    filter(rate>1e-10) %>% mutate(rate_display=pmin(rate,1), log_rate=log10(rate))
  g <- graph_from_data_frame(ed, directed=TRUE, vertices=node_pos %>% rename(name=state))
  ggraph(g, layout="manual", x=node_pos$x, y=node_pos$y) +
    geom_edge_arc(aes(width=rate_display, color=log_rate, label=sprintf("%.4f",rate)),
      strength=0.25, arrow=arrow(length=unit(3,"mm"),type="closed"),
      end_cap=circle(7,"mm"), start_cap=circle(7,"mm"),
      label_size=2.8, label_colour="grey20", angle_calc="along",
      label_dodge=unit(2.5,"mm"), show.legend=c(width=FALSE,color=TRUE)) +
    geom_node_point(aes(color=I(color)), size=15) +
    geom_node_text(aes(label=label), size=2.7, fontface="bold", color="white", lineheight=0.85) +
    scale_edge_width_continuous(range=c(0.5,4)) +
    scale_edge_color_gradient2(low="#74add1", mid="#ffffbf", high="#d73027",
      midpoint=median(ed$log_rate), name="log10(rate/My)",
      guide=guide_edge_colorbar(title.position="top")) +
    labs(title=title_str) + coord_cartesian(xlim=c(-1.7,1.7), ylim=c(-1.5,1.5)) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          legend.position="bottom", legend.key.width=unit(1.1,"cm"))
}

model_dm <- models  # named list, same specs used for fitting
# Show ARD_irrevH for ALL datasets: it's the fullest directional model
# (H irreversible; Sat->Trans and Trans->Sat as SEPARATE rates, so the cycle
# asymmetry is visible, unlike symST/SYM). It is best (Metazoa) or within
# dAICc ~0.4 (Full) / ~2.2 (Viridiplantae) of the top model.
best_per_ds <- data.frame(
  dataset = unname(datasets),
  model   = "ARD_irrevH",
  stringsAsFactors = FALSE
)
ds_id_map <- setNames(names(datasets), unname(datasets))

flows <- lapply(seq_len(nrow(best_per_ds)), function(i) {
  dsl <- best_per_ds$dataset[i]; mn <- best_per_ds$model[i]
  Q <- fit_Q(ds_id_map[[dsl]], model_dm[[mn]])
  flow_one(Q, sprintf("%s — %s", dsl, mn))
})
p_flow <- wrap_plots(flows, nrow=1, guides="collect") &
  theme(legend.position="bottom")
p_flow <- p_flow + plot_annotation(
  title="Centromere-architecture transition rates (best 3-state model per dataset, chronos-correlated)",
  subtitle="ARD_irrevH shown for all (H irreversible; Sat->Trans vs Trans->Sat as separate rates). Bidirectional Sat<->Trans = the a-o cycle.",
  theme=theme(plot.title=element_text(face="bold",size=12),
              plot.subtitle=element_text(size=9,colour="grey40")))
ggsave(file.path(out_dir,"flow_diagrams_3state_ARD_irrevH.pdf"), p_flow, width=15, height=6)
ggsave(file.path(out_dir,"flow_diagrams_3state_ARD_irrevH.png"), p_flow, width=15, height=6, dpi=300)

cat("\nOutputs in:", out_dir, "\n")
