#!/usr/bin/env Rscript
# calibrate_chronos_correlated_325sp.R
# Re-date the 325-sp FastSpeciesTree with chronos (penalised likelihood, correlated
# rates, smoothing lambda=0.1 — the model selected by TimeTree benchmarking) using a
# constraints TSV (label, tip_a, tip_b, age_min, age_max) of 62 calibration points.
# Usage: Rscript calibrate_chronos_correlated_325sp.R <constraints.tsv> <out.nwk> [model] [lambda]
# Publication run: Rscript calibrate_chronos_correlated_325sp.R over_calib.tsv \
#                    outputs/full_325sp_calibrated_correlatedlambda01.nwk correlated 0.1
suppressPackageStartupMessages({library(ape); library(phytools)})
args <- commandArgs(trailingOnly=TRUE)
CONSTR <- if (length(args)>=1) args[1] else "over_calib.tsv"
OUT    <- if (length(args)>=2) args[2] else "outputs/full_325sp_calibrated_correlatedlambda01.nwk"
MODEL  <- if (length(args)>=3) args[3] else "correlated"   # chronos model: correlated | relaxed | discrete
LAMBDA <- if (length(args)>=4) as.numeric(args[4]) else 0.1 # smoothing parameter (0.1 = selected model)

SP  <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/01_species_tree"
CTRL<- chronos.control(iter.max=1e6, eval.max=1e6, dual.iter.max=1e4, tol=1e-8)
NRESTART <- 12

tree <- read.tree(file.path(SP,"fast_species_tree_325sp_renamed.nwk"))
tree <- multi2di(midpoint.root(tree), random=FALSE)

tip <- function(base){
  h <- grep(paste0("^",gsub("\\.","\\\\.",base),"(\\.fa|\\.fasta)?$"), tree$tip.label, value=TRUE)
  if(length(h)==0) NA_character_ else h[1]
}
node_of <- function(a,b) if(is.na(tip(a))||is.na(tip(b))) NA_integer_ else as.integer(getMRCA(tree,c(tip(a),tip(b))))

cal <- read.delim(file.path(SP,CONSTR), stringsAsFactors=FALSE)
cal$node <- mapply(function(l,a,b){
  if(l=="Eukaryota_root") return(as.integer(Ntip(tree)+1L))
  node_of(a,b)
}, cal$label, cal$tip_a, cal$tip_b)
cat(sprintf("Resolved %d/%d constraints\n", sum(!is.na(cal$node)), nrow(cal)))
cal <- cal[!is.na(cal$node),]

# collapse duplicate nodes: intersect bounds (max of mins, min of maxs)
mn <- tapply(cal$age_min, cal$node, max)
mx <- tapply(cal$age_max, cal$node, min)
nodes <- as.integer(names(mn)); keep <- mn <= mx
if(any(!keep)) cat("Dropping", sum(!keep), "nodes with inconsistent merged bounds\n")
calib <- makeChronosCalib(tree, node=nodes[keep], age.min=mn[keep], age.max=mx[keep])
cat("Calibration nodes:", nrow(calib), "\n")

# run chronos with restarts; keep the CONVERGED run with best (lowest) PHIIC
best <- NULL; best_phiic <- Inf; best_conv <- FALSE
for (i in seq_len(NRESTART)) {
  conv <- TRUE
  fiti <- withCallingHandlers(
    tryCatch(chronos(tree, lambda=LAMBDA, model=MODEL, calibration=calib, control=CTRL),
             error=function(e){conv<<-FALSE; NULL}),
    warning=function(w){ if(grepl("without convergence", conditionMessage(w))) conv<<-FALSE
                         invokeRestart("muffleWarning") })
  if (is.null(fiti)) next
  ph <- attr(fiti, "PHIIC")$PHIIC; if (is.null(ph)) ph <- attr(fiti,"PHIIC")
  cat(sprintf("  restart %2d: PHIIC=%.2f  converged=%s\n", i, ph, conv))
  # prefer converged runs; among same class pick lowest PHIIC
  better <- (conv && !best_conv) || (conv==best_conv && ph < best_phiic)
  if (better) { best <- fiti; best_phiic <- ph; best_conv <- conv }
  if (best_conv) break   # first fully-converged run is sufficient
}
if (is.null(best)) stop("chronos failed on all restarts")
cat(sprintf("CHOSEN: PHIIC=%.2f  converged=%s\n", best_phiic, best_conv))
if (!best_conv) cat("WARNING: no fully-converged run found\n")
fit <- as.phylo(best)
fit$tip.label <- sub("\\.(fa|fasta)$","", fit$tip.label)   # strip .fa to match published tips
outpath <- if (startsWith(OUT,"/")) OUT else file.path(SP, OUT)
write.tree(fit, outpath)
# also write a .fa-tipped sibling (kept in sync for plot_tree_v5_chronos.R)
fit_fa <- fit; fit_fa$tip.label <- paste0(fit_fa$tip.label, ".fa")
write.tree(fit_fa, sub("\\.nwk$", "_fa.nwk", outpath))
cat("root age:", round(max(node.depth.edgelength(fit)),1), " | Saved:", outpath,
    "(+ _fa sibling)\n")
