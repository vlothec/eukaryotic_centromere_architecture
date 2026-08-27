#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(lmtest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 02_custom_models.R <input_dir> <output_dir>")
}

input_dir <- normalizePath(args[1])
output_dir <- normalizePath(args[2], mustWork = FALSE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load tree & annotation
 tree <- read.tree(file.path(input_dir, "tree_renamed.nw"))
 tree <- midpoint.root(tree)

 anno <- read.csv(file.path(input_dir, "branch_symbol_anno.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
 anno$tree_label <- gsub(" ", "", anno$Species)

# Build tip states
char_raw <- anno$architecture[match(tree$tip.label, anno$tree_label)]
names(char_raw) <- tree$tip.label

# Map to short state labels (Unknown omitted)
map_levels <- c("Holocentric"               = "H",
                "Satellite/transposon"      = "Mixed",
                "Satellite"                 = "Sat",
                "Transposon"                = "Trans")

char_short <- map_levels[char_raw]
names(char_short) <- tree$tip.label

state_levels <- c("H","Mixed","Sat","Trans")
char_factor <- factor(char_short, levels=state_levels)
names(char_factor) <- tree$tip.label
keep <- !is.na(char_factor)
tree_full <- drop.tip(tree, tree$tip.label[!keep])
char_factor <- droplevels(char_factor[keep])
names(char_factor) <- tree_full$tip.label

cat("Tree tips:", Ntip(tree_full), "\n")
cat("Char length:", length(char_factor), "\n")
cat("Unique states:", unique(as.character(char_factor)), "\n")

# ----- Helper: build ARD design matrix (all rates free) -----
make_dm_ard <- function() {
  dm <- matrix(0, nrow=length(state_levels), ncol=length(state_levels))
  rownames(dm) <- colnames(dm) <- state_levels
  k <- 1
  for (i in state_levels) {
    for (j in state_levels) {
      if (i != j) {
        dm[i, j] <- k
        k <- k + 1
      }
    }
  }
  dm
}

# ----- ARD-constrained custom models -----
make_dm_irrev_h <- function() {
  # Constraint: H is an endpoint: H -> Mixed/Sat/Trans = 0 (irreversible holocentricity)
  dm <- make_dm_ard()
  dm["H","Mixed"] <- 0
  dm["H","Sat"] <- 0
  dm["H","Trans"] <- 0
  dm
}

make_dm_no_direct_st <- function() {
  # Constraint: no direct Sat <-> Trans
  dm <- make_dm_ard()
  dm["Sat","Trans"] <- 0
  dm["Trans","Sat"] <- 0
  dm
}

make_dm_equal_source <- function() {
  # Constraint: Sat -> H equals Trans -> H
  dm <- make_dm_ard()
  dm["Trans","H"] <- dm["Sat","H"]
  dm
}

make_dm_irrev_h_noMixedToH <- function() {
  # Constraint: H -> Mixed/Sat/Trans = 0 AND Mixed -> H = 0
  dm <- make_dm_ard()
  dm["H","Mixed"] <- 0
  dm["H","Sat"] <- 0
  dm["H","Trans"] <- 0
  dm["Mixed","H"] <- 0
  dm
}

make_dm_irrevh_nodirectst <- function() {
  # Constraint: H is endpoint (H -> Mixed/Sat/Trans = 0) AND Sat <-> Trans forbidden
  dm <- make_dm_ard()
  dm["H","Mixed"] <- 0
  dm["H","Sat"] <- 0
  dm["H","Trans"] <- 0
  dm["Sat","Trans"] <- 0
  dm["Trans","Sat"] <- 0
  dm
}

fit_one <- function(name, model) {
  cat("Fitting", name, "...\n")
  fitMk(tree_full, char_factor, model=model, states=state_levels,
        control=list(maxit=2000))
}

fit_ER  <- fit_one("ER", "ER")
fit_SYM <- fit_one("SYM", "SYM")
fit_ARD <- fit_one("ARD", "ARD")

fit_irrev_h   <- fit_one("ARD_irrevH", make_dm_irrev_h())
fit_no_direct <- fit_one("ARD_noDirectST", make_dm_no_direct_st())
fit_equal_src <- fit_one("ARD_equalSource", make_dm_equal_source())
fit_irrev_h_noMixedToH <- fit_one("ARD_irrevH_noMixedToH", make_dm_irrev_h_noMixedToH())
fit_irrevh_nodirect <- fit_one("ARD_irrevH_noDirectST", make_dm_irrevh_nodirectst())

models <- list(ER=fit_ER, SYM=fit_SYM, ARD=fit_ARD,
               ARD_irrevH=fit_irrev_h,
               ARD_noDirectST=fit_no_direct,
               ARD_equalSource=fit_equal_src,
               ARD_irrevH_noMixedToH=fit_irrev_h_noMixedToH,
               ARD_irrevH_noDirectST=fit_irrevh_nodirect)

 aic_vec  <- sapply(models, function(x) if (!is.null(x$AIC)) x$AIC else AIC(x))
 logL_vec <- sapply(models, function(x) if (!is.null(x$logLik)) as.numeric(x$logLik) else as.numeric(logLik(x)))
 k_vec    <- sapply(models, function(x) if (!is.null(x$k)) x$k else length(x$rates))

 summary_df <- data.frame(k=k_vec, logL=logL_vec, AIC=aic_vec)
 write.table(summary_df, file=file.path(output_dir, "custom_model_aic.tsv"),
             sep="\t", quote=FALSE, col.names=NA)

 aic_weights <- aic.w(aic_vec)
 write.table(aic_weights, file=file.path(output_dir, "custom_model_aic_weights.tsv"),
             sep="\t", quote=FALSE, col.names=NA)

 lrt_out <- capture.output({
   print(lrtest(fit_ER, fit_SYM))
   print(lrtest(fit_ER, fit_ARD))
   print(lrtest(fit_SYM, fit_ARD))
 })
 writeLines(lrt_out, con=file.path(output_dir, "custom_model_lrtests.txt"))

png(file.path(output_dir, "Q_matrices_variants.png"), width=2600, height=2000, res=200)
par(mfrow=c(3,3))
safe_plot <- function(obj, title) {
  tryCatch({
    plot(obj, main=title, show.zeros=FALSE)
  }, error=function(e) {
    plot.new()
    text(0.5, 0.5, paste0(title, "\n(plot failed)"), cex=0.9)
  })
}
safe_plot(fit_ER,    "ER")
safe_plot(fit_SYM,   "SYM")
safe_plot(fit_ARD,   "ARD")
safe_plot(fit_irrev_h,   "ARD + H irreversible")
safe_plot(fit_no_direct, "ARD + no direct ST")
safe_plot(fit_equal_src, "ARD + equal sources")
safe_plot(fit_irrev_h_noMixedToH, "ARD + H irrev + no Mixed->H")
safe_plot(fit_irrevh_nodirect, "ARD + H irrevers + no ST")
dev.off()

 cat("Custom model comparison complete for:", input_dir, "\n")
 cat("Outputs in:", output_dir, "\n")
