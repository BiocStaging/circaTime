# Reproducible consensus-vs-single-engine benchmark. Not part of the package
# build (see .Rbuildignore). Run from the project root:
#   Rscript benchmarks/run_consensus_benchmark.R
#
# AGENTS.md Section 6 Step 3.1 subtask 4: benchmark method = "consensus" against
# each individual engine and report where consensus improves vs. underperforms
# the best single engine.
#
# Why this does NOT call estimatePhase(method = "consensus") in a loop: doing so
# would score the bundled consensus models on the very data they were trained on.
# Instead this reproduces the consensus rule inside a real donor-split
# cross-validation via the existing, tested benchmarkPhase() - the same design
# inst/scripts/train_reference_models.R Pass 1 uses - so every number below is a
# genuine out-of-fold prediction. The consensus weights are derived per tissue
# from the same CV run's per-engine MAE, exactly as the bundled models' are, so
# this measures the shipped rule rather than an idealised variant of it.
#
# Requires benchmarks/data/GSE54650_se.rds (run inst/scripts/ingest_GSE54650.R
# first). Slow: tauFisher's MetaCycle step dominates, as in run_benchmark.R.

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

dataPath <- "benchmarks/data/GSE54650_se.rds"
if (!file.exists(dataPath)) {
  stop("Missing ", dataPath, " - run inst/scripts/ingest_GSE54650.R first.")
}
se <- readRDS(dataPath)

# Only the engines that are actually bundled as pre-trained models can form a
# consensus a user could invoke (see estimatePhase()'s Consensus section).
methodsUsed <- c("timetable", "taufisher")
engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01,
                              nrep = 1, numbasis = 5))
)
nFolds <- 5
seed <- 1

rows <- list()
for (tis in unique(colData(se)$tissue)) {
  cat(sprintf("[%s] tissue=%s\n", format(Sys.time()), tis))
  seT <- se[, colData(se)$tissue == tis]
  bm <- benchmarkPhase(seT,
    truth_col = "phase", methods = methodsUsed, assay_name = "logExpr",
    donor_col = "donor", nFolds = nFolds, seed = seed, engineArgs = engineArgs
  )
  res <- bm$residuals[!bm$residuals$fold_failed, ]

  perEngineMAE <- vapply(methodsUsed, function(m) {
    mean(res$circular_error[res$method == m], na.rm = TRUE)
  }, numeric(1))

  wide <- stats::reshape(
    res[, c("sample", "method", "pred", "truth")],
    idvar = c("sample", "truth"), timevar = "method", direction = "wide"
  )
  predCols <- paste0("pred.", methodsUsed)
  if (!all(predCols %in% colnames(wide))) {
    cat("  skipping consensus for this tissue: not every engine produced predictions\n")
    next
  }
  phaseMat <- as.matrix(wide[, predCols, drop = FALSE])

  schemes <- list(
    consensus_equal = rep(1, length(methodsUsed)) / length(methodsUsed),
    consensus_accuracy = (1 / perEngineMAE) / sum(1 / perEngineMAE)
  )

  for (m in methodsUsed) {
    rows[[length(rows) + 1]] <- data.frame(
      tissue = tis, method = m, circular_mae = perEngineMAE[[m]],
      n = sum(!is.na(res$pred[res$method == m]))
    )
  }
  for (schemeName in names(schemes)) {
    combined <- .combineEnginePhases(phaseMat, unname(schemes[[schemeName]]))
    ok <- !is.na(combined$phase)
    rows[[length(rows) + 1]] <- data.frame(
      tissue = tis, method = schemeName,
      circular_mae = mean(circularError(combined$phase[ok], wide$truth[ok])),
      n = sum(ok)
    )
  }
}

summaryDf <- do.call(rbind, rows)
rownames(summaryDf) <- NULL

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(summaryDf, "benchmarks/results/phase3_consensus_GSE54650.csv", row.names = FALSE)

cat("\n== Phase 3 Step 3.1: consensus vs. single engines, GSE54650 ==\n")
print(summaryDf, row.names = FALSE)

cat("\nMean circular MAE by method (pooled across tissues):\n")
pooled <- stats::aggregate(circular_mae ~ method, data = summaryDf, FUN = mean)
print(pooled[order(pooled$circular_mae), ], row.names = FALSE)

# Per-tissue head-to-head: does consensus beat the best single engine, and by how
# much? Reported per tissue rather than only pooled, because the whole premise of
# consensus is that which engine wins varies by tissue.
cat("\nPer-tissue: best single engine vs. each consensus scheme (hours, negative = consensus better)\n")
cmp <- do.call(rbind, lapply(split(summaryDf, summaryDf$tissue), function(d) {
  singles <- d[d$method %in% methodsUsed, ]
  bestSingle <- singles$circular_mae[which.min(singles$circular_mae)]
  data.frame(
    tissue = d$tissue[1],
    best_single = round(bestSingle, 3),
    best_single_method = singles$method[which.min(singles$circular_mae)],
    consensus_equal = round(d$circular_mae[d$method == "consensus_equal"], 3),
    delta_equal = round(d$circular_mae[d$method == "consensus_equal"] - bestSingle, 3),
    consensus_accuracy = round(d$circular_mae[d$method == "consensus_accuracy"], 3),
    delta_accuracy = round(d$circular_mae[d$method == "consensus_accuracy"] - bestSingle, 3)
  )
}))
print(cmp, row.names = FALSE)

cat(sprintf(
  "\nTissues where consensus_accuracy beats the best single engine: %d/%d\n",
  sum(cmp$delta_accuracy < 0), nrow(cmp)
))
cat(sprintf(
  "Tissues where consensus_equal beats the best single engine: %d/%d\n",
  sum(cmp$delta_equal < 0), nrow(cmp)
))
cat(
  "\nNote: 'best single engine' is an oracle - it is the best engine chosen with",
  "hindsight per tissue, which no real user could pick in advance. The fair",
  "practical comparison is consensus vs. committing to one engine everywhere",
  "(the pooled table above).\n"
)
