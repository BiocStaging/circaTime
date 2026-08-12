# T11 (revision_and_fix.md): Table 1, Table 4, and both panels of Figure 1
# ("fig:engineconsensus") currently come from two different cross-validation
# runs (run_complete_case_benchmark.R: 4 engines, nFolds=6; run_consensus_
# benchmark.R: 2 engines, nFolds=5), which is why the figure caption has to
# apologise that per-tissue numbers for the same engine "are not expected to
# match exactly" between panels.
#
# This script does NOT re-run cross-validation. T6's already-committed
# revision_persample_residuals_GSE54650.csv IS a single CV run (nFolds=6,
# seed=1, all four engines scored together per tissue, same seed -> same donor
# folds regardless of which methods were requested) - the per-sample record T2/
# T6 already needed. Consensus is recomputed directly from that one file's
# Molecular Timetable + tauFisher predictions (the only two engines with
# bundled pre-trained models, so the only pair a real consensus call could
# combine - see estimatePhase()'s Consensus section), using
# .combineEnginePhases(), so single engines and consensus now share one fold
# assignment and there is nothing left to reconcile between panels.
# Run from the project root: Rscript benchmarks/run_unified_benchmark.R

suppressMessages(library(pkgload))
pkgload::load_all(".", quiet = TRUE)

resid <- utils::read.csv("benchmarks/results/revision_persample_residuals_GSE54650.csv",
                          stringsAsFactors = FALSE)
resid <- resid[!resid$fold_failed, ]

allEngines <- c("timetable", "zeitzeiger", "timesignatr", "taufisher")
consensusEngines <- c("timetable", "taufisher")

unifiedRows <- list()

# Complete-case sample set (same definition as run_complete_case_benchmark.R /
# Table 1): the samples every one of the four engines produced a prediction
# for. Table 4 / Panel B are restricted to this same set - not just the
# timetable+taufisher intersection - so that a given tissue's Molecular
# Timetable or tauFisher number is IDENTICAL in Table 1/Panel A and Table
# 4/Panel B, per this task's requirement to score consensus "on the same
# held-out samples as everything else".
okByEngine <- lapply(allEngines, function(m) {
  s <- resid[resid$method == m & !is.na(resid$pred), ]
  split(s$sample, s$tissue)
})
names(okByEngine) <- allEngines
completeCaseSamples <- lapply(unique(resid$tissue), function(tis) {
  Reduce(intersect, lapply(okByEngine, function(x) x[[tis]]))
})
names(completeCaseSamples) <- unique(resid$tissue)

# ---- Single engines: restricted to complete-case sample set ----------
# Each engine reports error only on the samples all four engines could score.
# This makes Table 1 and Table 6 single-engine values identical for shared tissue/engine.
for (m in allEngines) {
  rows <- list()
  for (tis in unique(resid$tissue)) {
    sub <- resid[resid$tissue == tis & resid$method == m & resid$sample %in% completeCaseSamples[[tis]] & !is.na(resid$pred), ]
    if (nrow(sub) > 0) {
      rows[[tis]] <- sub
    }
  }
  if (length(rows) > 0) {
    sub <- do.call(rbind, rows)
    unifiedRows[[m]] <- data.frame(
      dataset = "GSE54650", tissue = sub$tissue, fold = sub$fold, sample_id = sub$sample,
      engine_or_method = m, true_time = sub$truth, predicted_time = sub$pred,
      circular_error_h = sub$circular_error
    )
  }
}

# ---- Consensus schemes: recomputed from the SAME run's timetable+taufisher ----
# predictions, restricted to the complete-case sample set, per tissue, exactly
# as run_consensus_benchmark.R did but on this unified file instead of a
# separate nFolds=5 run.
consensusRows <- list()
perTissueMAE <- list()
for (tis in unique(resid$tissue)) {
  sub <- resid[resid$tissue == tis & resid$method %in% consensusEngines &
               resid$sample %in% completeCaseSamples[[tis]], ]
  perEngineMAE <- vapply(consensusEngines, function(m) {
    mean(sub$circular_error[sub$method == m], na.rm = TRUE)
  }, numeric(1))
  perTissueMAE[[tis]] <- perEngineMAE

  wide <- stats::reshape(
    sub[, c("sample", "method", "pred", "truth", "fold")],
    idvar = c("sample", "truth", "fold"), timevar = "method", direction = "wide"
  )
  predCols <- paste0("pred.", consensusEngines)
  if (!all(predCols %in% colnames(wide))) next
  phaseMat <- as.matrix(wide[, predCols, drop = FALSE])

  schemes <- list(
    consensus_equal = rep(1, length(consensusEngines)) / length(consensusEngines),
    consensus_accuracy = (1 / perEngineMAE) / sum(1 / perEngineMAE)
  )
  for (schemeName in names(schemes)) {
    combined <- .combineEnginePhases(phaseMat, unname(schemes[[schemeName]]))
    ok <- !is.na(combined$phase)
    consensusRows[[paste(tis, schemeName)]] <- data.frame(
      dataset = "GSE54650", tissue = tis, fold = wide$fold[ok], sample_id = wide$sample[ok],
      engine_or_method = schemeName, true_time = wide$truth[ok],
      predicted_time = combined$phase[ok],
      circular_error_h = circularError(combined$phase[ok], wide$truth[ok])
    )
  }
}
unifiedRows <- c(unifiedRows, consensusRows)

out <- do.call(rbind, unifiedRows)
rownames(out) <- NULL
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "benchmarks/results/revision_unified_benchmark.csv", row.names = FALSE)
cat("Wrote benchmarks/results/revision_unified_benchmark.csv (", nrow(out), "rows)\n")

# ---- Table 1: complete-case MAE per engine per tissue (matches T6 exactly, --
# same source file, restated here so it is visibly part of the one unified file)
cc <- utils::read.csv("benchmarks/results/revision_complete_case_MAE.csv")
cat("\n== Table 1 (complete-case MAE, from the same nFolds=6/seed=1 run) ==\n")
print(stats::aggregate(mae_complete_case ~ engine, data = cc, FUN = mean))

# ---- Table 4: consensus vs. single engines, pooled + per-tissue oracle count -
# Single-engine numbers restricted to the same complete-case sample set as
# consensus above (and as Table 1), so a tissue's Molecular Timetable/tauFisher
# number is identical across Table 1, Table 4, and both figure panels.
t4 <- do.call(rbind, lapply(consensusRows, function(d) {
  data.frame(tissue = d$tissue[1], method = d$engine_or_method[1], circular_mae = mean(d$circular_error_h))
}))
ccSingles <- out[out$engine_or_method %in% consensusEngines, ]
ccSingles <- ccSingles[mapply(function(s, t) s %in% completeCaseSamples[[t]],
                               ccSingles$sample_id, ccSingles$tissue), ]
singles4 <- stats::aggregate(circular_error_h ~ tissue + engine_or_method, data = ccSingles, FUN = mean)
colnames(singles4) <- c("tissue", "method", "circular_mae")
t4all <- rbind(t4, singles4)
rownames(t4all) <- NULL
utils::write.csv(t4all, "benchmarks/results/revision_consensus_GSE54650.csv", row.names = FALSE)

cat("\n== Table 4 (consensus vs. single engines, unified run) ==\n")
pooled4 <- stats::aggregate(circular_mae ~ method, data = t4all, FUN = mean)
print(pooled4[order(pooled4$circular_mae), ])

cmp <- do.call(rbind, lapply(split(t4all, t4all$tissue), function(d) {
  singles <- d[d$method %in% consensusEngines, ]
  bestSingle <- singles$circular_mae[which.min(singles$circular_mae)]
  data.frame(
    tissue = d$tissue[1], best_single = round(bestSingle, 3),
    best_single_method = singles$method[which.min(singles$circular_mae)],
    consensus_equal = round(d$circular_mae[d$method == "consensus_equal"], 3),
    delta_equal = round(d$circular_mae[d$method == "consensus_equal"] - bestSingle, 3),
    consensus_accuracy = round(d$circular_mae[d$method == "consensus_accuracy"], 3),
    delta_accuracy = round(d$circular_mae[d$method == "consensus_accuracy"] - bestSingle, 3)
  )
}))
cat("\nPer-tissue oracle comparison:\n")
print(cmp, row.names = FALSE)
cat(sprintf("\nTissues where consensus_accuracy beats the best single engine: %d/%d\n",
            sum(cmp$delta_accuracy < 0), nrow(cmp)))
cat(sprintf("Tissues where consensus_equal beats the best single engine: %d/%d\n",
            sum(cmp$delta_equal < 0), nrow(cmp)))

cat("\nWrote benchmarks/results/revision_consensus_GSE54650.csv\n")
