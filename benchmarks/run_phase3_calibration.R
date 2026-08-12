# Real Phase 3 gate deliverable (AGENTS.md Section 9): does the uncertainty layer's
# stated coverage actually hold, and do low-amplitude-flagged samples have
# demonstrably worse error, on real labelled data (GSE54650)?
#
# Scope of this run: the `timetable` engine only. Each timetable fit is ~0.5s, so a
# full 12-tissue x 6-fold x 50-bootstrap-replicate run costs minutes. tauFisher's
# fit is ~2.5min (real MetaCycle rhythmicity detection per Phase 2), so the same
# bootstrap design for tauFisher would cost order-of-days (nBoot x nFolds x
# nTissues fits) - that is a separate, explicitly-scoped decision, not something to
# spend unattended compute on here. Not part of the package build; see
# .Rbuildignore. Run from the project root: Rscript benchmarks/run_phase3_calibration.R
# Requires benchmarks/data/GSE54650_se.rds (run inst/scripts/ingest_GSE54650.R first).
#
# First attempt at this script used bootstrapPhaseCI() with no `residualPool`
# (training-resampling variance only) and found severe undercoverage: a stated 80%
# interval covered truth only ~48% of the time, because refitting `timetable` on
# resampled training donors barely changes the prediction - it does not capture
# real test-time/aleatoric error. Fixed by a two-pass design: pass 1 does a plain
# (non-bootstrap) cross-validated run per tissue to collect real held-out circular
# residuals; pass 2 bootstraps with `residualPool` set to the OTHER 11 tissues'
# pass-1 residuals (leave-tissue-out, so a tissue's own calibration is never
# evaluated using its own residuals - AGENTS.md 7.3's "don't validate against your
# own output" logic applied to residual pooling).

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

seFile <- "benchmarks/data/GSE54650_se.rds"
if (!file.exists(seFile)) {
  stop("Run inst/scripts/ingest_GSE54650.R first to produce ", seFile)
}
se <- readRDS(seFile)

nFolds <- 6
nBoot <- 50
level <- 0.8
seed <- 1
fitArgs <- list(nGenes = 100, minCor = 0.5)

tissues <- unique(colData(se)$tissue)

# ---- Pass 1: plain (non-bootstrap) cross-validated point estimates per tissue,
# to build real held-out residual pools for pass 2's residual-augmented bootstrap.
cat("== Pass 1: cross-validated point estimates (building residual pools) ==\n")
pass1 <- list()
for (tis in tissues) {
  seT <- se[, colData(se)$tissue == tis]
  exprMat <- assay(seT, "logExpr")
  truth <- colData(seT)$phase
  donors <- colData(seT)$donor
  fold <- .assignDonorFolds(donors, nFolds, seed = seed)

  for (k in seq_len(nFolds)) {
    trainIdx <- which(fold != k)
    testIdx <- which(fold == k)
    if (length(testIdx) == 0 || length(trainIdx) == 0) next

    fitObj <- fitMolecularTimetable(exprMat[, trainIdx, drop = FALSE], truth[trainIdx],
                                     nGenes = fitArgs$nGenes, minCor = fitArgs$minCor)
    predDf <- predictMolecularTimetable(fitObj, exprMat[, testIdx, drop = FALSE])
    pass1[[paste(tis, k)]] <- data.frame(tissue = tis, truth = truth[testIdx], pred = predDf$phase)
  }
}
pass1All <- do.call(rbind, pass1)
pass1All$residual <- circularDiff(pass1All$pred, pass1All$truth)

# ---- Pass 2: bootstrap confidence intervals, residual-augmented using the OTHER
# 11 tissues' pass-1 residuals (leave-tissue-out - a tissue's interval is never
# calibrated against its own residuals).
cat("== Pass 2: bootstrap confidence intervals (residual-augmented) ==\n")
allCI <- list()
allDraws <- list()
allAmp <- list()

for (tis in tissues) {
  cat(sprintf("[%s] tissue=%s\n", format(Sys.time()), tis))
  seT <- se[, colData(se)$tissue == tis]
  exprMat <- assay(seT, "logExpr")
  truth <- colData(seT)$phase
  donors <- colData(seT)$donor
  sampleNames <- colnames(seT)

  amp <- clockAmplitude(exprMat)
  allAmp[[tis]] <- data.frame(tissue = tis, sample = sampleNames, amp = as.numeric(amp))

  residualPool <- pass1All$residual[pass1All$tissue != tis]

  fold <- .assignDonorFolds(donors, nFolds, seed = seed)

  for (k in seq_len(nFolds)) {
    trainIdx <- which(fold != k)
    testIdx <- which(fold == k)
    if (length(testIdx) == 0 || length(trainIdx) == 0) next

    ciFold <- tryCatch(
      bootstrapPhaseCI(
        method = "timetable",
        trainMat = exprMat[, trainIdx, drop = FALSE],
        trainTime = truth[trainIdx],
        testMat = exprMat[, testIdx, drop = FALSE],
        donors = donors[trainIdx],
        nBoot = nBoot, level = level, seed = seed,
        fitArgs = fitArgs,
        residualPool = residualPool
      ),
      error = function(e) {
        message(sprintf("  FAILED tissue=%s fold=%d: %s", tis, k, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(ciFold)) next

    draws_k <- attr(ciFold, "draws") # capture before column assignment below
    ciFold$tissue <- tis
    ciFold$fold <- k
    ciFold$truth <- truth[testIdx]

    key <- paste(tis, k)
    allCI[[key]] <- ciFold
    allDraws[[key]] <- draws_k
  }
}

ciAll <- do.call(rbind, allCI)
rownames(ciAll) <- NULL
drawsAll <- do.call(rbind, allDraws)
attr(ciAll, "draws") <- drawsAll
ampAll <- do.call(rbind, allAmp)

cal <- calibrationCurve(ciAll, ciAll$truth, levels = c(0.5, 0.7, 0.8, 0.9, 0.95))
cat("\n== Calibration: nominal vs empirical coverage (pooled, all tissues) ==\n")
print(cal)

err <- circularError(ciAll$phase, ciAll$truth)
ampVec <- ampAll$amp[match(paste(ciAll$tissue, ciAll$sample), paste(ampAll$tissue, ampAll$sample))]

thr <- deriveAmplitudeThreshold(ampVec, err)
cat("\n== Amplitude threshold (derived empirically) ==\n")
cat(sprintf(
  "threshold=%.3f  low_frac=%.3f  low_mae=%.3fh  high_mae=%.3fh  effect_size=%.3f  p=%.4g\n",
  thr$threshold, thr$low_frac, thr$low_mae, thr$high_mae, thr$effect_size, thr$p_value
))

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)

write.csv(cal, "benchmarks/results/phase3_calibration_GSE54650_timetable.csv", row.names = FALSE)

thrDf <- as.data.frame(thr[c("threshold", "effect_size", "p_value", "low_mae", "high_mae", "low_frac")])
write.csv(thrDf, "benchmarks/results/phase3_amplitude_threshold_GSE54650_timetable.csv", row.names = FALSE)

resid <- data.frame(
  tissue = ciAll$tissue, sample = ciAll$sample, fold = ciAll$fold,
  truth = ciAll$truth, phase = ciAll$phase,
  phase_lo = ciAll$phase_lo, phase_hi = ciAll$phase_hi,
  amp = ampVec, error = err,
  low_amplitude = ampVec <= thr$threshold
)
write.csv(resid, "benchmarks/results/phase3_residuals_GSE54650_timetable.csv", row.names = FALSE)

grDevices::png("benchmarks/results/phase3_calibration_GSE54650_timetable.png", width = 600, height = 600)
plotCalibration(cal)
grDevices::dev.off()

cat("\nWrote benchmarks/results/phase3_calibration_GSE54650_timetable.{csv,png}, ",
    "phase3_amplitude_threshold_GSE54650_timetable.csv, ",
    "phase3_residuals_GSE54650_timetable.csv\n", sep = "")
