# Phase 3 gate deliverable (AGENTS.md Section 9), second dataset: does the
# clockAmplitude() low-amplitude flag actually predict worse phase error on a
# dataset where damping is plausible (aged tissue, AGENTS.md 7.4), rather than
# GSE54650's healthy-tissue-only atlas (which gave p=0.083, not significant -
# see phase3_amplitude_threshold_GSE54650_timetable.csv)?
#
# See inst/scripts/ingest_GSE84511.R header for why this dataset (young vs aged
# mouse epidermal stem cells, real ZT ground truth) and the honest caveat that its
# source paper (Solanas et al. 2017 Cell) reports core clock GENE expression is
# population-average unaltered by age - this run tests whether clockAmplitude()'s
# per-SAMPLE PCA amplitude nonetheless carries signal, not a guaranteed-positive
# hypothesis.
#
# Single tissue, so no leave-tissue-out residual pooling is needed here (unlike
# run_phase3_calibration.R) - this script only needs point-estimate CV error, since
# the amplitude-threshold gate metric depends on circularError(), not on bootstrap
# CI width. Not part of the package build; see .Rbuildignore.
# Run from the project root: Rscript benchmarks/run_phase3_amplitude_GSE84511.R
# Requires benchmarks/data/GSE84511_se.rds (run inst/scripts/ingest_GSE84511.R first).

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

seFile <- "benchmarks/data/GSE84511_se.rds"
if (!file.exists(seFile)) {
  stop("Run inst/scripts/ingest_GSE84511.R first to produce ", seFile)
}
se <- readRDS(seFile)

nFolds <- 6
seed <- 1
fitArgs <- list(nGenes = 100, minCor = 0.5)

exprMat <- assay(se, "logExpr")
truth <- colData(se)$phase
donors <- colData(se)$donor
ageGroup <- colData(se)$age_group
sampleNames <- colnames(se)

# ---- clock amplitude, computed once across all 48 samples (single tissue group,
# as a real user would - age is not known to clockAmplitude()).
amp <- clockAmplitude(exprMat)
cat(sprintf(
  "clockAmplitude(): %d/%d samples had a finite amplitude (genesUsed=%d: %s)\n",
  sum(!is.na(amp)), length(amp), length(attr(amp, "genesUsed")),
  paste(attr(amp, "genesUsed"), collapse = ", ")
))

# Sanity/interpretability check (not the gate metric itself): does the amplitude
# measure even separate young vs aged on average?
ampByAge <- tapply(amp, ageGroup, mean, na.rm = TRUE)
cat("\n== Mean clockAmplitude by age group (sanity check, not the gate metric) ==\n")
print(ampByAge)
wt <- stats::wilcox.test(amp ~ ageGroup, alternative = "greater") # H1: young > aged
cat(sprintf("Wilcoxon young > aged: W=%.1f p=%.4g\n", wt$statistic, wt$p.value))

# ---- Point-estimate cross-validated phase error (plain CV, no bootstrap needed
# for this gate metric).
cat("\n== Cross-validated point estimates (timetable engine) ==\n")
fold <- .assignDonorFolds(donors, nFolds, seed = seed)
predAll <- rep(NA_real_, length(truth))
for (k in seq_len(nFolds)) {
  trainIdx <- which(fold != k)
  testIdx <- which(fold == k)
  if (length(testIdx) == 0 || length(trainIdx) == 0) next
  fitObj <- fitMolecularTimetable(exprMat[, trainIdx, drop = FALSE], truth[trainIdx],
    nGenes = fitArgs$nGenes, minCor = fitArgs$minCor
  )
  predDf <- predictMolecularTimetable(fitObj, exprMat[, testIdx, drop = FALSE])
  predAll[testIdx] <- predDf$phase
}
stopifnot(!anyNA(predAll))

err <- circularError(predAll, truth)
cat(sprintf("Overall circularMAE: %.3fh\n", mean(err)))
errByAge <- tapply(err, ageGroup, mean)
cat("Mean circularError by age group:\n")
print(errByAge)

# ---- Amplitude threshold (derived empirically) - the actual Phase 3 gate metric.
thr <- deriveAmplitudeThreshold(amp, err)
cat("\n== Amplitude threshold (derived empirically) ==\n")
cat(sprintf(
  "threshold=%.3f  low_frac=%.3f  low_mae=%.3fh  high_mae=%.3fh  effect_size=%.3f  p=%.4g\n",
  thr$threshold, thr$low_frac, thr$low_mae, thr$high_mae, thr$effect_size, thr$p_value
))

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)

thrDf <- as.data.frame(thr[c("threshold", "effect_size", "p_value", "low_mae", "high_mae", "low_frac")])
write.csv(thrDf, "benchmarks/results/phase3_amplitude_threshold_GSE84511_timetable.csv", row.names = FALSE)

resid <- data.frame(
  sample = sampleNames, age_group = ageGroup, zt = colData(se)$zt,
  truth = truth, phase = predAll, amp = as.numeric(amp), error = err,
  low_amplitude = amp <= thr$threshold
)
write.csv(resid, "benchmarks/results/phase3_residuals_GSE84511_timetable.csv", row.names = FALSE)

cat("\nWrote benchmarks/results/phase3_amplitude_threshold_GSE84511_timetable.csv, ",
  "phase3_residuals_GSE84511_timetable.csv\n",
  sep = ""
)
