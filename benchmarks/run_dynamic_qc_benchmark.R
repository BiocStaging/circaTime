# Reproducible dynamic-QC-vs-legacy-QC benchmark. Not part of the package build
# (see .Rbuildignore). Run from the project root:
#   Rscript benchmarks/run_dynamic_qc_benchmark.R
#
# AGENTS.md Section 6 Step 3.2 subtask 3: statistically validate that the new QC
# flag distinguishes low-quality/damped samples from good ones - the specific
# weakness of the legacy fixed 16-gene amplitude flag, which reached only p=0.083
# (GSE54650) and p=0.20 (GSE84511). Subtask 4 wants the two side by side, so both
# scores are computed on exactly the same samples and scored the same way.
#
# Reuses the cached out-of-fold cross-validation predictions rather than re-running
# CV, so the phase estimates here are genuinely held out:
#   benchmarks/results/phase3_residuals_GSE54650_timetable.csv
#   benchmarks/results/phase3_residuals_GSE84511_timetable.csv
# (regenerate them with benchmarks/run_phase3_validation.R if they are missing).
#
# Note on what "phase" means here: clockCoherence() is given the *estimated*
# out-of-fold phase, exactly as estimatePhase() would at inference - not the true
# ZT. The ground truth is used only to score the flag afterwards. See the
# circularity section of ?clockCoherence.

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

resultsDir <- "benchmarks/results"

# AUC via the Mann-Whitney U identity - avoids a dependency just for one number.
aucFromRanks <- function(score, isBad) {
  ok <- !is.na(score) & !is.na(isBad)
  score <- score[ok]; isBad <- isBad[ok]
  if (length(unique(isBad)) < 2) return(NA_real_)
  r <- rank(score)
  nBad <- sum(isBad); nGood <- sum(!isBad)
  # AUC for "low score predicts bad", so flip: 1 - P(score_bad > score_good)
  1 - (sum(r[isBad]) - nBad * (nBad + 1) / 2) / (nBad * nGood)
}

# deriveAmplitudeThreshold() picks the best of ~9 candidate cutoffs by effect size
# and then reports that winner's uncorrected Wilcoxon p. That p is optimistic: the
# cutoff was chosen using the same errors it is tested against. This re-runs the
# whole selection on permuted errors, so the null distribution includes the
# selection step and the resulting p is honest. The legacy flag was originally
# judged by the uncorrected p, so both are reported side by side rather than
# quietly swapping one for the other.
selectionCorrectedP <- function(score, error, nPerm = 2000, seed = 1) {
  obs <- tryCatch(deriveAmplitudeThreshold(score, error), error = function(e) NULL)
  if (is.null(obs)) return(NA_real_)
  set.seed(seed)
  hits <- vapply(seq_len(nPerm), function(b) {
    perm <- tryCatch(deriveAmplitudeThreshold(score, sample(error)),
                     error = function(e) NULL)
    if (is.null(perm)) return(NA)
    perm$effect_size >= obs$effect_size
  }, logical(1))
  (sum(hits, na.rm = TRUE) + 1) / (sum(!is.na(hits)) + 1)
}

scoreOneFlag <- function(label, dataset, score, error, badQuantile = 0.75) {
  ok <- !is.na(score) & !is.na(error)
  score <- score[ok]; error <- error[ok]
  isBad <- error >= stats::quantile(error, badQuantile)

  rho <- suppressWarnings(stats::cor(score, error, method = "spearman"))
  auc <- aucFromRanks(score, isBad)
  thr <- tryCatch(deriveAmplitudeThreshold(score, error), error = function(e) NULL)

  data.frame(
    dataset = dataset, flag = label, n = length(score),
    spearman_rho = round(rho, 3),
    auc_worst_quartile = round(auc, 3),
    threshold = if (is.null(thr)) NA_real_ else round(thr$threshold, 3),
    frac_flagged = if (is.null(thr)) NA_real_ else round(thr$low_frac, 3),
    low_mae = if (is.null(thr)) NA_real_ else round(thr$low_mae, 3),
    high_mae = if (is.null(thr)) NA_real_ else round(thr$high_mae, 3),
    p_value = if (is.null(thr)) NA_real_ else signif(thr$p_value, 3),
    p_selection_corrected = signif(selectionCorrectedP(score, error), 3),
    stringsAsFactors = FALSE
  )
}

rows <- list()

# ---- GSE54650: 12-tissue healthy mouse atlas, per-tissue panels ----------------
resPath <- file.path(resultsDir, "phase3_residuals_GSE54650_timetable.csv")
sePath <- "benchmarks/data/GSE54650_se.rds"
if (file.exists(resPath) && file.exists(sePath)) {
  cat("[GSE54650] loading\n")
  res <- utils::read.csv(resPath, stringsAsFactors = FALSE)
  se <- readRDS(sePath)
  mat <- assay(se, "logExpr")[, res$sample, drop = FALSE]

  cat("[GSE54650] dynamic panel selection (per tissue)\n")
  coh <- clockCoherence(mat, phase = res$phase, tissue = res$tissue)
  cat("  panel sizes:\n")
  print(attr(coh, "panelSize"))
  if (any(attr(coh, "relaxed"), na.rm = TRUE)) {
    cat("  NOTE: relaxed (top-minPanel fallback) in: ",
        paste(names(which(attr(coh, "relaxed"))), collapse = ", "), "\n")
  }

  rows[[length(rows) + 1]] <- scoreOneFlag("legacy_amplitude", "GSE54650", res$amp, res$error)
  rows[[length(rows) + 1]] <- scoreOneFlag("dynamic_coherence", "GSE54650", coh, res$error)

  res$coherence <- coh
  utils::write.csv(res, file.path(resultsDir, "phase3_dynamic_qc_GSE54650_persample.csv"),
                   row.names = FALSE)
} else {
  cat("SKIP GSE54650: missing ", resPath, " or ", sePath, "\n")
}

# ---- GSE84511: young vs aged epidermal stem cells -----------------------------
# The dataset the legacy flag specifically failed on, and the only one here with a
# labelled biological damping contrast, so it is the real test of Step 3.2.
resPath <- file.path(resultsDir, "phase3_residuals_GSE84511_timetable.csv")
sePath <- "benchmarks/data/GSE84511_se.rds"
if (file.exists(resPath) && file.exists(sePath)) {
  cat("[GSE84511] loading\n")
  res <- utils::read.csv(resPath, stringsAsFactors = FALSE)
  se <- readRDS(sePath)
  mat <- assay(se, "logExpr")[, res$sample, drop = FALSE]

  # One panel for the whole cohort: splitting by age group would select a
  # different panel for each arm and make the two scores incomparable, which is
  # the opposite of what a QC flag should do.
  coh <- clockCoherence(mat, phase = res$phase)
  cat("  panel size: ", attr(coh, "panelSize"), " relaxed: ",
      attr(coh, "relaxed"), "\n", sep = "")

  rows[[length(rows) + 1]] <- scoreOneFlag("legacy_amplitude", "GSE84511", res$amp, res$error)
  rows[[length(rows) + 1]] <- scoreOneFlag("dynamic_coherence", "GSE84511", coh, res$error)

  res$coherence <- coh
  utils::write.csv(res, file.path(resultsDir, "phase3_dynamic_qc_GSE84511_persample.csv"),
                   row.names = FALSE)

  # Does either score recover the known biological contrast (aged clocks are the
  # damped ones)? The legacy flag did not: p=0.74 for amplitude by age group.
  cat("\nGSE84511 by age group (aged expected to be more damped):\n")
  for (nm in c("amp", "coherence")) {
    v <- res[[nm]]
    aged <- res$age_group == "aged"
    tt <- suppressWarnings(stats::wilcox.test(v[aged], v[!aged], alternative = "less"))
    cat(sprintf("  %-10s aged mean=%.3f young mean=%.3f  Wilcoxon p=%.3g\n",
                nm, mean(v[aged], na.rm = TRUE), mean(v[!aged], na.rm = TRUE), tt$p.value))
  }
} else {
  cat("SKIP GSE84511: missing ", resPath, " or ", sePath, "\n")
}

summaryDf <- do.call(rbind, rows)
rownames(summaryDf) <- NULL
utils::write.csv(summaryDf, file.path(resultsDir, "phase3_dynamic_qc_vs_legacy.csv"),
                 row.names = FALSE)

cat("\n== Phase 3 Step 3.2: dynamic coherence vs. legacy amplitude ==\n")
print(summaryDf, row.names = FALSE)
cat(
  "\nspearman_rho: correlation of the score with absolute circular error.",
  "\n  Both flags are 'low score = bad', so a NEGATIVE rho is the desired direction.",
  "\nauc_worst_quartile: probability the score ranks a worst-quartile-error sample",
  "\n  below a better one. 0.5 = no discrimination.",
  "\np_value: one-sided Wilcoxon that samples below the empirically derived",
  "\n  threshold have higher error - the same test the legacy flag was judged by.",
  "\np_selection_corrected: the same statistic with the threshold SEARCH inside the",
  "\n  permutation null. This is the number to believe; p_value is optimistic.\n"
)
