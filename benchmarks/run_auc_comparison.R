# T10 (revision_and_fix.md): the fixed-panel vs. dynamic-coherence AUCs
# (Table~\ref{tab:dynamicqc}, GSE54650 and GSE84511) are two AUCs on the SAME
# samples and were never compared to each other directly - one clearing p<0.05
# and the other not is not evidence the two differ. This computes the paired
# comparison: DeLong's test (pROC::roc.test) plus a paired bootstrap over
# samples for a CI on the AUC difference.
#
# Reuses the already-committed per-sample files from run_dynamic_qc_benchmark.R
# (amp, coherence, error columns for the same held-out samples) rather than
# rerunning CV, per ground rule 1.
# Run from the project root: Rscript benchmarks/run_auc_comparison.R

suppressMessages(library(pROC))

set.seed(1)
nBoot <- 10000

pairedAucComparison <- function(dataset, df, badQuantile = 0.75) {
  ok <- !is.na(df$amp) & !is.na(df$coherence) & !is.na(df$error)
  df <- df[ok, ]
  isBad <- df$error >= stats::quantile(df$error, badQuantile)

  # Both flags are "low score = bad": controls (good) score higher than cases
  # (bad), so direction ">" (matches the aucFromRanks() convention already used
  # in run_dynamic_qc_benchmark.R / Table~\ref{tab:dynamicqc}).
  rocFixed   <- pROC::roc(response = isBad, predictor = df$amp,       direction = ">", quiet = TRUE)
  rocDynamic <- pROC::roc(response = isBad, predictor = df$coherence, direction = ">", quiet = TRUE)

  dl <- pROC::roc.test(rocFixed, rocDynamic, method = "delong", paired = TRUE)

  aucFixed <- as.numeric(pROC::auc(rocFixed))
  aucDynamic <- as.numeric(pROC::auc(rocDynamic))
  aucDiffObs <- aucDynamic - aucFixed

  n <- nrow(df)
  boot <- vapply(seq_len(nBoot), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    isBadB <- isBad[idx]
    if (length(unique(isBadB)) < 2) return(NA_real_)
    aFixed <- tryCatch(as.numeric(pROC::auc(pROC::roc(
      response = isBadB, predictor = df$amp[idx], direction = ">", quiet = TRUE))),
      error = function(e) NA_real_)
    aDynamic <- tryCatch(as.numeric(pROC::auc(pROC::roc(
      response = isBadB, predictor = df$coherence[idx], direction = ">", quiet = TRUE))),
      error = function(e) NA_real_)
    aDynamic - aFixed
  }, numeric(1))
  boot <- boot[!is.na(boot)]

  ciLo <- stats::quantile(boot, 0.025, names = FALSE)
  ciHi <- stats::quantile(boot, 0.975, names = FALSE)
  # Two-sided bootstrap p-value: proportion of resamples on the other side of 0
  # from the observed difference, doubled (percentile-based, symmetric convention).
  bootP <- min(1, 2 * min(mean(boot <= 0), mean(boot >= 0)))

  data.frame(
    dataset = dataset, auc_fixed = round(aucFixed, 4), auc_dynamic = round(aucDynamic, 4),
    auc_difference = round(aucDiffObs, 4), ci_lo = round(ciLo, 4), ci_hi = round(ciHi, 4),
    delong_p = signif(dl$p.value, 3), bootstrap_p = signif(bootP, 3), n_samples = n
  )
}

d54650 <- utils::read.csv("benchmarks/results/phase3_dynamic_qc_GSE54650_persample.csv")
d84511 <- utils::read.csv("benchmarks/results/phase3_dynamic_qc_GSE84511_persample.csv")

out <- rbind(
  pairedAucComparison("GSE54650", d54650),
  pairedAucComparison("GSE84511", d84511)
)

utils::write.csv(out, "benchmarks/results/revision_auc_comparison.csv", row.names = FALSE)
print(out, row.names = FALSE)
cat("\nWrote benchmarks/results/revision_auc_comparison.csv\n")
