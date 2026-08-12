# Reproducible Bayesian-vs-bootstrap CI benchmark. Not part of the package build
# (see .Rbuildignore). Run from the project root:
#   Rscript benchmarks/run_bayesian_ci_benchmark.R
#
# AGENTS.md Section 6 Step 3.3 / Section 7 DoD item 5: "each new algorithm has a
# benchmark comparison against the pre-existing baseline it replaces or augments
# ... Bayesian CI vs bootstrap CI". test-confidence.R already unit-tests
# consistency on synthetic data; this is the real-data companion, using exactly
# the bundled reference models' own residual pools (PhaseReferenceModel@residualPool)
# and the real held-out cross-validated predictions already on disk from
# benchmarks/run_phase3_validation.R, so both intervals are compared on the
# identical genuine out-of-fold predictions rather than a fresh simulation.
#
# Requires benchmarks/results/phase3_residuals_GSE54650_timetable.csv (regenerate
# with benchmarks/run_phase3_validation.R if missing) and the bundled reference
# models (R/sysdata.rda, built by inst/scripts/train_reference_models.R).

suppressMessages(library(pkgload))
pkgload::load_all(".", quiet = TRUE)

resPath <- "benchmarks/results/phase3_residuals_GSE54650_timetable.csv"
if (!file.exists(resPath)) {
  stop("Missing ", resPath, " - run benchmarks/run_phase3_validation.R first.")
}
res <- utils::read.csv(resPath, stringsAsFactors = FALSE)

level <- 0.8
rows <- list()
densityRows <- list()

for (tis in unique(res$tissue)) {
  refModel <- circaTime:::.referenceModels[[paste("Mmusculus", tis, "timetable", "microarray", sep = ".")]]
  if (is.null(refModel)) {
    cat("skipping", tis, ": no bundled timetable model\n")
    next
  }
  sub <- res[res$tissue == tis, ]

  post <- bayesianPhasePosterior(point = sub$phase, residualPool = refModel@residualPool, level = level)

  bayesCovers <- circularIntervalCovers(sub$truth, post$credible_lo, post$credible_hi)
  bootCovers <- circularIntervalCovers(sub$truth, sub$phase_lo, sub$phase_hi)
  bayesWidth <- circularDiff(post$credible_hi, post$credible_lo)
  bootWidth <- circularDiff(sub$phase_hi, sub$phase_lo)
  # both intervals are typically well under half the circle, but guard against the
  # circularDiff sign-fold the same way the unit tests had to (see
  # test-circular-stats.R): a "negative width" means the true width exceeded 12h.
  bayesWidth <- ifelse(bayesWidth < 0, bayesWidth + 24, bayesWidth)
  bootWidth <- ifelse(bootWidth < 0, bootWidth + 24, bootWidth)

  rows[[length(rows) + 1]] <- data.frame(
    tissue = tis, n = nrow(sub), kappa = post$kappa, concentration = post$concentration,
    bayes_coverage = mean(bayesCovers), boot_coverage = mean(bootCovers),
    bayes_mean_width = mean(bayesWidth), boot_mean_width = mean(bootWidth),
    mean_width_diff_h = mean(bayesWidth - bootWidth)
  )
}

summaryDf <- do.call(rbind, rows)
rownames(summaryDf) <- NULL

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(summaryDf, "benchmarks/results/phase3_bayesian_vs_bootstrap_ci_GSE54650.csv", row.names = FALSE)

cat("\n== Phase 3 Step 3.3: Bayesian credible interval vs. bootstrap CI, GSE54650 ==\n")
print(summaryDf, row.names = FALSE)

cat(sprintf(
  "\nNominal level: %.0f%%. Pooled empirical coverage: Bayesian %.1f%%, bootstrap %.1f%%.\n",
  level * 100, 100 * mean(rep(summaryDf$bayes_coverage, summaryDf$n)),
  100 * mean(rep(summaryDf$boot_coverage, summaryDf$n))
))
cat(sprintf(
  "Pooled mean width: Bayesian %.2fh, bootstrap %.2fh (difference %.2fh).\n",
  weighted.mean(summaryDf$bayes_mean_width, summaryDf$n),
  weighted.mean(summaryDf$boot_mean_width, summaryDf$n),
  weighted.mean(summaryDf$mean_width_diff_h, summaryDf$n)
))
cat(
  "\nInterpretation: the two intervals are built from the SAME residual pool by",
  "\nconstruction, so close agreement here is expected, not a coincidence - this",
  "\nchecks that the parametric (von Mises) approximation does not distort the",
  "\nempirically-validated bootstrap interval's calibration, not that the two are",
  "\nindependent evidence of anything. See ?bayesianPhasePosterior's Relationship",
  "\nsection for when the two are expected to diverge (a non-von-Mises-shaped pool).\n"
)
