# T12 (revision_and_fix.md): proper Results-section treatment for the human
# (GSE56931) and pseudobulk (GSE145197) cohorts, replacing the "unvalidated
# preview" parentheticals in Section 4 with real, properly-scoped analyses.
#
# GSE56931 (human whole blood): restricted to the 130 unconfounded
# normal-routine baseline samples (protocol_phase == "baseline"; the sleep-
# deprivation and recovery arms are a genuine circadian confound, see
# inst/scripts/ingest_GSE56931.R), donor-split CV over the 14 subjects (this is
# the one cohort in this package where donor != sample, so donor-splitting is
# load-bearing - see Section 2.3/T4). All four engines scored, same style as
# run_complete_case_benchmark.R.
#
# GSE145197 (mouse liver pseudobulk, single-cell aggregated): only 10 samples
# across 4 timepoints - not enough for a defensible CV split, so this is NOT a
# cross-validation benchmark. Instead it evaluates the bundled Molecular
# Timetable RNA-seq liver reference model (trained on GSE54651, see
# inst/scripts/train_reference_models.R) out-of-the-box via estimatePhase(),
# exactly as a real user would call it - zero additional training on this
# script's part.
# Run from the project root: Rscript benchmarks/run_human_pseudobulk_benchmark.R

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

resultsDir <- "benchmarks/results"
dir.create(resultsDir, showWarnings = FALSE, recursive = TRUE)

bootCI <- function(x, nBoot = 10000, seed = 1) {
  set.seed(seed)
  n <- length(x)
  boots <- vapply(seq_len(nBoot), function(i) mean(sample(x, n, replace = TRUE)), numeric(1))
  stats::quantile(boots, c(0.025, 0.975), names = FALSE)
}

## ---- GSE56931: human whole blood, baseline-only, donor-split CV --------------
sePath <- "benchmarks/data/GSE56931_se.rds"
se <- readRDS(sePath)
seBaseline <- se[, colData(se)$protocol_phase == "baseline"]
cat(sprintf("[GSE56931] %d baseline samples, %d donors\n",
            ncol(seBaseline), length(unique(colData(seBaseline)$donor))))
stopifnot(ncol(seBaseline) == 130, length(unique(colData(seBaseline)$donor)) == 14)

engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01,
                              nrep = 1, numbasis = 5)),
  zeitzeiger = list(),
  timesignatr = list(fit = list(seed = 1))
)

humanRows <- list()
for (m in names(engineArgs)) {
  cat(sprintf("[%s] GSE56931 method=%s\n", format(Sys.time()), m))
  bm <- tryCatch(
    benchmarkPhase(
      seBaseline, truth_col = "phase", methods = m, assay_name = "logExpr",
      donor_col = "donor", nFolds = 5, seed = 1, engineArgs = engineArgs[m]
    ),
    error = function(e) {
      message(sprintf("  FAILED (%s): %s", m, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(bm)) humanRows[[m]] <- cbind(method = m, bm$residuals)
}
humanResid <- do.call(rbind, humanRows)
rownames(humanResid) <- NULL

# Save per-sample residuals for significance testing (R6)
utils::write.csv(humanResid, file.path(resultsDir, "revision_human_GSE56931_persample.csv"),
                 row.names = FALSE)

# Summary: per-method pooled MAE with subject-level bootstrap CI (14 subjects)
humanSummary <- do.call(rbind, lapply(split(humanResid, humanResid$method), function(d) {
  ok <- d[!d$fold_failed & !is.na(d$pred), ]
  # Aggregate to subject level (mean error per subject)
  subjMeans <- tapply(ok$circular_error, ok$donor, mean)
  # Subject-level bootstrap CI (resample 14 subjects with replacement, 10k draws)
  ci <- if (length(subjMeans) > 1) {
    set.seed(1)
    boots <- vapply(seq_len(10000), function(i) mean(sample(subjMeans, length(subjMeans), replace = TRUE)), numeric(1))
    stats::quantile(boots, c(0.025, 0.975), names = FALSE)
  } else c(NA, NA)
  data.frame(
    dataset = "GSE56931", method = d$method[1], n = nrow(ok),
    circular_mae = mean(ok$circular_error), ci_lo = ci[1], ci_hi = ci[2]
  )
}))
utils::write.csv(humanSummary, file.path(resultsDir, "revision_human_GSE56931_benchmark.csv"),
                 row.names = FALSE)
cat("\n== GSE56931 (human whole blood, baseline-only, donor-split CV) ==\n")
print(humanSummary, row.names = FALSE)

## ---- GSE145197: mouse liver pseudobulk, out-of-the-box transfer test ---------
sePath2 <- "benchmarks/data/GSE145197_se.rds"
se2 <- readRDS(sePath2)
cat(sprintf("\n[GSE145197] %d samples\n", ncol(se2)))
stopifnot(ncol(se2) == 10)

# estimatePhase() overwrites colData(se)$phase with its OWN prediction (see
# R/estimatePhase.R), so the true ZT must be saved before calling it.
trueZt <- colData(se2)$phase
se2Scored <- estimatePhase(se2, method = "timetable", tissue = "Liv",
                           species = "Mmusculus", platform = "rnaseq",
                           assay_name = "logExpr")
pseudoResid <- data.frame(
  sample = colnames(se2Scored),
  truth = trueZt,
  pred = colData(se2Scored)$phase,
  flag = as.character(colData(se2Scored)$phase_flag)
)
pseudoResid$circular_error <- circularError(pseudoResid$pred, pseudoResid$truth)

pseudoSummary <- data.frame(
  dataset = "GSE145197", method = "timetable (bundled RNA-seq liver model, no training)",
  n = nrow(pseudoResid), circular_mae = mean(pseudoResid$circular_error, na.rm = TRUE),
  min_error = min(pseudoResid$circular_error, na.rm = TRUE),
  max_error = max(pseudoResid$circular_error, na.rm = TRUE)
)
utils::write.csv(pseudoSummary, file.path(resultsDir, "revision_pseudobulk_GSE145197_transfer.csv"),
                 row.names = FALSE)
utils::write.csv(pseudoResid, file.path(resultsDir, "revision_pseudobulk_GSE145197_persample.csv"),
                 row.names = FALSE)
cat("\n== GSE145197 (mouse liver pseudobulk, bundled model, no CV, n=10) ==\n")
print(pseudoSummary, row.names = FALSE)
print(pseudoResid, row.names = FALSE)

cat("\nWrote benchmarks/results/revision_human_GSE56931_benchmark.csv\n")
cat("Wrote benchmarks/results/revision_pseudobulk_GSE145197_transfer.csv\n")
cat("Wrote benchmarks/results/revision_pseudobulk_GSE145197_persample.csv\n")
