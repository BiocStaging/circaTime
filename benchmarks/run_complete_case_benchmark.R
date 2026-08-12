# T6 (revision_and_fix.md): complete-case MAE across all four engines.
#
# run_benchmark.R (which produced benchmarks/results/phase2_benchmark_GSE54650.csv)
# only kept bm$summary and discarded bm$residuals, so there is no saved per-sample
# record of exactly which samples tauFisher's failed folds excluded. Per this
# revision's ground rule 1 ("re-run only if per-sample errors were not saved"),
# this script re-runs the identical benchmark (same engineArgs, nFolds, seed as
# run_benchmark.R) and this time keeps bm$residuals, so the complete-case
# intersection can be computed exactly rather than estimated.
#
# Same runtime profile as run_benchmark.R: order of hours, entirely dominated by
# tauFisher's real MetaCycle rhythmicity detection (the other three engines are
# low single-digit seconds per fold).
#
# Run from the project root: Rscript benchmarks/run_complete_case_benchmark.R

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

# Identical to run_benchmark.R's engineArgs - see that script for the rationale
# behind each setting (each engine's own authors' recommended defaults).
engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01,
                              nrep = 1, numbasis = 5)),
  zeitzeiger = list(),
  timesignatr = list(fit = list(seed = 1))
)

summaryResults <- list()
residualResults <- list()
for (tis in unique(colData(se)$tissue)) {
  seT <- se[, colData(se)$tissue == tis]
  for (m in names(engineArgs)) {
    cat(sprintf("[%s] tissue=%s method=%s\n", format(Sys.time()), tis, m))
    bm <- tryCatch(
      benchmarkPhase(
        seT,
        truth_col = "phase",
        methods = m,
        assay_name = "logExpr",
        nFolds = 6,
        seed = 1,
        engineArgs = engineArgs[m]
      ),
      error = function(e) {
        message(sprintf("  FAILED (%s/%s): %s", tis, m, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(bm)) {
      summaryResults[[paste(tis, m)]] <- cbind(tissue = tis, bm$summary)
      residualResults[[paste(tis, m)]] <- cbind(tissue = tis, bm$residuals)
    }
  }
}

summaryDf <- do.call(rbind, summaryResults)
rownames(summaryDf) <- NULL
residDf <- do.call(rbind, residualResults)
rownames(residDf) <- NULL

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(residDf, "benchmarks/results/revision_persample_residuals_GSE54650.csv",
                 row.names = FALSE)

# Sanity check against the already-committed summary file: all-available MAE here
# should match phase2_benchmark_GSE54650.csv (same seed/engineArgs/nFolds), since
# this is the same benchmark re-run, just with residuals kept this time.
allAvail <- stats::aggregate(circular_error ~ tissue + method, data = residDf[!residDf$fold_failed, ],
                              FUN = function(x) c(mae = mean(x), n = length(x)))
allAvail <- do.call(data.frame, allAvail)
names(allAvail) <- c("tissue", "engine", "mae_all_available", "n_all_available")

# Complete case: for each tissue, the sample IDs every one of the four engines
# scored (fold not failed, prediction not NA).
engines <- names(engineArgs)
okByEngine <- lapply(engines, function(m) {
  sub <- residDf[residDf$method == m & !residDf$fold_failed & !is.na(residDf$pred), ]
  split(sub$sample, sub$tissue)
})
names(okByEngine) <- engines

ccRows <- list()
for (tis in unique(residDf$tissue)) {
  common <- Reduce(intersect, lapply(okByEngine, function(x) x[[tis]]))
  for (m in engines) {
    sub <- residDf[residDf$method == m & residDf$tissue == tis & residDf$sample %in% common, ]
    ccRows[[paste(tis, m)]] <- data.frame(
      tissue = tis, engine = m,
      mae_complete_case = mean(sub$circular_error),
      n_complete_case = nrow(sub)
    )
  }
}
ccDf <- do.call(rbind, ccRows)
rownames(ccDf) <- NULL

out <- merge(allAvail, ccDf, by = c("tissue", "engine"))
out <- out[order(out$tissue, out$engine), ]
utils::write.csv(out, "benchmarks/results/revision_complete_case_MAE.csv", row.names = FALSE)

cat("\n== T6: all-available vs. complete-case MAE ==\n")
print(out, row.names = FALSE)
cat("\nPooled means:\n")
print(stats::aggregate(cbind(mae_all_available, mae_complete_case) ~ engine, data = out, FUN = mean))
cat("\nWrote benchmarks/results/revision_persample_residuals_GSE54650.csv\n")
cat("Wrote benchmarks/results/revision_complete_case_MAE.csv\n")
