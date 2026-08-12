# Reproducible benchmark run. Not part of the package build (see .Rbuildignore).
# Run from the project root: Rscript benchmarks/run_benchmark.R
#
# Phase 2 scope (AGENTS.md Section 9): every wrapped/implemented engine, per-tissue
# on GSE54650 (the benchmark backbone, Section 8), each run with its own authors'
# recommended settings (Phase 2 task 3) - no engine is tuned worse than its paper.
# Requires benchmarks/data/GSE54650_se.rds - if missing, run
# inst/scripts/ingest_GSE54650.R first (it downloads from GEO).
#
# tauFisher's fit step runs real MetaCycle rhythmicity detection on ~22k genes per
# fold (Phase 0 licensing note: this is the actual method, not a stand-in), so this
# script is slow (order of hours for 12 tissues x 6 folds, dominated entirely by
# tauFisher - the other three engines each take low single-digit seconds per fold) -
# that cost is inherent to running the real methods rather than reimplementing them.

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

engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  # tauFisher's own SingleDataset vignette defaults: both MetaCycle methods, no
  # p-value filtering (thres = 1.01), nrep = 1 (GSE54650 has one array per
  # (tissue, CT) cell - no replicate structure to encode), numbasis = 5.
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01,
                              nrep = 1, numbasis = 5)),
  # zeitzeiger's own vignette defaults (nKnots = 3, nTime = 10, sumabsv = 1,
  # orth = TRUE) and its worked example's nSpc = 2 - see R/engines-zeitzeiger.R.
  zeitzeiger = list(),
  # TimeSignatR's own package defaults (a = 0.5, recalib = FALSE) - not the paper's
  # Fig 1 hand-tuned penalty (s = exp(-1.42)) or fixed foldid, which are specific
  # to reproducing that one figure, not a general recommended setting - see
  # R/engines-timesignatr.R.
  timesignatr = list(fit = list(seed = 1))
)

results <- list()
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
    if (!is.null(bm)) results[[paste(tis, m)]] <- cbind(tissue = tis, bm$summary)
  }
}

summaryDf <- do.call(rbind, results)
rownames(summaryDf) <- NULL

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(summaryDf, "benchmarks/results/phase2_benchmark_GSE54650.csv",
          row.names = FALSE)

cat("\n== Phase 2 benchmark: GSE54650, per tissue x method ==\n")
print(summaryDf, row.names = FALSE)
cat("\nMean circular MAE by method:\n")
print(stats::aggregate(circular_mae ~ method, data = summaryDf, FUN = mean))
cat("Random-guess baseline on a 24h circle: 6h\n")
cat(
  "\nNote: circular_cor is unreliable here because GSE54650's ground truth is\n",
  "evenly spaced across the full 24h circle - see the \"Known instability\"\n",
  "section of circularCor()'s documentation (R/circular-stats.R). circular_mae\n",
  "is the metric that matters for this dataset.\n",
  sep = ""
)
