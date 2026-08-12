# Reproducible cross-platform transfer benchmark. Not part of the package build
# (see .Rbuildignore). Run from the project root: Rscript benchmarks/run_transfer_benchmark.R
#
# AGENTS.md Section 5 Step 2.2: measures circular MAE when a model trained on
# microarray data (GSE54650) is applied to RNA-seq data from the SAME study/tissues
# (GSE54651, its RNA-seq SubSeries companion - notes/decisions.md Section 10), and
# vice versa, per tissue, for every wrapped engine. Uses transferPhase()
# (R/benchmark.R) rather than a new harness - it reuses the exact same engine
# registry/metrics benchmarkPhase() uses, just with a single fixed train/test split
# (the two platforms) instead of within-dataset CV folds.
#
# Requires benchmarks/data/GSE54650_se.rds and benchmarks/data/GSE54651_se.rds - if
# missing, run inst/scripts/ingest_GSE54650.R and inst/scripts/ingest_GSE54651.R
# first (both download from GEO).
#
# Same caveat as run_benchmark.R: tauFisher's fit step runs real MetaCycle
# rhythmicity detection on tens of thousands of genes, so this script is slow.
# Unlike run_benchmark.R (6-fold CV per tissue), transferPhase() only fits once per
# (tissue, direction, method), so this is roughly 1/6 the cost of the within-dataset
# benchmark for the same engine set - still on the order of tens of minutes to a few
# hours end to end, dominated by tauFisher.

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

microPath <- "benchmarks/data/GSE54650_se.rds"
rnaseqPath <- "benchmarks/data/GSE54651_se.rds"
if (!file.exists(microPath)) {
  stop("Missing ", microPath, " - run inst/scripts/ingest_GSE54650.R first.")
}
if (!file.exists(rnaseqPath)) {
  stop("Missing ", rnaseqPath, " - run inst/scripts/ingest_GSE54651.R first.")
}
seMicro <- readRDS(microPath)
seRnaseq <- readRDS(rnaseqPath)

tissues <- intersect(unique(colData(seMicro)$tissue), unique(colData(seRnaseq)$tissue))
message(sprintf("%d tissues common to both platforms: %s", length(tissues), paste(tissues, collapse = ", ")))

# Same recommended per-engine settings as run_benchmark.R, for apples-to-apples
# comparison against the within-dataset (non-transfer) numbers in
# benchmarks/results/phase2_benchmark_GSE54650.csv.
engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01,
                              nrep = 1, numbasis = 5)),
  zeitzeiger = list(),
  timesignatr = list(fit = list(seed = 1))
)

directions <- list(
  list(name = "microarray_to_rnaseq", trainSe = seMicro, testSe = seRnaseq),
  list(name = "rnaseq_to_microarray", trainSe = seRnaseq, testSe = seMicro)
)

results <- list()
for (dir in directions) {
  for (tis in tissues) {
    trainT <- dir$trainSe[, colData(dir$trainSe)$tissue == tis]
    testT <- dir$testSe[, colData(dir$testSe)$tissue == tis]
    for (m in names(engineArgs)) {
      cat(sprintf("[%s] direction=%s tissue=%s method=%s\n", format(Sys.time()), dir$name, tis, m))
      tb <- tryCatch(
        transferPhase(
          trainT, testT,
          truth_col = "phase", methods = m,
          train_assay_name = "logExpr", test_assay_name = "logExpr",
          engineArgs = engineArgs[m]
        ),
        error = function(e) {
          message(sprintf("  FAILED (%s/%s/%s): %s", dir$name, tis, m, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(tb)) {
        results[[paste(dir$name, tis, m)]] <- cbind(direction = dir$name, tissue = tis, tb$summary)
      }
    }
  }
}

summaryDf <- do.call(rbind, results)
rownames(summaryDf) <- NULL

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(summaryDf, "benchmarks/results/phase2_transfer_benchmark_GSE54650_GSE54651.csv",
          row.names = FALSE)

cat("\n== Phase 2 transfer benchmark: GSE54650 (microarray) <-> GSE54651 (RNA-seq) ==\n")
print(summaryDf, row.names = FALSE)

cat("\nMean circular MAE by direction x method (transfer-error matrix):\n")
print(stats::aggregate(circular_mae ~ direction + method, data = summaryDf, FUN = mean))
cat("Random-guess baseline on a 24h circle: 6h\n")
cat(
  "\nNote: zeitzeiger/timesignatr/taufisher all require their ENTIRE fitted gene",
  "set to be present in the other platform's genes and silently return all-NA",
  "predictions otherwise (only timetable tolerates partial gene-panel overlap -",
  "see notes/decisions.md Section 11) - this shows up here as n == 0 with",
  "n_failed == 0 (nothing errored) and circular_mae == NaN, not as a",
  "benchmarking bug.\n",
  sep = " "
)
