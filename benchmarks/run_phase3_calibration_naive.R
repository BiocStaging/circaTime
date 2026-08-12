# Naive-bootstrap counterpart to run_phase3_calibration.R: same design (6-fold
# donor CV per tissue, nBoot=50, level=0.8), but bootstrapPhaseCI() called with
# NO residualPool - training-resampling variance only. run_phase3_calibration.R's
# header comment already reports one measured point from an earlier version of
# this run (~48% empirical coverage at 80% nominal); this script regenerates the
# full 5-level curve for that variant so it isn't a single anecdote anymore.
# Run from the project root: Rscript benchmarks/run_phase3_calibration_naive.R

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

cat("== Naive bootstrap (no residualPool) confidence intervals ==\n")
allCI <- list()
allDraws <- list()

for (tis in tissues) {
  cat(sprintf("[%s] tissue=%s\n", format(Sys.time()), tis))
  seT <- se[, colData(se)$tissue == tis]
  exprMat <- assay(seT, "logExpr")
  truth <- colData(seT)$phase
  donors <- colData(seT)$donor

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
        fitArgs = fitArgs
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

cal <- calibrationCurve(ciAll, ciAll$truth, levels = c(0.5, 0.7, 0.8, 0.9, 0.95))
cat("\n== Calibration: nominal vs empirical coverage (pooled, all tissues, NO residual pool) ==\n")
print(cal)

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(cal, "benchmarks/results/phase3_calibration_GSE54650_timetable_naive.csv", row.names = FALSE)
cat("\nWrote benchmarks/results/phase3_calibration_GSE54650_timetable_naive.csv\n")
