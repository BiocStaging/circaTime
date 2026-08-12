# Benchmark harness. Cross-validates a phase-inference engine against labelled
# ground truth, split by donor/animal (never by sample - AGENTS.md Section 9 Phase 1
# task 5), and reports circular error metrics via R/circular-stats.R.
#
# Engine registry: each entry is a list(fit = function(trainMat, trainTime, ...),
# predict = function(fit, testMat, ...) -> data.frame(sample, phase, confidence),
# supervised = TRUE/FALSE). Unsupervised engines get their held-out predictions
# aligned to truth via circularAlign() before scoring,
# per AGENTS.md 7.2; supervised engines (Phase 1: timetable) do not need this, since
# they are trained directly on absolute time and have no rotation ambiguity.

.circaTimeEngines <- list(
  timetable = list(
    fit = function(trainMat, trainTime, ...) fitMolecularTimetable(trainMat, trainTime, ...),
    predict = function(fit, testMat, ...) predictMolecularTimetable(fit, testMat, ...),
    supervised = TRUE
  ),
  taufisher = list(
    fit = function(trainMat, trainTime, ...) fitTauFisher(trainMat, trainTime, ...),
    predict = function(fit, testMat, ...) predictTauFisher(fit, testMat, ...),
    supervised = TRUE
  ),
  zeitzeiger = list(
    fit = function(trainMat, trainTime, ...) fitZeitZeiger(trainMat, trainTime, ...),
    predict = function(fit, testMat, ...) predictZeitZeiger(fit, testMat, ...),
    supervised = TRUE
  ),
  timesignatr = list(
    fit = function(trainMat, trainTime, ...) fitTimeSignatR(trainMat, trainTime, ...),
    predict = function(fit, testMat, ...) predictTimeSignatR(fit, testMat, ...),
    supervised = TRUE
  )
)

#' Assign donors to cross-validation folds
#'
#' Splits by unique donor, not by sample (AGENTS.md Section 9 Phase 1 task 5), so
#' that repeated measures from the same donor never straddle train/test. If `seed`
#' is supplied, the global RNG state is saved before seeding and restored on exit,
#' so this does not affect the caller's subsequent random draws (AGENTS.md Section
#' 10: no `set.seed()` inside package functions without the caller asking for it).
#'
#' @param donors Character/factor vector, one donor ID per sample.
#' @param nFolds Number of folds.
#' @param seed Optional integer seed for the fold assignment.
#' @return A vector the same length as `donors`, giving each sample's fold (1..nFolds).
#' @keywords internal
.assignDonorFolds <- function(donors, nFolds, seed = NULL) {
  uniqueDonors <- unique(donors)
  if (length(uniqueDonors) < nFolds) {
    stop("Fewer unique donors (", length(uniqueDonors), ") than nFolds (", nFolds, ")")
  }

  drawFolds <- function() sample(rep(seq_len(nFolds), length.out = length(uniqueDonors)))

  if (!is.null(seed)) withr::local_seed(seed)

  foldByDonor <- stats::setNames(drawFolds(), uniqueDonors)
  unname(foldByDonor[as.character(donors)])
}

#' Cross-validated benchmark of phase-inference engines against ground truth
#'
#' @param se A `SummarizedExperiment` with an expression assay (genes x samples)
#'   and a ground-truth phase column in `colData`.
#' @param truth_col Name of the `colData(se)` column holding ground-truth phase, in
#'   `[0, period)`.
#' @param methods Character vector of engine names to benchmark: `"timetable"`
#'   (Phase 1), `"taufisher"`, `"zeitzeiger"`, and `"timesignatr"` (Phase 2).
#' @param assay_name Name of the assay to use; defaults to the first assay.
#' @param donor_col Name of the `colData(se)` column identifying donor/animal, used
#'   for held-out splitting. If `NULL` (default), each sample is treated as its own
#'   donor - only appropriate for datasets with no repeated-measures structure (see
#'   `inst/scripts/ingest_GSE54650.R` for a worked justification).
#' @param nFolds Number of cross-validation folds, split by donor (default 5).
#' @param period Length of the cycle (default 24).
#' @param seed Optional integer seed for reproducible fold assignment; see
#'   [.assignDonorFolds()] for how the global RNG state is protected.
#' @param engineArgs Named list, keyed by method, of extra arguments to pass to
#'   that engine's `fit` and `predict` functions separately (they generally take
#'   different arguments), e.g.
#'   `list(timetable = list(fit = list(nGenes = 100), predict = list(minGenes = 5)))`.
#' @examples
#' set.seed(1)
#' nGenes <- 15; nSamples <- 24
#' time <- rep(seq(0, 23, length.out = 6), 4)
#' theta <- time / 24 * 2 * pi
#' truePhase <- runif(nGenes, 0, 2 * pi)
#' mat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - truePhase[i]) + rnorm(nSamples, sd = 0.3)
#' }))
#' rownames(mat) <- paste0("gene", seq_len(nGenes))
#' colnames(mat) <- paste0("s", seq_len(nSamples))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
#' SummarizedExperiment::colData(se)$ZT <- time
#' bm <- benchmarkPhase(se, truth_col = "ZT", methods = "timetable", nFolds = 4, seed = 1)
#' bm$summary
#' @return A list with two elements: `summary`, a `data.frame` with one row per
#'   method (`method`, `circular_mae`, `circular_cor`, `n`, `n_failed`); and
#'   `residuals`, a `data.frame` with one row per (method, held-out sample)
#'   (`method`, `sample`, `donor`, `fold`, `truth`, `pred`, `circular_error`,
#'   `aligned`, `fold_failed`).
#'
#'   If an engine's `fit`/`predict` errors on a given fold (e.g. because that
#'   fold's training data happens to violate an engine-specific requirement, such
#'   as tauFisher's rhythmicity-detection step needing every clock hour
#'   represented), that fold's test samples get `pred = NA`, `fold_failed = TRUE`,
#'   and a `message()` is emitted - the rest of that method's folds still run and
#'   count normally (`n`/`n_failed` in `summary` report how many held-out samples
#'   got a real prediction vs. how many were lost to a fold failure). This never
#'   fabricates a result for a failed sample (AGENTS.md Section 0.3).
#'
#'   `circular_mae` is the primary accuracy number. For evenly-sampled time-course
#'   ground truth (e.g. GSE54650), `circular_cor` can be numerically unstable and
#'   even strongly negative for accurate predictions - see the "Known instability"
#'   section of [circularCor()] before relying on it as a benchmark gate.
#' @export
benchmarkPhase <- function(se, truth_col, methods = "timetable", assay_name = NULL,
                            donor_col = NULL, nFolds = 5, period = 24, seed = NULL,
                            engineArgs = list()) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment")
  }
  cd <- SummarizedExperiment::colData(se)
  if (!truth_col %in% colnames(cd)) {
    stop("truth_col '", truth_col, "' not found in colData(se)")
  }
  unknown <- setdiff(methods, names(.circaTimeEngines))
  if (length(unknown) > 0) {
    stop("Unknown method(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(.circaTimeEngines), collapse = ", "))
  }

  exprMat <- if (is.null(assay_name)) {
    SummarizedExperiment::assay(se)
  } else {
    SummarizedExperiment::assay(se, assay_name)
  }
  truth <- cd[[truth_col]]
  donors <- if (is.null(donor_col)) colnames(se) else cd[[donor_col]]
  sampleNames <- colnames(se)
  if (is.null(sampleNames)) sampleNames <- as.character(seq_len(ncol(se)))

  fold <- .assignDonorFolds(donors, nFolds, seed = seed)

  residualsList <- list()

  for (m in methods) {
    engine <- .circaTimeEngines[[m]]
    fitArgs <- engineArgs[[m]]$fit
    if (is.null(fitArgs)) fitArgs <- list()
    predictArgs <- engineArgs[[m]]$predict
    if (is.null(predictArgs)) predictArgs <- list()

    for (k in seq_len(nFolds)) {
      trainIdx <- which(fold != k)
      testIdx <- which(fold == k)
      if (length(testIdx) == 0 || length(trainIdx) == 0) next

      testTruth <- truth[testIdx]

      # A single fold's data can violate an engine's own requirements in ways that
      # have nothing to do with its accuracy (e.g. tauFisher's JTK_CYCLE step
      # hard-errors if a donor-based CV split happens to remove every replicate of
      # one clock hour from the training set). Per AGENTS.md 0.3/13: never fabricate
      # a result, but a fold failure is not the whole method failing - isolate it so
      # the other folds' genuine results still count, and report the failure rather
      # than silently dropping it or crashing the whole benchmark.
      foldResult <- tryCatch(
        {
          fitObj <- do.call(engine$fit, c(
            list(trainMat = exprMat[, trainIdx, drop = FALSE], trainTime = truth[trainIdx]),
            fitArgs
          ))
          predDf <- do.call(engine$predict, c(
            list(fit = fitObj, testMat = exprMat[, testIdx, drop = FALSE]),
            predictArgs
          ))
          list(phase = predDf$phase, error = NA_character_)
        },
        error = function(e) {
          message(sprintf(
            "benchmarkPhase: method '%s' fold %d/%d failed (%d test sample(s) set to NA): %s",
            m, k, nFolds, length(testIdx), conditionMessage(e)
          ))
          list(phase = rep(NA_real_, length(testIdx)), error = conditionMessage(e))
        }
      )

      predPhase <- foldResult$phase
      foldFailed <- !is.na(foldResult$error)
      aligned <- FALSE

      if (isFALSE(engine$supervised) && !all(is.na(predPhase))) {
        alignRes <- circularAlign(predPhase, testTruth, period = period)
        predPhase <- alignRes$aligned
        aligned <- TRUE
      }

      residualsList[[length(residualsList) + 1]] <- data.frame(
        method = m,
        sample = sampleNames[testIdx],
        donor = as.character(donors[testIdx]),
        fold = k,
        truth = testTruth,
        pred = predPhase,
        circular_error = circularError(predPhase, testTruth),
        aligned = aligned,
        fold_failed = foldFailed
      )
    }
  }

  residuals <- do.call(rbind, residualsList)
  rownames(residuals) <- NULL

  summaryList <- lapply(methods, function(m) {
    sub <- residuals[residuals$method == m, ]
    data.frame(
      method = m,
      circular_mae = mean(sub$circular_error, na.rm = TRUE),
      circular_cor = circularCor(sub$pred, sub$truth, period = period),
      n = sum(!is.na(sub$pred)),
      n_failed = sum(sub$fold_failed)
    )
  })
  summaryDf <- do.call(rbind, summaryList)

  list(summary = summaryDf, residuals = residuals)
}

#' Cross-dataset transfer benchmark of phase-inference engines
#'
#' Fits each engine on one `SummarizedExperiment` and scores it on a *different*
#' one (AGENTS.md Section 5 Step 2.2), e.g. training on microarray data and
#' predicting on RNA-seq data from the same tissue/study. This is a distinct
#' question from [benchmarkPhase()]'s cross-validation (how well an engine
#' predicts held-out samples from the *same* dataset it was trained on) and
#' [benchmarkPhase()] has no way to ask it - there is exactly one train/test split
#' here, fixed by which dataset a sample belongs to, not `nFolds` random ones. This
#' function is not a parallel benchmarking harness, though: it reuses the same
#' engine registry, gene-name-based fit/predict matching (each engine's own
#' `predict` function already intersects `fit$genes` against `rownames(testMat)`,
#' so `trainSe`/`testSe` do not need to share an identical gene panel), and
#' circular-error metrics that [benchmarkPhase()] uses internally.
#'
#' @param trainSe,testSe `SummarizedExperiment` objects to fit on and score
#'   against, respectively (e.g. a microarray cohort and an RNA-seq cohort of the
#'   same tissue). Their gene panels need not match exactly.
#' @param truth_col Name of the `colData` column holding ground-truth phase, in
#'   `[0, period)`. Must exist in both `trainSe` and `testSe`.
#' @param methods Character vector of engine names; see [benchmarkPhase()].
#' @param train_assay_name,test_assay_name Assay to use in `trainSe`/`testSe`
#'   respectively; each defaults to that object's first assay.
#' @param period Length of the cycle (default 24).
#' @param engineArgs Named list, keyed by method, of extra arguments to pass to
#'   that engine's `fit` and `predict` functions; see [benchmarkPhase()].
#' @return A list with two elements: `summary`, a `data.frame` with one row per
#'   method (`method`, `circular_mae`, `circular_cor`, `n`, `n_failed`); and
#'   `residuals`, a `data.frame` with one row per (method, `testSe` sample)
#'   (`method`, `sample`, `truth`, `pred`, `circular_error`, `aligned`, `failed`).
#'
#'   If an engine's `fit`/`predict` call actually *errors*, that method's
#'   predictions are all `NA` and `n_failed` records how many test samples were
#'   lost, exactly as [benchmarkPhase()] does for a failed fold. Separately, three
#'   of the four wrapped engines ([predictZeitZeiger()], [predictTimeSignatR()],
#'   and [predictTauFisher()]) require the *entire* fitted gene set to be present
#'   in `testMat` and silently return all-`NA` phase (no error at all) when it is
#'   not - only [predictMolecularTimetable()] tolerates partial gene-panel overlap
#'   (see their own docs). That case shows up here as `n = 0` with `n_failed = 0`
#'   (not counted as a "failure" - nothing errored), not `n_failed = n`; check `n`
#'   as well as `n_failed` when auditing a transfer result, and see
#'   `notes/decisions.md` for a worked example against real cross-platform data.
#'   This never fabricates a result for a failed sample (AGENTS.md Section 0.3).
#' @examples
#' set.seed(1)
#' makeSe <- function(n) {
#'   time <- seq(0, 23, length.out = n)
#'   theta <- time / 24 * 2 * pi
#'   mat <- rbind(
#'     g1 = 5 + 2 * cos(theta) + rnorm(n, sd = 0.2),
#'     g2 = 5 + 2 * sin(theta) + rnorm(n, sd = 0.2)
#'   )
#'   colnames(mat) <- paste0("s", seq_len(n))
#'   se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
#'   SummarizedExperiment::colData(se)$ZT <- time
#'   se
#' }
#' tb <- transferPhase(makeSe(12), makeSe(8), truth_col = "ZT", methods = "timetable")
#' tb$summary
#' @export
transferPhase <- function(trainSe, testSe, truth_col, methods = "timetable",
                           train_assay_name = NULL, test_assay_name = NULL,
                           period = 24, engineArgs = list()) {
  if (!methods::is(trainSe, "SummarizedExperiment") ||
      !methods::is(testSe, "SummarizedExperiment")) {
    stop("trainSe and testSe must both be SummarizedExperiment objects")
  }
  trainCd <- SummarizedExperiment::colData(trainSe)
  testCd <- SummarizedExperiment::colData(testSe)
  if (!truth_col %in% colnames(trainCd) || !truth_col %in% colnames(testCd)) {
    stop("truth_col '", truth_col, "' not found in colData of trainSe and/or testSe")
  }
  unknown <- setdiff(methods, names(.circaTimeEngines))
  if (length(unknown) > 0) {
    stop("Unknown method(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(.circaTimeEngines), collapse = ", "))
  }

  trainMat <- if (is.null(train_assay_name)) {
    SummarizedExperiment::assay(trainSe)
  } else {
    SummarizedExperiment::assay(trainSe, train_assay_name)
  }
  testMat <- if (is.null(test_assay_name)) {
    SummarizedExperiment::assay(testSe)
  } else {
    SummarizedExperiment::assay(testSe, test_assay_name)
  }
  trainTime <- trainCd[[truth_col]]
  testTruth <- testCd[[truth_col]]
  sampleNames <- colnames(testSe)
  if (is.null(sampleNames)) sampleNames <- as.character(seq_len(ncol(testSe)))

  residualsList <- list()
  for (m in methods) {
    engine <- .circaTimeEngines[[m]]
    fitArgs <- engineArgs[[m]]$fit
    if (is.null(fitArgs)) fitArgs <- list()
    predictArgs <- engineArgs[[m]]$predict
    if (is.null(predictArgs)) predictArgs <- list()

    fitResult <- tryCatch(
      {
        fitObj <- do.call(engine$fit, c(list(trainMat = trainMat, trainTime = trainTime), fitArgs))
        predDf <- do.call(engine$predict, c(list(fit = fitObj, testMat = testMat), predictArgs))
        list(phase = predDf$phase, error = NA_character_)
      },
      error = function(e) {
        message(sprintf(
          "transferPhase: method '%s' failed (%d test sample(s) set to NA): %s",
          m, ncol(testMat), conditionMessage(e)
        ))
        list(phase = rep(NA_real_, ncol(testMat)), error = conditionMessage(e))
      }
    )

    predPhase <- fitResult$phase
    failed <- !is.na(fitResult$error)
    aligned <- FALSE

    if (isFALSE(engine$supervised) && !all(is.na(predPhase))) {
      alignRes <- circularAlign(predPhase, testTruth, period = period)
      predPhase <- alignRes$aligned
      aligned <- TRUE
    }

    residualsList[[length(residualsList) + 1]] <- data.frame(
      method = m,
      sample = sampleNames,
      truth = testTruth,
      pred = predPhase,
      circular_error = circularError(predPhase, testTruth),
      aligned = aligned,
      failed = failed
    )
  }

  residuals <- do.call(rbind, residualsList)
  rownames(residuals) <- NULL

  summaryList <- lapply(methods, function(m) {
    sub <- residuals[residuals$method == m, ]
    data.frame(
      method = m,
      circular_mae = mean(sub$circular_error, na.rm = TRUE),
      circular_cor = circularCor(sub$pred, sub$truth, period = period),
      n = sum(!is.na(sub$pred)),
      n_failed = if (isTRUE(sub$failed[1])) nrow(sub) else 0L
    )
  })
  summaryDf <- do.call(rbind, summaryList)

  list(summary = summaryDf, residuals = residuals)
}
