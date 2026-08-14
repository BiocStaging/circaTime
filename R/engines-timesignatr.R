# TimeSignatR engine (Braun et al. 2018, PNAS: 10.1073/pnas.1800314115). Wraps the
# authors' own R package (github.com/braunr/TimeSignatR, GPL-2) rather than
# reimplementing it, per AGENTS.md Section 4/9 Phase 2. Method: elastic-net regression
# (glmnet, family = "mgaussian") of each training sample's expression onto the
# (sin, cos) unit-circle coordinates of its known clock time, via the authors' own
# `trainTimeStamp()`/`predTimeStamp()`.
#
# Two adaptations at the wrapper boundary, not inside the wrapped algorithm itself:
#   1. `trainTimeStamp()` always internally holds out `1 - trainFrac` of the training
#      subjects for its own diagnostic prediction (`trainFrac` defaults to 0.5 in the
#      authors' own API). circaTime's `benchmarkPhase()` already does its own held-out
#      cross-validation, so this wrapper always calls it with `trainFrac = 1` - using
#      every training sample for the fit rather than silently discarding half of an
#      already donor-limited training fold to an internal holdout it never uses.
#   2. `trainTimeStamp()`/`XY2dectime()` hardcode a 24-hour period. This wrapper
#      rescales `trainTime`/predictions to/from hours at the boundary so the public
#      API still honours `period`, matching the other engines.

#' Fit a TimeSignatR model from labelled training data
#'
#' Calls the authors' own `TimeSignatR::trainTimeStamp()` rather than reimplementing
#' it. Unlike [fitMolecularTimetable()]/[fitTauFisher()], TimeSignatR does no gene
#' selection of its own - the fitted `cv.glmnet` object's coefficient matrix spans
#' every gene in `trainMat`, with the elastic-net penalty (not a pre-filtering step)
#' doing the shrinkage - so [predictTimeSignatR()] requires every one of those genes
#' to be present in `testMat`, same all-or-nothing rule as [predictTauFisher()].
#'
#' @param trainMat Numeric matrix of expression, genes x samples.
#' @param trainTime Numeric vector of known circadian phase in `[0, period)`, one
#'   value per column of `trainMat`.
#' @param subjIDs Character/factor vector of subject/donor IDs, one per column of
#'   `trainMat`, passed to `trainTimeStamp()` (used for its `recalib` option and
#'   its internal subject sampling, which is a no-op here since `trainFrac = 1`
#'   selects every subject regardless). Default: `colnames(trainMat)`, i.e. each
#'   sample treated as its own subject, if `trainMat` has column names; otherwise
#'   `seq_len(ncol(trainMat))`.
#' @param period Length of the cycle (default 24, for hours). `trainTimeStamp()`
#'   itself is hardcoded to a 24-hour period; this wrapper rescales `trainTime` to
#'   hours before calling it (see file header).
#' @param alpha Elastic-net mixing parameter, passed to `trainTimeStamp()`'s `a`
#'   argument (`1` = lasso, `0` = ridge). Default 0.5, the package's own default.
#' @param recalib Logical, passed to `trainTimeStamp()`'s `recalib` argument:
#'   whether to per-subject mean-centre expression (via
#'   `TimeSignatR::recalibrateExprs()`) before fitting. Default `FALSE`, the
#'   package's own default; only meaningful with repeated-measures `subjIDs`.
#' @param seed Optional integer seed. `trainTimeStamp()` draws its (here, no-op)
#'   subject sample and `cv.glmnet()`'s cross-validation fold assignment from the
#'   global RNG; if `seed` is supplied, the global RNG state is saved and restored
#'   (AGENTS.md Section 10), so this call does not affect the caller's subsequent
#'   random draws.
#' @param ... Additional arguments passed to `trainTimeStamp()` (and on to
#'   `glmnet::cv.glmnet()`).
#' @return An object of class `"timeSignatRFit"`: a list with the training gene
#'   names (`genes`), the raw `trainTimeStamp()` result (`timestamp`, containing the
#'   fitted `cv.glmnet` object `predTimeStamp()` needs), and `period`.
#' @examples
#' if (requireNamespace("TimeSignatR", quietly = TRUE)) {
#'   set.seed(1)
#'   nGenes <- 40; time <- seq(0, 23, by = 1)
#'   theta <- time / 24 * 2 * pi
#'   trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'     5 + 2 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(length(time), sd = 0.3)
#'   }))
#'   rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#'   colnames(trainMat) <- paste0("s", seq_len(length(time)))
#'   fit <- fitTimeSignatR(trainMat, time, seed = 1)
#'   length(fit$genes)
#' }
#' @export
fitTimeSignatR <- function(trainMat, trainTime, subjIDs = NULL, period = 24,
                            alpha = 0.5, recalib = FALSE, seed = NULL, ...) {
  .requireEngine("TimeSignatR", 'remotes::install_github("braunr/TimeSignatR")')
  stopifnot(is.matrix(trainMat), ncol(trainMat) == length(trainTime))
  if (is.null(rownames(trainMat))) {
    stop("trainMat must have rownames (gene identifiers)")
  }
  if (is.null(subjIDs)) {
    subjIDs <- colnames(trainMat)
    if (is.null(subjIDs)) subjIDs <- seq_len(ncol(trainMat))
  }
  stopifnot(length(subjIDs) == ncol(trainMat))

  timeHours <- (trainTime %% period) * (24 / period)

  if (!is.null(seed)) withr::local_seed(seed) else withr::local_preserve_seed()

  timestamp <- TimeSignatR::trainTimeStamp(
    expr = trainMat, subjIDs = subjIDs, times = timeHours,
    trainFrac = 1, a = alpha, recalib = recalib, ...
  )

  structure(
    list(genes = rownames(trainMat), timestamp = timestamp, period = period),
    class = "timeSignatRFit"
  )
}

#' Predict circadian phase with a fitted TimeSignatR model
#'
#' Calls the authors' own `TimeSignatR::predTimeStamp()`. `cv.glmnet`'s prediction
#' matches features by column position, not by name, so - like [predictZeitZeiger()]
#' - this requires every one of `fit$genes` to be present in `testMat`'s rownames;
#' if any is missing, every sample gets `NA` rather than a prediction computed on a
#' reordered or shrunken gene set.
#'
#' @param fit A `"timeSignatRFit"` object from [fitTimeSignatR()].
#' @param testMat Numeric matrix of expression, genes x samples, same scale as the
#'   training matrix.
#' @return A `data.frame` with one row per column of `testMat`: `sample`, `phase`
#'   (numeric, `[0, period)`, or `NA` if `testMat` was missing any fitted gene), and
#'   `confidence` (`NA_real_` for every row - `predTimeStamp()` returns a point
#'   estimate with no native per-sample confidence signal; calibrated uncertainty is
#'   Phase 3 scope, AGENTS.md Section 9).
#' @examples
#' if (requireNamespace("TimeSignatR", quietly = TRUE)) {
#'   set.seed(1)
#'   nGenes <- 40; time <- seq(0, 23, by = 1)
#'   theta <- time / 24 * 2 * pi
#'   trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'     5 + 2 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(length(time), sd = 0.3)
#'   }))
#'   rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#'   colnames(trainMat) <- paste0("s", seq_len(length(time)))
#'   fit <- fitTimeSignatR(trainMat, time, seed = 1)
#'   pred <- predictTimeSignatR(fit, trainMat)
#'   head(pred)
#' }
#' @export
predictTimeSignatR <- function(fit, testMat) {
  .requireEngine("TimeSignatR", 'remotes::install_github("braunr/TimeSignatR")')
  stopifnot(inherits(fit, "timeSignatRFit"))
  if (is.null(rownames(testMat))) {
    stop("testMat must have rownames (gene identifiers)")
  }

  sampleNames <- colnames(testMat)
  nSamples <- ncol(testMat)
  if (is.null(sampleNames)) sampleNames <- seq_len(nSamples)

  outPhase <- rep(NA_real_, nSamples)

  if (all(fit$genes %in% rownames(testMat))) {
    testSub <- testMat[fit$genes, , drop = FALSE]
    predHours <- TimeSignatR::predTimeStamp(fit$timestamp, newx = testSub)
    outPhase <- as.numeric(predHours) * (fit$period / 24)
  }

  data.frame(
    sample = sampleNames,
    phase = outPhase,
    confidence = NA_real_,
    row.names = NULL
  )
}
