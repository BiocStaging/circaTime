# ZeitZeiger engine (Hughey, Hastie & Butte 2016, Nucleic Acids Res: 10.1093/nar/gkw030).
# Wraps the authors' own R package (github.com/hugheylab/zeitzeiger, GPL-2) rather than
# reimplementing it, per AGENTS.md Section 4/9 Phase 2. Method: fit a periodic smoothing
# spline per gene against a labelled training time course (zeitzeigerFit()), reduce the
# time-dependent variation to a handful of sparse principal components (zeitzeigerSpc()),
# then predict a new sample's phase by maximum-likelihood optimisation over those
# components (zeitzeigerPredict()).
#
# zeitzeiger's own functions use observations-in-rows/features-in-columns matrices and
# a time scale normalised to [0, 1); this wrapper transposes to/from circaTime's
# genes-x-samples convention and circaTime's [0, period) time scale at the boundary only.

#' Fit a ZeitZeiger model from labelled training data
#'
#' Calls the authors' own pipeline (`zeitzeiger::zeitzeigerFit()` +
#' `zeitzeiger::zeitzeigerSpc()`) rather than reimplementing it. Unlike
#' [fitTauFisher()], zeitzeiger has no separate "trained model" object - the fitted
#' periodic splines are re-derived from `xTrain`/`timeTrain` at prediction time
#' (`zeitzeiger::zeitzeigerPredict()`'s own design) - so the training matrix and
#' times are carried in the returned object as well as the SPC result.
#'
#' @param trainMat Numeric matrix of expression, genes x samples.
#' @param trainTime Numeric vector of known circadian phase in `[0, period)`, one
#'   value per column of `trainMat`.
#' @param period Length of the cycle (default 24, for hours).
#' @param nKnots Number of internal knots for the periodic smoothing spline, passed
#'   to `zeitzeigerFit()`. Default 3, the package's own default.
#' @param nTime Number of time-points used to discretise each gene's fitted curve
#'   before computing sparse principal components, passed to `zeitzeigerSpc()`.
#'   Default 10, the package's own default.
#' @param sumabsv L1-constraint on the sparse principal components, passed to
#'   `zeitzeigerSpc()`. Default 1, the package's own default.
#' @param orth Logical; require the components be orthogonal, passed to
#'   `zeitzeigerSpc()`. Default `TRUE`, the package's own default.
#' @param useSpc Logical; use `PMA::SPC()` for the sparse components (default) or
#'   fall back to `base::svd()`, passed to `zeitzeigerSpc()`.
#' @return An object of class `"zeitZeigerFit"`: a list with the training gene names
#'   (`genes`), the transposed training matrix (`xTrain`) and fractional training
#'   time (`timeTrain`) that `zeitzeigerPredict()` needs at prediction time, the SPC
#'   result (`spcResult`), `nKnots`, and `period`.
#' @examples
#' set.seed(1)
#' nGenes <- 40; time <- seq(0, 23, by = 1)
#' theta <- time / 24 * 2 * pi
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(length(time), sd = 0.3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' fit <- fitZeitZeiger(trainMat, time)
#' length(fit$genes)
#' @export
fitZeitZeiger <- function(trainMat, trainTime, period = 24, nKnots = 3, nTime = 10,
                           sumabsv = 1, orth = TRUE, useSpc = TRUE) {
  .requireEngine("zeitzeiger", 'remotes::install_github("hugheylab/zeitzeiger")')
  stopifnot(is.matrix(trainMat), ncol(trainMat) == length(trainTime))
  if (is.null(rownames(trainMat))) {
    stop("trainMat must have rownames (gene identifiers)")
  }

  xTrain <- t(trainMat)
  timeTrain <- (trainTime %% period) / period

  fitRes <- zeitzeiger::zeitzeigerFit(xTrain, timeTrain, nKnots = nKnots)
  spcResult <- zeitzeiger::zeitzeigerSpc(
    fitRes$xFitMean, fitRes$xFitResid,
    nTime = nTime, useSpc = useSpc, sumabsv = sumabsv, orth = orth
  )

  structure(
    list(
      genes = rownames(trainMat),
      xTrain = xTrain,
      timeTrain = timeTrain,
      spcResult = spcResult,
      nKnots = nKnots,
      period = period
    ),
    class = "zeitZeigerFit"
  )
}

#' Predict circadian phase with a fitted ZeitZeiger model
#'
#' Calls the authors' own `zeitzeiger::zeitzeigerPredict()`. zeitzeiger matches
#' training and test features by column position, not by name, so - unlike
#' [predictMolecularTimetable()], which silently uses whatever fitted genes overlap
#' `testMat` - this requires every one of `fit$genes` to be present in `testMat`'s
#' rownames; if any is missing, every sample gets `NA` rather than a prediction
#' computed on a reordered or shrunken gene set (same all-or-nothing rule as
#' [predictTauFisher()]).
#'
#' @param fit A `"zeitZeigerFit"` object from [fitZeitZeiger()].
#' @param testMat Numeric matrix of expression, genes x samples, same scale as the
#'   training matrix.
#' @param nSpc Number of sparse principal components to use for prediction, passed
#'   to `zeitzeigerPredict()`'s `nSpc` argument as a single value (a single
#'   prediction per sample, not zeitzeiger's default sweep over `1:K`). Default 2.
#'   Clamped down to the number of components `zeitzeigerSpc()` actually returned
#'   (`length(fit$spcResult$d)`) if `nSpc` exceeds it, e.g. for small gene sets.
#' @param timeRange Vector of fractional-time grid points searched for the
#'   maximum-likelihood initial value, passed through to `zeitzeigerPredict()`.
#'   Default matches the package's own default resolution (100 points over one
#'   cycle).
#' @return A `data.frame` with one row per column of `testMat`: `sample`, `phase`
#'   (numeric, `[0, period)`, or `NA` if `testMat` was missing any fitted gene), and
#'   `confidence` (`NA_real_` for every row - zeitzeiger's per-sample likelihood
#'   surface is not converted to a calibrated confidence score here; calibrated
#'   uncertainty is Phase 3 scope, AGENTS.md Section 9).
#' @examples
#' set.seed(1)
#' nGenes <- 40; time <- seq(0, 23, by = 1)
#' theta <- time / 24 * 2 * pi
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(length(time), sd = 0.3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' fit <- fitZeitZeiger(trainMat, time)
#' pred <- predictZeitZeiger(fit, trainMat)
#' head(pred)
#' @export
predictZeitZeiger <- function(fit, testMat, nSpc = 2,
                               timeRange = seq(0, 1 - 0.01, 0.01)) {
  .requireEngine("zeitzeiger", 'remotes::install_github("hugheylab/zeitzeiger")')
  stopifnot(inherits(fit, "zeitZeigerFit"))
  if (is.null(rownames(testMat))) {
    stop("testMat must have rownames (gene identifiers)")
  }

  sampleNames <- colnames(testMat)
  nSamples <- ncol(testMat)
  if (is.null(sampleNames)) sampleNames <- seq_len(nSamples)

  outPhase <- rep(NA_real_, nSamples)

  if (all(fit$genes %in% rownames(testMat))) {
    xTest <- t(testMat[fit$genes, , drop = FALSE])
    nSpcUse <- min(nSpc, length(fit$spcResult$d))

    predRes <- zeitzeiger::zeitzeigerPredict(
      fit$xTrain, fit$timeTrain, xTest, fit$spcResult,
      nKnots = fit$nKnots, nSpc = nSpcUse, timeRange = timeRange
    )
    outPhase <- (predRes$timePred[, 1] %% 1) * fit$period
  }

  data.frame(
    sample = sampleNames,
    phase = outPhase,
    confidence = NA_real_,
    row.names = NULL
  )
}
