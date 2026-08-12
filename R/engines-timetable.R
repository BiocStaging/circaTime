# Molecular Timetable engine (Ueda et al. 2004, PNAS: 10.1073/pnas.0401882101).
# Historical baseline / floor engine per AGENTS.md Section 2 - implemented from the
# paper's description, since no maintained third-party software release exists to
# wrap. Method: fit a harmonic (cosine) curve per gene against a labelled training
# time course, keep the genes whose fit is good ("time-indicating genes"), then for
# a new single-timepoint sample find the phase whose per-gene reference template
# best correlates with the sample's observed expression of those genes.

#' Fit a Molecular Timetable model from labelled training data
#'
#' @param trainMat Numeric matrix of expression, genes x samples (log-scale
#'   expression is assumed, as for the harmonic-regression fit to be meaningful).
#' @param trainTime Numeric vector of known circadian phase (same units as
#'   `period`, one value per column of `trainMat`).
#' @param period Length of the cycle (default 24, for hours).
#' @param nGenes Maximum number of time-indicating genes to keep, ranked by fit
#'   quality (Pearson correlation between fitted and observed expression). Default
#'   100, following the general magnitude used in the original method.
#' @param minCor Minimum per-gene fit correlation required to be eligible as a
#'   time-indicating gene, evaluated before the `nGenes` cap. Default 0.5.
#' @return An object of class `"molecularTimetableFit"`: a list with the selected
#'   gene names (`genes`), their harmonic-fit coefficients (`mesor`, `b`, `c`), fit
#'   correlation (`r`), and `period`.
#' @examples
#' set.seed(1)
#' nGenes <- 20; nSamples <- 24
#' trainTime <- runif(nSamples, 0, 24)
#' theta <- trainTime / 24 * 2 * pi
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(nSamples, sd = 0.3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' fit <- fitMolecularTimetable(trainMat, trainTime, nGenes = 10)
#' length(fit$genes)
#' @export
fitMolecularTimetable <- function(trainMat, trainTime, period = 24, nGenes = 100,
                                   minCor = 0.5) {
  stopifnot(is.matrix(trainMat), ncol(trainMat) == length(trainTime))
  if (is.null(rownames(trainMat))) {
    stop("trainMat must have rownames (gene identifiers)")
  }

  theta <- trainTime / period * 2 * pi
  design <- cbind(1, cos(theta), sin(theta))

  # Vectorised least-squares fit of every gene at once: coefs is genes x 3
  # (mesor, b, c), solving trainMat ~= coefs %*% t(design) per gene.
  qrDesign <- qr(design)
  coefs <- t(qr.coef(qrDesign, t(trainMat)))
  colnames(coefs) <- c("mesor", "b", "c")

  fittedVals <- design %*% t(coefs) # samples x genes
  r <- vapply(seq_len(nrow(trainMat)), function(i) {
    obs <- trainMat[i, ]
    fit <- fittedVals[, i]
    if (stats::sd(obs) == 0 || stats::sd(fit) == 0) return(0)
    stats::cor(obs, fit)
  }, numeric(1))

  keep <- which(r >= minCor)
  if (length(keep) == 0) {
    stop("No genes met minCor = ", minCor, "; lower minCor or check trainMat/trainTime")
  }
  keep <- keep[order(r[keep], decreasing = TRUE)]
  if (length(keep) > nGenes) keep <- keep[seq_len(nGenes)]

  structure(
    list(
      genes = rownames(trainMat)[keep],
      mesor = coefs[keep, "mesor"],
      b = coefs[keep, "b"],
      c = coefs[keep, "c"],
      r = r[keep],
      period = period
    ),
    class = "molecularTimetableFit"
  )
}

#' Predict circadian phase with a fitted Molecular Timetable model
#'
#' For each test sample, finds the phase (on a fine grid over one cycle) whose
#' per-gene harmonic-fit template best correlates with the sample's observed
#' expression of the fitted model's time-indicating genes.
#'
#' @param fit A `"molecularTimetableFit"` object from [fitMolecularTimetable()].
#' @param testMat Numeric matrix of expression, genes x samples, on the same scale
#'   as the training matrix. Only genes overlapping `fit$genes` are used per sample;
#'   genes absent from `testMat` are silently dropped from that computation.
#' @param gridN Number of grid points spanning one cycle to search over (default
#'   288, i.e. 5-minute resolution over 24h).
#' @param minGenes Minimum number of overlapping time-indicating genes required to
#'   attempt a prediction for a sample; below this, the sample gets `NA` (default 5).
#' @return A `data.frame` with one row per column of `testMat`: `phase` (numeric,
#'   `[0, period)`, or `NA` if prediction failed), and `confidence` (the maximum
#'   per-sample correlation achieved across the grid, or `NA`).
#' @examples
#' set.seed(1)
#' nGenes <- 20; nSamples <- 24
#' trainTime <- runif(nSamples, 0, 24)
#' theta <- trainTime / 24 * 2 * pi
#' truePhase <- runif(nGenes, 0, 2 * pi)
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - truePhase[i]) + rnorm(nSamples, sd = 0.1)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' fit <- fitMolecularTimetable(trainMat, trainTime, nGenes = 10, minCor = 0.3)
#' pred <- predictMolecularTimetable(fit, trainMat)
#' head(pred)
#' @export
predictMolecularTimetable <- function(fit, testMat, gridN = 288, minGenes = 5) {
  stopifnot(inherits(fit, "molecularTimetableFit"))
  if (is.null(rownames(testMat))) {
    stop("testMat must have rownames (gene identifiers)")
  }

  period <- fit$period
  common <- intersect(fit$genes, rownames(testMat))

  grid <- seq(0, period, length.out = gridN + 1)[seq_len(gridN)]
  thetaGrid <- grid / period * 2 * pi

  nSamples <- ncol(testMat)
  outPhase <- rep(NA_real_, nSamples)
  outConf <- rep(NA_real_, nSamples)

  if (length(common) >= minGenes) {
    idx <- match(common, fit$genes)
    mesor <- fit$mesor[idx]; b <- fit$b[idx]; c <- fit$c[idx]
    # template: genes(common) x gridN
    template <- mesor + outer(b, cos(thetaGrid)) + outer(c, sin(thetaGrid))
    testSub <- testMat[common, , drop = FALSE]

    for (j in seq_len(nSamples)) {
      obs <- testSub[, j]
      valid <- is.finite(obs)
      if (sum(valid) < minGenes) next
      cors <- stats::cor(obs[valid], template[valid, , drop = FALSE])
      best <- which.max(cors)
      outPhase[j] <- grid[best]
      outConf[j] <- cors[best]
    }
  }

  sampleNames <- colnames(testMat)
  if (is.null(sampleNames)) sampleNames <- seq_len(nSamples)

  data.frame(
    sample = sampleNames,
    phase = outPhase,
    confidence = outConf,
    row.names = NULL
  )
}
