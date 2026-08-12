# Tests for R/engines-zeitzeiger.R.
#
# zeitzeiger ships no bundled real/example dataset in its installed package - only a
# simulated-data vignette (github.com/hugheylab/zeitzeiger/vignettes/introduction.Rmd),
# not real biological data. The AGENTS.md Section 9 Phase 2 mandatory fidelity test
# below reproduces that vignette's simulated dataset exactly (same seed, same
# generative code) and compares this wrapper's output to the authors' own functions
# called directly on the identical data - a genuine reproduction of "the original
# implementation's output on the original authors' own example data," just not real
# biological data, because the authors' own example isn't either. Every other test in
# this file uses its own smaller synthetic data for speed.

testthat::skip_if_not_installed("zeitzeiger")

makeSyntheticZeitZeigerData <- function(nGenes = 40, seed = 1) {
  hasOldSeed <- exists(".Random.seed", envir = .GlobalEnv)
  oldSeed <- if (hasOldSeed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (hasOldSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)

  time <- seq(0, 23, by = 1)
  theta <- time / 24 * 2 * pi
  mat <- t(vapply(seq_len(nGenes), function(i) {
    5 + 2 * cos(theta - stats::runif(1, 0, 2 * pi)) + stats::rnorm(length(time), sd = 0.3)
  }, numeric(length(time))))
  rownames(mat) <- paste0("gene", seq_len(nGenes))
  colnames(mat) <- paste0("s", seq_len(length(time)))

  list(mat = mat, time = time)
}

test_that("fitZeitZeiger returns expected structure on synthetic oscillatory data", {
  d <- makeSyntheticZeitZeigerData()
  fit <- fitZeitZeiger(d$mat, d$time)

  expect_s3_class(fit, "zeitZeigerFit")
  expect_equal(fit$genes, rownames(d$mat))
  expect_equal(dim(fit$xTrain), c(ncol(d$mat), nrow(d$mat)))
  expect_true(all(fit$timeTrain >= 0 & fit$timeTrain < 1))
  expect_equal(fit$period, 24)
  expect_true(length(fit$spcResult$d) >= 1)
})

test_that("fitZeitZeiger requires rownames on trainMat", {
  d <- makeSyntheticZeitZeigerData()
  matNoNames <- d$mat
  rownames(matNoNames) <- NULL
  expect_error(fitZeitZeiger(matNoNames, d$time), "rownames")
})

test_that("fitZeitZeiger wraps trainTime into [0, period) before scaling to [0, 1)", {
  d <- makeSyntheticZeitZeigerData()
  fitBase <- fitZeitZeiger(d$mat, d$time)
  fitWrapped <- fitZeitZeiger(d$mat, d$time + 24)
  expect_equal(fitBase$timeTrain, fitWrapped$timeTrain)
})

test_that("predictZeitZeiger returns expected structure and phases within range", {
  d <- makeSyntheticZeitZeigerData()
  fit <- fitZeitZeiger(d$mat, d$time)
  pred <- predictZeitZeiger(fit, d$mat)

  expect_equal(nrow(pred), ncol(d$mat))
  expect_true(all(c("sample", "phase", "confidence") %in% names(pred)))
  expect_true(all(is.na(pred$confidence))) # no native confidence signal - see docs
  expect_true(all(pred$phase >= 0 & pred$phase < 24))
})

test_that("predictZeitZeiger returns NA for every sample when fitted genes are missing from testMat", {
  d <- makeSyntheticZeitZeigerData()
  fit <- fitZeitZeiger(d$mat, d$time)

  testMat <- d$mat[1:5, , drop = FALSE] # far fewer genes than fit$genes
  pred <- predictZeitZeiger(fit, testMat)
  expect_true(all(is.na(pred$phase)))
})

test_that("predictZeitZeiger requires rownames on testMat", {
  d <- makeSyntheticZeitZeigerData()
  fit <- fitZeitZeiger(d$mat, d$time)
  testNoNames <- d$mat
  rownames(testNoNames) <- NULL
  expect_error(predictZeitZeiger(fit, testNoNames), "rownames")
})

test_that("predictZeitZeiger clamps nSpc to the number of available SPCs", {
  d <- makeSyntheticZeitZeigerData(nGenes = 3) # fewer genes than default nTime = 10
  fit <- fitZeitZeiger(d$mat, d$time)
  expect_true(length(fit$spcResult$d) < 10)

  pred <- predictZeitZeiger(fit, d$mat, nSpc = 10) # should not error despite nSpc > K
  expect_equal(nrow(pred), ncol(d$mat))
})

test_that("fitZeitZeiger + predictZeitZeiger reproduce zeitzeiger's own vignette example (mandatory fidelity test, AGENTS.md 9)", {
  # Reproduces github.com/hugheylab/zeitzeiger vignettes/introduction.Rmd exactly:
  # same nObs/nFeatures/seed and the same three-signal generative code, so xTrain/
  # xTest are bit-for-bit what the vignette itself works with. zeitzeigerFit/Spc/
  # Predict are deterministic given their inputs (no RNG use of their own), so gold
  # values from calling them directly are reproducible without a stored fixture.
  nObs <- 100; nFeatures <- 200; amp <- 0.5
  time <- rep((0:9) / 10, length.out = nObs)
  set.seed(42)
  x <- matrix(stats::rnorm(nObs * nFeatures, sd = 0.7), nrow = nObs)

  baseSignal <- sin(time * 2 * pi)
  x[, 1] <- x[, 1] + 3 * baseSignal
  x[, 2] <- x[, 2] - 0.5 * baseSignal
  for (ii in 3:10) x[, ii] <- x[, ii] + stats::runif(1, -amp, amp) * baseSignal

  baseSignal <- cos(time * 2 * pi)
  x[, 11] <- x[, 11] + 2 * baseSignal
  x[, 12] <- x[, 12] - baseSignal
  for (ii in 13:20) x[, ii] <- x[, ii] + stats::runif(1, -amp, amp) * baseSignal

  baseSignal <- cos(time * 4 * pi + pi / 6)
  x[, 21] <- x[, 21] + baseSignal
  x[, 22] <- x[, 22] - baseSignal
  for (ii in 23:30) x[, ii] <- x[, ii] + stats::runif(1, -amp, amp) * baseSignal

  idxTrain <- 1:round(nObs * 0.6)
  xTrain <- x[idxTrain, , drop = FALSE]; timeTrain <- time[idxTrain]
  xTest <- x[-idxTrain, , drop = FALSE]; timeTest <- time[-idxTrain]

  # Gold: the vignette's own calls, bypassing our wrapper entirely.
  fitResult <- zeitzeiger::zeitzeigerFit(xTrain, timeTrain)
  spcResult <- zeitzeiger::zeitzeigerSpc(fitResult$xFitMean, fitResult$xFitResid, sumabsv = 1)
  predResult <- zeitzeiger::zeitzeigerPredict(xTrain, timeTrain, xTest, spcResult, nSpc = 2)
  goldPred <- as.numeric(predResult$timePred[, 1])

  # circaTime convention: genes(features) x samples(observations), so transpose the
  # vignette's obs x features matrices; period = 1 since the vignette's time is
  # already fractional [0, 1), matching zeitzeiger's own native scale.
  trainMat <- t(xTrain)
  rownames(trainMat) <- paste0("feature", seq_len(nFeatures))
  colnames(trainMat) <- paste0("s", idxTrain)
  testMat <- t(xTest)
  rownames(testMat) <- paste0("feature", seq_len(nFeatures))
  colnames(testMat) <- paste0("s", setdiff(seq_len(nObs), idxTrain))

  fit <- fitZeitZeiger(trainMat, timeTrain, period = 1, sumabsv = 1)
  pred <- predictZeitZeiger(fit, testMat, nSpc = 2)

  expect_equal(pred$phase, goldPred)

  # Sanity check this is a real, working predictor, not a reproduction of a broken
  # one: circular MAE (as a fraction of the cycle) should be well under a
  # coin-flip's 0.25 on this vignette's own (fairly easy) simulated test set.
  circErr <- pmin(abs(pred$phase - timeTest), 1 - abs(pred$phase - timeTest))
  expect_true(mean(circErr) < 0.1)
})
