# Tests for R/engines-timetable.R (Molecular Timetable floor baseline).

makeSyntheticTimetableData <- function(nGenes = 30, nSamples = 40, noiseSd = 0.2, seed = 1) {
  hasOldSeed <- exists(".Random.seed", envir = .GlobalEnv)
  oldSeed <- if (hasOldSeed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (hasOldSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)

  time <- stats::runif(nSamples, 0, 24)
  theta <- time / 24 * 2 * pi
  truePhase <- stats::runif(nGenes, 0, 2 * pi)
  trueAmp <- stats::runif(nGenes, 0.5, 3)

  mat <- t(vapply(seq_len(nGenes), function(i) {
    5 + trueAmp[i] * cos(theta - truePhase[i]) + stats::rnorm(nSamples, sd = noiseSd)
  }, numeric(nSamples)))
  rownames(mat) <- paste0("gene", seq_len(nGenes))

  list(mat = mat, time = time)
}

test_that("fitMolecularTimetable selects genes and returns expected structure", {
  d <- makeSyntheticTimetableData()
  fit <- fitMolecularTimetable(d$mat, d$time, nGenes = 10, minCor = 0.3)

  expect_s3_class(fit, "molecularTimetableFit")
  expect_true(length(fit$genes) <= 10)
  expect_true(all(fit$genes %in% rownames(d$mat)))
  expect_true(all(fit$r >= 0.3))
  expect_equal(fit$period, 24)
})

test_that("fitMolecularTimetable requires rownames on trainMat", {
  d <- makeSyntheticTimetableData()
  matNoNames <- d$mat
  rownames(matNoNames) <- NULL
  expect_error(fitMolecularTimetable(matNoNames, d$time), "rownames")
})

test_that("fitMolecularTimetable errors when no gene meets minCor", {
  d <- makeSyntheticTimetableData(noiseSd = 50) # pure noise, no oscillation signal
  expect_error(fitMolecularTimetable(d$mat, d$time, minCor = 0.99), "minCor")
})

test_that("predictMolecularTimetable recovers phase accurately on low-noise synthetic data", {
  d <- makeSyntheticTimetableData(nGenes = 40, nSamples = 60, noiseSd = 0.1)
  fit <- fitMolecularTimetable(d$mat, d$time, nGenes = 20, minCor = 0.5)
  pred <- predictMolecularTimetable(fit, d$mat)

  expect_equal(nrow(pred), ncol(d$mat))
  expect_true(all(c("sample", "phase", "confidence") %in% names(pred)))
  mae <- circularMAE(pred$phase, d$time)
  expect_lt(mae, 1) # well under the ~6h random-guess baseline
})

test_that("predictMolecularTimetable returns NA for samples with too few overlapping genes", {
  d <- makeSyntheticTimetableData()
  fit <- fitMolecularTimetable(d$mat, d$time, nGenes = 10, minCor = 0.3)

  testMat <- d$mat[1:2, , drop = FALSE] # far fewer genes than minGenes default (5)
  pred <- predictMolecularTimetable(fit, testMat, minGenes = 5)
  expect_true(all(is.na(pred$phase)))
})

test_that("predictMolecularTimetable requires rownames on testMat", {
  d <- makeSyntheticTimetableData()
  fit <- fitMolecularTimetable(d$mat, d$time, nGenes = 10, minCor = 0.3)
  testNoNames <- d$mat
  rownames(testNoNames) <- NULL
  expect_error(predictMolecularTimetable(fit, testNoNames), "rownames")
})
