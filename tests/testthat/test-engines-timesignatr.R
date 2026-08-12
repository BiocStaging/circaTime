# Tests for R/engines-timesignatr.R.
#
# Unlike tauFisher (whose example data ships inside the installed package), the
# authors' own worked example - real Moller/Archer/Arnardottir human blood data -
# lives only in the GitHub repo's example/DATA/TSexampleData.Rdata (56 MB, so it is
# downloaded here rather than bundled: AGENTS.md caps in-package data objects at
# 5 MB). The AGENTS.md Section 9 Phase 2 mandatory fidelity test is below
# ("... reproduces TimeSignatR's own real training/test split ..."); it downloads
# that file into the session temp dir (skipping, not failing, if the network/GitHub
# is unreachable) and is the one real-data check for this engine. Every other test
# in this file uses synthetic data.

testthat::skip_if_not_installed("TimeSignatR")

makeSyntheticTimeSignatRData <- function(nGenes = 40, seed = 1) {
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

test_that("fitTimeSignatR returns expected structure on synthetic oscillatory data", {
  d <- makeSyntheticTimeSignatRData()
  fit <- fitTimeSignatR(d$mat, d$time, seed = 1)

  expect_s3_class(fit, "timeSignatRFit")
  expect_equal(fit$genes, rownames(d$mat))
  expect_equal(fit$period, 24)
  expect_true(all(fit$timestamp$train)) # trainFrac = 1 - no internal holdout
  expect_s3_class(fit$timestamp$cv.fit, "cv.glmnet")
})

test_that("fitTimeSignatR requires rownames on trainMat", {
  d <- makeSyntheticTimeSignatRData()
  matNoNames <- d$mat
  rownames(matNoNames) <- NULL
  expect_error(fitTimeSignatR(matNoNames, d$time), "rownames")
})

test_that("fitTimeSignatR with a seed does not leak global RNG state", {
  d <- makeSyntheticTimeSignatRData()
  set.seed(123)
  before <- runif(1)

  set.seed(123)
  invisible(fitTimeSignatR(d$mat, d$time, seed = 999))
  after <- runif(1)

  expect_equal(before, after)
})

test_that("predictTimeSignatR returns expected structure and phases within range", {
  d <- makeSyntheticTimeSignatRData()
  fit <- fitTimeSignatR(d$mat, d$time, seed = 1)
  pred <- predictTimeSignatR(fit, d$mat)

  expect_equal(nrow(pred), ncol(d$mat))
  expect_true(all(c("sample", "phase", "confidence") %in% names(pred)))
  expect_true(all(is.na(pred$confidence))) # no native confidence signal - see docs
  expect_true(all(pred$phase >= 0 & pred$phase < 24))
})

test_that("predictTimeSignatR returns NA for every sample when fitted genes are missing from testMat", {
  d <- makeSyntheticTimeSignatRData()
  fit <- fitTimeSignatR(d$mat, d$time, seed = 1)

  testMat <- d$mat[1:5, , drop = FALSE] # far fewer genes than fit$genes
  pred <- predictTimeSignatR(fit, testMat)
  expect_true(all(is.na(pred$phase)))
})

test_that("predictTimeSignatR requires rownames on testMat", {
  d <- makeSyntheticTimeSignatRData()
  fit <- fitTimeSignatR(d$mat, d$time, seed = 1)
  testNoNames <- d$mat
  rownames(testNoNames) <- NULL
  expect_error(predictTimeSignatR(fit, testNoNames), "rownames")
})

test_that("fitTimeSignatR rescales a non-24 period to hours before fitting", {
  d <- makeSyntheticTimeSignatRData()
  timeFrac <- d$time / 24 # same phases expressed on a [0, 1) period

  fit24 <- fitTimeSignatR(d$mat, d$time, period = 24, seed = 1)
  fit1 <- fitTimeSignatR(d$mat, timeFrac, period = 1, seed = 1)

  pred24 <- predictTimeSignatR(fit24, d$mat)
  pred1 <- predictTimeSignatR(fit1, d$mat)

  expect_equal(pred24$phase / 24, pred1$phase, tolerance = 1e-6)
})

test_that("fitTimeSignatR + predictTimeSignatR reproduce TimeSignatR's own real training/test split (mandatory fidelity test, AGENTS.md 9)", {
  # Not skip_on_cran(): this is the mandatory Phase 2 fidelity gate (AGENTS.md
  # Section 9). It IS allowed to skip when the data can't be downloaded (offline
  # sandbox, GitHub unreachable) - that's an environment limitation, not a reason
  # to treat the gate as optional; run it with network access before release.
  cacheFile <- file.path(tools::R_user_dir("circaTime", "cache"), "TSexampleData.Rdata")
  if (!file.exists(cacheFile)) {
    dir.create(dirname(cacheFile), recursive = TRUE, showWarnings = FALSE)
    ok <- tryCatch({
      utils::download.file(
        "https://raw.githubusercontent.com/braunr/TimeSignatR/master/example/DATA/TSexampleData.Rdata",
        destfile = cacheFile, mode = "wb", quiet = TRUE
      )
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(ok)) {
      unlink(cacheFile)
      skip("TSexampleData.Rdata could not be downloaded (no network/GitHub access)")
    }
  }

  envData <- new.env()
  load(cacheFile, envir = envData)
  allExpr <- envData$all.expr
  allMeta <- envData$all.meta

  # The Moller et al. train/test/subset ("TrTe" study) split is pre-recorded in the
  # data itself (all.meta$train) - the authors' own script re-derives an identical
  # split from a fixed seed for illustration, but using the recorded column directly
  # is simpler and exactly as faithful.
  studyIdx <- allMeta$study == "TrTe"
  trainIdx <- studyIdx & allMeta$train == 1
  testIdx <- studyIdx & allMeta$train == 0

  trainMat <- allExpr[, trainIdx]
  testMat <- allExpr[, testIdx]
  trainTime <- allMeta$LocalTime[trainIdx]
  testTime <- allMeta$LocalTime[testIdx]
  subjIDs <- allMeta$ID[trainIdx]

  # Gold values: calling TimeSignatR's own trainTimeStamp()/predTimeStamp() directly,
  # bypassing our wrapper entirely, with the package's own default settings
  # (a = 0.5, recalib = FALSE) and trainFrac = 1 (the split above is already the
  # full intended training set, so there is nothing left for TimeSignatR's own
  # internal holdout to hold out).
  hasOldSeed <- exists(".Random.seed", envir = .GlobalEnv)
  oldSeed <- if (hasOldSeed) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(42)
  goldFit <- TimeSignatR::trainTimeStamp(
    expr = trainMat, subjIDs = subjIDs, times = trainTime,
    trainFrac = 1, a = 0.5, recalib = FALSE
  )
  goldPred <- as.numeric(TimeSignatR::predTimeStamp(goldFit, newx = testMat))
  if (hasOldSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)

  fit <- fitTimeSignatR(trainMat, trainTime, subjIDs = subjIDs, seed = 42)
  pred <- predictTimeSignatR(fit, testMat)

  expect_equal(pred$phase, goldPred)

  # Sanity check that this is a real, working predictor on real data, not just a
  # reproduction of a broken one: circular MAE against the true local time should
  # be well under a 12h coin-flip on this real held-out test set.
  circErr <- pmin(abs(pred$phase - testTime), 24 - abs(pred$phase - testTime))
  expect_true(mean(circErr) < 6)
})
