# Tests for R/benchmark.R.

makeSyntheticSE <- function(nGenes = 40, nSamples = 60, noiseSd = 0.15, seed = 1) {
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
  colnames(mat) <- paste0("sample", seq_len(nSamples))

  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = mat),
    colData = S4Vectors::DataFrame(phase = time, row.names = colnames(mat))
  )
}

# For transferPhase(): unlike makeSyntheticSE() (used for within-dataset CV, where
# each call's own random gene-phase relationships are all that matters), a transfer
# test needs train and test sets that share the SAME underlying gene-to-true-phase
# relationship, generated once, with independent noise/sample draws standing in for
# two different "platforms" measuring the same biology - exactly what a real
# cross-platform transfer benchmark assumes. Two independent makeSyntheticSE() calls
# do not do this (each draws its own unrelated truePhase per gene), so training on
# one and testing on the other has no real signal to transfer.
makeTransferPair <- function(nGenes = 40, nTrain = 60, nTest = 30, noiseSd = 0.15, seed = 1) {
  hasOldSeed <- exists(".Random.seed", envir = .GlobalEnv)
  oldSeed <- if (hasOldSeed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (hasOldSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)

  truePhase <- stats::runif(nGenes, 0, 2 * pi)
  trueAmp <- stats::runif(nGenes, 0.5, 3)
  geneNames <- paste0("gene", seq_len(nGenes))

  makeOne <- function(nSamples, prefix) {
    time <- stats::runif(nSamples, 0, 24)
    theta <- time / 24 * 2 * pi
    mat <- t(vapply(seq_len(nGenes), function(i) {
      5 + trueAmp[i] * cos(theta - truePhase[i]) + stats::rnorm(nSamples, sd = noiseSd)
    }, numeric(nSamples)))
    rownames(mat) <- geneNames
    colnames(mat) <- paste0(prefix, seq_len(nSamples))
    SummarizedExperiment::SummarizedExperiment(
      assays = list(expr = mat),
      colData = S4Vectors::DataFrame(phase = time, row.names = colnames(mat))
    )
  }

  list(train = makeOne(nTrain, "train"), test = makeOne(nTest, "test"))
}

test_that(".assignDonorFolds splits by donor, never straddling a donor across folds", {
  donors <- rep(c("d1", "d2", "d3", "d4", "d5", "d6"), each = 3)
  fold <- .assignDonorFolds(donors, nFolds = 3, seed = 42)
  expect_equal(length(fold), length(donors))
  # every sample from the same donor must land in the same fold
  byDonor <- split(fold, donors)
  expect_true(all(vapply(byDonor, function(x) length(unique(x)) == 1, logical(1))))
})

test_that(".assignDonorFolds errors when there are fewer donors than folds", {
  expect_error(.assignDonorFolds(c("a", "b"), nFolds = 5), "Fewer unique donors")
})

test_that(".assignDonorFolds with a seed does not leak global RNG state", {
  set.seed(123)
  before <- runif(1)

  set.seed(123)
  invisible(.assignDonorFolds(rep(letters[1:10], each = 2), nFolds = 4, seed = 999))
  after <- runif(1)

  expect_equal(before, after)
})

test_that("benchmarkPhase runs end-to-end with the timetable engine and beats the random-guess baseline", {
  se <- makeSyntheticSE()
  bm <- benchmarkPhase(se, truth_col = "phase", methods = "timetable",
                        assay_name = "expr", nFolds = 4, seed = 1,
                        engineArgs = list(timetable = list(
                          fit = list(nGenes = 15, minCor = 0.3)
                        )))

  expect_true(all(c("summary", "residuals") %in% names(bm)))
  expect_equal(nrow(bm$summary), 1)
  expect_equal(bm$summary$method, "timetable")
  # random guessing on a 24h circle gives an expected circular MAE of 6h
  expect_lt(bm$summary$circular_mae, 6)
  expect_equal(nrow(bm$residuals), ncol(se))
  expect_false(any(bm$residuals$aligned)) # timetable is supervised, no alignment
})

test_that("benchmarkPhase isolates a single fold's engine failure instead of losing the whole method", {
  # Regression test: a fold whose training data happens to violate an engine's own
  # requirements (e.g. tauFisher's JTK_CYCLE erroring when a CV split removes every
  # replicate of one clock hour) must not discard the other folds' genuine results.
  se <- makeSyntheticSE()
  callCount <- 0
  flakyEngine <- list(
    fit = function(trainMat, trainTime, ...) {
      callCount <<- callCount + 1
      if (callCount == 2) stop("synthetic fold failure")
      list(meanTime = mean(trainTime))
    },
    predict = function(fit, testMat, ...) {
      data.frame(sample = colnames(testMat), phase = rep(fit$meanTime, ncol(testMat)),
                 confidence = NA_real_)
    },
    supervised = TRUE
  )

  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(".circaTimeEngines", c(oldRegistry, list(flaky = flakyEngine)), ns = "circaTime")

  expect_message(
    bm <- benchmarkPhase(se, truth_col = "phase", methods = "flaky",
                          assay_name = "expr", nFolds = 4, seed = 1),
    "fold 2/4 failed"
  )

  expect_true(any(bm$residuals$fold_failed))
  expect_true(any(!bm$residuals$fold_failed))
  expect_true(all(is.na(bm$residuals$pred[bm$residuals$fold_failed])))
  expect_true(all(!is.na(bm$residuals$pred[!bm$residuals$fold_failed])))
  expect_equal(bm$summary$n_failed, sum(bm$residuals$fold_failed))
  expect_equal(bm$summary$n, sum(!is.na(bm$residuals$pred)))
})

test_that("benchmarkPhase errors on an unknown method", {
  se <- makeSyntheticSE()
  expect_error(benchmarkPhase(se, truth_col = "phase", methods = "not_a_method"),
               "Unknown method")
})

test_that("benchmarkPhase errors on a missing truth_col", {
  se <- makeSyntheticSE()
  expect_error(benchmarkPhase(se, truth_col = "not_a_column", methods = "timetable"),
               "not found in colData")
})

test_that("benchmarkPhase threads `period` through to the engine's fit (12h cycle)", {
  set.seed(1)
  nGenes <- 40; nSamples <- 48
  time <- stats::runif(nSamples, 0, 12)
  theta <- time / 12 * 2 * pi
  truePhase <- stats::runif(nGenes, 0, 2 * pi)
  mat <- t(vapply(seq_len(nGenes), function(i) {
    5 + 2 * cos(theta - truePhase[i]) + stats::rnorm(nSamples, sd = 0.15)
  }, numeric(nSamples)))
  rownames(mat) <- paste0("gene", seq_len(nGenes))
  colnames(mat) <- paste0("s", seq_len(nSamples))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = mat),
    colData = S4Vectors::DataFrame(phase = time, row.names = colnames(mat))
  )

  bm <- benchmarkPhase(se, truth_col = "phase", methods = "timetable",
                       assay_name = "expr", nFolds = 4, seed = 1, period = 12,
                       engineArgs = list(timetable = list(fit = list(nGenes = 15, minCor = 0.3))))

  # random guessing on a 12h circle has an expected circular MAE of 3h
  expect_lt(bm$summary$circular_mae, 3)
  expect_true(all(bm$residuals$pred < 12, na.rm = TRUE))
})

test_that("benchmarkPhase warns when a 24h-only engine is asked for a non-24h period", {
  se <- makeSyntheticSE()
  engine24 <- list(
    fit = function(trainMat, trainTime, ...) list(genes = rownames(trainMat)[1:5]),
    predict = function(fit, testMat, ...) {
      data.frame(sample = colnames(testMat), phase = rep(12, ncol(testMat)),
                 confidence = NA_real_)
    },
    supervised = TRUE,
    supportsPeriod = FALSE
  )
  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(".circaTimeEngines", c(oldRegistry, list(engine24 = engine24)), ns = "circaTime")

  expect_warning(
    benchmarkPhase(se, truth_col = "phase", methods = "engine24", assay_name = "expr",
                   nFolds = 2, period = 12),
    "inherently 24-hour"
  )
})

# ---- transferPhase() (AGENTS.md Section 5 Step 2.2) ----

test_that("transferPhase fits on one dataset, scores on another, and beats random guessing", {
  pair <- makeTransferPair(nTrain = 60, nTest = 30, seed = 1)
  trainSe <- pair$train
  testSe <- pair$test

  tb <- transferPhase(trainSe, testSe, truth_col = "phase", methods = "timetable",
                       train_assay_name = "expr", test_assay_name = "expr",
                       engineArgs = list(timetable = list(fit = list(nGenes = 15, minCor = 0.3))))

  expect_true(all(c("summary", "residuals") %in% names(tb)))
  expect_equal(nrow(tb$summary), 1)
  expect_equal(tb$summary$method, "timetable")
  expect_lt(tb$summary$circular_mae, 6) # random-guess baseline on a 24h circle
  expect_equal(nrow(tb$residuals), ncol(testSe))
  expect_equal(tb$summary$n, ncol(testSe))
  expect_equal(tb$summary$n_failed, 0L)
  expect_false(any(tb$residuals$aligned)) # timetable is supervised
})

test_that("transferPhase tolerates a train/test gene panel that only partially overlaps", {
  pair <- makeTransferPair(nGenes = 40, nTrain = 60, nTest = 30, seed = 1)
  trainSe <- pair$train
  testSe <- pair$test
  # Simulate a platform with a smaller gene panel (e.g. RNA-seq vs microarray): keep
  # only half the genes, so trainSe and testSe do NOT share an identical panel.
  testSe <- testSe[rownames(testSe)[1:20], ]

  tb <- transferPhase(trainSe, testSe, truth_col = "phase", methods = "timetable",
                       train_assay_name = "expr", test_assay_name = "expr",
                       engineArgs = list(timetable = list(fit = list(nGenes = 15, minCor = 0.3))))
  expect_equal(tb$summary$n, ncol(testSe))
  expect_true(all(!is.na(tb$residuals$pred)))
})

test_that("transferPhase records a whole-method failure as all-NA predictions with n_failed set", {
  trainSe <- makeSyntheticSE(seed = 1)
  testSe <- makeSyntheticSE(nSamples = 10, seed = 2)
  brokenEngine <- list(
    fit = function(trainMat, trainTime, ...) stop("synthetic fit failure"),
    predict = function(fit, testMat, ...) stop("should not be called"),
    supervised = TRUE
  )
  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(".circaTimeEngines", c(oldRegistry, list(broken = brokenEngine)), ns = "circaTime")

  expect_message(
    tb <- transferPhase(trainSe, testSe, truth_col = "phase", methods = "broken",
                         train_assay_name = "expr", test_assay_name = "expr"),
    "failed"
  )
  expect_true(all(is.na(tb$residuals$pred)))
  expect_equal(tb$summary$n, 0L)
  expect_equal(tb$summary$n_failed, ncol(testSe))
  expect_true(is.nan(tb$summary$circular_mae) || is.na(tb$summary$circular_mae))
})

test_that("transferPhase errors on an unknown method or a missing truth_col", {
  trainSe <- makeSyntheticSE(seed = 1)
  testSe <- makeSyntheticSE(nSamples = 10, seed = 2)
  expect_error(transferPhase(trainSe, testSe, truth_col = "phase", methods = "not_a_method"),
               "Unknown method")
  expect_error(transferPhase(trainSe, testSe, truth_col = "not_a_column", methods = "timetable"),
               "not found in colData")
})
