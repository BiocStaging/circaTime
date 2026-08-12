# Tests for R/engines-taufisher.R.
#
# Includes the AGENTS.md Section 9 Phase 2 mandatory fidelity test: the wrapped
# engine must reproduce tauFisher's own SingleDataset vignette result on the
# authors' own bundled GSE157077 example, using the authors' own recommended
# settings. This test genuinely invokes MetaCycle (via fitTauFisher()) on ~24,800
# genes and takes a few minutes - that cost is inherent to testing the real engine,
# not a synthetic stand-in.

testthat::skip_if_not_installed("tauFisher")

makeSyntheticTauFisherData <- function(nGenes = 60, seed = 1) {
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
    100 + 40 * cos(theta - stats::runif(1, 0, 2 * pi)) + stats::rnorm(length(time), sd = 3)
  }, numeric(length(time))))
  # Uppercase, since train_tauFisher() uppercases rownames(train_df) internally
  # (tauFisher assumes gene-symbol conventions) - lowercase names here would make
  # fit$Genes fail to match rownames(testMat) purely on casing, not biology.
  rownames(mat) <- paste0("GENE", seq_len(nGenes))
  colnames(mat) <- paste0("s", seq_len(length(time)))

  list(mat = mat, time = time)
}

test_that("fitTauFisher returns expected structure on synthetic oscillatory data", {
  d <- makeSyntheticTauFisherData()
  fit <- fitTauFisher(d$mat, d$time, seed = 1)

  expect_s3_class(fit, "tauFisherFit")
  expect_true(length(fit$Genes) >= 2)
  expect_true(all(fit$Genes %in% rownames(d$mat)))
  expect_true(fit$method %in% c("JTK", "LS"))
  expect_equal(fit$period, 24)
})

test_that("fitTauFisher requires rownames on trainMat", {
  d <- makeSyntheticTauFisherData()
  matNoNames <- d$mat
  rownames(matNoNames) <- NULL
  expect_error(fitTauFisher(matNoNames, d$time), "rownames")
})

test_that("fitTauFisher with a seed does not leak global RNG state", {
  d <- makeSyntheticTauFisherData()
  set.seed(123)
  before <- runif(1)

  set.seed(123)
  invisible(fitTauFisher(d$mat, d$time, seed = 999))
  after <- runif(1)

  expect_equal(before, after)
})

test_that("predictTauFisher returns expected structure and predicts within-training-range hours", {
  d <- makeSyntheticTauFisherData()
  fit <- fitTauFisher(d$mat, d$time, seed = 1)
  pred <- predictTauFisher(fit, d$mat)

  expect_equal(nrow(pred), ncol(d$mat))
  expect_true(all(c("sample", "phase", "confidence") %in% names(pred)))
  expect_true(all(is.na(pred$confidence))) # no native confidence signal - see docs
  expect_true(all(pred$phase >= 0 & pred$phase < 24))
})

test_that("predictTauFisher returns NA for every sample when fitted genes are missing from testMat", {
  d <- makeSyntheticTauFisherData()
  fit <- fitTauFisher(d$mat, d$time, seed = 1)

  testMat <- d$mat[1, , drop = FALSE] # far fewer genes than fit$Genes
  pred <- predictTauFisher(fit, testMat)
  expect_true(all(is.na(pred$phase)))
})

test_that("predictTauFisher requires rownames on testMat", {
  d <- makeSyntheticTauFisherData()
  fit <- fitTauFisher(d$mat, d$time, seed = 1)
  testNoNames <- d$mat
  rownames(testNoNames) <- NULL
  expect_error(predictTauFisher(fit, testNoNames), "rownames")
})

test_that("predictTauFisher works in a fresh session that never called fitTauFisher (regression test)", {
  # Real bug found while validating Phase 3 Step 3.1's consensus vignette chunk:
  # tauFisher's own NAMESPACE does not import nnet, so nnet's S3 method table
  # (predict.multinom, which the fitted $Model dispatches to) is only ever
  # populated as a side effect of fitTauFisher() -> nnet::multinom() having run
  # IN THIS SESSION. Every other test in this file calls fitTauFisher() first,
  # which masks the bug completely - the real-world failure mode is a
  # prediction-only session against a bundled reference model (exactly what
  # estimatePhase() does), which never trains anything itself. Reproducing that
  # honestly needs a genuinely fresh subprocess; running inline here would just
  # inherit nnet already loaded by an earlier test in this file/session.
  testthat::skip_if_not_installed("callr")
  testthat::skip_if_not_installed("pkgload")
  d <- makeSyntheticTauFisherData()
  fitPath <- tempfile(fileext = ".rds")
  saveRDS(fitTauFisher(d$mat, d$time, seed = 1), fitPath)

  pkgPath <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  result <- callr::r(
    function(pkgPath, fitPath, testMat) {
      if (!requireNamespace("circaTime", quietly = TRUE)) {
        pkgload::load_all(pkgPath, quiet = TRUE)
      } else {
        library(circaTime)
      }
      fit <- readRDS(fitPath)
      predictTauFisher(fit, testMat)
    },
    args = list(pkgPath = pkgPath, fitPath = fitPath, testMat = d$mat)
  )
  expect_true(all(result$phase >= 0 & result$phase < 24))
})

test_that("fitTauFisher + predictTauFisher reproduce tauFisher's own SingleDataset vignette result (mandatory fidelity test, AGENTS.md 9)", {
  # Deliberately NOT skip_on_cran()/skip_if(): this is the mandatory Phase 2
  # fidelity gate (AGENTS.md Section 9) and must always run, slow as it is.

  datafile <- system.file("extdata", "GSE157077_mouse_scn_control.tsv",
                           package = "tauFisher", mustWork = TRUE)
  raw <- utils::read.delim(file = datafile, stringsAsFactors = FALSE)

  # Collapse non-unique gene symbols by mean expression, as the vignette does.
  raw$ID <- toupper(raw$ID)
  raw <- raw[stats::complete.cases(raw[, -1]), ]
  df <- stats::aggregate(raw[, -1], by = list(ID = raw$ID), FUN = mean)
  rownames(df) <- df$ID
  df <- df[, -1]

  time <- as.numeric(vapply(strsplit(colnames(df), "_"), `[`, 2, FUN.VALUE = character(1)))
  replicate <- as.numeric(vapply(strsplit(colnames(df), "_"), `[`, 4, FUN.VALUE = character(1)))

  # Train on replicates 1-2, test on replicate 3, exactly as the vignette does.
  trainIdx <- which(replicate != 3)
  trainTimeAdj <- time[trainIdx] + 24 * (replicate[trainIdx] - 1)
  trainAdj <- df[, trainIdx]
  colnames(trainAdj) <- trainTimeAdj
  trainAdj <- trainAdj[, order(as.numeric(colnames(trainAdj)))]
  trainTimeAdj <- sort(trainTimeAdj)

  testDf <- df[, -trainIdx]
  testTime <- time[-trainIdx]

  fit <- fitTauFisher(as.matrix(trainAdj), trainTimeAdj,
                       method = c("JTK", "LS"), thres = 1.01,
                       nrep = 2, numbasis = 5, seed = 123)
  pred <- predictTauFisher(fit, as.matrix(testDf))

  # Gold values: running the vignette's own direct train_tauFisher()/test_tauFisher()
  # calls (bypassing our wrapper entirely) on this exact train/test split with the
  # same seed = 123 immediately before the training call.
  goldPhase <- c(22, 0, 9, 13, 14, 21)
  goldTruth <- c(0, 4, 8, 12, 16, 20)

  expect_equal(testTime, goldTruth)
  expect_equal(pred$phase, goldPhase)
  expect_equal(fit$method, "JTK") # LS finds nothing significant on this dataset
})
