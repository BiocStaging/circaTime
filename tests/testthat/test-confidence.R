# Tests for R/confidence.R (Phase 3 uncertainty layer).

test_that("coreClockGenes returns the bundled reference table with the AGENTS.md 7.7 gene set", {
  tbl <- coreClockGenes()
  expect_true(all(c("symbol", "aliases", "source") %in% colnames(tbl)))
  expect_setequal(tbl$symbol, c(
    "ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
    "NR1D1", "NR1D2", "DBP", "TEF", "HLF", "CIART", "BHLHE40", "BHLHE41"
  ))
  expect_true(all(nzchar(tbl$source)))
})

makeClockLikeExpr <- function(nGenes = 8, nSamples = 30, tissue = NULL, seed = 1) {
  set.seed(seed)
  genes <- coreClockGenes()$symbol[seq_len(nGenes)]
  phase <- stats::runif(nSamples, 0, 2 * pi)
  genePhase <- stats::runif(nGenes, 0, 2 * pi)
  mat <- t(vapply(seq_len(nGenes), function(i) {
    cos(phase - genePhase[i]) + stats::rnorm(nSamples, sd = 0.1)
  }, numeric(nSamples)))
  rownames(mat) <- genes
  colnames(mat) <- paste0("s", seq_len(nSamples))
  mat
}

test_that("clockAmplitude matches core clock genes case-insensitively", {
  mat <- makeClockLikeExpr()
  rownames(mat) <- tolower(rownames(mat))
  amp <- clockAmplitude(mat)
  expect_length(amp, ncol(mat))
  expect_true(all(amp >= 0))
  expect_equal(length(attr(amp, "genesUsed")), nrow(mat))
})

test_that("clockAmplitude errors with too few matched genes", {
  mat <- matrix(rnorm(2 * 20), nrow = 2, dimnames = list(c("ARNTL", "CLOCK"), paste0("s", 1:20)))
  expect_error(clockAmplitude(mat), "core clock genes")
})

test_that("clockAmplitude is higher for rhythmic signal than for pure noise", {
  rhythmic <- makeClockLikeExpr(seed = 2)
  noise <- rhythmic
  set.seed(3)
  noise[] <- rnorm(length(noise), sd = 1)
  ampRhythmic <- clockAmplitude(rhythmic)
  ampNoise <- clockAmplitude(noise)
  expect_gt(mean(ampRhythmic), mean(ampNoise))
})

test_that("clockAmplitude z-scores within tissue groups separately and returns NA for tiny groups", {
  mat1 <- makeClockLikeExpr(nSamples = 30, seed = 4)
  mat2 <- makeClockLikeExpr(nSamples = 2, seed = 5)
  colnames(mat2) <- paste0("t2_", colnames(mat2))
  mat <- cbind(mat1, mat2)
  tissue <- c(rep("liver", 30), rep("heart", 2))

  amp <- clockAmplitude(mat, tissue = tissue)
  expect_false(any(is.na(amp[tissue == "liver"])))
  expect_true(all(is.na(amp[tissue == "heart"]))) # heart group has only 2 samples
})

# A cheap synthetic engine for bootstrap/calibration tests: point estimate is the
# circular mean of trainTime, jittered by circular noise whose scale depends on
# training-set size (fewer donors -> a shakier estimate -> wider bootstrap spread) -
# a proxy for a real engine's sampling variability without any real computation cost.
makeJitterEngine <- function(sd = 2, failOnReplicate = NULL) {
  callCount <- 0
  list(
    fit = function(trainMat, trainTime, ...) {
      callCount <<- callCount + 1
      if (!is.null(failOnReplicate) && callCount == failOnReplicate) {
        stop("synthetic engine failure")
      }
      list(center = circularMean(trainTime), n = length(trainTime))
    },
    predict = function(fit, testMat, ...) {
      jitter <- wrapPhase(stats::rnorm(ncol(testMat), sd = sd / sqrt(fit$n)))
      data.frame(
        sample = colnames(testMat),
        phase = wrapPhase(fit$center + jitter),
        confidence = NA_real_
      )
    },
    supervised = TRUE
  )
}

test_that("bootstrapPhaseCI returns a point estimate and a widening interval with fewer donors", {
  set.seed(10)
  trainTime <- rep(c(2, 6, 10, 14, 18, 22), each = 4)
  trainMat <- matrix(rnorm(length(trainTime)), nrow = 1, dimnames = list("g", NULL))
  testMat <- matrix(rnorm(3), nrow = 1, dimnames = list("g", paste0("t", 1:3)))
  donors <- paste0("d", seq_along(trainTime))

  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(
    ".circaTimeEngines",
    c(oldRegistry, list(jitter = makeJitterEngine(sd = 3))),
    ns = "circaTime"
  )

  ci <- bootstrapPhaseCI("jitter", trainMat, trainTime, testMat, donors = donors,
                          nBoot = 40, level = 0.8, seed = 1)

  expect_equal(nrow(ci), 3)
  expect_true(all(c("sample", "phase", "phase_lo", "phase_hi", "level",
                     "n_boot_ok", "n_boot_failed") %in% colnames(ci)))
  expect_equal(ci$n_boot_ok, rep(40, 3))
  expect_false(any(is.na(ci$phase_lo)))
  expect_equal(dim(attr(ci, "draws")), c(3, 40))
})

test_that("bootstrapPhaseCI isolates a failed replicate instead of aborting the whole bootstrap", {
  set.seed(11)
  trainTime <- rep(c(2, 6, 10, 14, 18, 22), each = 4)
  trainMat <- matrix(rnorm(length(trainTime)), nrow = 1, dimnames = list("g", NULL))
  testMat <- matrix(rnorm(2), nrow = 1, dimnames = list("g", paste0("t", 1:2)))
  donors <- paste0("d", seq_along(trainTime))

  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  # replicate 3 = the 3rd call to fit(): call 1 is the point-estimate fit, so this
  # fails the 2nd bootstrap replicate specifically.
  utils::assignInNamespace(
    ".circaTimeEngines",
    c(oldRegistry, list(flakyJitter = makeJitterEngine(sd = 3, failOnReplicate = 3))),
    ns = "circaTime"
  )

  expect_message(
    ci <- bootstrapPhaseCI("flakyJitter", trainMat, trainTime, testMat, donors = donors,
                            nBoot = 10, level = 0.8, seed = 1),
    "replicate 2/10 failed"
  )

  expect_equal(ci$n_boot_ok, rep(9, 2))
  expect_equal(ci$n_boot_failed, rep(1, 2))
  expect_false(any(is.na(ci$phase))) # point estimate (replicate 0) always succeeds
})

test_that("bootstrapPhaseCI with a seed does not leak global RNG state", {
  set.seed(123)
  before <- runif(1)

  trainTime <- rep(c(2, 6, 10, 14, 18, 22), each = 4)
  trainMat <- matrix(rnorm(length(trainTime)), nrow = 1, dimnames = list("g", NULL))
  testMat <- matrix(rnorm(2), nrow = 1, dimnames = list("g", paste0("t", 1:2)))

  set.seed(123)
  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(
    ".circaTimeEngines", c(oldRegistry, list(jitter2 = makeJitterEngine(sd = 1))), ns = "circaTime"
  )
  invisible(bootstrapPhaseCI("jitter2", trainMat, trainTime, testMat, nBoot = 5, seed = 999))
  after <- runif(1)

  expect_equal(before, after)
})

test_that("bootstrapPhaseCI errors on an unknown method", {
  expect_error(
    bootstrapPhaseCI("not_a_method", matrix(1, 1, 1), 1, matrix(1, 1, 1)),
    "Unknown method"
  )
})

test_that("residualPool fixes severe undercoverage from a refit-insensitive engine (regression test for the real GSE54650 finding)", {
  # A training-resampling-only bootstrap on a real engine (timetable, on GSE54650)
  # produced a stated 80% interval that actually covered truth only ~48% of the
  # time: refitting on resampled training donors barely changes the prediction, so
  # the bootstrap interval collapses to near-zero width regardless of how noisy the
  # real test-time error is. This engine reproduces that failure mode exactly
  # (fit() ignores the training data, always predicting 12) so the fix
  # (`residualPool`) can be verified without a slow real engine.
  deterministicEngine <- list(
    fit = function(trainMat, trainTime, ...) list(),
    predict = function(fit, testMat, ...) {
      data.frame(sample = colnames(testMat), phase = rep(12, ncol(testMat)), confidence = NA_real_)
    },
    supervised = TRUE
  )
  oldRegistry <- circaTime:::.circaTimeEngines
  on.exit(utils::assignInNamespace(".circaTimeEngines", oldRegistry, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(".circaTimeEngines", c(oldRegistry, list(det = deterministicEngine)), ns = "circaTime")

  nTrials <- 300
  trainMat <- matrix(1, nrow = 1, ncol = 10, dimnames = list("g", paste0("d", 1:10)))
  trainTime <- rep(12, 10)
  testMat <- matrix(1, nrow = 1, ncol = nTrials, dimnames = list("g", paste0("t", 1:nTrials)))

  set.seed(41)
  truth <- wrapPhase(12 + stats::rnorm(nTrials, sd = 2))
  # an independent pool from the same noise-generating process, as a prior
  # cross-validated benchmark run would provide
  residualPool <- stats::rnorm(500, sd = 2)

  ciPlain <- bootstrapPhaseCI("det", trainMat, trainTime, testMat, nBoot = 20, level = 0.8, seed = 1)
  calPlain <- calibrationCurve(ciPlain, truth, levels = 0.8)

  ciResid <- bootstrapPhaseCI("det", trainMat, trainTime, testMat, nBoot = 20, level = 0.8, seed = 1,
                               residualPool = residualPool)
  calResid <- calibrationCurve(ciResid, truth, levels = 0.8)

  expect_lt(calPlain$empirical, 0.3) # severely undercovered without residualPool
  expect_gt(calResid$empirical, 0.6) # residualPool recovers much closer to nominal
})

test_that("calibrationCurve recovers approximately correct empirical coverage from known circular noise", {
  set.seed(20)
  nSamples <- 800
  nBoot <- 60
  sd <- 2
  truth <- stats::runif(nSamples, 0, 24)
  # A valid percentile-bootstrap calibration check needs two *independent* noise
  # draws per sample: the point estimate's own offset from truth, and each
  # bootstrap replicate's offset from that point estimate (not from truth
  # directly - conflating the two makes the replicate spread hugely overstate
  # how far the point estimate itself typically lands from truth, since B=60
  # replicates average out much tighter than any single one, and produces
  # spuriously ~100% coverage regardless of nominal level). This mirrors the
  # real bootstrapPhaseCI(): replicates are refits whose spread approximates the
  # point estimate's own sampling variability, not repeated noisy copies of truth.
  pointEst <- wrapPhase(truth + stats::rnorm(nSamples, sd = sd))
  draws <- t(vapply(pointEst, function(pe) {
    wrapPhase(pe + stats::rnorm(nBoot, sd = sd))
  }, numeric(nBoot)))

  ci <- data.frame(sample = paste0("s", seq_len(nSamples)), phase = pointEst)
  attr(ci, "draws") <- draws

  cal <- calibrationCurve(ci, truth, levels = c(0.5, 0.8, 0.95))
  expect_equal(nrow(cal), 3)
  # empirical coverage should track nominal reasonably closely for a genuinely
  # well-calibrated bootstrap (replicate spread ~= point estimate's own spread)
  expect_lt(max(abs(cal$nominal - cal$empirical)), 0.12)
})

test_that("calibrationCurve errors without a draws attribute", {
  expect_error(calibrationCurve(data.frame(phase = 1), truth = 1), "draws")
})

test_that("plotCalibration runs without error and returns its input invisibly", {
  cal <- data.frame(nominal = c(0.5, 0.8), empirical = c(0.48, 0.81), n = c(100, 100))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  result <- plotCalibration(cal)
  expect_equal(result, cal)
})

test_that("deriveAmplitudeThreshold finds a threshold that separates a known low/high split", {
  set.seed(30)
  n <- 200
  amp <- c(stats::runif(n / 2, 0, 1), stats::runif(n / 2, 2, 5))
  error <- c(stats::rnorm(n / 2, mean = 5, sd = 1), stats::rnorm(n / 2, mean = 1, sd = 1))

  res <- deriveAmplitudeThreshold(amp, error)
  expect_true(res$threshold > 0 && res$threshold < 5)
  expect_gt(res$low_mae, res$high_mae)
  expect_lt(res$p_value, 0.05)
  expect_gt(res$effect_size, 0.5)
})

test_that("deriveAmplitudeThreshold errors with too few pairs", {
  expect_error(deriveAmplitudeThreshold(1:5, 1:5), "at least 20")
})

# clockCoherence() - Phase 3 Step 3.2 dynamic QC. nGenes/nSamples kept small and
# nPerm kept low (these tests only need the RANK of coherence to be right, not a
# precisely calibrated p-value) so the suite stays fast.
makeCoherenceExpr <- function(nGenes = 60, nSamples = 24, dampIdx = integer(0),
                              dampFrac = 0.1, seed = 1) {
  set.seed(seed)
  t <- stats::runif(nSamples, 0, 24)
  truePhase <- stats::runif(nGenes, 0, 24)
  amp <- c(rep(2, nGenes * 0.6), rep(0, nGenes - round(nGenes * 0.6))) # only 60% truly rhythmic
  mat <- t(vapply(seq_len(nGenes), function(g) {
    a <- if (length(dampIdx) == 0) amp[g] else amp[g]
    theta <- 2 * pi * (t - truePhase[g]) / 24
    base <- a * cos(theta) + stats::rnorm(nSamples, sd = 1)
    if (length(dampIdx) > 0) base[dampIdx] <- dampFrac * a * cos(theta[dampIdx]) + stats::rnorm(length(dampIdx), sd = 1)
    base
  }, numeric(nSamples)))
  dimnames(mat) <- list(paste0("g", seq_len(nGenes)), paste0("s", seq_len(nSamples)))
  list(mat = mat, phase = t)
}

test_that("clockCoherence scores damped samples lower than healthy ones", {
  d <- makeCoherenceExpr(nGenes = 150, nSamples = 40, dampIdx = 1:8, dampFrac = 0.1, seed = 42)
  coh <- clockCoherence(d$mat, phase = d$phase, nPerm = 50, seed = 1)
  expect_gt(mean(coh[-(1:8)], na.rm = TRUE), mean(coh[1:8], na.rm = TRUE))
  expect_gte(attr(coh, "panelSize")[["all"]], 20)
  expect_false(isTRUE(attr(coh, "relaxed")[["all"]]))
})

test_that("clockCoherence groups by tissue and selects a panel per group", {
  d1 <- makeCoherenceExpr(nGenes = 100, nSamples = 20, seed = 5)
  d2 <- makeCoherenceExpr(nGenes = 100, nSamples = 20, seed = 6)
  colnames(d2$mat) <- paste0("t2_", colnames(d2$mat))
  mat <- cbind(d1$mat, d2$mat)
  phase <- c(d1$phase, d2$phase)
  tissue <- c(rep("liver", 20), rep("heart", 20))

  coh <- clockCoherence(mat, phase = phase, tissue = tissue, nPerm = 50, seed = 1)
  expect_length(coh, 40)
  expect_setequal(names(attr(coh, "panelSize")), c("liver", "heart"))
  expect_true(all(attr(coh, "panelSize") > 0))
})

test_that("clockCoherence returns NA for groups below minSamples", {
  d <- makeCoherenceExpr(nGenes = 60, nSamples = 8, seed = 7)
  coh <- clockCoherence(d$mat, phase = d$phase, minSamples = 10, nPerm = 20)
  expect_true(all(is.na(coh)))
})

test_that("clockCoherence falls back to top-minPanel when no gene clears the null (relaxed = TRUE)", {
  set.seed(8)
  nGenes <- 30; nSamples <- 20
  mat <- matrix(rnorm(nGenes * nSamples), nrow = nGenes,
                dimnames = list(paste0("g", 1:nGenes), paste0("s", 1:nSamples)))
  phase <- stats::runif(nSamples, 0, 24)
  coh <- clockCoherence(mat, phase = phase, nPerm = 50, minPanel = 15, seed = 1)
  expect_true(attr(coh, "relaxed")[["all"]])
  expect_equal(attr(coh, "panelSize")[["all"]], 15)
  expect_true(is.na(attr(coh, "r2Cutoff")[["all"]]))
})

test_that("clockCoherence validates inputs", {
  mat <- matrix(rnorm(20), nrow = 4, dimnames = list(paste0("g", 1:4), paste0("s", 1:5)))
  expect_error(clockCoherence(mat, phase = 1:3), "length ncol")
  expect_error(clockCoherence(mat, phase = 1:5, tissue = 1:2), "tissue must be")
  expect_error(clockCoherence(mat, phase = 1:5, nPerm = 0), "nPerm")
  expect_error(clockCoherence(mat, phase = 1:5, fdr = 1), "fdr must be")
  noNames <- mat; rownames(noNames) <- NULL
  expect_error(clockCoherence(noNames, phase = 1:5), "rownames")
})

# bayesianPhasePosterior() - Phase 3 Step 3.3.
test_that("bayesianPhasePosterior's credible interval tracks a directly-computed one from the same pool", {
  set.seed(11)
  pool <- circularDiff(rnorm(300, sd = 1.5), 0)
  R <- circularConcentration(pool)
  kappa <- circaTime:::.vonMisesKappaFromR(R)
  direct <- circularCredibleInterval(mu = 6, kappa = kappa, level = 0.8)

  post <- bayesianPhasePosterior(point = 6, residualPool = pool, level = 0.8)
  expect_equal(post$kappa, kappa)
  expect_equal(post$concentration, R)
  expect_equal(post$credible_lo, direct$lo)
  expect_equal(post$credible_hi, direct$hi)
})

test_that("bayesianPhasePosterior's credible interval is close to (not wildly different from) the empirical bootstrap interval on the same pool (AGENTS.md Step 3.3.4)", {
  set.seed(12)
  pool <- circularDiff(rnorm(400, sd = 1.2), 0) # roughly von Mises-shaped by construction
  mu <- 10
  post <- bayesianPhasePosterior(point = mu, residualPool = pool, level = 0.8)
  empirical <- circularQuantileInterval(wrapPhase(mu + pool), level = 0.8)

  bayesWidth <- circularDiff(post$credible_hi, post$credible_lo)
  empWidth <- circularDiff(empirical$hi, empirical$lo)
  expect_lt(abs(bayesWidth - empWidth), 1) # within an hour of each other
  # centres should also roughly agree
  expect_lt(circularError(mean(c(post$credible_lo, post$credible_hi)), mu), 1)
})

test_that("bayesianPhasePosterior returns a density that integrates to ~1 per sample and peaks at the point estimate", {
  set.seed(13)
  pool <- circularDiff(rnorm(300, sd = 1), 0)
  post <- bayesianPhasePosterior(point = c(6, 18), residualPool = pool, nGrid = 241)
  expect_length(post$density, 2)
  for (i in 1:2) {
    d <- post$density[[i]]
    expect_equal(sum(d$density) * (24 / nrow(d)), 1, tolerance = 1e-3)
    peakPhase <- d$phase[which.max(d$density)]
    expect_lt(circularError(peakPhase, c(6, 18)[i]), 24 / nrow(d) * 2)
  }
})

test_that("bayesianPhasePosterior's kappa is a single shared value across multiple points", {
  set.seed(14)
  pool <- circularDiff(rnorm(200, sd = 2), 0)
  post <- bayesianPhasePosterior(point = c(1, 8, 15), residualPool = pool)
  expect_length(post$kappa, 1)
  expect_length(post$credible_lo, 3)
})

test_that("bayesianPhasePosterior errors with too few residual pool values", {
  expect_error(bayesianPhasePosterior(point = 6, residualPool = 1), "at least 2")
})

test_that("clockCoherence is reproducible given the same seed", {
  d <- makeCoherenceExpr(nGenes = 80, nSamples = 20, seed = 9)
  coh1 <- clockCoherence(d$mat, phase = d$phase, nPerm = 30, seed = 42)
  coh2 <- clockCoherence(d$mat, phase = d$phase, nPerm = 30, seed = 42)
  expect_identical(as.numeric(coh1), as.numeric(coh2))
})
