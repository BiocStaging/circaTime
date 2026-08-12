# T3 (revision_and_fix.md): fidelity anchor for the Molecular Timetable engine.
#
# Molecular Timetable is a reimplementation from the Ueda et al. 2004 paper text,
# not a wrapper around the original authors' software (no maintained third-party
# release exists to wrap - see R/engines-timetable.R). Every other engine in this
# package carries a mandatory fidelity test against an external reference; this is
# that anchor for Molecular Timetable.
#
# Validation tier used: (c), the synthetic ground-truth test (the task's own
# fallback option). (a) was not attempted: Ueda et al. 2004's own mouse liver time
# course is not retrievable in a form this package's harness can reproduce (a
# different, decades-old microarray platform with no public GEO accession this
# package's ingest tooling supports). (b), an independently-coded second
# implementation of the same algorithm, was also not used: the current
# `fitMolecularTimetable()` training step already *is* a direct linear
# least-squares cosinor fit (see R/engines-timetable.R lines 43-50), which is the
# same closed-form solution (b) would ask a second implementation to converge on -
# a second implementation of an already-closed-form linear regression would not be
# independent evidence, it would just be the same normal equations solved twice.
# (c) instead validates the thing that is not already provably correct by
# construction: whether the discretized grid-search *prediction* step recovers a
# known, simulated phase, and whether its accuracy degrades in the expected
# direction (worse) as per-gene noise increases. This is what an implementation
# bug (a sign error, a template misalignment, an off-by-one in the grid) would be
# expected to break.
#
# Design: fit on one sample set, predict on a disjoint held-out sample set drawn
# from the same known-phase generative model (so this is a genuine train/test
# split, not self-recovery on the training data already covered by
# test-engines-timetable.R). Repeated at three noise levels; asserts (1) low noise
# recovers phase tightly, (2) error is monotonically non-decreasing as noise rises,
# (3) even the high-noise condition stays well under the ~6h random-guess floor.

makeFidelityData <- function(nGenes, nSamples, noiseSd, seed) {
  hasOldSeed <- exists(".Random.seed", envir = .GlobalEnv)
  oldSeed <- if (hasOldSeed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (hasOldSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)

  truePhase <- stats::runif(nGenes, 0, 2 * pi)
  trueAmp <- stats::runif(nGenes, 0.5, 3)
  time <- stats::runif(nSamples, 0, 24)
  theta <- time / 24 * 2 * pi

  mat <- t(vapply(seq_len(nGenes), function(i) {
    5 + trueAmp[i] * cos(theta - truePhase[i]) + stats::rnorm(nSamples, sd = noiseSd)
  }, numeric(nSamples)))
  rownames(mat) <- paste0("gene", seq_len(nGenes))
  colnames(mat) <- paste0("s", seq_len(nSamples))

  list(mat = mat, time = time, truePhase = truePhase, trueAmp = trueAmp)
}

# One shared set of true per-gene phase/amplitude so train and test noise draws
# are independent perturbations of the same underlying rhythm, not two unrelated
# simulations.
makeTrainTestPair <- function(nGenes, nTrain, nTest, noiseSd, seed) {
  hasOldSeed <- exists(".Random.seed", envir = .GlobalEnv)
  oldSeed <- if (hasOldSeed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (hasOldSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)

  truePhase <- stats::runif(nGenes, 0, 2 * pi)
  trueAmp <- stats::runif(nGenes, 0.5, 3)

  simulate <- function(nSamples) {
    time <- stats::runif(nSamples, 0, 24)
    theta <- time / 24 * 2 * pi
    mat <- t(vapply(seq_len(nGenes), function(i) {
      5 + trueAmp[i] * cos(theta - truePhase[i]) + stats::rnorm(nSamples, sd = noiseSd)
    }, numeric(nSamples)))
    rownames(mat) <- paste0("gene", seq_len(nGenes))
    list(mat = mat, time = time)
  }

  list(train = simulate(nTrain), test = simulate(nTest))
}

test_that("Molecular Timetable recovers known phase on held-out synthetic samples (low noise)", {
  d <- makeTrainTestPair(nGenes = 60, nTrain = 60, nTest = 40, noiseSd = 0.1, seed = 101)
  fit <- fitMolecularTimetable(d$train$mat, d$train$time, nGenes = 30, minCor = 0.5)
  pred <- predictMolecularTimetable(fit, d$test$mat)

  mae <- circularMAE(pred$phase, d$test$time)
  expect_lt(mae, 1) # tight recovery expected at low noise on held-out samples
})

test_that("Molecular Timetable fidelity degrades in the expected direction as noise rises", {
  noiseLevels <- c(low = 0.1, medium = 0.6, high = 2.0)
  maes <- vapply(noiseLevels, function(sd) {
    d <- makeTrainTestPair(nGenes = 60, nTrain = 60, nTest = 40, noiseSd = sd, seed = 202)
    fit <- fitMolecularTimetable(d$train$mat, d$train$time, nGenes = 30, minCor = 0.3)
    pred <- predictMolecularTimetable(fit, d$test$mat)
    circularMAE(pred$phase, d$test$time)
  }, numeric(1))

  # Monotonic (non-decreasing) degradation is the signature of a correctly wired
  # phase-recovery pipeline; a bug in the grid search or template construction
  # would typically produce an error close to the ~6h random-guess floor
  # regardless of noise, not a graded response to it.
  expect_true(maes["medium"] >= maes["low"] - 1e-8)
  expect_true(maes["high"] >= maes["medium"] - 1e-8)
  # Even the high-noise condition should stay well clear of random guessing.
  expect_lt(maes["high"], 6)
})

test_that("Molecular Timetable fidelity anchor is a genuine train/test split, not self-recovery", {
  # Guards against the anchor accidentally regressing into what
  # test-engines-timetable.R already covers (predicting on the training matrix).
  d <- makeTrainTestPair(nGenes = 40, nTrain = 50, nTest = 30, noiseSd = 0.2, seed = 303)
  expect_false(isTRUE(all.equal(d$train$mat[, seq_len(min(ncol(d$train$mat), ncol(d$test$mat)))],
                                 d$test$mat[, seq_len(min(ncol(d$train$mat), ncol(d$test$mat)))])))
  expect_equal(length(intersect(d$train$time, d$test$time)), 0)
})
