# Tests for R/phase-confounding.R.
#
# The behaviour that actually matters here is discrimination: a design where
# collection time tracks the condition must come out significant, and one where phase
# is balanced across groups must not - a test that only checked "returns a data.frame"
# would pass on a function that always reported confounding. The wraparound case
# (groups at ZT23 vs ZT1) is the circular-correctness check AGENTS.md 7.1 demands:
# those groups are 2h apart, not 22h, and any Euclidean implementation gets it wrong.

makeConfoundingSE <- function(phase, ...) {
  n <- length(phase)
  mat <- matrix(
    stats::rnorm(50 * n), nrow = 50,
    dimnames = list(paste0("g", seq_len(50)), paste0("s", seq_len(n)))
  )
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  SummarizedExperiment::colData(se)$phase <- wrapPhase(phase)
  extras <- list(...)
  for (nm in names(extras)) {
    SummarizedExperiment::colData(se)[[nm]] <- extras[[nm]]
  }
  se
}

test_that("a design where collection time tracks condition is flagged", {
  set.seed(1)
  se <- makeConfoundingSE(
    c(stats::rnorm(10, 6, 1), stats::rnorm(10, 14, 1)),
    condition = rep(c("case", "control"), each = 10)
  )
  res <- phaseConfounding(se, ~condition, nPerm = 999, seed = 1)

  expect_equal(nrow(res), 1L)
  expect_equal(res$term, "condition")
  expect_equal(res$type, "categorical")
  expect_equal(res$n, 20L)
  expect_equal(res$n_groups, 2L)
  expect_lt(res$p_value, 0.01)
  # cases ~ZT6, controls ~ZT14: the reported gap should recover that ~8h separation
  expect_gt(res$effect_size, 6)
  expect_lt(res$effect_size, 10)
})

test_that("a design with phase balanced across groups is not flagged", {
  set.seed(1)
  se <- makeConfoundingSE(
    rep(seq(0, 22, by = 2), 2),
    condition = rep(c("case", "control"), each = 12)
  )
  res <- phaseConfounding(se, ~condition, nPerm = 999, seed = 1)

  expect_gt(res$p_value, 0.05)
  expect_lt(res$effect_size, 1)
})

test_that("group separation is measured circularly across the 0/24 boundary", {
  set.seed(1)
  # ZT23 vs ZT1: 2h apart the short way round, 22h the wrong way round.
  se <- makeConfoundingSE(
    c(stats::rnorm(10, 23, 0.5), stats::rnorm(10, 25, 0.5)),
    condition = rep(c("a", "b"), each = 10)
  )
  res <- phaseConfounding(se, ~condition, nPerm = 999, seed = 1)

  expect_lt(res$effect_size, 4)
  gm <- attr(res, "groupMeans")
  expect_equal(nrow(gm), 2L)
  expect_true(all(gm$circular_mean >= 0 & gm$circular_mean < 24))
  # the ZT25 group must wrap to ~1, not stay at 25
  expect_lt(min(gm$circular_mean), 3)
})

test_that("numeric terms are tested by circular-linear correlation", {
  set.seed(1)
  n <- 30
  ph <- stats::runif(n, 0, 24)
  se <- makeConfoundingSE(
    ph,
    age = 50 + 10 * cos(ph / 24 * 2 * pi) + stats::rnorm(n, sd = 0.5),
    noise = stats::rnorm(n)
  )
  res <- phaseConfounding(se, ~ age + noise, nPerm = 999, seed = 1)

  expect_equal(res$type, c("numeric", "numeric"))
  expect_true(all(is.na(res$n_groups)))
  # age was constructed as a function of phase; noise was not
  expect_lt(res$p_value[res$term == "age"], 0.01)
  expect_gt(res$statistic[res$term == "age"], 0.8)
  expect_gt(res$p_value[res$term == "noise"], 0.05)
  # no categorical terms, so no group means to report
  expect_null(attr(res, "groupMeans"))
})

test_that("only 'ok'-flagged samples are used, and the drop count is reported", {
  set.seed(1)
  se <- makeConfoundingSE(
    c(stats::rnorm(10, 6, 1), stats::rnorm(10, 14, 1)),
    condition = rep(c("case", "control"), each = 10),
    phase_flag = rep(c("ok", "low_amplitude"), times = 10)
  )
  res <- phaseConfounding(se, ~condition, nPerm = 199, seed = 1)

  expect_equal(res$n, 10L)
  expect_equal(attr(res, "nDropped"), 10L)

  # flag_col = NULL opts out of the filtering entirely
  resAll <- phaseConfounding(se, ~condition, flag_col = NULL, nPerm = 199, seed = 1)
  expect_equal(resAll$n, 20L)
  expect_equal(attr(resAll, "nDropped"), 0L)
})

test_that("interaction terms are skipped with a warning, main effects kept", {
  set.seed(1)
  n <- 30
  se <- makeConfoundingSE(
    stats::runif(n, 0, 24),
    condition = rep(c("a", "b", "c"), each = 10),
    age = stats::rnorm(n)
  )
  expect_warning(
    res <- phaseConfounding(se, ~ condition * age, nPerm = 99, seed = 1),
    "Skipping interaction term"
  )
  expect_equal(res$term, c("condition", "age"))
})

test_that("phaseConfounding does not leak global RNG state", {
  set.seed(1)
  se <- makeConfoundingSE(
    stats::runif(20, 0, 24),
    condition = rep(c("a", "b"), each = 10)
  )
  set.seed(42)
  before <- .Random.seed
  invisible(phaseConfounding(se, ~condition, nPerm = 99, seed = 1))
  expect_identical(before, .Random.seed)

  # and without a seed argument either
  set.seed(42)
  before <- .Random.seed
  invisible(phaseConfounding(se, ~condition, nPerm = 99))
  expect_identical(before, .Random.seed)
})

test_that("seeded runs are reproducible", {
  set.seed(1)
  se <- makeConfoundingSE(
    stats::runif(20, 0, 24),
    condition = rep(c("a", "b"), each = 10)
  )
  a <- phaseConfounding(se, ~condition, nPerm = 199, seed = 7)
  b <- phaseConfounding(se, ~condition, nPerm = 199, seed = 7)
  expect_equal(a$p_value, b$p_value)
})

test_that("malformed input errors informatively", {
  set.seed(1)
  se <- makeConfoundingSE(
    stats::runif(20, 0, 24),
    condition = rep(c("a", "b"), each = 10)
  )
  expect_error(phaseConfounding(se, "not a formula"), "must be a formula")
  expect_error(phaseConfounding(se, ~nope, nPerm = 9), "not found in colData")
  expect_error(phaseConfounding(se, ~1, nPerm = 9), "no terms to test")
  expect_error(phaseConfounding(matrix(1:4, 2), ~condition), "SummarizedExperiment")

  # no phase column and no way to compute one
  bare <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logExpr = matrix(
      stats::rnorm(20), nrow = 10,
      dimnames = list(paste0("g", 1:10), c("s1", "s2"))
    ))
  )
  expect_error(phaseConfounding(bare, ~condition), "requires both 'method' and 'tissue'")
})

test_that("a term with a single level yields NA rather than aborting", {
  set.seed(1)
  se <- makeConfoundingSE(
    stats::runif(20, 0, 24),
    condition = rep("only", 20)
  )
  res <- phaseConfounding(se, ~condition, nPerm = 99, seed = 1)
  expect_true(is.na(res$statistic))
  expect_true(is.na(res$p_value))
})

test_that("phase is estimated on the fly when se has no phase column", {
  set.seed(1)
  # Deliberately random genes: estimatePhase()'s own contract is to flag rather than
  # error, so this exercises the integration path and the resulting all-dropped guard.
  mat <- matrix(
    stats::rnorm(100 * 4), nrow = 100,
    dimnames = list(paste0("g", seq_len(100)), paste0("s", seq_len(4)))
  )
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  SummarizedExperiment::colData(se)$condition <- rep(c("a", "b"), each = 2)

  expect_error(
    phaseConfounding(se, ~condition, method = "timetable", tissue = "Liv", nPerm = 9),
    "Fewer than 3 usable phase values"
  )
})
