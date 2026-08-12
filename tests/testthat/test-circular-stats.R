# Tests for R/circular-stats.R. The wraparound-at-the-0/24-boundary tests here are
# an explicit Phase 1 gate condition (AGENTS.md Section 9).

test_that("wrapPhase wraps onto [0, period)", {
  expect_equal(wrapPhase(-1), 23)
  expect_equal(wrapPhase(24), 0)
  expect_equal(wrapPhase(25.5), 1.5)
  expect_equal(wrapPhase(c(-1, 0, 23.5, 24, 25.5)), c(23, 0, 23.5, 0, 1.5))
})

test_that("circularDiff treats 23.5 and 0.5 as 1 hour apart, not 23", {
  expect_equal(circularDiff(0.5, 23.5), 1)
  expect_equal(circularDiff(23.5, 0.5), -1)
})

test_that("circularDiff is antisymmetric and zero for identical inputs", {
  expect_equal(circularDiff(5, 5), 0)
  expect_equal(circularDiff(10, 3), -circularDiff(3, 10))
})

test_that("circularDiff stays within (-12, 12] at every point around the circle", {
  truth <- seq(0, 23.9, by = 0.1)
  pred <- wrapPhase(truth + 0.5)
  d <- circularDiff(pred, truth)
  expect_true(all(d > -12 & d <= 12))
})

test_that("circularMAE handles the 0/24 wraparound correctly", {
  # naive |pred - truth| would give 23h error here; true circular error is 1h
  expect_equal(circularMAE(0.5, 23.5), 1)
  expect_equal(circularMAE(23.5, 0.5), 1)
  # pair 1: 0.5 vs 23.5 -> 1h error; pair 2: 1 vs 23 -> 2h error; mean = 1.5
  expect_equal(circularMAE(c(0.5, 1), c(23.5, 23)), 1.5)
})

test_that("circularMAE is 0 for identical vectors and symmetric bounded near 12", {
  x <- c(0, 6, 12, 18, 23.9)
  expect_equal(circularMAE(x, x), 0)
  expect_equal(circularMAE(0, 12), 12) # antipodal: max possible error
})

test_that("circularMAE respects na.rm", {
  pred <- c(0.5, NA, 5)
  truth <- c(23.5, 3, 5)
  expect_equal(circularMAE(pred, truth, na.rm = TRUE), mean(c(1, 0)))
  expect_true(is.na(circularMAE(pred, truth, na.rm = FALSE)))
})

test_that("circularError is the per-sample vector version of circularMAE", {
  pred <- c(0.5, 5)
  truth <- c(23.5, 5)
  expect_equal(circularError(pred, truth), c(1, 0))
  expect_equal(mean(circularError(pred, truth)), circularMAE(pred, truth))
})

test_that("circularMean wraps around correctly (23 and 1 average to ~0, not 12)", {
  expect_equal(round(circularMean(c(23, 1)), 6), 0)
  expect_equal(round(circularMean(c(0, 0, 0)), 6), 0)
})

test_that("circularMean of evenly spaced points around the whole circle is undefined-ish (near 0 magnitude) but doesn't error", {
  x <- seq(0, 23, length.out = 12)
  expect_true(is.finite(circularMean(x)))
})

test_that("circularMean drops NA by default and returns NA if nothing left", {
  expect_equal(round(circularMean(c(23, 1, NA)), 6), 0)
  expect_true(is.na(circularMean(c(NA_real_, NA_real_))))
})

test_that("circularWeightedMean with equal weights matches circularMean", {
  x <- c(1, 5, 23, 14)
  expect_equal(
    as.numeric(circularWeightedMean(x, rep(1, 4))),
    circularMean(x)
  )
  # weight scale must not matter, only relative size
  expect_equal(
    as.numeric(circularWeightedMean(x, rep(7, 4))),
    as.numeric(circularWeightedMean(x, rep(1, 4)))
  )
})

test_that("circularWeightedMean wraps around the midnight boundary", {
  # 23 and 1 should average to ~0, never 12
  expect_equal(round(as.numeric(circularWeightedMean(c(23, 1), c(1, 1))), 6), 0)
})

test_that("circularWeightedMean pulls toward the more heavily weighted phase", {
  # 6 and 8, second trusted twice as much -> should land past the midpoint (7)
  out <- as.numeric(circularWeightedMean(c(6, 8), c(1, 2)))
  expect_gt(out, 7)
  expect_lt(out, 8)
  # a zero weight removes that input entirely
  expect_equal(as.numeric(circularWeightedMean(c(6, 18), c(1, 0))), 6)
})

test_that("circularWeightedMean reports concentration as agreement in [0, 1]", {
  # identical inputs agree perfectly
  expect_equal(attr(circularWeightedMean(c(6, 6), c(1, 1)), "concentration"), 1)
  # antipodal inputs with equal weight cancel out completely
  expect_equal(
    round(attr(circularWeightedMean(c(0, 12), c(1, 1)), "concentration"), 12), 0
  )
  # partial disagreement lands strictly between
  conc <- attr(circularWeightedMean(c(4, 8), c(1, 1)), "concentration")
  expect_gt(conc, 0)
  expect_lt(conc, 1)
})

test_that("circularWeightedMean drops NA phases and renormalises the rest", {
  # the NA entry must not drag the result toward 0 or produce NA
  expect_equal(as.numeric(circularWeightedMean(c(6, NA), c(1, 1))), 6)
  out <- circularWeightedMean(c(NA_real_, NA_real_), c(1, 1))
  expect_true(is.na(out))
  expect_true(is.na(attr(out, "concentration")))
})

test_that("circularWeightedMean returns NA when all weights are zero, and rejects negatives", {
  out <- circularWeightedMean(c(6, 8), c(0, 0))
  expect_true(is.na(out))
  expect_error(circularWeightedMean(c(6, 8), c(1, -1)), "non-negative")
  expect_error(circularWeightedMean(c(6, 8), c(1, 1, 1)))
})

test_that("circularCor is ~1 for perfectly matched circular data", {
  x <- c(0, 4, 8, 12, 16, 20)
  y <- x
  expect_equal(round(circularCor(x, y), 6), 1)
})

test_that("circularCor is ~1 for a pure rotation (shift) of the same pattern", {
  # NB: points evenly spaced around the *whole* circle (e.g. c(0,4,8,12,16,20))
  # have a degenerate (zero-length) resultant vector, which makes the circular
  # mean - and hence this test - numerically unstable. Use an irregularly spaced,
  # clustered set instead, as real phase data would look like.
  x <- c(1, 2, 3, 9, 10, 14)
  y <- wrapPhase(x + 6)
  expect_gt(circularCor(x, y), 0.99)
})

test_that("circularCor is ~-1 for a pure reflection", {
  x <- c(1, 2, 3, 9, 10, 14)
  y <- wrapPhase(-x)
  expect_lt(circularCor(x, y), -0.99)
})

test_that("circularCor handles NA pairwise and short input", {
  x <- c(0, 4, 8, NA, 16, 20)
  y <- c(0, 4, 8, 12, 16, 20)
  expect_false(is.na(circularCor(x, y)))
  expect_true(is.na(circularCor(c(1), c(1))))
})

test_that("circularAlign recovers a known rotation with no reflection", {
  truth <- c(0, 4, 8, 12, 16, 20)
  offset <- 5
  est <- wrapPhase(truth + offset)
  # est = truth + offset, so to align est back to truth we need rotation ~ -offset
  res <- circularAlign(est, truth)
  expect_false(res$reflected)
  expect_lt(res$mae, 1e-6)
  expect_equal(round(wrapPhase(res$rotation + offset), 4), 0)
})

test_that("circularAlign recovers a known rotation + reflection", {
  truth <- c(0, 4, 8, 12, 16, 20)
  est <- wrapPhase(-truth + 5) # reflected direction, then rotated
  res <- circularAlign(est, truth)
  expect_true(res$reflected)
  expect_lt(res$mae, 1e-6)
})

test_that("circularAlign without allow_reflection cannot fix a reflected input", {
  truth <- c(0, 4, 8, 12, 16, 20)
  est <- wrapPhase(-truth + 5)
  res <- circularAlign(est, truth, allow_reflection = FALSE)
  expect_false(res$reflected)
  expect_gt(res$mae, 1) # cannot be aligned well without reflection
})

test_that("circularAlign errors with fewer than 2 usable pairs", {
  expect_error(circularAlign(c(1, NA), c(NA, 2)))
})

test_that("circularQuantileInterval brackets a tight cluster near its center", {
  x <- c(11.5, 11.8, 12, 12.2, 12.5)
  ci <- circularQuantileInterval(x, level = 0.8)
  expect_equal(round(ci$center, 1), 12)
  expect_lt(ci$lo, ci$center)
  expect_gt(ci$hi, ci$center)
})

test_that("circularQuantileInterval handles the 0/24 wraparound without inflating width", {
  # draws clustered right across midnight; a naive (non-circular) quantile
  # interval on raw values would come out enormous (spanning ~23h)
  x <- c(23.5, 23.8, 0, 0.2, 0.5)
  ci <- circularQuantileInterval(x, level = 0.8)
  expect_equal(round(ci$center, 1), 0)
  width <- (ci$hi - ci$lo) %% 24
  expect_lt(width, 3)
})

test_that("circularQuantileInterval returns NA with fewer than 2 non-NA draws", {
  ci <- circularQuantileInterval(c(5, NA))
  expect_true(is.na(ci$center))
  expect_true(is.na(ci$lo))
  expect_true(is.na(ci$hi))
})

test_that("circularQuantileInterval errors on an out-of-range level", {
  expect_error(circularQuantileInterval(c(1, 2, 3), level = 1))
  expect_error(circularQuantileInterval(c(1, 2, 3), level = 0))
})

test_that("circularIntervalCovers handles the 0/24 wraparound arc correctly", {
  # arc from 23 to 1, wrapping through midnight
  expect_true(circularIntervalCovers(23.9, 23, 1))
  expect_true(circularIntervalCovers(0, 23, 1))
  expect_true(circularIntervalCovers(0.9, 23, 1))
  expect_false(circularIntervalCovers(12, 23, 1))
  expect_false(circularIntervalCovers(22, 23, 1))
})

test_that("circularIntervalCovers handles a non-wrapping arc correctly", {
  expect_true(circularIntervalCovers(10, 8, 12))
  expect_false(circularIntervalCovers(14, 8, 12))
})

test_that("circularIntervalCovers is TRUE only at the exact point for a zero-width interval", {
  expect_true(circularIntervalCovers(5, 5, 5))
  expect_false(circularIntervalCovers(5.1, 5, 5))
})

test_that("circularIntervalCovers propagates NA", {
  expect_true(is.na(circularIntervalCovers(NA, 1, 2)))
  expect_true(is.na(circularIntervalCovers(1, NA, 2)))
})

# circularConcentration() / circularCredibleInterval() - Phase 3 Step 3.3.
test_that("circularConcentration is 1 for identical values and ~0 for a uniform spread", {
  expect_equal(circularConcentration(c(5, 5, 5)), 1)
  expect_lt(circularConcentration(c(0, 6, 12, 18)), 0.01)
})

test_that("circularConcentration decreases as spread increases", {
  set.seed(1)
  tight <- circularConcentration(rnorm(500, mean = 0, sd = 1))
  loose <- circularConcentration(rnorm(500, mean = 0, sd = 6))
  expect_gt(tight, loose)
})

test_that("circularConcentration drops NA and returns NA_real_ for all-NA/empty input", {
  expect_equal(circularConcentration(c(1, NA, 1, NA)), circularConcentration(c(1, 1)))
  expect_true(is.na(circularConcentration(c(NA, NA))))
  expect_true(is.na(circularConcentration(numeric(0))))
})

test_that("circularCredibleInterval widens as kappa decreases and spans the whole circle at kappa = 0", {
  narrow <- circularCredibleInterval(mu = 12, kappa = 20, level = 0.8)
  wide <- circularCredibleInterval(mu = 12, kappa = 1, level = 0.8)
  # mu = 12 keeps both intervals away from the 0/24 boundary, so raw hi - lo is a
  # safe width measure here (circularDiff's shortest-rotation fold breaks once an
  # interval exceeds half the circle, which the wide/kappa=1 case does).
  width <- function(ci) ci$hi - ci$lo
  expect_lt(width(narrow), width(wide))

  # width > half the circle here (19.2h), so circularDiff's shortest-rotation
  # convention would fold it to its -4.8h complement - use raw hi - lo instead,
  # which is safe since this interval (2.4-21.6) does not itself wrap past 0/24.
  uniform <- circularCredibleInterval(mu = 12, kappa = 0, level = 0.8)
  expect_equal(uniform$hi - uniform$lo, 24 * 0.8, tolerance = 1e-6)
})

test_that("circularCredibleInterval is centred on mu (symmetric deviation either side)", {
  ci <- circularCredibleInterval(mu = 12, kappa = 5, level = 0.8)
  loDev <- circularDiff(12, ci$lo)
  hiDev <- circularDiff(ci$hi, 12)
  # grid-based quantile (n = 2000 over a 24h period -> ~0.012h/cell) introduces
  # sub-grid-cell rounding, not a real asymmetry
  expect_equal(loDev, hiDev, tolerance = 0.02)
})

test_that("circularCredibleInterval vectorises over mu with a shared kappa", {
  ci <- circularCredibleInterval(mu = c(0, 12), kappa = 5, level = 0.8)
  expect_length(ci$lo, 2)
  expect_length(ci$hi, 2)
  expect_equal(circularDiff(ci$hi[1], ci$lo[1]), circularDiff(ci$hi[2], ci$lo[2]), tolerance = 0.02)
})

test_that("circularCredibleInterval errors on invalid level or negative kappa", {
  expect_error(circularCredibleInterval(12, 5, level = 0), "level")
  expect_error(circularCredibleInterval(12, -1, level = 0.8), "kappa")
})
