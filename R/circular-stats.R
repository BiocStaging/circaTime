# Circular statistics primitives. Per AGENTS.md 7.1, all circular maths lives here;
# every other file must call into these functions rather than compute on raw phase
# values directly.

#' Wrap a phase value onto \verb{[0, period)}
#'
#' @param x Numeric vector of phase values, in the same units as `period`.
#' @param period Length of the cycle (default 24, for hours).
#' @return Numeric vector, each element wrapped into `[0, period)`.
#' @examples
#' wrapPhase(c(-1, 0, 23.5, 24, 25.5))
#' @export
wrapPhase <- function(x, period = 24) {
  x %% period
}

#' Signed circular difference between two phases
#'
#' @param pred,truth Numeric vectors of phase values (same units, same `period`).
#' @param period Length of the cycle (default 24).
#' @return Numeric vector of signed differences in `(-period/2, period/2]`, i.e. the
#'   shortest signed rotation from `truth` to `pred`.
#' @examples
#' circularDiff(0.5, 23.5) # 1, not -23
#' @export
circularDiff <- function(pred, truth, period = 24) {
  d <- (pred - truth) %% period
  ifelse(d > period / 2, d - period, d)
}

#' Circular mean absolute error
#'
#' The primary error metric for this package (AGENTS.md 7.1). Computed as
#' `pmin(abs(d), period - abs(d))` where `d = pred - truth`, so that e.g. a
#' prediction of 23.5 against a truth of 0.5 counts as 1 hour of error, not 23.
#'
#' @param pred,truth Numeric vectors of phase values in `[0, period)`. `NA` entries
#'   are dropped pairwise before averaging.
#' @param period Length of the cycle (default 24, for hours).
#' @param na.rm Logical; if `FALSE` (default `TRUE`), `NA` in either input makes the
#'   result `NA` instead of being dropped.
#' @return A single numeric value: the mean circular absolute error.
#' @examples
#' circularMAE(c(0.5, 12), c(23.5, 0))
#' @export
circularMAE <- function(pred, truth, period = 24, na.rm = TRUE) {
  stopifnot(length(pred) == length(truth))
  d <- pred - truth
  err <- pmin(abs(d), period - abs(d))
  mean(err, na.rm = na.rm)
}

#' Per-sample circular absolute error
#'
#' Vectorised version of the error term inside [circularMAE()], useful for residual
#' analysis (e.g. `benchmarkPhase()` per-sample output).
#'
#' @inheritParams circularMAE
#' @return Numeric vector of per-sample circular absolute errors.
#' @examples
#' circularError(c(0.5, 12), c(23.5, 0))
#' @export
circularError <- function(pred, truth, period = 24) {
  stopifnot(length(pred) == length(truth))
  d <- pred - truth
  pmin(abs(d), period - abs(d))
}

#' Circular mean direction
#'
#' @param x Numeric vector of phase values.
#' @param period Length of the cycle (default 24).
#' @param na.rm Logical; drop `NA` values before averaging (default `TRUE`).
#' @return A single numeric value in `[0, period)`: the mean direction of `x`.
#' @examples
#' circularMean(c(23, 1)) # ~0, not 12
#' @export
circularMean <- function(x, period = 24, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  theta <- x / period * 2 * pi
  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta)))
  wrapPhase(mean_theta / (2 * pi) * period, period)
}

#' Weighted circular mean direction
#'
#' The weighted generalisation of [circularMean()]: each phase is converted to a
#' unit vector, the vectors are summed with the supplied weights, and the
#' direction of the resultant is returned. Used by [estimatePhase()]'s
#' `method = "consensus"` to combine several engines' phase estimates (AGENTS.md
#' Section 6 Step 3.1). `circularWeightedMean(x, rep(1, length(x)))` is identical
#' to `circularMean(x)`.
#'
#' The resultant vector's *length* is also returned (as attribute
#' `"concentration"`, in `[0, 1]`): 1 means every input phase agreed exactly, 0
#' means they cancelled out completely. This is the standard circular measure of
#' agreement, and for a consensus estimate it is a directly interpretable
#' "did the engines agree?" signal rather than a hidden intermediate.
#'
#' @param x Numeric vector of phase values.
#' @param w Numeric vector of non-negative weights, same length as `x`. Weights
#'   are normalised internally, so only their relative sizes matter.
#' @param period Length of the cycle (default 24).
#' @param na.rm Logical; drop pairs where `x` is `NA` before averaging
#'   (default `TRUE`). Weights of dropped entries are dropped with them, and the
#'   remainder renormalised.
#' @return A single numeric value in `[0, period)`, with attribute
#'   `"concentration"` (numeric, `[0, 1]`). `NA_real_` (concentration `NA`) if no
#'   usable input remains or all weights are zero.
#' @examples
#' # two engines disagreeing by 2h, the second trusted twice as much
#' circularWeightedMean(c(6, 8), c(1, 2))
#' # perfect agreement -> concentration 1
#' attr(circularWeightedMean(c(6, 6), c(1, 1)), "concentration")
#' @export
circularWeightedMean <- function(x, w, period = 24, na.rm = TRUE) {
  stopifnot(length(x) == length(w))
  if (any(w < 0, na.rm = TRUE)) stop("weights must be non-negative", call. = FALSE)
  if (na.rm) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }
  naOut <- structure(NA_real_, concentration = NA_real_)
  if (length(x) == 0 || sum(w) == 0) return(naOut)

  w <- w / sum(w)
  theta <- x / period * 2 * pi
  sumSin <- sum(w * sin(theta))
  sumCos <- sum(w * cos(theta))
  structure(
    wrapPhase(atan2(sumSin, sumCos) / (2 * pi) * period, period),
    concentration = sqrt(sumSin^2 + sumCos^2)
  )
}

#' Circular-circular correlation coefficient (Jammalamadaka & SenGupta, 2001)
#'
#' The circular analogue of Pearson correlation. Must be used instead of `cor()`
#' whenever either input is a phase (AGENTS.md 7.1).
#'
#' @section Known instability on uniformly-sampled time courses:
#' This coefficient centres each variable on its own circular mean before
#' comparing them. When one variable's values are (near-)evenly spaced around the
#' *entire* circle - e.g. ground truth from a time course sampled every 2h across
#' a full 24h day, such as GSE54650 - that variable's resultant vector length is
#' ~0, its circular mean becomes numerically arbitrary (driven by floating-point
#' noise, not signal), and the resulting correlation can come out strongly
#' negative *even for near-perfect predictions*. This is a property of the
#' Jammalamadaka-SenGupta coefficient itself, not a bug: verified directly on
#' GSE54650-shaped ground truth (24 evenly-spaced timepoints), a prediction set
#' with 0.2h mean circular error still yields a coefficient around -0.46.
#' [circularMAE()] does not share this failure mode and should be treated as the
#' primary accuracy metric for evenly-sampled benchmark data; treat `circularCor`
#' as a secondary/diagnostic number on such data, not a reliability gate.
#'
#' @param x,y Numeric vectors of phase values (same units, same `period`).
#' @param period Length of the cycle (default 24).
#' @return A single numeric value in `[-1, 1]`.
#' @examples
#' circularCor(c(0, 6, 12, 18), c(1, 7, 13, 19))
#' @export
circularCor <- function(x, y, period = 24) {
  stopifnot(length(x) == length(y))
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 2) return(NA_real_)

  thetaX <- x / period * 2 * pi
  thetaY <- y / period * 2 * pi
  xBar <- atan2(mean(sin(thetaX)), mean(cos(thetaX)))
  yBar <- atan2(mean(sin(thetaY)), mean(cos(thetaY)))

  sinX <- sin(thetaX - xBar)
  sinY <- sin(thetaY - yBar)

  num <- sum(sinX * sinY)
  den <- sqrt(sum(sinX^2) * sum(sinY^2))
  if (den == 0) return(NA_real_)
  num / den
}

#' Circular quantile interval around a set of phase draws
#'
#' Used by [bootstrapPhaseCI()] to turn a sample's bootstrap replicate
#' predictions into a `phase_lo`/`phase_hi` interval. Centres the draws on
#' their circular mean, takes the empirical `(1-level)/2` and `1-(1-level)/2`
#' quantiles of the *signed* circular deviation from that centre (via
#' [circularDiff()], so a deviation is always in `(-period/2, period/2]` and
#' wraparound never inflates the spread), then maps those quantile deviations
#' back onto the circle. This assumes the draws are unimodal around their
#' circular mean, which holds for a bootstrap distribution of predictions
#' clustered near the true phase, but would not hold for, say, a bimodal
#' rotation/reflection-ambiguous set (AGENTS.md 7.2) - resolve that ambiguity
#' before calling this.
#'
#' @param x Numeric vector of phase draws, in `[0, period)`; `NA` entries are
#'   dropped.
#' @param level Nominal interval coverage, in `(0, 1)` (default 0.8).
#' @param period Length of the cycle (default 24).
#' @return A list with `center` (circular mean of `x`), `lo`, `hi` (interval
#'   bounds in `[0, period)`). All `NA_real_` if fewer than 2 non-`NA` draws
#'   remain.
#' @examples
#' circularQuantileInterval(c(23, 0, 1, 23.5, 0.5))
#' @export
circularQuantileInterval <- function(x, level = 0.8, period = 24) {
  stopifnot(level > 0, level < 1)
  x <- x[!is.na(x)]
  if (length(x) < 2) {
    return(list(center = NA_real_, lo = NA_real_, hi = NA_real_))
  }

  center <- circularMean(x, period = period)
  dev <- circularDiff(x, center, period = period)
  alpha <- (1 - level) / 2
  devLo <- stats::quantile(dev, probs = alpha, type = 7, names = FALSE)
  devHi <- stats::quantile(dev, probs = 1 - alpha, type = 7, names = FALSE)

  list(
    center = center,
    lo = wrapPhase(center + devLo, period = period),
    hi = wrapPhase(center + devHi, period = period)
  )
}

#' Mean resultant length of a set of phases (circular concentration)
#'
#' The standard circular-statistics measure of how tightly a sample of phases
#' clusters around its own circular mean: 1 for perfect agreement, 0 for phases
#' spread uniformly around the whole circle. Equivalent to `1 - circular
#' variance`. This is the primitive [circularWeightedMean()]'s `"concentration"`
#' attribute already computes internally for a *weighted combination of engine
#' estimates*; this function computes the same quantity directly for any raw
#' sample of phases (or phase deviations), which is what Phase 3 Step 3.3's
#' Bayesian posterior (see [bayesianPhasePosterior()]) needs applied to a
#' reference model's residual pool.
#'
#' @param x Numeric vector of phase values (or circular deviations, e.g. from
#'   [circularDiff()]).
#' @param period Length of the cycle (default 24).
#' @param na.rm Logical; drop `NA` values before computing (default `TRUE`).
#' @return A single numeric value in `[0, 1]`. `NA_real_` if no usable input
#'   remains.
#' @examples
#' circularConcentration(c(0, 0, 0)) # 1: no spread at all
#' circularConcentration(c(0, 6, 12, 18)) # ~0: spread evenly around the circle
#' @export
circularConcentration <- function(x, period = 24, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  theta <- x / period * 2 * pi
  sqrt(mean(cos(theta))^2 + mean(sin(theta))^2)
}

# Converts a mean resultant length R to the von Mises concentration parameter
# kappa via the standard closed-form approximation (Fisher 1993, "Statistical
# Analysis of Circular Data", Eq. 4.40; equivalently Mardia & Jupp 2000), rather
# than numerically inverting A(kappa) = R with a root-finder - accurate to
# ~O(1e-3) across the whole R range and avoids pulling in a circular-stats
# package (AGENTS.md: no heavy dependency for something this size).
.vonMisesKappaFromR <- function(R) {
  ifelse(R < 0.53, 2 * R + R^3 + 5 * R^5 / 6,
  ifelse(R < 0.85, -0.4 + 1.39 * R + 0.43 / (1 - R),
         1 / (R^3 - 4 * R^2 + 3 * R)))
}

# Quantile of the SIGNED deviation from a von Mises(0, kappa) distribution, on a
# fine grid via numerical integration of the density - there is no closed-form
# von Mises quantile function. kappa = 0 (uniform) is handled directly since the
# grid approach is numerically unstable there. n controls grid resolution, not
# precision-critical for a QC/uncertainty display, so a moderate default is fine.
.vonMisesQuantile <- function(p, kappa, period = 24, n = 2000) {
  if (kappa <= 0) return((p - 0.5) * period) # uniform: every deviation equally likely
  half <- period / 2
  dev <- seq(-half, half, length.out = n)
  theta <- dev / period * 2 * pi
  dens <- exp(kappa * (cos(theta) - 1)) # shift for numerical stability; shape unchanged
  cdf <- cumsum(dens)
  cdf <- cdf / cdf[length(cdf)]
  # At very high kappa the density concentrates into a handful of grid cells, so
  # cdf plateaus at 0/1 and repeats - approx() requires strictly increasing x.
  # Keep only the first occurrence of each distinct cdf value rather than passing
  # ties through (harmless numerically, but approx() warns on every call).
  keep <- !duplicated(cdf)
  stats::approx(cdf[keep], dev[keep], xout = p, rule = 2)$y
}

#' Bayesian (von Mises) credible interval around a point phase estimate
#'
#' The parametric counterpart to [circularQuantileInterval()]: instead of taking
#' empirical quantiles of a finite set of bootstrap draws, this fits a von Mises
#' distribution - the circular analogue of the normal distribution, and the
#' standard choice for Bayesian circular models (Fisher 1993) - to the
#' concentration already measured from real data, and returns the analytic
#' credible interval of that distribution centred on `mu`. Because a von Mises
#' distribution is unimodal and symmetric, the equal-tailed interval computed
#' here is also the shortest (highest-density) interval at the same level - no
#' separate HDI computation is needed.
#'
#' @param mu Single numeric point phase estimate (or a vector, recycled per
#'   sample), in `[0, period)`.
#' @param kappa Von Mises concentration parameter (non-negative). See
#'   [circularConcentration()] + the internal R-to-kappa conversion, or supply a
#'   value fit elsewhere. `kappa = 0` is the uniform distribution (maximally
#'   uninformative posterior - the interval spans the whole circle).
#' @param level Nominal credible level (default 0.8, matching
#'   [circularQuantileInterval()]'s default so the two are directly comparable).
#' @param period Length of the cycle (default 24).
#' @return A list with `lo`, `hi` (interval bounds in `[0, period)`, same length
#'   as `mu`).
#' @examples
#' circularCredibleInterval(mu = 12, kappa = 5, level = 0.8)
#' @export
circularCredibleInterval <- function(mu, kappa, level = 0.8, period = 24) {
  stopifnot(level > 0, level < 1, all(kappa >= 0))
  alpha <- (1 - level) / 2
  devLo <- .vonMisesQuantile(alpha, kappa, period = period)
  devHi <- .vonMisesQuantile(1 - alpha, kappa, period = period)
  list(lo = wrapPhase(mu + devLo, period = period), hi = wrapPhase(mu + devHi, period = period))
}

#' Does a circular interval cover a truth value?
#'
#' Tests whether `truth` lies on the arc walked from `lo` to `hi` in the
#' increasing direction (wrapping past `period` if needed) - i.e. the same
#' direction [circularQuantileInterval()] builds its interval in. A
#' zero-width interval (`lo == hi`) only covers an exact match, not the whole
#' circle.
#'
#' @param truth Numeric vector of ground-truth phase values.
#' @param lo,hi Numeric vectors of interval bounds, same length as `truth`
#'   (as returned by [circularQuantileInterval()] per sample).
#' @param period Length of the cycle (default 24).
#' @return Logical vector, `NA` wherever `truth`, `lo`, or `hi` is `NA`.
#' @examples
#' circularIntervalCovers(23.9, 23, 1) # TRUE: 23.9 is on the 23 -> 1 arc
#' circularIntervalCovers(12, 23, 1)   # FALSE
#' @export
circularIntervalCovers <- function(truth, lo, hi, period = 24) {
  stopifnot(length(truth) == length(lo), length(lo) == length(hi))
  arcLen <- (hi - lo) %% period
  pos <- (truth - lo) %% period
  covers <- pos <= arcLen
  covers[is.na(truth) | is.na(lo) | is.na(hi)] <- NA
  covers
}

#' Align unsupervised phase estimates to a reference by rotation and reflection
#'
#' Unsupervised engines (e.g. CYCLOPS-family, CHIRAL) recover sample ordering around
#' the circle only up to an arbitrary rotation and reflection (AGENTS.md 7.2). This
#' finds the rotation offset and reflection (out of the two possible orientations)
#' that best aligns `est` to `ref`, then returns the aligned estimate.
#'
#' The rotation for each orientation is estimated as the circular mean of the
#' pairwise differences to `ref`; this is the standard closed-form alignment used in
#' the CYCLOPS/CHIRAL literature and is a close practical approximation to the
#' rotation that minimises circular MAE.
#'
#' @param est Numeric vector of unsupervised phase estimates, in `[0, period)`.
#' @param ref Numeric vector of reference/ground-truth phases, same length as `est`.
#'   `NA` entries are excluded from fitting the alignment but the alignment is still
#'   applied to every element of `est`.
#' @param period Length of the cycle (default 24).
#' @param allow_reflection Logical; if `TRUE` (default), also try reversing the
#'   direction of `est` and keep whichever orientation gives lower circular MAE
#'   against `ref`.
#' @return A list with elements `aligned` (numeric vector, same length as `est`),
#'   `rotation` (the applied additive offset), `reflected` (logical), and `mae` (the
#'   circular MAE of `aligned` vs `ref`, over the non-`NA` fitting subset).
#' @examples
#' truth <- c(0, 4, 8, 12, 16, 20)
#' est <- wrapPhase(-truth + 5) # reflected and rotated
#' circularAlign(est, truth)$rotation
#' @export
circularAlign <- function(est, ref, period = 24, allow_reflection = TRUE) {
  stopifnot(length(est) == length(ref))
  keep <- !is.na(est) & !is.na(ref)
  if (sum(keep) < 2) {
    stop("circularAlign() needs at least 2 non-NA (est, ref) pairs to fit an alignment")
  }

  fitOrientation <- function(candidate) {
    rotation <- circularMean(ref[keep] - candidate[keep], period = period)
    aligned <- wrapPhase(candidate + rotation, period = period)
    mae <- circularMAE(aligned[keep], ref[keep], period = period)
    list(aligned = aligned, rotation = rotation, mae = mae)
  }

  fwd <- fitOrientation(est)
  best <- c(fwd, reflected = FALSE)

  if (allow_reflection) {
    rev <- fitOrientation(wrapPhase(-est, period = period))
    if (rev$mae < best$mae) {
      best <- c(rev, reflected = TRUE)
    }
  }

  list(aligned = best$aligned, rotation = best$rotation,
       reflected = best$reflected, mae = best$mae)
}
