# phaseConfounding() - the AGENTS.md Section 6 diagnostic aimed at users who are not
# circadian researchers: "could unmodelled collection time be confounding my
# differential expression design?"
#
# The question is deliberately narrow. Circadian phase confounds a `~ condition`
# contrast only when phase is UNBALANCED across the levels of that term: if collection
# time is spread evenly within every group, rhythmic genes add variance but do not bias
# the condition coefficient. If the cases were collected in the morning and the controls
# in the afternoon, every rhythmic gene shows a spurious condition effect. So what is
# tested here is the association between (inferred) phase and each design term - not
# per-gene rhythmicity, which is `limorhyde`/`dryR` territory and an explicit non-goal
# (AGENTS.md Section 3).
#
# Phase is circular, so every test here is circular (AGENTS.md 7.1): a Watson-Williams
# statistic for categorical terms, a Mardia/Johnson-Wehrly circular-linear correlation
# for numeric ones. Both get their null distribution by permutation rather than from the
# usual parametric approximations, which assume von Mises phases and (for
# Watson-Williams) equal concentration across groups - assumptions real phase estimates
# routinely violate.

# Watson-Williams between-group statistic, in its F form. Under label permutation the
# pooled resultant R, N and k are all invariant, so the ranking depends only on
# sum(Ri) - the usual (1 + 3/(8*kappa)) correction is likewise permutation-invariant
# and is therefore omitted rather than estimated.
.watsonWilliamsStat <- function(theta, grp) {
  n <- length(theta)
  groups <- unique(grp)
  k <- length(groups)
  if (k < 2) return(NA_real_)

  resultant <- function(th) sqrt(sum(cos(th))^2 + sum(sin(th))^2)
  sumRi <- sum(vapply(groups, function(g) resultant(theta[grp == g]), numeric(1)))
  R <- resultant(theta)

  denom <- n - sumRi
  # Every observation identical within its group: no residual dispersion for the
  # statistic to scale against. Not a failure of the data, just undefined here.
  if (denom <= .Machine$double.eps) return(NA_real_)
  (n - k) * (sumRi - R) / ((k - 1) * denom)
}

# Circular-linear correlation (Mardia 1976; Johnson & Wehrly 1977): how much of a
# linear covariate's variation is explained by the (cos, sin) pair of the phase.
.circularLinearStat <- function(theta, x) {
  cosT <- cos(theta); sinT <- sin(theta)
  if (stats::sd(x) == 0 || stats::sd(cosT) == 0 || stats::sd(sinT) == 0) {
    return(NA_real_)
  }
  rxc <- stats::cor(x, cosT)
  rxs <- stats::cor(x, sinT)
  rcs <- stats::cor(cosT, sinT)
  denom <- 1 - rcs^2
  if (denom <= .Machine$double.eps) return(NA_real_)
  r2 <- (rxc^2 + rxs^2 - 2 * rxc * rxs * rcs) / denom
  sqrt(max(r2, 0))
}

#' Is circadian phase confounding a differential expression design?
#'
#' The diagnostic for users who are not circadian researchers (AGENTS.md Section 6):
#' given a design such as `~ condition`, is the samples' circadian phase distributed
#' unevenly across that design's terms? If it is, rhythmic genes will show apparent
#' condition effects that are really collection-time effects.
#'
#' Phase is only a confounder when it is *unbalanced* across a term's levels. Phase
#' that varies randomly within every group inflates residual variance but does not
#' bias the term's coefficient, so a non-significant result here is genuinely
#' reassuring about bias (though not about power).
#'
#' @section What is tested, and what is not:
#' This tests the association between phase and each design term. It does **not**
#' test which genes are rhythmic - rhythmicity testing on time-labelled data is a
#' solved problem and an explicit non-goal (AGENTS.md Section 3); hand the phased
#' object to `limorhyde`/`dryR` for that. It also does not adjust or correct the
#' design; it reports the imbalance so the analyst can decide whether to add phase
#' as a covariate, which is a modelling decision, not something to apply silently.
#'
#' @section Statistics:
#' Both tests are circular (AGENTS.md 7.1) and both take their null distribution
#' from permuting phase against the term, which avoids the von Mises and
#' equal-concentration assumptions the parametric versions rely on:
#' * **Categorical terms** (factor, character, logical): the Watson-Williams
#'   between-group statistic in its F form.
#' * **Numeric terms**: the Mardia/Johnson-Wehrly circular-linear correlation.
#'
#' P-values use the add-one estimator `(1 + #{perm >= observed}) / (nPerm + 1)`, so
#' the smallest reportable value is `1 / (nPerm + 1)` and a p-value is never
#' reported as exactly zero.
#'
#' @param se A `SummarizedExperiment`.
#' @param design A one-sided formula naming `colData(se)` columns, e.g.
#'   `~ condition` or `~ condition + age`. Each term is tested separately.
#'   Interaction terms (`a:b`) are skipped with a warning - a term-by-term
#'   imbalance report is what is actionable here, and a full multivariable
#'   circular regression is out of scope.
#' @param phase_col Name of the `colData(se)` column holding phase. If absent,
#'   [estimatePhase()] is called first, which requires `method` and `tissue`.
#' @param flag_col Name of the `colData(se)` column holding [estimatePhase()]'s
#'   quality flag. When present, only samples flagged `"ok"` are tested -
#'   `"low_amplitude"`/`"gene_coverage"`/`"failed"` phases are unreliable and would
#'   only add noise. Set to `NULL` to use every non-`NA` phase regardless.
#' @param method,tissue,species Passed to [estimatePhase()], and used only when
#'   `se` does not already carry `phase_col`.
#' @param nPerm Number of permutations (default 1000).
#' @param period Length of the cycle (default 24).
#' @param seed Optional integer seed; the global RNG state is left unchanged either
#'   way (AGENTS.md Section 10).
#' @return A `data.frame` with one row per tested term: `term`, `type`
#'   (`"categorical"`/`"numeric"`), `n` (samples used), `n_groups` (`NA` for numeric
#'   terms), `statistic` (Watson-Williams F, or circular-linear correlation),
#'   `effect_size` (for categorical terms the largest pairwise circular distance
#'   between group mean phases, in the same units as `period`, i.e. directly
#'   readable as hours; for numeric terms the circular-linear correlation itself),
#'   and `p_value`.
#'
#'   Two attributes carry the supporting detail: `"groupMeans"`, a `data.frame`
#'   (`term`, `group`, `n`, `circular_mean`) of per-group mean phases for the
#'   categorical terms, and `"nDropped"`, the number of samples excluded for having
#'   an `NA` or non-`"ok"` phase.
#' @examples
#' set.seed(1)
#' # 20 samples: cases collected around ZT6, controls around ZT14 - a design where
#' # collection time tracks the condition almost perfectly.
#' phase <- c(rnorm(10, mean = 6, sd = 1), rnorm(10, mean = 14, sd = 1))
#' mat <- matrix(rnorm(100 * 20), nrow = 100)
#' rownames(mat) <- paste0("g", seq_len(100))
#' colnames(mat) <- paste0("s", seq_len(20))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
#' SummarizedExperiment::colData(se)$phase <- wrapPhase(phase)
#' SummarizedExperiment::colData(se)$condition <- rep(c("case", "control"), each = 10)
#' phaseConfounding(se, ~ condition, nPerm = 199, seed = 1)
#' @export
phaseConfounding <- function(se, design, phase_col = "phase",
                              flag_col = "phase_flag", method = NULL,
                              tissue = NULL, species = "Mmusculus",
                              nPerm = 1000, period = 24, seed = NULL) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment")
  }
  if (!inherits(design, "formula")) {
    stop("design must be a formula, e.g. ~ condition")
  }

  cd <- SummarizedExperiment::colData(se)

  # Phase either already exists (the user ran estimatePhase(), or has real recorded
  # collection times) or is computed here. Computing it needs the same explicit
  # tissue/method estimatePhase() itself demands - never guessed (AGENTS.md 7.5).
  if (!phase_col %in% colnames(cd)) {
    if (is.null(method) || is.null(tissue)) {
      stop("se has no '", phase_col, "' column, so phase must be estimated first, ",
           "which requires both 'method' and 'tissue'. Either call estimatePhase() ",
           "yourself, or pass method= and tissue= here. ",
           "See availableReferenceModels().")
    }
    se <- estimatePhase(se, method = method, tissue = tissue, species = species)
    cd <- SummarizedExperiment::colData(se)
  }

  phase <- as.numeric(cd[[phase_col]])
  usable <- !is.na(phase)
  if (!is.null(flag_col) && flag_col %in% colnames(cd)) {
    usable <- usable & as.character(cd[[flag_col]]) == "ok"
  }
  nDropped <- sum(!usable)
  if (sum(usable) < 3) {
    stop("Fewer than 3 usable phase values (", sum(usable), " of ", length(phase),
         "); nothing to test. Check '", phase_col, "' and '", flag_col, "'.")
  }

  termLabels <- attr(stats::terms(design), "term.labels")
  if (length(termLabels) == 0) {
    stop("design has no terms to test, e.g. ~ condition")
  }
  isInteraction <- grepl(":", termLabels, fixed = TRUE)
  if (any(isInteraction)) {
    warning("Skipping interaction term(s): ",
            paste(termLabels[isInteraction], collapse = ", "),
            ". Only main effects are tested.")
    termLabels <- termLabels[!isInteraction]
  }
  missingTerms <- setdiff(termLabels, colnames(cd))
  if (length(missingTerms) > 0) {
    stop("design term(s) not found in colData(se): ",
         paste(missingTerms, collapse = ", "))
  }
  if (length(termLabels) == 0) {
    stop("No testable main-effect terms left in design")
  }

  if (!is.null(seed)) withr::local_seed(seed) else withr::local_preserve_seed()

  rows <- list()
  groupMeansList <- list()

  for (term in termLabels) {
    values <- cd[[term]][usable]
    theta <- phase[usable] / period * 2 * pi

    # A term can be unusable on its own (all one level, all NA) without the rest of
    # the design being unusable - report that row as NA rather than aborting the
    # whole diagnostic.
    keep <- !is.na(values)
    thetaT <- theta[keep]
    valuesT <- values[keep]
    isNumericTerm <- is.numeric(valuesT)

    statFun <- if (isNumericTerm) {
      function(th) .circularLinearStat(th, as.numeric(valuesT))
    } else {
      grp <- as.character(valuesT)
      function(th) .watsonWilliamsStat(th, grp)
    }

    nUsed <- length(thetaT)
    nGroups <- if (isNumericTerm) NA_integer_ else length(unique(as.character(valuesT)))

    observed <- if (nUsed >= 3 && (isNumericTerm || nGroups >= 2)) statFun(thetaT) else NA_real_

    pValue <- NA_real_
    if (!is.na(observed)) {
      permStats <- vapply(seq_len(nPerm), function(i) {
        s <- statFun(sample(thetaT))
        if (is.na(s)) -Inf else s
      }, numeric(1))
      pValue <- (1 + sum(permStats >= observed)) / (nPerm + 1)
    }

    effectSize <- NA_real_
    if (isNumericTerm) {
      effectSize <- observed
    } else if (!is.na(observed)) {
      groups <- sort(unique(as.character(valuesT)))
      means <- vapply(groups, function(g) {
        circularMean(phase[usable][keep][as.character(valuesT) == g], period = period)
      }, numeric(1))
      groupMeansList[[length(groupMeansList) + 1]] <- data.frame(
        term = term,
        group = groups,
        n = vapply(groups, function(g) sum(as.character(valuesT) == g), integer(1)),
        circular_mean = unname(means),
        row.names = NULL
      )
      # Largest pairwise gap between group mean phases: the plain-language version of
      # "how far apart were these groups collected", in hours.
      pairs <- utils::combn(length(groups), 2)
      effectSize <- max(abs(circularDiff(
        means[pairs[1, ]], means[pairs[2, ]], period = period
      )))
    }

    rows[[length(rows) + 1]] <- data.frame(
      term = term,
      type = if (isNumericTerm) "numeric" else "categorical",
      n = nUsed,
      n_groups = nGroups,
      statistic = observed,
      effect_size = effectSize,
      p_value = pValue,
      row.names = NULL
    )
  }

  out <- do.call(rbind, rows)
  attr(out, "groupMeans") <- if (length(groupMeansList) > 0) {
    do.call(rbind, groupMeansList)
  } else {
    NULL
  }
  attr(out, "nDropped") <- nDropped
  out
}
