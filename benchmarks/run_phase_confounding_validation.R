# T2 (revision_and_fix.md): validate phaseConfounding()'s type-I error and power
# with a real simulation study, calling the exported function itself (not a
# reimplementation of its statistics) so the validation covers the actual code
# path a user would run.
#
# Grid, per the task spec:
#  - Type-I error: n/group in (10,20,50,100), groups in (2,3), phase distribution
#    in (uniform, von Mises kappa=2, von Mises kappa=5); >=2000 reps/cell.
#  - Power: same grid, plus a mean-phase offset (0,1,2,3,4,6 h) applied to one
#    group; >=2000 reps/cell.
#  - Robustness: null simulation with unequal concentration between two groups
#    (kappa=5 vs kappa=1) across the same n/group grid; >=2000 reps/cell.
#  - Numeric-term branch (Mardia/Johnson-Wehrly): null (phase independent of a
#    continuous covariate) and graded association, across the same n grid;
#    >=2000 reps/cell.
#
# Runtime: dominated by phaseConfounding()'s own permutation loop (nPerm=199
# here, enough resolution near alpha=0.05); ~50-70 min total for the full grid.
# Run from the project root: Rscript benchmarks/run_phase_confounding_validation.R

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

set.seed(1)
nReps <- 2000
nPerm <- 199
alpha <- 0.05

# Best & Fisher (1979) rejection algorithm for von Mises deviates. kappa=0 means
# "no concentration parameter": sample uniform on the circle directly.
rvonmises <- function(n, mu, kappa) {
  if (kappa <= 1e-8) return((stats::runif(n, -pi, pi) + mu) %% (2 * pi))
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a)) / (2 * kappa)
  r <- (1 + b^2) / (2 * b)
  out <- numeric(n)
  for (i in seq_len(n)) {
    repeat {
      u1 <- stats::runif(1); z <- cos(pi * u1); f <- (1 + r * z) / (r + z); c <- kappa * (r - f)
      u2 <- stats::runif(1)
      if (c * (2 - c) - u2 > 0 || log(c / u2) + 1 - c >= 0) break
    }
    u3 <- stats::runif(1)
    theta <- if (u3 > 0.5) acos(f) else -acos(f)
    out[i] <- theta + mu
  }
  out %% (2 * pi)
}

hoursToRad <- function(h) h / 24 * 2 * pi
radToHours <- function(r) r / (2 * pi) * 24

# Cache one minimal SE per total sample size so the assay isn't reallocated
# every replicate; only colData is replaced per call.
seCache <- new.env()
getSE <- function(nTotal) {
  key <- as.character(nTotal)
  if (is.null(seCache[[key]])) {
    mat <- matrix(0, nrow = 1, ncol = nTotal)
    rownames(mat) <- "g1"; colnames(mat) <- paste0("s", seq_len(nTotal))
    seCache[[key]] <- SummarizedExperiment(assays = list(logExpr = mat))
  }
  seCache[[key]]
}

runCategorical <- function(nPerGroup, nGroups, kappaVec, muOffsetHours = NULL) {
  nTotal <- nPerGroup * nGroups
  se <- getSE(nTotal)
  grp <- rep(paste0("g", seq_len(nGroups)), each = nPerGroup)
  muHours <- rep(12, nGroups)
  if (!is.null(muOffsetHours)) muHours[nGroups] <- muHours[nGroups] + muOffsetHours
  kappaVec <- rep(kappaVec, length.out = nGroups)

  rejections <- 0L
  for (i in seq_len(nReps)) {
    theta <- unlist(lapply(seq_len(nGroups), function(g) {
      rvonmises(nPerGroup, hoursToRad(muHours[g]), kappaVec[g])
    }))
    colData(se)$phase <- wrapPhase(radToHours(theta))
    colData(se)$grp <- grp
    r <- phaseConfounding(se, ~grp, flag_col = NULL, nPerm = nPerm)
    if (!is.na(r$p_value[1]) && r$p_value[1] < alpha) rejections <- rejections + 1L
  }
  rejections / nReps
}

runNumeric <- function(n, kappa, assocHoursPerUnit) {
  se <- getSE(n)
  x <- stats::runif(n, 0, 1)
  muHours <- 12 + assocHoursPerUnit * x
  rejections <- 0L
  for (i in seq_len(nReps)) {
    theta <- vapply(muHours, function(m) rvonmises(1, hoursToRad(m), kappa), numeric(1))
    colData(se)$phase <- wrapPhase(radToHours(theta))
    colData(se)$x <- x
    r <- phaseConfounding(se, ~x, flag_col = NULL, nPerm = nPerm)
    if (!is.na(r$p_value[1]) && r$p_value[1] < alpha) rejections <- rejections + 1L
  }
  rejections / nReps
}

nPerGroupGrid <- c(10, 20, 50, 100)
groupsGrid <- c(2, 3)
distGrid <- list(uniform = 0, kappa2 = 2, kappa5 = 5)
offsetGrid <- c(0, 1, 2, 3, 4, 6)

rows <- list()

## ---- Type-I error --------------------------------------------------------
cat("== Type-I error grid ==\n")
for (nGrp in groupsGrid) for (nPer in nPerGroupGrid) for (distName in names(distGrid)) {
  kappa <- distGrid[[distName]]
  rate <- runCategorical(nPer, nGrp, kappa)
  cat(sprintf("[%s] groups=%d n/group=%d dist=%s rejection_rate=%.4f\n",
              format(Sys.time()), nGrp, nPer, distName, rate))
  rows[[length(rows) + 1]] <- data.frame(
    study = "type1_error", n_per_group = nPer, n_groups = nGrp,
    phase_dist = distName, offset_hours = NA, rejection_rate = rate, n_reps = nReps
  )
}

## ---- Power ----------------------------------------------------------------
cat("== Power grid ==\n")
for (nGrp in groupsGrid) for (nPer in nPerGroupGrid) for (distName in names(distGrid)) {
  kappa <- distGrid[[distName]]
  for (offset in offsetGrid) {
    rate <- runCategorical(nPer, nGrp, kappa, muOffsetHours = offset)
    cat(sprintf("[%s] groups=%d n/group=%d dist=%s offset=%dh rejection_rate=%.4f\n",
                format(Sys.time()), nGrp, nPer, distName, offset, rate))
    rows[[length(rows) + 1]] <- data.frame(
      study = "power", n_per_group = nPer, n_groups = nGrp,
      phase_dist = distName, offset_hours = offset, rejection_rate = rate, n_reps = nReps
    )
  }
}

## ---- Robustness: unequal concentration between two groups ----------------
cat("== Robustness grid (unequal concentration, kappa 5 vs 1) ==\n")
for (nPer in nPerGroupGrid) {
  rate <- runCategorical(nPer, 2, kappaVec = c(5, 1))
  cat(sprintf("[%s] n/group=%d rejection_rate=%.4f\n", format(Sys.time()), nPer, rate))
  rows[[length(rows) + 1]] <- data.frame(
    study = "robustness_unequal_kappa", n_per_group = nPer, n_groups = 2,
    phase_dist = "kappa5_vs_kappa1", offset_hours = NA, rejection_rate = rate, n_reps = nReps
  )
}

## ---- Numeric term (Mardia / Johnson-Wehrly) -------------------------------
cat("== Numeric-term grid ==\n")
assocGrid <- c(0, 2, 4, 8)
for (n in nPerGroupGrid) for (assoc in assocGrid) {
  rate <- runNumeric(n, kappa = 2, assocHoursPerUnit = assoc)
  cat(sprintf("[%s] n=%d assoc_hours_per_unit=%d rejection_rate=%.4f\n",
              format(Sys.time()), n, assoc, rate))
  rows[[length(rows) + 1]] <- data.frame(
    study = if (assoc == 0) "numeric_type1_error" else "numeric_power",
    n_per_group = n, n_groups = NA, phase_dist = paste0("kappa2_assoc", assoc),
    offset_hours = assoc, rejection_rate = rate, n_reps = nReps
  )
}

out <- do.call(rbind, rows)
rownames(out) <- NULL
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "benchmarks/results/revision_phase_confounding_validation.csv", row.names = FALSE)
cat("\nWrote benchmarks/results/revision_phase_confounding_validation.csv\n")
