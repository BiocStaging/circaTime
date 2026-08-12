# T5 (revision_and_fix.md): document each engine's prediction resolution and
# show sub-2h MAE on GSE54650 (12 atlas timepoints, 2h apart) is not a grid
# artefact. Reuses T6's already-committed per-sample predictions
# (benchmarks/results/revision_persample_residuals_GSE54650.csv) rather than
# rerunning the benchmark - those predictions are real script output already
# checked in, per ground rule 1 (never invent numbers, everything traces to a
# committed CSV).
# Run from the project root: Rscript benchmarks/run_grid_resolution_check.R

suppressMessages(library(pkgload))
pkgload::load_all(".", quiet = TRUE)

set.seed(1)

df <- utils::read.csv("benchmarks/results/revision_persample_residuals_GSE54650.csv")
df <- df[!df$fold_failed, ]

# Each engine's own prediction step size (see R/engines-*.R): Molecular
# Timetable searches a 288-point grid over 24h (gridN=288, R/engines-timetable.R);
# tauFisher classifies onto an integer hour via nnet::multinom
# (R/engines-taufisher.R); ZeitZeiger's timeRange grid is only the initial value
# for a continuous likelihood optimisation, not the final output
# (R/engines-zeitzeiger.R); TimeSignatR's predTimeStamp() is a continuous glmnet
# regression output (R/engines-timesignatr.R). Confirmed empirically below via
# unique-value counts, not just read off the code.
stepHours <- c(timetable = 24 / 288, taufisher = 1, zeitzeiger = NA, timesignatr = NA)

atlasHours <- seq(0, 22, by = 2) # GSE54650's own 12 collection timepoints

onAtlasGrid <- function(pred, tol = 1e-6) {
  vapply(pred, function(p) {
    d <- abs(((p - atlasHours + 12) %% 24) - 12)
    any(d < tol)
  }, logical(1))
}

rows <- list()
for (m in unique(df$method)) {
  sub <- df[df$method == m, ]
  step <- stepHours[[m]]

  nUnique <- length(unique(sub$pred))
  fracOnAtlas <- mean(onAtlasGrid(sub$pred))
  maeObserved <- circularMAE(sub$pred, sub$truth)

  if (!is.na(step)) {
    jitter <- stats::runif(nrow(sub), -step / 2, step / 2)
    predJit <- (sub$pred + jitter) %% 24
    maeJittered <- circularMAE(predJit, sub$truth)
  } else {
    maeJittered <- NA_real_
  }

  rows[[length(rows) + 1]] <- data.frame(
    engine = m, n_predictions = nrow(sub), n_unique_pred_values = nUnique,
    prediction_step_hours = step, frac_on_atlas_grid = fracOnAtlas,
    mae_observed = maeObserved, mae_after_grid_jitter = maeJittered
  )
}

out <- do.call(rbind, rows)
rownames(out) <- NULL
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "benchmarks/results/revision_grid_resolution.csv", row.names = FALSE)

cat("\nChance baseline: expected circular error of a uniform random guess on a 24h\n")
cat("circle. E[min(|X|, 24-|X|)] for X ~ Uniform(-12, 12) = 24/4 = 6h.\n")

print(out)
cat("\nWrote benchmarks/results/revision_grid_resolution.csv\n")
