# R3 (revision_and_fix_v2.md): SINGLE unified benchmark script.
#
# This script is the **sole source** for Tables 2, 3, 6, Figure 1 Panels A & B,
# and every downstream derived quantity (winner counts, oracle benefit, bootstrap
# CIs, significance tests). No other benchmark driver script should be needed for
# these outputs.
#
# Modes:
#   Default    : reads existing per-sample predictions from T6's committed
#                revision_persample_residuals_GSE54650.csv (fast, < 10 seconds)
#   --rerun    : re-runs the full 6-fold CV from scratch (~4 hours, dominated by
#                tauFisher's MetaCycle). Use for reproducibility verification.
#
# Seeds:
#   CV seed      = 1  (donor-split fold assignment; recorded in output)
#   Bootstrap seed = 42 (all CIs in this script; recorded in output)
#
# Outputs (all under benchmarks/results/):
#   revision2_unified_per_sample.csv  — the ONE tidy per-sample CSV
#   revision2_summary_per_tissue.csv  — per-tissue MAE for all methods (Tables 2 & 6)
#   revision2_significance_tests.csv  — bootstrap CIs, Friedman, Wilcoxon (Tables 2 & 6)
#   revision2_winner_counts.csv       — per-tissue best engine (Table 2 "Tissues won")
#   revision2_oracle_benefit.csv      — oracle vs Timetable-everywhere (Sections 3.2, 4.1)
#
# Pooled-mean definition (stated here and in every caption):
#   Unweighted mean of 12 per-tissue circular MAE values. Because n = 20 per tissue
#   is constant under complete-case scoring, this equals the mean of 240 per-sample
#   circular errors.
#
# Run from the project root:
#   Rscript benchmarks/run_unified_benchmark_v2.R            # fast, reads existing
#   Rscript benchmarks/run_unified_benchmark_v2.R --rerun    # full CV re-run
# ==============================================================================

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

# ==============================================================================
# 0. Configuration
# ==============================================================================
CV_SEED    <- 1L
BOOT_SEED  <- 42L
N_BOOT     <- 10000L
RERUN      <- "--rerun" %in% commandArgs(trailingOnly = TRUE)
OUTDIR     <- "benchmarks/results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

allEngines       <- c("timetable", "zeitzeiger", "timesignatr", "taufisher")
consensusEngines <- c("timetable", "taufisher")  # engines with bundled pretrained models

# Engine-specific arguments — identical to run_complete_case_benchmark.R / run_benchmark.R
# so that --rerun produces bit-identical results.
engineArgs <- list(
  timetable   = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher   = list(fit = list(method = c("JTK", "LS"), thres = 1.01,
                                nrep = 1, numbasis = 5)),
  zeitzeiger  = list(),
  timesignatr = list(fit = list(seed = 1))
)

# ==============================================================================
# 1. Per-sample predictions: run CV or read existing
# ==============================================================================
if (RERUN) {
  cat("=== MODE: full CV re-run (seed =", CV_SEED, ", nFolds = 6) ===\n")
  dataPath <- "benchmarks/data/GSE54650_se.rds"
  if (!file.exists(dataPath)) {
    stop("Missing ", dataPath, " - run inst/scripts/ingest_GSE54650.R first.")
  }
  se <- readRDS(dataPath)

  residualResults <- list()
  for (tis in unique(colData(se)$tissue)) {
    seT <- se[, colData(se)$tissue == tis]
    for (m in names(engineArgs)) {
      cat(sprintf("[%s] tissue=%s method=%s\n", format(Sys.time()), tis, m))
      bm <- tryCatch(
        benchmarkPhase(
          seT,
          truth_col   = "phase",
          methods     = m,
          assay_name  = "logExpr",
          nFolds      = 6,
          seed        = CV_SEED,
          engineArgs  = engineArgs[m]
        ),
        error = function(e) {
          message(sprintf("  FAILED (%s/%s): %s", tis, m, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(bm)) {
        residualResults[[paste(tis, m)]] <- cbind(tissue = tis, bm$residuals)
      }
    }
  }
  residDf <- do.call(rbind, residualResults)
  rownames(residDf) <- NULL
} else {
  cat("=== MODE: reading existing per-sample predictions ===\n")
  existingFile <- file.path(OUTDIR, "revision_persample_residuals_GSE54650.csv")
  if (!file.exists(existingFile)) {
    stop("Missing ", existingFile,
         ". Run with --rerun, or run benchmarks/run_complete_case_benchmark.R first.")
  }
  residDf <- utils::read.csv(existingFile, stringsAsFactors = FALSE)
  cat("Read", nrow(residDf), "rows from", existingFile, "\n")
}

# ==============================================================================
# 2. Complete-case filtering
# ==============================================================================
# Remove failed folds (tauFisher JTK_CYCLE failures)
residOk <- residDf[!residDf$fold_failed, ]

# Complete-case sample set: samples that every engine produced a valid prediction
# for. This is the common denominator across all four engines.
okByEngine <- lapply(allEngines, function(m) {
  s <- residOk[residOk$method == m & !is.na(residOk$pred), ]
  split(s$sample, s$tissue)
})
names(okByEngine) <- allEngines

tissues <- sort(unique(residOk$tissue))
completeCaseSamples <- lapply(tissues, function(tis) {
  Reduce(intersect, lapply(okByEngine, function(x) x[[tis]]))
})
names(completeCaseSamples) <- tissues

nCC <- sum(vapply(completeCaseSamples, length, integer(1)))
cat("Complete-case sample set:", nCC, "samples across", length(tissues), "tissues\n")
stopifnot(nCC == 240)  # 20 per tissue × 12 tissues

# ==============================================================================
# 3. Build single-engine rows (complete-case only)
# ==============================================================================
singleEngineRows <- list()
for (m in allEngines) {
  for (tis in tissues) {
    sub <- residOk[residOk$method == m & residOk$tissue == tis &
                     residOk$sample %in% completeCaseSamples[[tis]], ]
    if (nrow(sub) > 0) {
      singleEngineRows[[paste(tis, m)]] <- data.frame(
        dataset         = "GSE54650",
        tissue          = tis,
        fold            = sub$fold,
        sample_id       = sub$sample,
        method          = m,
        true_time       = sub$truth,
        predicted_time  = sub$pred,
        circular_error_h = sub$circular_error,
        stringsAsFactors = FALSE
      )
    }
  }
}

# ==============================================================================
# 4. Consensus computation (on complete-case samples)
# ==============================================================================
# Consensus combines Molecular Timetable + tauFisher predictions using
# .combineEnginePhases() — the only two engines with bundled pretrained models.
consensusRows <- list()
perTissueEngineMAE <- list()  # needed for accuracy weighting

for (tis in tissues) {
  sub <- residOk[residOk$tissue == tis & residOk$method %in% consensusEngines &
                   residOk$sample %in% completeCaseSamples[[tis]], ]

  # Per-engine MAE for accuracy weighting
  perEngineMAE <- vapply(consensusEngines, function(m) {
    mean(sub$circular_error[sub$method == m], na.rm = TRUE)
  }, numeric(1))
  perTissueEngineMAE[[tis]] <- perEngineMAE

  # Reshape to wide: one row per sample, columns for each engine's prediction
  wide <- stats::reshape(
    sub[, c("sample", "method", "pred", "truth", "fold")],
    idvar = c("sample", "truth", "fold"), timevar = "method", direction = "wide"
  )
  predCols <- paste0("pred.", consensusEngines)
  if (!all(predCols %in% colnames(wide))) {
    warning("Missing consensus engine predictions for tissue ", tis)
    next
  }
  phaseMat <- as.matrix(wide[, predCols, drop = FALSE])

  # Two consensus schemes
  schemes <- list(
    consensus_equal    = rep(1, length(consensusEngines)) / length(consensusEngines),
    consensus_accuracy = (1 / perEngineMAE) / sum(1 / perEngineMAE)
  )

  for (schemeName in names(schemes)) {
    combined <- .combineEnginePhases(phaseMat, unname(schemes[[schemeName]]))
    ok <- !is.na(combined$phase)
    consensusRows[[paste(tis, schemeName)]] <- data.frame(
      dataset         = "GSE54650",
      tissue          = tis,
      fold            = wide$fold[ok],
      sample_id       = wide$sample[ok],
      method          = schemeName,
      true_time       = wide$truth[ok],
      predicted_time  = combined$phase[ok],
      circular_error_h = circularError(combined$phase[ok], wide$truth[ok]),
      stringsAsFactors = FALSE
    )
  }
}

# ==============================================================================
# 5. Output: revision2_unified_per_sample.csv
# ==============================================================================
allRows <- c(singleEngineRows, consensusRows)
perSample <- do.call(rbind, allRows)
perSample$seed <- CV_SEED
rownames(perSample) <- NULL

# Sort for reproducibility
perSample <- perSample[order(perSample$dataset, perSample$tissue,
                              perSample$method, perSample$fold,
                              perSample$sample_id), ]

outFile <- file.path(OUTDIR, "revision2_unified_per_sample.csv")
utils::write.csv(perSample, outFile, row.names = FALSE)
cat("\nWrote", outFile, "(", nrow(perSample), "rows )\n")

# ==============================================================================
# 6. Aggregations — all derived from perSample, nothing else
# ==============================================================================

# --- 6a. Per-tissue MAE for every method ----------------------------------
perTissueMAE <- stats::aggregate(
  circular_error_h ~ tissue + method,
  data = perSample,
  FUN  = mean
)
colnames(perTissueMAE)[3] <- "circular_mae"
perTissueMAE <- perTissueMAE[order(perTissueMAE$method, perTissueMAE$tissue), ]

# Also record n per cell (should be 20 for all complete-case)
perTissueN <- stats::aggregate(
  circular_error_h ~ tissue + method,
  data = perSample,
  FUN  = length
)
colnames(perTissueN)[3] <- "n"
perTissueSummary <- merge(perTissueMAE, perTissueN, by = c("tissue", "method"))

outSummary <- file.path(OUTDIR, "revision2_summary_per_tissue.csv")
utils::write.csv(perTissueSummary, outSummary, row.names = FALSE)
cat("Wrote", outSummary, "\n")

# --- 6b. Bootstrap CIs (tissue-level resampling, single seed) -------------
# Resampling unit: 12 tissues. This is stated in every caption.
# We resample the 12 per-tissue MAE values, not per-sample errors.
bootCI <- function(x, nBoot = N_BOOT, seed = BOOT_SEED) {
  set.seed(seed)
  n <- length(x)
  boots <- vapply(seq_len(nBoot), function(i) {
    mean(sample(x, n, replace = TRUE))
  }, numeric(1))
  stats::quantile(boots, c(0.025, 0.975), names = FALSE)
}

allMethods <- sort(unique(perSample$method))

# --- 6c. Significance tests -----------------------------------------------
sigRows <- list()

# -- Table 2 block: four single engines ------------------------------------
mat_engines <- stats::reshape(
  perTissueMAE[perTissueMAE$method %in% allEngines,
               c("tissue", "method", "circular_mae")],
  idvar = "tissue", timevar = "method", direction = "wide"
)
colnames(mat_engines) <- sub("^circular_mae\\.", "", colnames(mat_engines))
mat_engines <- mat_engines[, c("tissue", allEngines)]
stopifnot(nrow(mat_engines) == 12, !anyNA(mat_engines))

# Friedman omnibus
fried_engines <- stats::friedman.test(as.matrix(mat_engines[, allEngines]))
sigRows[[length(sigRows) + 1]] <- data.frame(
  table_source = "Table2", comparison = "all_single_engines", test = "friedman",
  statistic = unname(fried_engines$statistic), p_raw = fried_engines$p.value,
  p_holm = NA, n_tissues = 12, effect_size_h = NA, ci_lo_h = NA, ci_hi_h = NA,
  bootstrap_seed = NA, bootstrap_n = NA, resampling_unit = NA
)

# Pairwise Wilcoxon (all 6 contrasts, Holm-corrected over 6)
pairs_engines <- utils::combn(allEngines, 2, simplify = FALSE)
pvals_engines <- vapply(pairs_engines, function(p) {
  stats::wilcox.test(mat_engines[[p[1]]], mat_engines[[p[2]]],
                     paired = TRUE, exact = FALSE)$p.value
}, numeric(1))
pholm_engines <- stats::p.adjust(pvals_engines, method = "holm")

for (i in seq_along(pairs_engines)) {
  p <- pairs_engines[[i]]
  eff <- mean(mat_engines[[p[1]]] - mat_engines[[p[2]]])
  sigRows[[length(sigRows) + 1]] <- data.frame(
    table_source = "Table2",
    comparison = paste0(p[1], "_vs_", p[2]),
    test = "wilcoxon_paired", statistic = NA,
    p_raw = pvals_engines[i], p_holm = pholm_engines[i],
    n_tissues = 12, effect_size_h = eff, ci_lo_h = NA, ci_hi_h = NA,
    bootstrap_seed = NA, bootstrap_n = NA, resampling_unit = NA
  )
}

# Bootstrap CIs for single engines
for (eng in allEngines) {
  ci <- bootCI(mat_engines[[eng]])
  sigRows[[length(sigRows) + 1]] <- data.frame(
    table_source = "Table2",
    comparison = paste0("pooled_mean_", eng),
    test = "bootstrap_tissue", statistic = NA, p_raw = NA, p_holm = NA,
    n_tissues = 12, effect_size_h = mean(mat_engines[[eng]]),
    ci_lo_h = ci[1], ci_hi_h = ci[2],
    bootstrap_seed = BOOT_SEED, bootstrap_n = N_BOOT,
    resampling_unit = "tissue"
  )
}

# -- Table 6 block: consensus + timetable + taufisher ----------------------
methods6 <- c("timetable", "taufisher", "consensus_equal", "consensus_accuracy")
mat_consensus <- stats::reshape(
  perTissueMAE[perTissueMAE$method %in% methods6,
               c("tissue", "method", "circular_mae")],
  idvar = "tissue", timevar = "method", direction = "wide"
)
colnames(mat_consensus) <- sub("^circular_mae\\.", "", colnames(mat_consensus))
mat_consensus <- mat_consensus[, c("tissue", methods6)]
stopifnot(nrow(mat_consensus) == 12, !anyNA(mat_consensus))

# Friedman omnibus
fried_cons <- stats::friedman.test(as.matrix(mat_consensus[, methods6]))
sigRows[[length(sigRows) + 1]] <- data.frame(
  table_source = "Table6", comparison = "all_methods", test = "friedman",
  statistic = unname(fried_cons$statistic), p_raw = fried_cons$p.value,
  p_holm = NA, n_tissues = 12, effect_size_h = NA, ci_lo_h = NA, ci_hi_h = NA,
  bootstrap_seed = NA, bootstrap_n = NA, resampling_unit = NA
)

# Pairwise contrasts for Table 6: every contrast against the best single engine
# (timetable) plus the two consensus schemes against each other = 4 contrasts.
# Pre-specification: "every contrast against the point-estimate leader, plus
# the two weighting schemes against each other."
contrasts6 <- list(
  c("timetable", "consensus_accuracy"),
  c("timetable", "consensus_equal"),
  c("timetable", "taufisher"),
  c("consensus_accuracy", "consensus_equal")
)
pvals6 <- vapply(contrasts6, function(p) {
  stats::wilcox.test(mat_consensus[[p[1]]], mat_consensus[[p[2]]],
                     paired = TRUE, exact = FALSE)$p.value
}, numeric(1))
pholm6 <- stats::p.adjust(pvals6, method = "holm")

for (i in seq_along(contrasts6)) {
  p <- contrasts6[[i]]
  eff <- mean(mat_consensus[[p[1]]] - mat_consensus[[p[2]]])
  sigRows[[length(sigRows) + 1]] <- data.frame(
    table_source = "Table6",
    comparison = paste0(p[1], "_vs_", p[2]),
    test = "wilcoxon_paired", statistic = NA,
    p_raw = pvals6[i], p_holm = pholm6[i],
    n_tissues = 12, effect_size_h = eff, ci_lo_h = NA, ci_hi_h = NA,
    bootstrap_seed = NA, bootstrap_n = NA, resampling_unit = NA
  )
}

# Bootstrap CIs for Table 6 methods
for (m in methods6) {
  ci <- bootCI(mat_consensus[[m]])
  sigRows[[length(sigRows) + 1]] <- data.frame(
    table_source = "Table6",
    comparison = paste0("pooled_mean_", m),
    test = "bootstrap_tissue", statistic = NA, p_raw = NA, p_holm = NA,
    n_tissues = 12, effect_size_h = mean(mat_consensus[[m]]),
    ci_lo_h = ci[1], ci_hi_h = ci[2],
    bootstrap_seed = BOOT_SEED, bootstrap_n = N_BOOT,
    resampling_unit = "tissue"
  )
}

sigDf <- do.call(rbind, sigRows)
rownames(sigDf) <- NULL
outSig <- file.path(OUTDIR, "revision2_significance_tests.csv")
utils::write.csv(sigDf, outSig, row.names = FALSE)
cat("Wrote", outSig, "\n")

# --- 6d. Winner counts (Table 2 "Tissues won") ----------------------------
# Per-tissue winner: engine with lowest complete-case MAE.
# Tie rule: if two engines are equal at full (double) precision, report as tie.
# At display precision (2 decimal) ties that are not ties at full precision are
# resolved by reporting the true winner and noting the displayed rounding.
winnerRows <- list()
for (tis in tissues) {
  tisMAE <- mat_engines[mat_engines$tissue == tis, allEngines, drop = FALSE]
  vals <- unlist(tisMAE)
  bestVal <- min(vals)
  winners <- allEngines[vals == bestVal]
  winnerRows[[tis]] <- data.frame(
    tissue = tis,
    best_mae = bestVal,
    best_mae_display = round(bestVal, 2),
    winner = paste(winners, collapse = " / "),
    is_tie  = length(winners) > 1,
    stringsAsFactors = FALSE
  )
}
winnerDf <- do.call(rbind, winnerRows)
rownames(winnerDf) <- NULL

# Summary counts
winnerCounts <- table(unlist(strsplit(winnerDf$winner, " / ")))
tieCounts <- sum(winnerDf$is_tie)
cat("\n== Tissues won (complete-case, full precision) ==\n")
for (eng in allEngines) {
  nWon <- sum(winnerDf$winner == eng)  # outright wins
  nTied <- sum(grepl(eng, winnerDf$winner) & winnerDf$is_tie)
  if (nTied > 0) {
    cat(sprintf("  %s: %d (+%d tie)\n", eng, nWon, nTied))
  } else {
    cat(sprintf("  %s: %d\n", eng, nWon))
  }
}

# Add summary row to winner CSV
summaryRow <- data.frame(
  tissue = "SUMMARY",
  best_mae = NA,
  best_mae_display = NA,
  winner = paste(
    vapply(allEngines, function(eng) {
      nWon <- sum(winnerDf$winner == eng)
      nTied <- sum(grepl(eng, winnerDf$winner) & winnerDf$is_tie)
      if (nTied > 0) paste0(eng, "=", nWon, "(+", nTied, " tie)")
      else paste0(eng, "=", nWon)
    }, character(1)),
    collapse = "; "
  ),
  is_tie = NA,
  stringsAsFactors = FALSE
)
winnerDf <- rbind(winnerDf, summaryRow)

outWinners <- file.path(OUTDIR, "revision2_winner_counts.csv")
utils::write.csv(winnerDf, outWinners, row.names = FALSE)
cat("Wrote", outWinners, "\n")

# --- 6e. Oracle benefit (Sections 3.2, 4.1) --------------------------------
# Oracle: for each tissue, pick the engine with the lowest MAE.
# Compare against Timetable-everywhere.
oraclePerTissue <- vapply(tissues, function(tis) {
  vals <- unlist(mat_engines[mat_engines$tissue == tis, allEngines])
  min(vals)
}, numeric(1))
timetablePerTissue <- mat_engines$timetable

oraclePooled    <- mean(oraclePerTissue)
timetablePooled <- mean(timetablePerTissue)
oracleBenefit   <- timetablePooled - oraclePooled

cat(sprintf("\n== Oracle benefit ==\n"))
cat(sprintf("  Oracle pooled MAE:    %.4f h\n", oraclePooled))
cat(sprintf("  Timetable pooled MAE: %.4f h\n", timetablePooled))
cat(sprintf("  Benefit:              %.4f h (%.1f min)\n",
            oracleBenefit, oracleBenefit * 60))

# Bootstrap CI on the oracle benefit (resample 12 tissues, 10,000 draws)
set.seed(BOOT_SEED)
bootBenefit <- vapply(seq_len(N_BOOT), function(i) {
  idx <- sample(12, 12, replace = TRUE)
  mean(timetablePerTissue[idx]) - mean(oraclePerTissue[idx])
}, numeric(1))
benefitCI <- stats::quantile(bootBenefit, c(0.025, 0.975), names = FALSE)
cat(sprintf("  95%% CI:               [%.4f, %.4f]\n", benefitCI[1], benefitCI[2]))
cat(sprintf("  CI spans zero:         %s\n", benefitCI[1] <= 0))

# Per-engine marginal contribution:
# For each of ZeitZeiger and TimeSignature, the benefit of adding that engine
# to the {Timetable, tauFisher} pool.
basePool <- c("timetable", "taufisher")
baseOracle <- vapply(tissues, function(tis) {
  vals <- unlist(mat_engines[mat_engines$tissue == tis, basePool])
  min(vals)
}, numeric(1))
basePooledMAE <- mean(baseOracle)

marginalBenefit <- list()
for (extraEng in c("zeitzeiger", "timesignatr")) {
  extPool <- c(basePool, extraEng)
  extOracle <- vapply(tissues, function(tis) {
    vals <- unlist(mat_engines[mat_engines$tissue == tis, extPool])
    min(vals)
  }, numeric(1))
  extPooledMAE <- mean(extOracle)
  mb <- basePooledMAE - extPooledMAE
  marginalBenefit[[extraEng]] <- mb
  cat(sprintf("  Adding %s to {Timetable,tauFisher}: buys %.4f h (%.1f min)\n",
              extraEng, mb, mb * 60))
}

oracleOut <- data.frame(
  quantity = c("oracle_pooled_mae", "timetable_pooled_mae", "oracle_benefit",
               "oracle_benefit_ci_lo", "oracle_benefit_ci_hi",
               "oracle_benefit_ci_spans_zero",
               "base_pool_mae_timetable_taufisher",
               "marginal_benefit_zeitzeiger", "marginal_benefit_timesignatr"),
  value = c(oraclePooled, timetablePooled, oracleBenefit,
            benefitCI[1], benefitCI[2],
            as.numeric(benefitCI[1] <= 0),
            basePooledMAE,
            marginalBenefit[["zeitzeiger"]],
            marginalBenefit[["timesignatr"]]),
  bootstrap_seed = c(NA, NA, NA, BOOT_SEED, BOOT_SEED, NA, NA, NA, NA),
  bootstrap_n    = c(NA, NA, NA, N_BOOT, N_BOOT, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)
outOracle <- file.path(OUTDIR, "revision2_oracle_benefit.csv")
utils::write.csv(oracleOut, outOracle, row.names = FALSE)
cat("Wrote", outOracle, "\n")

# ==============================================================================
# 7. Console summary for verification
# ==============================================================================
cat("\n========================================\n")
cat("UNIFIED BENCHMARK SUMMARY (R3)\n")
cat("========================================\n")

cat("\n-- Table 2: Complete-case MAE per engine --\n")
pooledEngines <- stats::aggregate(circular_mae ~ method, data = perTissueMAE[perTissueMAE$method %in% allEngines, ], FUN = mean)
pooledEngines <- pooledEngines[order(pooledEngines$circular_mae), ]
for (i in seq_len(nrow(pooledEngines))) {
  eng <- pooledEngines$method[i]
  mae <- pooledEngines$circular_mae[i]
  ci_row <- sigDf[sigDf$comparison == paste0("pooled_mean_", eng) & sigDf$table_source == "Table2", ]
  cat(sprintf("  %-12s  %.3f  [%.3f, %.3f]\n", eng, mae, ci_row$ci_lo_h, ci_row$ci_hi_h))
}

cat("\n-- Table 6: Consensus vs single engines --\n")
pooledAll <- stats::aggregate(circular_mae ~ method, data = perTissueMAE[perTissueMAE$method %in% methods6, ], FUN = mean)
pooledAll <- pooledAll[order(pooledAll$circular_mae), ]
for (i in seq_len(nrow(pooledAll))) {
  m <- pooledAll$method[i]
  mae <- pooledAll$circular_mae[i]
  ci_row <- sigDf[sigDf$comparison == paste0("pooled_mean_", m) & sigDf$table_source == "Table6", ]
  cat(sprintf("  %-22s  %.3f  [%.3f, %.3f]\n", m, mae, ci_row$ci_lo_h, ci_row$ci_hi_h))
}

cat("\n-- Consistency check --\n")
# Table 2 and Table 6 should report identical Timetable and tauFisher values
for (eng in consensusEngines) {
  t2 <- sigDf[sigDf$comparison == paste0("pooled_mean_", eng) & sigDf$table_source == "Table2", ]
  t6 <- sigDf[sigDf$comparison == paste0("pooled_mean_", eng) & sigDf$table_source == "Table6", ]
  mae_match <- t2$effect_size_h == t6$effect_size_h
  ci_match  <- t2$ci_lo_h == t6$ci_lo_h && t2$ci_hi_h == t6$ci_hi_h
  cat(sprintf("  %s: MAE match=%s, CI match=%s\n", eng, mae_match, ci_match))
}

cat("\n-- Output files --\n")
cat("  ", outFile, "\n")
cat("  ", outSummary, "\n")
cat("  ", outSig, "\n")
cat("  ", outWinners, "\n")
cat("  ", outOracle, "\n")
cat("\nDone. CV seed =", CV_SEED, ", bootstrap seed =", BOOT_SEED, "\n")
