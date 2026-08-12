# EXPERIMENTAL PREVIEW - NOT part of the shipped package, does NOT touch
# R/sysdata.rda, and its output is NOT loaded by estimatePhase(). Requested as a
# quick look at what a human (GSE56931) and a pseudobulk (GSE145197) reference
# model would look like, using the exact same two-pass procedure as
# inst/scripts/train_reference_models.R (the real, shipped training script - read
# that one first, this one only differs in which cohorts it points at and where
# it writes its output).
#
# Why these two are NOT bundled for real (notes/decisions.md Section 10/12,
# paper/manuscript.tex Section "Discussion and limitations"): GSE56931 is a
# single tissue (whole blood) with only 130/249 unconfounded baseline samples
# from 14 donors, and GSE145197 is a single 10-sample, 4-Zeitgeber-timepoint
# pseudobulk cohort. Both residual pools and gene-selection steps below will be
# thin by construction - this script exists to SHOW that concretely, not to
# produce something fit to ship. Treat every number this prints as a rough,
# unvalidated preview, not a benchmark result.
#
# Run from the project root: Rscript benchmarks/train_reference_models_experimental_preview.R
# Requires benchmarks/data/GSE56931_se.rds and benchmarks/data/GSE145197_se.rds
# (run inst/scripts/ingest_GSE56931.R / ingest_GSE145197.R first) and
# benchmarks/results/phase3_amplitude_threshold_GSE54650_timetable.csv (same
# reused-not-rederived amplitude threshold as the real training script, an even
# larger simplifying assumption here since it was derived from mouse liver
# microarray, not human blood or mouse hepatocyte pseudobulk data).

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

thrFile <- "benchmarks/results/phase3_amplitude_threshold_GSE54650_timetable.csv"
if (!file.exists(thrFile)) {
  stop("Run benchmarks/run_phase3_calibration.R first to produce ", thrFile)
}
ampThreshold <- utils::read.csv(thrFile)$threshold[1]

engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01, nrep = 1, numbasis = 5))
)
seed <- 1

# Same trainPlatform() logic as inst/scripts/train_reference_models.R, with two
# differences suited to a single-cohort preview rather than a 12-tissue atlas:
# `nFolds` is a parameter (donor counts here are much smaller than GSE54650's 12
# tissues x lots-of-donors), and it prints each fitted model's in-memory size so
# an oversized tauFisher fit (as happened for GSE54651/rnaseq, notes/decisions.md
# Section 12) is visible immediately rather than discovered at a package-size
# check that does not even apply here.
trainPreview <- function(se, species, platform, tissueLabel, citationNote,
                          methods = c("timetable", "taufisher"), nFolds = 5) {
  seT <- se
  colData(seT)$tissue <- tissueLabel # trainPlatform()'s real counterpart splits by
  # tissue; these cohorts are one group, so this is a constant label, not a split.

  cat(sprintf("== [%s/%s/%s] Pass 1: cross-validated residuals (nFolds=%d) ==\n",
              species, platform, tissueLabel, nFolds))
  bm <- benchmarkPhase(seT,
    truth_col = "phase", methods = methods,
    donor_col = "donor", nFolds = nFolds, seed = seed, engineArgs = engineArgs[methods]
  )
  res <- bm$residuals
  res <- res[!res$fold_failed, ]
  res$residual <- circularDiff(res$pred, res$truth)
  cat(sprintf(
    "Held-out circular MAE: %s\n",
    paste(sprintf("%s=%.2fh (n=%d)", bm$summary$method, bm$summary$circular_mae, bm$summary$n), collapse = ", ")
  ))

  cat(sprintf("\n== [%s/%s/%s] Pass 2: full fit + amplitude PCA basis ==\n", species, platform, tissueLabel))
  exprMat <- assay(seT, "logExpr")
  truth <- colData(seT)$phase

  amp <- clockAmplitude(exprMat)
  ampGenes <- attr(amp, "genesUsed")
  ampCenter <- attr(amp, "center")[["all"]]
  ampScale <- attr(amp, "scale")[["all"]]
  ampRotation <- attr(amp, "rotation")[["all"]]

  models <- list()
  if ("timetable" %in% methods) {
    fit <- do.call(fitMolecularTimetable, c(list(trainMat = exprMat, trainTime = truth), engineArgs$timetable$fit))
    models$timetable <- fit
  }
  if ("taufisher" %in% methods) {
    fit <- do.call(fitTauFisher, c(list(trainMat = exprMat, trainTime = truth, seed = seed), engineArgs$taufisher$fit))
    environment(fit$Model$terms) <- baseenv() # same stripping fix as the real script
    fit$Model$call <- NULL
    cat(sprintf("tauFisher selected %d periodic genes (watch for the GSE54651/rnaseq size blowup)\n",
                length(fit$Genes)))
    models$taufisher <- fit
  }

  referenceModels <- list()
  for (m in names(models)) {
    residualPool <- res$residual[res$method == m]
    key <- paste(species, tissueLabel, m, platform, sep = ".")
    rm <- methods::new("PhaseReferenceModel",
      species = species, tissue = tissueLabel, method = m, platform = platform,
      fit = models[[m]],
      # NOT leave-tissue-out here (there is only one group) - this residual pool
      # is the SAME cohort's own CV residuals, a materially weaker guarantee than
      # every other bundled model's leave-tissue-out pooling. Said explicitly in
      # trainedOn below, not left implicit.
      residualPool = residualPool,
      ampCenter = ampCenter, ampScale = ampScale, ampRotation = ampRotation,
      ampGenes = ampGenes, ampThreshold = ampThreshold,
      trainedOn = paste0(
        citationNote, ", EXPERIMENTAL PREVIEW - not bundled, not shipped. ",
        "residual pool is this SAME cohort's own within-dataset CV residuals ",
        "(no leave-tissue-out - there is only one group here), n=", length(residualPool),
        "; amplitude threshold inherited from Phase 3's GSE54650 mouse-liver-",
        "microarray derivation, not re-derived for this cohort - see notes/decisions.md"
      )
    )
    cat(sprintf("%s: %.1f KB, %d held-out residuals\n", key, as.numeric(object.size(rm)) / 1024, length(residualPool)))
    referenceModels[[key]] <- rm
  }
  referenceModels
}

humanFile <- "benchmarks/data/GSE56931_se.rds"
if (!file.exists(humanFile)) stop("Run inst/scripts/ingest_GSE56931.R first to produce ", humanFile)
seHuman <- readRDS(humanFile)
seHumanBaseline <- seHuman[, colData(seHuman)$protocol_phase == "baseline"]
cat(sprintf("GSE56931 baseline subset: %d samples, %d donors\n",
            ncol(seHumanBaseline), length(unique(colData(seHumanBaseline)$donor))))

pseudobulkFile <- "benchmarks/data/GSE145197_se.rds"
if (!file.exists(pseudobulkFile)) stop("Run inst/scripts/ingest_GSE145197.R first to produce ", pseudobulkFile)
sePseudobulk <- readRDS(pseudobulkFile)
cat(sprintf("GSE145197: %d samples, %d donors, ZT range %s\n",
            ncol(sePseudobulk), length(unique(colData(sePseudobulk)$donor)),
            paste(range(colData(sePseudobulk)$phase), collapse = "-")))

cat("\n############ Human (GSE56931, whole blood, baseline only) ############\n")
humanModels <- tryCatch(
  trainPreview(
    seHumanBaseline, species = "Hsapiens", platform = "microarray", tissueLabel = "Blood",
    citationNote = "GSE56931 (Arnardottir et al. 2014 Sleep), baseline samples only",
    nFolds = 5
  ),
  error = function(e) {
    message("Human preview FAILED: ", conditionMessage(e))
    list()
  }
)

cat("\n############ Pseudobulk (GSE145197, mouse liver hepatocytes) ############\n")
pseudobulkModels <- tryCatch(
  trainPreview(
    sePseudobulk, species = "Mmusculus", platform = "pseudobulk_scrnaseq", tissueLabel = "LivHep",
    citationNote = "GSE145197 (Droin et al. 2021 Nat Metab), pseudobulk hepatocytes",
    nFolds = 5
  ),
  error = function(e) {
    message("Pseudobulk preview FAILED: ", conditionMessage(e))
    list()
  }
)

previewModels <- c(humanModels, pseudobulkModels)
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(previewModels, "benchmarks/results/experimental_preview_reference_models.rds")
cat(sprintf(
  "\nSaved %d EXPERIMENTAL PREVIEW model(s) to benchmarks/results/experimental_preview_reference_models.rds (NOT R/sysdata.rda, NOT used by estimatePhase()).\n",
  length(previewModels)
))
