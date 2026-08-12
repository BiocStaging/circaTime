# Data-prep script (AGENTS.md Section 5) that builds R/sysdata.rda, the bundled
# `.referenceModels` used by estimatePhase() (R/estimatePhase.R). Unlike the
# ingest_*.R scripts, whose *output* is gitignored benchmark data, this script's
# output IS shipped package data - rerun it and commit the resulting R/sysdata.rda
# whenever the training data or engine settings change.
#
# AGENTS.md Section 5 Step 2.3: trains one model set PER PLATFORM (microarray from
# GSE54650, rnaseq from GSE54651 - its same-study, same-tissue RNA-seq counterpart,
# notes/decisions.md Section 10), keyed as species.tissue.method.platform
# (PhaseReferenceModel@platform, R/AllClasses.R) so estimatePhase(platform = ...)
# can pick the model matching a user's own data's assay technology explicitly,
# rather than always silently falling back to the microarray model. Both platforms
# use IDENTICAL engine settings and the identical two-pass procedure below, run
# once per platform.
#
# Trains both engines (timetable, tauFisher) on all 12 tissues, using the exact
# settings already used/validated in benchmarks/run_benchmark.R (Phase 2 gate) and
# benchmarks/run_phase3_calibration.R (Phase 3 gate):
#   - Full per-tissue fit (all 24 samples) -> the model estimatePhase() actually
#     calls predict() on.
#   - A plain per-tissue 5-fold CV pass (via the existing, tested benchmarkPhase(),
#     not a bespoke loop) -> real held-out circular residuals, pooled leave-tissue-
#     out (a tissue's own reference model never uses its own residuals) for
#     phase_lo/phase_hi, same principle as run_phase3_calibration.R.
#   - clockAmplitude() per tissue (single group, matching how run_phase3_calibration.R
#     computed it) -> the PCA center/scale/rotation estimatePhase() projects new
#     samples onto (see .projectClockAmplitude() in R/estimatePhase.R).
#   - Per-(tissue, method) held-out CV circular MAE -> PhaseReferenceModel@cvMAE,
#     which is what makes method = "consensus"'s accuracy weighting tissue-aware
#     (AGENTS.md Section 6 Step 3.1).
#   - Consensus entries (one per weighting scheme: consensus_accuracy,
#     consensus_equal) whenever a platform bundles >= 2 engines. Their residual
#     pools are REAL held-out consensus predictions: the weighting rule is applied
#     to the same out-of-fold per-engine predictions Pass 1 already produced, not
#     to a refit. See the trainedOn string for the one honest caveat (the accuracy
#     weights themselves come from that same CV run's aggregate MAE).
#   - The amplitude threshold is NOT re-derived here - it reuses the value already
#     produced and honestly documented in Phase 3
#     (benchmarks/results/phase3_amplitude_threshold_GSE54650_timetable.csv, see
#     notes/decisions.md for the p=0.083 caveat) rather than
#     silently recomputing a fresh number (AGENTS.md 0.3).
#
# Run from the project root: Rscript inst/scripts/train_reference_models.R
# Requires benchmarks/data/GSE54650_se.rds (run inst/scripts/ingest_GSE54650.R
# first) and benchmarks/results/phase3_amplitude_threshold_GSE54650_timetable.csv
# (run benchmarks/run_phase3_calibration.R first, or use the version already
# committed to that path).

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
cat(sprintf("Reusing Phase 3's derived amplitude threshold: %.4f\n", ampThreshold))
# NOT re-derived per platform (AGENTS.md 0.3: this is a stated simplifying
# assumption, not silently presented as freshly computed for GSE54651/GSE56931's
# own clock_amp scale) - see notes/decisions.md Section 12.

engineArgs <- list(
  timetable = list(fit = list(nGenes = 100, minCor = 0.5)),
  taufisher = list(fit = list(method = c("JTK", "LS"), thres = 1.01, nrep = 1, numbasis = 5))
)
nFolds <- 5
seed <- 1

# Trains one platform's full model set (both passes) and returns a named list of
# PhaseReferenceModel, keyed species.tissue.method.platform - the unit
# inst/scripts/train_reference_models.R produces per call, one call per platform.
# `methods` lets a platform bundle fewer engines than the default two - see the
# GSE54651/rnaseq call below for why.
trainPlatform <- function(se, species, platform, citationNote, methods = c("timetable", "taufisher")) {
  tissues <- unique(colData(se)$tissue)

  # ---- Pass 1: plain per-tissue CV (real benchmarkPhase(), not a bespoke loop) to
  # build real held-out residual pools, pooled leave-tissue-out below.
  cat(sprintf("== [%s/%s] Pass 1: cross-validated residuals per tissue ==\n", species, platform))
  pass1 <- list()
  for (tis in tissues) {
    cat(sprintf("[%s] CV tissue=%s\n", format(Sys.time()), tis))
    seT <- se[, colData(se)$tissue == tis]
    bm <- benchmarkPhase(seT,
      truth_col = "phase", methods = methods,
      donor_col = "donor", nFolds = nFolds, seed = seed, engineArgs = engineArgs[methods]
    )
    res <- bm$residuals
    res <- res[!res$fold_failed, ]
    res$tissue <- tis
    res$residual <- circularDiff(res$pred, res$truth)
    pass1[[tis]] <- res
  }
  pass1All <- do.call(rbind, pass1)
  rownames(pass1All) <- NULL

  # Per-(tissue, method) held-out CV circular MAE -> PhaseReferenceModel@cvMAE,
  # and the tissue-aware prior behind method = "consensus" (AGENTS.md Section 6
  # Step 3.1). This is deliberately NOT mean(abs(residualPool)): the residual
  # pool is leave-tissue-out, so it describes other tissues, whereas the
  # consensus weighting needs this tissue's own number.
  # benchmarkPhase() already computed circular_error per held-out sample; reuse it
  # rather than recomputing the same quantity from pred/truth.
  cvMAEByTissueMethod <- stats::aggregate(
    circular_error ~ tissue + method, data = pass1All,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  lookupCvMAE <- function(tis, m) {
    hit <- cvMAEByTissueMethod$tissue == tis & cvMAEByTissueMethod$method == m
    if (!any(hit)) NA_real_ else cvMAEByTissueMethod$circular_error[hit][1]
  }

  # Real held-out consensus residuals: apply each weighting scheme to the SAME
  # per-fold held-out predictions the individual engines produced above, so the
  # consensus CI rests on genuine out-of-fold predictions rather than a refit.
  # Honest caveat recorded in trainedOn: the accuracy weights themselves are
  # derived from this same CV run's aggregate MAE, so that one hyperparameter is
  # not itself held out - each individual prediction being combined still is.
  consensusWeightsFor <- function(tis, scheme, ms) {
    if (identical(scheme, "equal")) {
      w <- rep(1, length(ms))
    } else {
      maes <- vapply(ms, function(m) lookupCvMAE(tis, m), numeric(1))
      if (anyNA(maes) || any(maes <= 0)) return(NULL)
      w <- 1 / maes
    }
    stats::setNames(w / sum(w), ms)
  }

  consensusResiduals <- list() # scheme -> data.frame(tissue, residual)
  if (length(methods) >= 2) {
    for (scheme in c("accuracy", "equal")) {
      perTissue <- list()
      for (tis in tissues) {
        sub <- pass1All[pass1All$tissue == tis, ]
        w <- consensusWeightsFor(tis, scheme, methods)
        if (is.null(w)) next
        wide <- stats::reshape(
          sub[, c("sample", "method", "pred", "truth")],
          idvar = c("sample", "truth"), timevar = "method", direction = "wide"
        )
        predCols <- paste0("pred.", methods)
        if (!all(predCols %in% colnames(wide))) next
        phaseMat <- as.matrix(wide[, predCols, drop = FALSE])
        combined <- .combineEnginePhases(phaseMat, unname(w[methods]))
        keep <- !is.na(combined$phase)
        if (!any(keep)) next
        perTissue[[tis]] <- data.frame(
          tissue = tis,
          residual = circularDiff(combined$phase[keep], wide$truth[keep])
        )
      }
      if (length(perTissue) > 0) {
        consensusResiduals[[scheme]] <- do.call(rbind, perTissue)
      }
    }
  }

  # ---- Pass 2: full per-tissue fit (what estimatePhase() actually predicts with),
  # plus clockAmplitude()'s PCA basis for that tissue.
  cat(sprintf("\n== [%s/%s] Pass 2: full per-tissue fits + amplitude PCA basis ==\n", species, platform))
  referenceModels <- list()
  for (tis in tissues) {
    cat(sprintf("[%s] fitting tissue=%s\n", format(Sys.time()), tis))
    seT <- se[, colData(se)$tissue == tis]
    exprMat <- assay(seT, "logExpr")
    truth <- colData(seT)$phase

    amp <- clockAmplitude(exprMat) # single group ("all") - matches run_phase3_calibration.R's usage
    ampGenes <- attr(amp, "genesUsed")
    ampCenter <- attr(amp, "center")[["all"]]
    ampScale <- attr(amp, "scale")[["all"]]
    ampRotation <- attr(amp, "rotation")[["all"]]

    fits <- list()
    if ("timetable" %in% methods) {
      fits$timetable <- do.call(fitMolecularTimetable, c(list(trainMat = exprMat, trainTime = truth), engineArgs$timetable$fit))
    }
    if ("taufisher" %in% methods) {
      fits$taufisher <- do.call(fitTauFisher, c(list(trainMat = exprMat, trainTime = truth, seed = seed), engineArgs$taufisher$fit))

      # tauFisher's fit$Model is an nnet::multinom object whose $terms carries a
      # reference to the R environment train_tauFisher() fit it in - which, via R's
      # normal lexical scoping, transitively holds that call's entire local training
      # data (train_df: the full genes x samples matrix, several MB per tissue, not
      # just the periodic genes actually used). Left attached, this alone blew
      # R/sysdata.rda past the 5 MB in-package limit even though every genuinely-
      # needed piece of the fit is tiny. Stripping the environment (and the now-
      # dangling $call) does not change predictions - verified identical before/after
      # on real GSE54650 data - since predictTauFisher() supplies newdata directly
      # and nnet's predict.multinom() does not need to re-evaluate the original
      # training environment.
      environment(fits$taufisher$Model$terms) <- baseenv()
      fits$taufisher$Model$call <- NULL
    }

    for (m in names(fits)) {
      residualPool <- pass1All$residual[pass1All$tissue != tis & pass1All$method == m]
      key <- paste(species, tis, m, platform, sep = ".")
      referenceModels[[key]] <- methods::new("PhaseReferenceModel",
        species = species, tissue = tis, method = m, platform = platform,
        fit = fits[[m]], # keep its molecularTimetableFit/tauFisherFit S3 class -
        # predictMolecularTimetable()/predictTauFisher() require it via inherits()
        residualPool = residualPool,
        ampCenter = ampCenter, ampScale = ampScale, ampRotation = ampRotation,
        ampGenes = ampGenes, ampThreshold = ampThreshold,
        cvMAE = lookupCvMAE(tis, m),
        trainedOn = paste0(
          citationNote, ", tissue=", tis, ", platform=", platform,
          "; residual pool leave-tissue-out from other ", length(tissues) - 1, " tissues (n=",
          length(residualPool), "); amplitude threshold from Phase 3's GSE54650/",
          "timetable pooled derivation, reused across platforms (not re-derived) - ",
          "see notes/decisions.md"
        )
      )
    }

    # Consensus entries (AGENTS.md Section 6 Step 3.1), one per weighting scheme.
    # These carry no engine fit of their own - `fit` instead records which engines
    # to combine and with what weights, so estimatePhase() applies exactly the
    # combination whose residuals were measured here, and a reviewer can read the
    # weights straight off the stored model.
    for (scheme in names(consensusResiduals)) {
      w <- consensusWeightsFor(tis, scheme, methods)
      if (is.null(w)) next
      pool <- consensusResiduals[[scheme]]
      residualPool <- pool$residual[pool$tissue != tis]
      if (length(residualPool) == 0) next
      key <- paste(species, tis, paste0("consensus_", scheme), platform, sep = ".")
      referenceModels[[key]] <- methods::new("PhaseReferenceModel",
        species = species, tissue = tis, method = paste0("consensus_", scheme),
        platform = platform,
        fit = list(engines = methods, weights = w),
        residualPool = residualPool,
        ampCenter = ampCenter, ampScale = ampScale, ampRotation = ampRotation,
        ampGenes = ampGenes, ampThreshold = ampThreshold,
        cvMAE = {
          own <- pool$residual[pool$tissue == tis]
          if (length(own) == 0) NA_real_ else mean(abs(own))
        },
        trainedOn = paste0(
          citationNote, ", tissue=", tis, ", platform=", platform,
          "; CONSENSUS of ", paste(methods, collapse = "+"), " weighted by '", scheme,
          "' (", paste(sprintf("%s=%.3f", names(w), w), collapse = ", "), ")",
          "; residual pool is real held-out consensus predictions, leave-tissue-out ",
          "from other ", length(tissues) - 1, " tissues (n=", length(residualPool), ")",
          if (identical(scheme, "accuracy")) {
            paste0(
              ". Caveat: the accuracy weights are derived from the same ", nFolds,
              "-fold CV run whose out-of-fold predictions form this pool, so the ",
              "weighting hyperparameter itself is not held out (each combined ",
              "prediction still is) - see notes/decisions.md"
            )
          } else {
            ""
          }
        )
      )
    }
  }
  referenceModels
}

microarrayFile <- "benchmarks/data/GSE54650_se.rds"
if (!file.exists(microarrayFile)) {
  stop("Run inst/scripts/ingest_GSE54650.R first to produce ", microarrayFile)
}
rnaseqFile <- "benchmarks/data/GSE54651_se.rds"
if (!file.exists(rnaseqFile)) {
  stop("Run inst/scripts/ingest_GSE54651.R first to produce ", rnaseqFile)
}

microarrayModels <- trainPlatform(
  readRDS(microarrayFile), species = "Mmusculus", platform = "microarray",
  citationNote = "GSE54650 (Zhang et al. 2014 PNAS)"
)
# tauFisher excluded here: MetaCycle's rhythmicity detection at the same thres =
# 1.01 setting used for microarray selects 65-373 "periodic" genes per tissue on
# this RNA-seq dataset (vs. ~15 on GSE54650), inflating tauFisher's pairwise-
# gene-difference PCA_embedding to 30-47 MB PER TISSUE even after the
# environment-stripping fix above - verified directly (notes/decisions.md Section
# 12), not assumed. That is a property of MetaCycle's detection behaving
# differently on this platform/dataset, not a bug in the stripping fix or in
# fitTauFisher() itself, and re-tuning thres specifically for this platform would
# break the "identical settings across platforms" comparison the whole point of
# this benchmark depends on. Molecular Timetable alone is bundled for rnaseq;
# tauFisher on RNA-seq data remains usable via fitTauFisher()/benchmarkPhase()
# directly, just not as a bundled pre-trained model.
rnaseqModels <- trainPlatform(
  readRDS(rnaseqFile), species = "Mmusculus", platform = "rnaseq",
  citationNote = "GSE54651 (Zhang et al. 2014 PNAS, RNA-seq SubSeries of GSE54650's study)",
  methods = "timetable"
)

referenceModels <- c(microarrayModels, rnaseqModels)
.referenceModels <- referenceModels
save(.referenceModels, file = "R/sysdata.rda", compress = "xz")
sizeKb <- file.size("R/sysdata.rda") / 1024
cat(sprintf(
  "\nSaved R/sysdata.rda: %d reference models, %.1f KB\n",
  length(referenceModels), sizeKb
))
if (sizeKb > 5 * 1024) {
  stop("R/sysdata.rda exceeds 5 MB (AGENTS.md Section 9 Phase 4 task 4) - stop and ",
       "report rather than silently shipping an oversized package (AGENTS.md Section 4).")
}
