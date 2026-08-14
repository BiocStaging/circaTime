# Main user entry point (AGENTS.md Section 6): SummarizedExperiment in, phased
# SummarizedExperiment out. Wires together pieces already built in Phases 1-3
# (engines in R/engines-*.R via the .circaTimeEngines registry in R/benchmark.R,
# clockAmplitude()/circularQuantileInterval() in R/confidence.R and
# R/circular-stats.R) against a bundled PhaseReferenceModel (R/AllClasses.R) rather
# than inventing new inference logic (AGENTS.md 0.3).

#' List the bundled reference models available to [estimatePhase()]
#'
#' @return A `data.frame`: `species`, `tissue`, `method`, `platform`, `cv_mae`
#'   (that model's own held-out cross-validated circular MAE in hours; `NA` if
#'   not measured), `n_amp_genes`, `n_residual_draws`. Rows whose `method` starts
#'   with `consensus_` are the combined models reachable via
#'   `estimatePhase(method = "consensus", consensus_weights = ...)`.
#' @details This list is the exhaustive set of models [estimatePhase()] can use.
#'   It is currently `species = "Mmusculus"` only, covering all 12 GSE54650
#'   tissues: Molecular Timetable and
#'   tauFisher, plus their `consensus_accuracy`/`consensus_equal` combinations, on
#'   `platform = "microarray"`, and Molecular Timetable only on
#'   `platform = "rnaseq"`. ZeitZeiger and TimeSignatR are wrapped engines but are
#'   deliberately absent here: they are reachable through [fitZeitZeiger()]/
#'   [predictZeitZeiger()], [fitTimeSignatR()]/[predictTimeSignatR()],
#'   [benchmarkPhase()], and [transferPhase()] rather than as pre-trained
#'   `estimatePhase()` models.
#' @examples
#' availableReferenceModels()
#' @export
availableReferenceModels <- function() {
  rows <- lapply(.referenceModels, function(rm) {
    data.frame(
      species = rm@species, tissue = rm@tissue, method = rm@method, platform = rm@platform,
      cv_mae = rm@cvMAE,
      n_amp_genes = length(rm@ampGenes), n_residual_draws = length(rm@residualPool)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Projects a *new* sample's core clock genes onto an already-fitted reference
# cohort's PCA basis (rather than re-fitting PCA on just this call's samples, as
# clockAmplitude() does), so a single-sample estimatePhase() call still gets a
# real amplitude instead of the NA that clockAmplitude()'s own <3-samples-per-
# group rule would produce. Genes in ampGenes missing from exprMat contribute 0 to
# the projection (same convention clockAmplitude() uses for zero-variance genes).
# Returns NA for every sample if fewer than 4 of ampGenes are present at all (below
# clockAmplitude()'s own minGenes default) - not enough to be a meaningful summary.
.projectClockAmplitude <- function(exprMat, refModel) {
  nSamples <- ncol(exprMat)
  matchIdx <- match(toupper(refModel@ampGenes), toupper(rownames(exprMat)))
  present <- !is.na(matchIdx)
  if (sum(present) < 4) {
    return(rep(NA_real_, nSamples))
  }

  block <- t(exprMat[matchIdx[present], , drop = FALSE]) # samples x genes(present)
  centered <- sweep(block, 2, refModel@ampCenter[present], "-")
  scaled <- sweep(centered, 2, refModel@ampScale[present], "/")
  scaled[is.na(scaled)] <- 0

  proj <- scaled %*% refModel@ampRotation[present, , drop = FALSE]
  sqrt(rowSums(proj^2))
}

# Runs one engine's predict() against one PhaseReferenceModel, including the
# gene-name normalisation and the gene-coverage pre-check that used to be inline
# in estimatePhase(). Factored out so method = "consensus" can call it once per
# constituent engine without duplicating any of that logic (AGENTS.md 0.3: the
# consensus path must reuse the single-engine path, not reimplement it).
.predictWithReferenceModel <- function(exprMat, refModel, species, map_aliases) {
  method <- refModel@method
  fitObj <- refModel@fit

  fittedGenes <- switch(method,
    timetable = fitObj$genes,
    taufisher = fitObj$Genes,
    stop("No gene-coverage rule defined for method '", method, "'")
  )
  targetGenes <- unique(c(fittedGenes, refModel@ampGenes))
  mapping <- .prepareExpressionForModel(
    exprMat, targetGenes = targetGenes, species = species, map_aliases = map_aliases
  )
  mappedMat <- mapping$exprMat
  overlap <- sum(toupper(fittedGenes) %in% toupper(rownames(mappedMat)))
  # Mirrors each predict function's own internal gating exactly (see
  # predictMolecularTimetable()'s default minGenes and predictTauFisher()'s
  # "every fitted gene must be present" rule), so this pre-check agrees with
  # what predict() is actually about to do.
  geneCoverageOK <- switch(method,
    timetable = overlap >= 5,
    taufisher = overlap == length(fittedGenes)
  )

  engine <- .circaTimeEngines[[method]]
  predDf <- do.call(engine$predict, list(fit = fitObj, testMat = mappedMat))

  list(
    phase = predDf$phase, geneCoverageOK = geneCoverageOK,
    exprMat = mappedMat, report = mapping$report
  )
}

# Combines several engines' per-sample phase estimates into one, per AGENTS.md
# Section 6 Step 3.1. `weights` is the stored per-engine prior (already
# tissue-specific - see PhaseReferenceModel@cvMAE); engines that failed their
# gene-coverage check, or that returned NA for a given sample, are dropped for
# that sample and the remaining weights renormalised, rather than being folded in
# as if they had contributed. Returns per-sample phase, the concentration
# (engine-agreement) measure from circularWeightedMean(), and how many engines
# actually contributed - none of which is inferred later from the phase alone.
.combineEnginePhases <- function(phaseMat, weights) {
  nSamples <- nrow(phaseMat)
  out <- data.frame(
    phase = rep(NA_real_, nSamples),
    concentration = rep(NA_real_, nSamples),
    n_engines = rep(0L, nSamples)
  )
  for (i in seq_len(nSamples)) {
    p <- phaseMat[i, ]
    usable <- !is.na(p) & weights > 0
    if (!any(usable)) next
    combined <- circularWeightedMean(p[usable], weights[usable])
    out$phase[i] <- as.numeric(combined)
    out$concentration[i] <- attr(combined, "concentration")
    out$n_engines[i] <- sum(usable)
  }
  out
}

#' Estimate circadian phase for new samples against a bundled reference model
#'
#' The package's main entry point (AGENTS.md Section 6). Looks up a pre-trained
#' [PhaseReferenceModel-class] for the requested `(species, tissue, method, platform)` -
#' see [availableReferenceModels()] for what is currently bundled - predicts phase
#' with that engine's own `predict` function, attaches a circular confidence
#' interval built from the reference model's real held-out residual pool (no live
#' bootstrap refit at call time - see [PhaseReferenceModel-class]), a core-clock
#' amplitude projected onto the reference cohort's PCA basis (works even for a
#' single-sample call, unlike a fresh [clockAmplitude()] fit), and a QC flag.
#'
#' Per the Section 6 contract, this never hard-errors on a per-sample basis - a
#' sample that can't be scored gets `phase = NA` and an informative `phase_flag`.
#' Hard errors are reserved for malformed input: an unrecognised
#' species/tissue/method combination, or an `se` missing gene rownames.
#'
#' @section What this does NOT do (yet):
#' Only `species = "Mmusculus"` is currently bundled, trained from GSE54650
#' (`platform = "microarray"`) and GSE54651 (`platform = "rnaseq"`, its RNA-seq
#' counterpart - see `inst/scripts/train_reference_models.R`) - calling with any
#' other species/tissue/method/platform combination is a hard error listing what
#' IS available, not a silent guess (AGENTS.md 7.5/0.3). `platform` is never
#' inferred from `se` itself; picking the wrong one does not error, it just scores
#' against a model trained on a different assay technology than your data. Gene matching is case-insensitive symbol matching,
#' the same convention [clockAmplitude()]/[coreClockGenes()] already use.
#' `data_species`/`map_orthologs` (AGENTS.md Section 4 Step 1.2) let `se` carry
#' genes in a *different* species' symbol convention than the bundled model - see
#' those arguments below - but this remains one-to-one symbol/ortholog matching,
#' not full Ensembl-ID mapping within a single species (see [standardizeGeneIds()]
#' for that). A larger multi-species/tissue reference cohort (the `circaTimeData`
#' ExperimentHub package, AGENTS.md Section 9 Phase 4 task 4) is future work, not
#' silently substituted here.
#'
#' @section Bayesian credible interval (Phase 3 Step 3.3):
#' Alongside the empirical `phase_lo`/`phase_hi` interval (residual-pool
#' bootstrap quantiles), every call also
#' returns `phase_credible_lo`/`phase_credible_hi`: a von Mises credible interval
#' at the same `level`, fit from the same reference model's residual pool via
#' [bayesianPhasePosterior()] rather than a second, independent computation. The
#' two intervals answer related but distinct questions - the empirical one is
#' "what interval did the observed residuals actually span", the Bayesian one is
#' "what is a smooth probability distribution over the true phase" - and are
#' expected to be close but not identical (see `?bayesianPhasePosterior`'s
#' Relationship section for when and why they can diverge).
#' `metadata(se)$circaTime_bayesian` records the fitted `kappa` and
#' `concentration`; call [bayesianPhasePosterior()] directly for the full
#' posterior density rather than just the interval endpoints.
#'
#' @section Consensus (`method = "consensus"`):
#' Rather than making the user pick one engine, `method = "consensus"` runs every
#' engine bundled for that (species, tissue, platform) and combines their phase
#' estimates on the circle via [circularWeightedMean()] (AGENTS.md Section 6
#' Step 3.1). Two weighting schemes are available via `consensus_weights`:
#'
#' * `"accuracy"` (default) - **tissue-aware**: each engine is weighted by
#'   `1/cvMAE`, where `cvMAE` is that engine's own held-out cross-validated
#'   circular MAE *on this specific tissue*, measured during training and stored
#'   on the model (see [PhaseReferenceModel-class]). An engine that is more
#'   accurate on liver than on brainstem is therefore weighted differently in
#'   each.
#' * `"equal"` - a plain circular mean of the engines' estimates, with no
#'   accuracy prior at all. Useful as the honest baseline for whether the
#'   accuracy weighting is buying anything.
#'
#' The weights are not hidden: `metadata(se)$circaTime_consensus` records the
#' engines used, the `stored_weights` (the trained prior) and the
#' `effective_weights` actually applied to this call, each engine's own
#' per-sample `per_engine_phase`, the per-sample `concentration` (how strongly
#' the engines agreed, in `[0, 1]`) and `n_engines_used`. An engine that fails
#' its gene-coverage check on your data has its weight set to zero and the
#' remainder renormalised, rather than contributing an unusable phase; the
#' consensus is flagged `"gene_coverage"` only if *no* engine had coverage.
#'
#' What consensus is available over is limited by what is bundled: mouse
#' microarray has both Molecular Timetable and tauFisher, so consensus there
#' combines two engines. ZeitZeiger and TimeSignature are deliberately not
#' bundled as pre-trained models (their fitted objects are two orders of
#' magnitude larger for no measured accuracy benefit), and the mouse RNA-seq
#' platform bundles Molecular Timetable only, so no consensus model exists for
#' it - requesting one is an error listing what does exist, not a silent
#' fallback to a single engine.
#'
#' @param se A `SummarizedExperiment` with an expression assay, genes x samples.
#' @param method Name of a bundled engine for the requested tissue/species (see
#'   [availableReferenceModels()]); currently `"timetable"` or `"taufisher"`, or
#'   `"consensus"` to combine every engine bundled for that
#'   tissue/species/platform (see the Consensus section below).
#' @param tissue Required, no default (AGENTS.md 7.5) - a tissue code from
#'   [availableReferenceModels()] for the requested species (GSE54650's mouse
#'   tissue abbreviations, e.g. `"Liv"` for liver - see
#'   `inst/scripts/ingest_GSE54650.R`).
#' @param species Default `"Mmusculus"` - species of the bundled reference model
#'   to use. See the section above.
#' @param platform Default `"microarray"` - the data type the bundled reference
#'   model was trained on (AGENTS.md Section 5 Step 2.3), e.g. `"microarray"`
#'   (GSE54650) or `"rnaseq"` (GSE54651, the RNA-seq counterpart of the same
#'   study/tissues). Not inferred from `se` or from `data_species`/`map_orthologs` -
#'   pick the platform whose reference model matches (or is closest to) how your
#'   own data was generated; see [availableReferenceModels()] for what is bundled.
#' @param consensus_weights Only used when `method = "consensus"`; one of
#'   `"accuracy"` (default) or `"equal"`. See the Consensus section below. These
#'   are two separately-trained models, not one model with a runtime switch: each
#'   carries its own held-out residual pool measured under that weighting, since
#'   an accuracy-weighted and an equally-weighted consensus do not have the same
#'   error distribution and must not share a confidence interval.
#' @param assay_name Name of the assay to use. If `NULL` (default), uses the
#'   assay named `"logExpr"` if present (every bundled reference model is
#'   trained on log-scale expression), else the first assay.
#' @param level Nominal coverage of `phase_lo`/`phase_hi` (default 0.8).
#' @param qc_method `"legacy"` (default) uses the fixed 16-gene core-clock
#'   amplitude ([clockAmplitude()]) projected onto the bundled reference model's
#'   training-cohort PCA basis, with the threshold fit once during training
#'   ([deriveAmplitudeThreshold()]). `"dynamic"` uses [clockCoherence()] instead:
#'   the informative gene panel is selected from `se` itself rather than a fixed
#'   panel (Phase 3 Step 3.2), which needs a real cohort of samples in this call
#'   (see `qc_min_samples`) rather than a stored basis, since the whole point is
#'   that the panel is dataset-specific. See [clockCoherence()] for why this
#'   scores self-consistency with the (usually estimated, not true) `phase`
#'   column, and `benchmarks/run_dynamic_qc_benchmark.R` for validation against
#'   the legacy flag on real labelled data.
#' @param qc_threshold Only used when `qc_method = "dynamic"`. Coherence values at
#'   or below this are flagged `"low_amplitude"`. Default `0`: a coherence of 0
#'   means the fitted cohort rhythm explains no more of this sample's panel than
#'   its own per-gene mean would - a model-free zero point, not fit to any one
#'   dataset. `benchmarks/run_dynamic_qc_benchmark.R` shows tightening this
#'   (e.g. to `0.1`) sharpens separation on GSE84511; left as a scale-appropriate
#'   default here rather than a per-tissue value snooped from one benchmark.
#' @param qc_min_samples Only used when `qc_method = "dynamic"`; minimum samples
#'   in `se` for [clockCoherence()]'s per-dataset panel selection to be
#'   meaningful (default 10, matching [clockCoherence()]'s own default). Below
#'   this, `qc_method` silently falls back to `"legacy"` for this call and
#'   `metadata(se)$circaTime_qc$method_used` records the fallback - a cohort this
#'   small cannot support fitting its own rhythm, so "dynamic" would just be
#'   noise dressed up as a QC score.
#' @param map_aliases Logical; if `TRUE` (default), normalize common clock-gene
#'   aliases (for example `BMAL1` -> `ARNTL`/`Arntl`) and strip Ensembl version
#'   suffixes before matching genes to the reference model. A mapping summary is
#'   stored in `metadata(se)$circaTime_gene_mapping`.
#' @param data_species Species of `rownames(se)` itself, default equal to
#'   `species` (the common case: data and model are the same species). Set this
#'   when `se`'s genes are in a *different* species' symbol convention than the
#'   requested model - e.g. human input data (`data_species = "Hsapiens"`)
#'   against the bundled mouse model (`species = "Mmusculus"`, the only one
#'   currently bundled).
#' @param map_orthologs Logical, default `TRUE` whenever `data_species` differs
#'   from `species`. Ortholog-maps `se`'s genes from `data_species` to `species`
#'   via [mapOrthologs()] (one-to-one only; genes with no unambiguous ortholog
#'   are dropped, not guessed at) before the usual gene-matching step. A mapping
#'   summary is stored in `metadata(se)$circaTime_ortholog_mapping`. Setting this
#'   `FALSE` while `data_species != species` is a hard error rather than silently
#'   comparing mismatched-species gene symbols.
#' @return `se` with `colData` columns added: `phase`, `phase_lo`, `phase_hi`
#'   (empirical residual-pool interval), `phase_credible_lo`, `phase_credible_hi`
#'   (the Phase 3 Step 3.3 von Mises credible interval built on the same
#'   residual pool - see [bayesianPhasePosterior()]; `metadata(se)$circaTime_bayesian`
#'   records the fitted concentration), `clock_amp`, `phase_flag` (factor:
#'   `"ok"` / `"low_amplitude"` / `"gene_coverage"` / `"failed"`). `clock_amp`
#'   holds whichever score `qc_method` produced (amplitude or coherence -
#'   `metadata(se)$circaTime_qc` records which). For `method = "consensus"`,
#'   `metadata(se)$circaTime_consensus` additionally records the per-engine
#'   estimates and the weights used (see the Consensus section).
#' @examples
#' mat <- matrix(rnorm(200), nrow = 100, dimnames = list(paste0("g", 1:100), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
#' # random genes won't overlap the bundled model's trained gene set - this shows
#' # the "never hard-error on a sample" contract: phase_flag is set, not an error.
#' se <- estimatePhase(se, method = "timetable", tissue = "Liv")
#' SummarizedExperiment::colData(se)$phase_flag
#' @export
estimatePhase <- function(se, method, tissue, species = "Mmusculus",
                           platform = "microarray",
                           consensus_weights = c("accuracy", "equal"),
                           assay_name = NULL, level = 0.8,
                           qc_method = c("legacy", "dynamic"),
                           qc_threshold = 0, qc_min_samples = 10,
                           map_aliases = TRUE,
                           data_species = species,
                           map_orthologs = !identical(data_species, species)) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment")
  }
  consensus_weights <- match.arg(consensus_weights)
  qc_method <- match.arg(qc_method)
  if (missing(tissue) || length(tissue) != 1 || is.na(tissue)) {
    stop("tissue is required and must be a single value (AGENTS.md 7.5 - ",
         "there is no default). See availableReferenceModels().")
  }
  if (missing(method) || length(method) != 1) {
    stop("method must be a single engine name. See availableReferenceModels().")
  }
  if (!identical(data_species, species) && !isTRUE(map_orthologs)) {
    stop(
      "data_species ('", data_species, "') differs from species ('", species,
      "') but map_orthologs = FALSE. Either set map_orthologs = TRUE (the ",
      "default whenever they differ) or align data_species with species.",
      call. = FALSE
    )
  }

  # "consensus" is a user-facing alias: the actual stored model is per weighting
  # scheme, because each scheme has its OWN real held-out residual pool measured
  # under that scheme (an accuracy-weighted consensus and an equally-weighted one
  # do not have the same error distribution, so they must not share a CI).
  isConsensus <- identical(method, "consensus")
  lookupMethod <- if (isConsensus) paste0("consensus_", consensus_weights) else method

  key <- paste(species, tissue, lookupMethod, platform, sep = ".")
  refModel <- .referenceModels[[key]]
  if (is.null(refModel)) {
    avail <- availableReferenceModels()
    stop(
      "No bundled reference model for species = '", species, "', tissue = '",
      tissue, "', method = '", lookupMethod, "', platform = '", platform,
      "'. Available combinations:\n",
      paste(
        utils::capture.output(utils::write.table(avail, row.names = FALSE, quote = FALSE)),
        collapse = "\n"
      )
    )
  }

  exprMat <- if (is.null(assay_name)) {
    # Every engine is trained on log-scale expression (see fitMolecularTimetable()/
    # fitTauFisher() and inst/scripts/train_reference_models.R); prefer an assay
    # literally named "logExpr" when present (the convention this package's own
    # ingest scripts use) over silently taking whichever assay happens to be
    # first - a raw-scale first assay produces materially worse predictions
    # without erroring, which is a real trap for an SE built the same way
    # inst/scripts/ingest_GSE54650.R builds one (assays = list(expr = ..., logExpr = ...)).
    if ("logExpr" %in% SummarizedExperiment::assayNames(se)) {
      SummarizedExperiment::assay(se, "logExpr")
    } else {
      SummarizedExperiment::assay(se)
    }
  } else {
    SummarizedExperiment::assay(se, assay_name)
  }
  if (is.null(rownames(exprMat))) {
    stop("se's assay must have rownames (gene identifiers)")
  }

  orthoReport <- NULL
  if (!identical(data_species, species)) {
    orthoReport <- mapOrthologs(rownames(exprMat), fromSpecies = data_species, toSpecies = species)
    idx <- match(rownames(exprMat), orthoReport$input_gene)
    keepRow <- !is.na(idx) & orthoReport$status[idx] == "ok"
    exprMat <- exprMat[keepRow, , drop = FALSE]
    rownames(exprMat) <- orthoReport$ortholog[idx[keepRow]]
  }

  nSamples <- ncol(se)
  consensusReport <- NULL

  if (isConsensus) {
    # Each constituent engine gets its own gene-name normalisation pass (their
    # target gene sets differ), so predict per engine, then combine.
    engines <- refModel@fit$engines
    storedWeights <- refModel@fit$weights
    if (!all(engines %in% names(storedWeights))) {
      stop(
        "Consensus model '", key, "' stores weights for (",
        paste(names(storedWeights), collapse = ", "), ") but lists engines (",
        paste(engines, collapse = ", "), "). The stored model is malformed - ",
        "rebuild it with inst/scripts/train_reference_models.R.",
        call. = FALSE
      )
    }
    perEngine <- lapply(engines, function(e) {
      constituentKey <- paste(species, tissue, e, platform, sep = ".")
      constituent <- .referenceModels[[constituentKey]]
      if (is.null(constituent)) {
        stop(
          "Consensus model for tissue '", tissue, "' lists engine '", e,
          "' but no bundled model '", constituentKey, "' exists. The bundled ",
          "reference models and consensus definitions are out of sync - ",
          "rebuild them with inst/scripts/train_reference_models.R.",
          call. = FALSE
        )
      }
      .predictWithReferenceModel(exprMat, constituent, species, map_aliases)
    })
    names(perEngine) <- engines

    phaseMat <- do.call(cbind, lapply(perEngine, `[[`, "phase"))
    colnames(phaseMat) <- engines
    coverageOK <- vapply(perEngine, `[[`, logical(1), "geneCoverageOK")

    # An engine whose genes are not adequately covered contributes nothing rather
    # than contributing a garbage phase; its weight goes to zero and the rest are
    # renormalised inside circularWeightedMean().
    effectiveWeights <- storedWeights[engines]
    effectiveWeights[!coverageOK] <- 0

    combined <- .combineEnginePhases(phaseMat, effectiveWeights)
    point <- combined$phase
    geneCoverageOK <- any(coverageOK)

    # Amplitude basis is identical across a tissue's engines (one clockAmplitude()
    # fit per tissue in training), so the consensus model's own copy is the same
    # basis the constituents carry - use it directly.
    legacyAmpVec <- .projectClockAmplitude(perEngine[[1]]$exprMat, refModel)
    qcExprMat <- perEngine[[1]]$exprMat

    mapping <- perEngine[[1]]
    consensusReport <- list(
      engines = engines,
      stored_weights = storedWeights[engines],
      effective_weights = effectiveWeights / sum(effectiveWeights),
      gene_coverage_ok = coverageOK,
      per_engine_phase = phaseMat,
      concentration = combined$concentration,
      n_engines_used = combined$n_engines
    )
  } else {
    single <- .predictWithReferenceModel(exprMat, refModel, species, map_aliases)
    point <- single$phase
    geneCoverageOK <- single$geneCoverageOK
    legacyAmpVec <- .projectClockAmplitude(single$exprMat, refModel)
    qcExprMat <- single$exprMat
    mapping <- single
  }

  # Dynamic QC needs a real cohort in THIS call to select its own gene panel
  # (that is the whole point of Step 3.2 - see ?clockCoherence); below
  # qc_min_samples it cannot do that meaningfully, so it silently falls back to
  # the legacy fixed-panel flag rather than returning a noisy score.
  qcMethodUsed <- qc_method
  if (identical(qc_method, "dynamic") && nSamples >= qc_min_samples) {
    ampVec <- clockCoherence(qcExprMat, phase = point, minSamples = qc_min_samples)
    ampThreshold <- qc_threshold
    qcAttrs <- list(panel = attr(ampVec, "panel"), panelSize = attr(ampVec, "panelSize"),
                    relaxed = attr(ampVec, "relaxed"), r2Cutoff = attr(ampVec, "r2Cutoff"))
  } else {
    if (identical(qc_method, "dynamic")) qcMethodUsed <- "legacy"
    ampVec <- legacyAmpVec
    ampThreshold <- refModel@ampThreshold
    qcAttrs <- NULL
  }

  ciList <- lapply(seq_len(nSamples), function(i) {
    if (is.na(point[i])) {
      return(list(lo = NA_real_, hi = NA_real_))
    }
    circularQuantileInterval(wrapPhase(point[i] + refModel@residualPool), level = level)
  })
  loVec <- vapply(ciList, `[[`, numeric(1), "lo")
  hiVec <- vapply(ciList, `[[`, numeric(1), "hi")

  # Phase 3 Step 3.3: a parametric (von Mises) posterior built on the SAME
  # residual pool as the empirical interval above, not a second live bootstrap -
  # see ?bayesianPhasePosterior. kappa is a single number fit once from the pool
  # and shared by every sample in this call, so this is cheap regardless of
  # nSamples; only the per-sample credible bounds (a shift of that one
  # distribution by each sample's own point estimate) vary per sample.
  bayesPost <- if (length(refModel@residualPool) >= 2) {
    bayesianPhasePosterior(
      point = ifelse(is.na(point), 0, point), # placeholder mu, masked out below
      residualPool = refModel@residualPool, level = level
    )
  } else {
    NULL
  }
  credLoVec <- if (is.null(bayesPost)) rep(NA_real_, nSamples) else bayesPost$credible_lo
  credHiVec <- if (is.null(bayesPost)) rep(NA_real_, nSamples) else bayesPost$credible_hi
  credLoVec[is.na(point)] <- NA_real_
  credHiVec[is.na(point)] <- NA_real_

  flag <- rep("ok", nSamples)
  if (!geneCoverageOK) {
    flag[] <- "gene_coverage"
  } else {
    flag[is.na(point)] <- "failed"
    lowAmpIdx <- flag == "ok" & !is.na(ampVec) & ampVec <= ampThreshold
    flag[lowAmpIdx] <- "low_amplitude"
  }
  flag <- factor(flag, levels = c("ok", "low_amplitude", "gene_coverage", "failed"))

  cd <- SummarizedExperiment::colData(se)
  cd$phase <- point
  cd$phase_lo <- loVec
  cd$phase_hi <- hiVec
  cd$phase_credible_lo <- credLoVec
  cd$phase_credible_hi <- credHiVec
  cd$clock_amp <- ampVec
  cd$phase_flag <- flag
  SummarizedExperiment::colData(se) <- cd
  oldMeta <- S4Vectors::metadata(se)
  oldMeta$circaTime_gene_mapping <- mapping$report
  oldMeta$circaTime_qc <- c(
    list(method_requested = qc_method, method_used = qcMethodUsed, threshold = ampThreshold),
    qcAttrs
  )
  oldMeta$circaTime_bayesian <- list(
    kappa = if (is.null(bayesPost)) NA_real_ else bayesPost$kappa,
    concentration = if (is.null(bayesPost)) NA_real_ else bayesPost$concentration,
    level = level
  )
  if (!is.null(orthoReport)) oldMeta$circaTime_ortholog_mapping <- orthoReport
  if (!is.null(consensusReport)) oldMeta$circaTime_consensus <- consensusReport
  S4Vectors::metadata(se) <- oldMeta

  se
}
