# S4 class definitions (AGENTS.md Section 5/10: S4 over S3 for the public API).
# The engines' own fit objects (molecularTimetableFit, tauFisherFit) stay S3 -
# they are internal implementation detail behind the fit()/predict() functions in
# R/engines-*.R, not part of the Section 6 public contract, and converting them
# would be a large, low-value refactor of already-tested code. PhaseReferenceModel
# below is the one new concept estimatePhase() actually introduces.

#' A trained circadian phase-inference reference model
#'
#' Bundles everything [estimatePhase()] needs to predict phase for new samples
#' against one (species, tissue, method, platform) reference cohort: the engine's
#' own fit object, a pool of real held-out circular residuals for
#' `phase_lo`/`phase_hi` (AGENTS.md Section 9 Phase 3 - the same residual-pool
#' mechanism [bootstrapPhaseCI()] uses, but precomputed once at training time
#' rather than live-bootstrapped per call, since refitting e.g. tauFisher dozens of
#' times per user call would be prohibitively slow), and the core-clock PCA basis +
#' derived threshold for `clock_amp`/`phase_flag` (AGENTS.md 7.4). Built by
#' `inst/scripts/train_reference_models.R` from GSE54650/GSE54651 and bundled
#' internally as `.referenceModels` (see [availableReferenceModels()] to list what
#' is covered).
#'
#' @slot species,tissue,method Character scalars identifying this model. Only
#'   `species = "Mmusculus"` is currently bundled - see [estimatePhase()]'s
#'   documentation for why (AGENTS.md 7.6 gap).
#' @slot platform Character scalar: the data type this model was trained on, e.g.
#'   `"microarray"` (GSE54650) or `"rnaseq"` (GSE54651, AGENTS.md Section 5 Step
#'   2.3). Defaults to `"microarray"` for models built before this slot existed,
#'   so old `new("PhaseReferenceModel", ...)` calls that omit it keep working.
#' @slot fit The engine's own fit object (e.g. from [fitMolecularTimetable()] or
#'   [fitTauFisher()]), passed as-is to that engine's `predict` function. Typed
#'   `ANY` rather than `list`: both fit objects are S3-classed lists
#'   (`"molecularTimetableFit"`/`"tauFisherFit"`), and S4's `validObject()`
#'   rejects an S3-classed object against a `list`-typed slot even though
#'   `is.list()` on it is `TRUE` - the engines' own `predict` functions require
#'   that S3 class via `inherits()`, so it must survive storage in this slot
#'   unchanged.
#' @slot residualPool Numeric vector of real held-out circular residuals
#'   (`circularDiff(pred, truth)`) pooled from tissues OTHER than this model's own
#'   (leave-tissue-out, same principle as `benchmarks/run_phase3_calibration.R`).
#' @slot ampCenter,ampScale Named numeric vectors: the per-gene z-score
#'   center/scale from [clockAmplitude()]'s PCA fit on this cohort.
#' @slot ampRotation Numeric matrix (genes x nPC): the PCA rotation from
#'   [clockAmplitude()] on this cohort.
#' @slot ampGenes Character vector, the core clock genes (row order matching
#'   `ampCenter`/`ampScale`/`ampRotation`) used for this cohort's amplitude PCA.
#' @slot ampThreshold Numeric scalar: the [clockAmplitude()] cutoff below which
#'   `phase_flag` becomes `"low_amplitude"` (from [deriveAmplitudeThreshold()]).
#' @slot cvMAE Numeric scalar: this exact (species, tissue, method, platform)
#'   combination's OWN held-out cross-validated circular MAE, in hours, measured
#'   during training. Distinct from `mean(abs(residualPool))`, which summarises
#'   *other* tissues' residuals (the pool is leave-tissue-out) and therefore
#'   cannot say how well this engine does on this particular tissue. That
#'   distinction is the whole point of the slot: [estimatePhase()]'s
#'   `method = "consensus"` weights each engine by `1/cvMAE`, which is only
#'   tissue-aware because this number is. `NA_real_` if not measured.
#' @slot trainedOn Character scalar, a citation/provenance string.
#' @return An S4 object of class `PhaseReferenceModel`.
#' @examples
#' rm <- methods::new("PhaseReferenceModel",
#'   species = "Mmusculus", tissue = "Liv", method = "timetable",
#'   fit = list(genes = c("Arntl", "Per1")), residualPool = rnorm(50, sd = 1.5),
#'   ampCenter = c(Arntl = 0, Per1 = 0), ampScale = c(Arntl = 1, Per1 = 1),
#'   ampRotation = matrix(c(0.7, 0.7, 0.7, -0.7), nrow = 2),
#'   ampGenes = c("Arntl", "Per1"), ampThreshold = 2,
#'   trainedOn = "toy example, not a real fit"
#' )
#' rm
#' @exportClass PhaseReferenceModel
methods::setClass("PhaseReferenceModel", representation = methods::representation(
  species = "character", tissue = "character", method = "character",
  platform = "character",
  fit = "ANY",
  residualPool = "numeric",
  ampCenter = "numeric", ampScale = "numeric", ampRotation = "matrix",
  ampGenes = "character", ampThreshold = "numeric",
  cvMAE = "numeric",
  trainedOn = "character"
), prototype = methods::prototype(platform = "microarray", cvMAE = NA_real_))
