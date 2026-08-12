# Phase 3 uncertainty layer (AGENTS.md Section 9). Three pieces, in the order
# AGENTS.md lists them:
#   1. clockAmplitude()        - per-sample core-clock amplitude (7.4), the fixed
#                                 16-gene "legacy" QC flag.
#   1b. clockCoherence()        - Phase 3 Step 3.2's dynamic replacement for it:
#                                 the informative gene panel is selected from the
#                                 input dataset itself, not from a fixed list.
#   2. bootstrapPhaseCI()       - phase_lo/phase_hi via bootstrap refitting, engine-
#                                 agnostic (wraps whichever engine is named, does not
#                                 invent a new inference method - Section 0.3).
#   3. calibrationCurve()/plotCalibration() - does the stated coverage hold on held-out
#                                 labelled data.
#   4. deriveAmplitudeThreshold() - the low_amplitude cutoff, fit empirically per the
#                                 Phase 3 gate, not guessed.
#   5. bayesianPhasePosterior()  - Phase 3 Step 3.3: a von Mises posterior density +
#                                 credible interval built directly from the same
#                                 residual pool bootstrapPhaseCI()/estimatePhase()
#                                 already use, not a new/independent inference method.
# All circular maths (quantile intervals, coverage) lives in R/circular-stats.R and is
# only called into here, per AGENTS.md 7.1.

#' The bundled core clock gene reference table
#'
#' Returns the AGENTS.md 7.7 core clock gene set, one row per gene, with common
#' aliases and a source citation. Matching against expression data rownames is
#' always case-insensitive (see [clockAmplitude()]) - the mouse orthologs of every
#' gene in this set share the same letters as the human symbol, differing only in
#' case, so a single table serves both species without a separate mouse column.
#'
#' @return A `data.frame` with columns `symbol` (canonical human HGNC symbol),
#'   `aliases` (semicolon-separated), `source` (citation).
#' @examples
#' head(coreClockGenes())
#' @export
coreClockGenes <- function() {
  path <- system.file("extdata", "core_clock_genes.tsv", package = "circaTime")
  if (!nzchar(path)) {
    stop("core_clock_genes.tsv not found; is circaTime installed correctly?")
  }
  utils::read.delim(path, stringsAsFactors = FALSE)
}

#' Per-sample core-clock amplitude
#'
#' Tissue-normalised amplitude of the AGENTS.md 7.7 core clock gene set, used to
#' flag samples where the circadian clock is too damped for a phase estimate to be
#' meaningful (AGENTS.md 7.4). Each matched gene is z-scored within its tissue
#' group (so amplitude is comparable across tissues with different baseline
#' expression/variance), then the z-scored block is reduced with PCA; amplitude is
#' the Euclidean norm of a sample's first two PC scores. This mirrors the
#' "eigengene" amplitude used in the CYCLOPS/circadian-transcriptomics literature
#' (Anafi et al. 2017) for the same QC purpose: PC1/PC2 capture the near-circular
#' manifold that a co-oscillating gene set traces out across a cohort, so their
#' combined magnitude reflects how strongly a given sample's clock genes
#' co-oscillate, independent of what phase they point to. This is a QC/amplitude
#' summary, not a new phase-inference algorithm (AGENTS.md 0.3).
#'
#' @section Known limitation — validated on real data, effect not yet confirmed:
#' AGENTS.md 7.4 expects this amplitude to be low "in aged tissue, tumours, and
#' critically ill donors," and Phase 3's gate requires showing flagged samples have
#' demonstrably worse phase error. On the two labelled real datasets tested so far
#' this effect was **not** statistically significant: GSE54650 (12-tissue healthy
#' mouse atlas, no aged/damped-clock samples; p=0.083) and GSE84511 (young vs aged
#' mouse epidermal stem cells, real ZT ground truth; p=0.20 for the derived
#' threshold, and mean amplitude barely differs by age group, p=0.74 - see
#' `benchmarks/results/phase3_amplitude_threshold_GSE84511_timetable.csv` and
#' `notes/decisions.md`). The likely reason, per GSE84511's source
#' paper (Solanas et al. 2017 Cell) and its companion liver study (Sato et
#' al. 2017 Cell): aging leaves core clock GENE expression population-average
#' "virtually unaltered" and instead reprograms which downstream genes oscillate -
#' so a core-clock-gene PCA amplitude may not be the right proxy for the kind of
#' age-related damping the literature actually reports, even though aged samples
#' in GSE84511 did have ~2.5x higher real cross-validated phase error than young
#' ones. Kept as documented, honest scope for now rather than silently dropped or
#' rounded up (AGENTS.md 0.3); revisit if a dataset or amplitude definition
#' surfaces that separates the two groups.
#'
#' @param exprMat Genes x samples numeric matrix, log-scale, with gene symbol
#'   rownames (any case; matched case-insensitively against [coreClockGenes()]).
#' @param tissue Optional character/factor vector, one tissue label per sample
#'   (`ncol(exprMat)` long). If supplied, genes are z-scored within each tissue
#'   group separately before PCA. If `NULL` (default), all samples are treated as
#'   one group.
#' @param minGenes Minimum number of matched core clock genes required (default
#'   4); errors below this since amplitude from too few genes is not a meaningful
#'   summary.
#' @return A numeric vector, `ncol(exprMat)` long, of non-negative amplitude
#'   values, with attributes `genesUsed` (character vector of matched rownames),
#'   `varExplained` (named numeric, fraction of variance captured by PC1+PC2 per
#'   tissue group), and `center`/`scale`/`rotation` (named lists, one element per
#'   tissue group: the [scale()] center/scale vectors and [prcomp()] rotation
#'   matrix used for that group's PCA) - the latter three let a caller reproduce
#'   this exact projection for a *new* sample against an already-fitted reference
#'   group, which is what [estimatePhase()] does for small (often n=1) inference
#'   calls where there aren't enough samples in the call itself to fit a fresh PCA.
#'   `NA` for samples in a tissue group with fewer than 3 members (too few for
#'   PCA/z-scoring to be meaningful).
#' @examples
#' set.seed(1)
#' genes <- c("Arntl", "Per1", "Per2", "Cry1", "Nr1d1", "Dbp")
#' mat <- matrix(rnorm(6 * 20), nrow = 6, dimnames = list(genes, paste0("s", 1:20)))
#' amp <- clockAmplitude(mat)
#' @export
clockAmplitude <- function(exprMat, tissue = NULL, minGenes = 4) {
  stopifnot(is.matrix(exprMat))
  if (is.null(rownames(exprMat))) {
    stop("exprMat must have rownames (gene identifiers)")
  }
  nSamples <- ncol(exprMat)
  if (!is.null(tissue) && length(tissue) != nSamples) {
    stop("tissue must be NULL or length ncol(exprMat)")
  }
  grp <- if (is.null(tissue)) rep("all", nSamples) else as.character(tissue)

  clockGenes <- coreClockGenes()$symbol
  matchIdx <- match(toupper(clockGenes), toupper(rownames(exprMat)))
  genesUsed <- rownames(exprMat)[matchIdx[!is.na(matchIdx)]]
  if (length(genesUsed) < minGenes) {
    stop("Only ", length(genesUsed), " of ", length(clockGenes),
         " core clock genes found in exprMat rownames (need >= ", minGenes, ")")
  }

  sub <- exprMat[genesUsed, , drop = FALSE]
  amp <- rep(NA_real_, nSamples)
  groups <- unique(grp)
  varExplained <- stats::setNames(rep(NA_real_, length(groups)), groups)
  center <- stats::setNames(vector("list", length(groups)), groups)
  scaleVec <- stats::setNames(vector("list", length(groups)), groups)
  rotation <- stats::setNames(vector("list", length(groups)), groups)

  for (g in groups) {
    idx <- which(grp == g)
    if (length(idx) < 3) next # too few samples in this group for PCA to be meaningful

    block <- t(sub[, idx, drop = FALSE]) # samples x genes
    z <- scale(block) # z-score each gene within this tissue group
    z[is.na(z)] <- 0 # zero-variance gene in this group carries no signal, not missingness

    pc <- stats::prcomp(z, center = FALSE, scale. = FALSE)
    nPC <- min(2, ncol(pc$x))
    amp[idx] <- sqrt(rowSums(pc$x[, seq_len(nPC), drop = FALSE]^2))
    varExplained[g] <- sum(pc$sdev[seq_len(nPC)]^2) / sum(pc$sdev^2)

    grpCenter <- attr(z, "scaled:center"); grpScale <- attr(z, "scaled:scale")
    grpScale[is.na(grpScale) | grpScale == 0] <- 1 # matches the z[is.na(z)] <- 0 fallback above
    center[[g]] <- grpCenter
    scaleVec[[g]] <- grpScale
    rotation[[g]] <- pc$rotation[, seq_len(nPC), drop = FALSE]
  }

  structure(
    amp,
    genesUsed = genesUsed, varExplained = varExplained,
    center = center, scale = scaleVec, rotation = rotation
  )
}

# Selects the rhythmic gene panel for one already-z-scored block and returns the
# per-sample coherence on it. `z` is genes x samples (rows z-scored across the
# samples of this group), `Q` an orthonormal basis of the centred cosinor design.
# Split out from clockCoherence() only so the per-tissue loop stays readable.
.coherenceOnBlock <- function(z, Q, nPerm, fdr, minPanel, seed) {
  ssTot <- rowSums(z^2)
  usable <- ssTot > .Machine$double.eps
  z <- z[usable, , drop = FALSE]
  ssTot <- ssTot[usable]
  if (nrow(z) < minPanel) {
    return(list(coherence = rep(NA_real_, ncol(z)), panel = character(0),
                r2 = numeric(0), cutoff = NA_real_, relaxed = FALSE))
  }

  # R^2 of the cosinor fit per gene. Projecting onto the orthonormal Q rather than
  # forming the n x n hat matrix keeps the permutation loop O(G n) instead of
  # O(G n^2), which is what makes a real permutation null affordable on a
  # transcriptome-sized matrix.
  r2 <- rowSums((z %*% Q)^2) / ssTot

  # Data-driven cutoff: the panel is whatever genes beat a phase-permuted null at
  # the requested FDR, so panel size adapts to how much circadian signal this
  # particular dataset actually carries (Step 3.2 subtask 1) instead of being a
  # fixed count or a fixed variance quantile.
  oldSeed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  set.seed(seed)
  nullR2 <- unlist(lapply(seq_len(nPerm), function(b) {
    Qp <- Q[sample.int(nrow(Q)), , drop = FALSE]
    rowSums((z %*% Qp)^2) / ssTot
  }), use.names = FALSE)
  if (is.null(oldSeed)) {
    if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", oldSeed, envir = globalenv())
  }

  # Empirical p from the pooled null, then BH. Pooling the null across genes is
  # valid here because every gene shares the same design and the same z-scaling,
  # so its null R^2 distribution depends only on n - and pooling buys nPerm * G
  # null draws instead of nPerm per gene.
  nullSorted <- sort(nullR2)
  pEmp <- (length(nullSorted) - findInterval(r2, nullSorted) + 1) / (length(nullSorted) + 1)
  q <- stats::p.adjust(pEmp, method = "BH")
  panelIdx <- which(q < fdr)

  relaxed <- FALSE
  if (length(panelIdx) < minPanel) {
    # Not enough genes clear the null - fall back to the top minPanel by R^2 so QC
    # still returns a score, but mark it so the caller can report that this
    # dataset's circadian signal was too weak for the null-based selection.
    relaxed <- TRUE
    panelIdx <- utils::head(order(r2, decreasing = TRUE), minPanel)
  }

  zPanel <- z[panelIdx, , drop = FALSE]
  fitted <- (zPanel %*% Q) %*% t(Q)
  resid <- zPanel - fitted
  ssTotSample <- colSums(zPanel^2)
  coh <- 1 - colSums(resid^2) / ssTotSample
  coh[!is.finite(coh) | ssTotSample <= .Machine$double.eps] <- NA_real_

  list(coherence = coh, panel = rownames(zPanel), r2 = r2[panelIdx],
       cutoff = if (relaxed) NA_real_ else min(r2[panelIdx]), relaxed = relaxed)
}

#' Dynamic per-sample circadian coherence (Phase 3 QC)
#'
#' The Phase 3 Step 3.2 replacement for the fixed 16-gene amplitude flag
#' ([clockAmplitude()]). Rather than scoring every dataset against the same
#' core-clock panel, this selects the informative gene panel *from the input
#' dataset itself*: each gene is fitted with a cosinor against the cohort's
#' (estimated) phases, and the panel is the set of genes whose fit beats a
#' phase-permuted null at a chosen FDR. A sample's coherence is then the fraction
#' of its own deviation-from-gene-mean that the cohort's fitted rhythm explains -
#' 1 for a sample sitting exactly on the cohort's circadian manifold, near 0 (or
#' negative) for a sample whose clock genes are damped, desynchronised, or simply
#' inconsistent with the phase it was assigned.
#'
#' @section Why replace the fixed panel:
#' The legacy flag scores amplitude on 16 core clock genes. Its documented failure
#' (see [clockAmplitude()]) is that aging and similar insults leave core clock
#' gene *expression* close to unaltered while reprogramming which downstream genes
#' oscillate - so a fixed core-clock panel is measuring the wrong thing. Selecting
#' the panel per dataset lets those reprogrammed, dataset-specific oscillators
#' carry the QC signal.
#'
#' @section Circularity - read before interpreting:
#' `phase` is normally the *estimated* phase from [estimatePhase()], not ground
#' truth, so a sample whose phase was estimated badly will fit the cohort cosinor
#' badly and score low. That is the intended behaviour (the flag is meant to catch
#' exactly those samples) but it does mean coherence is **not** an independent
#' measurement of data quality - it is a self-consistency score. Validation
#' against held-out true phase (`benchmarks/run_dynamic_qc_benchmark.R`) is
#' therefore the only claim made for it here; do not read it as a general
#' sample-quality metric. Supplying known sampling times as `phase` removes the
#' circularity and turns it into a pure rhythm-strength score.
#'
#' @section Cohort-level by construction:
#' Unlike the legacy flag, this needs a cohort: the panel and the fitted rhythm
#' are both estimated across samples. Below `minSamples` (default 10) it returns
#' all `NA` rather than a number computed from too little data, and
#' [estimatePhase()] falls back to the legacy flag in that case. This is a
#' deliberate design choice over storing a fixed dynamic basis in the bundled
#' reference models - a stored panel would be "dynamic" only with respect to the
#' training set, which is the very thing Step 3.2 sets out to remove.
#'
#' @param exprMat Genes x samples numeric matrix, log-scale, with gene rownames.
#' @param phase Numeric vector of phases in hours, `ncol(exprMat)` long, usually
#'   [estimatePhase()]'s output. `NA` phases drop those samples from panel
#'   selection and receive `NA` coherence.
#' @param tissue Optional character/factor vector, one label per sample. If
#'   supplied, genes are z-scored and the panel is selected within each tissue
#'   group separately (as in [clockAmplitude()]).
#' @param period Rhythm period in hours (default 24).
#' @param nPerm Number of phase permutations forming the null (default 200).
#' @param fdr Benjamini-Hochberg FDR for panel selection (default 0.05).
#' @param minPanel Minimum panel size (default 20). If fewer genes clear the null,
#'   the top `minPanel` genes by cosinor R-squared are used instead and the
#'   `relaxed` attribute is set for that group.
#' @param minSamples Minimum samples per group (default 10); smaller groups get
#'   `NA`.
#' @param seed Seed for the permutation null, so the panel is reproducible.
#' @return A numeric vector, `ncol(exprMat)` long, of coherence scores (1 = perfect
#'   fit to the cohort rhythm, lower = worse; can be negative). Attributes:
#'   `panel` (named list, selected gene symbols per group), `panelSize` (named
#'   integer), `relaxed` (named logical, `TRUE` where the null-based selection
#'   fell through to the top-`minPanel` fallback), `r2Cutoff` (named numeric, the
#'   smallest R-squared admitted, `NA` when relaxed).
#' @seealso [clockAmplitude()] for the legacy fixed-panel flag,
#'   [deriveAmplitudeThreshold()] for turning either score into a cutoff.
#' @examples
#' set.seed(1)
#' nGenes <- 200; nSamples <- 24
#' t <- seq(0, 23, length.out = nSamples)
#' truePhase <- runif(nGenes, 0, 24)
#' mat <- t(vapply(seq_len(nGenes), function(g) {
#'   2 * cos(2 * pi * (t - truePhase[g]) / 24) + rnorm(nSamples, sd = 1)
#' }, numeric(nSamples)))
#' dimnames(mat) <- list(paste0("g", seq_len(nGenes)), paste0("s", seq_len(nSamples)))
#' coh <- clockCoherence(mat, phase = t)
#' summary(coh)
#' @export
clockCoherence <- function(exprMat, phase, tissue = NULL, period = 24,
                           nPerm = 200, fdr = 0.05, minPanel = 20,
                           minSamples = 10, seed = 1) {
  stopifnot(is.matrix(exprMat))
  if (is.null(rownames(exprMat))) {
    stop("exprMat must have rownames (gene identifiers)")
  }
  nSamples <- ncol(exprMat)
  if (length(phase) != nSamples) {
    stop("phase must be length ncol(exprMat)")
  }
  if (!is.null(tissue) && length(tissue) != nSamples) {
    stop("tissue must be NULL or length ncol(exprMat)")
  }
  if (nPerm < 1) stop("nPerm must be >= 1")
  if (fdr <= 0 || fdr >= 1) stop("fdr must be in (0, 1)")

  grp <- if (is.null(tissue)) rep("all", nSamples) else as.character(tissue)
  groups <- unique(grp)
  coh <- rep(NA_real_, nSamples)
  panel <- stats::setNames(vector("list", length(groups)), groups)
  panelSize <- stats::setNames(rep(0L, length(groups)), groups)
  relaxed <- stats::setNames(rep(NA, length(groups)), groups)
  r2Cutoff <- stats::setNames(rep(NA_real_, length(groups)), groups)

  for (g in groups) {
    idx <- which(grp == g & !is.na(phase))
    if (length(idx) < minSamples) next

    block <- exprMat[, idx, drop = FALSE]
    z <- t(scale(t(block)))
    z[is.na(z)] <- 0 # a gene with no variance in this group carries no signal

    theta <- 2 * pi * phase[idx] / period
    D <- scale(cbind(cos(theta), sin(theta)), center = TRUE, scale = FALSE)
    qrD <- qr(D)
    if (qrD$rank < 2) next # phases all identical/collinear - no rhythm to fit
    Q <- qr.Q(qrD)[, seq_len(2), drop = FALSE]

    res <- .coherenceOnBlock(z, Q, nPerm = nPerm, fdr = fdr,
                             minPanel = minPanel, seed = seed)
    coh[idx] <- res$coherence
    panel[[g]] <- res$panel
    panelSize[g] <- length(res$panel)
    relaxed[g] <- res$relaxed
    r2Cutoff[g] <- res$cutoff
  }

  structure(coh, panel = panel, panelSize = panelSize,
            relaxed = relaxed, r2Cutoff = r2Cutoff)
}

#' Bootstrap confidence interval for a phase-inference engine's predictions
#'
#' Refits the named engine on `nBoot` donor-resampled-with-replacement copies of
#' the training data and predicts `testMat` each time; the spread of those
#' bootstrap predictions gives each test sample a circular confidence interval
#' (via [circularQuantileInterval()]), while a single fit on the original,
#' non-resampled training data gives the reported point estimate. This wraps
#' whichever registered engine's `fit`/`predict` functions (the same
#' `.circaTimeEngines` registry `benchmarkPhase()` uses); it does not invent a new
#' phase-inference method (AGENTS.md 0.3).
#'
#' A bootstrap replicate's `fit`/`predict` call is allowed to fail (e.g.
#' tauFisher's MetaCycle step erroring when a resample happens to drop every
#' replicate of one clock hour - the same failure mode [benchmarkPhase()] isolates
#' per-fold). Failed replicates are dropped from that sample's interval rather
#' than aborting the whole bootstrap (AGENTS.md 0.3: report the failure, do not
#' fabricate or crash).
#'
#' @section Why `residualPool` matters:
#' Refitting on resampled training data only captures how sensitive the fit is to
#' *which training donors were used* - it says nothing about the error a genuinely
#' new test sample's own biological noise or the model's irreducible bias
#' contributes, which on real data is usually the dominant source of error.
#' Verified directly on GSE54650 (mouse 12-tissue atlas): a training-resampling-only
#' bootstrap for the `timetable` engine gave a stated 80% interval that actually
#' covered truth only ~48% of the time - badly overconfident, exactly the failure
#' mode AGENTS.md Phase 3 warns about ("uncalibrated intervals are worse than
#' none"). Supplying `residualPool` - an independent sample of real circular
#' residuals (e.g. `circularDiff(pred, truth)` from a prior cross-validated
#' benchmark run on data the current `testMat` was not scored against) adds one
#' resampled residual to every bootstrap draw, folding that test-time noise back
#' into the interval. This is a standard prediction-interval technique (residual/
#' bagging bootstrap, as used e.g. for random forest prediction intervals), not a
#' new phase-inference method (AGENTS.md 0.3).
#'
#' @param method Name of a registered engine.
#' @param trainMat,trainTime,testMat As in the engine's own `fit`/`predict`.
#' @param donors Optional character/factor vector, one donor ID per training
#'   column (`ncol(trainMat)` long); resampling draws whole donors with
#'   replacement so repeated-measures structure isn't broken. If `NULL` (default),
#'   each training column is its own donor.
#' @param nBoot Number of bootstrap replicates (default 50).
#' @param level Nominal coverage of the returned interval, in `(0, 1)` (default
#'   0.8).
#' @param seed Optional integer seed; RNG state is saved/restored (AGENTS.md
#'   Section 10).
#' @param fitArgs,predictArgs Extra named arguments passed to the engine's
#'   `fit`/`predict` respectively.
#' @param period Length of the cycle (default 24).
#' @param residualPool Optional numeric vector of independent circular residuals
#'   (each in `(-period/2, period/2]`, e.g. from [circularDiff()]) to fold into
#'   the bootstrap draws - see the section above. `NULL` (default) reproduces the
#'   training-resampling-only interval, which AGENTS.md Phase 3's own calibration
#'   check found to be overconfident on real data; supplying a residual pool from
#'   an independent benchmark run is strongly recommended for any interval that
#'   will actually be reported to a user.
#' @return A `data.frame` with one row per test sample: `sample`, `phase` (point
#'   estimate from the original unresampled fit), `phase_lo`, `phase_hi`, `level`,
#'   `n_boot_ok`, `n_boot_failed`. The raw per-replicate draws are attached as
#'   `attr(., "draws")`, an `ncol(testMat) x nBoot` matrix (`NA` for failed
#'   replicates), so [calibrationCurve()] can recompute intervals at other levels
#'   without refitting anything.
#' @examples
#' set.seed(1)
#' nGenes <- 15; nSamples <- 20
#' trainTime <- runif(nSamples, 0, 24)
#' theta <- trainTime / 24 * 2 * pi
#' truePhase <- runif(nGenes, 0, 2 * pi)
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - truePhase[i]) + rnorm(nSamples, sd = 0.3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' ci <- bootstrapPhaseCI("timetable", trainMat, trainTime, trainMat[, 1:5, drop = FALSE],
#'   nBoot = 5, seed = 1
#' )
#' ci
#' @export
bootstrapPhaseCI <- function(method, trainMat, trainTime, testMat, donors = NULL,
                              nBoot = 50, level = 0.8, seed = NULL,
                              fitArgs = list(), predictArgs = list(), period = 24,
                              residualPool = NULL) {
  engine <- .circaTimeEngines[[method]]
  if (is.null(engine)) {
    stop("Unknown method '", method, "'. Available: ",
         paste(names(.circaTimeEngines), collapse = ", "))
  }
  stopifnot(level > 0, level < 1)

  nTrain <- ncol(trainMat)
  if (is.null(donors)) donors <- seq_len(nTrain)
  if (length(donors) != nTrain) stop("donors must be length ncol(trainMat)")
  uniqueDonors <- unique(donors)

  sampleNames <- colnames(testMat)
  if (is.null(sampleNames)) sampleNames <- as.character(seq_len(ncol(testMat)))
  nTest <- ncol(testMat)

  runOnce <- function(trIdx) {
    fitObj <- do.call(engine$fit, c(
      list(trainMat = trainMat[, trIdx, drop = FALSE], trainTime = trainTime[trIdx]),
      fitArgs
    ))
    predDf <- do.call(engine$predict, c(list(fit = fitObj, testMat = testMat), predictArgs))
    predDf$phase
  }

  # Protect RNG state around the *whole* function body, not just the bootstrap
  # loop - the point-estimate fit/predict call below can itself consume randomness
  # (e.g. a stochastic engine), and that consumption must not leak into the
  # caller's subsequent draws either (AGENTS.md Section 10).
  if (!is.null(seed)) withr::local_seed(seed) else withr::local_preserve_seed()

  point <- runOnce(seq_len(nTrain))

  draws <- matrix(NA_real_, nrow = nTest, ncol = nBoot)
  for (b in seq_len(nBoot)) {
    bootDonors <- sample(uniqueDonors, length(uniqueDonors), replace = TRUE)
    trIdx <- unlist(lapply(bootDonors, function(d) which(donors == d)), use.names = FALSE)
    rep <- tryCatch(
      runOnce(trIdx),
      error = function(e) {
        message(sprintf(
          "bootstrapPhaseCI: method '%s' replicate %d/%d failed: %s",
          method, b, nBoot, conditionMessage(e)
        ))
        rep(NA_real_, nTest)
      }
    )
    if (!is.null(residualPool)) {
      rep <- wrapPhase(rep + sample(residualPool, nTest, replace = TRUE), period = period)
    }
    draws[, b] <- rep
  }

  intervals <- lapply(seq_len(nTest), function(i) {
    circularQuantileInterval(draws[i, ], level = level, period = period)
  })
  loVec <- vapply(intervals, `[[`, numeric(1), "lo")
  hiVec <- vapply(intervals, `[[`, numeric(1), "hi")
  nOk <- rowSums(!is.na(draws))

  out <- data.frame(
    sample = sampleNames,
    phase = point,
    phase_lo = loVec,
    phase_hi = hiVec,
    level = level,
    n_boot_ok = nOk,
    n_boot_failed = nBoot - nOk,
    row.names = NULL
  )
  attr(out, "draws") <- draws
  out
}

#' Empirical coverage of bootstrap intervals across nominal levels
#'
#' Recomputes [circularQuantileInterval()] at each of `levels` from the raw
#' bootstrap draws already computed by [bootstrapPhaseCI()] (its `"draws"`
#' attribute), and reports what fraction of the time the interval actually
#' contains the known ground truth. This is the calibration check AGENTS.md
#' Phase 3 task 3 requires before a stated interval can be trusted at all:
#' "uncalibrated intervals are worse than none."
#'
#' @param ci A `data.frame` returned by [bootstrapPhaseCI()] (or `rbind()` of
#'   several such calls made with the same `nBoot`), carrying a `"draws"`
#'   attribute.
#' @param truth Numeric vector of ground-truth phases, same length and order as
#'   `nrow(ci)`.
#' @param levels Numeric vector of nominal coverage levels to evaluate (default
#'   `c(0.5, 0.7, 0.8, 0.9, 0.95)`).
#' @param period Length of the cycle (default 24).
#' @return A `data.frame` with one row per level: `nominal`, `empirical`
#'   (observed coverage fraction), `n` (non-`NA` intervals evaluated).
#' @examples
#' set.seed(1)
#' nGenes <- 15; nSamples <- 20
#' trainTime <- runif(nSamples, 0, 24)
#' theta <- trainTime / 24 * 2 * pi
#' truePhase <- runif(nGenes, 0, 2 * pi)
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   5 + 2 * cos(theta - truePhase[i]) + rnorm(nSamples, sd = 0.3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' ci <- bootstrapPhaseCI("timetable", trainMat, trainTime, trainMat[, 1:5, drop = FALSE],
#'   nBoot = 5, seed = 1
#' )
#' calibrationCurve(ci, truth = trainTime[1:5], levels = c(0.5, 0.8))
#' @export
calibrationCurve <- function(ci, truth, levels = c(0.5, 0.7, 0.8, 0.9, 0.95), period = 24) {
  draws <- attr(ci, "draws")
  if (is.null(draws)) {
    stop("ci must carry a 'draws' attribute, as produced by bootstrapPhaseCI()")
  }
  if (nrow(draws) != length(truth)) {
    stop("truth must be length nrow(ci)")
  }

  rows <- lapply(levels, function(lv) {
    intervals <- lapply(seq_len(nrow(draws)), function(i) {
      circularQuantileInterval(draws[i, ], level = lv, period = period)
    })
    lo <- vapply(intervals, `[[`, numeric(1), "lo")
    hi <- vapply(intervals, `[[`, numeric(1), "hi")
    covers <- circularIntervalCovers(truth, lo, hi, period = period)
    data.frame(nominal = lv, empirical = mean(covers, na.rm = TRUE), n = sum(!is.na(covers)))
  })
  do.call(rbind, rows)
}

#' Plot a calibration curve
#'
#' Base-graphics plot of empirical vs. nominal coverage (AGENTS.md Phase 3 task 3:
#' "Plot it"), with the identity line marking perfect calibration. Points below
#' the line mean the stated interval is too narrow (overconfident); points above
#' mean it is too wide (underconfident, but honest).
#'
#' @param calCurve Output of [calibrationCurve()].
#' @return `calCurve`, invisibly. Called for its plotting side effect.
#' @examples
#' calCurve <- data.frame(nominal = c(0.5, 0.8, 0.95), empirical = c(0.53, 0.81, 0.96))
#' plotCalibration(calCurve)
#' @export
plotCalibration <- function(calCurve) {
  graphics::plot(
    calCurve$nominal, calCurve$empirical,
    xlim = c(0, 1), ylim = c(0, 1), pch = 19,
    xlab = "Nominal coverage", ylab = "Empirical coverage",
    main = "Calibration"
  )
  graphics::abline(0, 1, lty = 2, col = "grey50")
  invisible(calCurve)
}

#' Empirically derive a low-amplitude flag threshold
#'
#' Finds the [clockAmplitude()] cutoff below which phase error is meaningfully
#' worse, from held-out labelled data - AGENTS.md 7.4/Phase 3 task 4 requires this
#' be derived empirically, not from intuition.
#'
#' Candidate thresholds are deciles of `amp` (10th-50th percentile). For each
#' candidate, samples split into `low` (`amp <= candidate`) and `high`
#' (`amp > candidate`) groups; the candidate is scored by a one-sided
#' Wilcoxon rank-sum test of `error(low) > error(high)`, converted to a
#' common-language effect size (probability a random low-group error exceeds a
#' random high-group error, in `[0, 1]`), subject to the `low` group holding at
#' least `minLowFrac` of samples - too small a `low` group is not a reliable flag,
#' however well it happens to separate this particular dataset. The candidate
#' with the highest effect size among those satisfying the fraction constraint is
#' returned.
#'
#' @param amp Numeric vector of [clockAmplitude()] values.
#' @param error Numeric vector of per-sample circular error (e.g. from
#'   [circularError()] against known ground truth), same length as `amp`.
#' @param minLowFrac Minimum fraction of samples the `low` group must contain
#'   (default 0.05).
#' @return A list: `threshold` (chosen amplitude cutoff), `effect_size`
#'   (common-language effect size at that cutoff), `p_value` (Wilcoxon test
#'   p-value), `low_mae`/`high_mae` (mean error in each group), `low_frac`
#'   (fraction of samples flagged `low_amplitude`).
#' @examples
#' set.seed(1)
#' amp <- c(runif(15, 0, 1), runif(15, 2, 4))
#' error <- c(runif(15, 3, 5), runif(15, 0, 1))
#' deriveAmplitudeThreshold(amp, error)
#' @export
deriveAmplitudeThreshold <- function(amp, error, minLowFrac = 0.05) {
  keep <- !is.na(amp) & !is.na(error)
  amp <- amp[keep]; error <- error[keep]
  if (length(amp) < 20) {
    stop("Need at least 20 non-NA (amp, err) pairs to derive a threshold reliably")
  }

  candidates <- unique(stats::quantile(amp, probs = seq(0.1, 0.5, by = 0.05), names = FALSE))

  best <- NULL
  for (cand in candidates) {
    lowIdx <- amp <= cand
    lowFrac <- mean(lowIdx)
    if (lowFrac < minLowFrac || lowFrac > 0.9) next
    if (sum(lowIdx) < 5 || sum(!lowIdx) < 5) next

    test <- stats::wilcox.test(error[lowIdx], error[!lowIdx], alternative = "greater")
    effectSize <- unname(test$statistic) / (sum(lowIdx) * sum(!lowIdx))
    if (is.null(best) || effectSize > best$effect_size) {
      best <- list(
        threshold = cand, effect_size = effectSize, p_value = test$p.value,
        low_mae = mean(error[lowIdx]), high_mae = mean(error[!lowIdx]),
        low_frac = lowFrac
      )
    }
  }

  if (is.null(best)) {
    stop("No candidate threshold satisfied minLowFrac constraints; check amp/err inputs")
  }
  best
}

#' Bayesian phase posterior built on the residual-pool bootstrap
#'
#' Phase 3 Step 3.3. This does not replace [bootstrapPhaseCI()]/
#' [circularQuantileInterval()]'s empirical interval - it derives a genuinely
#' different, complementary object from the *same* residual pool: a von Mises
#' posterior predictive density for a sample's true phase, centred on the point
#' estimate, with concentration ([circularConcentration()]) measured directly
#' from that pool rather than assumed. [circularQuantileInterval()] answers "what
#' interval did the observed residuals fall in"; this answers "what is the
#' probability the true phase is any given value" - useful whenever a caller
#' needs more than an interval endpoint, e.g. `P(|error| < 2h)` for a specific
#' downstream decision threshold, or a smooth density to plot.
#'
#' @section Relationship to the residual pool's own coverage:
#' Because the credible interval here is the equal-tailed interval of a
#' *parametric* (unimodal, symmetric) fit to the pool's concentration, while
#' [circularQuantileInterval()]'s interval is the *empirical* quantile interval
#' of the pool itself, the two intervals agree exactly only if the pool really is
#' von Mises-distributed. Real residual pools are approximately but not exactly
#' von Mises (see `benchmarks/run_bayesian_ci_benchmark.R` for a direct
#' comparison on real data) - the intended use is as a complementary posterior
#' summary, not a replacement calibration target for the already-validated
#' empirical interval.
#'
#' @param point Single numeric point phase estimate (or a vector, one per
#'   sample), in `[0, period)` - typically `estimatePhase()`'s own point estimate.
#' @param residualPool Numeric vector of real circular residuals (signed
#'   deviations from truth, as stored in `PhaseReferenceModel@residualPool`).
#'   Concentration is measured from this pool once and shared across every
#'   `point` (matching how the reference model's residual pool is
#'   tissue/model-wide, not per-sample).
#' @param level Nominal credible level for the returned interval (default 0.8).
#' @param period Length of the cycle (default 24).
#' @param nGrid Number of points in the returned density grid per sample
#'   (default 121, i.e. every 12 minutes over a 24h period).
#' @return A list with `kappa` (fitted von Mises concentration, shared across all
#'   `point`s), `concentration` (the underlying mean resultant length, `[0, 1]`),
#'   `credible_lo`/`credible_hi` (numeric vectors, one per `point`), and `density`
#'   (a list of `data.frame(phase, density)`, one per `point`, each integrating to
#'   1 over `[0, period)`).
#' @examples
#' pool <- circularDiff(rnorm(200, sd = 1.5), 0) # a tight, roughly von Mises pool
#' post <- bayesianPhasePosterior(point = c(6, 18), residualPool = pool)
#' post$kappa
#' post$credible_lo
#' @export
bayesianPhasePosterior <- function(point, residualPool, level = 0.8, period = 24,
                                   nGrid = 121) {
  residualPool <- residualPool[!is.na(residualPool)]
  if (length(residualPool) < 2) {
    stop("residualPool must have at least 2 non-NA values")
  }
  R <- circularConcentration(residualPool, period = period)
  kappa <- .vonMisesKappaFromR(R)

  ci <- circularCredibleInterval(point, kappa, level = level, period = period)

  grid <- seq(0, period, length.out = nGrid)[-nGrid] # exclude the duplicate endpoint
  densityAt <- function(mu) {
    theta <- (grid - mu) / period * 2 * pi
    dens <- exp(kappa * cos(theta))
    if (kappa == 0) {
      dens <- rep(1, length(grid))
    }
    dens <- dens / (sum(dens) * (period / length(grid))) # normalise to integrate to 1
    data.frame(phase = grid, density = dens)
  }

  list(
    kappa = kappa, concentration = R,
    credible_lo = ci$lo, credible_hi = ci$hi,
    density = lapply(point, densityAt)
  )
}
