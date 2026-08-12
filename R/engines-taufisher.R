# tauFisher engine (Duan & Ngo et al. 2024, Nat Commun: 10.1038/s41467-024-48041-6).
# Wraps the authors' own R package (github.com/micnngo/tauFisher, GPL-3.0) rather than
# reimplementing it, per AGENTS.md Section 4/9 Phase 2. Method: find within-tissue
# rhythmic genes via MetaCycle (JTK_CYCLE / Lomb-Scargle), smooth each gene's training
# curve onto an hourly grid with functional data analysis, build a scaled matrix of
# pairwise gene-expression differences, reduce with PCA, and classify a new sample's
# hour via multinomial logistic regression on the first two PCs.

#' Fit a tauFisher model from labelled training data
#'
#' Calls the authors' own pipeline (`tauFisher::run_metacycle()` +
#' `find_periodic_genes()` + `train_tauFisher()`) rather than reimplementing it.
#' `run_metacycle()` writes intermediate files; these are written to and cleaned up
#' from a temporary directory, never the working directory.
#'
#' @param trainMat Numeric matrix of expression, genes x samples (linear or
#'   log-scale; `train_tauFisher()` log2-transforms internally).
#' @param trainTime Numeric vector of known circadian phase in `[0, 24)`, one value
#'   per column of `trainMat`.
#' @param method Character vector of MetaCycle method(s) to detect rhythmic genes,
#'   passed to `run_metacycle()`/`find_periodic_genes()`. Default `c("JTK", "LS")`,
#'   the package's own default.
#' @param thres p-value threshold passed to `find_periodic_genes()`. Default 1.01
#'   (i.e. no filtering), matching the authors' own `SingleDataset` vignette, which
#'   ranks candidate genes by p-value and takes the top 10 rather than thresholding.
#' @param nrep Number of training replicates, passed to `train_tauFisher()`. Default
#'   1, appropriate whenever there is exactly one array per (tissue, timepoint) cell
#'   (as for the GSE54650 benchmark backbone); set higher only if `trainMat` truly
#'   contains stacked biological replicates as separate columns per timepoint.
#' @param numbasis Number of Fourier basis functions for the FDA smoothing step,
#'   passed to `train_tauFisher()`. Default 5, the package's own default.
#' @param seed Optional integer seed. `train_tauFisher()` fits a `nnet::multinom`
#'   model with random initial weights; if `seed` is supplied, the global RNG state
#'   is saved and restored (AGENTS.md Section 10), so this call does not affect the
#'   caller's subsequent random draws.
#' @return An object of class `"tauFisherFit"`: a list with the chosen periodic genes
#'   (`Genes`), the PCA embedding (`PCA_embedding`) and trained multinomial model
#'   (`Model`) returned by `train_tauFisher()`, which MetaCycle method supplied the
#'   genes (`method`), and `period` (always 24).
#' @examples
#' \donttest{
#' set.seed(1)
#' nGenes <- 40; time <- seq(0, 23, by = 1)
#' theta <- time / 24 * 2 * pi
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   100 + 40 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(length(time), sd = 3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' fit <- fitTauFisher(trainMat, time, seed = 1)
#' length(fit$Genes)
#' }
#' @export
fitTauFisher <- function(trainMat, trainTime, method = c("JTK", "LS"), thres = 1.01,
                          nrep = 1, numbasis = 5, seed = NULL) {
  .requireEngine("tauFisher", 'remotes::install_github("micnngo/tauFisher")')
  stopifnot(is.matrix(trainMat), ncol(trainMat) == length(trainTime))
  if (is.null(rownames(trainMat))) {
    stop("trainMat must have rownames (gene identifiers)")
  }

  workDir <- tempfile("taufisher_metacycle_")
  dir.create(workDir)
  on.exit(unlink(workDir, recursive = TRUE), add = TRUE)
  metaFile <- file.path(workDir, "expr_for_meta2d.txt")

  utils::capture.output(
    tauFisher::run_metacycle(
      df = trainMat, timepoints = trainTime,
      out_file = metaFile, out_dir = workDir, method = method
    )
  )
  meta2dFile <- file.path(workDir, paste0("meta2d_", basename(metaFile)))
  genesByMethod <- tauFisher::find_periodic_genes(
    input_file = meta2dFile, test_genes = rownames(trainMat),
    per_method = method, thres = thres
  )

  # The authors' own vignette picks whichever method yields usable genes (in their
  # worked example, Lomb-Scargle found none and JTK_CYCLE was used); we follow the
  # same rule, trying `method` in the order given.
  periodicGenes <- NULL
  usedMethod <- NA_character_
  for (m in method) {
    cand <- genesByMethod[[m]]
    if (!is.null(cand) && length(cand) >= 2) {
      periodicGenes <- cand
      usedMethod <- m
      break
    }
  }
  if (is.null(periodicGenes)) {
    stop("tauFisher found no usable periodic genes via ", paste(method, collapse = "/"),
         " at thres = ", thres)
  }

  if (!is.null(seed)) withr::local_seed(seed) else withr::local_preserve_seed()

  utils::capture.output(
    trained <- tauFisher::train_tauFisher(
      periodic_genes = periodicGenes, train_df = trainMat, train_time = trainTime,
      nrep = nrep, numbasis = numbasis
    )
  )

  structure(
    list(
      Genes = trained$Genes,
      PCA_embedding = trained$PCA_embedding,
      Model = trained$Model,
      method = usedMethod,
      period = 24
    ),
    class = "tauFisherFit"
  )
}

#' Predict circadian phase with a fitted tauFisher model
#'
#' Calls the authors' own `tauFisher::test_tauFisher()`. Unlike `train_tauFisher()`,
#' this step is deterministic given the fitted model (no RNG use).
#'
#' @param fit A `"tauFisherFit"` object from [fitTauFisher()].
#' @param testMat Numeric matrix of expression, genes x samples, same scale as the
#'   training matrix. Prediction requires every one of `fit$Genes` to be present in
#'   `testMat`'s rownames (tauFisher's gene-pair-difference features need the full
#'   set); samples are only ever predicted as a whole model, so if any fitted gene is
#'   missing, every sample gets `NA` rather than a prediction computed on a silently
#'   shrunken gene set.
#' @param minGenes Minimum number of `fit$Genes` that must be present in `testMat`
#'   to attempt prediction at all (default 2, the minimum tauFisher needs to form
#'   even one gene pair).
#' @return A `data.frame` with one row per column of `testMat`: `sample`, `phase`
#'   (integer hour in `[0, 24)`, from `nnet::multinom`'s discrete class prediction,
#'   or `NA` if prediction was not attempted), and `confidence` (`NA_real_` for
#'   every row - tauFisher has no native per-sample confidence signal; calibrated
#'   uncertainty is Phase 3 scope, AGENTS.md Section 9).
#' @examples
#' \donttest{
#' set.seed(1)
#' nGenes <- 40; time <- seq(0, 23, by = 1)
#' theta <- time / 24 * 2 * pi
#' trainMat <- t(sapply(seq_len(nGenes), function(i) {
#'   100 + 40 * cos(theta - runif(1, 0, 2 * pi)) + rnorm(length(time), sd = 3)
#' }))
#' rownames(trainMat) <- paste0("gene", seq_len(nGenes))
#' fit <- fitTauFisher(trainMat, time, seed = 1)
#' pred <- predictTauFisher(fit, trainMat)
#' head(pred)
#' }
#' @export
predictTauFisher <- function(fit, testMat, minGenes = 2) {
  .requireEngine("tauFisher", 'remotes::install_github("micnngo/tauFisher")')
  # tauFisher's own NAMESPACE does not import nnet, so nnet's namespace - and with
  # it the S3 method predict.multinom the fitted $Model dispatches to - is only
  # ever loaded as a side effect of tauFisher::train_tauFisher() having called
  # nnet::multinom() during TRAINING. A prediction-only session that never trained
  # anything itself (the normal case: predict against a bundled reference model)
  # never triggers that load, so predict() on the trained multinom object fails
  # with "no applicable method for 'predict'" - not a tauFisher/testMat problem at
  # all. nnet ships with every R installation (a base "Recommended" package), so
  # this is always available; it just isn't loaded yet.
  .requireEngine("nnet", "install.packages(\"nnet\")")
  stopifnot(inherits(fit, "tauFisherFit"))
  if (is.null(rownames(testMat))) {
    stop("testMat must have rownames (gene identifiers)")
  }

  sampleNames <- colnames(testMat)
  nSamples <- ncol(testMat)
  if (is.null(sampleNames)) sampleNames <- seq_len(nSamples)

  # fit$Genes is always uppercase (train_tauFisher() uppercases internally); mirror
  # test_tauFisher()'s own toupper(rownames(test_df)) before checking presence, so
  # this check tests exactly what test_tauFisher() will actually match on.
  present <- fit$Genes %in% toupper(rownames(testMat))
  outPhase <- rep(NA_real_, nSamples)

  if (sum(present) >= minGenes && sum(present) == length(fit$Genes)) {
    result <- tauFisher::test_tauFisher(
      chosen_genes = fit$Genes, test_df = testMat,
      test_time = rep(NA_real_, nSamples),
      pca_embedding = fit$PCA_embedding, trained_model = fit$Model
    )
    outPhase <- as.numeric(result$predicted)
  }

  data.frame(
    sample = sampleNames,
    phase = outPhase,
    confidence = NA_real_,
    row.names = NULL
  )
}
