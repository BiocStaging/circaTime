# Tests for R/estimatePhase.R (AGENTS.md Section 6 main entry point) and
# R/AllClasses.R's PhaseReferenceModel. The happy-path tests use the REAL bundled
# .referenceModels (R/sysdata.rda, built by inst/scripts/train_reference_models.R
# from GSE54650) - synthetic testMat rows are constructed from that real fit's own
# template coefficients so the test is deterministic without depending on
# benchmarks/data/GSE54650_se.rds (gitignored, not guaranteed present in CI). Flag-
# specific edge cases use a temporarily injected fake PhaseReferenceModel for full
# control, restored via on.exit - the same pattern test-confidence.R's residualPool
# regression test uses for .circaTimeEngines.

test_that("PhaseReferenceModel can be constructed and shown", {
  rm <- methods::new("PhaseReferenceModel",
    species = "Mmusculus", tissue = "Liv", method = "timetable",
    fit = list(genes = c("Arntl", "Per1")), residualPool = rnorm(10, sd = 1),
    ampCenter = c(Arntl = 0, Per1 = 0), ampScale = c(Arntl = 1, Per1 = 1),
    ampRotation = matrix(c(0.7, 0.7), nrow = 2), ampGenes = c("Arntl", "Per1"),
    ampThreshold = 2, trainedOn = "unit test fixture"
  )
  out <- utils::capture.output(rm)
  expect_true(any(grepl("Mmusculus", out)))
  expect_true(any(grepl("Liv", out)))
})

test_that("availableReferenceModels lists real bundled models", {
  avail <- availableReferenceModels()
  expect_true(all(c("species", "tissue", "method", "platform", "n_amp_genes", "n_residual_draws") %in% colnames(avail)))
  expect_gt(nrow(avail), 0)
  expect_true(all(avail$species == "Mmusculus"))
  expect_true("timetable" %in% avail$method)
  expect_true("microarray" %in% avail$platform)
})

test_that("estimatePhase requires tissue with no default", {
  mat <- matrix(1, nrow = 5, dimnames = list(paste0("g", 1:5), "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  expect_error(estimatePhase(se, method = "timetable"), "tissue")
})

test_that("estimatePhase errors clearly on an unsupported combination, listing available ones", {
  mat <- matrix(1, nrow = 5, dimnames = list(paste0("g", 1:5), "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  expect_error(
    estimatePhase(se, method = "timetable", tissue = "NotARealTissue"),
    "No bundled reference model"
  )
})

# Builds a synthetic expression matrix that round-trips through a REAL bundled
# timetable model's own harmonic template at a chosen target phase, plus that
# model's own ampGenes at a mild offset from the training center - deterministic,
# real-model-shaped input without depending on the original GSE54650 file.
makeRoundTripMat <- function(refModel, targetPhase, nSamples = 1) {
  fit <- refModel@fit
  theta <- targetPhase / fit$period * 2 * pi
  fitVals <- fit$mesor + fit$b * cos(theta) + fit$c * sin(theta)
  fitMat <- matrix(rep(fitVals, nSamples), nrow = length(fitVals))
  rownames(fitMat) <- fit$genes

  ampVals <- refModel@ampCenter + refModel@ampScale # one SD above this tissue's mean
  ampMat <- matrix(rep(ampVals, nSamples), nrow = length(ampVals))
  rownames(ampMat) <- names(refModel@ampCenter)

  mat <- rbind(fitMat, ampMat)
  colnames(mat) <- paste0("s", seq_len(nSamples))
  mat
}

test_that("estimatePhase happy path: real bundled timetable model recovers its own template phase", {
  avail <- availableReferenceModels()
  tis <- avail$tissue[avail$method == "timetable"][1]
  refModel <- circaTime:::.referenceModels[[paste("Mmusculus", tis, "timetable", "microarray", sep = ".")]]

  mat <- makeRoundTripMat(refModel, targetPhase = 6, nSamples = 3)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  se <- estimatePhase(se, method = "timetable", tissue = tis)
  cd <- SummarizedExperiment::colData(se)

  expect_true(all(!is.na(cd$phase)))
  expect_true(all(circularError(cd$phase, rep(6, 3)) < 1)) # exact template match, should be near-perfect
  expect_true(all(!is.na(cd$phase_lo) & !is.na(cd$phase_hi)))
  expect_true(all(!is.na(cd$clock_amp))) # proves reference-projected amplitude works, not just NA
  expect_s3_class(cd$phase_flag, "factor")
})

test_that("estimatePhase computes a real clock_amp for a single-sample (n=1) call", {
  avail <- availableReferenceModels()
  tis <- avail$tissue[avail$method == "timetable"][1]
  refModel <- circaTime:::.referenceModels[[paste("Mmusculus", tis, "timetable", "microarray", sep = ".")]]

  mat <- makeRoundTripMat(refModel, targetPhase = 12, nSamples = 1)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  se <- estimatePhase(se, method = "timetable", tissue = tis)
  cd <- SummarizedExperiment::colData(se)

  expect_equal(nrow(cd), 1)
  expect_false(is.na(cd$clock_amp))
})

test_that("estimatePhase attaches phase_credible_lo/hi close to phase_lo/hi and records kappa in metadata (Phase 3 Step 3.3)", {
  avail <- availableReferenceModels()
  tis <- avail$tissue[avail$method == "timetable"][1]
  refModel <- circaTime:::.referenceModels[[paste("Mmusculus", tis, "timetable", "microarray", sep = ".")]]

  mat <- makeRoundTripMat(refModel, targetPhase = 6, nSamples = 3)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  se <- estimatePhase(se, method = "timetable", tissue = tis)
  cd <- SummarizedExperiment::colData(se)

  expect_true(all(!is.na(cd$phase_credible_lo) & !is.na(cd$phase_credible_hi)))
  widthBayes <- circularDiff(cd$phase_credible_hi, cd$phase_credible_lo)
  widthEmp <- circularDiff(cd$phase_hi, cd$phase_lo)
  expect_true(all(abs(widthBayes - widthEmp) < 2))

  bayesMeta <- S4Vectors::metadata(se)$circaTime_bayesian
  expect_gt(bayesMeta$kappa, 0)
  expect_equal(bayesMeta$level, 0.8)
})

test_that("estimatePhase defaults to the logExpr assay, not whichever assay is first (regression test)", {
  # Real bug: an SE built the way inst/scripts/ingest_GSE54650.R builds one has
  # assays = list(expr = <linear>, logExpr = <log2>) in that order.
  # SummarizedExperiment::assay(se) with no name picks the FIRST assay, which
  # silently fed every bundled (log-scale-trained) reference model raw linear
  # values and produced materially worse predictions with no error at all.
  avail <- availableReferenceModels()
  tis <- avail$tissue[avail$method == "timetable"][1]
  refModel <- circaTime:::.referenceModels[[paste("Mmusculus", tis, "timetable", "microarray", sep = ".")]]

  logMat <- makeRoundTripMat(refModel, targetPhase = 6, nSamples = 3)
  wrongMat <- 2^logMat # plausible-looking "linear scale" nonsense relative to the log-fitted template
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = wrongMat, logExpr = logMat))

  se <- estimatePhase(se, method = "timetable", tissue = tis) # assay_name left NULL
  cd <- SummarizedExperiment::colData(se)
  expect_true(all(circularError(cd$phase, rep(6, 3)) < 1)) # matches the logExpr-only happy path test
})

# Injected fake model gives full control over gene overlap and amplitude, isolated
# from real bundled data and restored on exit.
makeFakeRefModel <- function(ampThreshold = 2) {
  # Genes need *differing* b/c coefficients (real time-indicating genes would)
  # so both observed and template vectors have nonzero variance across genes -
  # a constant template/observation pair gives an undefined (NaN) correlation in
  # predictMolecularTimetable()'s grid search.
  fakeFit <- structure(
    list(genes = paste0("G", 1:5), mesor = rep(5, 5),
         b = c(1, 0, -1, 0, 0.7), c = c(0, 1, 0, -1, 0.7),
         r = rep(0.9, 5), period = 24),
    class = "molecularTimetableFit"
  )
  methods::new("PhaseReferenceModel",
    species = "Test", tissue = "T1", method = "timetable",
    fit = fakeFit, residualPool = stats::rnorm(50, sd = 1),
    ampCenter = stats::setNames(rep(0, 4), paste0("A", 1:4)),
    ampScale = stats::setNames(rep(1, 4), paste0("A", 1:4)),
    ampRotation = matrix(rep(0.5, 4), nrow = 4),
    ampGenes = paste0("A", 1:4), ampThreshold = ampThreshold,
    trainedOn = "unit test fixture"
  )
}

withFakeRefModel <- function(refModel, code) {
  oldModels <- circaTime:::.referenceModels
  on.exit(utils::assignInNamespace(".referenceModels", oldModels, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(
    ".referenceModels",
    c(oldModels, stats::setNames(list(refModel), "Test.T1.timetable.microarray")),
    ns = "circaTime"
  )
  force(code)
}

test_that("estimatePhase flags low_amplitude vs ok using a controlled injected model", {
  fake <- makeFakeRefModel(ampThreshold = 2)
  withFakeRefModel(fake, {
    fakeFit <- fake@fit
    gVals <- fakeFit$mesor + fakeFit$b * cos(0) + fakeFit$c * sin(0) # theta = 0 template
    gMat <- matrix(rep(gVals, 2), nrow = 5, dimnames = list(fakeFit$genes, c("lo", "hi")))
    ampRows <- matrix(c(0, 0, 0, 0, 10, 10, 10, 10), nrow = 4, dimnames = list(paste0("A", 1:4), c("lo", "hi")))
    mat <- rbind(gMat, ampRows)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
    se <- estimatePhase(se, method = "timetable", tissue = "T1", species = "Test")
    cd <- SummarizedExperiment::colData(se)

    expect_false(any(is.na(cd$phase))) # exact template match at theta=0 for both samples
    expect_equal(as.character(cd$phase_flag[colnames(se) == "lo"]), "low_amplitude")
    expect_equal(as.character(cd$phase_flag[colnames(se) == "hi"]), "ok")
  })
})

test_that("estimatePhase flags gene_coverage (and never errors) when almost no genes overlap", {
  fake <- makeFakeRefModel()
  withFakeRefModel(fake, {
    mat <- matrix(1, nrow = 3, ncol = 2, dimnames = list(paste0("X", 1:3), c("s1", "s2")))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
    se <- estimatePhase(se, method = "timetable", tissue = "T1", species = "Test")
    cd <- SummarizedExperiment::colData(se)

    expect_true(all(is.na(cd$phase)))
    expect_true(all(cd$phase_flag == "gene_coverage"))
  })
})

# Builds a multi-sample SE against the fake T1 model: the 5 fit genes and 4 amp
# genes exact-template-match a per-sample target phase (so the point estimate is
# always near-perfect, isolating what qc_method actually changes), plus nExtra
# extra genes carrying real rhythmic signal at those same target phases - the raw
# material clockCoherence() needs to select its own panel from. dampIdx samples
# get that extra-gene signal replaced with noise, i.e. a damped clock that the
# fixed 9-gene fit/amp panel above cannot see (mirrors the real GSE84511 finding
# motivating Step 3.2: legacy panel unchanged, downstream genes desynchronised).
makeDynamicQcSE <- function(fake, targetPhase, nExtra = 60, dampIdx = integer(0), seed = 1) {
  set.seed(seed)
  fakeFit <- fake@fit
  theta <- targetPhase / fakeFit$period * 2 * pi
  gMat <- vapply(theta, function(th) fakeFit$mesor + fakeFit$b * cos(th) + fakeFit$c * sin(th), numeric(5))
  rownames(gMat) <- fakeFit$genes

  ampVals <- fake@ampCenter + fake@ampScale
  ampMat <- matrix(rep(ampVals, length(targetPhase)), nrow = 4,
                   dimnames = list(names(fake@ampCenter), NULL))

  extraPhase <- stats::runif(nExtra, 0, 24)
  extraTheta <- outer(targetPhase, extraPhase, function(tp, ep) 2 * pi * (tp - ep) / 24)
  extraMat <- t(2 * cos(extraTheta)) + stats::rnorm(nExtra * length(targetPhase), sd = 0.3)
  if (length(dampIdx) > 0) {
    extraMat[, dampIdx] <- stats::rnorm(nExtra * length(dampIdx), sd = 0.3)
  }
  rownames(extraMat) <- paste0("Extra", seq_len(nExtra))

  mat <- rbind(gMat, ampMat, extraMat)
  colnames(mat) <- paste0("s", seq_along(targetPhase))
  SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
}

test_that("qc_method defaults to legacy (backward compatible)", {
  fake <- makeFakeRefModel(ampThreshold = 2)
  withFakeRefModel(fake, {
    se <- makeDynamicQcSE(fake, targetPhase = rep(6, 12), seed = 1)
    se <- estimatePhase(se, method = "timetable", tissue = "T1", species = "Test")
    expect_equal(S4Vectors::metadata(se)$circaTime_qc$method_used, "legacy")
    expect_equal(S4Vectors::metadata(se)$circaTime_qc$threshold, fake@ampThreshold)
  })
})

test_that("qc_method = 'dynamic' scores damped samples lower and flags them", {
  fake <- makeFakeRefModel(ampThreshold = 2)
  withFakeRefModel(fake, {
    targetPhase <- rep(seq(0, 22, by = 2), 2) # 24 samples, real phase spread
    dampIdx <- 1:4
    se <- makeDynamicQcSE(fake, targetPhase, nExtra = 80, dampIdx = dampIdx, seed = 7)
    se <- estimatePhase(se, method = "timetable", tissue = "T1", species = "Test",
                        qc_method = "dynamic", qc_threshold = 0.3)
    cd <- SummarizedExperiment::colData(se)
    qc <- S4Vectors::metadata(se)$circaTime_qc

    expect_equal(qc$method_used, "dynamic")
    expect_gt(mean(cd$clock_amp[-dampIdx], na.rm = TRUE), mean(cd$clock_amp[dampIdx], na.rm = TRUE))
    expect_true(all(as.character(cd$phase_flag[dampIdx]) == "low_amplitude"))
    expect_gt(qc$panelSize[["all"]], 0)
  })
})

test_that("qc_method = 'dynamic' falls back to legacy below qc_min_samples", {
  fake <- makeFakeRefModel(ampThreshold = 2)
  withFakeRefModel(fake, {
    se <- makeDynamicQcSE(fake, targetPhase = c(2, 8, 14), nExtra = 20, seed = 3)
    seDynamic <- estimatePhase(se, method = "timetable", tissue = "T1", species = "Test",
                               qc_method = "dynamic", qc_min_samples = 10)
    seLegacy <- estimatePhase(se, method = "timetable", tissue = "T1", species = "Test",
                              qc_method = "legacy")
    qc <- S4Vectors::metadata(seDynamic)$circaTime_qc

    expect_equal(qc$method_requested, "dynamic")
    expect_equal(qc$method_used, "legacy")
    expect_equal(qc$threshold, fake@ampThreshold)
    expect_equal(SummarizedExperiment::colData(seDynamic)$clock_amp,
                SummarizedExperiment::colData(seLegacy)$clock_amp)
  })
})

test_that("estimatePhase rescues common clock-gene aliases and reports mapping", {
  fakeFit <- structure(
    list(genes = c("Arntl", "Per1", "Per2", "Cry1", "Cry2"), mesor = rep(5, 5),
         b = c(1, 0, -1, 0, 0.7), c = c(0, 1, 0, -1, 0.7),
         r = rep(0.9, 5), period = 24),
    class = "molecularTimetableFit"
  )
  fake <- methods::new("PhaseReferenceModel",
    species = "Mmusculus", tissue = "AliasTissue", method = "timetable",
    fit = fakeFit, residualPool = stats::rnorm(50, sd = 1),
    ampCenter = stats::setNames(rep(0, 4), c("Arntl", "Per1", "Per2", "Cry1")),
    ampScale = stats::setNames(rep(1, 4), c("Arntl", "Per1", "Per2", "Cry1")),
    ampRotation = matrix(rep(0.5, 4), nrow = 4),
    ampGenes = c("Arntl", "Per1", "Per2", "Cry1"), ampThreshold = 0,
    trainedOn = "unit test fixture"
  )

  oldModels <- circaTime:::.referenceModels
  on.exit(utils::assignInNamespace(".referenceModels", oldModels, ns = "circaTime"), add = TRUE)
  utils::assignInNamespace(
    ".referenceModels",
    c(oldModels, stats::setNames(list(fake), "Mmusculus.AliasTissue.timetable.microarray")),
    ns = "circaTime"
  )

  gVals <- fakeFit$mesor + fakeFit$b * cos(0) + fakeFit$c * sin(0)
  mat <- matrix(gVals, ncol = 1, dimnames = list(c("BMAL1", "Per1", "Per2", "Cry1", "Cry2"), "s1"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  se <- estimatePhase(se, method = "timetable", tissue = "AliasTissue", species = "Mmusculus")
  cd <- SummarizedExperiment::colData(se)
  report <- S4Vectors::metadata(se)$circaTime_gene_mapping

  expect_false(is.na(cd$phase))
  expect_equal(as.character(cd$phase_flag), "ok")
  expect_equal(report$n_target_overlap_before, 4)
  expect_equal(report$n_target_overlap_after, 5)
  expect_equal(report$n_alias_or_id_rescued, 1)
})

test_that("estimatePhase(data_species=) ortholog-maps human-symbol input onto the bundled mouse model", {
  # AGENTS.md Section 4 Step 1.2.2: full round trip through mapOrthologs(), not a
  # mocked reference model - a real bundled model's own genes, relabelled to
  # their human orthologs, should still recover the model's own template phase.
  testthat::skip_if_not_installed("babelgene")
  avail <- availableReferenceModels()
  tis <- avail$tissue[avail$method == "timetable"][1]
  refModel <- circaTime:::.referenceModels[[paste("Mmusculus", tis, "timetable", "microarray", sep = ".")]]

  mat <- makeRoundTripMat(refModel, targetPhase = 6, nSamples = 1)
  ortho <- mapOrthologs(rownames(mat), fromSpecies = "Mmusculus", toSpecies = "Hsapiens")
  keep <- ortho$status[match(rownames(mat), ortho$input_gene)] == "ok"
  matHuman <- mat[keep, , drop = FALSE]
  rownames(matHuman) <- ortho$ortholog[match(rownames(matHuman), ortho$input_gene)]

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = matHuman))
  se <- estimatePhase(
    se, method = "timetable", tissue = tis, species = "Mmusculus", data_species = "Hsapiens"
  )
  cd <- SummarizedExperiment::colData(se)
  orthoRep <- S4Vectors::metadata(se)$circaTime_ortholog_mapping

  expect_false(is.na(cd$phase))
  expect_true(circularError(cd$phase, 6) < 1)
  expect_true(!is.null(orthoRep) && nrow(orthoRep) > 0)
})

# ---- method = "consensus" (AGENTS.md Section 6 Step 3.1) ----
#
# Uses injected fake models rather than the real bundled ones so the two
# constituent engines' phases can be set to known, deliberately-disagreeing
# values - with real models you cannot tell a correct weighted combination from a
# lucky one, because both engines agree closely on real data.

# Two fake timetable-shaped fits whose templates peak at different phases, so the
# engines genuinely disagree and the weighting has something to do.
makeConsensusFixture <- function(weights = c(timetable = 0.5, taufisher = 0.5)) {
  # `genes` drives the (patched, timetable-based) predict path; `Genes` is what
  # .predictWithReferenceModel()'s coverage rule reads for method "taufisher".
  # Both are set so the fake taufisher reproduces the real one's all-or-nothing
  # coverage behaviour instead of trivially passing on a NULL gene set.
  mkFit <- function(genes) {
    structure(
      list(genes = genes, Genes = toupper(genes), mesor = rep(5, 5),
           b = c(1, 0, -1, 0, 0.7), c = c(0, 1, 0, -1, 0.7),
           r = rep(0.9, 5), period = 24),
      class = "molecularTimetableFit"
    )
  }
  ampArgs <- list(
    ampCenter = stats::setNames(rep(0, 4), paste0("A", 1:4)),
    ampScale = stats::setNames(rep(1, 4), paste0("A", 1:4)),
    ampRotation = matrix(rep(0.5, 4), nrow = 4),
    ampGenes = paste0("A", 1:4), ampThreshold = 0
  )
  mk <- function(method, fit, cvMAE) {
    do.call(methods::new, c(
      list("PhaseReferenceModel", species = "Test", tissue = "T2", method = method,
           platform = "fake", fit = fit, residualPool = stats::rnorm(50, sd = 1),
           cvMAE = cvMAE, trainedOn = "unit test fixture"),
      ampArgs
    ))
  }
  list(
    # constituent engines are registered under real engine names so
    # .predictWithReferenceModel()'s switch() finds a gene-coverage rule
    "Test.T2.timetable.fake" = mk("timetable", mkFit(paste0("G", 1:5)), cvMAE = 1),
    "Test.T2.taufisher.fake" = mk("taufisher", mkFit(paste0("H", 1:5)), cvMAE = 3),
    "Test.T2.consensus_accuracy.fake" = mk(
      "consensus_accuracy",
      list(engines = c("timetable", "taufisher"), weights = weights), cvMAE = 1.5
    ),
    "Test.T2.consensus_equal.fake" = mk(
      "consensus_equal",
      list(engines = c("timetable", "taufisher"), weights = c(timetable = 0.5, taufisher = 0.5)),
      cvMAE = 2
    )
  )
}

# The taufisher fixture above carries a molecularTimetableFit, so point the
# taufisher engine slot at the timetable predict function for the duration of a
# consensus test - this keeps the fixture free of the real tauFisher package
# while still exercising the two-engine consensus code path.
withConsensusFixture <- function(models, code) {
  oldModels <- circaTime:::.referenceModels
  oldEngines <- circaTime:::.circaTimeEngines
  on.exit({
    utils::assignInNamespace(".referenceModels", oldModels, ns = "circaTime")
    utils::assignInNamespace(".circaTimeEngines", oldEngines, ns = "circaTime")
  }, add = TRUE)

  patchedEngines <- oldEngines
  patchedEngines$taufisher <- oldEngines$timetable
  utils::assignInNamespace(".circaTimeEngines", patchedEngines, ns = "circaTime")
  utils::assignInNamespace(".referenceModels", c(oldModels, models), ns = "circaTime")
  force(code)
}

# Expression matrix where engine A's genes (G*) encode phase 0 and engine B's
# genes (H*) encode phase 6, so the two engines return 0 and 6 respectively.
makeDisagreeingMat <- function() {
  tmpl <- function(theta) 5 + c(1, 0, -1, 0, 0.7) * cos(theta) + c(0, 1, 0, -1, 0.7) * sin(theta)
  mat <- rbind(
    matrix(tmpl(0), nrow = 5, dimnames = list(paste0("G", 1:5), "s1")),
    matrix(tmpl(6 / 24 * 2 * pi), nrow = 5, dimnames = list(paste0("H", 1:5), "s1")),
    matrix(10, nrow = 4, dimnames = list(paste0("A", 1:4), "s1"))
  )
  mat
}

test_that("consensus with equal weights lands on the circular mean of the engines", {
  withConsensusFixture(makeConsensusFixture(), {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(logExpr = makeDisagreeingMat())
    )
    out <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                          platform = "fake", consensus_weights = "equal")
    cd <- SummarizedExperiment::colData(out)
    rep <- S4Vectors::metadata(out)$circaTime_consensus

    # engines report 0 and 6 -> equal-weight circular mean is 3
    expect_equal(as.numeric(rep$per_engine_phase[1, ]), c(0, 6))
    expect_equal(round(cd$phase, 4), 3)
    expect_equal(unname(rep$effective_weights), c(0.5, 0.5))
    expect_equal(rep$n_engines_used, 2L)
  })
})

test_that("consensus with accuracy weights pulls toward the lower-cvMAE engine", {
  # timetable cvMAE = 1, taufisher cvMAE = 3 -> weights 0.75 / 0.25, so the
  # combination of 0 and 6 must land nearer 0 than the equal-weight answer of 3.
  w <- c(timetable = 0.75, taufisher = 0.25)
  withConsensusFixture(makeConsensusFixture(weights = w), {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(logExpr = makeDisagreeingMat())
    )
    out <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                          platform = "fake", consensus_weights = "accuracy")
    cd <- SummarizedExperiment::colData(out)
    rep <- S4Vectors::metadata(out)$circaTime_consensus

    expect_lt(cd$phase, 3)
    expect_gt(cd$phase, 0)
    expect_equal(unname(rep$stored_weights), c(0.75, 0.25))
  })
})

test_that("consensus exposes the weights, per-engine phases and agreement it used", {
  withConsensusFixture(makeConsensusFixture(), {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(logExpr = makeDisagreeingMat())
    )
    out <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                          platform = "fake", consensus_weights = "equal")
    rep <- S4Vectors::metadata(out)$circaTime_consensus

    expect_setequal(
      names(rep),
      c("engines", "stored_weights", "effective_weights", "gene_coverage_ok",
        "per_engine_phase", "concentration", "n_engines_used")
    )
    expect_equal(rep$engines, c("timetable", "taufisher"))
    expect_equal(colnames(rep$per_engine_phase), c("timetable", "taufisher"))
    # 0 vs 6 on a 24h circle is a 90-degree disagreement -> concentration cos(45deg)
    expect_equal(round(rep$concentration, 6), round(cos(pi / 4), 6))
    expect_true(all(rep$gene_coverage_ok))
  })
})

test_that("consensus drops an engine that lacks gene coverage and renormalises", {
  withConsensusFixture(makeConsensusFixture(), {
    mat <- makeDisagreeingMat()
    mat <- mat[!rownames(mat) %in% paste0("H", 1:5), , drop = FALSE] # kill engine B
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
    out <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                          platform = "fake", consensus_weights = "equal")
    cd <- SummarizedExperiment::colData(out)
    rep <- S4Vectors::metadata(out)$circaTime_consensus

    expect_equal(unname(rep$gene_coverage_ok), c(TRUE, FALSE))
    expect_equal(unname(rep$effective_weights), c(1, 0))
    expect_equal(rep$n_engines_used, 1L)
    expect_equal(round(cd$phase, 4), 0) # falls back to engine A's own answer
    expect_equal(as.character(cd$phase_flag), "ok") # one good engine is still coverage
  })
})

test_that("consensus flags gene_coverage only when no engine has coverage", {
  withConsensusFixture(makeConsensusFixture(), {
    mat <- matrix(1, nrow = 3, ncol = 1, dimnames = list(paste0("X", 1:3), "s1"))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
    out <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                          platform = "fake", consensus_weights = "equal")
    cd <- SummarizedExperiment::colData(out)
    expect_equal(as.character(cd$phase_flag), "gene_coverage")
    expect_true(is.na(cd$phase))
  })
})

test_that("consensus_weights picks a genuinely different stored model, not a runtime switch", {
  withConsensusFixture(makeConsensusFixture(), {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(logExpr = makeDisagreeingMat())
    )
    acc <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                          platform = "fake", consensus_weights = "accuracy")
    eq <- estimatePhase(se, method = "consensus", tissue = "T2", species = "Test",
                         platform = "fake", consensus_weights = "equal")
    # different residual pools -> different CI widths, proving each scheme carries
    # its own calibration rather than sharing one
    accModel <- circaTime:::.referenceModels[["Test.T2.consensus_accuracy.fake"]]
    eqModel <- circaTime:::.referenceModels[["Test.T2.consensus_equal.fake"]]
    expect_false(identical(accModel@residualPool, eqModel@residualPool))
    expect_false(identical(accModel@cvMAE, eqModel@cvMAE))
  })
})

test_that("estimatePhase errors clearly when no consensus model exists for a combination", {
  expect_error(
    estimatePhase(
      SummarizedExperiment::SummarizedExperiment(
        assays = list(logExpr = matrix(1, 3, 1, dimnames = list(paste0("g", 1:3), "s1")))
      ),
      method = "consensus", tissue = "Liv", species = "Mmusculus", platform = "rnaseq"
    ),
    "consensus_accuracy"
  )
})

test_that("estimatePhase rejects an unknown consensus_weights value", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logExpr = matrix(1, 3, 1, dimnames = list(paste0("g", 1:3), "s1")))
  )
  expect_error(
    estimatePhase(se, method = "consensus", tissue = "Liv", consensus_weights = "nope"),
    "should be one of"
  )
})

test_that("estimatePhase errors if data_species != species and map_orthologs = FALSE", {
  mat <- matrix(rnorm(10), nrow = 5, dimnames = list(paste0("g", 1:5), c("s1", "s2")))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  expect_error(
    estimatePhase(
      se, method = "timetable", tissue = "Liv", species = "Mmusculus",
      data_species = "Hsapiens", map_orthologs = FALSE
    ),
    "map_orthologs"
  )
})
