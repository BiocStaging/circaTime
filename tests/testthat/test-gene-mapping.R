# Tests for R/gene-mapping.R: standardizeGeneIds() (AGENTS.md Section 4 Step 1.1)
# and its internal helpers. Ensembl/Entrez/RefSeq IDs used below (mouse Actb,
# human ACTB) are real, verified NCBI/Ensembl identifiers, not made up - see
# notes/decisions.md Section 8.

makeSE <- function(ids) {
  mat <- matrix(rnorm(length(ids) * 3), nrow = length(ids), dimnames = list(ids, paste0("s", 1:3)))
  SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
}

test_that(".detectIdType classifies all four ID types correctly", {
  ids <- c("Arntl", "ENSMUSG00000055116", "ENSG00000075624", "11461", "60", "NM_007393", "NM_007393.5")
  expect_equal(
    .detectIdType(ids),
    c("symbol", "ensembl", "ensembl", "entrez", "entrez", "refseq", "refseq")
  )
})

test_that("standardizeGeneIds() leaves already-correct symbols unchanged (mouse)", {
  se <- makeSE(c("Arntl", "Clock", "Per1"))
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"), "standardizeGeneIds")
  expect_equal(rownames(out), c("Arntl", "Clock", "Per1"))
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_true(all(rep$resolved_by == "unchanged"))
  expect_true(all(rep$id_type == "symbol"))
})

test_that("standardizeGeneIds() leaves already-correct symbols unchanged (human)", {
  se <- makeSE(c("ARNTL", "CLOCK", "PER1"))
  expect_message(out <- standardizeGeneIds(se, species = "Hsapiens"))
  expect_equal(rownames(out), c("ARNTL", "CLOCK", "PER1"))
})

test_that("standardizeGeneIds() resolves clock-gene aliases via the free layer", {
  se <- makeSE(c("BMAL1", "MOP3"))
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  expect_equal(rownames(out), c("Arntl", "Arntl"))
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_true(all(rep$resolved_by == "free_layer"))
})

test_that("standardizeGeneIds() resolves core-clock Ensembl IDs without AnnotationDbi", {
  se <- makeSE(c("ENSMUSG00000055116", "ENSMUSG00000029238")) # Arntl, Clock
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  expect_equal(rownames(out), c("Arntl", "Clock"))
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_true(all(rep$id_type == "ensembl"))
  expect_true(all(rep$resolved_by == "free_layer"))
})

test_that("standardizeGeneIds() Ensembl version suffixes are stripped before matching", {
  se <- makeSE(c("ENSMUSG00000055116.15"))
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  expect_equal(rownames(out), "Arntl")
})

test_that("standardizeGeneIds() resolves additional verified circadian aliases (NPAS2, CIART)", {
  # MOP4 -> NPAS2 and CHRONO/C1orf51/Gm129 -> CIART are real HGNC/MGI-documented
  # aliases/previous symbols (notes/decisions.md Section 8), not made up.
  se <- makeSE(c("MOP4", "CHRONO", "C1ORF51", "GM129"))
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  expect_equal(rownames(out), c("Npas2", "Ciart", "Ciart", "Ciart"))
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_true(all(rep$resolved_by == "free_layer"))
})

test_that("standardizeGeneIds() falls through to AnnotationDbi for a non-clock Ensembl ID", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Mm.eg.db")
  se <- makeSE("ENSMUSG00000029580") # mouse Actb, not in the hardcoded clock table
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  expect_equal(rownames(out), "Actb")
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_equal(rep$resolved_by, "annotationdbi")
})

test_that("standardizeGeneIds() falls through to AnnotationDbi for Entrez IDs (both species)", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Mm.eg.db")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  seMouse <- makeSE("11461") # mouse Actb Entrez ID
  expect_message(outMouse <- standardizeGeneIds(seMouse, species = "Mmusculus"))
  expect_equal(rownames(outMouse), "Actb")

  seHuman <- makeSE("60") # human ACTB Entrez ID
  expect_message(outHuman <- standardizeGeneIds(seHuman, species = "Hsapiens"))
  expect_equal(rownames(outHuman), "ACTB")
})

test_that("standardizeGeneIds() falls through to AnnotationDbi for RefSeq IDs", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Mm.eg.db")
  se <- makeSE(c("NM_007393", "NM_007393.5")) # mouse Actb, with and without version
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  expect_equal(rownames(out), c("Actb", "Actb"))
})

test_that("standardizeGeneIds() reports unresolved IDs, warns, and keeps the original", {
  se <- makeSE(c("Arntl", "9999999999")) # second ID: not a real Entrez gene
  expect_warning(
    expect_message(out <- standardizeGeneIds(se, species = "Mmusculus")),
    "could not be mapped"
  )
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_equal(rep$resolved_by[rep$input_id == "9999999999"], "unresolved")
  expect_equal(rownames(out)[2], "9999999999")
})

test_that("standardizeGeneIds() report has one row per input gene with expected columns", {
  se <- makeSE(c("Arntl", "BMAL1", "ENSMUSG00000055116"))
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus"))
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_equal(nrow(rep), 3)
  expect_equal(
    colnames(rep),
    c("input_id", "id_type", "mapped_symbol", "resolved_by")
  )
})

test_that("standardizeGeneIds() rejects non-SummarizedExperiment input", {
  expect_error(standardizeGeneIds(matrix(1:4, 2, 2), species = "Mmusculus"), "SummarizedExperiment")
})

test_that("standardizeGeneIds() rejects an se with no rownames", {
  mat <- matrix(rnorm(6), nrow = 2)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
  expect_error(standardizeGeneIds(se, species = "Mmusculus"), "rownames")
})

test_that("standardizeGeneIds() rejects an unsupported species", {
  se <- makeSE(c("11461")) # a non-symbol ID forces the species check to run
  expect_error(standardizeGeneIds(se, species = "Drosophila"), "Mmusculus.*Hsapiens")
})

# ---- mapOrthologs() (AGENTS.md Section 4 Step 1.2) ----

test_that("mapOrthologs() maps known mouse<->human circadian gene orthologs", {
  testthat::skip_if_not_installed("babelgene")
  out <- mapOrthologs(c("Arntl", "Clock", "Per1"), fromSpecies = "Mmusculus", toSpecies = "Hsapiens")
  expect_equal(out$status, c("ok", "ok", "ok"))
  expect_equal(out$ortholog, c("ARNTL", "CLOCK", "PER1"))
})

test_that("mapOrthologs() resolves Arntl<->ARNTL despite babelgene's legacy BMAL1 symbol", {
  # Empirically verified (notes/decisions.md Section 8): babelgene's own
  # human_symbol column for this gene's ortholog is "BMAL1", not the current HGNC
  # symbol "ARNTL" - querying "ARNTL" directly returns zero rows from babelgene.
  # This is exactly the AGENTS.md Section 4 Step 1.2.4 mandated test case.
  testthat::skip_if_not_installed("babelgene")
  mouseToHuman <- mapOrthologs("Arntl", fromSpecies = "Mmusculus", toSpecies = "Hsapiens")
  expect_equal(mouseToHuman$ortholog, "ARNTL")
  expect_equal(mouseToHuman$status, "ok")

  humanToMouse <- mapOrthologs("ARNTL", fromSpecies = "Hsapiens", toSpecies = "Mmusculus")
  expect_equal(humanToMouse$ortholog, "Arntl")
  expect_equal(humanToMouse$status, "ok")
})

test_that("mapOrthologs() drops genes with no ortholog and warns", {
  testthat::skip_if_not_installed("babelgene")
  expect_warning(
    out <- mapOrthologs(c("Arntl", "NotARealGeneSymbolXYZ"), fromSpecies = "Mmusculus", toSpecies = "Hsapiens"),
    "no one-to-one ortholog"
  )
  expect_equal(out$status, c("ok", "not_found"))
  expect_true(is.na(out$ortholog[2]))
})

test_that("mapOrthologs() returns one row per unique input gene, input order preserved", {
  testthat::skip_if_not_installed("babelgene")
  out <- mapOrthologs(c("Per1", "Arntl", "Per1"), fromSpecies = "Mmusculus", toSpecies = "Hsapiens")
  expect_equal(nrow(out), 2)
  expect_equal(out$input_gene, c("Per1", "Arntl"))
})

test_that("mapOrthologs() rejects same fromSpecies/toSpecies and unsupported species", {
  expect_error(mapOrthologs("Arntl", fromSpecies = "Mmusculus", toSpecies = "Mmusculus"), "differ")
  expect_error(mapOrthologs("Arntl", fromSpecies = "Mmusculus", toSpecies = "Drosophila"), "Mmusculus.*Hsapiens")
})

test_that("standardizeGeneIds() map_aliases = FALSE skips the alias table", {
  # BMAL1 is a symbol-type ID (not Ensembl/Entrez/RefSeq), so it never reaches the
  # AnnotationDbi fallback either way - with the alias layer off, it is left
  # exactly as given, silently (a symbol standardizeGeneIds() doesn't recognise as
  # a known alias is not necessarily wrong, just not translated).
  se <- makeSE("BMAL1")
  expect_message(out <- standardizeGeneIds(se, species = "Mmusculus", map_aliases = FALSE))
  expect_equal(rownames(out), "BMAL1")
  rep <- S4Vectors::metadata(out)$circaTime_id_mapping
  expect_equal(rep$resolved_by, "unchanged")
})
