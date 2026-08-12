# Gene identifier normalization helpers. Two layers:
#   1. .normalizeGeneIds() - deterministic, dependency-free cleanup (Ensembl version
#      stripping, hand-curated clock-gene alias table). Used internally by
#      estimatePhase()/.prepareExpressionForModel() on every call, so it must stay
#      free of any Suggests dependency.
#   2. standardizeGeneIds() (below, AGENTS.md Section 4 Step 1.1) - the full
#      Ensembl/Entrez/RefSeq/symbol -> symbol translation layer for a
#      SummarizedExperiment, exported for direct use. Runs the free layer first;
#      anything still unmapped falls through to AnnotationDbi + org.Mm.eg.db/
#      org.Hs.eg.db, which are Suggests (notes/decisions.md Section 8) and only
#      loaded if installed.

.clockGeneAliases <- c(
  BMAL1 = "ARNTL",
  MOP3 = "ARNTL",
  BHLHE5 = "BHLHE41",
  SHARP1 = "BHLHE41",
  DEC1 = "BHLHE40",
  DEC2 = "BHLHE41",
  E4BP4 = "NFIL3",
  REV_ERBA = "NR1D1",
  REV_ERB_ALPHA = "NR1D1",
  REV_ERBB = "NR1D2",
  REV_ERB_BETA = "NR1D2",
  MOP4 = "NPAS2",
  CHRONO = "CIART",
  C1ORF51 = "CIART",
  GM129 = "CIART"
)

.mouseClockEnsembl <- c(
  ENSMUSG00000055116 = "Arntl",
  ENSMUSG00000029238 = "Clock",
  ENSMUSG00000026077 = "Npas2",
  ENSMUSG00000020893 = "Per1",
  ENSMUSG00000055866 = "Per2",
  ENSMUSG00000028957 = "Per3",
  ENSMUSG00000020038 = "Cry1",
  ENSMUSG00000068742 = "Cry2",
  ENSMUSG00000020889 = "Nr1d1",
  ENSMUSG00000021775 = "Nr1d2",
  ENSMUSG00000001270 = "Dbp",
  ENSMUSG00000020608 = "Tef",
  ENSMUSG00000003949 = "Hlf",
  ENSMUSG00000038542 = "Ciart",
  ENSMUSG00000030103 = "Bhlhe40",
  ENSMUSG00000023661 = "Bhlhe41"
)

.humanClockEnsembl <- c(
  ENSG00000133794 = "ARNTL",
  ENSG00000134852 = "CLOCK",
  ENSG00000170485 = "NPAS2",
  ENSG00000179094 = "PER1",
  ENSG00000132326 = "PER2",
  ENSG00000049246 = "PER3",
  ENSG00000008405 = "CRY1",
  ENSG00000121671 = "CRY2",
  ENSG00000126368 = "NR1D1",
  ENSG00000174738 = "NR1D2",
  ENSG00000105516 = "DBP",
  ENSG00000167074 = "TEF",
  ENSG00000108924 = "HLF",
  ENSG00000159208 = "CIART",
  ENSG00000134107 = "BHLHE40",
  ENSG00000123095 = "BHLHE41"
)

.stripEnsemblVersion <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

.canonicalForSpecies <- function(x, species) {
  if (identical(species, "Mmusculus")) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  } else {
    toupper(x)
  }
}

.normalizeGeneIds <- function(ids, species = "Mmusculus", map_aliases = TRUE) {
  raw <- as.character(ids)
  stripped <- .stripEnsemblVersion(raw)
  key <- toupper(gsub("[- .]", "_", stripped))

  ensMap <- switch(species,
    Mmusculus = .mouseClockEnsembl,
    Hsapiens = .humanClockEnsembl,
    c(.mouseClockEnsembl, .humanClockEnsembl)
  )

  out <- stripped
  ensHit <- key %in% names(ensMap)
  out[ensHit] <- unname(ensMap[key[ensHit]])

  if (isTRUE(map_aliases)) {
    aliasHit <- key %in% names(.clockGeneAliases)
    out[aliasHit] <- .canonicalForSpecies(unname(.clockGeneAliases[key[aliasHit]]), species)
  }

  out
}

.collapseDuplicateRows <- function(mat) {
  ids <- rownames(mat)
  if (!any(duplicated(ids))) return(mat)

  # Built via rbind of explicit 1-row matrices rather than vapply(..., numeric(ncol(mat)))
  # + t(): vapply silently drops to a plain vector (not a matrix) whenever
  # ncol(mat) == 1, which t() then transposes into the wrong orientation - a real
  # bug hit by single-sample estimatePhase() calls with an alias-induced duplicate
  # gene (found while testing standardizeGeneIds(), AGENTS.md Section 4 Step 1.1).
  groups <- split(seq_len(nrow(mat)), ids)
  rows <- lapply(groups, function(idx) {
    if (length(idx) == 1) {
      mat[idx, , drop = FALSE]
    } else {
      matrix(
        colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
        nrow = 1, dimnames = list(NULL, colnames(mat))
      )
    }
  })
  out <- do.call(rbind, rows)
  rownames(out) <- names(groups)
  out[unique(ids), , drop = FALSE] # restore original first-appearance row order
}

.geneMappingReport <- function(before, after, targetGenes) {
  targetUpper <- toupper(targetGenes)
  beforeOverlap <- sum(toupper(before) %in% targetUpper)
  afterOverlap <- sum(toupper(after) %in% targetUpper)
  data.frame(
    n_input_genes = length(before),
    n_unique_mapped_genes = length(unique(after)),
    n_target_genes = length(targetGenes),
    n_target_overlap_before = beforeOverlap,
    n_target_overlap_after = afterOverlap,
    target_overlap_rate = afterOverlap / length(targetGenes),
    n_alias_or_id_rescued = max(afterOverlap - beforeOverlap, 0),
    stringsAsFactors = FALSE
  )
}

.prepareExpressionForModel <- function(exprMat, targetGenes, species = "Mmusculus",
                                       map_aliases = TRUE) {
  before <- rownames(exprMat)
  after <- .normalizeGeneIds(before, species = species, map_aliases = map_aliases)
  rownames(exprMat) <- after
  exprMat <- .collapseDuplicateRows(exprMat)
  list(
    exprMat = exprMat,
    report = .geneMappingReport(before, after, targetGenes)
  )
}

# ---- AGENTS.md Section 4 Step 1.1: full ID-type detection + AnnotationDbi ----

#' @keywords internal
.detectIdType <- function(ids) {
  ids <- .stripEnsemblVersion(as.character(ids))
  out <- rep("symbol", length(ids))
  out[grepl("^ENS[A-Z]{0,4}G[0-9]{6,}$", ids)] <- "ensembl"
  out[grepl("^[0-9]+$", ids)] <- "entrez"
  out[grepl("^(NM|NR|XM|XR)_[0-9]+$", ids)] <- "refseq"
  out
}

.orgDbPackage <- function(species) {
  switch(species,
    Mmusculus = "org.Mm.eg.db",
    Hsapiens = "org.Hs.eg.db",
    stop(
      "standardizeGeneIds() only supports species = 'Mmusculus' or 'Hsapiens', got '",
      species, "'", call. = FALSE
    )
  )
}

# Maps non-symbol IDs to gene symbol via AnnotationDbi::mapIds(). Only called for
# the subset of rownames .normalizeGeneIds() (the free layer) did not already
# resolve, so this is skipped entirely for the common case of already-symbol,
# already-alias-resolved input. Returns NA_character_ for anything that fails to
# map or for which the org.*.eg.db package is not installed (caller keeps the
# original, un-mapped ID in that case rather than losing the gene entirely).
.annotationDbiFallback <- function(ids, idType, species) {
  out <- rep(NA_character_, length(ids))
  needsLookup <- idType != "symbol"
  if (!any(needsLookup)) return(out)

  orgPkg <- .orgDbPackage(species)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace(orgPkg, quietly = TRUE)) {
    return(out)
  }
  orgDb <- getExportedValue(orgPkg, orgPkg)

  keytypeFor <- c(ensembl = "ENSEMBL", entrez = "ENTREZID", refseq = "REFSEQ")
  for (kt in names(keytypeFor)) {
    idx <- needsLookup & idType == kt
    if (!any(idx)) next
    mapped <- tryCatch(
      AnnotationDbi::mapIds(
        orgDb, keys = .stripEnsemblVersion(ids[idx]), column = "SYMBOL",
        keytype = keytypeFor[[kt]], multiVals = "first"
      ),
      error = function(e) rep(NA_character_, sum(idx))
    )
    out[idx] <- unname(mapped)
  }
  out
}

# ---- AGENTS.md Section 4 Step 1.2: one-to-one ortholog mapping (babelgene) ----

.babelgeneSpeciesName <- c(Mmusculus = "mouse", Hsapiens = "human")

# Reverse lookup on .clockGeneAliases (alias -> canonical): given a canonical
# symbol, which aliases resolve to it? Used to work around babelgene's bundled
# symbol table lagging current HGNC symbols for at least one gene circaTime cares
# about directly - see mapOrthologs() below and notes/decisions.md Section 8.
.reverseClockAliases <- function(symbol) {
  hit <- toupper(symbol) == unname(.clockGeneAliases)
  names(.clockGeneAliases)[hit]
}

# Forward direction: if `symbol` is itself a known alias (e.g. babelgene's own
# "BMAL1"), resolve it to circaTime's canonical symbol ("ARNTL") before returning
# it from mapOrthologs() - otherwise the alias fix above would round-trip an
# ARNTL query back out as "BMAL1" instead of "ARNTL".
.canonicalizeClockAlias <- function(symbol) {
  key <- toupper(symbol)
  if (key %in% names(.clockGeneAliases)) unname(.clockGeneAliases[[key]]) else symbol
}

#' Map gene symbols to their one-to-one ortholog in another species
#'
#' Wraps `babelgene::orthologs()` (AGENTS.md Section 4 Step 1.2), restricted to
#' unambiguous one-to-one orthologs: a gene with zero or more than one candidate
#' ortholog is dropped (`ortholog = NA`, `status` records why), never guessed at
#' via a many-to-many mapping.
#'
#' `babelgene`'s bundled symbol table is a HomoloGene-vintage snapshot and can lag
#' current HGNC symbols - verified directly for a gene this package cares about:
#' the human ortholog of mouse `Arntl` is stored there under the legacy symbol
#' `"BMAL1"`, not the current official `"ARNTL"`, so querying `"ARNTL"` on its own
#' returns zero rows. To avoid that silently dropping known circadian genes, every
#' query also tries the gene's known aliases from `circaTime`'s own alias table
#' (the same one [standardizeGeneIds()] uses) in a single batched lookup, and
#' results are relabelled back to the originally requested symbol.
#'
#' @param genes Character vector of gene symbols in `fromSpecies`'s convention
#'   (`Mmusculus`: `Arntl`; `Hsapiens`: `ARNTL` - the same convention
#'   [standardizeGeneIds()] and [estimatePhase()] use).
#' @param fromSpecies,toSpecies `"Mmusculus"` or `"Hsapiens"`, must differ.
#' @return A `data.frame`, one row per unique input gene (input order preserved):
#'   `input_gene`, `ortholog` (`NA` if none), `support_n` (number of supporting
#'   ortholog databases, `NA` if none), `status` (`"ok"`, `"not_found"`, or
#'   `"ambiguous"`). A `warning()` reports how many genes were dropped and why.
#' @examples
#' \donttest{
#' mapOrthologs(c("Arntl", "Clock", "Per1"), fromSpecies = "Mmusculus", toSpecies = "Hsapiens")
#' }
#' @export
mapOrthologs <- function(genes, fromSpecies, toSpecies) {
  if (!fromSpecies %in% names(.babelgeneSpeciesName) ||
      !toSpecies %in% names(.babelgeneSpeciesName)) {
    stop(
      "fromSpecies/toSpecies must be 'Mmusculus' or 'Hsapiens', got '",
      fromSpecies, "'/'", toSpecies, "'", call. = FALSE
    )
  }
  if (identical(fromSpecies, toSpecies)) {
    stop("fromSpecies and toSpecies must differ", call. = FALSE)
  }
  .requireEngine("babelgene", 'install.packages("babelgene")')

  genes <- as.character(genes)
  uniqueGenes <- unique(genes)
  humanQuery <- identical(fromSpecies, "Hsapiens")
  babelSpecies <- .babelgeneSpeciesName[[if (humanQuery) toSpecies else fromSpecies]]
  fromCol <- if (humanQuery) "human_symbol" else "symbol"
  toCol <- if (humanQuery) "symbol" else "human_symbol"

  aliasExtra <- unlist(lapply(uniqueGenes, .reverseClockAliases), use.names = FALSE)
  queryGenes <- unique(c(uniqueGenes, aliasExtra))

  raw <- tryCatch(
    babelgene::orthologs(genes = queryGenes, species = babelSpecies, human = humanQuery, top = FALSE),
    error = function(e) NULL # babelgene errors (not NA/empty) when nothing at all matches
  )

  resolveOne <- function(g) {
    candidates <- c(g, .reverseClockAliases(g))
    hit <- if (is.null(raw)) raw else raw[raw[[fromCol]] %in% candidates, , drop = FALSE]
    n <- if (is.null(hit)) 0L else nrow(hit)
    if (n == 0) {
      list(ortholog = NA_character_, support_n = NA_real_, status = "not_found")
    } else if (n > 1) {
      list(ortholog = NA_character_, support_n = NA_real_, status = "ambiguous")
    } else {
      list(
        ortholog = .canonicalForSpecies(.canonicalizeClockAlias(hit[[toCol]][1]), toSpecies),
        support_n = hit$support_n[1], status = "ok"
      )
    }
  }

  results <- lapply(uniqueGenes, resolveOne)
  report <- data.frame(
    input_gene = uniqueGenes,
    ortholog = vapply(results, `[[`, character(1), "ortholog"),
    support_n = vapply(results, `[[`, numeric(1), "support_n"),
    status = vapply(results, `[[`, character(1), "status"),
    stringsAsFactors = FALSE
  )

  nDropped <- sum(report$status != "ok")
  if (nDropped > 0) {
    warning(
      nDropped, " of ", nrow(report), " genes have no one-to-one ortholog from ",
      fromSpecies, " to ", toSpecies, " (", sum(report$status == "not_found"),
      " not found, ", sum(report$status == "ambiguous"), " ambiguous) and were dropped ",
      "(kept as NA rather than guessed via a many-to-many mapping).",
      call. = FALSE
    )
  }
  report
}

#' Standardize gene identifiers in a SummarizedExperiment to gene symbol
#'
#' Auto-detects the identifier type of `rownames(se)` per gene (Ensembl, Entrez,
#' RefSeq, or already gene symbol - AGENTS.md Section 4 Step 1.1) and maps
#' everything to gene symbol, the internal identifier convention every other
#' circaTime function uses (`estimatePhase()`, `phaseConfounding()`,
#' `clockAmplitude()`). Runs in two layers: a dependency-free pass first (Ensembl
#' version-suffix stripping and the hand-curated circadian clock-gene alias table
#' also used internally by [estimatePhase()]), then, only for whatever remains
#' unresolved, `AnnotationDbi::mapIds()` against `org.Mm.eg.db`/`org.Hs.eg.db` if
#' installed (`Suggests`, not required for the common already-symbol case - see
#' `notes/decisions.md` Section 8).
#'
#' This does not collapse rows that map to the same symbol (e.g. two RefSeq
#' transcripts of one gene) - `rownames(se)` may contain duplicates afterward,
#' same as an unstandardized `SummarizedExperiment` can already have. Functions
#' that need a unique gene set (e.g. [estimatePhase()]'s internal
#' `.prepareExpressionForModel()`) already collapse duplicates by row mean at
#' that point.
#'
#' @param se A `SummarizedExperiment` with gene identifiers as `rownames(se)`.
#' @param species `"Mmusculus"` or `"Hsapiens"`.
#' @param map_aliases Logical, default `TRUE`: resolve circadian clock-gene
#'   aliases (e.g. `BMAL1` -> `Arntl`/`ARNTL`) in the free layer before falling
#'   through to `AnnotationDbi`.
#' @return `se` with `rownames(se)` replaced by standardized gene symbols.
#'   `metadata(se)$circaTime_id_mapping` is set to a one-row-per-input-gene
#'   `data.frame`: `input_id`, `id_type` (`"symbol"`/`"ensembl"`/`"entrez"`/
#'   `"refseq"`), `mapped_symbol`, `resolved_by` (`"unchanged"`, `"free_layer"`
#'   -- Ensembl version-stripping, the hardcoded core-clock-gene Ensembl table, or
#'   the alias table --, `"annotationdbi"`, or `"unresolved"`). A summary is also
#'   emitted via `message()`: how many genes were resolved, and by which layer.
#' @examples
#' mat <- matrix(rnorm(9), nrow = 3, dimnames = list(c("Arntl", "Bmal1", "Clock"), paste0("s", 1:3)))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(logExpr = mat))
#' se <- standardizeGeneIds(se, species = "Mmusculus")
#' rownames(se)
#' @export
standardizeGeneIds <- function(se, species = "Mmusculus", map_aliases = TRUE) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment", call. = FALSE)
  }
  ids <- rownames(se)
  if (is.null(ids)) {
    stop("se must have rownames (gene identifiers)", call. = FALSE)
  }

  idType <- .detectIdType(ids)
  free <- .normalizeGeneIds(ids, species = species, map_aliases = map_aliases)

  resolvedBy <- ifelse(free != .stripEnsemblVersion(as.character(ids)), "free_layer", "unchanged")
  # The free layer only ever changes an ID via the hardcoded core-clock Ensembl
  # table or the alias table (unrecognised Ensembl IDs pass through unchanged,
  # stripped of version only), so anything still equal to its own (stripped)
  # input AND not already a symbol needs the AnnotationDbi fallback.
  stillNeedsLookup <- resolvedBy == "unchanged" & idType != "symbol"

  final <- free
  if (any(stillNeedsLookup)) {
    fallback <- .annotationDbiFallback(ids[stillNeedsLookup], idType[stillNeedsLookup], species)
    gotFallback <- !is.na(fallback)
    idxLookup <- which(stillNeedsLookup)
    final[idxLookup[gotFallback]] <- .canonicalForSpecies(fallback[gotFallback], species)
    resolvedBy[idxLookup[gotFallback]] <- "annotationdbi"
    resolvedBy[idxLookup[!gotFallback]] <- "unresolved"
  }

  report <- data.frame(
    input_id = as.character(ids), id_type = idType, mapped_symbol = final,
    resolved_by = resolvedBy, stringsAsFactors = FALSE
  )

  nUnresolved <- sum(resolvedBy == "unresolved")
  nFree <- sum(resolvedBy == "free_layer")
  nAnnotationDbi <- sum(resolvedBy == "annotationdbi")
  message(
    "standardizeGeneIds(): ", length(ids), " genes; ", nFree, " resolved by the free layer, ",
    nAnnotationDbi, " by AnnotationDbi, ", nUnresolved, " unresolved (", length(ids) - nUnresolved,
    "/", length(ids), " = ", sprintf("%.1f%%", 100 * (length(ids) - nUnresolved) / length(ids)),
    " mapped)."
  )
  if (nUnresolved > 0) {
    warning(
      nUnresolved, " of ", length(ids), " gene identifiers could not be mapped to a symbol",
      if (any(idType[resolvedBy == "unresolved"] != "symbol") &&
          !requireNamespace(.orgDbPackage(species), quietly = TRUE)) {
        paste0(" (install ", .orgDbPackage(species), " for Ensembl/Entrez/RefSeq lookup)")
      } else {
        ""
      },
      ". They are kept under their original identifier rather than dropped.",
      call. = FALSE
    )
  }

  rownames(se) <- final
  oldMeta <- S4Vectors::metadata(se)
  oldMeta$circaTime_id_mapping <- report
  S4Vectors::metadata(se) <- oldMeta
  se
}
