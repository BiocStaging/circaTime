# Data-prep script (not shipped code, AGENTS.md Section 5, Phase 2 Step 2.1.1).
# Downloads GSE54651 (Zhang et al. 2014, 12-tissue mouse circadian atlas, RNA-seq
# arm) from GEO, builds a SummarizedExperiment with true CT/ZT in colData, and
# writes it to benchmarks/data/GSE54651_se.rds (gitignored - regenerate by
# re-running this script, do not commit the data itself).
#
# GSE54651 is the RNA-seq SubSeries of the same SuperSeries (GSE54652) as the
# microarray dataset already ingested by ingest_GSE54650.R (GSE54650) - same 12
# mouse tissues, same underlying study, profiled by Illumina HiSeq 2000 (GPL13112)
# instead of microarray. Design: 12 tissues x 8 timepoints (CT22-CT64, i.e. every 6h
# over 2 days), 1 sample per (tissue, CT) cell, no biological replicates. Verified
# on GEO 2026-08-08 (notes/decisions.md Section 10) - do not assume this dataset's
# shape without re-checking if GEO's listing ever changes.
#
# Unlike GSE54650 (probe-level microarray, needs a GPL annotation lookup), GSE54651
# ships one supplementary file per tissue, each a gene (Ensembl ID) x sample
# ("<Tissue>_CT<n>") table of GEO-author-computed "quant" values. GEO's metadata
# does not state the normalization method (RPKM/FPKM/TPM); values are treated as
# already-normalized, non-negative expression estimates, consistent with how
# ingest_GSE54650.R treats microarray intensities.
#
# Same donor/sample note as GSE54650: one sample per (tissue, CT) cell, so
# donor == sample here (see ingest_GSE54650.R header for the full rationale on
# why a `donor` column is still added).
#
# Ensembl ID -> gene symbol translation is delegated to circaTime's own
# standardizeGeneIds() (Phase 1 Step 1.1) rather than a hand-rolled AnnotationDbi
# call: org.Mm.eg.db's own SYMBOL field lags current HGNC/MGI nomenclature for at
# least one core clock gene (Arntl's Ensembl ID maps to legacy "Bmal1" via plain
# AnnotationDbi::mapIds(), verified while writing this script), and
# standardizeGeneIds()'s free layer already resolves the core clock genes'
# Ensembl IDs directly (bypassing that legacy-symbol problem entirely) before
# ever falling through to AnnotationDbi for the rest. Using the package's own
# tested function also guarantees this dataset's gene symbols land on the exact
# same convention as everything else circaTime produces.

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

tissues <- c("Adr", "Aor", "BFat", "Bstm", "Cer", "Hrt", "Hyp", "Kid", "Liv", "Lun", "Mus", "WFat")
baseUrl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54651/suppl/"

cacheDir <- file.path(tempdir(), "GSE54651_suppl")
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)

message("Downloading 12 per-tissue quantification files from GEO (this hits the network)...")
tissueTables <- lapply(tissues, function(tis) {
  fname <- sprintf("GSE54651_%s_quant.txt.gz", tis)
  dest <- file.path(cacheDir, fname)
  if (!file.exists(dest)) {
    utils::download.file(paste0(baseUrl, fname), dest, quiet = TRUE, mode = "wb")
  }
  tab <- utils::read.delim(dest, row.names = 1, check.names = FALSE)
  stopifnot(all(grepl(paste0("^", tis, "_CT[0-9]+$"), colnames(tab))))
  as.matrix(tab)
})
names(tissueTables) <- tissues

# All 12 files share the same Ensembl gene ID universe (same annotation build) but
# not the same row order - fail loudly rather than silently misaligning if the sets
# themselves ever differ, then reindex every table to one canonical (sorted) order.
geneIds <- sort(rownames(tissueTables[[1]]))
stopifnot(all(vapply(tissueTables, function(m) setequal(rownames(m), geneIds), logical(1))))
tissueTables <- lapply(tissueTables, function(m) m[geneIds, , drop = FALSE])

exprMat <- do.call(cbind, tissueTables)
message(sprintf(
  "%d genes x %d samples across %d tissues", nrow(exprMat), ncol(exprMat), length(tissues)
))

titleParts <- regmatches(
  colnames(exprMat),
  regexec("^([A-Za-z]+)_CT([0-9]+)$", colnames(exprMat))
)
stopifnot(all(lengths(titleParts) == 3)) # fail loudly if the column-naming format ever changes

tissue <- vapply(titleParts, `[`, character(1), 2)
ct <- as.numeric(vapply(titleParts, `[`, character(1), 3))

colData <- DataFrame(
  sample_title = colnames(exprMat),
  tissue = tissue,
  ct = ct, # raw circadian time, may exceed 24 (2-day design)
  phase = ct %% 24, # ground truth for benchmarking, in [0, 24)
  donor = colnames(exprMat), # see header note: sample == donor for this dataset
  row.names = colnames(exprMat)
)

seRaw <- SummarizedExperiment(
  # `expr` is the vendor-supplied "quant" value as deposited (normalization method
  # not stated in GEO metadata); `logExpr` is log2(expr + 1) - the +1 pseudocount is
  # needed here (unlike GSE54650's microarray intensities, which are always > 0)
  # because RNA-seq quantifications can be exactly 0 for unexpressed genes.
  assays = list(expr = exprMat, logExpr = log2(exprMat + 1)),
  colData = colData
)

message("Mapping Ensembl gene IDs -> gene symbols via standardizeGeneIds()...")
seMapped <- standardizeGeneIds(seRaw, species = "Mmusculus")
idReport <- S4Vectors::metadata(seMapped)$circaTime_id_mapping
mappedRate <- mean(idReport$resolved_by != "unresolved")

# standardizeGeneIds() keeps unmapped IDs under their original (Ensembl) identifier
# rather than dropping them (by design - see its docs), but a raw Ensembl ID sitting
# in an otherwise-all-symbols rowname vector is not a usable "gene" for anything
# downstream; drop those, matching ingest_GSE54650.R's policy of only keeping
# probes/genes that resolved to a real symbol.
seMapped <- seMapped[idReport$resolved_by != "unresolved", ]

# standardizeGeneIds() does not collapse rows that land on the same symbol (see
# its own docs) - multiple Ensembl IDs can map to one gene symbol, so collapse by
# row mean using the same internal helper estimatePhase() itself relies on.
exprGene <- circaTime:::.collapseDuplicateRows(assay(seMapped, "expr"))
logExprGene <- circaTime:::.collapseDuplicateRows(assay(seMapped, "logExpr"))
message(sprintf(
  "%d Ensembl IDs -> %d unique gene symbols (%d samples)",
  nrow(exprMat), nrow(exprGene), ncol(exprGene)
))

se <- SummarizedExperiment(
  assays = list(expr = exprGene, logExpr = logExprGene),
  colData = colData[colnames(exprGene), ]
)
metadata(se)$source <- "GSE54651"
metadata(se)$platform <- "GPL13112"
metadata(se)$mapping_rate <- mappedRate
metadata(se)$note <- paste(
  "Ground truth colData$phase = CT %% 24. One sample per (tissue, CT) cell, no",
  "biological replicates. RNA-seq SubSeries companion to GSE54650 (same study,",
  "same tissues) - see notes/decisions.md Section 10. Gene symbols standardized",
  "via circaTime::standardizeGeneIds(); full per-gene resolution report was in",
  "metadata(seMapped)$circaTime_id_mapping before duplicate collapse."
)

dir.create("benchmarks/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(se, "benchmarks/data/GSE54651_se.rds")
message("Saved benchmarks/data/GSE54651_se.rds (", nrow(se), " genes x ", ncol(se), " samples)")
