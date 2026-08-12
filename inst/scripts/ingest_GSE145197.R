# Data-prep script (not shipped code, AGENTS.md Section 5, Phase 2 Step 2.1.3).
# Downloads GSE145197 (Droin et al. 2021, Nat Metab, PMID 33432202 - "Space-time
# logic of liver gene expression at sub-lobular scale") from GEO, aggregates its
# single-cell UMI counts into one PSEUDOBULK profile per sample (summing counts
# across all cells of that sample, not per-cell), and writes the result to
# benchmarks/data/GSE145197_se.rds (gitignored - regenerate by re-running this
# script, do not commit the data itself).
#
# GSE145197 design (verified on GEO 2026-08-08, notes/decisions.md Section 10): 10x
# Genomics scRNA-seq of mouse liver HEPATOCYTES ONLY (a sub-lobular zonation study,
# not a whole-liver atlas - there is exactly one cell type here, and that must stay
# explicit rather than implied), 10 samples (GSM4308343-GSM4308352, titles
# ZT00A/B/C, ZT06A/B, ZT12A/B/C, ZT18A/B) across 4 Zeitgeber timepoints (ZT0, ZT6,
# ZT12, ZT18; 2-3 replicate animals per timepoint), ~1400-3700 cells per sample.
# Each per-sample supplementary file is a gene x cell UMI count table, same 14812
# genes in the same order across all 10 files (verified directly, not assumed).
#
# Why pseudobulk, and why sum (not mean/median): AGENTS.md Step 2.1.3 asks for a
# pseudobulk dataset to test whether circaTime's engines - all of which are trained
# and validated on bulk RNA-seq/microarray - still recover circadian phase when
# applied to an aggregated single-cell profile standing in for a bulk sample.
# Summing UMI counts across a sample's cells is the standard pseudobulk
# construction (it reproduces what a real bulk RNA-seq library from the same cell
# population would approximate, unlike averaging normalized per-cell values, which
# does not correspond to any real sequencing quantity).
#
# Gene symbol casing: this series' gene names are deposited all-lowercase (e.g.
# "arntl", not "Arntl") - not a real naming convention, an artifact of the
# submitters' own processing. Fixed via circaTime's own species-canonical casing
# (the same Title-case Mmusculus convention standardizeGeneIds() and
# estimatePhase() already use) and then run through standardizeGeneIds() so any
# still-nonstandard clock-gene alias (e.g. "bmal1" -> "Bmal1" -> "Arntl") resolves
# the same way it would for any other dataset. Riken-clone-style identifiers (e.g.
# "0610007p14rik") cannot be case-recovered this way - only their first letter gets
# capitalized, not the internal mixed case MGI actually uses - but those are not
# genes this package's engines score against, so the residual case mismatch does
# not affect anything downstream.
#
# Terminal harvest, one GSM per animal (like GSE54650/GSE84511, unlike GSE56931) -
# donor and sample coincide.

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

gsmTitles <- c(
  GSM4308343 = "ZT00A", GSM4308344 = "ZT00B", GSM4308345 = "ZT00C",
  GSM4308346 = "ZT06A", GSM4308347 = "ZT06B",
  GSM4308348 = "ZT12A", GSM4308349 = "ZT12B", GSM4308350 = "ZT12C",
  GSM4308351 = "ZT18A", GSM4308352 = "ZT18B"
)
baseUrl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145197/suppl/"

cacheDir <- file.path(tempdir(), "GSE145197_suppl")
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
tarDest <- file.path(cacheDir, "GSE145197_RAW.tar")
if (!file.exists(tarDest)) {
  message("Downloading GSE145197_RAW.tar (this hits the network)...")
  utils::download.file(paste0(baseUrl, "GSE145197_RAW.tar"), tarDest, quiet = TRUE, mode = "wb")
}
utils::untar(tarDest, exdir = cacheDir)

message("Reading per-cell UMI tables and aggregating to pseudobulk (sum over cells)...")
pseudobulk <- lapply(names(gsmTitles), function(gsm) {
  matches <- list.files(cacheDir, pattern = paste0("^", gsm, "_"), full.names = TRUE)
  stopifnot(length(matches) == 1)
  # row.names = 1 fails here: a handful of gene names are duplicated within a single
  # file (e.g. "pisd" in GSM4308343) - read the gene column separately and collapse
  # duplicates by summing, rather than assuming row.names are unique like the other
  # ingest scripts' inputs are.
  tab <- utils::read.delim(matches, check.names = FALSE)
  geneCol <- tab[[1]]
  cellCounts <- as.matrix(tab[, -1, drop = FALSE])
  perCellRowSum <- rowSums(cellCounts)
  list(
    counts = unname(vapply(split(perCellRowSum, geneCol), sum, numeric(1))),
    genes = names(split(perCellRowSum, geneCol)),
    n_cells = ncol(cellCounts)
  )
})
names(pseudobulk) <- names(gsmTitles)

geneIds <- pseudobulk[[1]]$genes
stopifnot(all(vapply(pseudobulk, function(x) identical(x$genes, geneIds), logical(1))))

exprMat <- do.call(cbind, lapply(pseudobulk, `[[`, "counts"))
rownames(exprMat) <- geneIds
nCells <- vapply(pseudobulk, `[[`, integer(1), "n_cells")

message(sprintf(
  "%d genes x %d pseudobulk samples (%d-%d cells aggregated per sample)",
  nrow(exprMat), ncol(exprMat), min(nCells), max(nCells)
))

zt <- as.numeric(sub("^ZT([0-9]+)[A-Z]$", "\\1", gsmTitles))
stopifnot(!anyNA(zt))

colData <- DataFrame(
  sample_title = unname(gsmTitles),
  geo_accession = names(gsmTitles),
  zt = zt,
  phase = zt %% 24,
  cell_type = "hepatocyte",
  n_cells = nCells,
  donor = names(gsmTitles), # terminal harvest, one animal per GSM - see header note
  row.names = names(gsmTitles)
)

# Title-case gene names to circaTime's Mmusculus convention (see header note) so
# standardizeGeneIds()'s alias table can recognise anything case-flattened.
rownames(exprMat) <- circaTime:::.canonicalForSpecies(rownames(exprMat), "Mmusculus")

seRaw <- SummarizedExperiment(
  # `expr` is pseudobulk-summed raw UMI counts; `logExpr` is log2(expr + 1), same
  # pseudocount convention as ingest_GSE54651.R for the same reason (counts can be
  # exactly 0).
  assays = list(expr = exprMat, logExpr = log2(exprMat + 1)),
  colData = colData
)

message("Resolving any remaining clock-gene aliases via standardizeGeneIds()...")
seMapped <- suppressMessages(standardizeGeneIds(seRaw, species = "Mmusculus"))

exprGene <- circaTime:::.collapseDuplicateRows(assay(seMapped, "expr"))
logExprGene <- circaTime:::.collapseDuplicateRows(assay(seMapped, "logExpr"))
message(sprintf(
  "%d raw gene names -> %d unique gene symbols after alias resolution",
  nrow(exprMat), nrow(exprGene)
))

se <- SummarizedExperiment(
  assays = list(expr = exprGene, logExpr = logExprGene),
  colData = colData[colnames(exprGene), ]
)
metadata(se)$source <- "GSE145197"
metadata(se)$platform <- "10x Genomics (Illumina NextSeq 550)"
metadata(se)$note <- paste(
  "PSEUDOBULK profile per sample (summed UMI counts across all cells of that",
  "sample) - not per-cell data. Hepatocytes only (sub-lobular zonation study, not",
  "a whole-liver atlas) - see colData$cell_type. Ground truth colData$phase = ZT",
  "%% 24 (already < 24, no multi-day design here). 2-3 replicate animals per",
  "timepoint, terminal harvest (donor == sample). Only 4 ZT timepoints - coarser",
  "than the bulk cohorts, adequate for validating that pseudobulk aggregation",
  "preserves circadian signal but not for fine-grained phase regression - see",
  "notes/decisions.md Section 10."
)

dir.create("benchmarks/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(se, "benchmarks/data/GSE145197_se.rds")
message("Saved benchmarks/data/GSE145197_se.rds (", nrow(se), " genes x ", ncol(se), " samples)")
