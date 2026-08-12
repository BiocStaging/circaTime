# Data-prep script (not shipped code, AGENTS.md Section 5). Downloads GSE54650
# (Zhang et al. 2014, 12-tissue mouse circadian atlas) from GEO, builds a
# SummarizedExperiment with true CT/ZT in colData, and writes it to
# benchmarks/data/GSE54650_se.rds (gitignored - regenerate by re-running this
# script, do not commit the data itself).
#
# GSE54650 design: 12 mouse tissues x 24 timepoints (CT18-CT64, i.e. every 2h over
# 2 days), 1 array per (tissue, CT) cell, no biological replicates. Platform GPL6246
# (Affymetrix Mouse Gene 1.0 ST). Sample titles encode tissue + CT directly, e.g.
# "Adr_CT18". CT is truncated to [0, 24) to get circadian phase; the underlying CT
# value (which day) is kept as a separate colData column in case day-effects ever
# matter, but is not part of the phase-inference ground truth.
#
# Because each (tissue, CT) cell has exactly one array from a (presumably)
# independent animal, "donor" and "sample" coincide here - there is no repeated-
# measures structure within a tissue to violate AGENTS.md Section 9 Phase 1 task 5's
# "split by donor, never by sample" rule. A `donor` colData column is still added
# (one distinct donor ID per sample) so downstream benchmark code can apply the same
# donor-aware splitting logic uniformly across datasets that do have repeats.

suppressMessages({
  library(GEOquery)
  library(SummarizedExperiment)
  library(Biobase)
})

message("Downloading GSE54650 (this hits GEO over the network)...")
gse <- getGEO("GSE54650", GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]
stopifnot(is(gse, "ExpressionSet"))

exprMat <- Biobase::exprs(gse)
pd <- Biobase::pData(gse)

titleParts <- regmatches(
  pd$title,
  regexec("^([A-Za-z]+)_CT([0-9]+)$", pd$title)
)
stopifnot(all(lengths(titleParts) == 3)) # fail loudly if the title format ever changes

tissue <- vapply(titleParts, `[`, character(1), 2)
ct <- as.numeric(vapply(titleParts, `[`, character(1), 3))

colData <- DataFrame(
  sample_title = pd$title,
  geo_accession = pd$geo_accession,
  tissue = tissue,
  ct = ct, # raw circadian time, may exceed 24 (2-day design)
  phase = ct %% 24, # ground truth for benchmarking, in [0, 24)
  donor = pd$geo_accession, # see header note: sample == donor for this dataset
  row.names = pd$geo_accession
)

message("Fetching GPL6246 platform annotation for probe -> gene symbol mapping...")
gpl <- getGEO("GPL6246", AnnotGPL = TRUE)
annotTable <- Table(gpl)
probeToSymbol <- setNames(annotTable$`Gene symbol`, annotTable$ID)

geneSymbol <- probeToSymbol[rownames(exprMat)]
mappedRate <- mean(!is.na(geneSymbol) & geneSymbol != "")
message(sprintf("Probe -> gene symbol mapping rate: %.1f%%", 100 * mappedRate))

# Collapse multiple probes per gene by keeping the highest-variance probe
# (standard practice; avoids arbitrarily averaging away real oscillatory signal).
keep <- !is.na(geneSymbol) & geneSymbol != ""
exprMapped <- exprMat[keep, , drop = FALSE]
symMapped <- geneSymbol[keep]

probeVar <- apply(exprMapped, 1, stats::var)
ord <- order(symMapped, -probeVar)
exprMapped <- exprMapped[ord, , drop = FALSE]
symMapped <- symMapped[ord]
firstPerGene <- !duplicated(symMapped)
exprGene <- exprMapped[firstPerGene, , drop = FALSE]
rownames(exprGene) <- symMapped[firstPerGene]

message(sprintf(
  "%d probes -> %d unique gene symbols (%d samples)",
  nrow(exprMat), nrow(exprGene), ncol(exprGene)
))

se <- SummarizedExperiment(
  # `expr` is the raw (linear-scale) probeset intensity as deposited on GEO; `logExpr`
  # is log2(expr) and is what circaTime engines should train/predict on - harmonic
  # regression assumes roughly-sinusoidal behaviour on a scale where fold-changes are
  # additive, which holds in log-space, not on raw linear intensities.
  assays = list(expr = exprGene, logExpr = log2(exprGene)),
  colData = colData[colnames(exprGene), ]
)
metadata(se)$source <- "GSE54650"
metadata(se)$platform <- "GPL6246"
metadata(se)$mapping_rate <- mappedRate
metadata(se)$note <- paste(
  "Ground truth colData$phase = CT %% 24. One array per (tissue, CT) cell,",
  "no biological replicates - see script header for the donor/sample note."
)

dir.create("benchmarks/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(se, "benchmarks/data/GSE54650_se.rds")
message("Saved benchmarks/data/GSE54650_se.rds (", nrow(se), " genes x ", ncol(se), " samples)")
