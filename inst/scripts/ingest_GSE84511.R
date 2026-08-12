# Data-prep script (not shipped code, AGENTS.md Section 5). Downloads GSE84511
# (Solanas G, Peixoto FO, et al. "Aged Stem Cells Reprogram Their Daily Rhythmic
# Functions to Adapt to Stress." Cell. 2017;170(4):678-692. PMID 28802040) from GEO
# and builds a SummarizedExperiment with true ZT in colData, written to
# benchmarks/data/GSE84511_se.rds (gitignored - regenerate by re-running this
# script, do not commit the data itself).
#
# Why this dataset: AGENTS.md Section 9's Phase 3 gate requires showing that
# low-amplitude-flagged samples have demonstrably worse phase error, and Section
# 7.4 expects this to matter "in aged tissue, tumours, and critically ill donors."
# GSE54650 (healthy adult mouse tissue atlas, no aged/damped-clock samples) could
# not test this - see benchmarks/results/phase3_amplitude_threshold_GSE54650_timetable.csv
# (p=0.083, not significant). GSE84511 was chosen after checking several
# aging/tumour/critical-illness candidates for (a) real ground-truth sample
# collection time and (b) public GEO availability - most either lacked one of
# these, or the source paper explicitly reports core clock GENE expression is
# "virtually unaltered" by age (only which downstream genes oscillate changes,
# not the core clock itself - true of both this study and its companion liver
# paper, Sato et al. 2017 Cell, GSE93903). This is a real, honest caveat: this
# run may reproduce GSE54650's weak result rather than resolve it. It is being
# tried because clockAmplitude() computes a per-SAMPLE PCA-based amplitude
# (capturing inter-individual heterogeneity), which is a different quantity from
# the per-gene population-average cosinor rhythmicity test the paper reports - so
# the paper's null result at the population level does not necessarily rule out
# a per-sample signal.
#
# GSE84511 design: mouse epidermal stem cells (FACS-sorted, pseudobulk), young (8
# weeks old) vs aged (102-116 weeks old), 6 ZT timepoints (ZT0,4,8,12,16,20), 4
# replicate animals per (age, ZT) cell, no repeated sampling of the same animal
# (terminal FACS harvest) - so donor == sample, same as GSE54650. Platform GPL11180
# (Affymetrix HT MG-430 PM), values already log2 (RMA-normalized) per GEO.

suppressMessages({
  library(GEOquery)
  library(SummarizedExperiment)
  library(Biobase)
})

message("Downloading GSE84511 (this hits GEO over the network)...")
gse <- getGEO("GSE84511", GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]
stopifnot(is(gse, "ExpressionSet"))

exprMat <- Biobase::exprs(gse) # already log2 (RMA) scale
pd <- Biobase::pData(gse)

stopifnot(all(c("age:ch1", "timepoint:ch1") %in% colnames(pd)))

ageGroup <- ifelse(pd[["age:ch1"]] == "8 week old", "young",
  ifelse(pd[["age:ch1"]] == "102-116 week old", "aged", NA_character_)
)
stopifnot(!anyNA(ageGroup))

ztParts <- regmatches(pd[["timepoint:ch1"]], regexec("^T([0-9]+)$", pd[["timepoint:ch1"]]))
stopifnot(all(lengths(ztParts) == 2)) # fail loudly if the timepoint format ever changes
zt <- as.numeric(vapply(ztParts, `[`, character(1), 2))

colData <- DataFrame(
  sample_title = pd$title,
  geo_accession = pd$geo_accession,
  tissue = "epidermal_stem_cells", # single tissue in this dataset
  age_group = ageGroup,
  zt = zt,
  phase = zt %% 24, # ground truth for benchmarking, in [0, 24)
  donor = pd$geo_accession, # terminal FACS harvest, no repeats - sample == donor
  row.names = pd$geo_accession
)

message("Fetching GPL11180 platform annotation for probe -> gene symbol mapping...")
gpl <- getGEO("GPL11180", AnnotGPL = TRUE)
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
exprGeneLog <- exprMapped[firstPerGene, , drop = FALSE]
rownames(exprGeneLog) <- symMapped[firstPerGene]

message(sprintf(
  "%d probes -> %d unique gene symbols (%d samples)",
  nrow(exprMat), nrow(exprGeneLog), ncol(exprGeneLog)
))

se <- SummarizedExperiment(
  # values as deposited are already log2 (RMA-normalized); `expr` reconstructs a
  # linear scale for parity with the GSE54650 ingest script's assay convention,
  # `logExpr` (what circaTime engines train/predict on) is the values as-is.
  assays = list(expr = 2^exprGeneLog, logExpr = exprGeneLog),
  colData = colData[colnames(exprGeneLog), ]
)
metadata(se)$source <- "GSE84511"
metadata(se)$platform <- "GPL11180"
metadata(se)$mapping_rate <- mappedRate
metadata(se)$note <- paste(
  "Ground truth colData$phase = ZT %% 24. 4 replicate animals per (age_group, ZT)",
  "cell, terminal FACS harvest - see script header for the donor/sample and",
  "amplitude-caveat notes."
)

dir.create("benchmarks/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(se, "benchmarks/data/GSE84511_se.rds")
message("Saved benchmarks/data/GSE84511_se.rds (", nrow(se), " genes x ", ncol(se), " samples)")
