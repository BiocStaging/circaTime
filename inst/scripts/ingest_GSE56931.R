# Data-prep script (not shipped code, AGENTS.md Section 5, Phase 2 Step 2.1.2).
# Downloads GSE56931 (human whole-blood circadian/sleep-deprivation transcriptome,
# Rosetta/Merck HuRSTA-2a520709 microarray) from GEO, builds a SummarizedExperiment
# with true clock time in colData, and writes it to benchmarks/data/GSE56931_se.rds
# (gitignored - regenerate by re-running this script, do not commit the data itself).
#
# GSE56931 design (verified on GEO 2026-08-08, notes/decisions.md Section 10): 249
# human whole-blood samples from 14 subjects (Diekelmann/Van Dongen-type protocol,
# PMID 25197809), platform GPL10379. Each subject was sampled every 4 hours across
# three protocol phases: a 24h normal-routine baseline (day 1-2, 6 timepoints), 38h
# of continuous wakefulness ("sleepdep", day 2-3), and a recovery-sleep period (day
# 4). Sample characteristics (`day`, `hour`, `patient`, `responder`, `Sex`, and
# `biological group` = protocol phase) are already parsed fields in the series
# matrix, not something this script has to extract from a title string. `hour` is
# clock time of day directly (values 0/4/8/12/16/20, not a cumulative protocol
# hour), i.e. already the [0, 24) phase circaTime wants.
#
# IMPORTANT caveat (carried into metadata/vignette, see notes/decisions.md Section
# 10): only the 130 baseline samples are unconfounded normal-routine circadian
# data. The 79 "sleepdep" and 40 "recovery" samples were collected under sustained
# sleep deprivation, which blunts/phase-shifts clock gene expression - they are
# NOT dropped here (excluding them is a modelling choice, not a data-cleaning one),
# but colData$protocol_phase lets downstream code subset to `== "baseline"` for a
# clean ground-truth cohort.
#
# Unlike GSE54650/GSE84511/GSE54651 (terminal harvest, one sample per animal), each
# of the 14 human subjects here contributes many samples across days - donor and
# sample do NOT coincide. `colData$donor` is set to `patient`, not `geo_accession`,
# so any cross-validation splitting code that groups by donor (AGENTS.md Section 9
# Phase 1 task 5's "split by donor, never by sample" rule) actually holds out whole
# subjects here, which matters because repeated samples from the same person are
# far more similar to each other than to a different person's samples.

suppressMessages({
  library(GEOquery)
  library(SummarizedExperiment)
  library(Biobase)
})

message("Downloading GSE56931 (this hits GEO over the network)...")
gse <- getGEO("GSE56931", GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]
stopifnot(is(gse, "ExpressionSet"))

exprMat <- Biobase::exprs(gse) # already RMA-normalized (log2 scale) - see data_processing field
pd <- Biobase::pData(gse)

requiredCols <- c("day:ch1", "hour:ch1", "patient:ch1", "responder:ch1", "Sex:ch1", "biological group:ch1")
stopifnot(all(requiredCols %in% colnames(pd))) # fail loudly if GEO's characteristics fields ever change

hour <- as.numeric(trimws(pd[["hour:ch1"]]))
stopifnot(all(hour >= 0 & hour < 24)) # confirms `hour` is clock time, not cumulative protocol hour

colData <- DataFrame(
  sample_title = pd$title,
  geo_accession = pd$geo_accession,
  patient = pd[["patient:ch1"]],
  day = as.integer(pd[["day:ch1"]]),
  hour = hour,
  phase = hour, # already in [0, 24) - see header note
  protocol_phase = factor(
    pd[["biological group:ch1"]],
    levels = c("baseline", "sleep deprivation", "recovery")
  ),
  responder = pd[["responder:ch1"]],
  sex = pd[["Sex:ch1"]],
  donor = pd[["patient:ch1"]], # see header note: sample != donor for this dataset
  row.names = pd$geo_accession
)

message("Fetching GPL10379 platform annotation for probe -> gene symbol mapping...")
gpl <- getGEO("GPL10379", AnnotGPL = TRUE)
annotTable <- Table(gpl)
probeToSymbol <- setNames(annotTable$`Gene symbol`, annotTable$ID)

geneSymbol <- probeToSymbol[rownames(exprMat)]
mappedRate <- mean(!is.na(geneSymbol) & geneSymbol != "")
message(sprintf("Probe -> gene symbol mapping rate: %.1f%%", 100 * mappedRate))

# Collapse multiple probes per gene by keeping the highest-variance probe, same
# policy as ingest_GSE54650.R.
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
  # The deposited values are already RMA-normalized log2 intensities (verified via
  # the series' own data_processing field) - `logExpr` is that value as-is, and
  # `expr` is its back-transform to linear scale, kept only for the same two-assay
  # convention as the other cohorts (it is not a raw intensity the way GSE54650's
  # `expr` is).
  assays = list(expr = 2^exprGene, logExpr = exprGene),
  colData = colData[colnames(exprGene), ]
)
metadata(se)$source <- "GSE56931"
metadata(se)$platform <- "GPL10379"
metadata(se)$mapping_rate <- mappedRate
metadata(se)$note <- paste(
  "Ground truth colData$phase = hour of day, already in [0, 24). Repeated-measures",
  "design: colData$donor = patient, NOT geo_accession - see script header. Only",
  "colData$protocol_phase == 'baseline' samples (130/249) are unconfounded",
  "normal-routine circadian data; 'sleep deprivation' and 'recovery' samples were",
  "collected during/after sustained sleep deprivation - see notes/decisions.md",
  "Section 10 before pooling them with baseline for training or benchmarking."
)

dir.create("benchmarks/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(se, "benchmarks/data/GSE56931_se.rds")
message("Saved benchmarks/data/GSE56931_se.rds (", nrow(se), " genes x ", ncol(se), " samples)")
