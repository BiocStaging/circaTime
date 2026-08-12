# Data-prep script (AGENTS.md Section 5) that builds inst/extdata/example_se.rds,
# the self-contained example SummarizedExperiment vignettes/circaTime.Rmd runs
# against. Vignettes must build without a network download (Bioconductor's build
# system has none), so this is a real, once-off subset of GSE54650 bundled with the
# package rather than something the vignette re-downloads itself.
#
# Subset: liver ("Liv"), all 24 samples, full gene set (22105 genes) - realistic
# (a real analyst passes whatever genes they have, not a hand-curated list) and
# still small (~2.1 MB compressed), comfortably inside the < 5 MB inst/extdata
# budget (AGENTS.md Section 5) alongside the existing core_clock_genes.tsv.
#
# Run from the project root: Rscript inst/scripts/make_vignette_example.R
# Requires benchmarks/data/GSE54650_se.rds (run inst/scripts/ingest_GSE54650.R first).

suppressMessages(library(SummarizedExperiment))

seFile <- "benchmarks/data/GSE54650_se.rds"
if (!file.exists(seFile)) {
  stop("Run inst/scripts/ingest_GSE54650.R first to produce ", seFile)
}
se <- readRDS(seFile)
exampleSe <- se[, colData(se)$tissue == "Liv"]

metadata(exampleSe)$note <- paste(
  "Subset of GSE54650 (Zhang et al. 2014 PNAS) for vignette/example use:",
  "liver tissue, all 24 samples, full gene set. See",
  "inst/scripts/make_vignette_example.R."
)

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(exampleSe, "inst/extdata/example_se.rds", compress = "xz")
cat(sprintf(
  "Saved inst/extdata/example_se.rds (%d genes x %d samples, %.1f KB)\n",
  nrow(exampleSe), ncol(exampleSe), file.size("inst/extdata/example_se.rds") / 1024
))
