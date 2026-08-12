## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----libs, message = FALSE----------------------------------------------------
library(circaTime)
library(SummarizedExperiment)

## ----load-data----------------------------------------------------------------
se <- readRDS(system.file("extdata", "example_se.rds", package = "circaTime"))
se

## ----show-truth---------------------------------------------------------------
head(colData(se))

## ----avail--------------------------------------------------------------------
availableReferenceModels()

## ----estimate-----------------------------------------------------------------
phased <- estimatePhase(se, method = "timetable", tissue = "Liv")
colData(phased)[, c("phase", "phase_lo", "phase_hi", "clock_amp", "phase_flag")]

## ----consensus----------------------------------------------------------------
consensus <- estimatePhase(se, method = "consensus", tissue = "Liv")
colData(consensus)[, c("phase", "phase_lo", "phase_hi", "phase_flag")]

## ----consensus-report---------------------------------------------------------
report <- S4Vectors::metadata(consensus)$circaTime_consensus
report$stored_weights      # the trained prior, per engine
report$per_engine_phase    # what each engine said, per sample
round(report$concentration, 3)  # how strongly the engines agreed, in [0, 1]

## ----standardize, eval=requireNamespace("org.Mm.eg.db", quietly = TRUE)-------
mixedIds <- SummarizedExperiment(
  assays = list(logExpr = matrix(rnorm(9), nrow = 3,
    dimnames = list(c("Arntl", "ENSMUSG00000029238", "11461"), c("s1", "s2", "s3"))))
)
mixedIds <- standardizeGeneIds(mixedIds, species = "Mmusculus")
rownames(mixedIds) # Arntl (already a symbol), Clock (Ensembl), Actb (Entrez)

## ----orthologs, eval=requireNamespace("babelgene", quietly = TRUE)------------
mapOrthologs(c("Arntl", "Clock", "Per1"), fromSpecies = "Mmusculus", toSpecies = "Hsapiens")

## ----benchmark----------------------------------------------------------------
bm <- benchmarkPhase(se, truth_col = "phase", methods = "timetable", nFolds = 5, seed = 1)
bm$summary

## ----transfer, eval = FALSE---------------------------------------------------
# seMicroarray <- circaTimeData::circaTimeData("GSE54650") # same study, microarray
# seRnaseq <- circaTimeData::circaTimeData("GSE54651")     # same study, RNA-seq
# liverMicro <- seMicroarray[, SummarizedExperiment::colData(seMicroarray)$tissue == "Liv"]
# liverRnaseq <- seRnaseq[, SummarizedExperiment::colData(seRnaseq)$tissue == "Liv"]
# 
# tb <- transferPhase(liverMicro, liverRnaseq, truth_col = "phase", methods = "timetable")
# tb$summary

## ----bayesian-----------------------------------------------------------------
metadata(phased)$circaTime_bayesian
colData(phased)[, c("phase", "phase_lo", "phase_hi",
                    "phase_credible_lo", "phase_credible_hi")]

## ----dynamic-qc---------------------------------------------------------------
se_dyn <- estimatePhase(se, method = "timetable", tissue = "Liv",
                        qc_method = "dynamic")
table(colData(se_dyn)$phase_flag)
metadata(se_dyn)$circaTime_qc[c("method_used", "threshold", "panelSize")]

## ----downstream, eval = FALSE-------------------------------------------------
# # Illustrative only - not run as part of this vignette.
# library(limorhyde)
# colData(phased)$time_var <- limorhyde(colData(phased)$phase, colnames(colData(phased)))
# # ... proceed with limorhyde's own differential-rhythmicity workflow, e.g. via limma

## ----confounding--------------------------------------------------------------
# A deliberately confounded design: cases collected around ZT6, controls around ZT14.
set.seed(1)
n <- 20
confPhase <- c(rnorm(10, mean = 6, sd = 1), rnorm(10, mean = 14, sd = 1))
confMat <- matrix(rnorm(100 * n), nrow = 100,
                  dimnames = list(paste0("g", seq_len(100)), paste0("s", seq_len(n))))
conf <- SummarizedExperiment(assays = list(logExpr = confMat))
colData(conf)$phase <- wrapPhase(confPhase)
colData(conf)$condition <- rep(c("case", "control"), each = 10)

phaseConfounding(conf, ~ condition, nPerm = 999, seed = 1)

## ----confounding-balanced-----------------------------------------------------
balanced <- conf
colData(balanced)$phase <- rep(seq(0, 22, by = 2), length.out = n)

phaseConfounding(balanced, ~ condition, nPerm = 999, seed = 1)

## ----session-info-------------------------------------------------------------
sessionInfo()

