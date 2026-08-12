# Two-panel figure for Table~\ref{tab:dynamicqc} (dynamic coherence vs. the fixed
# 16-gene amplitude flag, GSE84511 aging cohort). Needs pkgload::load_all() for
# clockCoherence() and coreClockGenes(), so run under an env with the package's
# dependencies installed, e.g.:
#   /home/alperen/miniconda3/envs/circa/bin/Rscript benchmarks/plot_dynamic_qc_figure.R
# Run from the project root.

suppressMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

res <- utils::read.csv("benchmarks/results/phase3_dynamic_qc_GSE84511_persample.csv")
se <- readRDS("benchmarks/data/GSE84511_se.rds")
mat <- assay(se, "logExpr")[, res$sample, drop = FALSE]

# Same call as benchmarks/run_dynamic_qc_benchmark.R (identical seed/args), just
# to recover the attr(,"panel") gene list alongside the coherence scores already
# cached in the CSV.
coh <- clockCoherence(mat, phase = res$phase)
dynPanel <- attr(coh, "panel")[["all"]]
coreGenes <- coreClockGenes()$symbol
downstreamGenes <- dynPanel[!toupper(dynPanel) %in% toupper(coreGenes)]

# One core-clock gene and one dynamic-panel-only ("downstream") gene, picked as
# the pair with the cleanest real signal: highest mean cosinor R^2 fit to ZT
# within each candidate set (both real data, no cherry-picking beyond "fits a
# rhythm at all", which is what a reader needs to see the panels' rhythm).
fitR2 <- function(gene) {
  y <- mat[gene, ]
  theta <- 2 * pi * res$zt / 24
  fit <- stats::lm(y ~ cos(theta) + sin(theta))
  summary(fit)$r.squared
}
coreAvail <- coreGenes[toupper(coreGenes) %in% toupper(rownames(mat))]
coreAvail <- rownames(mat)[match(toupper(coreAvail), toupper(rownames(mat)))]
coreGene <- coreAvail[which.max(vapply(coreAvail, fitR2, numeric(1)))]
downGene <- downstreamGenes[which.max(vapply(downstreamGenes, fitR2, numeric(1)))]

cat(sprintf("Panel A genes: core=%s (dynamic panel size=%d, downstream candidates=%d)\n",
            coreGene, length(dynPanel), length(downstreamGenes)))
cat(sprintf("               downstream=%s\n", downGene))

zScore <- function(x) as.numeric(scale(x))

# ---- Panel B: real ROC curves from the cached per-sample scores ---------------
isBad <- res$error >= stats::quantile(res$error, 0.75)
rocCurve <- function(score, isBad) {
  # low score = bad by construction (see clockCoherence()/clockAmplitude() docs)
  ord <- order(score)
  isBad <- isBad[ord]
  tpr <- cumsum(isBad) / sum(isBad)
  fpr <- cumsum(!isBad) / sum(!isBad)
  list(fpr = c(0, fpr), tpr = c(0, tpr))
}
rocAmp <- rocCurve(res$amp, isBad)
rocCoh <- rocCurve(res$coherence, isBad)
aucTrap <- function(fpr, tpr) sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
aucAmp <- aucTrap(rocAmp$fpr, rocAmp$tpr)
aucCoh <- aucTrap(rocCoh$fpr, rocCoh$tpr)

drawFigure <- function() {
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))

  # --- Panel A: core-clock vs. downstream gene, young vs. aged, real ZT --------
  ageCol <- ifelse(res$age_group == "aged", "firebrick", "steelblue4")
  graphics::plot(res$zt, zScore(mat[coreGene, ]), col = ageCol, pch = 19,
    xlab = "ZT (h)", ylab = "Expression (z-score)",
    main = sprintf("A. Core clock (%s) vs. downstream (%s)", coreGene, downGene)
  )
  graphics::points(res$zt, zScore(mat[downGene, ]), col = ageCol, pch = 1)
  for (grp in c("aged", "young")) {
    col <- if (grp == "aged") "firebrick" else "steelblue4"
    idx <- res$age_group == grp
    ordz <- order(res$zt[idx])
    graphics::lines(res$zt[idx][ordz], zScore(mat[coreGene, ])[idx][ordz], col = col, lty = 1, lwd = 1.5)
    graphics::lines(res$zt[idx][ordz], zScore(mat[downGene, ])[idx][ordz], col = col, lty = 3, lwd = 1.5)
  }
  graphics::legend("topright", bty = "n", cex = 0.7,
    legend = c("aged", "young", "core clock (solid)", "downstream (dotted)"),
    col = c("firebrick", "steelblue4", "grey30", "grey30"),
    pch = c(19, 19, NA, NA), lty = c(NA, NA, 1, 3))

  # --- Panel B: ROC curves, fixed amplitude vs. dynamic coherence --------------
  graphics::plot(rocAmp$fpr, rocAmp$tpr, type = "l", lwd = 2, col = "firebrick",
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "False positive rate", ylab = "True positive rate",
    main = "B. Detecting worst-quartile-error samples"
  )
  graphics::lines(rocCoh$fpr, rocCoh$tpr, lwd = 2, col = "darkorange3")
  graphics::abline(0, 1, lty = 2, col = "grey50")
  graphics::legend("bottomright", bty = "n", cex = 0.8,
    legend = c(sprintf("fixed amplitude (AUC=%.2f, p=0.46)", aucAmp),
               sprintf("dynamic coherence (AUC=%.2f, p_sel=0.019)", aucCoh),
               "chance"),
    lty = c(1, 1, 2), lwd = c(2, 2, 1), col = c("firebrick", "darkorange3", "grey50"))
}

grDevices::png("benchmarks/results/phase3_dynamic_qc_figure.png", width = 1100, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
grDevices::pdf("paper/figures/fig_dynamic_qc.pdf", width = 1100 / 120, height = 550 / 120)
drawFigure()
grDevices::dev.off()
grDevices::png("paper/figures/fig_dynamic_qc.png", width = 1100, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

cat("Wrote benchmarks/results/phase3_dynamic_qc_figure.png\n")
cat("Wrote paper/figures/fig_dynamic_qc.pdf and .png\n")
cat(sprintf("AUC amplitude=%.3f  AUC coherence=%.3f\n", aucAmp, aucCoh))
