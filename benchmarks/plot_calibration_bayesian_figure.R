# Two-panel figure for the manuscript's calibration + Bayesian-density story.
# Reads only the CSVs already written by run_phase3_calibration.R and
# run_bayesian_ci_benchmark.R - no package load needed, base graphics only
# (matches plotCalibration()'s own style in R/confidence.R).
# Run from the project root: Rscript benchmarks/plot_calibration_bayesian_figure.R

calDf   <- utils::read.csv("benchmarks/results/phase3_calibration_GSE54650_timetable.csv")
calNaiveDf <- utils::read.csv("benchmarks/results/phase3_calibration_GSE54650_timetable_naive.csv")
tabDf   <- utils::read.csv("benchmarks/results/phase3_bayesian_vs_bootstrap_ci_GSE54650.csv")
residDf <- utils::read.csv("benchmarks/results/phase3_residuals_GSE54650_timetable.csv")

# Panel B example: one real Liver sample + its real bootstrap CI + the real
# von Mises kappa fitted to Liver's residual pool (Table 5's "Liv" row).
ex <- residDf[residDf$tissue == "Liv" & residDf$sample == "GSM1321202", ]
kappaLiv <- tabDf$kappa[tabDf$tissue == "Liv"]
period <- 24
point <- ex$phase

# von Mises density on a grid, exactly bayesianPhasePosterior()'s construction
# (R/confidence.R): f(t) ∝ exp(kappa * cos(2*pi*(t-point)/period)).
grid <- seq(0, period, length.out = 721)
dens <- exp(kappaLiv * cos(2 * pi * (grid - point) / period))
dens <- dens / (sum(dens) * mean(diff(grid))) # normalize to integrate to 1

# equal-tailed 80% credible interval from the same grid (nGrid-style ECDF inversion)
level <- 0.8
ord <- order(abs(grid - point))     # rank grid points by distance from the mode
cum <- cumsum(dens[ord] * mean(diff(grid)))
inSet <- ord[cum <= level]
bayesLo <- point - min(period / 2, max(abs(grid[inSet] - point)))
bayesHi <- point + min(period / 2, max(abs(grid[inSet] - point)))

bootLo <- ex$phase_lo
bootHi <- ex$phase_hi

# interval overlap (Jaccard) - both intervals are non-wrapping for this sample
interLo <- max(bayesLo, bootLo); interHi <- min(bayesHi, bootHi)
overlapLen <- max(0, interHi - interLo)
unionLen <- (bayesHi - bayesLo) + (bootHi - bootLo) - overlapLen
jaccard <- overlapLen / unionLen

drawFigure <- function() {
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))

  # --- Panel A: calibration curve ---
  graphics::plot(calDf$nominal, calDf$empirical,
    type = "b", pch = 19, col = "steelblue4", lwd = 2,
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "Nominal coverage", ylab = "Empirical coverage",
    main = "A. Calibration (GSE54650, timetable)"
  )
  graphics::abline(0, 1, lty = 2, col = "grey50")
  # Naive (no-residualPool) curve: full 5-level rerun of the same CV/bootstrap
  # design, via run_phase3_calibration_naive.R (benchmarks/results/*_naive.csv).
  graphics::lines(calNaiveDf$nominal, calNaiveDf$empirical, type = "b", pch = 4, col = "firebrick", lwd = 2)
  graphics::legend("bottomright", bty = "n", cex = 0.8,
    legend = c("residual-pool bootstrap (real)", "no residual pool (real)", "identity"),
    pch = c(19, 4, NA), lty = c(1, 1, 2), col = c("steelblue4", "firebrick", "grey50"))

  # --- Panel B: Bayesian posterior density with both intervals overlaid ---
  graphics::plot(grid, dens, type = "l", lwd = 2, col = "darkorange3",
    xlab = "Phase (h)", ylab = "Posterior density",
    main = sprintf("B. Liver sample %s (kappa=%.1f)", ex$sample, kappaLiv)
  )
  graphics::abline(v = point, lty = 3, col = "grey40")
  graphics::segments(bootLo, 0, bootHi, 0, lwd = 6, col = grDevices::adjustcolor("steelblue4", 0.6))
  graphics::segments(bayesLo, -graphics::par("usr")[4] * 0.02, bayesHi, -graphics::par("usr")[4] * 0.02,
    lwd = 6, col = grDevices::adjustcolor("darkorange3", 0.6), xpd = TRUE)
  graphics::legend("topright", bty = "n", cex = 0.8,
    legend = c("bootstrap 80% CI", "Bayesian 80% credible", sprintf("overlap (Jaccard) = %.2f", jaccard)),
    lty = c(1, 1, NA), lwd = c(6, 6, NA), col = c("steelblue4", "darkorange3", NA))
}

grDevices::png("benchmarks/results/phase3_calibration_bayesian_figure.png",
                width = 1100, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
grDevices::pdf("paper/figures/fig_calibration_bayesian.pdf", width = 1100 / 120, height = 550 / 120)
drawFigure()
grDevices::dev.off()
grDevices::png("paper/figures/fig_calibration_bayesian.png", width = 1100, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

cat("Wrote benchmarks/results/phase3_calibration_bayesian_figure.png\n")
cat("Wrote paper/figures/fig_calibration_bayesian.pdf and .png\n")
cat(sprintf("bootstrap CI: [%.2f, %.2f]h (width %.2fh)\n", bootLo, bootHi, bootHi - bootLo))
cat(sprintf("bayes CI:     [%.2f, %.2f]h (width %.2fh)\n", bayesLo, bayesHi, bayesHi - bayesLo))
cat(sprintf("Jaccard overlap: %.3f\n", jaccard))
