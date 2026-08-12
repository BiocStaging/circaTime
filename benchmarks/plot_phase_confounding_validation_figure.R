# Two-panel figure for Section 3.9 (validating phaseConfounding()'s type-I
# error and power via simulation). Reads only the CSV already written by
# run_phase_confounding_validation.R - no package load needed, base graphics
# only. Run from the project root:
#   Rscript benchmarks/plot_phase_confounding_validation_figure.R

df <- utils::read.csv("benchmarks/results/revision_phase_confounding_validation.csv")

nGrid <- c(10, 20, 50, 100)
distLevels <- c("uniform", "kappa2", "kappa5")
distLabels <- c("Uniform", "von Mises k=2", "von Mises k=5")
distCol <- c("steelblue4", "darkorange3", "firebrick")

t1 <- df[df$study == "type1_error" & df$n_groups == 2, ]

pw <- df[df$study == "power" & df$n_groups == 2 & df$phase_dist == "kappa2", ]
offsets <- sort(unique(pw$offset_hours))
nCol <- c("grey70", "steelblue3", "steelblue4", "black")

drawFigure <- function() {
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1.05))

  # --- Panel A: type-I error vs n/group, one line per phase distribution -----
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  graphics::plot(NA, xlim = range(nGrid), ylim = c(0, 0.09), log = "x",
    xlab = "n per group", ylab = "Rejection rate (null, alpha=0.05)",
    main = "A. Type-I error (2 groups)", cex.main = 0.95, xaxt = "n")
  graphics::axis(1, at = nGrid, labels = nGrid)
  graphics::abline(h = 0.05, lty = 2, col = "grey40")
  for (i in seq_along(distLevels)) {
    sub <- t1[t1$phase_dist == distLevels[i], ]
    sub <- sub[order(sub$n_per_group), ]
    graphics::lines(sub$n_per_group, sub$rejection_rate, col = distCol[i], lwd = 2, type = "b", pch = 16)
  }
  graphics::legend("topright", bty = "n", cex = 0.7, col = distCol, lwd = 2, pch = 16,
    legend = distLabels)

  # --- Panel B: power curves vs offset (h), one line per n, kappa=2 ----------
  graphics::par(mar = c(4.5, 4.2, 3, 1))
  graphics::plot(NA, xlim = range(offsets), ylim = c(0, 1),
    xlab = "Mean phase offset between groups (h)", ylab = "Rejection rate (power)",
    main = "B. Power (2 groups, von Mises k=2)", cex.main = 0.95)
  graphics::abline(h = 0.05, lty = 2, col = "grey40")
  for (i in seq_along(nGrid)) {
    sub <- pw[pw$n_per_group == nGrid[i], ]
    sub <- sub[order(sub$offset_hours), ]
    graphics::lines(sub$offset_hours, sub$rejection_rate, col = nCol[i], lwd = 2, type = "b", pch = 16)
  }
  graphics::legend("bottomright", bty = "n", cex = 0.7, col = nCol, lwd = 2, pch = 16,
    legend = paste0("n=", nGrid, "/group"))
}

grDevices::png("benchmarks/results/revision_phase_confounding_validation_figure.png", width = 1300, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
grDevices::pdf("paper/figures/fig_phase_confounding_validation.pdf", width = 1300 / 120, height = 550 / 120)
drawFigure()
grDevices::dev.off()
grDevices::png("paper/figures/fig_phase_confounding_validation.png", width = 1300, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

cat("Wrote benchmarks/results/revision_phase_confounding_validation_figure.png\n")
cat("Wrote paper/figures/fig_phase_confounding_validation.pdf and .png\n")
