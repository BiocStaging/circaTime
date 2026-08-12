# Two-panel figure for Table~\ref{tab:transfer} / Section 3.8 (cross-platform
# transfer, all-or-nothing gene-panel failure). Reads only the CSV already written
# by run_transfer_benchmark.R - no package load needed, base graphics only.
# Run from the project root: Rscript benchmarks/plot_transfer_figure.R

df <- utils::read.csv("benchmarks/results/phase2_transfer_benchmark_GSE54650_GSE54651.csv")

tissues <- c("Adr", "Aor", "Bstm", "BFat", "Cer", "Hrt", "Hyp", "Kid", "Liv", "Lun", "Mus", "WFat")
directions <- c("microarray_to_rnaseq", "rnaseq_to_microarray")
dirLabels <- c("Micro->RNA-seq", "RNA-seq->Micro")
methods <- c("timetable", "taufisher", "zeitzeiger", "timesignatr")
methodLabels <- c("Molecular\nTimetable", "tauFisher", "ZeitZeiger", "TimeSignature")

# ---- Panel A matrix: Table 7's real timetable circular MAE, tissue x direction --
maeMat <- matrix(NA_real_, nrow = length(tissues), ncol = length(directions),
                  dimnames = list(tissues, dirLabels))
for (d in seq_along(directions)) {
  sub <- df[df$method == "timetable" & df$direction == directions[d], ]
  maeMat[sub$tissue, d] <- sub$circular_mae
}

# ---- Panel B: real survival counts - tissues with a non-NA prediction, /12 -----
survMat <- matrix(0L, nrow = length(methods), ncol = length(directions),
                   dimnames = list(methodLabels, dirLabels))
for (m in seq_along(methods)) {
  for (d in seq_along(directions)) {
    sub <- df[df$method == methods[m] & df$direction == directions[d], ]
    survMat[m, d] <- sum(!is.na(sub$circular_mae))
  }
}

drawFigure <- function() {
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

  # --- Panel A: bidirectional transfer heatmap (Table 7, Molecular Timetable) --
  graphics::par(mar = c(5, 6, 3, 1))
  pal <- grDevices::colorRampPalette(c("steelblue4", "white", "firebrick"))(100)
  rng <- range(maeMat, na.rm = TRUE)
  graphics::image(seq_len(ncol(maeMat)), seq_len(nrow(maeMat)), t(maeMat[nrow(maeMat):1, ]),
    col = pal, zlim = rng, axes = FALSE,
    xlab = "", ylab = "", main = "A. Transfer MAE (h), Molecular Timetable", cex.main = 0.95
  )
  graphics::axis(1, at = seq_len(ncol(maeMat)), labels = dirLabels, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(nrow(maeMat)), labels = rev(tissues), las = 1, cex.axis = 0.75)
  for (i in seq_len(nrow(maeMat))) {
    for (j in seq_len(ncol(maeMat))) {
      val <- maeMat[nrow(maeMat) - i + 1, j]
      graphics::text(j, i, sprintf("%.1f", val), cex = 0.65)
    }
  }
  graphics::box()

  # --- Panel B: all-or-nothing robustness under gene-panel mismatch ------------
  graphics::par(mar = c(6, 4.2, 3, 1))
  bp <- graphics::barplot(t(survMat), beside = TRUE,
    col = c("steelblue4", "firebrick"), ylim = c(0, 12.9),
    ylab = "Tissues with a prediction (of 12)",
    main = "B. Robustness under train/score platform mismatch",
    names.arg = methodLabels, cex.names = 0.75, las = 1, cex.main = 0.95
  )
  graphics::text(bp, t(survMat) + 0.4, labels = t(survMat), cex = 0.75, xpd = TRUE)
  graphics::legend("topright", bty = "n", cex = 0.75, fill = c("steelblue4", "firebrick"),
    legend = dirLabels)
}

grDevices::png("benchmarks/results/phase2_transfer_figure.png", width = 1300, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
grDevices::pdf("paper/figures/fig_transfer.pdf", width = 1300 / 120, height = 550 / 120)
drawFigure()
grDevices::dev.off()
grDevices::png("paper/figures/fig_transfer.png", width = 1300, height = 550, res = 120)
drawFigure()
grDevices::dev.off()

cat("Wrote benchmarks/results/phase2_transfer_figure.png\n")
cat("Wrote paper/figures/fig_transfer.pdf and .png\n")
