# T1 (revision_and_fix.md): statistical testing for every accuracy comparison
# reported in the manuscript. Uses only already-committed per-tissue summary
# files (benchmarks/results/phase2_benchmark_GSE54650.csv,
# benchmarks/results/phase3_consensus_GSE54650.csv) - no CV is re-run, because
# the Friedman/Wilcoxon tests and the tissue-level bootstrap CIs both operate
# on one MAE value per (tissue, method) cell, which those files already have.
# Run from the project root: Rscript benchmarks/run_significance_tests.R

set.seed(1)
nBoot <- 10000
alpha <- 0.05

bootCI <- function(x, nBoot = 10000, seed = 1) {
  set.seed(seed)
  n <- length(x)
  boots <- vapply(seq_len(nBoot), function(i) mean(sample(x, n, replace = TRUE)), numeric(1))
  stats::quantile(boots, c(0.025, 0.975), names = FALSE)
}

rows <- list()

## ---- Table 1: four engines x 12 tissues (GSE54650) --------------------------
t1 <- utils::read.csv("benchmarks/results/phase2_benchmark_GSE54650.csv", stringsAsFactors = FALSE)
engines <- c("timetable", "zeitzeiger", "timesignatr", "taufisher")
mat1 <- reshape(t1[t1$method %in% engines, c("tissue", "method", "circular_mae")],
                idvar = "tissue", timevar = "method", direction = "wide")
colnames(mat1) <- sub("^circular_mae\\.", "", colnames(mat1))
mat1 <- mat1[, c("tissue", engines)]
stopifnot(nrow(mat1) == 12, !anyNA(mat1))

fried1 <- stats::friedman.test(as.matrix(mat1[, engines]))
rows[[length(rows) + 1]] <- data.frame(
  comparison = "Table1_all_engines", test = "friedman",
  statistic = unname(fried1$statistic), p_raw = fried1$p.value, p_holm = NA,
  n_tissues = 12, effect_size_h = NA, ci_lo_h = NA, ci_hi_h = NA
)

pairs1 <- utils::combn(engines, 2, simplify = FALSE)
pvals1 <- vapply(pairs1, function(p) {
  stats::wilcox.test(mat1[[p[1]]], mat1[[p[2]]], paired = TRUE, exact = FALSE)$p.value
}, numeric(1))
pholm1 <- stats::p.adjust(pvals1, method = "holm")
for (i in seq_along(pairs1)) {
  p <- pairs1[[i]]
  eff <- mean(mat1[[p[1]]] - mat1[[p[2]]])
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("Table1_", p[1], "_vs_", p[2]), test = "wilcoxon_paired",
    statistic = NA, p_raw = pvals1[i], p_holm = pholm1[i],
    n_tissues = 12, effect_size_h = eff, ci_lo_h = NA, ci_hi_h = NA
  )
}

for (eng in engines) {
  ci <- bootCI(mat1[[eng]], nBoot)
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("Table1_pooled_mean_", eng), test = "bootstrap_tissue_10000",
    statistic = NA, p_raw = NA, p_holm = NA, n_tissues = 12,
    effect_size_h = mean(mat1[[eng]]), ci_lo_h = ci[1], ci_hi_h = ci[2]
  )
}

## ---- Table 4: consensus vs single engines (GSE54650) -------------------------
t4 <- utils::read.csv("benchmarks/results/phase3_consensus_GSE54650.csv", stringsAsFactors = FALSE)
methods4 <- c("timetable", "taufisher", "consensus_equal", "consensus_accuracy")
mat4 <- reshape(t4[t4$method %in% methods4, c("tissue", "method", "circular_mae")],
                idvar = "tissue", timevar = "method", direction = "wide")
colnames(mat4) <- sub("^circular_mae\\.", "", colnames(mat4))
mat4 <- mat4[, c("tissue", methods4)]
stopifnot(nrow(mat4) == 12, !anyNA(mat4))

contrasts4 <- list(
  c("consensus_accuracy", "consensus_equal"),
  c("consensus_accuracy", "timetable"),
  c("consensus_accuracy", "taufisher")
)
pvals4 <- vapply(contrasts4, function(p) {
  stats::wilcox.test(mat4[[p[1]]], mat4[[p[2]]], paired = TRUE, exact = FALSE)$p.value
}, numeric(1))
pholm4 <- stats::p.adjust(pvals4, method = "holm")
for (i in seq_along(contrasts4)) {
  p <- contrasts4[[i]]
  eff <- mean(mat4[[p[1]]] - mat4[[p[2]]])
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("Table4_", p[1], "_vs_", p[2]), test = "wilcoxon_paired",
    statistic = NA, p_raw = pvals4[i], p_holm = pholm4[i],
    n_tissues = 12, effect_size_h = eff, ci_lo_h = NA, ci_hi_h = NA
  )
}

for (m in methods4) {
  ci <- bootCI(mat4[[m]], nBoot)
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("Table4_pooled_mean_", m), test = "bootstrap_tissue_10000",
    statistic = NA, p_raw = NA, p_holm = NA, n_tissues = 12,
    effect_size_h = mean(mat4[[m]]), ci_lo_h = ci[1], ci_hi_h = ci[2]
  )
}

out <- do.call(rbind, rows)
rownames(out) <- NULL
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "benchmarks/results/revision_significance_tests.csv", row.names = FALSE)

cat("== T1: significance tests ==\n")
print(out, row.names = FALSE, digits = 4)
cat("\nWrote benchmarks/results/revision_significance_tests.csv\n")
