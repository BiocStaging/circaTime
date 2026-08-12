# T6 (revision_and_fix.md): re-run T1's Table 1 significance tests (Friedman +
# Holm-corrected paired Wilcoxon + tissue-level bootstrap) on complete-case MAE
# instead of all-available MAE, using benchmarks/results/revision_complete_case_MAE.csv
# (produced by run_complete_case_benchmark.R). Mirrors run_significance_tests.R's
# Table 1 block exactly, just swapping the MAE column.
# Run from the project root: Rscript benchmarks/run_significance_tests_complete_case.R

set.seed(1)
nBoot <- 10000

bootCI <- function(x, nBoot = 10000, seed = 1) {
  set.seed(seed)
  n <- length(x)
  boots <- vapply(seq_len(nBoot), function(i) mean(sample(x, n, replace = TRUE)), numeric(1))
  stats::quantile(boots, c(0.025, 0.975), names = FALSE)
}

rows <- list()

t6 <- utils::read.csv("benchmarks/results/revision_complete_case_MAE.csv", stringsAsFactors = FALSE)
engines <- c("timetable", "zeitzeiger", "timesignatr", "taufisher")
mat1 <- reshape(t6[t6$engine %in% engines, c("tissue", "engine", "mae_complete_case")],
                idvar = "tissue", timevar = "engine", direction = "wide")
colnames(mat1) <- sub("^mae_complete_case\\.", "", colnames(mat1))
mat1 <- mat1[, c("tissue", engines)]
stopifnot(nrow(mat1) == 12, !anyNA(mat1))

fried1 <- stats::friedman.test(as.matrix(mat1[, engines]))
rows[[length(rows) + 1]] <- data.frame(
  comparison = "Table1_CC_all_engines", test = "friedman",
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
    comparison = paste0("Table1_CC_", p[1], "_vs_", p[2]), test = "wilcoxon_paired",
    statistic = NA, p_raw = pvals1[i], p_holm = pholm1[i],
    n_tissues = 12, effect_size_h = eff, ci_lo_h = NA, ci_hi_h = NA
  )
}

for (eng in engines) {
  ci <- bootCI(mat1[[eng]], nBoot)
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("Table1_CC_pooled_mean_", eng), test = "bootstrap_tissue_10000",
    statistic = NA, p_raw = NA, p_holm = NA, n_tissues = 12,
    effect_size_h = mean(mat1[[eng]]), ci_lo_h = ci[1], ci_hi_h = ci[2]
  )
}

out <- do.call(rbind, rows)
rownames(out) <- NULL
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "benchmarks/results/revision_significance_tests_complete_case.csv", row.names = FALSE)

cat("== T6: significance tests on complete-case MAE ==\n")
print(out, row.names = FALSE, digits = 4)
cat("\nWrote benchmarks/results/revision_significance_tests_complete_case.csv\n")
