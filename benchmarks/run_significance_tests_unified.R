# T11 (revision_and_fix.md): re-run T1's Table 4 significance-testing block on
# the unified-run consensus numbers (benchmarks/results/revision_consensus_GSE54650.csv,
# produced by run_unified_benchmark.R from the same nFolds=6/seed=1 CV run Table 1
# uses), since run_unified_benchmark.R changed the point estimates from the old
# separate nFolds=5 consensus run. Mirrors run_significance_tests.R's Table 4 block
# exactly, just pointed at the new source file.
# Run from the project root: Rscript benchmarks/run_significance_tests_unified.R

set.seed(1)
nBoot <- 10000

bootCI <- function(x, nBoot = 10000, seed = 1) {
  set.seed(seed)
  n <- length(x)
  boots <- vapply(seq_len(nBoot), function(i) mean(sample(x, n, replace = TRUE)), numeric(1))
  stats::quantile(boots, c(0.025, 0.975), names = FALSE)
}

rows <- list()

t4 <- utils::read.csv("benchmarks/results/revision_consensus_GSE54650.csv", stringsAsFactors = FALSE)
methods4 <- c("timetable", "taufisher", "consensus_equal", "consensus_accuracy")
mat4 <- stats::reshape(t4[t4$method %in% methods4, c("tissue", "method", "circular_mae")],
                       idvar = "tissue", timevar = "method", direction = "wide")
colnames(mat4) <- sub("^circular_mae\\.", "", colnames(mat4))
mat4 <- mat4[, c("tissue", methods4)]
stopifnot(nrow(mat4) == 12, !anyNA(mat4))

fried4 <- stats::friedman.test(as.matrix(mat4[, methods4]))
rows[[length(rows) + 1]] <- data.frame(
  comparison = "Table4_all_methods", test = "friedman",
  statistic = unname(fried4$statistic), p_raw = fried4$p.value, p_holm = NA,
  n_tissues = 12, effect_size_h = NA, ci_lo_h = NA, ci_hi_h = NA
)

# Best single point estimate under the unified run is Molecular Timetable, not
# accuracy-weighted consensus (see run_unified_benchmark.R output) - contrasts
# are against the best point-estimate performer, timetable, not assumed to still
# be consensus_accuracy.
contrasts4 <- list(
  c("timetable", "consensus_accuracy"),
  c("timetable", "consensus_equal"),
  c("timetable", "taufisher"),
  c("consensus_accuracy", "consensus_equal")
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
utils::write.csv(out, "benchmarks/results/revision_significance_tests_unified.csv", row.names = FALSE)

cat("== T11: Table 4 significance tests, unified run ==\n")
print(out, row.names = FALSE, digits = 4)
cat("\nWrote benchmarks/results/revision_significance_tests_unified.csv\n")
