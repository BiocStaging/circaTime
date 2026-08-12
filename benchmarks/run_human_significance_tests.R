# R6 (revision_and_fix_v2.md): significance testing for human cohort (GSE56931)
# Paired tests across the 14 subjects, each contributing its own mean circular error
# per engine. Friedman omnibus + pairwise Wilcoxon (Holm-corrected) + cluster bootstrap CIs
# (resample 14 subjects, compute mean per resample).
# Run from the project root: Rscript benchmarks/run_human_significance_tests.R

set.seed(1)
nBoot <- 10000

# Subject-level cluster bootstrap CI (resample 14 subjects with replacement)
clusterBootstrapCI <- function(subjMeans, nBoot = 10000, seed = 1) {
  set.seed(seed)
  boots <- vapply(seq_len(nBoot), function(i) {
    mean(sample(subjMeans, length(subjMeans), replace = TRUE))
  }, numeric(1))
  stats::quantile(boots, c(0.025, 0.975), names = FALSE)
}

# Read per-sample residuals and aggregate to subject level
resid <- utils::read.csv("benchmarks/results/revision_human_GSE56931_persample.csv",
                         stringsAsFactors = FALSE)
ok <- resid[!resid$fold_failed & !is.na(resid$pred), ]

methods <- c("timetable", "zeitzeiger", "timesignatr", "taufisher")
subjMat <- stats::reshape(
  ok[, c("donor", "method", "circular_error")],
  idvar = "donor", timevar = "method", direction = "wide"
)
colnames(subjMat) <- sub("^circular_error\\.", "", colnames(subjMat))
subjMat <- subjMat[, c("donor", methods)]
stopifnot(nrow(subjMat) == 14, !anyNA(subjMat[, methods]))

# Friedman test (blocking factor: subject; treatment: method)
fried <- stats::friedman.test(as.matrix(subjMat[, methods]))

# Pairwise Wilcoxon tests (paired, on the 14 subjects)
contrasts <- list(
  c("timetable", "zeitzeiger"),
  c("timetable", "timesignatr"),
  c("timetable", "taufisher"),
  c("zeitzeiger", "timesignatr"),
  c("zeitzeiger", "taufisher"),
  c("timesignatr", "taufisher")
)
pvals <- vapply(contrasts, function(p) {
  stats::wilcox.test(subjMat[[p[1]]], subjMat[[p[2]]], paired = TRUE, exact = FALSE)$p.value
}, numeric(1))
pholm <- stats::p.adjust(pvals, method = "holm")

# Output rows
rows <- list()
rows[[1]] <- data.frame(
  comparison = "Friedman_all_methods", test = "friedman",
  statistic = unname(fried$statistic), p_raw = fried$p.value, p_holm = NA,
  n_subjects = 14, effect_size_h = NA, ci_lo_h = NA, ci_hi_h = NA
)

for (i in seq_along(contrasts)) {
  p <- contrasts[[i]]
  eff <- mean(subjMat[[p[1]]] - subjMat[[p[2]]])
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("wilcoxon_", p[1], "_vs_", p[2]), test = "wilcoxon_paired",
    statistic = NA, p_raw = pvals[i], p_holm = pholm[i],
    n_subjects = 14, effect_size_h = eff, ci_lo_h = NA, ci_hi_h = NA
  )
}

# Cluster bootstrap CIs per method (on the 14 subjects)
for (m in methods) {
  ci <- clusterBootstrapCI(subjMat[[m]], nBoot)
  rows[[length(rows) + 1]] <- data.frame(
    comparison = paste0("bootstrap_", m), test = "cluster_bootstrap_subject_10000",
    statistic = NA, p_raw = NA, p_holm = NA,
    n_subjects = 14, effect_size_h = mean(subjMat[[m]]), ci_lo_h = ci[1], ci_hi_h = ci[2]
  )
}

out <- do.call(rbind, rows)
rownames(out) <- NULL
utils::write.csv(out, "benchmarks/results/revision_human_significance_tests.csv", row.names = FALSE)

cat("== R6: GSE56931 human significance tests (subject-level, n=14) ==\n")
print(out, row.names = FALSE, digits = 4)
cat("\nWrote benchmarks/results/revision_human_significance_tests.csv\n")
