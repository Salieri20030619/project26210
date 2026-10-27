#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 02_validate_phase103_nci60_two_level_analysis_v1.R <all_table> <renal_table> <analysis_dir>")
all_dt <- fread(normalizePath(args[[1]], mustWork = TRUE))
renal <- fread(normalizePath(args[[2]], mustWork = TRUE))
analysis_dir <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
results_dir <- file.path(analysis_dir, "results_v1")
qc_dir <- file.path(analysis_dir, "qc_v1")

primary <- fread(file.path(results_dir, "all_lineages_primary_associations_v1.tsv"))
adjusted <- fread(file.path(results_dir, "all_lineages_tissue_adjusted_v1.tsv"))
kendall <- fread(file.path(results_dir, "all_lineages_kendall_floor_sensitivity_v1.tsv"))
mutual <- fread(file.path(results_dir, "all_lineages_mutually_adjusted_rank_model_v1.tsv"))
above <- fread(file.path(results_dir, "all_lineages_above_floor_sensitivity_v1.tsv"))
tissue_out <- fread(file.path(results_dir, "tissue_specific_effects_v1.tsv"))
decisions <- fread(file.path(results_dir, "all_lineages_decisions_v1.tsv"))
renal_out <- fread(file.path(results_dir, "renal_exploratory_effects_v1.tsv"))
renal_loo_out <- fread(file.path(results_dir, "renal_leave_one_out_v1.tsv"))
renal_line_out <- fread(file.path(results_dir, "renal_line_level_v1.tsv"))

spec <- data.table(
  relation = c("TYMP-MTA", "MTAP-MTA"),
  x_col = c("tymp_log2_fpkm_plus1", "mtap_log2_fpkm_plus1"),
  i = 1:2
)

boot_plain <- function(x, y, seed, b = 20000L) {
  set.seed(seed)
  n <- length(x)
  draws <- numeric(b)
  for (k in seq_len(b)) {
    take <- sample.int(n, n, replace = TRUE)
    draws[k] <- suppressWarnings(cor(rank(x[take]), rank(y[take])))
  }
  unname(quantile(draws[is.finite(draws)], c(0.025, 0.975), type = 7))
}

within_effect <- function(x, y, g) {
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  rx <- rx - ave(rx, g, FUN = mean)
  ry <- ry - ave(ry, g, FUN = mean)
  cor(rx, ry)
}

boot_strata <- function(x, y, g, seed, b = 20000L) {
  set.seed(seed)
  strata <- split(seq_along(x), g)
  draws <- numeric(b)
  for (k in seq_len(b)) {
    take <- unlist(lapply(strata, function(v) sample(v, length(v), replace = TRUE)), use.names = FALSE)
    draws[k] <- within_effect(x[take], y[take], g[take])
  }
  unname(quantile(draws[is.finite(draws)], c(0.025, 0.975), type = 7))
}

permute_within <- function(x, y, g, seed, b = 100000L) {
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  rx <- rx - ave(rx, g, FUN = mean)
  ry <- ry - ave(ry, g, FUN = mean)
  obs <- cor(rx, ry)
  strata <- split(seq_along(x), g)
  set.seed(seed)
  extreme <- 0L
  for (k in seq_len(b)) {
    order_idx <- seq_along(x)
    for (v in strata) order_idx[v] <- sample(v)
    r <- cor(rx[order_idx], ry)
    if (abs(r) >= abs(obs)) extreme <- extreme + 1L
  }
  (extreme + 1) / (b + 1)
}

re_primary <- rbindlist(lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k]
  x <- all_dt[[s$x_col]]
  y <- all_dt$mta_mean_raw
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  ci <- boot_plain(x, y, 103100L + s$i)
  data.table(relation = s$relation, rho = unname(test$estimate), p = test$p.value, ci_low = ci[1], ci_high = ci[2])
}))
re_primary[, q := p.adjust(p, "BH")]

re_adjusted <- rbindlist(lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k]
  x <- all_dt[[s$x_col]]
  y <- all_dt$mta_mean_raw
  ci <- boot_strata(x, y, all_dt$tissue, 103200L + s$i)
  data.table(relation = s$relation, effect = within_effect(x, y, all_dt$tissue),
             ci_low = ci[1], ci_high = ci[2], p = permute_within(x, y, all_dt$tissue, 103300L + s$i))
}))
re_adjusted[, q := p.adjust(p, "BH")]

re_kendall <- rbindlist(lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k]
  test <- suppressWarnings(cor.test(all_dt[[s$x_col]], all_dt$mta_mean_raw, method = "kendall", exact = FALSE))
  data.table(relation = s$relation, tau = unname(test$estimate), p = test$p.value)
}))
re_kendall[, q := p.adjust(p, "BH")]

af <- all_dt[mta_mean_raw > min(all_dt$mta_mean_raw)]
re_above <- rbindlist(lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k]
  test <- suppressWarnings(cor.test(af[[s$x_col]], af$mta_mean_raw, method = "spearman", exact = FALSE))
  data.table(relation = s$relation, rho = unname(test$estimate), p = test$p.value)
}))
re_above[, q := p.adjust(p, "BH")]

model_dt <- copy(all_dt)
model_dt[, y := as.numeric(scale(rank(mta_mean_raw)))]
model_dt[, tx := as.numeric(scale(rank(tymp_log2_fpkm_plus1)))]
model_dt[, mx := as.numeric(scale(rank(mtap_log2_fpkm_plus1)))]
model_dt[, tissue := factor(tissue)]
fit <- lm(y ~ tx + mx + tissue, model_dt)
co <- summary(fit)$coefficients
re_mutual <- data.table(
  relation = c("TYMP-MTA", "MTAP-MTA"),
  beta = c(co["tx", "Estimate"], co["mx", "Estimate"]),
  p = c(co["tx", "Pr(>|t|)"], co["mx", "Pr(>|t|)"])
)
re_mutual[, q := p.adjust(p, "BH")]

re_tissue <- rbindlist(lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k]
  rbindlist(lapply(split(all_dt, all_dt$tissue), function(z) {
    if (nrow(z) < 3L) return(NULL)
    data.table(relation = s$relation, tissue = unique(z$tissue), n = nrow(z),
               rho = suppressWarnings(cor(z[[s$x_col]], z$mta_mean_raw, method = "spearman")))
  }))
}))

re_renal <- rbindlist(lapply(seq_len(nrow(spec)), function(k) {
  s <- spec[k]
  x <- renal[[s$x_col]]
  y <- renal$mta_mean_raw
  loo <- vapply(seq_along(x), function(j) suppressWarnings(cor(x[-j], y[-j], method = "spearman")), numeric(1))
  data.table(relation = s$relation, rho = suppressWarnings(cor(x, y, method = "spearman")),
             tau = suppressWarnings(cor(x, y, method = "kendall")), loo_min = min(loo), loo_max = max(loo),
             same_sign_n = sum(sign(loo) == sign(suppressWarnings(cor(x, y, method = "spearman")))))
}))

tol <- 1e-12
same <- function(a, b, tolerance = tol) all(abs(a - b) <= tolerance)
checks <- data.table(check_id = character(), pass = logical(), detail = character())
add <- function(id, pass, detail) checks <<- rbind(checks, data.table(check_id = id, pass = pass, detail = detail))

pjoin <- merge(primary, re_primary, by = "relation")
add("primary_rho_recomputed", same(pjoin$spearman_rho, pjoin$rho), paste(max(abs(pjoin$spearman_rho - pjoin$rho))))
add("primary_p_recomputed", same(pjoin$p_value, pjoin$p), paste(max(abs(pjoin$p_value - pjoin$p))))
add("primary_q_recomputed", same(pjoin$fdr_bh, pjoin$q), paste(max(abs(pjoin$fdr_bh - pjoin$q))))
add("primary_bootstrap_ci_recomputed", same(pjoin$bootstrap_ci_low, pjoin$ci_low) && same(pjoin$bootstrap_ci_high, pjoin$ci_high), "20,000 deterministic resamples")

ajoin <- merge(adjusted, re_adjusted, by = "relation")
add("adjusted_effect_recomputed", same(ajoin$tissue_adjusted_partial_spearman, ajoin$effect), paste(max(abs(ajoin$tissue_adjusted_partial_spearman - ajoin$effect))))
add("adjusted_permutation_p_recomputed", same(ajoin$permutation_p_value, ajoin$p), paste(max(abs(ajoin$permutation_p_value - ajoin$p))))
add("adjusted_q_recomputed", same(ajoin$fdr_bh, ajoin$q), paste(max(abs(ajoin$fdr_bh - ajoin$q))))
add("adjusted_bootstrap_ci_recomputed", same(ajoin$stratified_bootstrap_ci_low, ajoin$ci_low) && same(ajoin$stratified_bootstrap_ci_high, ajoin$ci_high), "20,000 tissue-stratified resamples")

kjoin <- merge(kendall, re_kendall, by = "relation")
add("kendall_recomputed", same(kjoin$kendall_tau_b, kjoin$tau) && same(kjoin$p_value, kjoin$p) && same(kjoin$fdr_bh, kjoin$q), "tau/P/q")
mjoin <- merge(mutual, re_mutual, by = "relation")
add("mutual_rank_model_recomputed", same(mjoin$standardized_rank_beta, mjoin$beta) && same(mjoin$p_value, mjoin$p) && same(mjoin$fdr_bh, mjoin$q), "beta/P/q")
bjoin <- merge(above, re_above, by = "relation")
add("above_floor_sensitivity_recomputed", same(bjoin$spearman_rho, bjoin$rho) && same(bjoin$p_value, bjoin$p) && same(bjoin$fdr_bh, bjoin$q), "rho/P/q")

tjoin <- merge(tissue_out, re_tissue, by = c("relation", "tissue"))
add("tissue_effects_recomputed", nrow(tjoin) == nrow(tissue_out) && same(tjoin$spearman_rho, tjoin$rho), paste0("rows=", nrow(tjoin)))
rjoin <- merge(renal_out, re_renal, by = "relation")
add("renal_effects_recomputed", same(rjoin$spearman_rho, rjoin$rho) && same(rjoin$kendall_tau_b, rjoin$tau), "rho/tau")
add("renal_loo_ranges_recomputed", same(rjoin$leave_one_out_min_rho, rjoin$loo_min) && same(rjoin$leave_one_out_max_rho, rjoin$loo_max) && all(rjoin$leave_one_out_same_sign_n == rjoin$same_sign_n), "range/sign")
add("renal_tables_have_no_inference", !any(grepl("p_value|fdr|q_value|ci_low|ci_high", c(names(renal_out), names(renal_loo_out)), ignore.case = TRUE)), paste(c(names(renal_out), names(renal_loo_out)), collapse = ";"))
add("renal_line_values_match_input", all(renal_line_out$cell_key %chin% renal$cell_key) && nrow(renal_line_out) == 8L, "8/8 cell lines")

gate <- merge(re_primary[, .(relation, primary_rho = rho, primary_fdr = q)], re_adjusted[, .(relation, adjusted_rho = effect, adjusted_fdr = q)], by = "relation")
gate[, expected := sign(primary_rho) == sign(adjusted_rho) & primary_fdr < 0.05 & adjusted_fdr < 0.05]
djoin <- merge(decisions, gate, by = "relation")
add("decision_gate_recomputed", all(djoin$support_gate_pass == djoin$expected), paste(djoin$relation, djoin$expected, collapse = ";"))
add("figure_png_present", file.info(file.path(analysis_dir, "figures_v1", "NCI60_TYMP_MTAP_MTA_TWO_LEVEL_v1.png"))$size > 100000, "PNG size gate")
add("figure_pdf_present", file.info(file.path(analysis_dir, "figures_v1", "NCI60_TYMP_MTAP_MTA_TWO_LEVEL_v1.pdf"))$size > 10000, "PDF size gate")
add("summary_boundary_present", grepl("do not establish a TYMP perturbation effect", paste(readLines(file.path(analysis_dir, "PHASE103_ANALYSIS_SUMMARY_v1.md")), collapse = " "), fixed = TRUE), "noncausal boundary")

fwrite(checks, file.path(qc_dir, "independent_validation_checks_v1.tsv"), sep = "\t", quote = FALSE)
status <- list(
  status = if (all(checks$pass)) "PASS_INDEPENDENT_VALIDATION" else "FAIL_INDEPENDENT_VALIDATION",
  validated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks_passed = sum(checks$pass), checks_total = nrow(checks),
  failed_checks = checks[pass == FALSE, check_id]
)
write_json(status, file.path(qc_dir, "independent_validation_status_v1.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(checks$pass)) stop("Validation failed: ", paste(status$failed_checks, collapse = ", "))
cat("PASS_INDEPENDENT_VALIDATION\n")
cat("checks=", status$checks_passed, "/", status$checks_total, "\n", sep = "")
