#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 03_validate_phase104_v1.R <ccl180_all> <ccle_all> <shared> <analysis_dir>")
ccl <- fread(normalizePath(args[[1]], mustWork = TRUE))
ccle <- fread(normalizePath(args[[2]], mustWork = TRUE))
shared <- fread(normalizePath(args[[3]], mustWork = TRUE))
root <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
results_dir <- file.path(root, "results_v1")
qc_dir <- file.path(root, "qc_v1")

primary <- fread(file.path(results_dir, "all_lineages_primary_adjacent_associations_v1.tsv"))
adjusted <- fread(file.path(results_dir, "all_lineages_tissue_adjusted_adjacent_v1.tsv"))
decisions <- fread(file.path(results_dir, "all_lineages_adjacent_decisions_v1.tsv"))
renal <- fread(file.path(results_dir, "renal_adjacent_exploratory_effects_v1.tsv"))
platform <- fread(file.path(results_dir, "shared_platform_concordance_v1.tsv"))
feasibility <- fread(file.path(root, "preprocessed_v1", "direct_mta_feasibility_v1.tsv"))

spec <- rbindlist(list(
  CJ(dataset = "CCL180", gene = c("TYMP", "MTAP"), metabolite = c("Methionine", "SAH", "N-Acetylputrescine")),
  CJ(dataset = "CCLE", gene = c("TYMP", "MTAP"), metabolite = c("Methionine", "SAH", "Putrescine"))
))
spec[, gene_col := ifelse(gene == "TYMP", "tymp_expression", "mtap_expression")]
spec[, metabolite_col := fifelse(metabolite == "Methionine", "methionine",
                           fifelse(metabolite == "SAH", "sah",
                           fifelse(metabolite == "N-Acetylputrescine", "n_acetylputrescine", "putrescine")))]
spec[, relation := paste(gene, metabolite, sep = "-")]
setorder(spec, dataset, gene, metabolite)
spec[, relation_index := seq_len(.N)]
data_for <- function(dataset) if (dataset == "CCL180") ccl else ccle

boot_rho <- function(x, y, seed, b = 10000L) {
  set.seed(seed)
  n <- length(x)
  values <- numeric(b)
  for (i in seq_len(b)) {
    idx <- sample.int(n, n, replace = TRUE)
    values[i] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }
  unname(quantile(values[is.finite(values)], c(0.025, 0.975), type = 7))
}

partial_rank <- function(x, y, tissue) {
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  cor(xr - ave(xr, tissue, FUN = mean), yr - ave(yr, tissue, FUN = mean))
}

boot_partial <- function(x, y, tissue, seed, b = 10000L) {
  set.seed(seed)
  groups <- split(seq_along(x), tissue)
  values <- numeric(b)
  for (i in seq_len(b)) {
    idx <- unlist(lapply(groups, function(v) sample(v, length(v), replace = TRUE)), use.names = FALSE)
    values[i] <- partial_rank(x[idx], y[idx], tissue[idx])
  }
  unname(quantile(values[is.finite(values)], c(0.025, 0.975), type = 7))
}

perm_p <- function(x, y, tissue, seed, b = 20000L) {
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  xres <- xr - ave(xr, tissue, FUN = mean)
  yres <- yr - ave(yr, tissue, FUN = mean)
  observed <- cor(xres, yres)
  groups <- split(seq_along(x), tissue)
  set.seed(seed)
  extreme <- 0L
  for (i in seq_len(b)) {
    idx <- seq_along(x)
    for (group in groups) idx[group] <- sample(group)
    if (abs(cor(xres[idx], yres)) >= abs(observed)) extreme <- extreme + 1L
  }
  (extreme + 1) / (b + 1)
}

re_primary <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  row <- spec[i]
  dat <- data_for(row$dataset)
  x <- dat[[row$gene_col]]
  y <- dat[[row$metabolite_col]]
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  ci <- boot_rho(x, y, 104100L + row$relation_index)
  data.table(dataset = row$dataset, relation = row$relation, rho = unname(test$estimate), p = test$p.value, ci_low = ci[1], ci_high = ci[2])
}))
re_primary[, q := p.adjust(p, "BH"), by = dataset]

re_adjusted <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  row <- spec[i]
  dat <- data_for(row$dataset)
  x <- dat[[row$gene_col]]
  y <- dat[[row$metabolite_col]]
  ci <- boot_partial(x, y, dat$tissue, 104200L + row$relation_index)
  data.table(dataset = row$dataset, relation = row$relation, effect = partial_rank(x, y, dat$tissue),
             ci_low = ci[1], ci_high = ci[2], p = perm_p(x, y, dat$tissue, 104300L + row$relation_index))
}))
re_adjusted[, q := p.adjust(p, "BH"), by = dataset]

re_renal <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  row <- spec[i]
  dat <- data_for(row$dataset)[is_renal == TRUE]
  x <- dat[[row$gene_col]]
  y <- dat[[row$metabolite_col]]
  rho <- cor(x, y, method = "spearman")
  loo <- vapply(seq_along(x), function(k) cor(x[-k], y[-k], method = "spearman"), numeric(1))
  data.table(dataset = row$dataset, relation = row$relation, rho = rho, tau = cor(x, y, method = "kendall"),
             loo_min = min(loo), loo_max = max(loo), same_sign_n = sum(sign(loo) == sign(rho)))
}))

re_platform <- rbindlist(lapply(c("methionine", "sah"), function(metabolite) {
  test <- suppressWarnings(cor.test(shared[[paste0(metabolite, "_ccl180")]], shared[[paste0(metabolite, "_ccle")]], method = "spearman", exact = FALSE))
  data.table(metabolite = ifelse(metabolite == "sah", "SAH", "Methionine"), rho = unname(test$estimate), p = test$p.value)
}))
re_platform[, q := p.adjust(p, "BH")]

tol <- 1e-12
same <- function(a, b) all(abs(a - b) <= tol)
checks <- data.table(check_id = character(), pass = logical(), detail = character())
add <- function(id, pass, detail) checks <<- rbind(checks, data.table(check_id = id, pass = pass, detail = detail))

pj <- merge(primary, re_primary, by = c("dataset", "relation"))
add("primary_rho_recomputed", same(pj$spearman_rho, pj$rho), paste(max(abs(pj$spearman_rho - pj$rho))))
add("primary_p_recomputed", same(pj$p_value, pj$p), paste(max(abs(pj$p_value - pj$p))))
add("primary_q_recomputed", same(pj$fdr_bh, pj$q), paste(max(abs(pj$fdr_bh - pj$q))))
add("primary_ci_recomputed", same(pj$bootstrap_ci_low, pj$ci_low) && same(pj$bootstrap_ci_high, pj$ci_high), "10,000 resamples")

aj <- merge(adjusted, re_adjusted, by = c("dataset", "relation"))
add("adjusted_effect_recomputed", same(aj$tissue_adjusted_partial_spearman, aj$effect), paste(max(abs(aj$tissue_adjusted_partial_spearman - aj$effect))))
add("adjusted_p_recomputed", same(aj$permutation_p_value, aj$p), paste(max(abs(aj$permutation_p_value - aj$p))))
add("adjusted_q_recomputed", same(aj$fdr_bh, aj$q), paste(max(abs(aj$fdr_bh - aj$q))))
add("adjusted_ci_recomputed", same(aj$bootstrap_ci_low, aj$ci_low) && same(aj$bootstrap_ci_high, aj$ci_high), "10,000 stratified resamples")

rj <- merge(renal, re_renal, by = c("dataset", "relation"))
add("renal_effects_recomputed", same(rj$spearman_rho, rj$rho) && same(rj$kendall_tau_b, rj$tau), "rho/tau")
add("renal_loo_recomputed", same(rj$leave_one_out_min_rho, rj$loo_min) && same(rj$leave_one_out_max_rho, rj$loo_max) && all(rj$leave_one_out_same_sign_n == rj$same_sign_n), "range/sign")
add("renal_no_inference_fields", !any(grepl("p_value|fdr|q_value|ci_low|ci_high", names(renal), ignore.case = TRUE)), paste(names(renal), collapse = ";"))

plj <- merge(platform, re_platform, by = "metabolite")
add("platform_recomputed", same(plj$spearman_rho, plj$rho) && same(plj$p_value, plj$p) && same(plj$fdr_bh, plj$q), "rho/P/q")

gate <- merge(re_primary[, .(dataset, relation, primary_rho = rho, primary_fdr = q)],
              re_adjusted[, .(dataset, relation, adjusted_rho = effect, adjusted_fdr = q)], by = c("dataset", "relation"))
gate[, expected := sign(primary_rho) == sign(adjusted_rho) & primary_fdr < 0.05 & adjusted_fdr < 0.05]
dj <- merge(decisions, gate, by = c("dataset", "relation"))
add("decision_gate_recomputed", all(dj$support_gate_pass == dj$expected), paste(sum(dj$expected), "/12"))
add("direct_mta_not_evaluable", nrow(feasibility) == 2L && all(feasibility$decision == "DIRECT_MTA_NOT_EVALUABLE"), paste(feasibility$dataset, collapse = ";"))
add("preprocessing_counts", nrow(ccl) == 155L && nrow(ccle) == 148L && nrow(shared) == 148L, "155/148/148")
add("figure_png_present", file.info(file.path(root, "figures_v1", "CCL180_CCLE_ADJACENT_ASSOCIATIONS_v1.png"))$size > 100000, "PNG size")
add("figure_pdf_present", file.info(file.path(root, "figures_v1", "CCL180_CCLE_ADJACENT_ASSOCIATIONS_v1.pdf"))$size > 10000, "PDF size")
add("summary_boundary_present", grepl("They do not test MTA", paste(readLines(file.path(root, "PHASE104_ANALYSIS_SUMMARY_v1.md")), collapse = " "), fixed = TRUE), "MTA boundary")

fwrite(checks, file.path(qc_dir, "independent_validation_checks_v1.tsv"), sep = "\t", quote = FALSE)
status <- list(
  status = if (all(checks$pass)) "PASS_INDEPENDENT_VALIDATION" else "FAIL_INDEPENDENT_VALIDATION",
  validated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks_passed = sum(checks$pass), checks_total = nrow(checks),
  failed_checks = checks[pass == FALSE, check_id]
)
write_json(status, file.path(qc_dir, "independent_validation_status_v1.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(checks$pass)) stop("Independent validation failed: ", paste(status$failed_checks, collapse = ", "))
cat("PASS_INDEPENDENT_VALIDATION\n")
cat("checks=", status$checks_passed, "/", status$checks_total, "\n", sep = "")
