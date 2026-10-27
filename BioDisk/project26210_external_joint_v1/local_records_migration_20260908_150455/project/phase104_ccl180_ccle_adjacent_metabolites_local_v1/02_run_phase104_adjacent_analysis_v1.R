#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 02_run_phase104_adjacent_analysis_v1.R <ccl180_all> <ccle_all> <shared> <output_dir>")
ccl <- fread(normalizePath(args[[1]], mustWork = TRUE))
ccle <- fread(normalizePath(args[[2]], mustWork = TRUE))
shared <- fread(normalizePath(args[[3]], mustWork = TRUE))
output_dir <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
results_dir <- file.path(output_dir, "results_v1")
figures_dir <- file.path(output_dir, "figures_v1")
qc_dir <- file.path(output_dir, "qc_v1")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(ccl) != 155L || nrow(ccle) != 148L || nrow(shared) != 148L) stop("Frozen Phase 104 sample count mismatch")

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

bootstrap_rho <- function(x, y, seed, b = 10000L) {
  set.seed(seed)
  n <- length(x)
  values <- replicate(b, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  values <- values[is.finite(values)]
  unname(quantile(values, c(0.025, 0.975), type = 7, names = FALSE))
}

partial_rank <- function(x, y, tissue) {
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  xres <- xr - ave(xr, tissue, FUN = mean)
  yres <- yr - ave(yr, tissue, FUN = mean)
  suppressWarnings(cor(xres, yres))
}

stratified_bootstrap_partial <- function(x, y, tissue, seed, b = 10000L) {
  set.seed(seed)
  groups <- split(seq_along(x), tissue)
  values <- replicate(b, {
    idx <- unlist(lapply(groups, function(v) sample(v, length(v), replace = TRUE)), use.names = FALSE)
    partial_rank(x[idx], y[idx], tissue[idx])
  })
  values <- values[is.finite(values)]
  unname(quantile(values, c(0.025, 0.975), type = 7, names = FALSE))
}

within_tissue_p <- function(x, y, tissue, seed, b = 20000L) {
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  xres <- xr - ave(xr, tissue, FUN = mean)
  yres <- yr - ave(yr, tissue, FUN = mean)
  observed <- suppressWarnings(cor(xres, yres))
  groups <- split(seq_along(x), tissue)
  set.seed(seed)
  extreme <- 0L
  for (i in seq_len(b)) {
    idx <- seq_along(x)
    for (group in groups) idx[group] <- sample(group)
    candidate <- suppressWarnings(cor(xres[idx], yres))
    if (is.finite(candidate) && abs(candidate) >= abs(observed)) extreme <- extreme + 1L
  }
  (extreme + 1) / (b + 1)
}

primary <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  row <- spec[i]
  dat <- data_for(row$dataset)
  x <- dat[[row$gene_col]]
  y <- dat[[row$metabolite_col]]
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  ci <- bootstrap_rho(x, y, 104100L + row$relation_index)
  data.table(
    dataset = row$dataset, relation = row$relation, gene = row$gene, metabolite = row$metabolite,
    n_cell_lines = length(x), spearman_rho = unname(test$estimate),
    bootstrap_ci_low = ci[1], bootstrap_ci_high = ci[2], p_value = test$p.value,
    bootstrap_resamples = 10000L
  )
}))
primary[, fdr_bh := p.adjust(p_value, method = "BH"), by = dataset]

adjusted <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  row <- spec[i]
  dat <- data_for(row$dataset)
  x <- dat[[row$gene_col]]
  y <- dat[[row$metabolite_col]]
  ci <- stratified_bootstrap_partial(x, y, dat$tissue, 104200L + row$relation_index)
  data.table(
    dataset = row$dataset, relation = row$relation, gene = row$gene, metabolite = row$metabolite,
    n_cell_lines = length(x), tissue_levels = uniqueN(dat$tissue),
    tissue_adjusted_partial_spearman = partial_rank(x, y, dat$tissue),
    bootstrap_ci_low = ci[1], bootstrap_ci_high = ci[2],
    permutation_p_value = within_tissue_p(x, y, dat$tissue, 104300L + row$relation_index),
    bootstrap_resamples = 10000L, within_tissue_permutations = 20000L
  )
}))
adjusted[, fdr_bh := p.adjust(permutation_p_value, method = "BH"), by = dataset]

decisions <- merge(
  primary[, .(dataset, relation, gene, metabolite, primary_rho = spearman_rho, primary_fdr = fdr_bh)],
  adjusted[, .(dataset, relation, adjusted_rho = tissue_adjusted_partial_spearman, adjusted_fdr = fdr_bh)],
  by = c("dataset", "relation")
)
decisions[, same_direction := sign(primary_rho) == sign(adjusted_rho)]
decisions[, support_gate_pass := same_direction & primary_fdr < 0.05 & adjusted_fdr < 0.05]
decisions[, decision := ifelse(support_gate_pass,
  "SUPPORTED_ADJACENT_METABOLITE_ASSOCIATION",
  "CURRENT_DATA_DO_NOT_SUPPORT_ROBUST_ADJACENT_ASSOCIATION")]
decisions[, interpretation_boundary := "adjacent metabolite; not MTA, perturbation, mediation, or causality"]

renal_effects <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  row <- spec[i]
  dat <- data_for(row$dataset)[is_renal == TRUE]
  x <- dat[[row$gene_col]]
  y <- dat[[row$metabolite_col]]
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  tau <- suppressWarnings(cor(x, y, method = "kendall"))
  loo <- vapply(seq_along(x), function(k) suppressWarnings(cor(x[-k], y[-k], method = "spearman")), numeric(1))
  stable <- is.finite(rho) && sign(rho) != 0 && all(sign(loo[is.finite(loo)]) == sign(rho)) && min(loo, na.rm = TRUE) * max(loo, na.rm = TRUE) > 0
  data.table(
    dataset = row$dataset, relation = row$relation, gene = row$gene, metabolite = row$metabolite,
    n_cell_lines = nrow(dat), spearman_rho = rho, kendall_tau_b = tau,
    leave_one_out_min_rho = min(loo, na.rm = TRUE), leave_one_out_max_rho = max(loo, na.rm = TRUE),
    leave_one_out_same_sign_n = sum(sign(loo) == sign(rho), na.rm = TRUE),
    leave_one_out_total = sum(is.finite(loo)),
    descriptive_stability = ifelse(stable, "DIRECTIONALLY_STABLE_EXPLORATORY", "DIRECTIONALLY_UNSTABLE_EXPLORATORY"),
    inferential_test_performed = FALSE
  )
}))

renal_lines <- rbindlist(list(
  ccl[is_renal == TRUE, .(dataset = "CCL180", cell_key, cell_line, tissue, tymp_expression, mtap_expression,
                          methionine, sah, polyamine_adjacent = n_acetylputrescine,
                          polyamine_adjacent_identity = "N-Acetylputrescine")],
  ccle[is_renal == TRUE, .(dataset = "CCLE", cell_key, cell_line, tissue, tymp_expression, mtap_expression,
                           methionine, sah, polyamine_adjacent = putrescine,
                           polyamine_adjacent_identity = "Putrescine")]
), fill = TRUE)

platform <- rbindlist(lapply(c("methionine", "sah"), function(metabolite) {
  x <- shared[[paste0(metabolite, "_ccl180")]]
  y <- shared[[paste0(metabolite, "_ccle")]]
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.table(
    metabolite = ifelse(metabolite == "sah", "SAH", "Methionine"),
    n_shared_cell_lines = length(x), spearman_rho = unname(test$estimate), p_value = test$p.value
  )
}))
platform[, fdr_bh := p.adjust(p_value, method = "BH")]
platform[, interpretation_boundary := "cross-platform consistency in overlapping cell lines; not independent replication"]

fwrite(primary, file.path(results_dir, "all_lineages_primary_adjacent_associations_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(adjusted, file.path(results_dir, "all_lineages_tissue_adjusted_adjacent_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(decisions, file.path(results_dir, "all_lineages_adjacent_decisions_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(renal_effects, file.path(results_dir, "renal_adjacent_exploratory_effects_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(renal_lines, file.path(results_dir, "renal_adjacent_line_level_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(platform, file.path(results_dir, "shared_platform_concordance_v1.tsv"), sep = "\t", quote = FALSE)

heat_primary <- primary[, .(dataset, analysis = "Primary", gene, metabolite, rho = spearman_rho, q = fdr_bh)]
heat_adjusted <- adjusted[, .(dataset, analysis = "Tissue-adjusted", gene, metabolite, rho = tissue_adjusted_partial_spearman, q = fdr_bh)]
heat <- rbind(heat_primary, heat_adjusted)
heat[, metabolite_display := factor(metabolite, levels = c("Methionine", "SAH", "N-Acetylputrescine", "Putrescine"))]
heat[, gene := factor(gene, levels = c("MTAP", "TYMP"))]

figure <- ggplot(heat, aes(x = metabolite_display, y = gene, fill = rho)) +
  geom_tile(aes(colour = q < 0.05), linewidth = 1.2) +
  geom_text(aes(label = sprintf("rho=%.2f\nq=%.3g", rho, q)), size = 3.1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "black"), guide = "none") +
  facet_grid(dataset ~ analysis, scales = "free_x", space = "free_x") +
  labs(
    title = "CCL180 and CCLE adjacent-metabolite associations",
    subtitle = "Black border: within-dataset BH q < 0.05; panels do not test MTA",
    x = NULL, y = NULL, fill = "Spearman effect",
    caption = "Independent unit is the cell line. Adjacent metabolites are not MTA substitutes; renal analyses are reported separately without inference."
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom", panel.grid = element_blank())
ggsave(file.path(figures_dir, "CCL180_CCLE_ADJACENT_ASSOCIATIONS_v1.png"), figure, width = 11, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "CCL180_CCLE_ADJACENT_ASSOCIATIONS_v1.pdf"), figure, width = 11, height = 7.5, device = cairo_pdf, bg = "white")

supported <- decisions[support_gate_pass == TRUE]
renal_stable <- renal_effects[descriptive_stability == "DIRECTIONALLY_STABLE_EXPLORATORY"]
summary_lines <- c(
  "# Phase 104 CCL180 and CCLE analysis summary",
  "",
  "## Direct MTA feasibility",
  "",
  "- Direct TYMP-MTA and MTAP-MTA analyses are not evaluable in CCL180 or CCLE because neither dataset contains an accepted MTA feature.",
  "- CCL180 near-mass ions retain conflicting formulas and identities and were not relabelled as MTA.",
  "",
  "## Adjacent-metabolite analyses",
  "",
  paste0("- CCL180 all-lineage n=155; CCLE paired all-lineage n=148. Supported dataset-specific gates: ", nrow(supported), "/12."),
  if (nrow(supported)) paste0("- Supported relations: ", paste(paste(supported$dataset, supported$relation, sep = ":"), collapse = "; "), ".") else "- No relation passed both the primary and tissue-adjusted BH-FDR gates.",
  paste0("- Cross-platform shared-cell-line concordance: Methionine rho=", sprintf("%.4f", platform[metabolite == "Methionine", spearman_rho]),
         ", BH q=", sprintf("%.4g", platform[metabolite == "Methionine", fdr_bh]),
         "; SAH rho=", sprintf("%.4f", platform[metabolite == "SAH", spearman_rho]),
         ", BH q=", sprintf("%.4g", platform[metabolite == "SAH", fdr_bh]), "."),
  "",
  "## Renal-restricted exploratory layer",
  "",
  paste0("- CCL180 renal n=", ccl[is_renal == TRUE, .N], "; CCLE renal n=", ccle[is_renal == TRUE, .N], ". No renal P value, FDR or confidence interval was calculated."),
  paste0("- Directionally stable renal relations: ", nrow(renal_stable), "/12", if (nrow(renal_stable)) paste0(" (", paste(paste(renal_stable$dataset, renal_stable$relation, sep = ":"), collapse = "; "), ")") else "."),
  "",
  "## Evidence boundary",
  "",
  "These results describe explicitly named adjacent metabolites. They do not test MTA and do not establish TYMP-driven MTA production, perturbation, mediation, or causality."
)
writeLines(summary_lines, file.path(output_dir, "PHASE104_ANALYSIS_SUMMARY_v1.md"), useBytes = TRUE)

checks <- data.table(
  check_id = c("ccl180_n_155", "ccle_n_148", "primary_twelve", "adjusted_twelve", "decisions_twelve",
               "renal_ccl180_below_10", "renal_ccle_below_10", "renal_no_p_fields", "platform_two",
               "direct_mta_not_claimed", "png_exists", "pdf_exists"),
  pass = c(nrow(ccl) == 155L, nrow(ccle) == 148L, nrow(primary) == 12L, nrow(adjusted) == 12L, nrow(decisions) == 12L,
           ccl[is_renal == TRUE, .N] < 10L, ccle[is_renal == TRUE, .N] < 10L,
           !any(grepl("p_value|fdr|q_value|ci_low|ci_high", names(renal_effects), ignore.case = TRUE)),
           nrow(platform) == 2L,
           all(decisions$interpretation_boundary == "adjacent metabolite; not MTA, perturbation, mediation, or causality"),
           file.exists(file.path(figures_dir, "CCL180_CCLE_ADJACENT_ASSOCIATIONS_v1.png")),
           file.exists(file.path(figures_dir, "CCL180_CCLE_ADJACENT_ASSOCIATIONS_v1.pdf")))
)
checks[, detail := c(nrow(ccl), nrow(ccle), nrow(primary), nrow(adjusted), nrow(decisions),
                     ccl[is_renal == TRUE, .N], ccle[is_renal == TRUE, .N], paste(names(renal_effects), collapse = ";"),
                     nrow(platform), "boundary", "PNG", "PDF")]
fwrite(checks, file.path(qc_dir, "analysis_checks_v1.tsv"), sep = "\t", quote = FALSE)
status <- list(
  status = if (all(checks$pass)) "PASS_ANALYSIS_COMPLETE" else "FAIL_ANALYSIS_CHECKS",
  direct_mta_evaluable = FALSE,
  ccl180_n = nrow(ccl), ccle_n = nrow(ccle),
  ccl180_renal_n = ccl[is_renal == TRUE, .N], ccle_renal_n = ccle[is_renal == TRUE, .N],
  supported_adjacent_relations = nrow(supported),
  stable_renal_relations = nrow(renal_stable),
  checks_passed = sum(checks$pass), checks_total = nrow(checks),
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(status, file.path(qc_dir, "analysis_status_v1.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(checks$pass)) stop("Analysis checks failed: ", paste(checks[pass == FALSE, check_id], collapse = ", "))
cat("PASS_ANALYSIS_COMPLETE\n")
cat("supported=", nrow(supported), "/12 stable_renal=", nrow(renal_stable), "/12\n", sep = "")
