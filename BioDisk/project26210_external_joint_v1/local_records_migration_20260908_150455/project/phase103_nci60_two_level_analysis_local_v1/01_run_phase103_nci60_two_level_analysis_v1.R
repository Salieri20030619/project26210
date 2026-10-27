#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 01_run_phase103_nci60_two_level_analysis_v1.R <all_table> <renal_table> <output_dir>")
all_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
renal_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
results_dir <- file.path(output_dir, "results_v1")
figures_dir <- file.path(output_dir, "figures_v1")
qc_dir <- file.path(output_dir, "qc_v1")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

all_dt <- fread(all_path)
renal <- fread(renal_path)
required <- c("cell_key", "metabolomics_cell_line", "tissue", "tymp_log2_fpkm_plus1", "mtap_log2_fpkm_plus1", "mta_mean_raw", "mta_all_repeats_at_global_floor")
if (!all(required %in% names(all_dt)) || !all(required %in% names(renal))) stop("Required Phase 102 columns are missing")
if (nrow(all_dt) != 58L || nrow(renal) != 8L) stop("Frozen sample count mismatch")
if (anyNA(all_dt[, ..required]) || anyNA(renal[, ..required])) stop("Missing value in frozen analysis variables")
if (anyDuplicated(all_dt$cell_key) || anyDuplicated(renal$cell_key)) stop("Duplicate independent cell-line key")

relations <- data.table(
  relation = c("TYMP-MTA", "MTAP-MTA"),
  gene = c("TYMP", "MTAP"),
  x_col = c("tymp_log2_fpkm_plus1", "mtap_log2_fpkm_plus1"),
  relation_order = 1:2
)
floor_value <- min(all_dt$mta_mean_raw)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

bootstrap_spearman <- function(x, y, b, seed) {
  set.seed(seed)
  n <- length(x)
  vals <- replicate(b, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  vals <- vals[is.finite(vals)]
  unname(quantile(vals, c(0.025, 0.975), names = FALSE, type = 7))
}

partial_rank_effect <- function(x, y, tissue) {
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  xres <- xr - ave(xr, tissue, FUN = mean)
  yres <- yr - ave(yr, tissue, FUN = mean)
  suppressWarnings(cor(xres, yres, method = "pearson"))
}

stratified_bootstrap_partial <- function(x, y, tissue, b, seed) {
  set.seed(seed)
  groups <- split(seq_along(x), tissue)
  vals <- replicate(b, {
    idx <- unlist(lapply(groups, function(z) sample(z, length(z), replace = TRUE)), use.names = FALSE)
    partial_rank_effect(x[idx], y[idx], tissue[idx])
  })
  vals <- vals[is.finite(vals)]
  unname(quantile(vals, c(0.025, 0.975), names = FALSE, type = 7))
}

within_tissue_permutation_p <- function(x, y, tissue, b, seed) {
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  xres <- xr - ave(xr, tissue, FUN = mean)
  yres <- yr - ave(yr, tissue, FUN = mean)
  observed <- suppressWarnings(cor(xres, yres))
  groups <- split(seq_along(x), tissue)
  set.seed(seed)
  exceed <- 0L
  for (i in seq_len(b)) {
    perm_idx <- seq_along(x)
    for (g in groups) perm_idx[g] <- sample(g, length(g), replace = FALSE)
    permuted <- suppressWarnings(cor(xres[perm_idx], yres))
    if (is.finite(permuted) && abs(permuted) >= abs(observed)) exceed <- exceed + 1L
  }
  (exceed + 1) / (b + 1)
}

primary_list <- vector("list", nrow(relations))
adjusted_list <- vector("list", nrow(relations))
kendall_list <- vector("list", nrow(relations))
above_floor_list <- vector("list", nrow(relations))

for (i in seq_len(nrow(relations))) {
  rel <- relations[i]
  x <- all_dt[[rel$x_col]]
  y <- all_dt$mta_mean_raw
  sp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE, alternative = "two.sided"))
  sp_ci <- bootstrap_spearman(x, y, 20000L, 103100L + i)
  primary_list[[i]] <- data.table(
    relation = rel$relation, gene = rel$gene, n_cell_lines = length(x),
    spearman_rho = unname(sp$estimate), bootstrap_ci_low = sp_ci[1], bootstrap_ci_high = sp_ci[2],
    p_value = sp$p.value, bootstrap_resamples = 20000L, mta_floor = floor_value,
    mta_lines_at_floor = sum(all_dt$mta_mean_raw == floor_value)
  )

  part <- partial_rank_effect(x, y, all_dt$tissue)
  part_ci <- stratified_bootstrap_partial(x, y, all_dt$tissue, 20000L, 103200L + i)
  perm_p <- within_tissue_permutation_p(x, y, all_dt$tissue, 100000L, 103300L + i)
  adjusted_list[[i]] <- data.table(
    relation = rel$relation, gene = rel$gene, n_cell_lines = length(x), tissue_levels = uniqueN(all_dt$tissue),
    tissue_adjusted_partial_spearman = part, stratified_bootstrap_ci_low = part_ci[1],
    stratified_bootstrap_ci_high = part_ci[2], permutation_p_value = perm_p,
    bootstrap_resamples = 20000L, within_tissue_permutations = 100000L
  )

  kt <- suppressWarnings(cor.test(x, y, method = "kendall", exact = FALSE, alternative = "two.sided"))
  kendall_list[[i]] <- data.table(
    relation = rel$relation, gene = rel$gene, n_cell_lines = length(x),
    kendall_tau_b = unname(kt$estimate), p_value = kt$p.value
  )

  af <- all_dt[mta_mean_raw > floor_value]
  af_test <- suppressWarnings(cor.test(af[[rel$x_col]], af$mta_mean_raw, method = "spearman", exact = FALSE))
  above_floor_list[[i]] <- data.table(
    relation = rel$relation, gene = rel$gene, n_cell_lines_above_floor = nrow(af),
    spearman_rho = unname(af_test$estimate), p_value = af_test$p.value,
    sensitivity_only = TRUE
  )
}

primary <- rbindlist(primary_list)
primary[, fdr_bh := p.adjust(p_value, method = "BH")]
adjusted <- rbindlist(adjusted_list)
adjusted[, fdr_bh := p.adjust(permutation_p_value, method = "BH")]
kendall <- rbindlist(kendall_list)
kendall[, fdr_bh := p.adjust(p_value, method = "BH")]
above_floor <- rbindlist(above_floor_list)
above_floor[, fdr_bh := p.adjust(p_value, method = "BH")]

all_model <- copy(all_dt)
all_model[, y_rank_z := as.numeric(scale(rank(mta_mean_raw, ties.method = "average")))]
all_model[, tymp_rank_z := as.numeric(scale(rank(tymp_log2_fpkm_plus1, ties.method = "average")))]
all_model[, mtap_rank_z := as.numeric(scale(rank(mtap_log2_fpkm_plus1, ties.method = "average")))]
all_model[, tissue := factor(tissue)]
fit <- lm(y_rank_z ~ tymp_rank_z + mtap_rank_z + tissue, data = all_model)
fit_sum <- summary(fit)$coefficients
fit_ci <- confint(fit)
mutual <- data.table(
  relation = c("TYMP-MTA", "MTAP-MTA"),
  gene = c("TYMP", "MTAP"),
  n_cell_lines = nrow(all_model),
  standardized_rank_beta = c(fit_sum["tymp_rank_z", "Estimate"], fit_sum["mtap_rank_z", "Estimate"]),
  ci_low = c(fit_ci["tymp_rank_z", 1], fit_ci["mtap_rank_z", 1]),
  ci_high = c(fit_ci["tymp_rank_z", 2], fit_ci["mtap_rank_z", 2]),
  p_value = c(fit_sum["tymp_rank_z", "Pr(>|t|)"], fit_sum["mtap_rank_z", "Pr(>|t|)"]),
  residual_df = df.residual(fit),
  covariates = "other_gene_rank+tissue_fixed_effects"
)
mutual[, fdr_bh := p.adjust(p_value, method = "BH")]

tissue_effects <- rbindlist(lapply(seq_len(nrow(relations)), function(i) {
  rel <- relations[i]
  rbindlist(lapply(split(all_dt, all_dt$tissue), function(z) {
    if (nrow(z) < 3L) return(NULL)
    data.table(
      relation = rel$relation, gene = rel$gene, tissue = unique(z$tissue), n_cell_lines = nrow(z),
      spearman_rho = suppressWarnings(cor(z[[rel$x_col]], z$mta_mean_raw, method = "spearman")),
      lines_at_mta_floor = sum(z$mta_mean_raw == floor_value), inferential_test_performed = FALSE
    )
  }))
}))
setorder(tissue_effects, relation, tissue)

decisions <- merge(primary[, .(relation, primary_rho = spearman_rho, primary_fdr = fdr_bh)],
                   adjusted[, .(relation, adjusted_rho = tissue_adjusted_partial_spearman, adjusted_fdr = fdr_bh)],
                   by = "relation")
decisions[, same_direction := sign(primary_rho) == sign(adjusted_rho)]
decisions[, support_gate_pass := same_direction & primary_fdr < 0.05 & adjusted_fdr < 0.05]
decisions[, decision := ifelse(support_gate_pass,
  "SUPPORTED_BASELINE_CROSS_LINEAGE_ASSOCIATION",
  "CURRENT_DATA_DO_NOT_SUPPORT_ROBUST_CROSS_LINEAGE_ASSOCIATION")]
decisions[, interpretation_boundary := "baseline association; not perturbation, mediation, or causality"]

renal_effects <- rbindlist(lapply(seq_len(nrow(relations)), function(i) {
  rel <- relations[i]
  x <- renal[[rel$x_col]]
  y <- renal$mta_mean_raw
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  tau <- suppressWarnings(cor(x, y, method = "kendall"))
  loo <- vapply(seq_along(x), function(k) suppressWarnings(cor(x[-k], y[-k], method = "spearman")), numeric(1))
  stable <- sign(rho) != 0 && all(sign(loo[is.finite(loo)]) == sign(rho)) && min(loo, na.rm = TRUE) * max(loo, na.rm = TRUE) > 0
  data.table(
    relation = rel$relation, gene = rel$gene, n_cell_lines = length(x), spearman_rho = rho,
    kendall_tau_b = tau, leave_one_out_min_rho = min(loo, na.rm = TRUE),
    leave_one_out_max_rho = max(loo, na.rm = TRUE),
    leave_one_out_same_sign_n = sum(sign(loo) == sign(rho), na.rm = TRUE),
    leave_one_out_total = sum(is.finite(loo)),
    descriptive_stability = ifelse(stable, "DIRECTIONALLY_STABLE_EXPLORATORY", "DIRECTIONALLY_UNSTABLE_EXPLORATORY"),
    inferential_test_performed = FALSE
  )
}))

renal_loo <- rbindlist(lapply(seq_len(nrow(relations)), function(i) {
  rel <- relations[i]
  x <- renal[[rel$x_col]]
  y <- renal$mta_mean_raw
  data.table(
    relation = rel$relation,
    omitted_cell_line = renal$metabolomics_cell_line,
    n_remaining = length(x) - 1L,
    spearman_rho = vapply(seq_along(x), function(k) suppressWarnings(cor(x[-k], y[-k], method = "spearman")), numeric(1)),
    inferential_test_performed = FALSE
  )
}))

renal_line <- copy(renal)
renal_line[, mta_rank_average := rank(mta_mean_raw, ties.method = "average")]
renal_line[, tymp_rank_average := rank(tymp_log2_fpkm_plus1, ties.method = "average")]
renal_line[, mtap_rank_average := rank(mtap_log2_fpkm_plus1, ties.method = "average")]
renal_line[, mta_at_global_floor := mta_mean_raw == floor_value]
setorder(renal_line, mta_rank_average, metabolomics_cell_line)

renal_floor <- rbindlist(lapply(seq_len(nrow(relations)), function(i) {
  rel <- relations[i]
  x <- renal[[rel$x_col]]
  at_floor <- renal$mta_mean_raw == floor_value
  data.table(
    relation = rel$relation, gene = rel$gene,
    n_floor = sum(at_floor), n_above_floor = sum(!at_floor),
    median_expression_floor = median(x[at_floor]),
    median_expression_above_floor = median(x[!at_floor]),
    median_difference_above_minus_floor = median(x[!at_floor]) - median(x[at_floor]),
    inferential_test_performed = FALSE
  )
}))

fwrite(primary, file.path(results_dir, "all_lineages_primary_associations_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(adjusted, file.path(results_dir, "all_lineages_tissue_adjusted_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(kendall, file.path(results_dir, "all_lineages_kendall_floor_sensitivity_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(mutual, file.path(results_dir, "all_lineages_mutually_adjusted_rank_model_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(above_floor, file.path(results_dir, "all_lineages_above_floor_sensitivity_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(tissue_effects, file.path(results_dir, "tissue_specific_effects_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(decisions, file.path(results_dir, "all_lineages_decisions_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(renal_effects, file.path(results_dir, "renal_exploratory_effects_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(renal_loo, file.path(results_dir, "renal_leave_one_out_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(renal_floor, file.path(results_dir, "renal_floor_group_descriptives_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(renal_line, file.path(results_dir, "renal_line_level_v1.tsv"), sep = "\t", quote = FALSE)

make_all_panel <- function(x_col, gene_label, rho_value) {
  ggplot(all_dt, aes(x = .data[[x_col]], y = mta_mean_raw, colour = tissue,
                     shape = mta_mean_raw == floor_value)) +
    geom_hline(yintercept = floor_value, colour = "grey55", linetype = "dotted", linewidth = 0.5) +
    geom_point(size = 2.6, alpha = 0.88) +
    scale_y_continuous(trans = "log2") +
    scale_colour_viridis_d(end = 0.9) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), name = "At MTA floor") +
    labs(x = paste0(gene_label, " log2(FPKM+1)"), y = "MTA mean (log2 scale)",
         title = paste0("All NCI60 lineages: ", gene_label, "–MTA"),
         subtitle = paste0("n=58; Spearman rho=", fmt_num(rho_value, 3))) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
}

make_renal_panel <- function(x_col, gene_label, rho_value) {
  ggplot(renal, aes(x = .data[[x_col]], y = mta_mean_raw, label = metabolomics_cell_line,
                   shape = mta_mean_raw == floor_value)) +
    geom_hline(yintercept = floor_value, colour = "grey55", linetype = "dotted", linewidth = 0.5) +
    geom_point(size = 3.0, colour = "#7A3E9D") +
    geom_text_repel(size = 3.0, min.segment.length = 0, max.overlaps = Inf, seed = 103400) +
    scale_y_continuous(trans = "log2") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), name = "At MTA floor") +
    labs(x = paste0(gene_label, " log2(FPKM+1)"), y = "MTA mean (log2 scale)",
         title = paste0("Renal only: ", gene_label, "–MTA"),
         subtitle = paste0("n=8 exploratory; Spearman rho=", fmt_num(rho_value, 3), "; no inferential test")) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
}

p1 <- make_all_panel("tymp_log2_fpkm_plus1", "TYMP", primary[relation == "TYMP-MTA", spearman_rho])
p2 <- make_all_panel("mtap_log2_fpkm_plus1", "MTAP", primary[relation == "MTAP-MTA", spearman_rho])
p3 <- make_renal_panel("tymp_log2_fpkm_plus1", "TYMP", renal_effects[relation == "TYMP-MTA", spearman_rho])
p4 <- make_renal_panel("mtap_log2_fpkm_plus1", "MTAP", renal_effects[relation == "MTAP-MTA", spearman_rho])
p1 <- p1 + labs(title = 'All NCI60 lineages: TYMP-MTA')
p2 <- p2 + labs(title = 'All NCI60 lineages: MTAP-MTA')
p3 <- p3 + labs(title = 'Renal only: TYMP-MTA')
p4 <- p4 + labs(title = 'Renal only: MTAP-MTA')
p1 <- p1 + guides(shape = 'none')
p2 <- p2 + guides(shape = 'none')
p3 <- p3 + guides(shape = 'none') + labs(subtitle = paste0('n=8; rho=', fmt_num(renal_effects[relation == 'TYMP-MTA', spearman_rho], 3), '; descriptive only'))
p4 <- p4 + guides(shape = 'none') + labs(subtitle = paste0('n=8; rho=', fmt_num(renal_effects[relation == 'MTAP-MTA', spearman_rho], 3), '; descriptive only'))
combined <- ((p1 | p2) / (p3 | p4)) + plot_layout(guides = 'collect') + plot_annotation(
  title = "NCI60 baseline TYMP/MTAP–MTA analysis",
  subtitle = "Points are independent cell lines; open points lie at the observed MTA floor of 0.086",
  caption = "Baseline association only. Technical repeats were aggregated per cell line; renal n=8 receives no inferential test."
)
combined <- combined + plot_annotation(
  title = 'NCI60 baseline TYMP/MTAP-MTA analysis',
  subtitle = 'Points are independent cell lines; open points lie at the observed MTA floor of 0.086',
  caption = 'Baseline association only. Technical repeats were aggregated per cell line; renal n=8 receives no inferential test.'
)
combined <- combined & theme(legend.position = 'bottom')
ggsave(file.path(figures_dir, "NCI60_TYMP_MTAP_MTA_TWO_LEVEL_v1.png"), combined, width = 12, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "NCI60_TYMP_MTAP_MTA_TWO_LEVEL_v1.pdf"), combined, width = 12, height = 8.5, device = cairo_pdf, bg = "white")

summary_lines <- c(
  "# Phase 103 NCI60 two-level TYMP/MTAP–MTA analysis",
  "",
  "## All NCI60 lineages",
  "",
  paste0("- TYMP–MTA: Spearman rho=", fmt_num(primary[relation == "TYMP-MTA", spearman_rho]),
         ", 95% bootstrap CI [", fmt_num(primary[relation == "TYMP-MTA", bootstrap_ci_low]), ", ",
         fmt_num(primary[relation == "TYMP-MTA", bootstrap_ci_high]), "], P=", fmt_num(primary[relation == "TYMP-MTA", p_value]),
         ", BH q=", fmt_num(primary[relation == "TYMP-MTA", fdr_bh]), ". Tissue-adjusted partial rho=",
         fmt_num(adjusted[relation == "TYMP-MTA", tissue_adjusted_partial_spearman]), ", permutation P=",
         fmt_num(adjusted[relation == "TYMP-MTA", permutation_p_value]), ", BH q=",
         fmt_num(adjusted[relation == "TYMP-MTA", fdr_bh]), ". Decision: ", decisions[relation == "TYMP-MTA", decision], "."),
  paste0("- MTAP–MTA: Spearman rho=", fmt_num(primary[relation == "MTAP-MTA", spearman_rho]),
         ", 95% bootstrap CI [", fmt_num(primary[relation == "MTAP-MTA", bootstrap_ci_low]), ", ",
         fmt_num(primary[relation == "MTAP-MTA", bootstrap_ci_high]), "], P=", fmt_num(primary[relation == "MTAP-MTA", p_value]),
         ", BH q=", fmt_num(primary[relation == "MTAP-MTA", fdr_bh]), ". Tissue-adjusted partial rho=",
         fmt_num(adjusted[relation == "MTAP-MTA", tissue_adjusted_partial_spearman]), ", permutation P=",
         fmt_num(adjusted[relation == "MTAP-MTA", permutation_p_value]), ", BH q=",
         fmt_num(adjusted[relation == "MTAP-MTA", fdr_bh]), ". Decision: ", decisions[relation == "MTAP-MTA", decision], "."),
  paste0("- MTA floor: ", floor_value, "; ", sum(all_dt$mta_mean_raw == floor_value), "/58 cell-line means are at the floor."),
  "",
  "## Renal-restricted exploratory layer",
  "",
  paste0("- TYMP–MTA: Spearman rho=", fmt_num(renal_effects[relation == "TYMP-MTA", spearman_rho]),
         "; Kendall tau-b=", fmt_num(renal_effects[relation == "TYMP-MTA", kendall_tau_b]),
         "; leave-one-out rho range [", fmt_num(renal_effects[relation == "TYMP-MTA", leave_one_out_min_rho]), ", ",
         fmt_num(renal_effects[relation == "TYMP-MTA", leave_one_out_max_rho]), "]; ",
         renal_effects[relation == "TYMP-MTA", descriptive_stability], "."),
  paste0("- MTAP–MTA: Spearman rho=", fmt_num(renal_effects[relation == "MTAP-MTA", spearman_rho]),
         "; Kendall tau-b=", fmt_num(renal_effects[relation == "MTAP-MTA", kendall_tau_b]),
         "; leave-one-out rho range [", fmt_num(renal_effects[relation == "MTAP-MTA", leave_one_out_min_rho]), ", ",
         fmt_num(renal_effects[relation == "MTAP-MTA", leave_one_out_max_rho]), "]; ",
         renal_effects[relation == "MTAP-MTA", descriptive_stability], "."),
  paste0("- Four of eight renal cell-line means are at the observed MTA floor; no renal P value, FDR or confidence interval was calculated."),
  "",
  "## Evidence boundary",
  "",
  "These are baseline cross-cell-line associations. They do not establish a TYMP perturbation effect, temporal mediation, or causality. Tissue adjustment reduces but cannot remove all lineage and cell-state confounding."
)
summary_lines <- gsub('\u2013', '-', enc2utf8(summary_lines), fixed = TRUE)
writeLines(summary_lines, file.path(output_dir, "PHASE103_ANALYSIS_SUMMARY_v1.md"), useBytes = TRUE)

checks <- data.table(
  check_id = c("all_n_58", "renal_n_8", "primary_two_relations", "adjusted_two_relations", "renal_no_p_columns",
               "renal_no_fdr_columns", "floor_preserved", "decision_gate_deterministic", "png_exists", "pdf_exists"),
  pass = c(nrow(all_dt) == 58L, nrow(renal) == 8L, nrow(primary) == 2L, nrow(adjusted) == 2L,
           !any(grepl("p_value", names(renal_effects))), !any(grepl("fdr|q_value", names(renal_effects), ignore.case = TRUE)),
           floor_value == 0.086 && sum(all_dt$mta_mean_raw == floor_value) == 17L,
           all(decisions$support_gate_pass == (decisions$same_direction & decisions$primary_fdr < 0.05 & decisions$adjusted_fdr < 0.05)),
           file.exists(file.path(figures_dir, "NCI60_TYMP_MTAP_MTA_TWO_LEVEL_v1.png")),
           file.exists(file.path(figures_dir, "NCI60_TYMP_MTAP_MTA_TWO_LEVEL_v1.pdf")))
)
checks[, detail := c(nrow(all_dt), nrow(renal), nrow(primary), nrow(adjusted), paste(names(renal_effects), collapse = ";"),
                     paste(names(renal_effects), collapse = ";"), paste0("floor=", floor_value, ";n=", sum(all_dt$mta_mean_raw == floor_value)),
                     paste(decisions$decision, collapse = ";"), "PNG", "PDF")]
fwrite(checks, file.path(qc_dir, "analysis_checks_v1.tsv"), sep = "\t", quote = FALSE)
status <- list(
  status = if (all(checks$pass)) "PASS_ANALYSIS_COMPLETE" else "FAIL_ANALYSIS_CHECKS",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  all_lineages_n = nrow(all_dt), renal_n = nrow(renal),
  primary_relations = nrow(primary), checks_passed = sum(checks$pass), checks_total = nrow(checks),
  decisions = setNames(as.list(decisions$decision), decisions$relation)
)
write_json(status, file.path(qc_dir, "analysis_status_v1.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(checks$pass)) stop("Analysis checks failed: ", paste(checks[pass == FALSE, check_id], collapse = ", "))
cat("PASS_ANALYSIS_COMPLETE\n")
cat("all_n=58 renal_n=8 checks=", sum(checks$pass), "/", nrow(checks), "\n", sep = "")
