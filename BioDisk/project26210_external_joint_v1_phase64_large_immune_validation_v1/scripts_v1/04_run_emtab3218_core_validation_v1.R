args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript 04_run_emtab3218_core_validation_v1.R <rma_matrix> <annotation> <sample_manifest> <patient_clinical> <signature_genes> <output_dir>")
}

rma_path <- args[[1L]]
annotation_path <- args[[2L]]
sample_manifest_path <- args[[3L]]
patient_clinical_path <- args[[4L]]
signature_path <- args[[5L]]
output_dir <- args[[6L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(frame, path) {
  write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

numeric_or_na <- function(value) suppressWarnings(as.numeric(value))
z_vector <- function(value) {
  value <- as.numeric(value)
  standard_deviation <- sd(value, na.rm = TRUE)
  if (!is.finite(standard_deviation) || standard_deviation == 0) return(rep(NA_real_, length(value)))
  (value - mean(value, na.rm = TRUE)) / standard_deviation
}

empty_result <- function() {
  data.frame(
    cohort = character(), model = character(), analysis_family = character(),
    timepoint = character(), outcome = character(), response_definition = character(),
    feature = character(), evidence_tier = character(), n = integer(), n_event = integer(),
    effect_type = character(), effect = numeric(), ci_low = numeric(), ci_high = numeric(),
    secondary_metric_type = character(), secondary_metric = numeric(),
    p_value = numeric(), q_value = numeric(), direction = character(), notes = character(),
    stringsAsFactors = FALSE
  )
}

make_result <- function(model, family, timepoint, outcome, response_definition, feature,
                        evidence_tier, n, n_event, effect_type, effect, ci_low = NA_real_,
                        ci_high = NA_real_, secondary_metric_type = "", secondary_metric = NA_real_,
                        p_value = NA_real_, direction = "", notes = "") {
  data.frame(
    cohort = "E-MTAB-3218", model = model, analysis_family = family,
    timepoint = timepoint, outcome = outcome, response_definition = response_definition,
    feature = feature, evidence_tier = evidence_tier, n = as.integer(n), n_event = as.integer(n_event),
    effect_type = effect_type, effect = as.numeric(effect), ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high), secondary_metric_type = secondary_metric_type,
    secondary_metric = as.numeric(secondary_metric), p_value = as.numeric(p_value),
    q_value = NA_real_, direction = direction, notes = notes,
    stringsAsFactors = FALSE
  )
}

auc_higher <- function(score, response) {
  keep <- is.finite(score) & response %in% c(0, 1)
  score <- score[keep]
  response <- response[keep]
  n_positive <- sum(response == 1)
  n_negative <- sum(response == 0)
  if (n_positive == 0 || n_negative == 0) return(NA_real_)
  rank_sum <- sum(rank(score, ties.method = "average")[response == 1])
  u_statistic <- rank_sum - n_positive * (n_positive + 1) / 2
  u_statistic / (n_positive * n_negative)
}

response_models <- function(frame, feature_name, response_name, timepoint_name, tier) {
  score <- numeric_or_na(frame[[feature_name]])
  response <- numeric_or_na(frame[[response_name]])
  keep <- is.finite(score) & response %in% c(0, 1)
  score <- score[keep]
  response <- response[keep]
  if (length(score) < 10L || length(unique(response)) != 2L || min(table(response)) < 3L) return(empty_result())
  score_z <- z_vector(score)
  fit <- try(glm(response ~ score_z, family = binomial()), silent = TRUE)
  logistic <- empty_result()
  if (!inherits(fit, "try-error") && is.finite(coef(summary(fit))[2L, 4L])) {
    coefficient <- coef(fit)[[2L]]
    standard_error <- coef(summary(fit))[2L, 2L]
    logistic <- make_result(
      "logistic_unadjusted", paste0("response_", response_name, "_", timepoint_name, "_logistic"),
      timepoint_name, response_name, response_name, feature_name, tier,
      length(score), sum(response == 1), "odds_ratio_per_sd", exp(coefficient),
      exp(coefficient - 1.96 * standard_error), exp(coefficient + 1.96 * standard_error),
      "auc_higher_is_responder", auc_higher(score, response), coef(summary(fit))[2L, 4L],
      ifelse(coefficient > 0, "higher_in_response", "lower_in_response")
    )
  }
  wilcox_fit <- suppressWarnings(wilcox.test(score[response == 1], score[response == 0], exact = FALSE))
  median_difference <- median(score[response == 1]) - median(score[response == 0])
  rank_test <- make_result(
    "wilcoxon_unadjusted", paste0("response_", response_name, "_", timepoint_name, "_wilcoxon"),
    timepoint_name, response_name, response_name, feature_name, tier,
    length(score), sum(response == 1), "median_difference_responder_minus_nonresponder",
    median_difference, secondary_metric_type = "auc_higher_is_responder",
    secondary_metric = auc_higher(score, response), p_value = wilcox_fit$p.value,
    direction = ifelse(median_difference > 0, "higher_in_response", "lower_in_response")
  )
  rbind(logistic, rank_test)
}

cox_model <- function(frame, feature_name, time_name, event_name, outcome_name, tier) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- numeric_or_na(frame[[feature_name]])
  time <- numeric_or_na(frame[[time_name]])
  event <- numeric_or_na(frame[[event_name]])
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  score <- score[keep]
  time <- time[keep]
  event <- event[keep]
  if (length(score) < 20L || sum(event) < 5L) return(empty_result())
  score_z <- z_vector(score)
  fit <- try(survival::coxph(survival::Surv(time, event) ~ score_z, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  fit_summary <- summary(fit)
  coefficient <- fit_summary$coef[1L, "coef"]
  standard_error <- fit_summary$coef[1L, "se(coef)"]
  make_result(
    "cox_unadjusted", paste0("survival_", outcome_name, "_baseline"), "baseline",
    outcome_name, "", feature_name, tier, length(score), sum(event), "hazard_ratio_per_sd",
    exp(coefficient), exp(coefficient - 1.96 * standard_error), exp(coefficient + 1.96 * standard_error),
    p_value = fit_summary$coef[1L, "Pr(>|z|)"],
    direction = ifelse(coefficient > 0, "higher_worse_survival", "higher_better_survival")
  )
}

expr <- as.matrix(read.delim(
  rma_path, header = TRUE, row.names = 1L, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = ""
))
storage.mode(expr) <- "double"
annotation <- read.delim(
  annotation_path, header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE
)
samples <- read.delim(sample_manifest_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
patients <- read.delim(patient_clinical_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
signatures <- read.delim(signature_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
signatures$weight <- signatures$raw_weight

samples$matrix_column <- paste0("ENTREZG-", sub("\\.cel$", "", samples$array_file, ignore.case = TRUE))
if (!all(samples$matrix_column %in% colnames(expr))) stop("Incomplete sample-to-expression mapping")
expr <- expr[, samples$matrix_column, drop = FALSE]

annotation_index <- match(rownames(expr), annotation$Probeset)
symbols <- trimws(annotation$Symbol[annotation_index])
keep_symbol <- !is.na(symbols) & nzchar(symbols) & symbols != "NA"
symbol_sums <- rowsum(expr[keep_symbol, , drop = FALSE], group = symbols[keep_symbol], reorder = TRUE)
symbol_counts <- table(symbols[keep_symbol])
gene_expr <- sweep(symbol_sums, 1L, as.numeric(symbol_counts[rownames(symbol_sums)]), "/")

f1 <- signatures[signatures$signature_id == "f1_top100_weighted", , drop = FALSE]
f1$weight <- as.numeric(f1$weight)
f1 <- f1[f1$gene %in% rownames(gene_expr) & is.finite(f1$weight) & f1$weight > 0, , drop = FALSE]
if (nrow(f1) < 70L) stop("F1 coverage below frozen threshold")
f1_expression <- gene_expr[f1$gene, , drop = FALSE]
f1_gene_z <- t(apply(f1_expression, 1L, z_vector))
f1_score_z <- colSums(f1_gene_z * f1$weight) / sum(f1$weight)
f1_score_raw <- colSums(f1_expression * f1$weight) / sum(f1$weight)
if (!"TYMP" %in% rownames(gene_expr)) stop("TYMP missing from collapsed gene matrix")
tymp_expression <- as.numeric(gene_expr["TYMP", ])
tymp_z <- z_vector(tymp_expression)
joint_axis_z <- (z_vector(f1_score_z) + tymp_z) / 2

samples$f1_score_gene_z_weighted <- as.numeric(f1_score_z)
samples$f1_score_rma_weighted <- as.numeric(f1_score_raw)
samples$tymp_expression_rma <- tymp_expression
samples$tymp_expression_z <- tymp_z
samples$f1_tymp_joint_axis_z <- joint_axis_z
samples$orr_primary <- numeric_or_na(samples$orr_response)
response20_lower <- tolower(trimws(samples$response20pct))
samples$response20_sensitivity <- ifelse(
  grepl("^20pct", response20_lower), 1,
  ifelse(response20_lower == "not20pct", 0, NA_real_)
)

features <- data.frame(
  feature = c("f1_score_gene_z_weighted", "tymp_expression_z", "f1_tymp_joint_axis_z"),
  tier = c("PRIMARY_F1", "PRIMARY_TYMP", "EXPLORATORY_JOINT"),
  stringsAsFactors = FALSE
)
results <- empty_result()

for (timepoint_name in c("baseline", "week4")) {
  current <- samples[samples$timepoint == timepoint_name, , drop = FALSE]
  for (feature_index in seq_len(nrow(features))) {
    feature_name <- features$feature[[feature_index]]
    tier <- features$tier[[feature_index]]
    results <- rbind(
      results,
      response_models(current, feature_name, "orr_primary", timepoint_name, tier),
      response_models(current, feature_name, "response20_sensitivity", timepoint_name, tier)
    )
  }
  correlation <- suppressWarnings(cor.test(
    current$f1_score_gene_z_weighted, current$tymp_expression_rma,
    method = "spearman", exact = FALSE
  ))
  results <- rbind(results, make_result(
    "spearman", paste0("f1_tymp_correlation_", timepoint_name), timepoint_name,
    "f1_tymp_correlation", "", "F1_vs_TYMP", "PRIMARY_RELATIONSHIP",
    nrow(current), 0L, "spearman_rho", unname(correlation$estimate),
    p_value = correlation$p.value,
    direction = ifelse(unname(correlation$estimate) > 0, "positive", "negative")
  ))
}

baseline <- samples[samples$timepoint == "baseline", , drop = FALSE]
for (feature_index in seq_len(nrow(features))) {
  feature_name <- features$feature[[feature_index]]
  tier <- features$tier[[feature_index]]
  results <- rbind(
    results,
    cox_model(baseline, feature_name, "pfs_days", "pfs_event", "PFS", tier),
    cox_model(baseline, feature_name, "os_days", "os_event", "OS", tier)
  )
}

paired <- merge(
  samples[samples$timepoint == "baseline", , drop = FALSE],
  samples[samples$timepoint == "week4", , drop = FALSE],
  by = "patient_id", suffixes = c("_baseline", "_week4"), all = FALSE
)
for (feature_index in seq_len(nrow(features))) {
  feature_name <- features$feature[[feature_index]]
  tier <- features$tier[[feature_index]]
  baseline_score <- numeric_or_na(paired[[paste0(feature_name, "_baseline")]])
  week4_score <- numeric_or_na(paired[[paste0(feature_name, "_week4")]])
  delta <- week4_score - baseline_score
  paired[[paste0("delta_", feature_name)]] <- delta
  keep <- is.finite(delta)
  paired_test <- suppressWarnings(wilcox.test(delta[keep], mu = 0, exact = FALSE))
  results <- rbind(results, make_result(
    "paired_wilcoxon", paste0("longitudinal_overall_", feature_name), "paired",
    "week4_minus_baseline", "", feature_name, tier, sum(keep), 0L,
    "median_paired_delta", median(delta[keep]), p_value = paired_test$p.value,
    direction = ifelse(median(delta[keep]) > 0, "increased_week4", "decreased_week4")
  ))
  for (response_suffix in c("orr_primary", "response20_sensitivity")) {
    response <- numeric_or_na(paired[[paste0(response_suffix, "_baseline")]])
    response_keep <- keep & response %in% c(0, 1)
    if (sum(response_keep) >= 10L && length(unique(response[response_keep])) == 2L && min(table(response[response_keep])) >= 3L) {
      delta_test <- suppressWarnings(wilcox.test(
        delta[response_keep & response == 1], delta[response_keep & response == 0], exact = FALSE
      ))
      difference <- median(delta[response_keep & response == 1]) - median(delta[response_keep & response == 0])
      results <- rbind(results, make_result(
        "delta_wilcoxon", paste0("longitudinal_delta_", response_suffix), "paired",
        "delta_by_response", response_suffix, feature_name, tier, sum(response_keep),
        sum(response[response_keep] == 1), "median_delta_difference_responder_minus_nonresponder",
        difference, secondary_metric_type = "auc_higher_delta_is_responder",
        secondary_metric = auc_higher(delta[response_keep], response[response_keep]),
        p_value = delta_test$p.value,
        direction = ifelse(difference > 0, "greater_in_response", "lower_in_response")
      ))
    }
  }
}

finite_p <- is.finite(results$p_value)
for (family_name in unique(results$analysis_family[finite_p])) {
  family_rows <- which(results$analysis_family == family_name & finite_p)
  results$q_value[family_rows] <- p.adjust(results$p_value[family_rows], method = "BH")
}
results$significance_class <- ifelse(
  is.finite(results$q_value) & results$q_value < 0.05, "FDR_LT_0.05",
  ifelse(is.finite(results$p_value) & results$p_value < 0.05, "NOMINAL_P_LT_0.05", "NOT_SIGNIFICANT")
)
results <- results[order(results$analysis_family, results$p_value, results$feature), , drop = FALSE]

score_summary <- do.call(rbind, lapply(c("baseline", "week4"), function(timepoint_name) {
  current <- samples[samples$timepoint == timepoint_name, , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(features)), function(feature_index) {
    value <- numeric_or_na(current[[features$feature[[feature_index]]]])
    data.frame(
      timepoint = timepoint_name, feature = features$feature[[feature_index]],
      evidence_tier = features$tier[[feature_index]], n = sum(is.finite(value)),
      mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE), q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE), stringsAsFactors = FALSE
    )
  }))
}))

write_tsv(samples, file.path(output_dir, "emtab3218_sample_scores_v1.tsv"))
write_tsv(paired, file.path(output_dir, "emtab3218_paired_scores_v1.tsv"))
write_tsv(score_summary, file.path(output_dir, "emtab3218_score_summary_v1.tsv"))
write_tsv(results, file.path(output_dir, "emtab3218_core_models_v1.tsv"))
write_tsv(
  results[results$significance_class != "NOT_SIGNIFICANT", , drop = FALSE],
  file.path(output_dir, "emtab3218_nominal_or_fdr_signals_v1.tsv")
)

status <- data.frame(
  gate = c(
    "sample_scores_complete", "paired_patients", "f1_genes_used", "core_model_rows",
    "fdr_signal_rows", "nominal_only_rows", "survival_package_available"
  ),
  observed = c(
    nrow(samples), nrow(paired), nrow(f1), nrow(results),
    sum(results$significance_class == "FDR_LT_0.05"),
    sum(results$significance_class == "NOMINAL_P_LT_0.05"),
    requireNamespace("survival", quietly = TRUE)
  ),
  expected_or_threshold = c(114, 42, 70, 20, 0, 0, TRUE),
  status = c(
    ifelse(nrow(samples) == 114L, "PASS", "FAIL"),
    ifelse(nrow(paired) == 42L, "PASS", "FAIL"),
    ifelse(nrow(f1) >= 70L, "PASS", "FAIL"),
    ifelse(nrow(results) >= 20L, "PASS", "FAIL"),
    "INFO", "INFO",
    ifelse(requireNamespace("survival", quietly = TRUE), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "emtab3218_core_validation_status_v1.tsv"))
print(status, row.names = FALSE)
print(results[results$significance_class != "NOT_SIGNIFICANT", c(
  "model", "timepoint", "outcome", "feature", "effect_type", "effect", "p_value", "q_value", "significance_class"
)], row.names = FALSE)
if (any(status$status == "FAIL")) quit(status = 2L)
cat("PASS_E_MTAB_3218_CORE_VALIDATION\n")
