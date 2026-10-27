args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript 10_run_checkmate_core_validation_v1.R <rna.tsv.gz> <clinical.tsv> <cibersort.tsv.gz> <signature.tsv> <output_dir> <cohort_name>")
}

rna_path <- args[[1L]]
clinical_path <- args[[2L]]
cibersort_path <- args[[3L]]
signature_path <- args[[4L]]
output_dir <- args[[5L]]
cohort_label <- args[[6L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(frame, path) {
  write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
z_vector <- function(value) {
  value <- as.numeric(value)
  standard_deviation <- sd(value, na.rm = TRUE)
  if (!is.finite(standard_deviation) || standard_deviation == 0) return(rep(NA_real_, length(value)))
  (value - mean(value, na.rm = TRUE)) / standard_deviation
}
z_by_group <- function(value, group) {
  output <- rep(NA_real_, length(value))
  for (level in unique(group)) {
    index <- which(group == level)
    output[index] <- z_vector(value[index])
  }
  output
}

empty_result <- function() data.frame(
  cohort = character(), model = character(), analysis_family = character(), stratum = character(),
  outcome = character(), feature = character(), evidence_tier = character(), n = integer(),
  n_event = integer(), effect_type = character(), effect = numeric(), ci_low = numeric(),
  ci_high = numeric(), p_value = numeric(), q_value = numeric(), direction = character(),
  notes = character(), stringsAsFactors = FALSE
)
make_result <- function(model, family, stratum, outcome, feature, tier, n, n_event,
                        effect_type, effect, ci_low = NA_real_, ci_high = NA_real_,
                        p_value = NA_real_, direction = "", notes = "") {
  data.frame(
    cohort = cohort_label, model = model, analysis_family = family, stratum = stratum,
    outcome = outcome, feature = feature, evidence_tier = tier, n = as.integer(n),
    n_event = as.integer(n_event), effect_type = effect_type, effect = as.numeric(effect),
    ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high), p_value = as.numeric(p_value),
    q_value = NA_real_, direction = direction, notes = notes, stringsAsFactors = FALSE
  )
}

logistic_model <- function(frame, feature, response, family, tier, adjusted = FALSE) {
  score <- suppressWarnings(as.numeric(frame[[feature]]))
  outcome <- suppressWarnings(as.numeric(frame[[response]]))
  keep <- is.finite(score) & outcome %in% c(0, 1)
  if (adjusted) keep <- keep & nzchar(frame$Cohort) & nzchar(frame$MSKCC)
  current <- frame[keep, , drop = FALSE]
  score <- z_vector(score[keep])
  outcome <- outcome[keep]
  if (length(outcome) < 30L || length(unique(outcome)) != 2L || min(table(outcome)) < 5L) return(empty_result())
  if (adjusted) {
    cohort <- factor(current$Cohort)
    risk <- factor(current$MSKCC)
    fit <- try(glm(outcome ~ score + cohort + risk, family = binomial()), silent = TRUE)
    model_name <- "logistic_adjusted_cohort_mskcc"
  } else {
    fit <- try(glm(outcome ~ score, family = binomial()), silent = TRUE)
    model_name <- "logistic_unadjusted"
  }
  if (inherits(fit, "try-error")) return(empty_result())
  fit_table <- coef(summary(fit))
  if (!"score" %in% rownames(fit_table) || !is.finite(fit_table["score", "Pr(>|z|)"])) return(empty_result())
  coefficient <- fit_table["score", "Estimate"]
  standard_error <- fit_table["score", "Std. Error"]
  make_result(
    model_name, paste0(family, ifelse(adjusted, "__adjusted", "__unadjusted")),
    "nivolumab", response, feature, tier, length(outcome), sum(outcome == 1),
    "odds_ratio_per_sd", exp(coefficient), exp(coefficient - 1.96 * standard_error),
    exp(coefficient + 1.96 * standard_error), fit_table["score", "Pr(>|z|)"],
    ifelse(coefficient > 0, "higher_in_response", "lower_in_response")
  )
}

cox_model <- function(frame, feature, time_column, event_column, outcome_name, family, tier, adjusted = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- suppressWarnings(as.numeric(frame[[feature]]))
  time <- suppressWarnings(as.numeric(frame[[time_column]]))
  event <- 1 - suppressWarnings(as.numeric(frame[[event_column]]))
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  if (adjusted) keep <- keep & nzchar(frame$Cohort) & nzchar(frame$MSKCC)
  current <- frame[keep, , drop = FALSE]
  score <- z_vector(score[keep])
  time <- time[keep]
  event <- event[keep]
  if (length(time) < 50L || sum(event) < 10L) return(empty_result())
  if (adjusted) {
    cohort <- factor(current$Cohort)
    risk <- factor(current$MSKCC)
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score + cohort + risk, ties = "efron"), silent = TRUE)
    model_name <- "cox_adjusted_cohort_mskcc"
  } else {
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score, ties = "efron"), silent = TRUE)
    model_name <- "cox_unadjusted"
  }
  if (inherits(fit, "try-error")) return(empty_result())
  fit_table <- summary(fit)$coef
  if (!"score" %in% rownames(fit_table)) return(empty_result())
  coefficient <- fit_table["score", "coef"]
  standard_error <- fit_table["score", "se(coef)"]
  make_result(
    model_name, paste0(family, ifelse(adjusted, "__adjusted", "__unadjusted")),
    "nivolumab", outcome_name, feature, tier, length(time), sum(event),
    "hazard_ratio_per_sd", exp(coefficient), exp(coefficient - 1.96 * standard_error),
    exp(coefficient + 1.96 * standard_error), fit_table["score", "Pr(>|z|)"],
    ifelse(coefficient > 0, "higher_worse_survival", "higher_better_survival")
  )
}

interaction_model <- function(frame, feature, outcome_name, tier) {
  score <- suppressWarnings(as.numeric(frame[[feature]]))
  arm_nivo <- as.integer(frame$Arm == "NIVOLUMAB")
  if (outcome_name == "ORR") {
    outcome <- frame$orr_binary
    keep <- is.finite(score) & outcome %in% c(0, 1) & arm_nivo %in% c(0, 1)
    score <- z_vector(score[keep])
    arm_nivo <- arm_nivo[keep]
    outcome <- outcome[keep]
    fit <- try(suppressWarnings(glm(outcome ~ score * arm_nivo, family = binomial())), silent = TRUE)
    coefficient_name <- "score:arm_nivo"
    if (inherits(fit, "try-error")) return(empty_result())
    fit_table <- coef(summary(fit))
    if (!coefficient_name %in% rownames(fit_table)) return(empty_result())
    coefficient <- fit_table[coefficient_name, "Estimate"]
    standard_error <- fit_table[coefficient_name, "Std. Error"]
    p_value <- fit_table[coefficient_name, "Pr(>|z|)"]
    effect_type <- "interaction_odds_ratio"
    n_event <- sum(outcome == 1)
    n <- length(outcome)
    direction <- ifelse(coefficient > 0, "higher_more_nivolumab_predictive", "higher_more_everolimus_predictive")
  } else {
    if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
    time_column <- ifelse(outcome_name == "PFS", "PFS", "OS")
    event_column <- ifelse(outcome_name == "PFS", "PFS_CNSR", "OS_CNSR")
    time <- suppressWarnings(as.numeric(frame[[time_column]]))
    event <- 1 - suppressWarnings(as.numeric(frame[[event_column]]))
    keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1) & arm_nivo %in% c(0, 1)
    score <- z_vector(score[keep])
    arm_nivo <- arm_nivo[keep]
    time <- time[keep]
    event <- event[keep]
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score * arm_nivo, ties = "efron"), silent = TRUE)
    coefficient_name <- "score:arm_nivo"
    if (inherits(fit, "try-error")) return(empty_result())
    fit_table <- summary(fit)$coef
    if (!coefficient_name %in% rownames(fit_table)) return(empty_result())
    coefficient <- fit_table[coefficient_name, "coef"]
    standard_error <- fit_table[coefficient_name, "se(coef)"]
    p_value <- fit_table[coefficient_name, "Pr(>|z|)"]
    effect_type <- "interaction_hazard_ratio"
    n_event <- sum(event)
    n <- length(time)
    direction <- ifelse(coefficient < 0, "higher_more_nivolumab_benefit", "higher_more_everolimus_benefit")
  }
  if (!all(is.finite(c(coefficient, standard_error, p_value))) || abs(coefficient) > 20 || standard_error > 20) {
    return(make_result(
      paste0(tolower(outcome_name), "_treatment_interaction"),
      paste0("cm025_treatment_interaction__", tolower(outcome_name)), "CM-025_randomized",
      outcome_name, feature, tier, n, n_event, effect_type, NA_real_, NA_real_, NA_real_,
      NA_real_, "not_estimable", "Complete/quasi-complete separation or unstable interaction estimate"
    ))
  }
  make_result(
    paste0(tolower(outcome_name), "_treatment_interaction"),
    paste0("cm025_treatment_interaction__", tolower(outcome_name)), "CM-025_randomized",
    outcome_name, feature, tier, n, n_event, effect_type, exp(coefficient),
    exp(coefficient - 1.96 * standard_error), exp(coefficient + 1.96 * standard_error),
    p_value, direction
  )
}

rna <- read.delim(gzfile(rna_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_names <- as.character(rna[[1L]])
expr <- as.matrix(rna[, -1L, drop = FALSE])
storage.mode(expr) <- "double"
rm(rna)
gene_sums <- rowsum(expr, group = gene_names, reorder = TRUE)
gene_counts <- table(gene_names)
gene_expr <- sweep(gene_sums, 1L, as.numeric(gene_counts[rownames(gene_sums)]), "/")
rm(expr, gene_sums)

clinical <- read.delim(clinical_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
clinical[is.na(clinical)] <- ""
clinical_index <- match(colnames(gene_expr), clinical$RNA_ID)
if (anyNA(clinical_index)) stop("RNA-to-clinical mapping incomplete")
samples <- clinical[clinical_index, , drop = FALSE]

signature <- read.delim(signature_path, header = TRUE, stringsAsFactors = FALSE)
f1 <- signature[signature$signature_id == "f1_top100_weighted", , drop = FALSE]
f1$weight <- as.numeric(f1$raw_weight)
f1 <- f1[f1$gene %in% rownames(gene_expr) & is.finite(f1$weight) & f1$weight > 0, , drop = FALSE]
if (nrow(f1) < 70L) stop("F1 gene coverage below threshold")

f1_expression <- gene_expr[f1$gene, , drop = FALSE]
global_gene_z <- t(apply(f1_expression, 1L, z_vector))
cohort_gene_z <- matrix(NA_real_, nrow = nrow(f1_expression), ncol = ncol(f1_expression))
for (cohort_name in unique(samples$Cohort)) {
  index <- which(samples$Cohort == cohort_name)
  cohort_gene_z[, index] <- t(apply(f1_expression[, index, drop = FALSE], 1L, z_vector))
}
weights <- f1$weight / sum(f1$weight)
samples$f1_score_global_gene_z <- as.numeric(colSums(global_gene_z * weights))
samples$f1_score_cohort_gene_z <- as.numeric(colSums(cohort_gene_z * weights))
samples$f1_score_cohort_z <- z_by_group(samples$f1_score_cohort_gene_z, samples$Cohort)

if (!"TYMP" %in% rownames(gene_expr)) stop("TYMP missing")
samples$tymp_expression <- as.numeric(gene_expr["TYMP", ])
samples$tymp_global_z <- z_vector(samples$tymp_expression)
samples$tymp_cohort_z <- z_by_group(samples$tymp_expression, samples$Cohort)
samples$f1_tymp_joint_cohort_z <- (samples$f1_score_cohort_z + samples$tymp_cohort_z) / 2

cibersort <- read.delim(gzfile(cibersort_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(cibersort) <- cibersort[[1L]]
cibersort <- as.matrix(cibersort[, -1L, drop = FALSE])
storage.mode(cibersort) <- "double"
if (!"t_cells_regulatory_tregs" %in% rownames(cibersort)) stop("Treg feature missing")
samples$cibersortx_treg <- as.numeric(cibersort["t_cells_regulatory_tregs", samples$RNA_ID])

samples$orr_binary <- ifelse(
  samples$ORR %in% c("CR", "PR", "CRPR"), 1,
  ifelse(samples$ORR %in% c("SD", "PD"), 0, NA_real_)
)
samples$benefit_binary <- ifelse(
  samples$Benefit == "CB", 1,
  ifelse(samples$Benefit == "NCB", 0, NA_real_)
)

features <- data.frame(
  feature = c(
    "f1_score_cohort_z", "tymp_cohort_z", "f1_tymp_joint_cohort_z", "cibersortx_treg",
    "f1_score_global_gene_z", "tymp_global_z"
  ),
  tier = c(
    "PRIMARY_F1", "PRIMARY_TYMP", "EXPLORATORY_JOINT", "PRIMARY_IMMUNE_TREG",
    "SENSITIVITY_GLOBAL_F1", "SENSITIVITY_GLOBAL_TYMP"
  ),
  stringsAsFactors = FALSE
)

nivolumab <- samples[samples$Arm == "NIVOLUMAB", , drop = FALSE]
cm025 <- samples[samples$Cohort == "CM-025", , drop = FALSE]
results <- empty_result()
for (feature_index in seq_len(nrow(features))) {
  feature <- features$feature[[feature_index]]
  tier <- features$tier[[feature_index]]
  results <- rbind(
    results,
    logistic_model(nivolumab, feature, "orr_binary", "nivolumab_orr", tier, FALSE),
    logistic_model(nivolumab, feature, "orr_binary", "nivolumab_orr", tier, TRUE),
    logistic_model(nivolumab, feature, "benefit_binary", "nivolumab_benefit", tier, FALSE),
    logistic_model(nivolumab, feature, "benefit_binary", "nivolumab_benefit", tier, TRUE),
    cox_model(nivolumab, feature, "PFS", "PFS_CNSR", "PFS", "nivolumab_pfs", tier, FALSE),
    cox_model(nivolumab, feature, "PFS", "PFS_CNSR", "PFS", "nivolumab_pfs", tier, TRUE),
    cox_model(nivolumab, feature, "OS", "OS_CNSR", "OS", "nivolumab_os", tier, FALSE),
    cox_model(nivolumab, feature, "OS", "OS_CNSR", "OS", "nivolumab_os", tier, TRUE),
    interaction_model(cm025, feature, "ORR", tier),
    interaction_model(cm025, feature, "PFS", tier),
    interaction_model(cm025, feature, "OS", tier)
  )
}

correlation_pairs <- data.frame(
  feature_x = c("f1_score_cohort_z", "f1_score_cohort_z", "tymp_cohort_z", "f1_tymp_joint_cohort_z"),
  feature_y = c("tymp_cohort_z", "cibersortx_treg", "cibersortx_treg", "cibersortx_treg"),
  tier = c("PRIMARY_RELATIONSHIP", "PRIMARY_IMMUNE_TREG", "PRIMARY_IMMUNE_TREG", "EXPLORATORY_JOINT"),
  stringsAsFactors = FALSE
)
strata <- list(all = samples, nivolumab = nivolumab, everolimus = samples[samples$Arm == "EVEROLIMUS", , drop = FALSE])
for (stratum_name in names(strata)) {
  current <- strata[[stratum_name]]
  for (pair_index in seq_len(nrow(correlation_pairs))) {
    x_name <- correlation_pairs$feature_x[[pair_index]]
    y_name <- correlation_pairs$feature_y[[pair_index]]
    x <- suppressWarnings(as.numeric(current[[x_name]]))
    y <- suppressWarnings(as.numeric(current[[y_name]]))
    keep <- is.finite(x) & is.finite(y)
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    rho <- unname(test$estimate)
    results <- rbind(results, make_result(
      "spearman", paste0("correlation__", stratum_name), stratum_name,
      paste(x_name, y_name, sep = "__vs__"), paste(x_name, y_name, sep = "__vs__"),
      correlation_pairs$tier[[pair_index]], sum(keep), 0L, "spearman_rho", rho,
      p_value = test$p.value, direction = ifelse(rho > 0, "positive", "negative")
    ))
  }
}

finite_p <- is.finite(results$p_value)
for (family in unique(results$analysis_family[finite_p])) {
  index <- which(results$analysis_family == family & finite_p)
  results$q_value[index] <- p.adjust(results$p_value[index], method = "BH")
}
results$significance_class <- ifelse(
  is.finite(results$q_value) & results$q_value < 0.05, "FDR_LT_0.05",
  ifelse(is.finite(results$p_value) & results$p_value < 0.05, "NOMINAL_P_LT_0.05", "NOT_SIGNIFICANT")
)
results <- results[order(results$analysis_family, results$p_value, results$feature), , drop = FALSE]

write_tsv(samples, file.path(output_dir, "checkmate_sample_scores_v1.tsv"))
write_tsv(results, file.path(output_dir, "checkmate_core_models_v1.tsv"))
write_tsv(results[results$significance_class != "NOT_SIGNIFICANT", , drop = FALSE],
          file.path(output_dir, "checkmate_core_signals_v1.tsv"))

status <- data.frame(
  gate = c(
    "rna_clinical_samples", "nivolumab_samples", "everolimus_samples", "cm025_randomized_samples",
    "f1_genes_used", "treg_values_finite", "model_rows", "survival_available"
  ),
  observed = c(
    nrow(samples), nrow(nivolumab), sum(samples$Arm == "EVEROLIMUS"), nrow(cm025),
    nrow(f1), sum(is.finite(samples$cibersortx_treg)), nrow(results),
    requireNamespace("survival", quietly = TRUE)
  ),
  expected_or_threshold = c(311, 181, 130, 250, 70, 311, 60, TRUE),
  status = c(
    ifelse(nrow(samples) == 311L, "PASS", "FAIL"),
    ifelse(nrow(nivolumab) == 181L, "PASS", "FAIL"),
    ifelse(sum(samples$Arm == "EVEROLIMUS") == 130L, "PASS", "FAIL"),
    ifelse(nrow(cm025) == 250L, "PASS", "FAIL"),
    ifelse(nrow(f1) >= 70L, "PASS", "FAIL"),
    ifelse(sum(is.finite(samples$cibersortx_treg)) == 311L, "PASS", "FAIL"),
    ifelse(nrow(results) >= 60L, "PASS", "FAIL"),
    ifelse(requireNamespace("survival", quietly = TRUE), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "checkmate_core_validation_status_v1.tsv"))
print(status, row.names = FALSE)
print(results[results$significance_class != "NOT_SIGNIFICANT", c(
  "model", "analysis_family", "stratum", "outcome", "feature", "effect", "ci_low", "ci_high", "p_value", "q_value", "significance_class"
)], row.names = FALSE)
if (any(status$status == "FAIL")) quit(status = 2L)
cat("PASS_CHECKMATE_CORE_VALIDATION\n")
