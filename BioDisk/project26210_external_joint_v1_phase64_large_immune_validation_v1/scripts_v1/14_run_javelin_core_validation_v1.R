args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript 14_run_javelin_core_validation_v1.R <f1_tymp.tsv.gz> <clinical.tsv> <pathway.tsv.gz> <signature.tsv> <output_dir> <cohort_label>")
}

expression_path <- args[[1L]]
clinical_path <- args[[2L]]
pathway_path <- args[[3L]]
signature_path <- args[[4L]]
output_dir <- args[[5L]]
cohort_label <- args[[6L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
z_vector <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
empty_result <- function() data.frame(
  cohort = character(), model = character(), analysis_family = character(), stratum = character(),
  outcome = character(), feature = character(), evidence_tier = character(), n = integer(),
  n_event = integer(), effect_type = character(), effect = numeric(), ci_low = numeric(),
  ci_high = numeric(), p_value = numeric(), q_value = numeric(), direction = character(),
  notes = character(), stringsAsFactors = FALSE
)
make_result <- function(model, family, stratum, feature, tier, n, n_event, effect_type,
                        effect, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
                        direction = "", notes = "") {
  data.frame(
    cohort = cohort_label, model = model, analysis_family = family, stratum = stratum,
    outcome = "PFS", feature = feature, evidence_tier = tier, n = as.integer(n),
    n_event = as.integer(n_event), effect_type = effect_type, effect = as.numeric(effect),
    ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high), p_value = as.numeric(p_value),
    q_value = NA_real_, direction = direction, notes = notes, stringsAsFactors = FALSE
  )
}

cox_within_arm <- function(frame, feature, tier, adjusted = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- z_vector(frame[[feature]])
  time <- suppressWarnings(as.numeric(frame$PFS_P))
  event <- 1 - suppressWarnings(as.numeric(frame$PFS_P_CNSR))
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  if (adjusted) keep <- keep & is.finite(as.numeric(frame$AGE)) & nzchar(frame$SEX)
  score <- score[keep]
  time <- time[keep]
  event <- event[keep]
  if (adjusted) {
    age <- z_vector(as.numeric(frame$AGE[keep]))
    sex <- factor(frame$SEX[keep])
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score + age + sex, ties = "efron"), silent = TRUE)
    model_name <- "cox_adjusted_age_sex"
  } else {
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score, ties = "efron"), silent = TRUE)
    model_name <- "cox_unadjusted"
  }
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  if (!"score" %in% rownames(tab)) return(empty_result())
  beta <- tab["score", "coef"]
  se <- tab["score", "se(coef)"]
  make_result(
    model_name, paste0("pfs_within_arm__", gsub("[^a-z0-9]+", "_", tolower(unique(frame$TRT01P))),
                       ifelse(adjusted, "__adjusted", "__unadjusted")),
    unique(frame$TRT01P), feature, tier, length(time), sum(event), "hazard_ratio_per_sd",
    exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), tab["score", "Pr(>|z|)"],
    ifelse(beta < 0, "higher_better_pfs", "higher_worse_pfs")
  )
}

cox_interaction <- function(frame, feature, tier, adjusted = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- z_vector(frame[[feature]])
  arm_combo <- as.integer(frame$TRT01P == "Avelumab+Axitinib")
  time <- suppressWarnings(as.numeric(frame$PFS_P))
  event <- 1 - suppressWarnings(as.numeric(frame$PFS_P_CNSR))
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1) & arm_combo %in% c(0, 1)
  if (adjusted) keep <- keep & is.finite(as.numeric(frame$AGE)) & nzchar(frame$SEX)
  score <- score[keep]
  arm_combo <- arm_combo[keep]
  time <- time[keep]
  event <- event[keep]
  if (adjusted) {
    age <- z_vector(as.numeric(frame$AGE[keep]))
    sex <- factor(frame$SEX[keep])
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score * arm_combo + age + sex, ties = "efron"), silent = TRUE)
    model_name <- "pfs_treatment_interaction_adjusted_age_sex"
  } else {
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score * arm_combo, ties = "efron"), silent = TRUE)
    model_name <- "pfs_treatment_interaction_unadjusted"
  }
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  coefficient_name <- "score:arm_combo"
  if (!coefficient_name %in% rownames(tab)) return(empty_result())
  beta <- tab[coefficient_name, "coef"]
  se <- tab[coefficient_name, "se(coef)"]
  make_result(
    model_name, paste0("pfs_treatment_interaction", ifelse(adjusted, "__adjusted", "__unadjusted")),
    "randomized_all", feature, tier, length(time), sum(event), "interaction_hazard_ratio",
    exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), tab[coefficient_name, "Pr(>|z|)"],
    ifelse(beta < 0, "higher_more_combo_benefit", "higher_more_sunitinib_benefit")
  )
}

cox_median <- function(frame, feature, tier, global_median) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- as.numeric(frame[[feature]])
  high <- as.integer(score >= global_median)
  time <- suppressWarnings(as.numeric(frame$PFS_P))
  event <- 1 - suppressWarnings(as.numeric(frame$PFS_P_CNSR))
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  high <- high[keep]
  time <- time[keep]
  event <- event[keep]
  fit <- try(survival::coxph(survival::Surv(time, event) ~ high, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  beta <- tab["high", "coef"]
  se <- tab["high", "se(coef)"]
  make_result(
    "cox_global_median_high_vs_low", paste0("pfs_median_sensitivity__", gsub("[^a-z0-9]+", "_", tolower(unique(frame$TRT01P)))),
    unique(frame$TRT01P), feature, tier, length(time), sum(event), "hazard_ratio_high_vs_low",
    exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), tab["high", "Pr(>|z|)"],
    ifelse(beta < 0, "high_better_pfs", "high_worse_pfs"),
    paste0("global_median=", format(global_median, digits = 8))
  )
}

expression <- read.delim(gzfile(expression_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_names <- as.character(expression[[1L]])
expr <- as.matrix(expression[, -1L, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- gene_names
if (anyDuplicated(rownames(expr))) stop("Duplicated target gene rows")

clinical <- read.delim(clinical_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
clinical_index <- match(colnames(expr), clinical$ID)
if (anyNA(clinical_index)) stop("Expression-to-clinical mapping incomplete")
samples <- clinical[clinical_index, , drop = FALSE]

pathway <- read.delim(gzfile(pathway_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
pathway_index <- match(samples$ID, pathway$row_id)
if (anyNA(pathway_index)) stop("Pathway-to-clinical mapping incomplete")
pathway <- pathway[pathway_index, , drop = FALSE]

signature <- read.delim(signature_path, header = TRUE, stringsAsFactors = FALSE)
f1 <- signature[signature$signature_id == "f1_top100_weighted" & signature$gene %in% rownames(expr), , drop = FALSE]
f1$raw_weight <- as.numeric(f1$raw_weight)
if (nrow(f1) < 70L) stop("F1 gene coverage below threshold")
f1_z <- t(apply(expr[f1$gene, , drop = FALSE], 1L, z_vector))
weight <- f1$raw_weight / sum(f1$raw_weight)
samples$f1_score <- as.numeric(colSums(f1_z * weight))
samples$tymp_expression <- as.numeric(expr["TYMP", ])
samples$tymp_z <- z_vector(samples$tymp_expression)
samples$f1_tymp_joint <- (z_vector(samples$f1_score) + samples$tymp_z) / 2
samples$treg_official <- as.numeric(pathway$LM22.T_cells_regulatory_Tregs)
samples$javelin26 <- as.numeric(pathway$B9991003_Javelin_Renal_101_genes26)
samples$teff_mcd <- as.numeric(pathway$MCDERMOTT_TEFF)
samples$ifng_response <- as.numeric(pathway$B9991003_c20_Interferon_gamma_response)
samples$angio_mcd <- as.numeric(pathway$MCDERMOTT_ANGIO)

features <- data.frame(
  feature = c("f1_score", "tymp_z", "f1_tymp_joint", "treg_official", "javelin26", "teff_mcd", "ifng_response", "angio_mcd"),
  tier = c("PRIMARY_F1", "PRIMARY_TYMP", "EXPLORATORY_JOINT", "PRIMARY_IMMUNE_TREG",
           "POSITIVE_CONTROL_JAVELIN26", "IMMUNE_CONTEXT_TEFF", "IMMUNE_CONTEXT_IFNG", "POSITIVE_CONTROL_ANGIO"),
  stringsAsFactors = FALSE
)

combo <- samples[samples$TRT01P == "Avelumab+Axitinib", , drop = FALSE]
sunitinib <- samples[samples$TRT01P == "Sunitinib", , drop = FALSE]
results <- empty_result()
for (i in seq_len(nrow(features))) {
  feature <- features$feature[[i]]
  tier <- features$tier[[i]]
  global_median <- median(samples[[feature]], na.rm = TRUE)
  results <- rbind(
    results,
    cox_interaction(samples, feature, tier, FALSE),
    cox_interaction(samples, feature, tier, TRUE),
    cox_within_arm(combo, feature, tier, FALSE),
    cox_within_arm(combo, feature, tier, TRUE),
    cox_within_arm(sunitinib, feature, tier, FALSE),
    cox_within_arm(sunitinib, feature, tier, TRUE),
    cox_median(combo, feature, tier, global_median),
    cox_median(sunitinib, feature, tier, global_median)
  )
}

correlation_pairs <- data.frame(
  x = c("f1_score", "f1_score", "tymp_z", "f1_tymp_joint"),
  y = c("tymp_z", "treg_official", "treg_official", "treg_official"),
  tier = c("PRIMARY_RELATIONSHIP", "PRIMARY_IMMUNE_TREG", "PRIMARY_IMMUNE_TREG", "EXPLORATORY_JOINT"),
  stringsAsFactors = FALSE
)
for (stratum_name in c("all", "combo", "sunitinib")) {
  current <- switch(stratum_name, all = samples, combo = combo, sunitinib = sunitinib)
  for (i in seq_len(nrow(correlation_pairs))) {
    x_name <- correlation_pairs$x[[i]]
    y_name <- correlation_pairs$y[[i]]
    x <- as.numeric(current[[x_name]])
    y <- as.numeric(current[[y_name]])
    keep <- is.finite(x) & is.finite(y)
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    rho <- unname(test$estimate)
    row <- make_result(
      "spearman", paste0("correlation__", stratum_name), stratum_name,
      paste(x_name, y_name, sep = "__vs__"), correlation_pairs$tier[[i]], sum(keep), 0L,
      "spearman_rho", rho, p_value = test$p.value, direction = ifelse(rho > 0, "positive", "negative")
    )
    row$outcome <- "correlation"
    results <- rbind(results, row)
  }
}

for (family in unique(results$analysis_family)) {
  idx <- which(results$analysis_family == family & is.finite(results$p_value))
  results$q_value[idx] <- p.adjust(results$p_value[idx], method = "BH")
}
results$significance_class <- ifelse(
  is.finite(results$q_value) & results$q_value < 0.05, "FDR_LT_0.05",
  ifelse(is.finite(results$p_value) & results$p_value < 0.05, "NOMINAL_P_LT_0.05", "NOT_SIGNIFICANT")
)
results <- results[order(results$analysis_family, results$p_value, results$feature), , drop = FALSE]

write_tsv(samples, file.path(output_dir, "javelin_sample_scores_v1.tsv"))
write_tsv(results, file.path(output_dir, "javelin_core_models_v1.tsv"))
write_tsv(results[results$significance_class != "NOT_SIGNIFICANT", , drop = FALSE],
          file.path(output_dir, "javelin_core_signals_v1.tsv"))

status <- data.frame(
  gate = c("analysis_samples", "combo_samples", "sunitinib_samples", "pfs_complete", "f1_genes_used", "model_rows", "survival_available"),
  observed = c(nrow(samples), nrow(combo), nrow(sunitinib), sum(is.finite(as.numeric(samples$PFS_P))), nrow(f1), nrow(results), requireNamespace("survival", quietly = TRUE)),
  expected_or_threshold = c(726, 300, 300, 726, 70, 70, TRUE),
  status = c(ifelse(nrow(samples) == 726L, "PASS", "FAIL"), ifelse(nrow(combo) >= 300L, "PASS", "FAIL"),
             ifelse(nrow(sunitinib) >= 300L, "PASS", "FAIL"), ifelse(sum(is.finite(as.numeric(samples$PFS_P))) == 726L, "PASS", "FAIL"),
             ifelse(nrow(f1) >= 70L, "PASS", "FAIL"), ifelse(nrow(results) >= 70L, "PASS", "FAIL"),
             ifelse(requireNamespace("survival", quietly = TRUE), "PASS", "FAIL")),
  stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "javelin_core_validation_status_v1.tsv"))
print(status, row.names = FALSE)
print(results[results$significance_class != "NOT_SIGNIFICANT", c(
  "model", "analysis_family", "stratum", "feature", "effect", "ci_low", "ci_high", "p_value", "q_value", "significance_class"
)], row.names = FALSE)
