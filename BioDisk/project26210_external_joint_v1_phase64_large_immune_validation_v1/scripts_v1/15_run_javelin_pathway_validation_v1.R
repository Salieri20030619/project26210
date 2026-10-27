args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript 15_run_javelin_pathway_validation_v1.R <pathway.tsv.gz> <sample_scores.tsv> <output_dir> <cohort_label>")
pathway_path <- args[[1L]]
sample_path <- args[[2L]]
output_dir <- args[[3L]]
cohort_label <- args[[4L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
z_vector <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
make_result <- function(model, family, stratum, pathway, n, n_event, effect_type, effect,
                        ci_low, ci_high, p_value, direction, notes = "") {
  data.frame(
    cohort = cohort_label, model = model, analysis_family = family, stratum = stratum,
    outcome = "PFS", pathway = pathway, n = as.integer(n), n_event = as.integer(n_event),
    effect_type = effect_type, effect = as.numeric(effect), ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high), p_value = as.numeric(p_value), q_value = NA_real_,
    direction = direction, theme_related = FALSE, retention_class = "", notes = notes,
    stringsAsFactors = FALSE
  )
}
empty_result <- function() data.frame(
  cohort = character(), model = character(), analysis_family = character(), stratum = character(),
  outcome = character(), pathway = character(), n = integer(), n_event = integer(),
  effect_type = character(), effect = numeric(), ci_low = numeric(), ci_high = numeric(),
  p_value = numeric(), q_value = numeric(), direction = character(), theme_related = logical(),
  retention_class = character(), notes = character(), stringsAsFactors = FALSE
)

within_arm <- function(frame, score, pathway, adjusted = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- z_vector(score)
  time <- as.numeric(frame$PFS_P)
  event <- 1 - as.numeric(frame$PFS_P_CNSR)
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  if (adjusted) keep <- keep & is.finite(as.numeric(frame$AGE)) & nzchar(frame$SEX)
  score <- score[keep]; time <- time[keep]; event <- event[keep]
  if (adjusted) {
    age <- z_vector(as.numeric(frame$AGE[keep])); sex <- factor(frame$SEX[keep])
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score + age + sex, ties = "efron"), silent = TRUE)
    model <- "cox_adjusted_age_sex"
  } else {
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score, ties = "efron"), silent = TRUE)
    model <- "cox_unadjusted"
  }
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef; beta <- tab["score", "coef"]; se <- tab["score", "se(coef)"]
  arm_slug <- gsub("[^a-z0-9]+", "_", tolower(unique(frame$TRT01P)))
  make_result(model, paste0("pathway_within_arm__", arm_slug, ifelse(adjusted, "__adjusted", "__unadjusted")),
              unique(frame$TRT01P), pathway, length(time), sum(event), "hazard_ratio_per_sd",
              exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), tab["score", "Pr(>|z|)"],
              ifelse(beta < 0, "higher_better_pfs", "higher_worse_pfs"))
}

interaction_test <- function(frame, score, pathway, adjusted = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  score <- z_vector(score)
  arm_combo <- as.integer(frame$TRT01P == "Avelumab+Axitinib")
  time <- as.numeric(frame$PFS_P)
  event <- 1 - as.numeric(frame$PFS_P_CNSR)
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  if (adjusted) keep <- keep & is.finite(as.numeric(frame$AGE)) & nzchar(frame$SEX)
  score <- score[keep]; arm_combo <- arm_combo[keep]; time <- time[keep]; event <- event[keep]
  if (adjusted) {
    age <- z_vector(as.numeric(frame$AGE[keep])); sex <- factor(frame$SEX[keep])
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score * arm_combo + age + sex, ties = "efron"), silent = TRUE)
    model <- "pfs_treatment_interaction_adjusted_age_sex"
  } else {
    fit <- try(survival::coxph(survival::Surv(time, event) ~ score * arm_combo, ties = "efron"), silent = TRUE)
    model <- "pfs_treatment_interaction_unadjusted"
  }
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef; beta <- tab["score:arm_combo", "coef"]; se <- tab["score:arm_combo", "se(coef)"]
  make_result(model, paste0("pathway_treatment_interaction", ifelse(adjusted, "__adjusted", "__unadjusted")),
              "randomized_all", pathway, length(time), sum(event), "interaction_hazard_ratio",
              exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), tab["score:arm_combo", "Pr(>|z|)"],
              ifelse(beta < 0, "higher_more_combo_benefit", "higher_more_sunitinib_benefit"))
}

median_test <- function(frame, score, pathway, global_median) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  high <- as.integer(score >= global_median)
  time <- as.numeric(frame$PFS_P)
  event <- 1 - as.numeric(frame$PFS_P_CNSR)
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  high <- high[keep]; time <- time[keep]; event <- event[keep]
  fit <- try(survival::coxph(survival::Surv(time, event) ~ high, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef; beta <- tab["high", "coef"]; se <- tab["high", "se(coef)"]
  arm_slug <- gsub("[^a-z0-9]+", "_", tolower(unique(frame$TRT01P)))
  make_result("cox_global_median_high_vs_low", paste0("pathway_median__", arm_slug), unique(frame$TRT01P),
              pathway, length(time), sum(event), "hazard_ratio_high_vs_low", exp(beta),
              exp(beta - 1.96 * se), exp(beta + 1.96 * se), tab["high", "Pr(>|z|)"],
              ifelse(beta < 0, "high_better_pfs", "high_worse_pfs"),
              paste0("global_median=", format(global_median, digits = 8)))
}

pathway <- read.delim(gzfile(pathway_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
samples <- read.delim(sample_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
index <- match(pathway$row_id, samples$ID)
if (anyNA(index)) stop("Pathway-to-sample mapping incomplete")
samples <- samples[index, , drop = FALSE]
pathway_names <- setdiff(colnames(pathway), "row_id")
combo_idx <- which(samples$TRT01P == "Avelumab+Axitinib")
suni_idx <- which(samples$TRT01P == "Sunitinib")

results <- vector("list", length(pathway_names) * 8L)
cursor <- 0L
for (pathway_name in pathway_names) {
  score <- as.numeric(pathway[[pathway_name]])
  global_median <- median(score, na.rm = TRUE)
  current <- list(
    interaction_test(samples, score, pathway_name, FALSE),
    interaction_test(samples, score, pathway_name, TRUE),
    within_arm(samples[combo_idx, , drop = FALSE], score[combo_idx], pathway_name, FALSE),
    within_arm(samples[combo_idx, , drop = FALSE], score[combo_idx], pathway_name, TRUE),
    within_arm(samples[suni_idx, , drop = FALSE], score[suni_idx], pathway_name, FALSE),
    within_arm(samples[suni_idx, , drop = FALSE], score[suni_idx], pathway_name, TRUE),
    median_test(samples[combo_idx, , drop = FALSE], score[combo_idx], pathway_name, global_median),
    median_test(samples[suni_idx, , drop = FALSE], score[suni_idx], pathway_name, global_median)
  )
  for (item in current) { cursor <- cursor + 1L; results[[cursor]] <- item }
}
results <- do.call(rbind, results[seq_len(cursor)])

theme_pattern <- paste(c("IMMUNE", "JAVELIN", "T_CELL", "TREG", "CD8", "CD4", "NK_CELL", "B_CELL", "PLASMA_CELL",
                         "MACROPHAGE", "MONOCYTE", "DENDRITIC", "INTERFERON", "TNFA", "NFKB", "TEFF", "MERCK",
                         "ANGIO", "MYELOID_INFLAMMATION", "NEUTROPHIL", "CYTOTOX"), collapse = "|")
results$theme_related <- grepl(theme_pattern, toupper(results$pathway))
for (family in unique(results$analysis_family)) {
  idx <- which(results$analysis_family == family & is.finite(results$p_value))
  results$q_value[idx] <- p.adjust(results$p_value[idx], method = "BH")
}
results$retention_class <- "NOT_RETAINED"
results$retention_class[is.finite(results$q_value) & results$q_value < 0.05] <- "C1_FDR_LT_0.05"
idx <- which(results$retention_class == "NOT_RETAINED" & results$theme_related & is.finite(results$p_value) & results$p_value < 0.05)
results$retention_class[idx] <- "C2_THEME_NOMINAL_P_LT_0.05"
for (family in unique(results$analysis_family)) {
  idx <- which(results$analysis_family == family & results$retention_class == "NOT_RETAINED" &
               is.finite(results$p_value) & results$p_value < 0.05)
  idx <- idx[order(results$p_value[idx], results$pathway[idx])]
  results$retention_class[head(idx, 15L)] <- "C3_TOP15_NOMINAL_P_LT_0.05"
}
results <- results[order(results$analysis_family, results$p_value, results$pathway), , drop = FALSE]
retained <- results[results$retention_class != "NOT_RETAINED", , drop = FALSE]
summary_rows <- do.call(rbind, lapply(unique(results$analysis_family), function(family) {
  x <- results[results$analysis_family == family, , drop = FALSE]
  data.frame(analysis_family = family, tested = nrow(x), nominal_p_lt_0.05 = sum(x$p_value < 0.05, na.rm = TRUE),
             fdr_lt_0.05 = sum(x$q_value < 0.05, na.rm = TRUE), theme_nominal = sum(x$theme_related & x$p_value < 0.05, na.rm = TRUE),
             retained = sum(x$retention_class != "NOT_RETAINED"), stringsAsFactors = FALSE)
}))

write_tsv(results, file.path(output_dir, "javelin_pathway_models_all_v1.tsv"))
write_tsv(retained, file.path(output_dir, "javelin_pathway_models_retained_v1.tsv"))
write_tsv(summary_rows, file.path(output_dir, "javelin_pathway_summary_v1.tsv"))

j26 <- results[results$analysis_family == "pathway_median__avelumab_axitinib" & results$pathway == "B9991003_Javelin_Renal_101_genes26", ]
angio <- results[results$analysis_family == "pathway_median__sunitinib" & results$pathway == "MCDERMOTT_ANGIO", ]
status <- data.frame(
  gate = c("analysis_samples", "pathways", "model_rows", "javelin26_positive_control", "angio_positive_control", "retention_rule_applied"),
  observed = c(nrow(samples), length(pathway_names), nrow(results), ifelse(nrow(j26), j26$p_value, NA), ifelse(nrow(angio), angio$p_value, NA), TRUE),
  expected_or_threshold = c(726, 58, 464, 0.01, 0.01, TRUE),
  status = c(ifelse(nrow(samples) == 726L, "PASS", "FAIL"), ifelse(length(pathway_names) == 58L, "PASS", "FAIL"),
             ifelse(nrow(results) == 464L, "PASS", "FAIL"), ifelse(nrow(j26) && j26$p_value < 0.01, "PASS", "FAIL"),
             ifelse(nrow(angio) && angio$p_value < 0.01, "PASS", "FAIL"), "PASS"), stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "javelin_pathway_validation_status_v1.tsv"))
print(status, row.names = FALSE)
print(summary_rows, row.names = FALSE)
print(retained[, c("analysis_family", "pathway", "effect", "p_value", "q_value", "retention_class")], row.names = FALSE)
