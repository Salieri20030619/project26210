args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 11_run_checkmate_pathway_validation_v1.R <ssgsea.tsv.gz> <sample_scores.tsv> <output_dir> <cohort_label>")
}

ssgsea_path <- args[[1L]]
sample_path <- args[[2L]]
output_dir <- args[[3L]]
cohort_label <- args[[4L]]
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
  outcome = character(), geneset_collection = character(), pathway = character(), n = integer(),
  n_event = integer(), effect_type = character(), effect = numeric(), ci_low = numeric(),
  ci_high = numeric(), p_value = numeric(), q_value = numeric(), direction = character(),
  theme_related = logical(), retention_class = character(), notes = character(), stringsAsFactors = FALSE
)

make_result <- function(model, family, outcome, collection, pathway, n, n_event, effect_type,
                        effect, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
                        direction = "", notes = "") {
  data.frame(
    cohort = cohort_label, model = model, analysis_family = family, stratum = "nivolumab",
    outcome = outcome, geneset_collection = collection, pathway = pathway, n = as.integer(n),
    n_event = as.integer(n_event), effect_type = effect_type, effect = as.numeric(effect),
    ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high), p_value = as.numeric(p_value),
    q_value = NA_real_, direction = direction, theme_related = FALSE,
    retention_class = "", notes = notes, stringsAsFactors = FALSE
  )
}

binary_test <- function(score, outcome, family, outcome_name, collection, pathway) {
  keep <- is.finite(score) & outcome %in% c(0, 1)
  score <- score[keep]
  outcome <- outcome[keep]
  if (length(score) < 30L || length(unique(outcome)) != 2L || min(table(outcome)) < 5L) return(empty_result())
  test <- suppressWarnings(wilcox.test(score[outcome == 1], score[outcome == 0], exact = FALSE))
  n1 <- sum(outcome == 1)
  n0 <- sum(outcome == 0)
  w <- unname(test$statistic)
  auc <- w / (n1 * n0)
  difference <- median(score[outcome == 1]) - median(score[outcome == 0])
  make_result(
    "wilcoxon", family, outcome_name, collection, pathway, length(score), n1,
    "median_difference_response_minus_nonresponse", difference, p_value = test$p.value,
    direction = ifelse(difference > 0, "higher_in_response", "lower_in_response"),
    notes = paste0("AUC=", format(auc, digits = 7))
  )
}

binary_adjusted_test <- function(score, outcome, cohort, risk, family, outcome_name, collection, pathway) {
  keep <- is.finite(score) & outcome %in% c(0, 1) & !is.na(cohort) & nzchar(cohort) & !is.na(risk) & nzchar(risk)
  score <- z_vector(score[keep])
  outcome <- outcome[keep]
  cohort <- factor(cohort[keep])
  risk <- factor(risk[keep])
  if (length(score) < 30L || length(unique(outcome)) != 2L || min(table(outcome)) < 5L) return(empty_result())
  fit <- try(suppressWarnings(glm(outcome ~ score + cohort + risk, family = binomial())), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- coef(summary(fit))
  if (!"score" %in% rownames(tab)) return(empty_result())
  beta <- tab["score", "Estimate"]
  se <- tab["score", "Std. Error"]
  p <- tab["score", "Pr(>|z|)"]
  if (!all(is.finite(c(beta, se, p))) || abs(beta) > 20 || se > 20) return(empty_result())
  make_result(
    "logistic_adjusted_cohort_mskcc", family, outcome_name, collection, pathway, length(score), sum(outcome == 1),
    "odds_ratio_per_sd", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p,
    ifelse(beta > 0, "higher_in_response", "lower_in_response")
  )
}

cox_test <- function(score, time, event, family, outcome_name, collection, pathway) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  score <- z_vector(score[keep])
  time <- time[keep]
  event <- event[keep]
  if (length(score) < 50L || sum(event) < 10L || !all(is.finite(score))) return(empty_result())
  fit <- try(survival::coxph(survival::Surv(time, event) ~ score, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  if (!"score" %in% rownames(tab)) return(empty_result())
  beta <- tab["score", "coef"]
  se <- tab["score", "se(coef)"]
  p <- tab["score", "Pr(>|z|)"]
  make_result(
    "cox_unadjusted", family, outcome_name, collection, pathway, length(score), sum(event),
    "hazard_ratio_per_sd", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p,
    ifelse(beta > 0, "higher_worse_survival", "higher_better_survival")
  )
}

cox_adjusted_test <- function(score, time, event, cohort, risk, family, outcome_name, collection, pathway) {
  if (!requireNamespace("survival", quietly = TRUE)) return(empty_result())
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1) &
    !is.na(cohort) & nzchar(cohort) & !is.na(risk) & nzchar(risk)
  score <- z_vector(score[keep])
  time <- time[keep]
  event <- event[keep]
  cohort <- factor(cohort[keep])
  risk <- factor(risk[keep])
  if (length(score) < 50L || sum(event) < 10L || !all(is.finite(score))) return(empty_result())
  fit <- try(survival::coxph(survival::Surv(time, event) ~ score + cohort + risk, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  if (!"score" %in% rownames(tab)) return(empty_result())
  beta <- tab["score", "coef"]
  se <- tab["score", "se(coef)"]
  p <- tab["score", "Pr(>|z|)"]
  make_result(
    "cox_adjusted_cohort_mskcc", family, outcome_name, collection, pathway, length(score), sum(event),
    "hazard_ratio_per_sd", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p,
    ifelse(beta > 0, "higher_worse_survival", "higher_better_survival")
  )
}

correlation_test <- function(score, target, family, outcome_name, collection, pathway) {
  keep <- is.finite(score) & is.finite(target)
  if (sum(keep) < 30L) return(empty_result())
  test <- suppressWarnings(cor.test(score[keep], target[keep], method = "spearman", exact = FALSE))
  rho <- unname(test$estimate)
  make_result(
    "spearman", family, outcome_name, collection, pathway, sum(keep), 0L,
    "spearman_rho", rho, p_value = test$p.value,
    direction = ifelse(rho > 0, "positive", "negative")
  )
}

ssgsea <- read.delim(gzfile(ssgsea_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ssgsea) < 300L || nrow(ssgsea) < 1600L) stop("ssGSEA matrix below expected dimensions")
collections <- as.character(ssgsea[[1L]])
pathways <- as.character(ssgsea[[2L]])
scores <- as.matrix(ssgsea[, -(1:2), drop = FALSE])
storage.mode(scores) <- "double"
rm(ssgsea)

samples <- read.delim(sample_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
sample_index <- match(colnames(scores), samples$RNA_ID)
if (anyNA(sample_index)) stop("ssGSEA-to-sample mapping incomplete")
samples <- samples[sample_index, , drop = FALSE]
nivo_index <- which(samples$Arm == "NIVOLUMAB")
nivo <- samples[nivo_index, , drop = FALSE]
nivo_scores <- scores[, nivo_index, drop = FALSE]
rm(scores)

orr <- suppressWarnings(as.numeric(nivo$orr_binary))
benefit <- suppressWarnings(as.numeric(nivo$benefit_binary))
pfs <- suppressWarnings(as.numeric(nivo$PFS))
pfs_event <- 1 - suppressWarnings(as.numeric(nivo$PFS_CNSR))
os <- suppressWarnings(as.numeric(nivo$OS))
os_event <- 1 - suppressWarnings(as.numeric(nivo$OS_CNSR))
f1 <- suppressWarnings(as.numeric(nivo$f1_score_cohort_z))
tymp <- suppressWarnings(as.numeric(nivo$tymp_cohort_z))
treg <- suppressWarnings(as.numeric(nivo$cibersortx_treg))
cohort <- as.character(nivo$Cohort)
risk <- as.character(nivo$MSKCC)

results <- vector("list", nrow(nivo_scores) * 11L)
cursor <- 0L
for (i in seq_len(nrow(nivo_scores))) {
  score <- as.numeric(nivo_scores[i, ])
  collection <- collections[[i]]
  pathway <- pathways[[i]]
  current <- list(
    binary_test(score, orr, "nivolumab_pathway_orr", "ORR", collection, pathway),
    binary_adjusted_test(score, orr, cohort, risk, "nivolumab_pathway_orr_adjusted", "ORR", collection, pathway),
    binary_test(score, benefit, "nivolumab_pathway_benefit", "Benefit", collection, pathway),
    binary_adjusted_test(score, benefit, cohort, risk, "nivolumab_pathway_benefit_adjusted", "Benefit", collection, pathway),
    cox_test(score, pfs, pfs_event, "nivolumab_pathway_pfs", "PFS", collection, pathway),
    cox_adjusted_test(score, pfs, pfs_event, cohort, risk, "nivolumab_pathway_pfs_adjusted", "PFS", collection, pathway),
    cox_test(score, os, os_event, "nivolumab_pathway_os", "OS", collection, pathway),
    cox_adjusted_test(score, os, os_event, cohort, risk, "nivolumab_pathway_os_adjusted", "OS", collection, pathway),
    correlation_test(score, f1, "pathway_correlation_f1", "F1", collection, pathway),
    correlation_test(score, tymp, "pathway_correlation_tymp", "TYMP", collection, pathway),
    correlation_test(score, treg, "pathway_correlation_treg", "CIBERSORTx_Treg", collection, pathway)
  )
  for (item in current) {
    cursor <- cursor + 1L
    results[[cursor]] <- item
  }
}
results <- do.call(rbind, results[seq_len(cursor)])

theme_pattern <- paste(c(
  "IMMUNE", "INTERFERON", "(^|_)IFN([A-Z0-9_]|$)", "(^|_)T_CELLS?(_|$)", "T_CELL_RECEPTOR",
  "(^|_)TCR(_|$)", "(^|_)CD8(_|$)", "(^|_)CD4(_|$)", "REGULATORY_T", "TREG", "CYTOTOX",
  "ANTIGEN", "(^|_)MHC(_|$)", "(^|_)HLA(_|$)", "PD_1", "PD1", "PD_L1", "PDL1", "CTLA4",
  "CYTOKINE", "CHEMOKINE", "INFLAM", "(^|_)TNF(_|$)", "NFKB", "JAK_STAT", "(^|_)IL2(_|$)",
  "(^|_)IL6(_|$)", "(^|_)IL10(_|$)", "(^|_)IL12(_|$)", "(^|_)IL17(_|$)", "TGFB", "TGF_BETA",
  "MACROPHAGE", "MONOCYTE", "DENDRITIC", "NK_CELL", "NATURAL_KILLER", "(^|_)B_CELLS?(_|$)",
  "LYMPHOCYTE", "LEUKOCYTE", "IMMUNOTHERAP", "NIVOLUMAB", "IPILIMUMAB", "AVELUMAB", "JAVELIN",
  "IMOTION", "TEFF"
), collapse = "|")
results$theme_related <- grepl(theme_pattern, toupper(paste(results$geneset_collection, results$pathway, sep = "__")))

for (family in unique(results$analysis_family)) {
  idx <- which(results$analysis_family == family & is.finite(results$p_value))
  results$q_value[idx] <- p.adjust(results$p_value[idx], method = "BH")
}

results$retention_class <- "NOT_RETAINED"
fdr_idx <- which(is.finite(results$q_value) & results$q_value < 0.05)
results$retention_class[fdr_idx] <- "C1_FDR_LT_0.05"
theme_idx <- which(results$retention_class == "NOT_RETAINED" & results$theme_related &
                   is.finite(results$p_value) & results$p_value < 0.05)
results$retention_class[theme_idx] <- "C2_THEME_NOMINAL_P_LT_0.05"
for (family in unique(results$analysis_family)) {
  idx <- which(results$analysis_family == family & results$retention_class == "NOT_RETAINED" &
               is.finite(results$p_value) & results$p_value < 0.05)
  if (length(idx)) {
    idx <- idx[order(results$p_value[idx], results$pathway[idx])]
    idx <- head(idx, 15L)
    results$retention_class[idx] <- "C3_TOP15_NOMINAL_P_LT_0.05"
  }
}
results <- results[order(results$analysis_family, results$p_value, results$pathway), , drop = FALSE]
retained <- results[results$retention_class != "NOT_RETAINED", , drop = FALSE]

summary_rows <- do.call(rbind, lapply(unique(results$analysis_family), function(family) {
  x <- results[results$analysis_family == family, , drop = FALSE]
  data.frame(
    analysis_family = family, tested = nrow(x), nominal_p_lt_0.05 = sum(x$p_value < 0.05, na.rm = TRUE),
    fdr_lt_0.05 = sum(x$q_value < 0.05, na.rm = TRUE), theme_nominal = sum(x$theme_related & x$p_value < 0.05, na.rm = TRUE),
    retained = sum(x$retention_class != "NOT_RETAINED"), stringsAsFactors = FALSE
  )
}))

write_tsv(results, file.path(output_dir, "checkmate_pathway_models_all_v1.tsv"))
write_tsv(retained, file.path(output_dir, "checkmate_pathway_models_retained_v1.tsv"))
write_tsv(summary_rows, file.path(output_dir, "checkmate_pathway_summary_v1.tsv"))

status <- data.frame(
  gate = c("ssgsea_pathways", "nivolumab_samples", "analysis_families", "model_rows", "retention_rule_applied"),
  observed = c(nrow(nivo_scores), nrow(nivo), length(unique(results$analysis_family)), nrow(results), TRUE),
  expected_or_threshold = c(1704, 181, 11, 18744, TRUE),
  status = c(
    ifelse(nrow(nivo_scores) == 1704L, "PASS", "FAIL"), ifelse(nrow(nivo) == 181L, "PASS", "FAIL"),
    ifelse(length(unique(results$analysis_family)) == 11L, "PASS", "FAIL"),
    ifelse(nrow(results) == 18744L, "PASS", "FAIL"), "PASS"
  ), stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "checkmate_pathway_validation_status_v1.tsv"))
print(status, row.names = FALSE)
print(summary_rows, row.names = FALSE)
print(head(retained[, c("analysis_family", "pathway", "effect", "p_value", "q_value", "retention_class")], 80L), row.names = FALSE)
