args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(paste(
    "Usage: Rscript 02_run_r2_composite_clinical_validation_v1.R",
    "<members.tsv> <emtab_rma.tsv> <emtab_annotation.tsv> <emtab_samples.tsv>",
    "<emtab_patients.tsv> <checkmate_expr.tsv.gz> <checkmate_clinical.tsv>",
    "<javelin_expr.tsv.gz> <javelin_clinical.tsv> <output_dir>"
  ))
}

members_path <- args[[1L]]
emtab_rma_path <- args[[2L]]
emtab_annotation_path <- args[[3L]]
emtab_samples_path <- args[[4L]]
emtab_patients_path <- args[[5L]]
checkmate_expr_path <- args[[6L]]
checkmate_clinical_path <- args[[7L]]
javelin_expr_path <- args[[8L]]
javelin_clinical_path <- args[[9L]]
output_dir <- args[[10L]]

if (dir.exists(output_dir)) stop("Refusing to overwrite an existing output directory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else path
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

z_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_group <- function(x, group) {
  out <- rep(NA_real_, length(x))
  for (level in unique(group)) {
    idx <- which(group == level)
    out[idx] <- z_vector(x[idx])
  }
  out
}

collapse_gene_rows <- function(expr, genes) {
  genes <- toupper(trimws(as.character(genes)))
  keep <- !is.na(genes) & nzchar(genes) & genes != "NA"
  expr <- expr[keep, , drop = FALSE]
  genes <- genes[keep]
  sums <- rowsum(expr, group = genes, reorder = TRUE)
  counts <- table(genes)
  sweep(sums, 1L, as.numeric(counts[rownames(sums)]), "/")
}

score_core <- function(gene_expr, groups, dataset, members) {
  required_genes <- unique(members$gene)
  present_required <- intersect(required_genes, rownames(gene_expr))
  expr <- gene_expr[present_required, , drop = FALSE]
  gene_z <- matrix(NA_real_, nrow = nrow(expr), ncol = ncol(expr), dimnames = dimnames(expr))
  for (level in unique(groups)) {
    idx <- which(groups == level)
    current <- expr[, idx, drop = FALSE]
    means <- rowMeans(current, na.rm = TRUE)
    sds <- apply(current, 1L, sd, na.rm = TRUE)
    good <- is.finite(sds) & sds > 0
    if (any(good)) {
      gene_z[good, idx] <- sweep(sweep(current[good, , drop = FALSE], 1L, means[good], "-"), 1L, sds[good], "/")
    }
  }

  signatures <- unique(members[, c(
    "core_order", "family", "signature_id", "label", "expected_phase64_direction", "orientation"
  )])
  signatures <- signatures[order(signatures$core_order), , drop = FALSE]
  score_matrix <- matrix(
    NA_real_, nrow = nrow(signatures), ncol = ncol(gene_expr),
    dimnames = list(signatures$signature_id, colnames(gene_expr))
  )
  coverage_rows <- list()
  for (i in seq_len(nrow(signatures))) {
    sig <- signatures[i, , drop = FALSE]
    frozen <- unique(members$gene[members$signature_id == sig$signature_id])
    present <- intersect(frozen, rownames(gene_z))
    fraction <- length(present) / length(frozen)
    passed <- length(present) >= 5L && fraction >= 0.50
    if (passed) {
      raw <- colMeans(gene_z[present, , drop = FALSE], na.rm = TRUE)
      oriented <- raw * as.numeric(sig$orientation)
      score_matrix[i, ] <- z_by_group(oriented, groups)
    }
    coverage_rows[[i]] <- data.frame(
      dataset = dataset,
      core_order = as.integer(sig$core_order),
      family = sig$family,
      signature_id = sig$signature_id,
      label = sig$label,
      expected_phase64_direction = sig$expected_phase64_direction,
      n_frozen_genes = length(frozen),
      n_present_genes = length(present),
      coverage_fraction = fraction,
      coverage_pass = passed,
      stringsAsFactors = FALSE
    )
  }
  coverage <- do.call(rbind, coverage_rows)
  metabolic <- which(signatures$family == "metabolic" & coverage$coverage_pass)
  regulon <- which(signatures$family == "regulon" & coverage$coverage_pass)
  if (length(metabolic) < 7L || length(regulon) < 7L) {
    stop(sprintf("Coverage family gate failed for %s: metabolic=%d regulon=%d", dataset, length(metabolic), length(regulon)))
  }
  metabolic_raw <- colMeans(score_matrix[metabolic, , drop = FALSE], na.rm = TRUE)
  regulon_raw <- colMeans(score_matrix[regulon, , drop = FALSE], na.rm = TRUE)
  metabolic_core <- z_by_group(metabolic_raw, groups)
  regulon_core <- z_by_group(regulon_raw, groups)
  combined_core <- z_by_group((metabolic_core + regulon_core) / 2, groups)
  combined_rank <- rep(NA_real_, length(groups))
  for (level in unique(groups)) {
    idx <- which(groups == level)
    meta_rank <- rank(metabolic_core[idx], ties.method = "average") / (length(idx) + 1)
    reg_rank <- rank(regulon_core[idx], ties.method = "average") / (length(idx) + 1)
    combined_rank[idx] <- z_vector((meta_rank + reg_rank) / 2)
  }
  tymp_z <- if ("TYMP" %in% rownames(gene_expr)) z_by_group(gene_expr["TYMP", ], groups) else rep(NA_real_, ncol(gene_expr))
  composite <- data.frame(
    sample_id = colnames(gene_expr),
    normalization_group = groups,
    r2_metabolic_core = metabolic_core,
    r2_regulon_core = regulon_core,
    r2_combined_core = combined_core,
    r2_combined_rank = combined_rank,
    tymp_expression_z = tymp_z,
    stringsAsFactors = FALSE
  )
  signature_long <- data.frame(
    dataset = dataset,
    sample_id = rep(colnames(score_matrix), each = nrow(score_matrix)),
    signature_id = rep(rownames(score_matrix), times = ncol(score_matrix)),
    oriented_signature_z = as.vector(score_matrix),
    stringsAsFactors = FALSE
  )
  list(composite = composite, coverage = coverage, signature_long = signature_long)
}

empty_result <- function() data.frame(
  cohort = character(), model_id = character(), analysis_tier = character(), endpoint = character(),
  stratum = character(), feature = character(), n = integer(), n_event = integer(),
  effect_type = character(), effect = numeric(), ci_low = numeric(), ci_high = numeric(),
  p_value = numeric(), q_value = numeric(), direction = character(), adjustment = character(),
  notes = character(), stringsAsFactors = FALSE
)

make_result <- function(cohort, model_id, tier, endpoint, stratum, feature, n, n_event,
                        effect_type, effect, ci_low, ci_high, p_value, direction,
                        adjustment = "none", notes = "") {
  data.frame(
    cohort = cohort, model_id = model_id, analysis_tier = tier, endpoint = endpoint,
    stratum = stratum, feature = feature, n = as.integer(n), n_event = as.integer(n_event),
    effect_type = effect_type, effect = as.numeric(effect), ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high), p_value = as.numeric(p_value), q_value = NA_real_,
    direction = direction, adjustment = adjustment, notes = notes, stringsAsFactors = FALSE
  )
}

fit_logistic <- function(frame, feature, outcome_column, cohort, model_id, tier, stratum,
                         covariate_columns = character(), covariate_terms = character(), adjustment = "none") {
  score <- suppressWarnings(as.numeric(frame[[feature]]))
  outcome <- suppressWarnings(as.numeric(frame[[outcome_column]]))
  keep <- is.finite(score) & outcome %in% c(0, 1)
  for (column in covariate_columns) keep <- keep & !is.na(frame[[column]]) & nzchar(as.character(frame[[column]]))
  current <- frame[keep, , drop = FALSE]
  current$.outcome <- outcome[keep]
  current$score <- z_vector(score[keep])
  if (nrow(current) < 10L || length(unique(current$.outcome)) != 2L || min(table(current$.outcome)) < 3L) return(empty_result())
  formula <- as.formula(paste(".outcome ~ score", if (length(covariate_terms)) paste("+", paste(covariate_terms, collapse = " + ")) else ""))
  fit <- try(suppressWarnings(glm(formula, data = current, family = binomial())), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- coef(summary(fit))
  if (!"score" %in% rownames(tab)) return(empty_result())
  beta <- tab["score", "Estimate"]
  se <- tab["score", "Std. Error"]
  p <- tab["score", "Pr(>|z|)"]
  if (!all(is.finite(c(beta, se, p))) || abs(beta) > 20 || se > 20) {
    return(make_result(cohort, model_id, tier, outcome_column, stratum, feature, nrow(current), sum(current$.outcome),
      "odds_ratio_per_sd", NA, NA, NA, NA, "not_estimable", adjustment, "Unstable logistic estimate"))
  }
  make_result(
    cohort, model_id, tier, outcome_column, stratum, feature, nrow(current), sum(current$.outcome),
    "odds_ratio_per_sd", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p,
    ifelse(beta > 0, "higher_in_response", "lower_in_response"), adjustment
  )
}

fit_cox <- function(frame, feature, time_column, event_column, event_is_censor, cohort, model_id,
                    tier, endpoint, stratum, covariate_columns = character(),
                    covariate_terms = character(), adjustment = "none") {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package unavailable")
  score <- suppressWarnings(as.numeric(frame[[feature]]))
  time <- suppressWarnings(as.numeric(frame[[time_column]]))
  raw_event <- suppressWarnings(as.numeric(frame[[event_column]]))
  event <- if (event_is_censor) 1 - raw_event else raw_event
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1)
  for (column in covariate_columns) keep <- keep & !is.na(frame[[column]]) & nzchar(as.character(frame[[column]]))
  current <- frame[keep, , drop = FALSE]
  current$.time <- time[keep]
  current$.event <- event[keep]
  current$score <- z_vector(score[keep])
  if (nrow(current) < 20L || sum(current$.event) < 5L) return(empty_result())
  formula <- as.formula(paste("survival::Surv(.time, .event) ~ score", if (length(covariate_terms)) paste("+", paste(covariate_terms, collapse = " + ")) else ""))
  fit <- try(survival::coxph(formula, data = current, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  if (!"score" %in% rownames(tab)) return(empty_result())
  beta <- tab["score", "coef"]
  se <- tab["score", "se(coef)"]
  p <- tab["score", "Pr(>|z|)"]
  make_result(
    cohort, model_id, tier, endpoint, stratum, feature, nrow(current), sum(current$.event),
    "hazard_ratio_per_sd", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p,
    ifelse(beta > 0, "higher_worse_survival", "higher_better_survival"), adjustment
  )
}

fit_cox_interaction <- function(frame, feature, time_column, event_column, event_is_censor,
                                arm_indicator, cohort, model_id, tier, endpoint,
                                covariate_columns = character(), covariate_terms = character(),
                                adjustment = "none", positive_direction, negative_direction) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package unavailable")
  score <- suppressWarnings(as.numeric(frame[[feature]]))
  time <- suppressWarnings(as.numeric(frame[[time_column]]))
  raw_event <- suppressWarnings(as.numeric(frame[[event_column]]))
  event <- if (event_is_censor) 1 - raw_event else raw_event
  arm <- as.integer(arm_indicator)
  keep <- is.finite(score) & is.finite(time) & time > 0 & event %in% c(0, 1) & arm %in% c(0, 1)
  for (column in covariate_columns) keep <- keep & !is.na(frame[[column]]) & nzchar(as.character(frame[[column]]))
  current <- frame[keep, , drop = FALSE]
  current$.time <- time[keep]
  current$.event <- event[keep]
  current$.arm <- arm[keep]
  current$score <- z_vector(score[keep])
  formula <- as.formula(paste("survival::Surv(.time, .event) ~ score * .arm", if (length(covariate_terms)) paste("+", paste(covariate_terms, collapse = " + ")) else ""))
  fit <- try(survival::coxph(formula, data = current, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(empty_result())
  tab <- summary(fit)$coef
  coefficient_name <- "score:.arm"
  if (!coefficient_name %in% rownames(tab)) return(empty_result())
  beta <- tab[coefficient_name, "coef"]
  se <- tab[coefficient_name, "se(coef)"]
  p <- tab[coefficient_name, "Pr(>|z|)"]
  make_result(
    cohort, model_id, tier, endpoint, "randomized_all", feature, nrow(current), sum(current$.event),
    "interaction_hazard_ratio", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se), p,
    ifelse(beta < 0, negative_direction, positive_direction), adjustment
  )
}

members <- read.delim(members_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
members$gene <- toupper(members$gene)
members$orientation <- as.numeric(members$orientation)
if (length(unique(members$signature_id)) != 18L || any(members$gene == "TYMP")) stop("Frozen member gate failed")

# E-MTAB-3218
emtab_probe <- as.matrix(read.delim(emtab_rma_path, header = TRUE, row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE))
storage.mode(emtab_probe) <- "double"
emtab_annotation <- read.delim(emtab_annotation_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
annotation_index <- match(rownames(emtab_probe), emtab_annotation$Probeset)
emtab_gene <- collapse_gene_rows(emtab_probe, emtab_annotation$Symbol[annotation_index])
rm(emtab_probe)
emtab_samples <- read.delim(emtab_samples_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
emtab_samples$matrix_column <- paste0("ENTREZG-", sub("\\.cel$", "", emtab_samples$array_file, ignore.case = TRUE))
if (anyNA(match(emtab_samples$matrix_column, colnames(emtab_gene)))) stop("E-MTAB sample mapping failed")
emtab_gene <- emtab_gene[, emtab_samples$matrix_column, drop = FALSE]
emtab_scored <- score_core(emtab_gene, rep("E-MTAB-3218", ncol(emtab_gene)), "E-MTAB-3218", members)
emtab <- cbind(emtab_samples, emtab_scored$composite[match(emtab_samples$matrix_column, emtab_scored$composite$sample_id), -1L, drop = FALSE])
emtab$dataset <- "E-MTAB-3218"
emtab$patient_key <- as.character(emtab$patient_id)
emtab$cohort_group <- "E-MTAB-3218"
emtab$arm_group <- as.character(emtab$trial_arm)
emtab$orr_primary <- suppressWarnings(as.numeric(emtab$orr_response))
emtab$response20 <- ifelse(grepl("^20pct", tolower(emtab$response20pct)), 1, ifelse(tolower(emtab$response20pct) == "not20pct", 0, NA))

# CheckMate
checkmate_raw <- read.delim(gzfile(checkmate_expr_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
checkmate_genes <- checkmate_raw[[1L]]
checkmate_matrix <- as.matrix(checkmate_raw[, -1L, drop = FALSE])
storage.mode(checkmate_matrix) <- "double"
rm(checkmate_raw)
checkmate_gene <- collapse_gene_rows(checkmate_matrix, checkmate_genes)
rm(checkmate_matrix)
checkmate_clinical <- read.delim(checkmate_clinical_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
checkmate_index <- match(colnames(checkmate_gene), checkmate_clinical$RNA_ID)
if (anyNA(checkmate_index)) stop("CheckMate RNA mapping failed")
checkmate <- checkmate_clinical[checkmate_index, , drop = FALSE]
checkmate_scored <- score_core(checkmate_gene, checkmate$Cohort, "Braun2020_CheckMate", members)
checkmate <- cbind(checkmate, checkmate_scored$composite[match(checkmate$RNA_ID, checkmate_scored$composite$sample_id), -1L, drop = FALSE])
checkmate$dataset <- "Braun2020_CheckMate"
checkmate$patient_key <- as.character(checkmate$SUBJID)
checkmate$cohort_group <- as.character(checkmate$Cohort)
checkmate$arm_group <- as.character(checkmate$Arm)
checkmate$orr_binary <- ifelse(checkmate$ORR %in% c("CR", "PR", "CRPR"), 1, ifelse(checkmate$ORR %in% c("SD", "PD"), 0, NA))
checkmate$benefit_binary <- ifelse(checkmate$Benefit == "CB", 1, ifelse(checkmate$Benefit == "NCB", 0, NA))
checkmate$age_z <- z_vector(checkmate$Age)

# JAVELIN
javelin_raw <- read.delim(gzfile(javelin_expr_path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
javelin_genes <- javelin_raw[[1L]]
javelin_matrix <- as.matrix(javelin_raw[, -1L, drop = FALSE])
storage.mode(javelin_matrix) <- "double"
rm(javelin_raw)
javelin_gene <- collapse_gene_rows(javelin_matrix, javelin_genes)
rm(javelin_matrix)
javelin_clinical <- read.delim(javelin_clinical_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
javelin_index <- match(colnames(javelin_gene), javelin_clinical$ID)
if (anyNA(javelin_index)) stop("JAVELIN expression mapping failed")
javelin <- javelin_clinical[javelin_index, , drop = FALSE]
javelin_scored <- score_core(javelin_gene, rep("JAVELIN_Renal101", ncol(javelin_gene)), "JAVELIN_Renal101", members)
javelin <- cbind(javelin, javelin_scored$composite[match(javelin$ID, javelin_scored$composite$sample_id), -1L, drop = FALSE])
javelin$dataset <- "JAVELIN_Renal101"
javelin$patient_key <- as.character(javelin$ID)
javelin$cohort_group <- "JAVELIN_Renal101"
javelin$arm_group <- as.character(javelin$TRT01P)
javelin$age_z <- z_vector(javelin$AGE)

coverage <- rbind(emtab_scored$coverage, checkmate_scored$coverage, javelin_scored$coverage)
signature_scores <- rbind(emtab_scored$signature_long, checkmate_scored$signature_long, javelin_scored$signature_long)
write_tsv(coverage, file.path(output_dir, "R2_SIGNATURE_COVERAGE_BY_COHORT_v1.tsv"))
write_tsv(signature_scores, file.path(output_dir, "R2_SIGNATURE_SCORES_LONG_v1.tsv.gz"))

score_columns <- c("r2_metabolic_core", "r2_regulon_core", "r2_combined_core", "r2_combined_rank", "tymp_expression_z")
sample_scores <- rbind(
  data.frame(dataset = emtab$dataset, sample_id = emtab$sample_id, patient_id = emtab$patient_key, timepoint = emtab$timepoint, cohort = emtab$cohort_group, arm = emtab$arm_group, emtab[, score_columns], stringsAsFactors = FALSE),
  data.frame(dataset = checkmate$dataset, sample_id = checkmate$RNA_ID, patient_id = checkmate$patient_key, timepoint = "baseline", cohort = checkmate$cohort_group, arm = checkmate$arm_group, checkmate[, score_columns], stringsAsFactors = FALSE),
  data.frame(dataset = javelin$dataset, sample_id = javelin$ID, patient_id = javelin$patient_key, timepoint = "baseline", cohort = javelin$cohort_group, arm = javelin$arm_group, javelin[, score_columns], stringsAsFactors = FALSE)
)
write_tsv(sample_scores, file.path(output_dir, "R2_COMPOSITE_SAMPLE_SCORES_v1.tsv"))

features <- c("r2_combined_core", "r2_metabolic_core", "r2_regulon_core")
primary_feature <- "r2_combined_core"
results <- list()
append_result <- function(x) if (nrow(x)) results[[length(results) + 1L]] <<- x

run_seven_models <- function(feature, tier) {
  emtab_baseline <- emtab[emtab$timepoint == "baseline", , drop = FALSE]
  nivo <- checkmate[checkmate$Arm == "NIVOLUMAB", , drop = FALSE]
  cm025 <- checkmate[checkmate$Cohort == "CM-025", , drop = FALSE]
  combo <- javelin[javelin$TRT01P == "Avelumab+Axitinib", , drop = FALSE]
  append_result(fit_logistic(emtab_baseline, feature, "orr_primary", "E-MTAB-3218", "emtab_baseline_orr", tier, "baseline"))
  append_result(fit_cox(emtab_baseline, feature, "pfs_days", "pfs_event", FALSE, "E-MTAB-3218", "emtab_baseline_pfs", tier, "PFS", "baseline"))
  append_result(fit_logistic(nivo, feature, "orr_binary", "Braun2020_CheckMate", "checkmate_nivolumab_orr_adjusted", tier, "nivolumab", c("Cohort", "MSKCC"), c("factor(Cohort)", "factor(MSKCC)"), "cohort+MSKCC"))
  append_result(fit_cox(nivo, feature, "PFS", "PFS_CNSR", TRUE, "Braun2020_CheckMate", "checkmate_nivolumab_pfs_adjusted", tier, "PFS", "nivolumab", c("Cohort", "MSKCC"), c("factor(Cohort)", "factor(MSKCC)"), "cohort+MSKCC"))
  append_result(fit_cox_interaction(cm025, feature, "PFS", "PFS_CNSR", TRUE, cm025$Arm == "NIVOLUMAB", "Braun2020_CheckMate", "checkmate_cm025_pfs_interaction", tier, "PFS", adjustment = "randomized_unadjusted", positive_direction = "higher_more_everolimus_benefit", negative_direction = "higher_more_nivolumab_benefit"))
  append_result(fit_cox_interaction(javelin, feature, "PFS_P", "PFS_P_CNSR", TRUE, javelin$TRT01P == "Avelumab+Axitinib", "JAVELIN_Renal101", "javelin_pfs_interaction_adjusted", tier, "PFS", c("AGE", "SEX"), c("age_z", "factor(SEX)"), "age+sex", positive_direction = "higher_more_sunitinib_benefit", negative_direction = "higher_more_combo_benefit"))
  append_result(fit_cox(combo, feature, "PFS_P", "PFS_P_CNSR", TRUE, "JAVELIN_Renal101", "javelin_combo_pfs_adjusted", tier, "PFS", "Avelumab+Axitinib", c("AGE", "SEX"), c("age_z", "factor(SEX)"), "age+sex"))
}

run_seven_models(primary_feature, "PRIMARY")
for (feature in setdiff(features, primary_feature)) run_seven_models(feature, "SECONDARY_FAMILY")

# Prespecified exploratory models
for (feature in features) {
  emtab_baseline <- emtab[emtab$timepoint == "baseline", , drop = FALSE]
  emtab_week4 <- emtab[emtab$timepoint == "week4", , drop = FALSE]
  append_result(fit_logistic(emtab_baseline, feature, "response20", "E-MTAB-3218", "emtab_baseline_response20", "EXPLORATORY", "baseline"))
  append_result(fit_logistic(emtab_week4, feature, "orr_primary", "E-MTAB-3218", "emtab_week4_orr", "EXPLORATORY", "week4"))
  append_result(fit_logistic(emtab_week4, feature, "response20", "E-MTAB-3218", "emtab_week4_response20", "EXPLORATORY", "week4"))
  append_result(fit_cox(emtab_baseline, feature, "os_days", "os_event", FALSE, "E-MTAB-3218", "emtab_baseline_os", "EXPLORATORY", "OS", "baseline"))

  paired <- merge(emtab_baseline, emtab_week4, by = "patient_id", suffixes = c("_baseline", "_week4"))
  paired[[paste0(feature, "_delta")]] <- paired[[paste0(feature, "_week4")]] - paired[[paste0(feature, "_baseline")]]
  append_result(fit_logistic(paired, paste0(feature, "_delta"), "orr_primary_baseline", "E-MTAB-3218", "emtab_delta_orr", "EXPLORATORY", "paired"))
  append_result(fit_logistic(paired, paste0(feature, "_delta"), "response20_baseline", "E-MTAB-3218", "emtab_delta_response20", "EXPLORATORY", "paired"))

  nivo <- checkmate[checkmate$Arm == "NIVOLUMAB", , drop = FALSE]
  cm025 <- checkmate[checkmate$Cohort == "CM-025", , drop = FALSE]
  append_result(fit_logistic(nivo, feature, "orr_binary", "Braun2020_CheckMate", "checkmate_nivolumab_orr_unadjusted", "EXPLORATORY", "nivolumab"))
  append_result(fit_logistic(nivo, feature, "benefit_binary", "Braun2020_CheckMate", "checkmate_nivolumab_benefit_adjusted", "EXPLORATORY", "nivolumab", c("Cohort", "MSKCC"), c("factor(Cohort)", "factor(MSKCC)"), "cohort+MSKCC"))
  append_result(fit_cox(nivo, feature, "PFS", "PFS_CNSR", TRUE, "Braun2020_CheckMate", "checkmate_nivolumab_pfs_unadjusted", "EXPLORATORY", "PFS", "nivolumab"))
  append_result(fit_cox(nivo, feature, "OS", "OS_CNSR", TRUE, "Braun2020_CheckMate", "checkmate_nivolumab_os_adjusted", "EXPLORATORY", "OS", "nivolumab", c("Cohort", "MSKCC"), c("factor(Cohort)", "factor(MSKCC)"), "cohort+MSKCC"))
  append_result(fit_cox_interaction(cm025, feature, "OS", "OS_CNSR", TRUE, cm025$Arm == "NIVOLUMAB", "Braun2020_CheckMate", "checkmate_cm025_os_interaction", "EXPLORATORY", "OS", adjustment = "randomized_unadjusted", positive_direction = "higher_more_everolimus_benefit", negative_direction = "higher_more_nivolumab_benefit"))

  combo <- javelin[javelin$TRT01P == "Avelumab+Axitinib", , drop = FALSE]
  sunitinib <- javelin[javelin$TRT01P == "Sunitinib", , drop = FALSE]
  append_result(fit_cox_interaction(javelin, feature, "PFS_P", "PFS_P_CNSR", TRUE, javelin$TRT01P == "Avelumab+Axitinib", "JAVELIN_Renal101", "javelin_pfs_interaction_unadjusted", "EXPLORATORY", "PFS", adjustment = "randomized_unadjusted", positive_direction = "higher_more_sunitinib_benefit", negative_direction = "higher_more_combo_benefit"))
  append_result(fit_cox(combo, feature, "PFS_P", "PFS_P_CNSR", TRUE, "JAVELIN_Renal101", "javelin_combo_pfs_unadjusted", "EXPLORATORY", "PFS", "Avelumab+Axitinib"))
  append_result(fit_cox(sunitinib, feature, "PFS_P", "PFS_P_CNSR", TRUE, "JAVELIN_Renal101", "javelin_sunitinib_pfs_adjusted", "EXPLORATORY", "PFS", "Sunitinib", c("AGE", "SEX"), c("age_z", "factor(SEX)"), "age+sex"))
}

# Rank-composite sensitivity repeats the seven model skeletons as exploratory.
run_seven_models("r2_combined_rank", "EXPLORATORY")

models <- if (length(results)) do.call(rbind, results) else empty_result()
for (tier in unique(models$analysis_tier)) {
  idx <- which(models$analysis_tier == tier & is.finite(models$p_value))
  models$q_value[idx] <- p.adjust(models$p_value[idx], method = "BH")
}
models$significance_class <- ifelse(
  is.finite(models$q_value) & models$q_value < 0.05, "FDR_LT_0.05",
  ifelse(is.finite(models$p_value) & models$p_value < 0.05, "NOMINAL_P_LT_0.05", "NOT_SIGNIFICANT")
)
models <- models[order(match(models$analysis_tier, c("PRIMARY", "SECONDARY_FAMILY", "EXPLORATORY")), models$p_value, models$cohort, models$model_id), , drop = FALSE]
write_tsv(models, file.path(output_dir, "R2_COMPOSITE_CLINICAL_MODELS_ALL_v1.tsv"))
write_tsv(models[models$analysis_tier == "PRIMARY", , drop = FALSE], file.path(output_dir, "R2_COMPOSITE_PRIMARY_MODELS_v1.tsv"))
write_tsv(models[models$analysis_tier == "SECONDARY_FAMILY", , drop = FALSE], file.path(output_dir, "R2_COMPOSITE_SECONDARY_MODELS_v1.tsv"))
write_tsv(models[models$analysis_tier == "EXPLORATORY", , drop = FALSE], file.path(output_dir, "R2_COMPOSITE_EXPLORATORY_MODELS_v1.tsv"))

primary <- models[models$analysis_tier == "PRIMARY", , drop = FALSE]
secondary <- models[models$analysis_tier == "SECONDARY_FAMILY", , drop = FALSE]
convergence <- data.frame(
  metric = c(
    "primary_models", "primary_fdr", "primary_nominal", "primary_cohorts_with_fdr",
    "secondary_models", "secondary_fdr", "exploratory_fdr", "clinical_upgrade_rule_met"
  ),
  value = c(
    nrow(primary), sum(primary$significance_class == "FDR_LT_0.05"),
    sum(is.finite(primary$p_value) & primary$p_value < 0.05),
    length(unique(primary$cohort[primary$significance_class == "FDR_LT_0.05"])),
    nrow(secondary), sum(secondary$significance_class == "FDR_LT_0.05"),
    sum(models$analysis_tier == "EXPLORATORY" & models$significance_class == "FDR_LT_0.05"),
    length(unique(primary$cohort[primary$significance_class == "FDR_LT_0.05"])) >= 2
  ),
  stringsAsFactors = FALSE
)
write_tsv(convergence, file.path(output_dir, "R2_COMPOSITE_CLINICAL_CONVERGENCE_v1.tsv"))

coverage_pass <- aggregate(coverage_pass ~ dataset + family, coverage, sum)
write_tsv(coverage_pass, file.path(output_dir, "R2_SIGNATURE_FAMILY_COVERAGE_SUMMARY_v1.tsv"))

primary_lines <- apply(primary, 1L, function(row) {
  sprintf("- %s / %s: %s=%.4g, 95%%CI %.4g–%.4g, P=%.4g, q=%.4g, %s.",
    row[["cohort"]], row[["model_id"]], row[["effect_type"]], as.numeric(row[["effect"]]),
    as.numeric(row[["ci_low"]]), as.numeric(row[["ci_high"]]), as.numeric(row[["p_value"]]),
    as.numeric(row[["q_value"]]), row[["direction"]])
})
upgrade <- convergence$value[convergence$metric == "clinical_upgrade_rule_met"] == "TRUE"
summary_text <- c(
  "# Phase69：R2髓系核心复合分数独立临床验证",
  "",
  "## 设计",
  "",
  "18个Phase67 R2 Myeloid/DC–TYMP签名被冻结为9个代谢签名和9个调控签名，并以等家族权重形成连续复合分数。TYMP不在任何签名成员中。临床验证仅使用独立的E-MTAB-3218、Braun2020 CheckMate和JAVELIN Renal 101。",
  "",
  "## 主模型",
  "",
  primary_lines,
  "",
  "## 判定",
  "",
  sprintf("七个主模型中FDR结果为%d项，涉及%d个独立临床队列。预注册的临床升级规则%s。",
    sum(primary$significance_class == "FDR_LT_0.05"),
    length(unique(primary$cohort[primary$significance_class == "FDR_LT_0.05"])),
    ifelse(upgrade, "达到", "未达到")),
  ifelse(upgrade,
    "复合分数获得跨队列临床支持；仍需按治疗背景解释，不等同于TYMP单基因预测。",
    "R2核心的机制复现不自动转化为response/PFS预测；F1/TYMP直接临床预测等级不升级。"),
  "",
  "## 边界",
  "",
  "bulk复合分数是组织级机制投影，不是纯Myeloid/DC细胞内读数。次级家族、治疗后变化、OS、benefit和未调整模型不能替代七个预注册主检验。"
)
writeLines(summary_text, file.path(output_dir, "PHASE69_R2_COMPOSITE_CLINICAL_SUMMARY_v1.md"), useBytes = TRUE)

status <- data.frame(
  gate = c(
    "frozen_signatures", "tymp_members", "emtab_samples", "checkmate_samples", "javelin_samples",
    "coverage_rows", "family_coverage_all_ge7", "primary_models", "secondary_models", "survival_available"
  ),
  observed = c(
    length(unique(members$signature_id)), sum(members$gene == "TYMP"), nrow(emtab), nrow(checkmate), nrow(javelin),
    nrow(coverage), all(coverage_pass$coverage_pass >= 7), nrow(primary), nrow(secondary), requireNamespace("survival", quietly = TRUE)
  ),
  expected = c(18, 0, 114, 311, 726, 54, TRUE, 7, 14, TRUE),
  status = c(
    ifelse(length(unique(members$signature_id)) == 18, "PASS", "FAIL"),
    ifelse(sum(members$gene == "TYMP") == 0, "PASS", "FAIL"),
    ifelse(nrow(emtab) == 114, "PASS", "FAIL"), ifelse(nrow(checkmate) == 311, "PASS", "FAIL"),
    ifelse(nrow(javelin) == 726, "PASS", "FAIL"), ifelse(nrow(coverage) == 54, "PASS", "FAIL"),
    ifelse(all(coverage_pass$coverage_pass >= 7), "PASS", "FAIL"), ifelse(nrow(primary) == 7, "PASS", "FAIL"),
    ifelse(nrow(secondary) == 14, "PASS", "FAIL"), ifelse(requireNamespace("survival", quietly = TRUE), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "PHASE69_RUN_STATUS_v1.tsv"))
print(status, row.names = FALSE)
print(primary[, c("cohort", "model_id", "effect", "ci_low", "ci_high", "p_value", "q_value", "significance_class")], row.names = FALSE)
if (any(status$status == "FAIL")) quit(status = 2L)
cat("PASS_PHASE69_R2_COMPOSITE_CLINICAL_VALIDATION\n")
