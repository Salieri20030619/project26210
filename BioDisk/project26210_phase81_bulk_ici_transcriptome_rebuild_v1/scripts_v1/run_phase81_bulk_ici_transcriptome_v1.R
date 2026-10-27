#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9L) {
  stop("Usage: run_phase81_bulk_ici_transcriptome_v1.R members checkmate_expr checkmate_clinical checkmate_ciber javelin_expr javelin_clinical javelin_deconv prereg outdir")
}
members_path <- args[[1L]]
checkmate_expr_path <- args[[2L]]
checkmate_clinical_path <- args[[3L]]
checkmate_ciber_path <- args[[4L]]
javelin_expr_path <- args[[5L]]
javelin_clinical_path <- args[[6L]]
javelin_deconv_path <- args[[7L]]
prereg_path <- args[[8L]]
outdir <- args[[9L]]
part <- paste0(outdir, ".part")
if (dir.exists(outdir) || file.exists(outdir) || dir.exists(part) || file.exists(part)) stop("Refusing existing output")
dir.create(part, recursive=TRUE)

z_vector <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm=TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm=TRUE)) / s
}

sha256_file <- function(path) {
  unname(tools::md5sum(path))
}

read_expression <- function(path) {
  table <- fread(path, check.names=FALSE)
  gene_column <- names(table)[[1L]]
  genes <- toupper(as.character(table[[gene_column]]))
  matrix <- as.matrix(table[, -1L, with=FALSE])
  storage.mode(matrix) <- "double"
  if (anyDuplicated(genes)) {
    matrix <- rowsum(matrix, group=genes, reorder=TRUE)
    counts <- table(genes)
    matrix <- matrix / as.numeric(counts[rownames(matrix)])
    genes <- rownames(matrix)
  }
  rownames(matrix) <- genes
  matrix
}

module_score <- function(expression, genes, label) {
  present <- intersect(genes, rownames(expression))
  if (length(present) != length(genes)) {
    stop(label, " missing genes: ", paste(setdiff(genes, present), collapse=";"))
  }
  matrix <- expression[present, , drop=FALSE]
  zmat <- t(apply(matrix, 1L, z_vector))
  colMeans(zmat)
}

score_expression <- function(expression, members, dataset) {
  modules <- split(members$gene, members$signature_id)
  required <- c("MTA_RESPONSE", "SPP1_MACROPHAGE_NO_SPP1", "APOE_TAM", "POSTN_CAF", "TYMP")
  if (!all(required %in% names(modules))) stop("Missing frozen modules")
  scores <- data.table(sample_id=colnames(expression), dataset=dataset)
  scores[, tymp_z := z_vector(as.numeric(expression["TYMP", ]))]
  scores[, mta_response := z_vector(module_score(expression, modules$MTA_RESPONSE, "MTA_RESPONSE"))]
  scores[, spp1_macrophage := z_vector(module_score(expression, modules$SPP1_MACROPHAGE_NO_SPP1, "SPP1"))]
  scores[, apoe_tam := z_vector(module_score(expression, modules$APOE_TAM, "APOE_TAM"))]
  scores[, postn_caf := z_vector(module_score(expression, modules$POSTN_CAF, "POSTN_CAF"))]
  scores[, full_chain_proxy := z_vector((tymp_z + mta_response + spp1_macrophage) / 3)]
  scores[, tymp_spp1_joint := z_vector((tymp_z + spp1_macrophage) / 2)]
  scores[, tib_paper_proxy := z_vector((apoe_tam + postn_caf) / 2)]
  scores[, extended_chain_tib := z_vector((full_chain_proxy + tib_paper_proxy) / 2)]
  scores[, tymp_x_spp1 := z_vector(tymp_z * spp1_macrophage)]
  scores
}

extract_checkmate_macro <- function(path, sample_ids) {
  table <- fread(path, check.names=FALSE)
  feature <- tolower(as.character(table[[1L]]))
  matrix <- as.matrix(table[, -1L, with=FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- feature
  targets <- c("macrophages_m0", "macrophages_m1", "macrophages_m2")
  if (!all(targets %in% rownames(matrix))) stop("CheckMate macrophage rows missing")
  values <- colSums(matrix[targets, sample_ids, drop=FALSE])
  data.table(sample_id=sample_ids, macro_abundance=z_vector(values))
}

extract_javelin_macro <- function(path) {
  table <- fread(path, check.names=FALSE)
  targets <- c("Macrophages_M0", "Macrophages_M1", "Macrophages_M2")
  if (!all(targets %in% names(table))) stop("JAVELIN macrophage columns missing")
  data.table(
    sample_id=as.character(table$row_id),
    macro_abundance=z_vector(rowSums(as.matrix(table[, ..targets])))
  )
}

empty_model <- function() {
  data.table(
    cohort=character(), model_id=character(), analysis_tier=character(),
    endpoint=character(), stratum=character(), feature=character(),
    n=integer(), n_event=integer(), effect_type=character(), effect=numeric(),
    ci_low=numeric(), ci_high=numeric(), p_value=numeric(), ph_p_value=numeric(),
    direction=character(), adjustment=character(), notes=character()
  )
}

model_row <- function(cohort, model_id, tier, endpoint, stratum, feature, n,
                      n_event, effect_type, coefficient, se, p, ph_p=NA_real_,
                      direction, adjustment, notes="") {
  data.table(
    cohort=cohort, model_id=model_id, analysis_tier=tier, endpoint=endpoint,
    stratum=stratum, feature=feature, n=as.integer(n), n_event=as.integer(n_event),
    effect_type=effect_type, effect=exp(coefficient),
    ci_low=exp(coefficient - 1.96 * se), ci_high=exp(coefficient + 1.96 * se),
    p_value=p, ph_p_value=ph_p, direction=direction,
    adjustment=adjustment, notes=notes
  )
}

cox_checkmate <- function(frame, feature, endpoint="PFS", tier="PRIMARY",
                          composition=FALSE, interaction=FALSE) {
  time_name <- ifelse(endpoint == "PFS", "PFS", "OS")
  event_name <- ifelse(endpoint == "PFS", "PFS_CNSR", "OS_CNSR")
  data <- copy(frame)
  data[, score := z_vector(get(feature))]
  data[, time := as.numeric(get(time_name))]
  data[, event := as.numeric(get(event_name))]
  if (interaction) {
    data <- data[Cohort == "CM-025"]
    data[, arm_nivo := as.integer(Arm == "NIVOLUMAB")]
    keep <- is.finite(data$score) & is.finite(data$time) & data$time > 0 &
      data$event %in% c(0, 1) & data$arm_nivo %in% c(0, 1) & nzchar(data$MSKCC)
    data <- data[keep]
    formula <- if (composition) {
      Surv(time, event) ~ score * arm_nivo + factor(MSKCC) + macro_abundance
    } else {
      Surv(time, event) ~ score * arm_nivo + factor(MSKCC)
    }
    fit <- coxph(formula, data=data, ties="efron")
    coefficient_name <- "score:arm_nivo"
    table <- summary(fit)$coef
    coefficient <- table[coefficient_name, "coef"]
    se <- table[coefficient_name, "se(coef)"]
    p <- table[coefficient_name, "Pr(>|z|)"]
    return(model_row(
      "Braun2020_CheckMate", paste0("checkmate_", tolower(endpoint), "_treatment_interaction"),
      tier, endpoint, "CM-025_randomized", feature, nrow(data), sum(data$event),
      "interaction_hazard_ratio", coefficient, se, p, NA_real_,
      ifelse(coefficient > 0, "higher_more_harmful_with_nivolumab", "higher_more_benefit_with_nivolumab"),
      ifelse(composition, "MSKCC+macrophage", "MSKCC")
    ))
  }
  data <- data[Arm == "NIVOLUMAB"]
  keep <- is.finite(data$score) & is.finite(data$time) & data$time > 0 &
    data$event %in% c(0, 1) & nzchar(data$Cohort) & nzchar(data$MSKCC)
  if (composition) keep <- keep & is.finite(data$macro_abundance)
  data <- data[keep]
  formula <- if (composition) {
    Surv(time, event) ~ score + factor(Cohort) + factor(MSKCC) + macro_abundance
  } else {
    Surv(time, event) ~ score + factor(Cohort) + factor(MSKCC)
  }
  fit <- coxph(formula, data=data, ties="efron")
  table <- summary(fit)$coef
  coefficient <- table["score", "coef"]
  se <- table["score", "se(coef)"]
  p <- table["score", "Pr(>|z|)"]
  ph <- try(cox.zph(fit), silent=TRUE)
  ph_p <- if (inherits(ph, "try-error")) NA_real_ else ph$table["score", "p"]
  model_row(
    "Braun2020_CheckMate", paste0("checkmate_nivolumab_", tolower(endpoint)),
    tier, endpoint, "nivolumab", feature, nrow(data), sum(data$event),
    "hazard_ratio_per_sd", coefficient, se, p, ph_p,
    ifelse(coefficient > 0, "higher_worse_survival", "higher_better_survival"),
    ifelse(composition, "cohort+MSKCC+macrophage", "cohort+MSKCC")
  )
}

cox_javelin <- function(frame, feature, tier="PRIMARY", composition=FALSE,
                        interaction=FALSE) {
  data <- copy(frame)
  data[, score := z_vector(get(feature))]
  data[, time := as.numeric(PFS_P)]
  data[, event := 1 - as.numeric(PFS_P_CNSR)]
  data[, combo := as.integer(TRT01P == "Avelumab+Axitinib")]
  if (interaction) {
    keep <- is.finite(data$score) & is.finite(data$time) & data$time > 0 &
      data$event %in% c(0, 1) & data$combo %in% c(0, 1) &
      is.finite(data$AGE) & nzchar(data$SEX)
    if (composition) keep <- keep & is.finite(data$macro_abundance)
    data <- data[keep]
    formula <- if (composition) {
      Surv(time, event) ~ score * combo + AGE + factor(SEX) + macro_abundance
    } else {
      Surv(time, event) ~ score * combo + AGE + factor(SEX)
    }
    fit <- coxph(formula, data=data, ties="efron")
    table <- summary(fit)$coef
    coefficient <- table["score:combo", "coef"]
    se <- table["score:combo", "se(coef)"]
    p <- table["score:combo", "Pr(>|z|)"]
    return(model_row(
      "JAVELIN_Renal101", "javelin_pfs_treatment_interaction", tier, "PFS",
      "randomized_all", feature, nrow(data), sum(data$event),
      "interaction_hazard_ratio", coefficient, se, p, NA_real_,
      ifelse(coefficient > 0, "higher_more_harmful_with_combo", "higher_more_benefit_with_combo"),
      ifelse(composition, "age+sex+macrophage", "age+sex")
    ))
  }
  data <- data[TRT01P == "Avelumab+Axitinib"]
  keep <- is.finite(data$score) & is.finite(data$time) & data$time > 0 &
    data$event %in% c(0, 1) & is.finite(data$AGE) & nzchar(data$SEX)
  if (composition) keep <- keep & is.finite(data$macro_abundance)
  data <- data[keep]
  formula <- if (composition) {
    Surv(time, event) ~ score + AGE + factor(SEX) + macro_abundance
  } else {
    Surv(time, event) ~ score + AGE + factor(SEX)
  }
  fit <- coxph(formula, data=data, ties="efron")
  table <- summary(fit)$coef
  coefficient <- table["score", "coef"]
  se <- table["score", "se(coef)"]
  p <- table["score", "Pr(>|z|)"]
  ph <- try(cox.zph(fit), silent=TRUE)
  ph_p <- if (inherits(ph, "try-error")) NA_real_ else ph$table["score", "p"]
  model_row(
    "JAVELIN_Renal101", "javelin_combo_pfs", tier, "PFS",
    "Avelumab+Axitinib", feature, nrow(data), sum(data$event),
    "hazard_ratio_per_sd", coefficient, se, p, ph_p,
    ifelse(coefficient > 0, "higher_worse_survival", "higher_better_survival"),
    ifelse(composition, "age+sex+macrophage", "age+sex")
  )
}

logistic_checkmate <- function(frame, feature, endpoint, tier) {
  data <- copy(frame[Arm == "NIVOLUMAB"])
  data[, score := z_vector(get(feature))]
  if (endpoint == "ORR") {
    data[, outcome := fifelse(ORR %in% c("CR", "PR", "CRPR"), 1,
      fifelse(ORR %in% c("SD", "PD"), 0, NA_real_))]
  } else {
    data[, outcome := fifelse(Benefit == "CB", 1,
      fifelse(Benefit == "NCB", 0, NA_real_))]
  }
  keep <- is.finite(data$score) & data$outcome %in% c(0, 1) &
    nzchar(data$Cohort) & nzchar(data$MSKCC)
  data <- data[keep]
  fit <- glm(outcome ~ score + factor(Cohort) + factor(MSKCC), data=data, family=binomial())
  table <- coef(summary(fit))
  coefficient <- table["score", "Estimate"]
  se <- table["score", "Std. Error"]
  p <- table["score", "Pr(>|z|)"]
  model_row(
    "Braun2020_CheckMate", paste0("checkmate_nivolumab_", tolower(endpoint)),
    tier, endpoint, "nivolumab", feature, nrow(data), sum(data$outcome),
    "odds_ratio_per_sd", coefficient, se, p, NA_real_,
    ifelse(coefficient > 0, "higher_in_response", "lower_in_response"),
    "cohort+MSKCC"
  )
}

members <- fread(members_path)
members[, gene := toupper(gene)]
if (nrow(members) != 30L || uniqueN(members$gene) != 30L) stop("Frozen signature membership failed")

checkmate_expression <- read_expression(checkmate_expr_path)
javelin_expression <- read_expression(javelin_expr_path)
required_genes <- unique(members$gene)
if (!all(required_genes %in% rownames(checkmate_expression))) {
  stop("CheckMate missing: ", paste(setdiff(required_genes, rownames(checkmate_expression)), collapse=";"))
}
if (!all(required_genes %in% rownames(javelin_expression))) {
  stop("JAVELIN missing: ", paste(setdiff(required_genes, rownames(javelin_expression)), collapse=";"))
}

checkmate_scores <- score_expression(checkmate_expression, members, "Braun2020_CheckMate")
javelin_scores <- score_expression(javelin_expression, members, "JAVELIN_Renal101")
checkmate_macro <- extract_checkmate_macro(checkmate_ciber_path, checkmate_scores$sample_id)
javelin_macro <- extract_javelin_macro(javelin_deconv_path)
checkmate_scores <- merge(checkmate_scores, checkmate_macro, by="sample_id", all.x=TRUE, sort=FALSE)
javelin_scores <- merge(javelin_scores, javelin_macro, by="sample_id", all.x=TRUE, sort=FALSE)

checkmate_clinical <- fread(checkmate_clinical_path, check.names=FALSE)
checkmate <- merge(
  checkmate_scores, checkmate_clinical,
  by.x="sample_id", by.y="RNA_ID", all.x=TRUE, sort=FALSE
)
if (nrow(checkmate) != 311L || anyNA(checkmate$Arm)) stop("CheckMate clinical merge failed")

javelin_clinical <- fread(javelin_clinical_path, check.names=FALSE)
javelin <- merge(
  javelin_scores, javelin_clinical,
  by.x="sample_id", by.y="ID", all.x=TRUE, sort=FALSE
)
if (nrow(javelin) != 726L || anyNA(javelin$TRT01P)) stop("JAVELIN clinical merge failed")

score_columns <- c(
  "tymp_z", "mta_response", "spp1_macrophage", "full_chain_proxy",
  "tymp_spp1_joint", "apoe_tam", "postn_caf", "tib_paper_proxy",
  "extended_chain_tib", "tymp_x_spp1", "macro_abundance"
)
sample_scores <- rbindlist(list(
  checkmate[, c(list(sample_id=sample_id, dataset=dataset, arm=Arm), .SD), .SDcols=score_columns],
  javelin[, c(list(sample_id=sample_id, dataset=dataset, arm=TRT01P), .SD), .SDcols=score_columns]
), use.names=TRUE)
fwrite(sample_scores, file.path(part, "phase81_sample_scores_v1.tsv"), sep="\t")

primary_features <- c("full_chain_proxy", "tymp_spp1_joint")
primary_rows <- list()
index <- 0L
for (feature in primary_features) {
  index <- index + 1L
  primary_rows[[index]] <- cox_checkmate(checkmate, feature, "PFS", "PRIMARY", FALSE, FALSE)
  index <- index + 1L
  primary_rows[[index]] <- cox_javelin(javelin, feature, "PRIMARY", FALSE, FALSE)
}
primary <- rbindlist(primary_rows)
primary[, q_value := p.adjust(p_value, method="BH")]
primary[, significance_class := fifelse(
  q_value < 0.05 & effect > 1, "COHORT_SUPPORTED",
  fifelse(p_value < 0.05 & effect > 1, "NOMINAL_EXPECTED_DIRECTION",
    fifelse(p_value < 0.05 & effect < 1, "NOMINAL_DIRECTION_CONFLICT", "NOT_SIGNIFICANT"))
)]

secondary_features <- c(
  "tymp_z", "mta_response", "spp1_macrophage", "apoe_tam", "postn_caf",
  "tib_paper_proxy", "extended_chain_tib", "tymp_x_spp1"
)
secondary_rows <- list()
index <- 0L
for (feature in secondary_features) {
  index <- index + 1L
  secondary_rows[[index]] <- cox_checkmate(checkmate, feature, "PFS", "SECONDARY", FALSE, FALSE)
  index <- index + 1L
  secondary_rows[[index]] <- cox_javelin(javelin, feature, "SECONDARY", FALSE, FALSE)
}
for (feature in primary_features) {
  index <- index + 1L
  secondary_rows[[index]] <- cox_checkmate(checkmate, feature, "OS", "SECONDARY", FALSE, FALSE)
  index <- index + 1L
  secondary_rows[[index]] <- logistic_checkmate(checkmate, feature, "ORR", "SECONDARY")
  index <- index + 1L
  secondary_rows[[index]] <- logistic_checkmate(checkmate, feature, "BENEFIT", "SECONDARY")
  index <- index + 1L
  secondary_rows[[index]] <- cox_checkmate(checkmate, feature, "PFS", "SECONDARY_INTERACTION", FALSE, TRUE)
  index <- index + 1L
  secondary_rows[[index]] <- cox_javelin(javelin, feature, "SECONDARY_INTERACTION", FALSE, TRUE)
}
secondary <- rbindlist(secondary_rows)
secondary[, q_value := p.adjust(p_value, method="BH")]
secondary[, significance_class := fifelse(
  q_value < 0.05, "FDR_LT_0.05",
  fifelse(p_value < 0.05, "NOMINAL_P_LT_0.05", "NOT_SIGNIFICANT")
)]

sensitivity_rows <- list()
index <- 0L
for (feature in primary_features) {
  index <- index + 1L
  sensitivity_rows[[index]] <- cox_checkmate(checkmate, feature, "PFS", "COMPOSITION_SENSITIVITY", TRUE, FALSE)
  index <- index + 1L
  sensitivity_rows[[index]] <- cox_javelin(javelin, feature, "COMPOSITION_SENSITIVITY", TRUE, FALSE)
}
sensitivity <- rbindlist(sensitivity_rows)
sensitivity[, q_value := p.adjust(p_value, method="BH")]
sensitivity[, significance_class := fifelse(
  q_value < 0.05, "FDR_LT_0.05",
  fifelse(p_value < 0.05, "NOMINAL_P_LT_0.05", "NOT_SIGNIFICANT")
)]

decision_rows <- list()
for (feature in primary_features) {
  target_feature <- feature
  current <- primary[feature == target_feature]
  checkmate_row <- current[cohort == "Braun2020_CheckMate"]
  javelin_row <- current[cohort == "JAVELIN_Renal101"]
  both_supported <- nrow(checkmate_row) == 1L && nrow(javelin_row) == 1L &&
    checkmate_row$effect > 1 && checkmate_row$q_value < 0.05 &&
    javelin_row$effect > 1 && javelin_row$q_value < 0.05
  one_plus_nominal <- (
    checkmate_row$effect > 1 && checkmate_row$q_value < 0.05 &&
      javelin_row$effect > 1 && javelin_row$p_value < 0.05
  ) || (
    javelin_row$effect > 1 && javelin_row$q_value < 0.05 &&
      checkmate_row$effect > 1 && checkmate_row$p_value < 0.05
  )
  same_direction <- checkmate_row$effect > 1 && javelin_row$effect > 1
  decision_rows[[length(decision_rows) + 1L]] <- data.table(
    feature=feature,
    checkmate_hr=checkmate_row$effect,
    checkmate_p=checkmate_row$p_value,
    checkmate_q=checkmate_row$q_value,
    javelin_hr=javelin_row$effect,
    javelin_p=javelin_row$p_value,
    javelin_q=javelin_row$q_value,
    dual_cohort_supported=both_supported,
    one_cohort_plus_external_nominal=one_plus_nominal,
    directional_consistency=same_direction,
    decision=fifelse(both_supported, "DUAL_COHORT_SUPPORTED",
      fifelse(one_plus_nominal, "ONE_COHORT_PLUS_EXTERNAL_NOMINAL",
        fifelse(same_direction, "DIRECTIONALLY_CONSISTENT_NOT_SUPPORTED", "DIRECTION_CONFLICT")))
  )
}
decisions <- rbindlist(decision_rows)

coverage <- rbindlist(lapply(c("Braun2020_CheckMate", "JAVELIN_Renal101"), function(dataset) {
  expression <- if (dataset == "Braun2020_CheckMate") checkmate_expression else javelin_expression
  members[, .(
    dataset=dataset,
    signature_id=signature_id,
    gene=gene,
    present=gene %in% rownames(expression)
  )]
}))

checks <- data.table(
  check=c(
    "frozen_members_30", "checkmate_samples_311", "javelin_samples_726",
    "checkmate_nivolumab_pfs_events_159", "javelin_combo_pfs_events_157",
    "primary_rows_4", "primary_effects_finite", "all_signature_genes_present"
  ),
  value=c(
    nrow(members), nrow(checkmate), nrow(javelin),
    sum(checkmate[Arm == "NIVOLUMAB"]$PFS_CNSR),
    sum(1 - javelin[TRT01P == "Avelumab+Axitinib"]$PFS_P_CNSR),
    nrow(primary), sum(is.finite(primary$effect)), sum(coverage$present)
  ),
  pass=c(
    nrow(members) == 30L,
    nrow(checkmate) == 311L,
    nrow(javelin) == 726L,
    sum(checkmate[Arm == "NIVOLUMAB"]$PFS_CNSR) == 159L,
    sum(1 - javelin[TRT01P == "Avelumab+Axitinib"]$PFS_P_CNSR) == 157L,
    nrow(primary) == 4L,
    all(is.finite(primary$effect)),
    all(coverage$present)
  )
)

fwrite(primary, file.path(part, "phase81_primary_models_v1.tsv"), sep="\t")
fwrite(secondary, file.path(part, "phase81_secondary_models_v1.tsv"), sep="\t")
fwrite(sensitivity, file.path(part, "phase81_composition_sensitivity_v1.tsv"), sep="\t")
fwrite(decisions, file.path(part, "phase81_cross_cohort_decisions_v1.tsv"), sep="\t")
fwrite(coverage, file.path(part, "phase81_signature_coverage_v1.tsv"), sep="\t")
fwrite(checks, file.path(part, "phase81_validation_checks_v1.tsv"), sep="\t")

summary_lines <- c(
  "# Phase 81 bulk ICI transcriptome summary v1",
  "",
  paste0("- Completed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "- CheckMate PFS event coding repaired to direct PFS_CNSR.",
  "- JAVELIN PFS event coding retained as 1-PFS_P_CNSR.",
  "",
  "## Cross-cohort decisions",
  ""
)
for (i in seq_len(nrow(decisions))) {
  row <- decisions[i]
  summary_lines <- c(summary_lines, sprintf(
    "- %s: CheckMate HR %.4f (q %.4g); JAVELIN HR %.4f (q %.4g); %s.",
    row$feature, row$checkmate_hr, row$checkmate_q,
    row$javelin_hr, row$javelin_q, row$decision
  ))
}
writeLines(summary_lines, file.path(part, "PHASE81_SUMMARY_v1.md"), useBytes=TRUE)
if (!all(checks$pass)) stop("Phase81 validation checks failed")
writeLines("status=PASS", file.path(part, "PHASE81_PASS_v1"))
if (!file.rename(part, outdir)) stop("Atomic output rename failed")
print(primary)
print(decisions)
