#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7L) stop("Usage: script phase81_scores tam.csv caf.csv checkmate_clinical javelin_clinical prereg outdir")
scores_path <- args[[1L]]
tam_path <- args[[2L]]
caf_path <- args[[3L]]
checkmate_clinical_path <- args[[4L]]
javelin_clinical_path <- args[[5L]]
prereg_path <- args[[6L]]
outdir <- args[[7L]]
part <- paste0(outdir, ".part")
if (dir.exists(outdir) || dir.exists(part) || file.exists(outdir) || file.exists(part)) stop("Refusing existing output")
dir.create(part, recursive=TRUE)

z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm=TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm=TRUE)) / s
}

fit_cox <- function(data, feature, cohort) {
  data <- copy(data)
  data[, score := z(get(feature))]
  if (cohort == "CheckMate") {
    data <- data[Arm == "NIVOLUMAB" & is.finite(PFS) & PFS > 0 &
      PFS_CNSR %in% c(0, 1) & nzchar(Cohort) & nzchar(MSKCC)]
    data[, event := as.numeric(PFS_CNSR)]
    data[, time := as.numeric(PFS)]
    fit <- coxph(Surv(time, event) ~ score + factor(Cohort) + factor(MSKCC), data=data)
    cohort_name <- "Braun2020_CheckMate"
    stratum <- "nivolumab"
    adjustment <- "cohort+MSKCC"
  } else {
    data <- data[TRT01P == "Avelumab+Axitinib" & is.finite(PFS_P) & PFS_P > 0 &
      PFS_P_CNSR %in% c(0, 1) & is.finite(AGE) & nzchar(SEX)]
    data[, event := 1 - as.numeric(PFS_P_CNSR)]
    data[, time := as.numeric(PFS_P)]
    fit <- coxph(Surv(time, event) ~ score + AGE + factor(SEX), data=data)
    cohort_name <- "JAVELIN_Renal101"
    stratum <- "Avelumab+Axitinib"
    adjustment <- "age+sex"
  }
  table <- summary(fit)$coef
  beta <- table["score", "coef"]
  se <- table["score", "se(coef)"]
  data.table(
    cohort=cohort_name, stratum=stratum, feature=feature,
    n=nrow(data), n_event=sum(data$event), effect_type="hazard_ratio_per_sd",
    effect=exp(beta), ci_low=exp(beta-1.96*se), ci_high=exp(beta+1.96*se),
    p_value=table["score", "Pr(>|z|)"], adjustment=adjustment
  )
}

scores <- fread(scores_path)
checkmate_scores <- scores[dataset == "Braun2020_CheckMate"]
javelin_scores <- scores[dataset == "JAVELIN_Renal101"]

tam <- fread(tam_path, check.names=FALSE)
setnames(tam, 1L, "sample_id")
caf <- fread(caf_path, check.names=FALSE)
setnames(caf, 1L, "sample_id")
required_tam <- c("C0_APOE_TAM", "C0_VEGFA_TAM", "C0_CD163_TAM")
if (!all(required_tam %in% names(tam)) || !"POSTN_CAF" %in% names(caf)) stop("Paper states missing")
paper <- merge(
  tam[, c("sample_id", required_tam), with=FALSE],
  caf[, .(sample_id, POSTN_CAF)],
  by="sample_id", all=FALSE
)
if (nrow(paper) != 311L) stop("Paper BayesPrism sample count mismatch")

checkmate <- merge(checkmate_scores, paper, by="sample_id", all.x=TRUE, sort=FALSE)
checkmate <- merge(
  checkmate, fread(checkmate_clinical_path),
  by.x="sample_id", by.y="RNA_ID", all.x=TRUE, sort=FALSE
)
if (nrow(checkmate) != 311L || anyNA(checkmate$C0_APOE_TAM)) stop("CheckMate merge failed")
checkmate[, apoe_state_z := z(C0_APOE_TAM)]
checkmate[, postn_state_z := z(POSTN_CAF)]
checkmate[, tymp_apoe_bridge := z((tymp_z + apoe_state_z)/2)]
checkmate[, tymp_postn_bridge := z((tymp_z + postn_state_z)/2)]
checkmate[, tymp_tib_bridge := z((tymp_z + apoe_state_z + postn_state_z)/3)]

javelin <- merge(
  javelin_scores, fread(javelin_clinical_path),
  by.x="sample_id", by.y="ID", all.x=TRUE, sort=FALSE
)
javelin[, apoe_state_z := z(apoe_tam)]
javelin[, postn_state_z := z(postn_caf)]
javelin[, tymp_apoe_bridge := z((tymp_z + apoe_state_z)/2)]
javelin[, tymp_postn_bridge := z((tymp_z + postn_state_z)/2)]
javelin[, tymp_tib_bridge := z((tymp_z + apoe_state_z + postn_state_z)/3)]

features <- c("tymp_apoe_bridge", "tymp_postn_bridge", "tymp_tib_bridge")
rows <- list()
for (feature in features) {
  rows[[length(rows)+1L]] <- fit_cox(checkmate, feature, "CheckMate")
  rows[[length(rows)+1L]] <- fit_cox(javelin, feature, "JAVELIN")
}
primary <- rbindlist(rows)
primary[, q_value := p.adjust(p_value, method="BH")]
primary[, significance_class := fifelse(
  q_value < 0.05 & effect > 1, "COHORT_SUPPORTED",
  fifelse(p_value < 0.05 & effect > 1, "NOMINAL_EXPECTED_DIRECTION",
    fifelse(p_value < 0.05 & effect < 1, "NOMINAL_DIRECTION_CONFLICT", "NOT_SIGNIFICANT"))
)]

secondary <- rbindlist(list(
  fit_cox(checkmate, "apoe_state_z", "CheckMate"),
  fit_cox(checkmate, "postn_state_z", "CheckMate"),
  fit_cox(javelin, "apoe_state_z", "JAVELIN"),
  fit_cox(javelin, "postn_state_z", "JAVELIN")
))
secondary[, q_value := p.adjust(p_value, method="BH")]

correlations <- rbindlist(lapply(
  c("C0_APOE_TAM", "POSTN_CAF", "C0_VEGFA_TAM", "C0_CD163_TAM"),
  function(target) {
    test <- cor.test(checkmate$tymp_z, checkmate[[target]], method="spearman", exact=FALSE)
    data.table(target=target, n=nrow(checkmate), rho=unname(test$estimate), p_value=test$p.value)
  }
))
correlations[, q_value := p.adjust(p_value, method="BH")]

decisions <- rbindlist(lapply(features, function(target_feature) {
  current <- primary[feature == target_feature]
  cm <- current[cohort == "Braun2020_CheckMate"]
  jav <- current[cohort == "JAVELIN_Renal101"]
  dual <- cm$effect > 1 && cm$q_value < 0.05 && jav$effect > 1 && jav$q_value < 0.05
  one_nominal <- (cm$effect > 1 && cm$q_value < 0.05 && jav$effect > 1 && jav$p_value < 0.05) ||
    (jav$effect > 1 && jav$q_value < 0.05 && cm$effect > 1 && cm$p_value < 0.05)
  same <- cm$effect > 1 && jav$effect > 1
  data.table(
    feature=target_feature, checkmate_hr=cm$effect, checkmate_p=cm$p_value,
    checkmate_q=cm$q_value, javelin_hr=jav$effect, javelin_p=jav$p_value,
    javelin_q=jav$q_value, decision=fifelse(dual, "DUAL_COHORT_SUPPORTED",
      fifelse(one_nominal, "ONE_COHORT_PLUS_EXTERNAL_NOMINAL",
        fifelse(same, "DIRECTIONALLY_CONSISTENT_NOT_SUPPORTED", "DIRECTION_CONFLICT")))
  )
}))

checks <- data.table(
  check=c("paper_samples_311", "phase82_primary_rows_6", "checkmate_events_159",
    "javelin_events_157", "all_effects_finite"),
  value=c(nrow(paper), nrow(primary),
    sum(checkmate[Arm=="NIVOLUMAB"]$PFS_CNSR),
    sum(1-javelin[TRT01P=="Avelumab+Axitinib"]$PFS_P_CNSR),
    sum(is.finite(primary$effect))),
  pass=c(nrow(paper)==311L, nrow(primary)==6L,
    sum(checkmate[Arm=="NIVOLUMAB"]$PFS_CNSR)==159L,
    sum(1-javelin[TRT01P=="Avelumab+Axitinib"]$PFS_P_CNSR)==157L,
    all(is.finite(primary$effect)))
)
fwrite(primary, file.path(part, "phase82_primary_models_v1.tsv"), sep="\t")
fwrite(secondary, file.path(part, "phase82_state_secondary_models_v1.tsv"), sep="\t")
fwrite(correlations, file.path(part, "phase82_tymp_state_correlations_v1.tsv"), sep="\t")
fwrite(decisions, file.path(part, "phase82_cross_cohort_decisions_v1.tsv"), sep="\t")
fwrite(checks, file.path(part, "phase82_validation_checks_v1.tsv"), sep="\t")
if (!all(checks$pass)) stop("Phase82 checks failed")
writeLines("status=PASS", file.path(part, "PHASE82_PASS_v1"))
if (!file.rename(part, outdir)) stop("Atomic rename failed")
print(primary)
print(decisions)
print(correlations)
