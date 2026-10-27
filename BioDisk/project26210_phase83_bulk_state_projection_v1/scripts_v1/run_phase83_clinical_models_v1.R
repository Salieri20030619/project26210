#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(survival)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9L) stop("Usage: weights diagnostics phase81_scores cm_clin cm_ciber jav_clin jav_deconv prereg outdir")
weights_path <- args[[1]]; diagnostics_path <- args[[2]]; scores81_path <- args[[3]]
cm_clin_path <- args[[4]]; cm_ciber_path <- args[[5]]; jav_clin_path <- args[[6]]
jav_deconv_path <- args[[7]]; prereg_path <- args[[8]]; outdir <- args[[9]]
part <- paste0(outdir, ".part")
if (dir.exists(outdir) || dir.exists(part) || file.exists(outdir) || file.exists(part)) stop("Refusing existing output")
dir.create(part, recursive=TRUE)

z <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm=TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm=TRUE)) / s
}

model_row <- function(cohort, id, tier, endpoint, stratum, feature, data, type,
                      beta, se, p, ph, adjustment, expected="positive_harmful") {
  data.table(cohort=cohort, model_id=id, analysis_tier=tier, endpoint=endpoint,
    stratum=stratum, feature=feature, n=nrow(data), n_event=sum(data$event),
    effect_type=type, effect=exp(beta), ci_low=exp(beta-1.96*se),
    ci_high=exp(beta+1.96*se), p_value=p, ph_p_value=ph,
    expected_direction=expected, observed_direction=ifelse(beta>0,"positive","negative"),
    adjustment=adjustment)
}

cox_term <- function(fit, term) {
  tab <- summary(fit)$coef
  if (!term %in% rownames(tab)) stop("Missing model term: ", term)
  ph <- try(cox.zph(fit), silent=TRUE)
  php <- if (inherits(ph,"try-error") || !term %in% rownames(ph$table)) NA_real_ else ph$table[term,"p"]
  list(beta=tab[term,"coef"], se=tab[term,"se(coef)"], p=tab[term,"Pr(>|z|)"], ph=php)
}

fit_cox <- function(frame, cohort, feature="spp1_macrophage_fraction",
                    biological=FALSE, endpoint="PFS", composition=FALSE, tier="PRIMARY") {
  d <- copy(frame)
  if (biological) {
    d[, spp1_z := z(spp1_macrophage_fraction)]
    d[, bio := tymp_z * spp1_z]
    term <- "bio"; feature_name <- "tymp_x_spp1_state"
  } else {
    d[, score := z(get(feature))]
    term <- "score"; feature_name <- feature
  }
  if (cohort == "Braun2020_CheckMate") {
    d <- d[Arm=="NIVOLUMAB"]; d[, time := as.numeric(get(endpoint))]
    d[, event := as.numeric(get(paste0(endpoint,"_CNSR")))]
    keep <- is.finite(d$time) & d$time>0 & d$event %in% c(0,1) & nzchar(d$Cohort) & nzchar(d$MSKCC)
    keep <- keep & if (biological) is.finite(d$tymp_z)&is.finite(d$spp1_z)&is.finite(d$bio) else is.finite(d$score)
    if (composition) keep <- keep & is.finite(d$official_macro_z)
    d <- d[keep]
    rhs <- if (biological) "tymp_z+spp1_z+bio" else "score"
    rhs <- paste(rhs,"factor(Cohort)+factor(MSKCC)",sep="+")
    if (composition) rhs <- paste(rhs,"official_macro_z",sep="+")
    fit <- coxph(as.formula(paste("Surv(time,event)~",rhs)), data=d, ties="efron")
    v <- cox_term(fit,term)
    return(model_row(cohort,paste0("cm_nivo_",tolower(endpoint),"_",feature_name),tier,
      endpoint,"nivolumab",feature_name,d,ifelse(biological,"interaction_hazard_ratio","hazard_ratio_per_sd"),
      v$beta,v$se,v$p,v$ph,ifelse(composition,"cohort+MSKCC+official_macrophage","cohort+MSKCC")))
  }
  d <- d[TRT01P=="Avelumab+Axitinib"]; d[, time:=as.numeric(PFS_P)]
  d[, event:=1-as.numeric(PFS_P_CNSR)]
  keep <- is.finite(d$time)&d$time>0&d$event %in% c(0,1)&is.finite(d$AGE)&nzchar(d$SEX)
  keep <- keep & if (biological) is.finite(d$tymp_z)&is.finite(d$spp1_z)&is.finite(d$bio) else is.finite(d$score)
  if (composition) keep <- keep & is.finite(d$official_macro_z)
  d <- d[keep]
  rhs <- if (biological) "tymp_z+spp1_z+bio" else "score"
  rhs <- paste(rhs,"AGE+factor(SEX)",sep="+")
  if (composition) rhs <- paste(rhs,"official_macro_z",sep="+")
  fit <- coxph(as.formula(paste("Surv(time,event)~",rhs)),data=d,ties="efron")
  v <- cox_term(fit,term)
  model_row(cohort,paste0("javelin_combo_pfs_",feature_name),tier,"PFS","avelumab_axitinib",
    feature_name,d,ifelse(biological,"interaction_hazard_ratio","hazard_ratio_per_sd"),
    v$beta,v$se,v$p,v$ph,ifelse(composition,"age+sex+official_macrophage","age+sex"))
}

fit_treatment <- function(frame, cohort) {
  d <- copy(frame); d[, score:=z(spp1_macrophage_fraction)]
  if (cohort=="Braun2020_CheckMate") {
    d <- d[Cohort=="CM-025"]; d[, time:=as.numeric(PFS)]; d[, event:=as.numeric(PFS_CNSR)]
    d[, treatment:=as.integer(Arm=="NIVOLUMAB")]
    d <- d[is.finite(score)&is.finite(time)&time>0&event %in% c(0,1)&nzchar(MSKCC)]
    fit <- coxph(Surv(time,event)~score*treatment+factor(MSKCC),d,ties="efron")
    adj <- "MSKCC"; stratum <- "CM-025_randomized"
  } else {
    d[, time:=as.numeric(PFS_P)]; d[, event:=1-as.numeric(PFS_P_CNSR)]
    d[, treatment:=as.integer(TRT01P=="Avelumab+Axitinib")]
    d <- d[is.finite(score)&is.finite(time)&time>0&event %in% c(0,1)&is.finite(AGE)&nzchar(SEX)]
    fit <- coxph(Surv(time,event)~score*treatment+AGE+factor(SEX),d,ties="efron")
    adj <- "age+sex"; stratum <- "randomized_all"
  }
  v <- cox_term(fit,"score:treatment")
  model_row(cohort,paste0(cohort,"_spp1_by_treatment"),"SECONDARY","PFS",stratum,
    "spp1_fraction_by_treatment",d,"interaction_hazard_ratio",v$beta,v$se,v$p,v$ph,
    adj,"negative_predictive_benefit")
}

fit_orr <- function(frame, biological=FALSE) {
  d <- copy(frame); d[, spp1_z:=z(spp1_macrophage_fraction)]; d[, bio:=tymp_z*spp1_z]
  d <- d[Arm=="NIVOLUMAB"&ORR %in% c("CR","PR","SD","PD")]
  d[, response:=as.integer(ORR %in% c("CR","PR"))]; d[, event:=response]
  if (biological) {
    fit <- glm(response~tymp_z+spp1_z+bio+factor(Cohort)+factor(MSKCC),d,family=binomial())
    term <- "bio"; feature <- "tymp_x_spp1_state"
  } else {
    fit <- glm(response~spp1_z+factor(Cohort)+factor(MSKCC),d,family=binomial())
    term <- "spp1_z"; feature <- "spp1_macrophage_fraction"
  }
  tab <- summary(fit)$coef
  model_row("Braun2020_CheckMate",paste0("cm_nivo_orr_",feature),"SECONDARY","ORR","nivolumab",
    feature,d,"odds_ratio",tab[term,"Estimate"],tab[term,"Std. Error"],tab[term,"Pr(>|z|)"],
    NA_real_,"cohort+MSKCC","negative_harmful")
}

weights <- fread(weights_path,check.names=FALSE); diag <- fread(diagnostics_path)
if (nrow(weights)!=1037L*13L || uniqueN(weights$state)!=13L) stop("Projection dimensions failed")
wide <- dcast(weights,dataset+sample_id~state,value.var="weight")
needed <- c("SPP1-state Macrophage","C1QC-state Macrophage","Other Macrophage")
if (!all(needed %in% names(wide))) stop("Macrophage states missing")
wide[, spp1_state_weight:=get("SPP1-state Macrophage")]
wide[, c1qc_state_weight:=get("C1QC-state Macrophage")]
wide[, other_macrophage_weight:=get("Other Macrophage")]
wide[, total_macrophage_weight:=spp1_state_weight+c1qc_state_weight+other_macrophage_weight]
wide[, spp1_macrophage_fraction:=fifelse(total_macrophage_weight>0,spp1_state_weight/total_macrophage_weight,NA_real_)]
wide[, c1qc_macrophage_fraction:=fifelse(total_macrophage_weight>0,c1qc_state_weight/total_macrophage_weight,NA_real_)]
wide[, spp1_c1qc_log_ratio:=log((spp1_state_weight+1e-6)/(c1qc_state_weight+1e-6))]

s81 <- fread(scores81_path)[,.(sample_id,tymp_z)]
if (anyDuplicated(s81$sample_id)) stop("Duplicated Phase81 sample IDs")
cm <- merge(wide[dataset=="CheckMate"],s81,by="sample_id",all.x=TRUE,sort=FALSE)
jav <- merge(wide[dataset=="JAVELIN"],s81,by="sample_id",all.x=TRUE,sort=FALSE)

cib <- fread(cm_ciber_path,check.names=FALSE)
cibm <- as.matrix(cib[,-1L,with=FALSE]); storage.mode(cibm)<-"double"
rownames(cibm) <- tolower(cib[[1L]])
cib_targets <- c("macrophages_m0","macrophages_m1","macrophages_m2")
if (!all(cib_targets %in% rownames(cibm)) || !all(cm$sample_id %in% colnames(cibm))) stop("CIBERSORTx merge failed")
cm_macro <- data.table(sample_id=cm$sample_id,official_macro=colSums(cibm[cib_targets,cm$sample_id,drop=FALSE]))
dec <- fread(jav_deconv_path,check.names=FALSE)
dec_targets <- c("Macrophages_M0","Macrophages_M1","Macrophages_M2")
if (!all(dec_targets %in% names(dec))) stop("ImmuneNet macrophage columns missing")
jav_macro <- dec[,.(sample_id=as.character(row_id),official_macro=rowSums(as.matrix(.SD))),.SDcols=dec_targets]
cm <- merge(cm,cm_macro,by="sample_id",all.x=TRUE,sort=FALSE)
jav <- merge(jav,jav_macro,by="sample_id",all.x=TRUE,sort=FALSE)
cm[,official_macro_z:=z(official_macro)]; jav[,official_macro_z:=z(official_macro)]
cm <- merge(cm,fread(cm_clin_path,check.names=FALSE),by.x="sample_id",by.y="RNA_ID",all.x=TRUE,sort=FALSE)
jav <- merge(jav,fread(jav_clin_path,check.names=FALSE),by.x="sample_id",by.y="ID",all.x=TRUE,sort=FALSE)
if (nrow(cm)!=311L || anyNA(cm$Arm) || anyNA(cm$tymp_z)) stop("CheckMate clinical merge failed")
if (nrow(jav)!=726L || anyNA(jav$TRT01P) || anyNA(jav$tymp_z)) stop("JAVELIN clinical merge failed")

cm_rho <- cor(cm$total_macrophage_weight,cm$official_macro,method="spearman",use="complete.obs")
jav_rho <- cor(jav$total_macrophage_weight,jav$official_macro,method="spearman",use="complete.obs")
checks <- rbindlist(list(
  data.table(check_id="projection_rows",observed=nrow(weights),threshold=1037L*13L,pass=nrow(weights)==1037L*13L),
  data.table(check_id="all_weights_finite",observed=sum(is.finite(weights$weight)),threshold=nrow(weights),pass=all(is.finite(weights$weight))),
  data.table(check_id="all_weights_nonnegative",observed=min(weights$weight),threshold=0,pass=min(weights$weight)>=-1e-12),
  data.table(check_id="weight_sum_error",observed=max(abs(diag$weight_sum-1)),threshold=1e-8,pass=max(abs(diag$weight_sum-1))<1e-8),
  data.table(check_id="cm_macro_nonzero_fraction",observed=mean(cm$total_macrophage_weight>0),threshold=.95,pass=mean(cm$total_macrophage_weight>0)>=.95),
  data.table(check_id="jav_macro_nonzero_fraction",observed=mean(jav$total_macrophage_weight>0),threshold=.95,pass=mean(jav$total_macrophage_weight>0)>=.95),
  data.table(check_id="cm_macro_external_spearman",observed=cm_rho,threshold=0,pass=is.finite(cm_rho)&&cm_rho>0),
  data.table(check_id="jav_macro_external_spearman",observed=jav_rho,threshold=0,pass=is.finite(jav_rho)&&jav_rho>0),
  data.table(check_id="cm_spp1_fraction_unique",observed=uniqueN(cm$spp1_macrophage_fraction),threshold=20,
    pass=uniqueN(cm$spp1_macrophage_fraction)>=20&&sd(cm$spp1_macrophage_fraction,na.rm=TRUE)>0),
  data.table(check_id="jav_spp1_fraction_unique",observed=uniqueN(jav$spp1_macrophage_fraction),threshold=20,
    pass=uniqueN(jav$spp1_macrophage_fraction)>=20&&sd(jav$spp1_macrophage_fraction,na.rm=TRUE)>0),
  data.table(check_id="cm_clinical_merge",observed=nrow(cm),threshold=311,pass=nrow(cm)==311L&&!anyNA(cm$Arm)),
  data.table(check_id="jav_clinical_merge",observed=nrow(jav),threshold=726,pass=nrow(jav)==726L&&!anyNA(jav$TRT01P))
))

primary <- rbindlist(list(
  fit_cox(cm,"Braun2020_CheckMate"),
  fit_cox(jav,"JAVELIN_Renal101"),
  fit_cox(cm,"Braun2020_CheckMate",biological=TRUE),
  fit_cox(jav,"JAVELIN_Renal101",biological=TRUE)
))
primary[,q_value:=p.adjust(p_value,"BH")]
primary[,significance_class:=fifelse(q_value<.05&effect>1,"COHORT_SUPPORTED_EXPECTED",
  fifelse(q_value<.05&effect<1,"COHORT_SUPPORTED_DIRECTION_CONFLICT",
    fifelse(p_value<.05&effect>1,"NOMINAL_EXPECTED",
      fifelse(p_value<.05&effect<1,"NOMINAL_DIRECTION_CONFLICT","NOT_SIGNIFICANT"))))]

secondary <- rbindlist(list(
  fit_cox(cm,"Braun2020_CheckMate","spp1_c1qc_log_ratio",tier="SECONDARY"),
  fit_cox(jav,"JAVELIN_Renal101","spp1_c1qc_log_ratio",tier="SECONDARY"),
  fit_cox(cm,"Braun2020_CheckMate","c1qc_macrophage_fraction",tier="SECONDARY"),
  fit_cox(jav,"JAVELIN_Renal101","c1qc_macrophage_fraction",tier="SECONDARY"),
  fit_cox(cm,"Braun2020_CheckMate","total_macrophage_weight",tier="SECONDARY"),
  fit_cox(jav,"JAVELIN_Renal101","total_macrophage_weight",tier="SECONDARY"),
  fit_cox(cm,"Braun2020_CheckMate",endpoint="OS",tier="SECONDARY"),
  fit_cox(cm,"Braun2020_CheckMate",biological=TRUE,endpoint="OS",tier="SECONDARY"),
  fit_orr(cm,FALSE),fit_orr(cm,TRUE),fit_treatment(cm,"Braun2020_CheckMate"),
  fit_treatment(jav,"JAVELIN_Renal101")
))
secondary[,q_value:=p.adjust(p_value,"BH")]
secondary[,significance_class:=fifelse(q_value<.05,"FDR_LT_0.05",fifelse(p_value<.05,"NOMINAL_P_LT_0.05","NOT_SIGNIFICANT"))]

sensitivity <- rbindlist(list(
  fit_cox(cm,"Braun2020_CheckMate",composition=TRUE,tier="COMPOSITION_SENSITIVITY"),
  fit_cox(jav,"JAVELIN_Renal101",composition=TRUE,tier="COMPOSITION_SENSITIVITY"),
  fit_cox(cm,"Braun2020_CheckMate",biological=TRUE,composition=TRUE,tier="COMPOSITION_SENSITIVITY"),
  fit_cox(jav,"JAVELIN_Renal101",biological=TRUE,composition=TRUE,tier="COMPOSITION_SENSITIVITY")
))
sensitivity[,q_value:=p.adjust(p_value,"BH")]

decisions <- rbindlist(lapply(c("spp1_macrophage_fraction","tymp_x_spp1_state"),function(f) {
  x <- primary[feature==f]; a <- x[cohort=="Braun2020_CheckMate"]; b <- x[cohort=="JAVELIN_Renal101"]
  dual <- a$effect>1&&a$q_value<.05&&b$effect>1&&b$q_value<.05
  one <- (a$effect>1&&a$q_value<.05&&b$effect>1&&b$p_value<.05) ||
    (b$effect>1&&b$q_value<.05&&a$effect>1&&a$p_value<.05)
  same <- a$effect>1&&b$effect>1
  data.table(feature=f,checkmate_effect=a$effect,checkmate_p=a$p_value,checkmate_q=a$q_value,
    javelin_effect=b$effect,javelin_p=b$p_value,javelin_q=b$q_value,
    dual_cohort_supported=dual,one_cohort_plus_external_nominal=one,
    expected_direction_both=same,decision=fifelse(dual,"DUAL_COHORT_SUPPORTED",
      fifelse(one,"ONE_COHORT_PLUS_EXTERNAL_NOMINAL",
        fifelse(same,"DIRECTIONALLY_CONSISTENT_NOT_SUPPORTED","DIRECTION_CONFLICT"))))
}))

checks <- rbind(checks,data.table(check_id="primary_model_count",observed=nrow(primary),
  threshold=4,pass=nrow(primary)==4L&&all(is.finite(primary$p_value))))
valid_ids <- c("all_weights_finite","all_weights_nonnegative","weight_sum_error",
  "cm_macro_nonzero_fraction","jav_macro_nonzero_fraction",
  "cm_macro_external_spearman","jav_macro_external_spearman",
  "cm_spp1_fraction_unique","jav_spp1_fraction_unique")
validity_pass <- all(checks[check_id %in% valid_ids,pass]); overall_pass <- all(checks$pass)

score_cols <- c("sample_id","dataset","tymp_z","spp1_state_weight","c1qc_state_weight",
  "other_macrophage_weight","total_macrophage_weight","spp1_macrophage_fraction",
  "c1qc_macrophage_fraction","spp1_c1qc_log_ratio","official_macro")
sample_scores <- rbindlist(list(cm[,..score_cols],jav[,..score_cols]))
fwrite(sample_scores,file.path(part,"phase83_sample_scores_v1.tsv"),sep="\t")
fwrite(diag,file.path(part,"phase83_projection_diagnostics_v1.tsv"),sep="\t")
fwrite(primary,file.path(part,"phase83_primary_models_v1.tsv"),sep="\t")
fwrite(secondary,file.path(part,"phase83_secondary_models_v1.tsv"),sep="\t")
fwrite(sensitivity,file.path(part,"phase83_composition_sensitivity_v1.tsv"),sep="\t")
fwrite(decisions,file.path(part,"phase83_cross_cohort_decisions_v1.tsv"),sep="\t")
fwrite(checks,file.path(part,"phase83_validation_checks_v1.tsv"),sep="\t")

summary_lines <- c("# Phase 83 bulk state projection summary v1","",
  paste0("- Projection validity gates: ",ifelse(validity_pass,"PASS","FAIL")),
  paste0("- Overall workflow checks: ",ifelse(overall_pass,"PASS","FAIL")),
  paste0("- CheckMate projected/official macrophage Spearman rho: ",signif(cm_rho,5)),
  paste0("- JAVELIN projected/official macrophage Spearman rho: ",signif(jav_rho,5)),
  paste0("- Median raw-NNLS relative reconstruction error: ",signif(median(diag$relative_reconstruction_error),5)),
  "","## Primary models","",
  paste(capture.output(print(primary[,.(cohort,feature,n,n_event,effect,ci_low,ci_high,p_value,q_value,ph_p_value,significance_class)])),collapse="\n"),
  "","## Cross-cohort decisions","",paste(capture.output(print(decisions)),collapse="\n"))
writeLines(summary_lines,file.path(part,"PHASE83_SUMMARY_v1.md"))
writeLines(ifelse(overall_pass,"PASS","FAIL"),file.path(part,
  ifelse(overall_pass,"PHASE83_PASS_v1","PHASE83_FAIL_v1")))
if (!file.rename(part,outdir)) stop("Failed to finalize output")
cat("Phase83 completed; validity=",validity_pass,"; overall=",overall_pass,"\n",sep="")
