#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=6L) stop("Usage: scores cm_clin jav_clin markers prereg outdir")
scores_path<-args[[1]]; cm_path<-args[[2]]; jav_path<-args[[3]]
markers_path<-args[[4]]; prereg_path<-args[[5]]; outdir<-args[[6]]
part<-paste0(outdir,".part")
if (dir.exists(outdir)||dir.exists(part)||file.exists(outdir)||file.exists(part)) stop("Refusing existing output")
dir.create(part,recursive=TRUE)

z <- function(x) {
  x<-as.numeric(x); s<-sd(x,na.rm=TRUE)
  if (!is.finite(s)||s<=0) return(rep(0,length(x)))
  (x-mean(x,na.rm=TRUE))/s
}

scores<-fread(scores_path)
markers<-fread(markers_path)
if (any(toupper(markers$gene) %in% c("TYMP","SPP1"))) stop("Target leakage in projection markers")
cm<-merge(scores[dataset=="CheckMate"],fread(cm_path,check.names=FALSE),
  by.x="sample_id",by.y="RNA_ID",all.x=TRUE,sort=FALSE)
jav<-merge(scores[dataset=="JAVELIN"],fread(jav_path,check.names=FALSE),
  by.x="sample_id",by.y="ID",all.x=TRUE,sort=FALSE)
if (nrow(cm)!=311L||nrow(jav)!=726L) stop("Clinical merge failed")
cm[,':='(study="CheckMate",trial=as.character(Cohort),arm=as.character(Arm),
  age=as.numeric(Age),sex=as.character(Sex))]
jav[,':='(study="JAVELIN",trial="JAVELIN",arm=as.character(TRT01P),
  age=as.numeric(AGE),sex=as.character(SEX))]

fit_assoc <- function(d,outcome,stratum="all",macro=TRUE,tier="PRIMARY") {
  x<-copy(d); x[,outcome_z:=z(get(outcome))]; x[,official_macro_z:=z(official_macro)]
  terms<-c("tymp_z","age")
  if (uniqueN(x$sex[!is.na(x$sex)&nzchar(x$sex)])>1) terms<-c(terms,"factor(sex)")
  if (uniqueN(x$trial[!is.na(x$trial)&nzchar(x$trial)])>1) terms<-c(terms,"factor(trial)")
  if (uniqueN(x$arm[!is.na(x$arm)&nzchar(x$arm)])>1) terms<-c(terms,"factor(arm)")
  if (macro) terms<-c(terms,"official_macro_z")
  fit<-lm(as.formula(paste("outcome_z~",paste(terms,collapse="+"))),data=x)
  tab<-summary(fit)$coef; used<-model.frame(fit)
  data.table(study=unique(x$study),model_id=paste(unique(x$study),outcome,stratum,
      ifelse(macro,"macro_adjusted","no_macro"),sep="_"),
    analysis_tier=tier,stratum=stratum,outcome=outcome,exposure="tymp_z",
    n=nrow(used),beta=tab["tymp_z","Estimate"],se=tab["tymp_z","Std. Error"],
    ci_low=tab["tymp_z","Estimate"]-1.96*tab["tymp_z","Std. Error"],
    ci_high=tab["tymp_z","Estimate"]+1.96*tab["tymp_z","Std. Error"],
    p_value=tab["tymp_z","Pr(>|t|)"],adjustment=paste(setdiff(terms,"tymp_z"),collapse="+"),
    expected_direction="positive")
}

primary_outcomes<-c("spp1_macrophage_fraction","spp1_c1qc_log_ratio")
primary<-rbindlist(c(
  lapply(primary_outcomes,function(o) fit_assoc(cm,o)),
  lapply(primary_outcomes,function(o) fit_assoc(jav,o))
))
primary[,q_value:=p.adjust(p_value,"BH")]
primary[,significance_class:=fifelse(q_value<.05&beta>0,"COHORT_SUPPORTED_EXPECTED",
  fifelse(q_value<.05&beta<0,"COHORT_SUPPORTED_DIRECTION_CONFLICT",
    fifelse(p_value<.05&beta>0,"NOMINAL_EXPECTED",
      fifelse(p_value<.05&beta<0,"NOMINAL_DIRECTION_CONFLICT","NOT_SIGNIFICANT"))))]

spearman<-rbindlist(c(lapply(list(cm,jav),function(d) rbindlist(lapply(primary_outcomes,function(o) {
  test<-cor.test(d$tymp_z,d[[o]],method="spearman",exact=FALSE)
  data.table(study=unique(d$study),outcome=o,n=sum(complete.cases(d$tymp_z,d[[o]])),
    rho=unname(test$estimate),p_value=test$p.value)
})))))
spearman[,q_value:=p.adjust(p_value,"BH")]

sensitivity<-rbindlist(list(
  fit_assoc(cm,"spp1_macrophage_fraction",macro=FALSE,tier="SENSITIVITY"),
  fit_assoc(cm,"spp1_c1qc_log_ratio",macro=FALSE,tier="SENSITIVITY"),
  fit_assoc(jav,"spp1_macrophage_fraction",macro=FALSE,tier="SENSITIVITY"),
  fit_assoc(jav,"spp1_c1qc_log_ratio",macro=FALSE,tier="SENSITIVITY"),
  fit_assoc(cm[arm=="NIVOLUMAB"],"spp1_macrophage_fraction","nivolumab",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(cm[arm=="NIVOLUMAB"],"spp1_c1qc_log_ratio","nivolumab",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(cm[arm=="EVEROLIMUS"],"spp1_macrophage_fraction","everolimus",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(cm[arm=="EVEROLIMUS"],"spp1_c1qc_log_ratio","everolimus",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(jav[arm=="Avelumab+Axitinib"],"spp1_macrophage_fraction","avelumab_axitinib",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(jav[arm=="Avelumab+Axitinib"],"spp1_c1qc_log_ratio","avelumab_axitinib",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(jav[arm=="Sunitinib"],"spp1_macrophage_fraction","sunitinib",TRUE,"ARM_SENSITIVITY"),
  fit_assoc(jav[arm=="Sunitinib"],"spp1_c1qc_log_ratio","sunitinib",TRUE,"ARM_SENSITIVITY")
))
sensitivity[,q_value:=p.adjust(p_value,"BH")]

secondary_outcomes<-c("spp1_state_weight","c1qc_macrophage_fraction","total_macrophage_weight")
secondary<-rbindlist(c(lapply(secondary_outcomes,function(o) fit_assoc(cm,o,tier="SECONDARY")),
  lapply(secondary_outcomes,function(o) fit_assoc(jav,o,tier="SECONDARY"))))
secondary[,q_value:=p.adjust(p_value,"BH")]

decisions<-rbindlist(lapply(primary_outcomes,function(o) {
  a<-primary[study=="CheckMate"&outcome==o]; b<-primary[study=="JAVELIN"&outcome==o]
  dual<-a$beta>0&&a$q_value<.05&&b$beta>0&&b$q_value<.05
  one<-(a$beta>0&&a$q_value<.05&&b$beta>0&&b$p_value<.05)||
    (b$beta>0&&b$q_value<.05&&a$beta>0&&a$p_value<.05)
  same<-a$beta>0&&b$beta>0
  data.table(outcome=o,checkmate_beta=a$beta,checkmate_p=a$p_value,checkmate_q=a$q_value,
    javelin_beta=b$beta,javelin_p=b$p_value,javelin_q=b$q_value,
    dual_cohort_supported=dual,one_cohort_plus_external_nominal=one,
    expected_direction_both=same,decision=fifelse(dual,"DUAL_COHORT_SUPPORTED",
      fifelse(one,"ONE_COHORT_PLUS_EXTERNAL_NOMINAL",
        fifelse(same,"DIRECTIONALLY_CONSISTENT_NOT_SUPPORTED","DIRECTION_CONFLICT"))))
}))

checks<-data.table(check_id=c("marker_target_exclusion","checkmate_rows","javelin_rows",
  "primary_model_count","primary_finite"),observed=c(1,nrow(cm),nrow(jav),nrow(primary),
  sum(is.finite(primary$p_value))),threshold=c(1,311,726,4,4),
  pass=c(!any(toupper(markers$gene)%in%c("TYMP","SPP1")),nrow(cm)==311L,nrow(jav)==726L,
    nrow(primary)==4L,all(is.finite(primary$p_value))))
fwrite(primary,file.path(part,"phase84_primary_models_v1.tsv"),sep="\t")
fwrite(spearman,file.path(part,"phase84_spearman_v1.tsv"),sep="\t")
fwrite(sensitivity,file.path(part,"phase84_sensitivity_v1.tsv"),sep="\t")
fwrite(secondary,file.path(part,"phase84_secondary_outcomes_v1.tsv"),sep="\t")
fwrite(decisions,file.path(part,"phase84_cross_cohort_decisions_v1.tsv"),sep="\t")
fwrite(checks,file.path(part,"phase84_validation_checks_v1.tsv"),sep="\t")
lines<-c("# Phase 84 TYMP-state coupling summary v1","",
  paste0("- Validation: ",ifelse(all(checks$pass),"PASS","FAIL")),"","## Primary models","",
  paste(capture.output(print(primary)),collapse="\n"),"","## Decisions","",
  paste(capture.output(print(decisions)),collapse="\n"))
writeLines(lines,file.path(part,"PHASE84_SUMMARY_v1.md"))
writeLines(ifelse(all(checks$pass),"PASS","FAIL"),file.path(part,
  ifelse(all(checks$pass),"PHASE84_PASS_v1","PHASE84_FAIL_v1")))
if (!file.rename(part,outdir)) stop("Failed to finalize output")
cat("Phase84 complete\n")
