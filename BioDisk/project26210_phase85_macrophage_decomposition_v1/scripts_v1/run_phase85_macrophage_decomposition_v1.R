#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=5L)stop("Usage: scores cm_clin jav_clin prereg outdir")
scores_path<-a[[1]];cm_path<-a[[2]];jav_path<-a[[3]];prereg_path<-a[[4]];outdir<-a[[5]]
part<-paste0(outdir,".part")
if(dir.exists(outdir)||dir.exists(part)||file.exists(outdir)||file.exists(part))stop("Refusing existing output")
dir.create(part,recursive=TRUE)
z<-function(x){x<-as.numeric(x);s<-sd(x,na.rm=TRUE);if(!is.finite(s)||s<=0)return(rep(0,length(x)));(x-mean(x,na.rm=TRUE))/s}
s<-fread(scores_path)
cm<-merge(s[dataset=="CheckMate"],fread(cm_path,check.names=FALSE),by.x="sample_id",by.y="RNA_ID",all.x=TRUE)
jav<-merge(s[dataset=="JAVELIN"],fread(jav_path,check.names=FALSE),by.x="sample_id",by.y="ID",all.x=TRUE)
cm[,':='(study="CheckMate",trial=as.character(Cohort),arm=as.character(Arm),age=as.numeric(Age),sex=as.character(Sex))]
jav[,':='(study="JAVELIN",trial="JAVELIN",arm=as.character(TRT01P),age=as.numeric(AGE),sex=as.character(SEX))]
if(nrow(cm)!=311L||nrow(jav)!=726L)stop("Merge failed")

fit_one<-function(d,outcome,condition_total=TRUE){
  x<-copy(d);x[,y:=z(get(outcome))];x[,total_z:=z(total_macrophage_weight)]
  x[,official_z:=z(official_macro)]
  terms<-c("tymp_z","age","factor(sex)")
  if(uniqueN(x$trial)>1)terms<-c(terms,"factor(trial)")
  if(uniqueN(x$arm)>1)terms<-c(terms,"factor(arm)")
  if(condition_total)terms<-c(terms,"total_z")
  terms<-c(terms,"official_z")
  fit<-lm(as.formula(paste("y~",paste(terms,collapse="+"))),data=x)
  tab<-summary(fit)$coef
  data.table(study=unique(x$study),outcome=outcome,
    model=ifelse(condition_total,"conditional_on_total_macro","coexpansion"),
    n=nrow(model.frame(fit)),beta=tab["tymp_z","Estimate"],se=tab["tymp_z","Std. Error"],
    ci_low=tab["tymp_z","Estimate"]-1.96*tab["tymp_z","Std. Error"],
    ci_high=tab["tymp_z","Estimate"]+1.96*tab["tymp_z","Std. Error"],
    p_value=tab["tymp_z","Pr(>|t|)"],adjustment=paste(setdiff(terms,"tymp_z"),collapse="+"))
}
components<-c("spp1_state_weight","c1qc_state_weight","other_macrophage_weight")
primary<-rbindlist(c(lapply(components,function(o)fit_one(cm,o,TRUE)),
  lapply(components,function(o)fit_one(jav,o,TRUE))))
primary[,q_value:=p.adjust(p_value,"BH")]
primary[,class:=fifelse(q_value<.05&beta>0,"FDR_POSITIVE",
  fifelse(q_value<.05&beta<0,"FDR_NEGATIVE",fifelse(p_value<.05,"NOMINAL","NOT_SIGNIFICANT")))]
secondary<-rbindlist(c(lapply(components,function(o)fit_one(cm,o,FALSE)),
  lapply(components,function(o)fit_one(jav,o,FALSE))))
secondary[,q_value:=p.adjust(p_value,"BH")]
decision<-rbindlist(lapply(components,function(o){
  x<-primary[outcome==o];u<-x[study=="CheckMate"];v<-x[study=="JAVELIN"]
  data.table(outcome=o,checkmate_beta=u$beta,checkmate_p=u$p_value,checkmate_q=u$q_value,
    javelin_beta=v$beta,javelin_p=v$p_value,javelin_q=v$q_value,
    positive_both=u$beta>0&&v$beta>0,supported_both=u$beta>0&&u$q_value<.05&&v$beta>0&&v$q_value<.05)
}))
checks<-data.table(check_id=c("cm_rows","jav_rows","primary_models","finite_p"),
  observed=c(nrow(cm),nrow(jav),nrow(primary),sum(is.finite(primary$p_value))),
  threshold=c(311,726,6,6),pass=c(nrow(cm)==311L,nrow(jav)==726L,nrow(primary)==6L,all(is.finite(primary$p_value))))
fwrite(primary,file.path(part,"phase85_primary_conditional_models_v1.tsv"),sep="\t")
fwrite(secondary,file.path(part,"phase85_secondary_coexpansion_models_v1.tsv"),sep="\t")
fwrite(decision,file.path(part,"phase85_cross_cohort_decisions_v1.tsv"),sep="\t")
fwrite(checks,file.path(part,"phase85_validation_checks_v1.tsv"),sep="\t")
writeLines(c("# Phase 85 macrophage decomposition summary v1","",
  paste0("- Validation: ",ifelse(all(checks$pass),"PASS","FAIL")),"","## Conditional models","",
  paste(capture.output(print(primary)),collapse="\n"),"","## Unconditioned co-expansion","",
  paste(capture.output(print(secondary)),collapse="\n")),file.path(part,"PHASE85_SUMMARY_v1.md"))
writeLines(ifelse(all(checks$pass),"PASS","FAIL"),file.path(part,ifelse(all(checks$pass),"PHASE85_PASS_v1","PHASE85_FAIL_v1")))
if(!file.rename(part,outdir))stop("Finalize failed")
cat("Phase85 complete\n")
