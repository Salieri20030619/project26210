library(data.table);library(survival)
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase271_C4_state_clinical_extension_20260928_v1"
mem<-fread(file.path(p,"tables/C4_MODEL_MEMBERSHIP.tsv"),na.strings=c("NA",""))
out<-list();hrs<-list()
for(id in c(9,12))for(sens in c("stratify_trial","CM025_only")){
 d<-mem[model_id==id];if(sens=="CM025_only")d<-d[trial=="CM-025"]
 d[,':='(sex=factor(sex),risk=factor(risk),trial=factor(trial))]
 cv<-c("age","sex","risk",if(id==12)c("lit","mac"),if(sens=="stratify_trial")"strata(trial)")
 fit<-coxph(as.formula(paste("Surv(time,event)~x+tt(x)+",paste(cv,collapse="+"))),data=d,tt=function(x,t,...)x*log1p(t))
 base<-coxph(as.formula(paste("Surv(time,event)~",paste(cv,collapse="+"))),data=d)
 cons<-coxph(as.formula(paste("Surv(time,event)~x+",paste(cv,collapse="+"))),data=d)
 b<-coef(fit);v<-vcov(fit);lr<-2*(as.numeric(logLik(fit))-as.numeric(logLik(base)));ltd<-2*(as.numeric(logLik(fit))-as.numeric(logLik(cons)))
 out[[length(out)+1]]<-data.table(parent_model_id=id,sensitivity=sens,n=nrow(d),events=sum(d$event),beta_0=b["x"],beta_log1p_month=b["tt(x)"],joint_P=pchisq(lr,2,lower.tail=FALSE),time_dependence_P=pchisq(ltd,1,lower.tail=FALSE))
 for(tm in c(3,6,12)){
  ll<-log1p(tm);be<-b["x"]+ll*b["tt(x)"];se<-sqrt(v["x","x"]+ll^2*v["tt(x)","tt(x)"]+2*ll*v["x","tt(x)"])
  hrs[[length(hrs)+1]]<-data.table(parent_model_id=id,sensitivity=sens,month=tm,HR=exp(be),CI_low=exp(be-1.96*se),CI_high=exp(be+1.96*se),at_risk=sum(d$time>=tm))
 }
}
z<-rbindlist(out);z[,joint_FDR:=p.adjust(joint_P,"BH"),by=parent_model_id]
fwrite(z,file.path(p,"tables/C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv"),sep="\t")
fwrite(rbindlist(hrs),file.path(p,"tables/C4_TIME_VARIATION_TRIAL_HR.tsv"),sep="\t");print(as.data.frame(z))
# Influence check on the fully adjusted time-dependent nivo model; no sample removal from final estimate.
d<-mem[model_id==12];d[,':='(sex=factor(sex),risk=factor(risk),trial=factor(trial))]
loo<-list()
for(i in seq_len(nrow(d))){
 z<-d[-i];fit<-coxph(Surv(time,event)~x+tt(x)+age+sex+trial+risk+lit+mac,data=z,tt=function(x,t,...)x*log1p(t))
 b<-coef(fit);loo[[i]]<-data.table(omitted_patient=d$sample_id[i],beta_0=b["x"],beta_log1p_month=b["tt(x)"],HR_month6=exp(b["x"]+log1p(6)*b["tt(x)"]),HR_month12=exp(b["x"]+log1p(12)*b["tt(x)"]))
}
fwrite(rbindlist(loo),file.path(p,"tables/C4_TIME_VARIATION_LEAVE_ONE_PATIENT.tsv"),sep="\t");print(as.data.frame(rbindlist(loo)[,lapply(.SD,range),.SDcols=c("beta_0","beta_log1p_month","HR_month6","HR_month12")]))
writeLines("SUCCESS",file.path(p,"qc/TIME_VARIATION_SENSITIVITY_COMPLETE.txt"))
