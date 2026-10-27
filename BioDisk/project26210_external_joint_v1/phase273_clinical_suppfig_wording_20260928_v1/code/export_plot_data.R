library(data.table);library(survival)
r<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
p<-file.path(r,"phase273_clinical_suppfig_wording_20260928_v1")
a<-file.path(r,"phase270_SPP1_literature_clinical_revalidation_20260927_v1")
c<-file.path(r,"phase271_C4_state_clinical_extension_20260928_v1")
cm<-fread(file.path(a,"tables/CHECKMATE_ANALYSIS_INPUT.tsv"))
cutoff<-median(cm$Jiang5);d<-copy(cm[arm=="NIVOLUMAB"])
d[,group:=factor(ifelse(Jiang5>=cutoff,"High","Low"),levels=c("Low","High"))]
stopifnot(nrow(d)==181,sum(d$group=="High")==84,sum(d$group=="Low")==97,sum(d$OS_event)==123)
fit<-survfit(Surv(OS_time,OS_event)~group,data=d,conf.type="log")
sm<-summary(fit,censored=TRUE)
km<-data.table(time=sm$time,survival=sm$surv,lower=sm$lower,upper=sm$upper,n_risk=sm$n.risk,n_event=sm$n.event,n_censor=sm$n.censor,group=sub("group=","",as.character(sm$strata)))
start<-data.table(time=0,survival=1,lower=1,upper=1,n_risk=c(97,84),n_event=0,n_censor=0,group=c("Low","High"))
km<-rbind(start,km);setorder(km,group,time);fwrite(km,file.path(p,"data/KM_POINTS.tsv"),sep="\t")
end<-ceiling(max(d$OS_time)/12)*12;times<-seq(0,end,12)
risk<-rbindlist(lapply(c("Low","High"),function(g)data.table(group=g,time=times,n_risk=sapply(times,function(t)sum(d$group==g&d$OS_time>=t)))))
fwrite(risk,file.path(p,"data/KM_RISK_TABLE.tsv"),sep="\t")
lr<-survdiff(Surv(OS_time,OS_event)~group,data=d);lp<-pchisq(lr$chisq,1,lower.tail=FALSE)
orig<-fread(file.path(a,"tables/ALL_SURVIVAL_MODELS.tsv"))[dataset=="CheckMate"&endpoint=="OS"&arm=="NIVOLUMAB"&score=="Jiang5"&model=="median_cohort"]
stopifnot(abs(lp-orig$logrank_P)<1e-10)
fwrite(data.table(cutoff=cutoff,n=181,events=123,high=84,low=97,xmax=end,logrank_P=lp,logrank_FDR=orig$logrank_FDR,HR=orig$HR,CI_low=orig$CI_low,CI_high=orig$CI_high),file.path(p,"data/KM_STATISTICS.tsv"),sep="\t")
im<-fread(file.path(a,"tables/TREATMENT_INTERACTIONS.tsv"))[score=="Jiang5"]
fwrite(im,file.path(p,"data/LITERATURE_INTERACTIONS.tsv"),sep="\t")
dd<-fread(file.path(c,"tables/C4_MODEL_MEMBERSHIP.tsv"),na.strings=c("NA",""))[model_id==12]
dd[,':='(sex=factor(sex),trial=factor(trial),risk=factor(risk))]
f<-coxph(Surv(time,event)~x+tt(x)+age+sex+trial+risk+lit+mac,data=dd,tt=function(x,t,...)x*log1p(t))
be<-coef(f);vv<-vcov(f);orig<-fread(file.path(c,"tables/C4_PFS_TIME_VARYING_MODELS.tsv"))[model_id==12]
stopifnot(abs(be["x"]-orig$beta_0)<1e-10,abs(be["tt(x)"]-orig$beta_log1p_month)<1e-10)
ts<-seq(0,24,by=.1);lt<-log1p(ts);est<-be["x"]+lt*be["tt(x)"];se<-sqrt(vv["x","x"]+lt^2*vv["tt(x)","tt(x)"]+2*lt*vv["x","tt(x)"])
curve<-data.table(month=ts,HR=exp(est),CI_low=exp(est-1.96*se),CI_high=exp(est+1.96*se))
fwrite(curve,file.path(p,"data/C4_TIME_HR_CURVE.tsv"),sep="\t")
fwrite(data.table(month=seq(0,24,6),n_risk=sapply(seq(0,24,6),function(t)sum(dd$time>=t))),file.path(p,"data/C4_TIME_RISK_TABLE.tsv"),sep="\t")
orighr<-fread(file.path(c,"tables/C4_PFS_TIME_SPECIFIC_HR.tsv"))[model_id==12]
for(tm in c(0,3,6,12,24)){old<-orighr[month==tm];new<-curve[month==tm];stopifnot(max(abs(as.numeric(old[,.(HR,CI_low,CI_high)])-as.numeric(new[,.(HR,CI_low,CI_high)])))<1e-9)}
sen<-fread(file.path(c,"tables/C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv"))[parent_model_id==12]
senhr<-fread(file.path(c,"tables/C4_TIME_VARIATION_TRIAL_HR.tsv"))[parent_model_id==12&month==12]
sen<-merge(sen,senhr,by=c("parent_model_id","sensitivity"))
sens<-rbind(data.table(model="Full adjustment",n=orig$n,events=orig$events,HR=orighr[month==12,HR],CI_low=orighr[month==12,CI_low],CI_high=orighr[month==12,CI_high],P=orig$joint_P,FDR=orig$joint_FDR),
sen[,.(model=ifelse(sensitivity=="CM025_only","CM025 only","Trial stratified"),n,events,HR,CI_low,CI_high,P=joint_P,FDR=joint_FDR)])
fwrite(sens,file.path(p,"data/C4_SENSITIVITY.tsv"),sep="\t");fwrite(orig,file.path(p,"data/C4_TIME_MODEL.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(p,"qc/R_SESSION_INFO.txt"))
writeLines("PASS: frozen KM and time-model values reproduced; original source tables unchanged",file.path(p,"qc/PLOT_DATA_VALIDATION.txt"))
print(as.data.frame(sens));print(data.frame(KM_n=181,logrank_P=lp,xmax=end))
