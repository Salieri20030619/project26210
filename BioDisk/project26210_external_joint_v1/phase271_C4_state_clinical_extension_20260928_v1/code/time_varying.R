library(data.table);library(survival)
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase271_C4_state_clinical_extension_20260928_v1"
res<-fread(file.path(p,"tables/C4_SURVIVAL_MODELS.tsv"));mem<-fread(file.path(p,"tables/C4_MODEL_MEMBERSHIP.tsv"),na.strings=c("NA",""))
out<-list();hr<-list();cc<-list();k<-0
for(id in res[dataset=="CheckMate"&endpoint=="PFS"&model %in% c("clinical","plus_Jiang5_macro"),model_id]){
 d<-mem[model_id==id];meta<-res[model_id==id];d[,':='(sex=factor(sex),trial=factor(trial),risk=factor(risk))]
 cv<-c("age","sex","trial","risk",if(meta$model=="plus_Jiang5_macro")c("lit","mac"));cv<-cv[vapply(cv,function(v)uniqueN(d[[v]])>1,logical(1))]
 f<-as.formula(paste("Surv(time,event)~x+tt(x)+",paste(cv,collapse="+")))
 f0<-as.formula(paste("Surv(time,event)~",paste(cv,collapse="+")))
 fc<-as.formula(paste("Surv(time,event)~x+",paste(cv,collapse="+")))
 fit<-coxph(f,data=d,tt=function(x,t,...)x*log1p(t));fbase<-coxph(f0,data=d);fconst<-coxph(fc,data=d)
 b<-coef(fit);v<-vcov(fit);sm<-summary(fit);lr<-2*(as.numeric(logLik(fit))-as.numeric(logLik(fbase)));td<-2*(as.numeric(logLik(fit))-as.numeric(logLik(fconst)))
 k<-k+1;out[[k]]<-data.table(model_id=id,arm=meta$arm,model=meta$model,n=nrow(d),events=sum(d$event),joint_LR=lr,joint_P=pchisq(lr,2,lower.tail=FALSE),time_dependence_LR=td,time_dependence_P=pchisq(td,1,lower.tail=FALSE),beta_0=b["x"],beta_log1p_month=b["tt(x)"])
 coeff<-as.data.table(sm$coefficients,keep.rownames="term");coeff[,model_id:=id];cc[[k]]<-coeff
 for(tm in c(0,3,6,12,24)){
  l<-log1p(tm);be<-b["x"]+l*b["tt(x)"];se<-sqrt(v["x","x"]+l^2*v["tt(x)","tt(x)"]+2*l*v["x","tt(x)"])
  hr[[length(hr)+1]]<-data.table(model_id=id,arm=meta$arm,model=meta$model,month=tm,HR=exp(be),CI_low=exp(be-1.96*se),CI_high=exp(be+1.96*se),at_risk=sum(d$time>=tm),pointwise_P=2*pnorm(-abs(be/se)))
 }
}
z<-rbindlist(out);z[,family:=paste("CheckMate_PFS",model,sep="|")];z[,joint_FDR:=p.adjust(joint_P,"BH"),by=family];z[,time_dependence_FDR:=p.adjust(time_dependence_P,"BH"),by=family]
fwrite(z,file.path(p,"tables/C4_PFS_TIME_VARYING_MODELS.tsv"),sep="\t")
fwrite(rbindlist(hr),file.path(p,"tables/C4_PFS_TIME_SPECIFIC_HR.tsv"),sep="\t")
fwrite(rbindlist(cc),file.path(p,"tables/C4_PFS_TIME_VARYING_COEFFICIENTS.tsv"),sep="\t")
print(as.data.frame(z));print(as.data.frame(rbindlist(hr)))
# Treatment comparison with time-dependent C4 effect, restricted randomized CM025.
im<-fread(file.path(p,"tables/C4_INTERACTION_MEMBERSHIP.tsv"),na.strings=c("NA",""))
irs<-fread(file.path(p,"tables/C4_TREATMENT_INTERACTIONS.tsv"))
io<-list();ih<-list()
for(id in irs[dataset=="CheckMate"&endpoint=="PFS",interaction_id]){
 d<-im[interaction_id==id];meta<-irs[interaction_id==id];d[,':='(sex=factor(sex),risk=factor(risk),xt=x*tx)]
 cv<-c("age","sex","risk",if(meta$model!="clinical")c("lit*tx","mac"))
 f<-as.formula(paste("Surv(time,event)~x+tx+xt+tt(x)+tt(xt)+",paste(cv,collapse="+")))
 f0<-as.formula(paste("Surv(time,event)~x+tx+tt(x)+",paste(cv,collapse="+")))
 fit<-coxph(f,data=d,tt=function(x,t,...)x*log1p(t));base<-coxph(f0,data=d,tt=function(x,t,...)x*log1p(t));b<-coef(fit);v<-vcov(fit)
 lr<-2*(as.numeric(logLik(fit))-as.numeric(logLik(base)))
 io[[length(io)+1]]<-data.table(interaction_id=id,model=meta$model,n=nrow(d),events=sum(d$event),joint_LR=lr,joint_interaction_P=pchisq(lr,2,lower.tail=FALSE),beta_interaction_0=b["xt"],beta_interaction_log1p_month=b["tt(xt)"])
 for(tm in c(0,3,6,12,24)){
  l<-log1p(tm);be<-b["xt"]+l*b["tt(xt)"];se<-sqrt(v["xt","xt"]+l^2*v["tt(xt)","tt(xt)"]+2*l*v["xt","tt(xt)"])
  ih[[length(ih)+1]]<-data.table(interaction_id=id,model=meta$model,month=tm,HR_ratio=exp(be),CI_low=exp(be-1.96*se),CI_high=exp(be+1.96*se),at_risk=sum(d$time>=tm))
 }
}
j<-rbindlist(io);j[,FDR:=p.adjust(joint_interaction_P,"BH")]
fwrite(j,file.path(p,"tables/C4_TIME_VARYING_TREATMENT_INTERACTIONS.tsv"),sep="\t")
fwrite(rbindlist(ih),file.path(p,"tables/C4_TIME_VARYING_INTERACTION_HR.tsv"),sep="\t")
print(as.data.frame(j));writeLines("SUCCESS",file.path(p,"qc/TIME_VARYING_COMPLETE.txt"))
