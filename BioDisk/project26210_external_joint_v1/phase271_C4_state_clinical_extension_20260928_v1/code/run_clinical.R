library(data.table);library(survival)
h<-"/mnt/BioDisk/wangmingyuan";r<-file.path(h,"project26210_external_joint_v1")
p<-file.path(r,"phase271_C4_state_clinical_extension_20260928_v1");old<-file.path(r,"phase270_SPP1_literature_clinical_revalidation_20260927_v1")
scores<-fread(file.path(p,"tables/C4_CLINICAL_SCORES.tsv"));macro<-fread(file.path(h,"project26210_phase83_bulk_state_projection_v1/results_v1/clinical_v1/phase83_sample_scores_v1.tsv"))
out<-list();members<-list();coeff<-list();ii<-0
write_tab<-function(x,n)fwrite(x,file.path(p,"tables",n),sep="\t",na="NA")
get_data<-function(ds){
 d<-fread(file.path(old,"tables",paste0(toupper(ds),"_ANALYSIS_INPUT.tsv")),na.strings=c("","NA"))
 d<-merge(d,scores[,.(sample_id,C4_state)],by="sample_id",all.x=TRUE)
 d<-merge(d,macro[,.(sample_id,official_macro)],by="sample_id",all.x=TRUE)
 stopifnot(!anyDuplicated(d$sample_id),all(is.finite(d$C4_state)))
 d[,sex:=factor(fifelse(toupper(as.character(sex)) %in% c("M","MALE"),"Male",fifelse(toupper(as.character(sex)) %in% c("F","FEMALE"),"Female",NA_character_)))];if(ds=="CheckMate"){d[,trial:=factor(trial)];d[,risk:=factor(risk)]}
 d
}
for(ds in c("CheckMate","JAVELIN")){
 all<-get_data(ds);write_tab(all,paste0(toupper(ds),"_ANALYSIS_INPUT.tsv"))
 for(ep in if(ds=="CheckMate")c("OS","PFS") else "PFS")for(a in unique(all$arm)){
  group<-copy(all[arm==a])
  group[,':='(x=as.numeric(scale(C4_state)),spp=as.numeric(scale(SPP1)),lit=as.numeric(scale(Jiang5)),mac=as.numeric(scale(official_macro)),time=get(paste0(ep,"_time")),event=get(paste0(ep,"_event")))]
  for(model in c("clinical","plus_SPP1","plus_Jiang5","plus_Jiang5_macro")){
   cv<-if(ds=="CheckMate")c("age","sex","trial","risk") else c("age","sex")
   ext<-switch(model,clinical=character(),plus_SPP1="spp",plus_Jiang5="lit",plus_Jiang5_macro=c("lit","mac"))
   vars<-c("x","time","event",cv,ext);d<-group[complete.cases(group[,..vars])&time>=0]
   cv<-cv[vapply(cv,function(v)uniqueN(d[[v]])>1,logical(1))]
   baseterms<-c(cv,ext);f<-as.formula(paste("Surv(time,event)~",paste(c("x",baseterms),collapse="+")))
   f0<-as.formula(paste("Surv(time,event)~",paste(baseterms,collapse="+")))
   fit<-coxph(f,data=d,x=TRUE,y=TRUE);fit0<-coxph(f0,data=d,x=TRUE,y=TRUE)
   sm<-summary(fit);ci<-sm$conf.int["x",];co<-sm$coefficients["x",]
   lr<-2*(as.numeric(logLik(fit))-as.numeric(logLik(fit0)));lp<-pchisq(lr,1,lower.tail=FALSE)
   ph<-cox.zph(fit);vif<-1/(1-summary(lm(as.formula(paste("x~",paste(baseterms,collapse="+"))),data=d))$r.squared)
   ii<-ii+1;out[[ii]]<-data.table(model_id=ii,dataset=ds,endpoint=ep,arm=a,model=model,n=nrow(d),events=sum(d$event),HR=ci["exp(coef)"],CI_low=ci["lower .95"],CI_high=ci["upper .95"],P=co["Pr(>|z|)"],increment_LR=lr,increment_P=lp,PH_P=ph$table["x","p"],VIF=vif,formula=deparse(f))
   members[[ii]]<-d[,.(model_id=ii,sample_id,x,spp,lit,mac,age,sex,time,event,trial=if(ds=="CheckMate")as.character(trial) else NA_character_,risk=if(ds=="CheckMate")as.character(risk) else NA_character_)]
   cc<-as.data.table(sm$coefficients,keep.rownames="term");cc[,model_id:=ii];coeff[[ii]]<-cc
  }
 }
}
res<-rbindlist(out);res[,family:=paste(dataset,endpoint,model,sep="|")];res[,FDR:=p.adjust(P,"BH"),by=family];res[,increment_FDR:=p.adjust(increment_P,"BH"),by=family]
write_tab(res,"C4_SURVIVAL_MODELS.tsv");write_tab(rbindlist(members),"C4_MODEL_MEMBERSHIP.tsv");write_tab(rbindlist(coeff),"C4_ALL_MODEL_COEFFICIENTS.tsv")
io<-list();im<-list();ip<-list();k<-0
for(ds in c("CheckMate","JAVELIN")){
 all<-get_data(ds);if(ds=="CheckMate")all<-all[Cohort=="CM-025"]
 all[,':='(x=as.numeric(scale(C4_state)),lit=as.numeric(scale(Jiang5)),mac=as.numeric(scale(official_macro)),tx=as.integer(if(ds=="CheckMate")arm=="NIVOLUMAB" else arm=="Avelumab+Axitinib"))]
 for(ep in if(ds=="CheckMate")c("OS","PFS") else "PFS")for(model in c("clinical","plus_Jiang5_interaction_macro")){
  d<-copy(all);d[,':='(time=get(paste0(ep,"_time")),event=get(paste0(ep,"_event")))]
  cv<-if(ds=="CheckMate")c("age","sex","risk") else c("age","sex")
  ext<-if(model=="clinical")character() else c("lit","mac")
  vars<-c("x","tx","time","event",cv,ext);d<-d[complete.cases(d[,..vars])&time>=0]
  terms<-c("x*tx",cv,if(model!="clinical")c("lit*tx","mac"))
  fit<-coxph(as.formula(paste("Surv(time,event)~",paste(terms,collapse="+"))),data=d,x=TRUE,y=TRUE)
  terms0<-c("x","tx",cv,if(model!="clinical")c("lit*tx","mac"))
  fit0<-coxph(as.formula(paste("Surv(time,event)~",paste(terms0,collapse="+"))),data=d)
  sm<-summary(fit);ci<-sm$conf.int["x:tx",];co<-sm$coefficients["x:tx",];ph<-cox.zph(fit)
  lr<-2*(as.numeric(logLik(fit))-as.numeric(logLik(fit0)))
  k<-k+1;io[[k]]<-data.table(interaction_id=k,dataset=ds,endpoint=ep,model=model,n=nrow(d),events=sum(d$event),HR_ratio=ci["exp(coef)"],CI_low=ci["lower .95"],CI_high=ci["upper .95"],P=co["Pr(>|z|)"],increment_P=pchisq(lr,1,lower.tail=FALSE),PH_P=ph$table["x:tx","p"],formula=deparse(formula(fit)))
  im[[k]]<-d[,.(interaction_id=k,sample_id,x,lit,mac,tx,age,sex,time,event,risk=if(ds=="CheckMate")as.character(risk) else NA_character_)]
  tt<-as.data.table(ph$table,keep.rownames="term");tt[,interaction_id:=k];ip[[k]]<-tt
 }
}
inter<-rbindlist(io);inter[,family:=paste(dataset,endpoint,sep="|")];inter[,FDR:=p.adjust(P,"BH"),by=family];inter[,increment_FDR:=p.adjust(increment_P,"BH"),by=family]
write_tab(inter,"C4_TREATMENT_INTERACTIONS.tsv");write_tab(rbindlist(im),"C4_INTERACTION_MEMBERSHIP.tsv");write_tab(rbindlist(ip),"C4_INTERACTION_PH.tsv")
capture.output(sessionInfo(),file=file.path(p,"qc/R_SESSION_INFO.txt"))
print(as.data.frame(res[,.(dataset,endpoint,arm,model,n,HR,P,FDR,increment_P,PH_P,VIF)]));print(as.data.frame(inter))
writeLines("SUCCESS",file.path(p,"qc/CLINICAL_COMPLETE.txt"))
