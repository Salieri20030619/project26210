library(data.table)
library(survival)
h <- "/mnt/BioDisk/wangmingyuan"
p <- file.path(h,"project26210_external_joint_v1/phase270_SPP1_literature_clinical_revalidation_20260927_v1")
b <- file.path(h,"project26210_external_joint_v1_phase64_large_immune_validation_v1/inputs_v1")
genes <- c("SPP1","FN1","CLEC5A","OLR1","IL1RN")
getexpr <- function(file,log=FALSE) {
 x<-fread(file); names(x)[1]<-"gene";stopifnot(!anyDuplicated(x$gene), all(c(genes,"TYMP") %in% x$gene))
 m<-as.matrix(x[,-1]);rownames(m)<-x$gene;storage.mode(m)<-"double"
 stopifnot(all(is.finite(m)),!anyDuplicated(colnames(m)))
 if(log) {stopifnot(all(m>=0));m<-log2(m+1)}
 data.table(sample_id=colnames(m),SPP1=as.numeric(m["SPP1",]),Jiang5=colMeans(m[genes,,drop=FALSE]),TYMP=as.numeric(m["TYMP",]))
}
cm<-fread(file.path(b,"CheckMate_Braun2020/processed_v1/checkmate_clinical_rna_subset_v1.tsv"),na.strings=c("","NA"))
stopifnot(!anyDuplicated(cm$SUBJID),!anyDuplicated(cm$RNA_ID))
cm<-merge(cm,getexpr(file.path(p,"inputs/checkmate_genes.tsv")),by.x="RNA_ID",by.y="sample_id",all=FALSE)
cm[,':='(sample_id=RNA_ID,dataset="CheckMate",age=as.numeric(Age),sex=factor(Sex),arm=Arm,trial=factor(Cohort),risk=factor(MSKCC),PFS_time=as.numeric(PFS),PFS_event=as.integer(PFS_CNSR),OS_time=as.numeric(OS),OS_event=as.integer(OS_CNSR))]
ja<-fread(file.path(b,"JAVELIN_Renal101/processed_v1/javelin_clinical_v1.tsv"),na.strings=c("","NA"))
stopifnot(!anyDuplicated(ja$ID))
ja<-merge(ja,getexpr(file.path(p,"inputs/javelin_genes.tsv"),TRUE),by.x="ID",by.y="sample_id",all=FALSE)
ja[,':='(sample_id=ID,dataset="JAVELIN",age=as.numeric(AGE),sex=factor(SEX),arm=TRT01P,PFS_time=as.numeric(PFS_P),PFS_event=1L-as.integer(PFS_P_CNSR))]
stopifnot(nrow(cm)==311,nrow(ja)==726,all(cm$PFS_event %in% 0:1),all(cm$OS_event %in% 0:1),all(ja$PFS_event %in% 0:1))
write_table<-function(x,n) fwrite(x,file.path(p,"tables",n),sep="\t",na="NA")
write_table(cm,"CHECKMATE_ANALYSIS_INPUT.tsv");write_table(ja,"JAVELIN_ANALYSIS_INPUT.tsv")
out<-list(); members<-list(); ii<-0
for (ds in c("CheckMate","JAVELIN")) {
 all<-if(ds=="CheckMate") cm else ja
 for(endpoint in if(ds=="CheckMate") c("PFS","OS") else "PFS") {
  for (sc in c("SPP1","Jiang5")) for(a in unique(all$arm)) {
   base<-copy(all[arm==a])
   for(mod in c("median_arm","median_cohort","continuous","continuous_adjusted")) {
    d<-copy(base);cut<-NA_real_
    if(grepl("median",mod)){
     cut<-median(if(mod=="median_arm") base[[sc]] else all[[sc]])
     d[,x:=as.numeric(get(sc)>=cut)]
    } else d[,x:=as.numeric(scale(get(sc)))]
    d[,':='(time=get(paste0(endpoint,"_time")),event=get(paste0(endpoint,"_event")))]
    cv<-if(mod=="continuous_adjusted") if(ds=="CheckMate") c("age","sex","trial","risk") else c("age","sex") else character()
    vars<-c("time","event","x",cv);d<-d[complete.cases(d[,..vars]) & time>=0]
    cv<-cv[vapply(cv,function(z) uniqueN(d[[z]])>1,logical(1))]
    f<-as.formula(paste("Surv(time,event)~",paste(c("x",cv),collapse="+")))
    fit<-coxph(f,data=d,x=TRUE,y=TRUE); sm<-summary(fit);co<-sm$coefficients["x",];ci<-sm$conf.int["x",]
    ph<-tryCatch(cox.zph(fit)$table["x","p"],error=function(e) NA_real_)
    lp<-NA_real_;nh<-NA_integer_; nl<-NA_integer_
    if(grepl("median",mod)){sd<-survdiff(Surv(time,event)~x,data=d);lp<-pchisq(sd$chisq,1,lower.tail=FALSE);nh<-sum(d$x==1);nl<-sum(d$x==0)}
    ii<-ii+1
    out[[ii]]<-data.table(model_id=ii,dataset=ds,endpoint=endpoint,score=sc,arm=a,model=mod,n=nrow(d),events=sum(d$event),high_n=nh,low_n=nl,cutoff=cut,HR=ci["exp(coef)"],CI_low=ci["lower .95"],CI_high=ci["upper .95"],P=co["Pr(>|z|)"],logrank_P=lp,PH_P=ph,formula=deparse(f))
    members[[ii]]<-d[,.(model_id=ii,sample_id,x,time,event)]
   }
  }
 }
}
res<-rbindlist(out);res[,family:=paste(dataset,endpoint,model,sep="|")]
res[,FDR:=p.adjust(P,"BH"),by=family];res[,logrank_FDR:=p.adjust(logrank_P,"BH"),by=family]
write_table(res,"ALL_SURVIVAL_MODELS.tsv");write_table(rbindlist(members),"MODEL_MEMBERSHIP.tsv")
inter<-list();k<-0
for(ds in c("CheckMate","JAVELIN")){
 all<-if(ds=="CheckMate") copy(cm[Cohort=="CM-025"]) else copy(ja)
 all[,tx:=as.integer(if(ds=="CheckMate") arm=="NIVOLUMAB" else arm=="Avelumab+Axitinib")]
 for(ep in if(ds=="CheckMate") c("PFS","OS") else "PFS") for(sc in c("SPP1","Jiang5")){
  d<-copy(all);d[,':='(x=as.numeric(scale(get(sc))),time=get(paste0(ep,"_time")),event=get(paste0(ep,"_event")))]
  cv<-if(ds=="CheckMate") c("age","sex","risk") else c("age","sex")
  vars<-c("time","event","x","tx",cv);d<-d[complete.cases(d[,..vars]) & time>=0]
  f<-as.formula(paste("Surv(time,event)~x*tx+",paste(cv,collapse="+")))
  fit<-coxph(f,data=d,x=TRUE,y=TRUE);sm<-summary(fit);ci<-sm$conf.int["x:tx",];co<-sm$coefficients["x:tx",]
  k<-k+1;inter[[k]]<-data.table(dataset=ds,endpoint=ep,score=sc,n=nrow(d),events=sum(d$event),HR_ratio=ci["exp(coef)"],CI_low=ci["lower .95"],CI_high=ci["upper .95"],P=co["Pr(>|z|)"],formula=deparse(f))
 }
}
it<-rbindlist(inter);it[,family:=paste(dataset,endpoint,"interaction",sep="|")];it[,FDR:=p.adjust(P,"BH"),by=family];write_table(it,"TREATMENT_INTERACTIONS.tsv")
capture.output(sessionInfo(),file=file.path(p,"qc/R_SESSION_INFO.txt"))
print(as.data.frame(res[model=="median_cohort",.(dataset,endpoint,score,arm,n,high_n,HR,P,logrank_P,FDR)]))
print(as.data.frame(it))
writeLines("SUCCESS",file.path(p,"qc/MODELS_COMPLETE.txt"))
