library(data.table);library(survival)
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase270_SPP1_literature_clinical_revalidation_20260927_v1"
d<-fread(file.path(p,"tables/JAVELIN_ANALYSIS_INPUT.tsv"))
out<-list();k<-0
for(cutmode in c("median_arm","median_cohort")){
 z<-copy(d[arm=="Avelumab+Axitinib"]);cut<-median(if(cutmode=="median_arm")z$SPP1 else d$SPP1);z[,x:=as.integer(SPP1>=cut)]
 for(eventmode in c("correct_1_minus_CNSR","INVALID_direct_CNSR")){
 z[,event:=if(eventmode=="correct_1_minus_CNSR")1-as.integer(PFS_P_CNSR) else as.integer(PFS_P_CNSR)]
 f<-coxph(Surv(PFS_P,event)~x,data=z);s<-summary(f)
 lr<-survdiff(Surv(PFS_P,event)~x,data=z);k<-k+1
 out[[k]]<-data.table(cutmode,eventmode,n=nrow(z),events=sum(z$event),HR=exp(coef(f)[1]),P=s$coefficients[1,"Pr(>|z|)"],logrank_P=pchisq(lr$chisq,1,lower.tail=FALSE),USE_FOR_CLINICAL_INFERENCE=eventmode=="correct_1_minus_CNSR")
 }}
fwrite(rbindlist(out),file.path(p,"qc/INVALID_EVENT_CODING_DIAGNOSTIC_ONLY.tsv"),sep="\t")
print(as.data.frame(rbindlist(out)))
