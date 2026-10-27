library(data.table);library(survival)
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase271_C4_state_clinical_extension_20260928_v1"
old<-fread(file.path(dirname(p),"phase270_SPP1_literature_clinical_revalidation_20260927_v1/tables/ALL_SURVIVAL_MODELS.tsv"))
d<-fread(file.path(p,"tables/CHECKMATE_ANALYSIS_INPUT.tsv"),na.strings=c("NA",""));d[,':='(sex=factor(sex),risk=factor(risk),trial=factor(trial))]
out<-list()
for(ep in c("OS","PFS"))for(a in unique(d$arm))for(sc in c("SPP1","Jiang5")){
 z<-copy(d[arm==a]);z[,':='(x=as.numeric(scale(get(sc))),time=get(paste0(ep,"_time")),event=get(paste0(ep,"_event")))]
 vars<-c("x","time","event","age","sex","trial","risk");z<-z[complete.cases(z[,..vars])]
 cv<-c("age","sex","trial","risk");cv<-cv[vapply(cv,function(n)uniqueN(z[[n]])>1,logical(1))]
 fit<-coxph(as.formula(paste("Surv(time,event)~x+",paste(cv,collapse="+"))),data=z);ss<-summary(fit);ci<-ss$conf.int["x",];co<-ss$coefficients["x",]
 orig<-old[dataset=="CheckMate"&endpoint==ep&arm==a&score==sc&model=="continuous_adjusted"]
 out[[length(out)+1]]<-data.table(endpoint=ep,arm=a,score=sc,n=nrow(z),events=sum(z$event),HR=ci["exp(coef)"],CI_low=ci["lower .95"],CI_high=ci["upper .95"],P=co["Pr(>|z|)"],previous_HR=orig$HR,previous_P=orig$P)
}
q<-rbindlist(out);q[,FDR:=p.adjust(P,"BH"),by=endpoint];fwrite(q,file.path(p,"tables/PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT.tsv"),sep="\t");print(as.data.frame(q))
