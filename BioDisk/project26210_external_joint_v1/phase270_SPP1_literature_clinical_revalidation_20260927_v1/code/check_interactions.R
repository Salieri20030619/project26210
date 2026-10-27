library(data.table);library(survival)
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase270_SPP1_literature_clinical_revalidation_20260927_v1"
o<-list();j<-0
for(ds in c("CHECKMATE","JAVELIN")){
 all<-fread(file.path(p,"tables",paste0(ds,"_ANALYSIS_INPUT.tsv")),na.strings=c("","NA"))
 if(ds=="CHECKMATE")all<-all[Cohort=="CM-025"]
 all[,sex:=factor(sex)];if(ds=="CHECKMATE")all[,risk:=factor(risk)]
 all[,tx:=as.integer(if(ds=="CHECKMATE")arm=="NIVOLUMAB" else arm=="Avelumab+Axitinib")]
 for(ep in if(ds=="CHECKMATE")c("PFS","OS") else "PFS")for(sc in c("SPP1","Jiang5")){
 d<-copy(all);d[,':='(x=as.numeric(scale(get(sc))),time=get(paste0(ep,"_time")),event=get(paste0(ep,"_event")))]
 cv<-if(ds=="CHECKMATE")c("age","sex","risk") else c("age","sex");vars<-c("time","event","x","tx",cv)
 d<-d[complete.cases(d[,..vars])];f<-as.formula(paste("Surv(time,event)~x*tx+",paste(cv,collapse="+")))
 fit<-coxph(f,data=d,x=TRUE,y=TRUE);zph<-cox.zph(fit);v<-vcov(fit);bb<-coef(fit)
 for(term in names(bb)){j<-j+1;o[[j]]<-data.table(dataset=ds,endpoint=ep,score=sc,item=term,estimate=bb[term],SE=sqrt(v[term,term]),P=summary(fit)$coefficients[term,"Pr(>|z|)"],PH_P=if(term %in% rownames(zph$table))zph$table[term,"p"] else NA_real_,n=nrow(d))}
 for(arm in 0:1){
  lv<-bb["x"]+arm*bb["x:tx"];se<-sqrt(v["x","x"]+arm^2*v["x:tx","x:tx"]+2*arm*v["x","x:tx"])
  j<-j+1;o[[j]]<-data.table(dataset=ds,endpoint=ep,score=sc,item=paste0("score_effect_tx",arm),estimate=lv,SE=se,P=2*pnorm(-abs(lv/se)),PH_P=NA_real_,n=nrow(d))
 }
 saveRDS(fit,file.path(p,"qc",paste(ds,ep,sc,"interaction.rds",sep="_")))
}}
z<-rbindlist(o);z[,':='(HR=exp(estimate),CI_low=exp(estimate-1.96*SE),CI_high=exp(estimate+1.96*SE))]
fwrite(z,file.path(p,"tables/INTERACTION_DIAGNOSTICS.tsv"),sep="\t",na="NA")
print(as.data.frame(z[endpoint=="OS"&item %in% c("x:tx","score_effect_tx0","score_effect_tx1")]))
