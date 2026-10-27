suppressPackageStartupMessages(library(nlme));suppressPackageStartupMessages(library(jsonlite))
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase319_SPP1_direct_barrier_paper_reproduction_20261005_v1"
d<-read.delim(file.path(p,"tables/ALL_OBSERVATIONS.tsv"),check.names=FALSE)
x<-subset(d,panel=="Ext Fig10H");x$genotype<-factor(ifelse(x$group=="wt","WT","KO"),levels=c("WT","KO"));stopifnot(nrow(x)==13,!anyDuplicated(x$unit))
f1<-lm(value~genotype,data=x);a<-coef(summary(f1))[2,]
y<-subset(d,panel=="Ext Fig10C");y$genotype<-factor(y$group,levels=c("WT","KO"));y$mouse<-factor(y$unit);stopifnot(nlevels(y$mouse)==10)
f2<-lme(value~genotype,random=~1|mouse,data=y,method="REML",control=lmeControl(opt="optim"));b<-summary(f2)$tTable[2,]
res<-data.frame(panel=c("Ext Fig10H","Ext Fig10C"),endpoint=c("PDGFRbeta_CAF_count_per_tissue_weight","Ki67_positive_fraction_in_CD8"),n_WT=c(6,5),n_KO=c(7,5),estimate_KO_minus_WT=c(a[1],b[1]),SE=c(a[2],b[2]),df=c(df.residual(f1),b[3]),t=c(a[3],b[4]),P=c(a[4],b[5]),model=c("published_raw_linear_model","reconstructed_mouse_random_intercept_REML"))
res$FDR_two_model_audit<-p.adjust(res$P,method="BH");res$lower95<-res$estimate_KO_minus_WT-qt(.975,res$df)*res$SE;res$upper95<-res$estimate_KO_minus_WT+qt(.975,res$df)*res$SE
write.table(res,file.path(p,"tables/QUALIFIED_MODEL_CHECKS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write_json(list(CAF_coefficients=coef(f1),CAF_sigma=summary(f1)$sigma,CD8_variance=as.data.frame(VarCorr(f2)[,c("Variance","StdDev")]),CD8_logLik=as.numeric(logLik(f2)),CD8_fixed_cov=as.matrix(vcov(f2)),session=capture.output(sessionInfo())),file.path(p,"qc/MODEL_PROVENANCE.json"),pretty=TRUE,auto_unbox=TRUE)
saveRDS(list(CAF=f1,CD8=f2),file.path(p,"qc/RECONSTRUCTED_MODELS.rds"));print(as.data.frame(res));cat("CD8 Variance\n");print(VarCorr(f2))
