suppressPackageStartupMessages(library(data.table))
P <- normalizePath('phase269_clinical_diagnosis_story_20260927_v1'); H <- dirname(getwd())
s <- fread(file.path(H,'project26210_phase83_bulk_state_projection_v1/results_v1/clinical_v1/phase83_sample_scores_v1.tsv'))
c <- fread(file.path(H,'project26210_external_joint_v1_phase64_large_immune_validation_v1/inputs_v1/CheckMate_Braun2020/processed_v1/checkmate_clinical_rna_subset_v1.tsv'))
s <- s[dataset=='CheckMate'];stopifnot(nrow(s)==311L)
d <- merge(s,c,by.x='sample_id',by.y='RNA_ID',sort=FALSE);stopifnot(nrow(d)==311L,!anyDuplicated(d$sample_id))
z <- function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
d[,spp1_z:=z(spp1_macrophage_fraction)];d[,bio:=tymp_z*spp1_z]
fwrite(d[Arm=='NIVOLUMAB',.N,by=.(Cohort,ORR)],file.path(P,'tables/RESPONSE_LABELS_BY_TRIAL.tsv'),sep='\t')
all <- list(); audit <- list(); memberships<-list();idx<-0L
for(mode in c('legacy','labels_only','labels_and_complete_cases')){
 labels<-if(mode=='legacy')c('CR','PR','SD','PD')else c('CR','PR','CRPR','SD','PD')
 a<-copy(d[Arm=='NIVOLUMAB'&ORR %in% labels]);a[,response:=as.integer(ORR %in% c('CR','PR','CRPR'))]
 if(mode=='labels_and_complete_cases')a<-a[!is.na(Cohort)&nzchar(trimws(Cohort))&!is.na(MSKCC)&nzchar(trimws(MSKCC))&is.finite(spp1_z)&is.finite(tymp_z)&is.finite(bio)]
 for(biological in c(FALSE,TRUE)){
  feature<-if(biological)'tymp_x_spp1_state' else 'spp1_macrophage_fraction'
  form<-if(biological)response~tymp_z+spp1_z+bio+factor(Cohort)+factor(MSKCC) else response~spp1_z+factor(Cohort)+factor(MSKCC)
  warnings<-character();fit<-withCallingHandlers(glm(form,data=a,family=binomial(),model=TRUE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  term<-if(biological)'bio' else 'spp1_z';tab<-summary(fit)$coef;beta<-tab[term,'Estimate'];se<-tab[term,'Std. Error'];p<-tab[term,'Pr(>|z|)'];idx<-idx+1L
  all[[idx]]<-data.table(mode=mode,feature=feature,input_n=nrow(a),actual_model_n=nobs(fit),responders=sum(model.response(model.frame(fit))),OR=exp(beta),CI_low=exp(beta-1.96*se),CI_high=exp(beta+1.96*se),p_value=p,converged=fit$converged,iterations=fit$iter,rank=fit$rank,coefficients=length(coef(fit)),max_abs_coefficient=max(abs(coef(fit)),na.rm=TRUE),warnings=paste(unique(warnings),collapse='; '))
  co<-as.data.table(tab,keep.rownames='term');co[,`:=`(mode=mode,feature=feature)];audit[[idx]]<-co
  members<-a[,.(sample_id,SUBJID,Cohort,MSKCC,ORR,response)];members[,`:=`(mode=mode,feature=feature)];memberships[[idx]]<-members
 }
}
a<-rbindlist(all);a[,FDR_ORR2:=p.adjust(p_value,'BH'),by=mode]
legacy<-fread(file.path(H,'project26210_phase83_bulk_state_projection_v1/results_v1/clinical_v1/phase83_secondary_models_v1.tsv'))
checks<-merge(a[mode=='legacy'],legacy[endpoint=='ORR'],by='feature',suffixes=c('_reproduced','_archived'))
checks[,`:=`(OR_abs_error=abs(OR-effect),P_abs_error=abs(p_value_reproduced-p_value_archived))];stopifnot(max(checks$OR_abs_error)<1e-8,max(checks$P_abs_error)<1e-8)
corrected<-copy(legacy)
for(feature_name in a[mode=='labels_and_complete_cases',feature]){
 row<-a[mode=='labels_and_complete_cases'&feature==feature_name]
 corrected[endpoint=='ORR'&feature==feature_name,`:=`(n=row$actual_model_n,n_event=row$responders,effect=row$OR,ci_low=row$CI_low,ci_high=row$CI_high,p_value=row$p_value)]
}
corrected[,q_value_historical12_repaired:=p.adjust(p_value,'BH')]
fwrite(a,file.path(P,'tables/ORR_REPAIR_COMPARISON.tsv'),sep='\t');fwrite(rbindlist(audit),file.path(P,'tables/ORR_ALL_COEFFICIENTS.tsv'),sep='\t')
fwrite(rbindlist(memberships),file.path(P,'tables/ORR_MODEL_INPUT_MEMBERSHIP.tsv'),sep='\t');fwrite(checks,file.path(P,'qc/LEGACY_ORR_REPRODUCTION.tsv'),sep='\t')
fwrite(corrected,file.path(P,'tables/PHASE83_SECONDARY_AUDIT_WITH_CORRECTED_ORR.tsv'),sep='\t')
print(as.data.frame(a));writeLines(capture.output(sessionInfo()),file.path(P,'qc/R_SESSION_INFO.txt'))
writeLines('ORR_REPAIR_EXIT_STATUS=0',file.path(P,'qc/ORR_REPAIR_COMPLETE.txt'))
