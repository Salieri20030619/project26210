#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=6L)stop("usage: 07_validate_m13_v1.R TRAJECTORY EXPRESSION_DIR METAVIPER_DIR M11_TARGET_STATUS ANALYSIS OUTPUT")
traj_path<-args[[1L]];expr_dir<-args[[2L]];mv_dir<-args[[3L]];target_path<-args[[4L]];analysis<-args[[5L]];output<-args[[6L]];part<-paste0(output,".part")
if(file.exists(output)||dir.exists(output)||file.exists(part)||dir.exists(part))stop("refusing existing output");dir.create(part,recursive=TRUE,showWarnings=FALSE)
checks<-list();add<-function(name,ok,detail)checks[[length(checks)+1L]]<<-data.table(check=name,status=ifelse(isTRUE(ok),"PASS","FAIL"),detail=as.character(detail))
st<-fromJSON(file.path(analysis,"status.json"));add("analysis_status",st$status=="PASS_M13_ANALYSIS_COMPLETE_WITH_SCOPE_LIMITATION",st$status)
traj<-fread(traj_path,colClasses=c(exact_cell_id="character",curated_patient_id="character"));expr<-fread(file.path(expr_dir,"m13_tymp_rrm2b_cell_expression_v1.tsv.gz"),colClasses=c(exact_cell_id="character"));activity<-readRDS(file.path(mv_dir,"rds","metaviper_activity_entrez_by_cell_v1.rds"));idx<-match(colnames(activity),traj$exact_cell_id);dmeta<-traj[idx];dexpr<-expr[idx]
add("source_scope",nrow(traj)==11950L&&nrow(expr)==11950L&&ncol(activity)==11921L&&uniqueN(dmeta$curated_patient_id)==20L&&identical(dmeta$exact_cell_id,colnames(activity)),paste(nrow(traj),ncol(activity)))
registry<-fread(file.path(analysis,"m13_feature_registry_v1.tsv"),colClasses=c(entrez_id="character"));add("feature_registry",nrow(registry)==7L&&setequal(registry$feature,c("CASP9","CHMP1A","PRSS8","PXN","HIF1A","TYMP_RNA","RRM2B_RNA")),paste(registry$feature,collapse=","))
methods<-fread(file.path(analysis,"m13_method_registry_v1.tsv"));add("method_registry",nrow(methods)==5L&&sum(methods$method_family=="primary_pseudotime")==3L,nrow(methods))
corr<-fread(file.path(analysis,"m13_patient_trajectory_correlations_v1.tsv"),na.strings="NA",colClasses=c(curated_patient_id="character"));terminal<-fread(file.path(analysis,"m13_patient_terminal_shifts_v1.tsv"),na.strings="NA",colClasses=c(curated_patient_id="character"))
add("correlation_rows",nrow(corr)==700L&&uniqueN(corr$curated_patient_id)==20L,nrow(corr));add("terminal_rows",nrow(terminal)==420L&&uniqueN(terminal$curated_patient_id)==20L,nrow(terminal))
prss<-registry[feature=="PRSS8"];score<-as.numeric(activity[prss$entrez_id,]);tv<-dmeta$scTour_pseudotime_patient_v1;sent<-data.table(pid=dmeta$curated_patient_id,score=score,tv=tv)[,.(effect_new=suppressWarnings(cor(score,tv,method="spearman"))),by=pid];obs<-corr[feature=="PRSS8"&method=="scTour",.(pid=curated_patient_id,effect_old=effect)];z<-merge(sent,obs,by="pid");add("source_sentinel_recomputed",nrow(z)==20L&&max(abs(z$effect_new-z$effect_old))<1e-12,max(abs(z$effect_new-z$effect_old)))
cs<-fread(file.path(analysis,"m13_trajectory_correlation_summary_v1.tsv"),na.strings="NA");ts<-fread(file.path(analysis,"m13_terminal_shift_summary_v1.tsv"),na.strings="NA");add("summary_rows",nrow(cs)==35L&&nrow(ts)==21L,paste(nrow(cs),nrow(ts)))
summary_ok<-TRUE
for(i in seq_len(nrow(cs))){r<-cs[i];v<-corr[feature==r$feature&method==r$method&is.finite(effect),effect];p<-if(all(v==0))1 else suppressWarnings(wilcox.test(v,mu=0,exact=NULL)$p.value);summary_ok<-summary_ok&&length(v)==r$n_eligible_patients&&abs(median(v)-r$median_patient_effect)<1e-12&&abs(p-r$wilcoxon_p_two_sided)<1e-12}
for(i in seq_len(nrow(ts))){r<-ts[i];v<-terminal[feature==r$feature&method==r$method&is.finite(effect),effect];p<-if(all(v==0))1 else suppressWarnings(wilcox.test(v,mu=0,exact=NULL)$p.value);summary_ok<-summary_ok&&length(v)==r$n_eligible_patients&&abs(median(v)-r$median_patient_effect)<1e-12&&abs(p-r$wilcoxon_p_two_sided)<1e-12}
add("summaries_recomputed",summary_ok,"counts medians Wilcoxon")
bins<-fread(file.path(analysis,"m13_patient_binned_curve_source_v1.tsv.gz"),colClasses=c(curated_patient_id="character"));add("binned_curve_keys",nrow(bins)>5000L&&!anyDuplicated(bins[,.(feature,method,curated_patient_id,bin)]),nrow(bins))
pairs<-fread(file.path(analysis,"m13_pairwise_terminal_coactivation_v1.tsv"));add("pairwise_rows",nrow(pairs)==63L,nrow(pairs));rr<-pairs[feature_a=="TYMP_RNA"&feature_b=="RRM2B_RNA"];add("requested_rna_pair",nrow(rr)==3L, nrow(rr))
add("scope_limitation_recorded",identical(st$TYMP_RRM2B_activity_available,FALSE)&&isTRUE(st$TYMP_RRM2B_RNA_sensitivity),"activity unavailable; RNA separate")
add("no_figures",length(list.files(analysis,pattern="[.](png|pdf)$",recursive=TRUE))==0L,"tables only")
check_dt<-rbindlist(checks);fwrite(check_dt,file.path(part,"verification_checks_v1.tsv"),sep="\t");failed<-check_dt[status!="PASS",check];out<-list(status=if(length(failed)==0L)"PASS_INDEPENDENT_VALIDATION" else "FAIL_INDEPENDENT_VALIDATION",checks_passed=sum(check_dt$status=="PASS"),checks_total=nrow(check_dt),failed_checks=failed);write(toJSON(out,auto_unbox=TRUE,pretty=TRUE),file.path(part,"status.json"));file.rename(part,output);if(length(failed))stop("validation failures: ",paste(failed,collapse=","))
