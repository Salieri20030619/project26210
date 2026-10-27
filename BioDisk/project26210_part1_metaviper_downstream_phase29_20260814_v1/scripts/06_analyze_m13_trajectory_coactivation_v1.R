#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=5L) stop("usage: 06_analyze_m13_trajectory_coactivation_v1.R TRAJECTORY_SCORES EXPRESSION_DIR METAVIPER_DIR M11_TARGET_STATUS OUTPUT")
traj_path<-args[[1L]]; expr_dir<-args[[2L]]; mv_dir<-args[[3L]]; target_path<-args[[4L]]; output<-args[[5L]]; part<-paste0(output,".part")
stopifnot(file.exists(traj_path),dir.exists(expr_dir),dir.exists(mv_dir),file.exists(target_path))
if(file.exists(output)||dir.exists(output)||file.exists(part)||dir.exists(part))stop("refusing existing output")
dir.create(part,recursive=TRUE,showWarnings=FALSE)

traj<-fread(traj_path,colClasses=c(exact_cell_id="character",curated_patient_id="character"))
expr<-fread(file.path(expr_dir,"m13_tymp_rrm2b_cell_expression_v1.tsv.gz"),colClasses=c(exact_cell_id="character"))
activity<-readRDS(file.path(mv_dir,"rds","metaviper_activity_entrez_by_cell_v1.rds")); mv_meta<-as.data.table(readRDS(file.path(mv_dir,"rds","metaviper_cell_metadata_v1.rds")))
map<-fread(file.path(mv_dir,"tables","metaviper_regulator_id_map_v1.tsv"),colClasses="character")
if(nrow(traj)!=11950L||uniqueN(traj$exact_cell_id)!=11950L||nrow(expr)!=11950L||!identical(expr$exact_cell_id,traj$exact_cell_id))stop("trajectory/expression gate failed")
idx<-match(colnames(activity),traj$exact_cell_id); if(anyNA(idx)||ncol(activity)!=11921L||!identical(mv_meta$exact_cell_id,colnames(activity)))stop("common cell gate failed")
dmeta<-copy(traj[idx]); dexpr<-copy(expr[idx]); if(!identical(dmeta$exact_cell_id,colnames(activity))||uniqueN(dmeta$curated_patient_id)!=20L)stop("common metadata gate failed")

targets<-fread(target_path,colClasses=c(entrez_id="character"))[activity_row_present==TRUE,.(feature=protein,entrez_id,evidence_family="final_MR_inferable",feature_type="metaVIPER_activity")]
hif<-map[symbol=="HIF1A",.(feature="HIF1A",entrez_id,evidence_family="contextual_HIF1A",feature_type="metaVIPER_activity")]
rna<-data.table(feature=c("TYMP_RNA","RRM2B_RNA"),entrez_id=NA_character_,evidence_family="RNA_sensitivity",feature_type="RNA_log1p_cp10k")
features<-rbindlist(list(targets,hif,rna),fill=TRUE); if(nrow(features)!=7L)stop("feature gate failed")
fwrite(features,file.path(part,"m13_feature_registry_v1.tsv"),sep="\t",na="NA")
methods<-data.table(
  method=c("Monocle3","Slingshot","scTour","CytoTRACE2_potency","CellRank_TYMP_fate"),
  column=c("Monocle3_pseudotime_patient_v1","Slingshot_pseudotime_patient_v1","scTour_pseudotime_patient_v1","CytoTRACE2_Score","CellRank_fate_TYMP_state_v1"),
  method_family=c("primary_pseudotime","primary_pseudotime","primary_pseudotime","potency_sensitivity","fate_probability_sensitivity")
)
if(!all(methods$column%chin%names(dmeta)))stop("trajectory method columns missing")
fwrite(methods,file.path(part,"m13_method_registry_v1.tsv"),sep="\t")

feature_values<-list()
for(i in seq_len(nrow(features))){
  f<-features$feature[[i]]
  value<-if(features$feature_type[[i]]=="metaVIPER_activity") as.numeric(activity[features$entrez_id[[i]],]) else if(f=="TYMP_RNA") dexpr$TYMP_log1p_cp10k else dexpr$RRM2B_log1p_cp10k
  feature_values[[f]]<-value
}

corr_rows<-list(); terminal_rows<-list(); bin_rows<-list()
for(mi in seq_len(nrow(methods))){
  method<-methods$method[[mi]]; tvalue<-as.numeric(dmeta[[methods$column[[mi]]]])
  if(any(!is.finite(tvalue)))stop("nonfinite trajectory values")
  for(fi in seq_len(nrow(features))){
    f<-features$feature[[fi]]; score<-feature_values[[f]]
    d<-data.table(curated_patient_id=dmeta$curated_patient_id,trajectory=tvalue,score=score)
    for(pid in unique(d$curated_patient_id)){
      x<-d[curated_patient_id==pid]
      rho<-if(nrow(x)>=20L&&uniqueN(x$trajectory)>1L&&uniqueN(x$score)>1L)suppressWarnings(cor(x$score,x$trajectory,method="spearman"))else NA_real_
      corr_rows[[length(corr_rows)+1L]]<-data.table(feature=f,evidence_family=features$evidence_family[[fi]],feature_type=features$feature_type[[fi]],method=method,method_family=methods$method_family[[mi]],curated_patient_id=pid,n_cells=nrow(x),effect=rho)
      if(methods$method_family[[mi]]=="primary_pseudotime"){
        q20<-as.numeric(quantile(x$trajectory,.2,names=FALSE)); q80<-as.numeric(quantile(x$trajectory,.8,names=FALSE)); early<-x$trajectory<=q20; late<-x$trajectory>=q80
        delta<-if(q80>q20&&sum(early)>=5L&&sum(late)>=5L)median(x$score[late])-median(x$score[early])else NA_real_
        terminal_rows[[length(terminal_rows)+1L]]<-data.table(feature=f,evidence_family=features$evidence_family[[fi]],feature_type=features$feature_type[[fi]],method=method,curated_patient_id=pid,n_early=sum(early),n_late=sum(late),effect=delta)
        ranks<-frank(x$trajectory,ties.method="average")/nrow(x); bins<-pmax(1L,pmin(20L,ceiling(ranks*20)))
        b<-data.table(bin=bins,score=x$score)[,.(n_cells=.N,mean_score=mean(score),median_score=median(score)),by=bin]
        b[,`:=`(feature=f,evidence_family=features$evidence_family[[fi]],feature_type=features$feature_type[[fi]],method=method,curated_patient_id=pid,bin_midpoint=(bin-.5)/20)]
        bin_rows[[length(bin_rows)+1L]]<-b
      }
    }
  }
}
corr<-rbindlist(corr_rows); terminal<-rbindlist(terminal_rows); bins<-rbindlist(bin_rows)
fwrite(corr,file.path(part,"m13_patient_trajectory_correlations_v1.tsv"),sep="\t",na="NA")
fwrite(terminal,file.path(part,"m13_patient_terminal_shifts_v1.tsv"),sep="\t",na="NA")
fwrite(bins,file.path(part,"m13_patient_binned_curve_source_v1.tsv.gz"),sep="\t",compress="gzip")

summarize<-function(v){v<-v[is.finite(v)];p<-if(length(v)>=3L){if(all(v==0))1 else suppressWarnings(wilcox.test(v,mu=0,exact=NULL)$p.value)}else NA_real_;list(n_eligible_patients=length(v),n_positive=sum(v>0),n_negative=sum(v<0),n_zero=sum(v==0),median_patient_effect=if(length(v))median(v)else NA_real_,wilcoxon_p_two_sided=p)}
corr_summary<-corr[,summarize(effect),by=.(feature,evidence_family,feature_type,method,method_family)]
corr_summary[,BH_FDR:=p.adjust(wilcoxon_p_two_sided,method="BH"),by=.(evidence_family,method)]
corr_summary[,positive_support:=median_patient_effect>0&n_eligible_patients>=10L&BH_FDR<.05]
terminal_summary<-terminal[,summarize(effect),by=.(feature,evidence_family,feature_type,method)]
terminal_summary[,BH_FDR:=p.adjust(wilcoxon_p_two_sided,method="BH"),by=.(evidence_family,method)]
terminal_summary[,positive_support:=median_patient_effect>0&n_eligible_patients>=10L&BH_FDR<.05]
fwrite(corr_summary,file.path(part,"m13_trajectory_correlation_summary_v1.tsv"),sep="\t",na="NA")
fwrite(terminal_summary,file.path(part,"m13_terminal_shift_summary_v1.tsv"),sep="\t",na="NA")

curve_summary<-bins[,.(n_patients=uniqueN(curated_patient_id),patient_equal_mean=mean(mean_score),patient_median=median(mean_score)),by=.(feature,evidence_family,feature_type,method,bin,bin_midpoint)]
fwrite(curve_summary,file.path(part,"m13_patient_balanced_curve_summary_v1.tsv"),sep="\t")

pairs<-t(combn(features$feature,2)); pair_rows<-list()
for(method_name in methods[method_family=="primary_pseudotime",method]){
  d<-terminal[method==method_name]
  for(i in seq_len(nrow(pairs))){
    a<-d[feature==pairs[i,1],.(curated_patient_id,delta_a=effect)]; b<-d[feature==pairs[i,2],.(curated_patient_id,delta_b=effect)]; z<-merge(a,b,by="curated_patient_id"); z<-z[is.finite(delta_a)&is.finite(delta_b)]
    pair_rows[[length(pair_rows)+1L]]<-data.table(method=method_name,feature_a=pairs[i,1],feature_b=pairs[i,2],n_patients=nrow(z),n_both_positive=sum(z$delta_a>0&z$delta_b>0),fraction_both_positive=mean(z$delta_a>0&z$delta_b>0),median_joint_min_shift=median(pmin(z$delta_a,z$delta_b)),spearman_patient_shift=suppressWarnings(cor(z$delta_a,z$delta_b,method="spearman")))
  }
}
pair_dt<-rbindlist(pair_rows); fwrite(pair_dt,file.path(part,"m13_pairwise_terminal_coactivation_v1.tsv"),sep="\t",na="NA")
status<-list(status="PASS_M13_ANALYSIS_COMPLETE_WITH_SCOPE_LIMITATION",n_cells=ncol(activity),n_patients=uniqueN(dmeta$curated_patient_id),n_features=nrow(features),primary_pseudotime_methods=3,sensitivity_methods=2,TYMP_RRM2B_activity_available=FALSE,TYMP_RRM2B_RNA_sensitivity=TRUE,figures_created=FALSE)
write(toJSON(status,auto_unbox=TRUE,pretty=TRUE),file.path(part,"status.json"));file.rename(part,output)
