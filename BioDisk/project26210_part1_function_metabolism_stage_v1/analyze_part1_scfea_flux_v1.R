#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set.seed(26210)
project <- "/mnt/BioDisk/ludingyu/project26210"
root <- file.path(project,"results","part1_rebuild_v1")
out <- file.path(root,"09_metabolism_trajectory","metabolism_v1")
flux_file <- file.path(out,"source_data","scFEA_M168_cell_flux.tsv.gz")
meta_file <- file.path(out,"objects","metabolism_cell_metadata_four_auc_v1.rds")
stopifnot(file.exists(flux_file),file.exists(meta_file),file.exists(file.path(out,"AUC_COMPLETE")))
gate <- system2(file.path(project,"check_memory_gate.sh"),"part1_scfea_analysis_v1",stdout=file.path(out,"qc","memory_gate_scfea_analysis_v1.txt"),stderr=file.path(out,"qc","memory_gate_scfea_analysis_v1.txt"))
if(!identical(gate,0L)) stop("Memory gate refused scFEA analysis")
meta <- as.data.table(readRDS(meta_file)); cells <- meta$exact_cell_id
flux <- fread(flux_file); setnames(flux,1L,"exact_cell_id"); setkey(flux,exact_cell_id); flux <- flux[J(cells)]
stopifnot(identical(flux$exact_cell_id,cells),nrow(flux)==11950L,!anyNA(flux))
modules <- setdiff(names(flux),"exact_cell_id"); mat <- as.matrix(flux[,..modules]); storage.mode(mat)<-"double"
safe_z <- function(x) {s<-sd(x); if(!is.finite(s)||s==0) rep(0,length(x)) else (x-mean(x))/s}
axes <- list(Pyrimidine_flux=c("M_36","M_45","M_87",paste0("M_",150:166),"M_171"), CysMet_flux=c("M_17",paste0("M_",28:31),"M_70","M_80"), ArgPro_flux=c("M_42",paste0("M_",61:70),"M_100","M_102","M_103"), Choline_proxy=c("M_21","M_83"))
axes <- lapply(axes,intersect,y=modules); stopifnot(all(lengths(axes)>=2L))
zmat <- apply(mat,2,safe_z); scores <- as.data.table(vapply(axes,function(x) rowMeans(zmat[,x,drop=FALSE]),numeric(nrow(zmat))))
score_cols <- names(axes); meta <- cbind(meta,scores)
fwrite(meta[,c("exact_cell_id","study","curated_patient_id","TYMP_log1p_CP10K_v1","F1","clone_group","TYMP_group",score_cols),with=FALSE],file.path(out,"source_data","scFEA_four_axis_cell_scores_v1.tsv.gz"),sep="\t")
safe_wilcox <- function(x) if(length(x)<2L) NA_real_ else if(all(abs(x)<1e-12)) 1 else suppressWarnings(wilcox.test(x,mu=0,exact=FALSE)$p.value)
paired <- function(group_col,low,high,contrast,min_cells=20L){
  long<-melt(meta[!is.na(get(group_col))],id.vars=c("curated_patient_id",group_col),measure.vars=score_cols,variable.name="pathway",value.name="score")
  agg<-long[,.(median_score=median(score),n_cells=.N),by=c("curated_patient_id","pathway",group_col)]; setnames(agg,group_col,"group")
  wide<-dcast(agg[n_cells>=min_cells],curated_patient_id+pathway~group,value.var="median_score")
  if(!all(c(low,high)%in%names(wide))) return(list(detail=data.table(),stat=data.table()))
  detail<-wide[is.finite(get(low))&is.finite(get(high)),.(curated_patient_id,pathway,low_score=get(low),high_score=get(high),delta=get(high)-get(low))]
  stat<-detail[,.(contrast=contrast,n_patients=.N,median_patient_delta=median(delta),fraction_patients_positive=mean(delta>0),paired_wilcox_PValue=safe_wilcox(delta)),by=pathway]; stat[,FDR:=p.adjust(paired_wilcox_PValue,"BH")]
  list(detail=detail[,contrast:=contrast],stat=stat)
}
a<-paired("TYMP_group","TYMP_negative","TYMP_positive","TYMP-positive_vs_negative"); b<-paired("clone_group","other_primary_clones","F1_high_primary_clone","F1-high_clone_vs_other")
detail<-rbindlist(list(a$detail,b$detail),fill=TRUE); stat<-rbindlist(list(a$stat,b$stat),fill=TRUE)
fwrite(detail,file.path(out,"tables","scFEA_four_axis_patient_deltas_v1.tsv"),sep="\t"); fwrite(stat,file.path(out,"tables","scFEA_four_axis_paired_statistics_v1.tsv"),sep="\t")
safe_cor <- function(x,y){ok<-is.finite(x)&is.finite(y); if(sum(ok)<30||sd(x[ok])==0||sd(y[ok])==0) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman"))}
cors<-rbindlist(lapply(c("TYMP_log1p_CP10K_v1","F1"),function(pred) rbindlist(lapply(score_cols,function(sc) meta[,.(n_cells=.N,rho=safe_cor(get(sc),get(pred))),by=curated_patient_id][,`:=`(pathway=sc,predictor=pred)]))))[is.finite(rho)]
summ<-cors[,.(n_patients=.N,median_patient_rho=median(rho),fraction_positive=mean(rho>0),signed_rank_PValue=safe_wilcox(rho)),by=.(pathway,predictor)]; summ[,FDR:=p.adjust(signed_rank_PValue,"BH"),by=predictor]
fwrite(cors,file.path(out,"tables","scFEA_four_axis_correlations_by_patient_v1.tsv"),sep="\t"); fwrite(summ,file.path(out,"tables","scFEA_four_axis_correlation_summary_v1.tsv"),sep="\t")
definitions<-rbindlist(lapply(names(axes),function(n) data.table(pathway=n,module=axes[[n]]))); fwrite(definitions,file.path(out,"source_data","scFEA_four_axis_module_definitions_v1.tsv"),sep="\t")
limits<-data.table(requested_axis=c("MTA/polyamine","Glycerophospholipid/Kennedy","Matched metabolomics"),support=c("Indirect Cys/Met and Arg/Pro proxy; MTA is not explicit","Partial choline proxy M21/M83; not a complete Kennedy pathway","Unavailable in current inputs"),interpretation=c("Transcriptome-inferred proxy only","Do not label as full Kennedy flux","Cannot validate metabolite fold changes")); fwrite(limits,file.path(out,"tables","scFEA_scope_limitations_v1.tsv"),sep="\t")
stat[,label:=sprintf("%.3g%s",median_patient_delta,fifelse(FDR<0.05,"*",""))]
p<-ggplot(stat,aes(contrast,pathway,fill=median_patient_delta))+geom_tile(colour="white")+geom_text(aes(label=label),size=3)+scale_fill_gradient2(low="#355F9D",mid="white",high="#C94755",midpoint=0)+theme_minimal(base_size=9)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=25,hjust=1))+labs(x=NULL,y=NULL,fill="Median delta",title="Patient-paired scFEA pathway differences")
ggsave(file.path(out,"figures","scFEA_four_axis_paired_heatmap_v1.png"),p,width=6.8,height=3.8,dpi=300,bg="white"); ggsave(file.path(out,"figures","scFEA_four_axis_paired_heatmap_v1.pdf"),p,width=6.8,height=3.8)
metadata<-list(schema_version="part1_metabolism_scfea_v1",completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),n_cells=nrow(mat),n_modules=ncol(mat),n_axes=4,method="scFEA M168 masked implementation, 100 epochs",statistical_unit="patient-paired medians",limitations="MTA and Kennedy are proxy axes; no matched metabolomics",dependency_changes=FALSE)
jsonlite::write_json(metadata,file.path(out,"run_metadata_v1.json"),auto_unbox=TRUE,pretty=TRUE); writeLines(capture.output(sessionInfo()),file.path(out,"session_info_v1.txt"))
files<-list.files(out,recursive=TRUE,full.names=TRUE)
relative_files<-substring(files,nchar(out)+2)
files<-files[
  file.info(files)$isdir==FALSE &
  !basename(files)%in%c("RUN_COMPLETE","file_manifest_v1.tsv","background_status_v1.tsv") &
  !grepl("^logs/",relative_files)
]
sha<-function(x) strsplit(system2("sha256sum",x,stdout=TRUE),"[[:space:]]+")[[1]][1]
manifest<-data.table(file=substring(files,nchar(out)+2),size_bytes=file.info(files)$size,sha256=vapply(files,sha,character(1))); fwrite(manifest,file.path(out,"file_manifest_v1.tsv"),sep="\t")
writeLines(format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),file.path(out,"RUN_COMPLETE")); cat("PART1_SCFEA_ANALYSIS_V1_PASS\n")
