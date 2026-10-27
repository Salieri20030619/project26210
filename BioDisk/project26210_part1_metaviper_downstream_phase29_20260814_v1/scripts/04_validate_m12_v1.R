#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5L) stop("usage: 04_validate_m12_v1.R METADATA METAVIPER_DIR M11_TARGET_STATUS ANALYSIS OUTPUT")
metadata_path <- args[[1L]]; mv_dir <- args[[2L]]; target_path <- args[[3L]]; analysis <- args[[4L]]; output <- args[[5L]]
part <- paste0(output, ".part")
if (file.exists(output) || dir.exists(output) || file.exists(part) || dir.exists(part)) stop("refusing existing output")
dir.create(part, recursive=TRUE, showWarnings=FALSE)
checks <- list(); add <- function(name, ok, detail) checks[[length(checks)+1L]] <<- data.table(check=name,status=ifelse(isTRUE(ok),"PASS","FAIL"),detail=as.character(detail))

st <- fromJSON(file.path(analysis,"status.json")); add("analysis_status", st$status=="PASS_M12_ANALYSIS_COMPLETE", st$status)
all_meta <- fread(metadata_path, colClasses=c(exact_cell_id="character",curated_patient_id="character"))
activity <- readRDS(file.path(mv_dir,"rds","metaviper_activity_entrez_by_cell_v1.rds")); idx <- match(colnames(activity),all_meta$exact_cell_id); meta <- all_meta[idx]
add("source_scope", ncol(activity)==11921L && uniqueN(meta$curated_patient_id)==20L && identical(meta$exact_cell_id,colnames(activity)), paste(ncol(activity),uniqueN(meta$curated_patient_id)))
registry <- fread(file.path(analysis,"m12_feature_registry_v1.tsv"),colClasses=c(entrez_id="character"))
add("feature_registry", nrow(registry)==5L && setequal(registry$protein,c("CASP9","CHMP1A","PRSS8","PXN","HIF1A")),paste(registry$protein,collapse=","))
add("activity_rows", all(registry$entrez_id %chin% rownames(activity)), "all five present")
effects <- fread(file.path(analysis,"m12_patient_activity_effects_v1.tsv"),na.strings="NA",colClasses=c(entrez_id="character",curated_patient_id="character"))
add("effect_rows", nrow(effects)==300L && uniqueN(effects$curated_patient_id)==20L && uniqueN(effects$endpoint)==3L, nrow(effects))

sentinel <- registry[protein=="CASP9"]
vals <- as.numeric(activity[sentinel$entrez_id,]); d <- data.table(patient=meta$curated_patient_id,F1=as.numeric(meta$F1),dominant=meta$dominant_factor_rank7_v1,score=vals)
reproduced <- list()
for (pid in unique(d$patient)) {
  x <- d[patient==pid]
  rho <- suppressWarnings(cor(x$score,x$F1,method="spearman"))
  reproduced[[length(reproduced)+1L]] <- data.table(curated_patient_id=pid,effect=rho)
}
reproduced <- rbindlist(reproduced); observed <- effects[protein=="CASP9" & endpoint=="continuous_F1_spearman",.(curated_patient_id,effect)]
joined <- merge(reproduced,observed,by="curated_patient_id",suffixes=c("_new","_old"))
add("source_sentinel_recomputed", nrow(joined)==20L && max(abs(joined$effect_new-joined$effect_old))<1e-12, max(abs(joined$effect_new-joined$effect_old)))

summary <- fread(file.path(analysis,"m12_f1_activity_patient_summary_v1.tsv"),na.strings="NA",colClasses=c(entrez_id="character"))
add("summary_rows", nrow(summary)==15L, nrow(summary))
summary_ok <- TRUE
for (i in seq_len(nrow(summary))) {
  r <- summary[i]; v <- effects[protein==r$protein & endpoint==r$endpoint & is.finite(effect),effect]
  p <- if (length(v)>=3L) {if(all(v==0)) 1 else suppressWarnings(wilcox.test(v,mu=0,exact=NULL)$p.value)} else NA_real_
  summary_ok <- summary_ok && length(v)==r$n_eligible_patients && abs(median(v)-r$median_patient_effect)<1e-12 && abs(p-r$wilcoxon_p_two_sided)<1e-12
}
add("summaries_recomputed", summary_ok, "counts medians Wilcoxon")
add("dominant_scope_expected", all(summary[endpoint=="F1_dominant_median_difference",n_eligible_patients]==9L), paste(summary[endpoint=="F1_dominant_median_difference",n_eligible_patients],collapse=","))
add("no_figures", length(list.files(analysis,pattern="[.](png|pdf)$",recursive=TRUE))==0L, "tables only")
check_dt <- rbindlist(checks); fwrite(check_dt,file.path(part,"verification_checks_v1.tsv"),sep="\t")
failed <- check_dt[status!="PASS",check]; out <- list(status=if(length(failed)==0L)"PASS_INDEPENDENT_VALIDATION" else "FAIL_INDEPENDENT_VALIDATION",checks_passed=sum(check_dt$status=="PASS"),checks_total=nrow(check_dt),failed_checks=failed)
write(toJSON(out,auto_unbox=TRUE,pretty=TRUE),file.path(part,"status.json")); file.rename(part,output)
if(length(failed)) stop("validation failures: ",paste(failed,collapse=","))
