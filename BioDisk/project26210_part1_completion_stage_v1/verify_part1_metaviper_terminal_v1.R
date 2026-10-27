#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT <- "/mnt/BioDisk/ludingyu/project26210"
ROOT <- file.path(PROJECT, "results/part1_rebuild_v1")
INPUT <- file.path(ROOT, "07_pyscenic_v1/input")
MV <- file.path(ROOT, "08_metaviper_v1")
OUT <- file.path(MV, "terminal_verification_v1")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

gate <- system2(file.path(PROJECT,"check_memory_gate.sh"), "part1_metaviper_terminal_verification_v1", stdout=TRUE, stderr=TRUE)
gate_status <- attr(gate,"status"); if (is.null(gate_status)) gate_status <- 0L
writeLines(gate, file.path(OUT,"memory_gate_v1.txt"))
if (gate_status != 0L || !any(grepl("PASS",gate,fixed=TRUE))) stop("Memory gate failed")

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks)+1L]] <<- data.table(check=name,status=ifelse(isTRUE(passed),"PASS","FAIL"),detail=as.character(detail))
}

cells <- fread(file.path(INPUT,"cells_v1.tsv"), colClasses="character")[[1]]
meta <- fread(file.path(INPUT,"cell_metadata_v1.tsv.gz"), colClasses=c(exact_cell_id="character",curated_patient_id="character"))
patient_counts <- meta[,.(n_cells=.N),by=curated_patient_id][order(curated_patient_id)]
eligible_patients <- patient_counts[n_cells>=30L,curated_patient_id]
eligible_cells <- meta[curated_patient_id %chin% eligible_patients,exact_cell_id]
add_check("frozen_input_cells", length(cells)==11950L && uniqueN(cells)==11950L && identical(meta$exact_cell_id,cells), paste0("n=",length(cells),";unique=",uniqueN(cells)))
add_check("eligible_scope", length(eligible_patients)==20L && length(eligible_cells)==11921L, paste0("patients=",length(eligible_patients),";cells=",length(eligible_cells)))
add_check("run_complete_marker", file.exists(file.path(MV,"RUN_COMPLETE")), file.path(MV,"RUN_COMPLETE"))

activity <- readRDS(file.path(MV,"rds/metaviper_activity_entrez_by_cell_v1.rds"))
add_check("combined_activity_shape", identical(dim(activity),c(5566L,11921L)), paste(dim(activity),collapse="x"))
add_check("combined_activity_cell_order", identical(colnames(activity),eligible_cells), paste0("columns=",ncol(activity)))
add_check("combined_activity_unique_ids", uniqueN(rownames(activity))==nrow(activity) && uniqueN(colnames(activity))==ncol(activity), paste0("regulators=",uniqueN(rownames(activity)),";cells=",uniqueN(colnames(activity))))
add_check("combined_activity_finite", all(is.finite(activity)) && !anyNA(activity), paste0("finite=",sum(is.finite(activity)),";total=",length(activity)))
reg_sd <- apply(activity,1L,sd)
add_check("combined_activity_nonconstant", all(is.finite(reg_sd)) && all(reg_sd>0), paste0("nonconstant=",sum(reg_sd>0),"/",length(reg_sd)))

meta_rds <- as.data.table(readRDS(file.path(MV,"rds/metaviper_cell_metadata_v1.rds")))
add_check("metadata_rds_cell_order", identical(meta_rds$exact_cell_id,eligible_cells), paste0("rows=",nrow(meta_rds)))
add_check("metadata_rds_patient_scope", setequal(unique(meta_rds$curated_patient_id),eligible_patients), paste0("patients=",uniqueN(meta_rds$curated_patient_id)))

run_qc <- fread(file.path(MV,"tables/metaviper_patient_run_qc_v1.tsv"), colClasses=c(curated_patient_id="character"))
add_check("patient_run_qc_rows", nrow(run_qc)==20L && identical(run_qc$curated_patient_id,eligible_patients), paste0("rows=",nrow(run_qc)))
add_check("patient_run_qc_totals", sum(run_qc$n_cells)==11921L && all(run_qc$n_regulators==5566L), paste0("cells=",sum(run_qc$n_cells),";regulators=",paste(unique(run_qc$n_regulators),collapse=",")))

patient_verification <- list()
for (i in seq_along(eligible_patients)) {
  patient <- eligible_patients[[i]]
  patient_cells <- meta[curated_patient_id==patient,exact_cell_id]
  path <- file.path(MV,"rds/patient_activity",paste0(gsub("[^A-Za-z0-9_.-]","_",patient),".rds"))
  checkpoint <- readRDS(path)
  same <- identical(dim(checkpoint),c(5566L,length(patient_cells))) && setequal(rownames(checkpoint),rownames(activity)) && identical(colnames(checkpoint),patient_cells) && all(is.finite(checkpoint))
  max_diff <- if (same) max(abs(checkpoint-activity[rownames(checkpoint),patient_cells,drop=FALSE])) else Inf
  passed <- same && is.finite(max_diff) && max_diff==0
  add_check(paste0("patient_checkpoint_",i), passed, paste0(patient,";cells=",length(patient_cells),";regulator_set_equal=",setequal(rownames(checkpoint),rownames(activity)),";max_abs_diff=",format(max_diff,scientific=TRUE)))
  hash <- strsplit(system2("sha256sum",path,stdout=TRUE)," +")[[1]][1]
  patient_verification[[i]] <- data.table(curated_patient_id=patient,n_cells=length(patient_cells),n_regulators=nrow(checkpoint),max_abs_diff=max_diff,size_bytes=file.info(path)$size,sha256=hash,status=ifelse(passed,"PASS","FAIL"))
  rm(checkpoint); gc(verbose=FALSE)
}
fwrite(rbindlist(patient_verification),file.path(OUT,"patient_checkpoint_verification_v1.tsv"),sep="\t")

mapping <- fread(file.path(MV,"tables/metaviper_regulator_id_map_v1.tsv"),colClasses="character")
add_check("regulator_map_rows", nrow(mapping)==5566L && identical(mapping$entrez_id,rownames(activity)), paste0("rows=",nrow(mapping)))
add_check("regulator_symbol_coverage", sum(!is.na(mapping$symbol) & nzchar(mapping$symbol))>5000L, paste0("mapped_symbols=",sum(!is.na(mapping$symbol) & nzchar(mapping$symbol))))
network_qc <- fread(file.path(MV,"tables/kidney_interactome_qc_v1.tsv"))
add_check("kidney_network_qc", nrow(network_qc)==2L && setequal(network_qc$network,c("KIRC","KIRP")) && all(network_qc$n_regulators_minsize25>4000L), paste(capture.output(print(network_qc)),collapse=" "))

regulator_qc <- data.table(entrez_id=rownames(activity),mean_activity=rowMeans(activity),sd_activity=reg_sd,min_activity=apply(activity,1L,min),max_activity=apply(activity,1L,max))
regulator_qc <- merge(regulator_qc,mapping,by="entrez_id",all.x=TRUE,sort=FALSE)
fwrite(regulator_qc,file.path(OUT,"metaviper_regulator_activity_qc_v1.tsv"),sep="\t")

checks_dt <- rbindlist(checks)
fwrite(checks_dt,file.path(OUT,"verification_checks_v1.tsv"),sep="\t")
failed <- checks_dt[status!="PASS",check]
metadata_out <- list(
  schema_version="part1_metaviper_terminal_verification_v1",
  completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  n_checks=nrow(checks_dt),
  n_failed=length(failed),
  failed_checks=failed,
  n_regulators=nrow(activity),
  n_cells=ncol(activity),
  n_patients=length(eligible_patients),
  activity_min=min(activity),
  activity_max=max(activity),
  median_regulator_sd=median(reg_sd),
  interpretation="rank-based inferred regulator activity; not protein abundance or causal TF activity"
)
jsonlite::write_json(metadata_out,file.path(OUT,"verification_metadata_v1.json"),auto_unbox=TRUE,pretty=TRUE)
if (length(failed)>0L) stop("metaVIPER verification failures: ",paste(failed,collapse=","))
writeLines(metadata_out$completed_at,file.path(OUT,"RUN_COMPLETE"))
cat("PART1_METAVIPER_TERMINAL_VERIFICATION_V1_PASS checks=",nrow(checks_dt)," regulators=",nrow(activity)," cells=",ncol(activity),"\n",sep="")
