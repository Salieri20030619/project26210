#!/usr/bin/env Rscript

.libPaths(c(
  "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/08_metaviper_v1/r_library",
  .libPaths()
))

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(aracne.networks)
  library(data.table)
  library(Matrix)
  library(org.Hs.eg.db)
  library(viper)
})

PROJECT <- "/mnt/BioDisk/ludingyu/project26210"
ROOT <- file.path(PROJECT, "results/part1_rebuild_v1")
INPUT <- file.path(ROOT, "07_pyscenic_v1/input")
OUT <- file.path(ROOT, "08_metaviper_v1")
TABLES <- file.path(OUT, "tables")
RDS <- file.path(OUT, "rds")
CHECKPOINTS <- file.path(RDS, "patient_activity")
GATE <- file.path(PROJECT, "check_memory_gate.sh")
MIN_PATIENT_CELLS <- 30L
MIN_REGULON_SIZE <- 25L
SEED <- 26210L

dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINTS, recursive = TRUE, showWarnings = FALSE)

gate_records <- list()
gate <- function(stage) {
  output <- system2(GATE, stage, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status"); if (is.null(status)) status <- 0L
  text <- paste(output, collapse = " ")
  message(text)
  gate_records[[length(gate_records) + 1L]] <<- data.table(stage=stage, exit_status=status, gate_output=text)
  fwrite(rbindlist(gate_records), file.path(TABLES, "memory_gate_log_v1.tsv"), sep="\t")
  if (status != 0L || !grepl("PASS", text, fixed=TRUE)) stop("Memory gate blocked ", stage)
}

gate("metaviper_load_input_v1")
meta <- fread(file.path(INPUT, "cell_metadata_v1.tsv.gz"), colClasses=c(exact_cell_id="character", curated_patient_id="character"))
genes <- fread(file.path(INPUT, "genes_v1.tsv"), colClasses="character")$gene
cells <- fread(file.path(INPUT, "cells_v1.tsv"), colClasses="character")$exact_cell_id
if (nrow(meta) != 11950L || length(genes) != 37795L || !identical(meta$exact_cell_id, cells)) stop("Frozen input mismatch")

con <- gzfile(file.path(INPUT, "raw_counts_genes_by_cells_v1.mtx.gz"), open="rb")
raw_counts <- readMM(con); close(con)
raw_counts <- as(raw_counts, "dgCMatrix")
dimnames(raw_counts) <- list(genes, cells)
if (!identical(dim(raw_counts), c(37795L,11950L)) || length(raw_counts@x) != 35713877L) stop("Raw matrix mismatch")
if (any(!is.finite(raw_counts@x)) || any(raw_counts@x < 0) || any(abs(raw_counts@x-round(raw_counts@x)) > 1e-8)) stop("Invalid counts")

patient_qc <- meta[, .(n_cells=.N, n_samples=uniqueN(sample_id), n_studies=uniqueN(study)), by=curated_patient_id]
patient_qc[, eligible := n_cells >= MIN_PATIENT_CELLS]
setorder(patient_qc, curated_patient_id)
fwrite(patient_qc, file.path(TABLES, "metaviper_patient_input_qc_v1.tsv"), sep="\t")
eligible_patients <- patient_qc[eligible == TRUE, curated_patient_id]
eligible_cells <- meta[curated_patient_id %chin% eligible_patients, exact_cell_id]
if (length(eligible_patients) < 10L || length(eligible_cells) < 10000L) stop("Unexpected eligible metaVIPER cohort")

gate("metaviper_map_and_normalize_v1")
mapped <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=genes, keytype="SYMBOL", column="ENTREZID", multiVals="first")
map_dt <- data.table(gene=genes, entrez_id=as.character(mapped[genes]))
map_dt[, mapped := !is.na(entrez_id)]
fwrite(map_dt, file.path(TABLES, "expression_symbol_to_entrez_v1.tsv"), sep="\t")
idx <- which(map_dt$mapped)
entrez <- unique(map_dt$entrez_id[idx])
aggregation <- sparseMatrix(i=match(map_dt$entrez_id[idx], entrez), j=seq_along(idx), x=1, dims=c(length(entrez),length(idx)))
counts_entrez <- aggregation %*% raw_counts[idx, eligible_cells, drop=FALSE]
dimnames(counts_entrez) <- list(entrez, eligible_cells)
lib_size <- Matrix::colSums(raw_counts[, eligible_cells, drop=FALSE])
log_cp10k <- counts_entrez %*% Diagonal(x=1e4/lib_size)
log_cp10k@x <- log1p(log_cp10k@x)
dimnames(log_cp10k) <- list(entrez, eligible_cells)

data(regulonkirc)
data(regulonkirp)
filter_regulon <- function(regulon, observed) {
  result <- lapply(regulon, function(entry) {
    keep <- names(entry$tfmode) %in% observed
    likelihood <- entry$likelihood
    if (length(likelihood) == length(keep)) likelihood <- likelihood[keep]
    list(tfmode=entry$tfmode[keep], likelihood=likelihood)
  })
  class(result) <- class(regulon)
  result
}
kirc <- filter_regulon(regulonkirc, rownames(log_cp10k))
kirp <- filter_regulon(regulonkirp, rownames(log_cp10k))
network_qc <- data.table(
  network=c("KIRC","KIRP"),
  n_regulators=c(length(kirc),length(kirp)),
  n_regulators_minsize25=c(
    sum(vapply(kirc,function(x) length(x$tfmode)>=MIN_REGULON_SIZE,logical(1))),
    sum(vapply(kirp,function(x) length(x$tfmode)>=MIN_REGULON_SIZE,logical(1)))
  )
)
fwrite(network_qc, file.path(TABLES, "kidney_interactome_qc_v1.tsv"), sep="\t")

activity_by_patient <- list(); run_qc <- list()
for (i in seq_along(eligible_patients)) {
  gate(paste0("metaviper_patient_", i, "_v1"))
  patient <- eligible_patients[i]
  patient_cells <- meta[curated_patient_id == patient & exact_cell_id %chin% eligible_cells, exact_cell_id]
  checkpoint <- file.path(CHECKPOINTS, paste0(gsub("[^A-Za-z0-9_.-]","_",patient),".rds"))
  started <- proc.time()[["elapsed"]]
  if (file.exists(checkpoint)) {
    activity <- readRDS(checkpoint); reused <- TRUE
  } else {
    set.seed(SEED+i)
    expression_dense <- as.matrix(log_cp10k[,patient_cells,drop=FALSE])
    activity <- viper::viper(expression_dense, list(KIRC=kirc,KIRP=kirp), method="rank", minsize=MIN_REGULON_SIZE, pleiotropy=FALSE, cores=2, verbose=FALSE)
    rm(expression_dense); gc(verbose=FALSE)
    saveRDS(activity, checkpoint, compress="gzip"); reused <- FALSE
  }
  if (!identical(colnames(activity),patient_cells) || any(!is.finite(activity))) stop("Invalid patient output: ",patient)
  activity_by_patient[[patient]] <- activity
  run_qc[[i]] <- data.table(curated_patient_id=patient,n_cells=length(patient_cells),n_regulators=nrow(activity),elapsed_seconds=proc.time()[["elapsed"]]-started,checkpoint_reused=reused)
  fwrite(rbindlist(run_qc), file.path(TABLES,"metaviper_patient_run_qc_v1.tsv"), sep="\t")
  message("[metaVIPER] ",i,"/",length(eligible_patients)," ",patient," cells=",length(patient_cells)," regulators=",nrow(activity))
}

regulators <- sort(unique(unlist(lapply(activity_by_patient,rownames))))
activity_matrix <- matrix(NA_real_,nrow=length(regulators),ncol=length(eligible_cells),dimnames=list(regulators,eligible_cells))
for (patient in names(activity_by_patient)) {
  x <- activity_by_patient[[patient]]
  activity_matrix[rownames(x),colnames(x)] <- x
}
if (anyNA(activity_matrix) || any(!is.finite(activity_matrix))) stop("Combined activity contains NA/nonfinite")
saveRDS(activity_matrix,file.path(RDS,"metaviper_activity_entrez_by_cell_v1.rds"),compress="gzip")
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(activity_matrix),keytype="ENTREZID",column="SYMBOL",multiVals="first")
fwrite(data.table(entrez_id=rownames(activity_matrix),symbol=as.character(symbols[rownames(activity_matrix)])),file.path(TABLES,"metaviper_regulator_id_map_v1.tsv"),sep="\t")
saveRDS(meta[exact_cell_id %chin% eligible_cells],file.path(RDS,"metaviper_cell_metadata_v1.rds"),compress="gzip")
writeLines(capture.output(sessionInfo()),file.path(OUT,"session_info_v1.txt"))
writeLines(format(Sys.time(),tz="",usetz=TRUE),file.path(OUT,"RUN_COMPLETE"))
message("PART1_METAVIPER_V1_PASS regulators=",nrow(activity_matrix)," cells=",ncol(activity_matrix))
