#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 01_audit_authoritative_rds_v1.R <project_root> <output_dir>")

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (dir.exists(output_dir) || file.exists(output_dir)) stop("Refusing existing output: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

safe_unique <- function(meta, candidates) {
  hit <- candidates[candidates %in% names(meta)]
  if (!length(hit)) return(NA_integer_)
  values <- as.character(meta[[hit[[1L]]]])
  values <- values[!is.na(values) & nzchar(trimws(values))]
  uniqueN(values)
}

score_rows <- list()
annotated_rows <- list()

for (dataset_current in datasets) {
  score_path <- file.path(
    project_root, "phase14_patient_program_validation_local_v1", "runs",
    "20260812_082000_provisional_v3", dataset_current, "objects",
    paste0(dataset_current, "_cell_scores_compact_v1.rds")
  )
  annotated_path <- file.path(
    project_root, "phase7_annotation_local_v2", "runs", "20260810_234538_376",
    "objects_annotated", paste0(dataset_current, "_annotated_local_v2.rds")
  )
  if (!file.exists(score_path) || !file.exists(annotated_path)) {
    stop("Missing authoritative RDS for ", dataset_current)
  }

  score_started <- Sys.time()
  score <- as.data.table(readRDS(score_path))
  score_finished <- Sys.time()
  required <- c("dataset", "cell_id", "patient_key", "sample_id", "phase7_fine_annotation")
  missing_required <- setdiff(required, names(score))
  if (length(missing_required)) stop(dataset_current, " score RDS missing: ", paste(missing_required, collapse = ";"))
  if (nrow(score) < 1L || anyDuplicated(score$cell_id) || !all(score$dataset == dataset_current)) {
    stop("Score RDS identity failure: ", dataset_current)
  }
  score_info <- file.info(score_path)
  score_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    rds_role = "phase14_compact_score",
    rds_path = normalizePath(score_path, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(score_info$size),
    md5 = unname(tools::md5sum(score_path)),
    object_class = paste(class(score), collapse = ";"),
    cells = nrow(score),
    features = NA_integer_,
    metadata_columns = ncol(score),
    patients = uniqueN(score$patient_key),
    samples = uniqueN(score$sample_id),
    annotations = uniqueN(score$phase7_fine_annotation),
    score_columns = length(grep("^(primary|sensitivity)__", names(score))),
    duplicate_cell_ids = anyDuplicated(score$cell_id),
    read_started = format(score_started, "%Y-%m-%dT%H:%M:%S%z"),
    read_finished = format(score_finished, "%Y-%m-%dT%H:%M:%S%z"),
    status = "PASS"
  )
  rm(score)
  invisible(gc())

  annotated_started <- Sys.time()
  annotated <- readRDS(annotated_path)
  annotated_finished <- Sys.time()
  if (!inherits(annotated, "Seurat")) stop("Annotated RDS is not Seurat: ", dataset_current)
  meta <- annotated@meta.data
  cell_ids <- colnames(annotated)
  if (length(cell_ids) < 1L || anyDuplicated(cell_ids) || nrow(meta) != length(cell_ids)) {
    stop("Annotated RDS identity failure: ", dataset_current)
  }
  annotation_col <- c("phase7_fine_annotation", "phase7_broad_annotation")
  annotation_col <- annotation_col[annotation_col %in% names(meta)]
  annotated_info <- file.info(annotated_path)
  annotated_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    rds_role = "phase7_annotated_seurat",
    rds_path = normalizePath(annotated_path, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(annotated_info$size),
    md5 = unname(tools::md5sum(annotated_path)),
    object_class = paste(class(annotated), collapse = ";"),
    cells = ncol(annotated),
    features = nrow(annotated),
    metadata_columns = ncol(meta),
    patients = safe_unique(meta, c("patient_key", "patient_id", "patient")),
    samples = safe_unique(meta, c("sample_id", "orig.ident", "sample")),
    annotations = if (length(annotation_col)) uniqueN(as.character(meta[[annotation_col[[1L]]]])) else NA_integer_,
    score_columns = length(grep("^(primary|sensitivity)__", names(meta))),
    duplicate_cell_ids = anyDuplicated(cell_ids),
    read_started = format(annotated_started, "%Y-%m-%dT%H:%M:%S%z"),
    read_finished = format(annotated_finished, "%Y-%m-%dT%H:%M:%S%z"),
    status = "PASS"
  )
  rm(meta, annotated)
  invisible(gc())
}

audit <- rbindlist(c(score_rows, annotated_rows), use.names = TRUE)
atomic_fwrite(audit, file.path(output_dir, "authoritative_rds_audit_v1.tsv"))

summary <- audit[, .(
  n_rds = .N,
  total_cells = sum(cells),
  all_pass = all(status == "PASS"),
  duplicate_cell_ids_total = sum(duplicate_cell_ids)
), by = rds_role]
atomic_fwrite(summary, file.path(output_dir, "authoritative_rds_audit_summary_v1.tsv"))

software <- data.table(
  item = c("R", "platform", "data.table"),
  version = c(paste(R.version$major, R.version$minor, sep = "."), R.version$platform, as.character(packageVersion("data.table")))
)
atomic_fwrite(software, file.path(output_dir, "software_versions_v1.tsv"))

status <- data.table(
  status = "PASS_ALL_AUTHORITATIVE_RDS_READ",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = length(datasets),
  rds_files = nrow(audit),
  all_readable = all(audit$status == "PASS"),
  note = "Five Phase 7 annotated Seurat RDS files and five Phase 14 compact score RDS files were read serially; no input was modified."
)
atomic_fwrite(status, file.path(output_dir, "status.tsv"))
cat("PASS_ALL_AUTHORITATIVE_RDS_READ\n")
