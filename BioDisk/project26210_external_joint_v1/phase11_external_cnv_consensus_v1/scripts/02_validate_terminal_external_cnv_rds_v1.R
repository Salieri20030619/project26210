#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(infercnv)
  library(jsonlite)
})

parse_args <- function(x) {
  if (length(x) %% 2 != 0L) stop("Arguments must be --key value pairs")
  out <- list()
  for (i in seq(1L, length(x), by = 2L)) out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
phase10_root <- normalizePath(args[["phase10-root"]], mustWork = TRUE)
out_dir <- args[["out-dir"]]
if (is.null(out_dir) || dir.exists(out_dir) || file.exists(out_dir)) {
  stop("No-overwrite gate: output directory missing or exists")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

atomic_json <- function(x, path) {
  part <- paste0(path, ".part")
  write_json(x, part, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

read_input_mapping <- function(method, unit) {
  path <- file.path(phase10_root, "full_inputs_v2", method, unit, "cell_id_mapping_v1.tsv.gz")
  mapping <- fread(path, colClasses = "character")
  if (anyDuplicated(mapping$tool_id) || anyDuplicated(mapping$exact_cell_id)) {
    stop("Non-unique input mapping: ", method, "/", unit)
  }
  mapping
}

records <- list()
record_index <- 0L

infer_root <- file.path(phase10_root, "results", "full_v2", "infercnv")
infer_dirs <- if (dir.exists(infer_root)) list.dirs(infer_root, recursive = FALSE, full.names = TRUE) else character()
infer_dirs <- sort(infer_dirs[file.exists(file.path(infer_dirs, "RUN_COMPLETE"))])
for (unit_dir in infer_dirs) {
  unit <- basename(unit_dir)
  mapping <- read_input_mapping("infercnv", unit)
  output_mapping <- fread(file.path(unit_dir, "cell_id_mapping_v2.tsv"), colClasses = "character")
  if (!identical(mapping$tool_id, output_mapping$tool_id) ||
      !identical(mapping$exact_cell_id, output_mapping$exact_cell_id)) {
    stop("inferCNV mapping order mismatch: ", unit)
  }
  status <- fromJSON(file.path(unit_dir, "status.json"), simplifyVector = TRUE)
  rds_path <- file.path(unit_dir, "infercnv_object_post_run_v2.rds")
  object <- readRDS(rds_path)
  cells <- colnames(object@expr.data)
  passed <- inherits(object, "infercnv") && length(cells) == nrow(mapping) &&
    setequal(cells, mapping$tool_id) && identical(status$status, "completed")
  record_index <- record_index + 1L
  records[[record_index]] <- data.table(
    method = "infercnv", unit_id = unit, terminal_state = "completed",
    n_input_cells = nrow(mapping), n_object_cells = length(cells),
    n_query = sum(mapping$role == "query"), n_reference = sum(mapping$role == "reference"),
    rds_path = rds_path, rds_size_bytes = file.info(rds_path)$size,
    object_class = paste(class(object), collapse = ";"), validation_status = ifelse(passed, "PASS", "FAIL")
  )
  rm(object, mapping, output_mapping)
  gc(verbose = FALSE)
}

copy_root <- file.path(phase10_root, "results", "full_v2", "copykat")
copy_dirs <- if (dir.exists(copy_root)) list.dirs(copy_root, recursive = FALSE, full.names = TRUE) else character()
copy_dirs <- sort(copy_dirs[file.exists(file.path(copy_dirs, "RUN_COMPLETE"))])
for (unit_dir in copy_dirs) {
  unit <- basename(unit_dir)
  mapping <- read_input_mapping("copykat", unit)
  status <- fromJSON(file.path(unit_dir, "status.json"), simplifyVector = TRUE)
  rds_path <- file.path(unit_dir, "copykat_result_v2.rds")
  object <- readRDS(rds_path)
  mapped <- fread(file.path(unit_dir, "copykat_predictions_mapped_v2.tsv"), colClasses = "character")
  prediction <- as.data.table(object$prediction)
  prediction_ids <- as.character(prediction$cell.names)
  passed <- is.list(object) && all(c("prediction", "CNAmat", "hclustering") %in% names(object)) &&
    nrow(prediction) == nrow(mapped) && setequal(prediction_ids, mapped$tool_id) &&
    all(mapped$exact_cell_id %in% mapping$exact_cell_id) && !anyDuplicated(mapped$exact_cell_id) &&
    identical(status$status, "completed")
  record_index <- record_index + 1L
  records[[record_index]] <- data.table(
    method = "copykat", unit_id = unit, terminal_state = "completed",
    n_input_cells = nrow(mapping), n_object_cells = nrow(prediction),
    n_query = sum(mapping$role == "query"), n_reference = sum(mapping$role == "reference"),
    rds_path = rds_path, rds_size_bytes = file.info(rds_path)$size,
    object_class = paste(class(object), collapse = ";"), validation_status = ifelse(passed, "PASS", "FAIL")
  )
  rm(object, mapping, mapped, prediction)
  gc(verbose = FALSE)
}

if (!length(records)) stop("No completed RDS outputs found")
result <- rbindlist(records, use.names = TRUE, fill = TRUE)
failed <- result[validation_status != "PASS"]
atomic_fwrite(result, file.path(out_dir, "terminal_rds_validation_v1.tsv"))
payload <- list(
  schema_version = "external_cnv_terminal_rds_validation_v1",
  status = if (nrow(failed)) "FAIL" else "PASS",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  phase10_root = phase10_root,
  n_infercnv_rds = sum(result$method == "infercnv"),
  n_copykat_rds = sum(result$method == "copykat"),
  n_passed = sum(result$validation_status == "PASS"),
  n_failed = nrow(failed),
  software = list(R = R.version.string, data.table = as.character(packageVersion("data.table")),
                  infercnv = as.character(packageVersion("infercnv")), jsonlite = as.character(packageVersion("jsonlite")))
)
atomic_json(payload, file.path(out_dir, "TERMINAL_RDS_VALIDATION_STATUS.json"))
marker <- if (nrow(failed)) "RUN_FAILED" else "RUN_COMPLETE"
part <- file.path(out_dir, paste0(marker, ".part"))
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), part)
if (!file.rename(part, file.path(out_dir, marker))) stop("Terminal marker rename failed")
if (nrow(failed)) stop("Terminal RDS validation failed for ", nrow(failed), " units")
cat("EXTERNAL_CNV_TERMINAL_RDS_VALIDATION_PASS infercnv=", payload$n_infercnv_rds,
    " copykat=", payload$n_copykat_rds, "\n", sep = "")
