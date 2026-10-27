#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(data.table)
})
if (!requireNamespace("digest", quietly = TRUE)) stop("Required package digest is unavailable")

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
project_root <- if (length(args) >= 1L) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else {
  normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
}

input_dir <- file.path(
  project_root, "phase8_spatial_alignment_local_v1", "annotation_spatial_crosswalk_v1",
  "runs", "20260811_152500_full_v1", "objects_annotated"
)
stage_root <- file.path(project_root, "phase9_cnv_alignment_local_v1")
run_dir <- file.path(stage_root, "runs", "20260811_exact_cell_alignment_v1")
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

final_crosswalk <- file.path(table_dir, "phase8_all_cells_cnv_alignment_crosswalk_v1.tsv.gz")
part_crosswalk <- file.path(table_dir, "phase8_all_cells_cnv_alignment_crosswalk_v1.part.tsv.gz")
final_summary <- file.path(table_dir, "phase8_crosswalk_export_summary_v1.tsv")
part_summary <- file.path(table_dir, "phase8_crosswalk_export_summary_v1.part.tsv")
status_file <- file.path(run_dir, "PHASE9_CNV_CROSSWALK_EXPORT_STATUS")

stop_if_exists <- c(final_crosswalk, final_summary, status_file)
if (any(file.exists(stop_if_exists))) {
  stop("No-overwrite gate: one or more final outputs already exist")
}
unlink(c(part_crosswalk, part_summary), force = TRUE)

input_files <- sort(list.files(input_dir, pattern = "_spatial_annotated_v1[.]rds$", full.names = TRUE))
if (length(input_files) != 5L) stop("Expected exactly five Phase 8 input objects")

keep_core <- c(
  "cell_id", "barcode_original", "dataset", "object_name", "sample_id", "patient_id",
  "patient_key", "tissue_status", "disease_status", "tissue_confirmed",
  "phase7_broad_annotation", "phase7_fine_annotation", "phase7_annotation_state",
  "phase7_annotation_confidence", "phase8_spatial_reference_state",
  "phase8_spatial_display_label", "phase8_spatial_mapping_status",
  "phase8_spatial_annotation_confidence", "phase8_f1_projection_score",
  "phase8_f1_percentile_patient", "phase8_f1_eligible_patient",
  "phase8_malignancy_claim", "phase8_annotation_version"
)
keep_optional <- c(
  "src__sample_id", "src__gsm", "src__patient_id", "src__source_dataset",
  "src__input_file", "src__original_barcode"
)
keep <- c(keep_core, keep_optional)

pieces <- vector("list", length(input_files))
summary_rows <- vector("list", length(input_files))
for (i in seq_along(input_files)) {
  input_path <- input_files[[i]]
  message(sprintf("[%d/5] Reading %s", i, basename(input_path)))
  obj <- readRDS(input_path)
  md <- as.data.table(obj[[]], keep.rownames = "seurat_cell_name")
  missing_cols <- setdiff(keep_core, names(md))
  if (length(missing_cols)) stop("Missing required fields: ", paste(missing_cols, collapse = ", "))
  missing_optional <- setdiff(keep_optional, names(md))
  if (length(missing_optional)) {
    for (field in missing_optional) set(md, j = field, value = NA_character_)
  }
  out <- md[, c("seurat_cell_name", keep), with = FALSE]
  dataset_value <- unique(out$dataset)
  if (length(dataset_value) != 1L) stop("Dataset field is not constant within ", basename(input_path))
  if (nrow(out) != ncol(obj)) stop("Metadata/cell count mismatch in ", basename(input_path))
  summary_rows[[i]] <- data.table(
    dataset = dataset_value,
    input_file = normalizePath(input_path, winslash = "/", mustWork = TRUE),
    input_bytes = file.info(input_path)$size,
    n_cells = nrow(out),
    n_unique_cell_id = uniqueN(out$cell_id),
    n_missing_cell_id = sum(is.na(out$cell_id) | out$cell_id == ""),
    n_missing_barcode = sum(is.na(out$barcode_original) | out$barcode_original == ""),
    n_phase8_malignancy_claim_true = sum(out$phase8_malignancy_claim %in% TRUE, na.rm = TRUE)
  )
  pieces[[i]] <- out
  rm(obj, md, out)
  invisible(gc())
}

crosswalk <- rbindlist(pieces, use.names = TRUE, fill = FALSE)
summary_dt <- rbindlist(summary_rows)
if (nrow(crosswalk) != 585119L) stop("Expected 585119 cells, found ", nrow(crosswalk))
if (uniqueN(crosswalk$cell_id) != nrow(crosswalk)) stop("cell_id is not globally unique")
if (anyNA(crosswalk$cell_id) || any(crosswalk$cell_id == "")) stop("Missing cell_id detected")
if (anyNA(crosswalk$barcode_original) || any(crosswalk$barcode_original == "")) stop("Missing barcode detected")
if (any(crosswalk$phase8_malignancy_claim %in% TRUE, na.rm = TRUE)) stop("Unexpected Phase 8 malignancy claim")

setorder(crosswalk, dataset, object_name, cell_id)
fwrite(crosswalk, part_crosswalk, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
if (!file.rename(part_crosswalk, final_crosswalk)) stop("Atomic rename failed for crosswalk")

summary_dt <- rbind(
  summary_dt,
  data.table(
    dataset = "ALL",
    input_file = NA_character_,
    input_bytes = sum(summary_dt$input_bytes),
    n_cells = sum(summary_dt$n_cells),
    n_unique_cell_id = uniqueN(crosswalk$cell_id),
    n_missing_cell_id = sum(is.na(crosswalk$cell_id) | crosswalk$cell_id == ""),
    n_missing_barcode = sum(is.na(crosswalk$barcode_original) | crosswalk$barcode_original == ""),
    n_phase8_malignancy_claim_true = sum(crosswalk$phase8_malignancy_claim %in% TRUE, na.rm = TRUE)
  )
)
fwrite(summary_dt, part_summary, sep = "\t", quote = FALSE, na = "NA")
if (!file.rename(part_summary, final_summary)) stop("Atomic rename failed for summary")

sha <- digest::digest(file = final_crosswalk, algo = "sha256", serialize = FALSE)
status_lines <- c(
  "PASS",
  sprintf("completed_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("cells=%d", nrow(crosswalk)),
  sprintf("datasets=%d", uniqueN(crosswalk$dataset)),
  sprintf("unique_cell_ids=%d", uniqueN(crosswalk$cell_id)),
  sprintf("malignancy_claim_true=%d", sum(crosswalk$phase8_malignancy_claim %in% TRUE, na.rm = TRUE)),
  sprintf("crosswalk_sha256=%s", sha),
  sprintf("R=%s", R.version.string),
  sprintf("SeuratObject=%s", as.character(packageVersion("SeuratObject"))),
  sprintf("data.table=%s", as.character(packageVersion("data.table"))),
  sprintf("digest=%s", as.character(packageVersion("digest")))
)
writeLines(status_lines, paste0(status_file, ".part"), useBytes = TRUE)
if (!file.rename(paste0(status_file, ".part"), status_file)) stop("Atomic rename failed for status")
message("PASS: ", final_crosswalk)
