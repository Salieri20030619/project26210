#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: 01_validate_metaviper_aligned_v1.R CELLS_TSV METADATA_TSV_GZ METAVIPER_DIR EXPECTED_CELLS EXPECTED_PATIENTS OUTPUT_DIR")
}
cells_path <- args[[1L]]
metadata_path <- args[[2L]]
mv_dir <- args[[3L]]
expected_cells <- as.integer(args[[4L]])
expected_patients <- as.integer(args[[5L]])
output <- args[[6L]]
part <- paste0(output, ".part")

stopifnot(file.exists(cells_path), file.exists(metadata_path), dir.exists(mv_dir))
if (file.exists(output) || dir.exists(output)) stop("refusing existing output: ", output)
if (file.exists(part) || dir.exists(part)) stop("refusing stale partial output: ", part)
dir.create(part, recursive = TRUE, showWarnings = FALSE)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    detail = as.character(detail)
  )
}

cells <- fread(cells_path, colClasses = "character")[[1L]]
meta <- fread(metadata_path, colClasses = c(exact_cell_id = "character", curated_patient_id = "character"))
patient_counts <- meta[, .(n_cells = .N), by = curated_patient_id][order(curated_patient_id)]
eligible_patients <- patient_counts[n_cells >= 30L, curated_patient_id]
eligible_cells <- meta[curated_patient_id %chin% eligible_patients, exact_cell_id]
add_check("frozen_input_cells", length(cells) == 11950L && uniqueN(cells) == 11950L && identical(meta$exact_cell_id, cells), paste0("n=", length(cells)))
add_check("eligible_scope", length(eligible_patients) == expected_patients && length(eligible_cells) == expected_cells, paste0("patients=", length(eligible_patients), ";cells=", length(eligible_cells)))
add_check("upstream_run_complete", file.exists(file.path(mv_dir, "RUN_COMPLETE")), file.path(mv_dir, "RUN_COMPLETE"))

activity_path <- file.path(mv_dir, "rds", "metaviper_activity_entrez_by_cell_v1.rds")
activity <- readRDS(activity_path)
add_check("combined_activity_shape", identical(dim(activity), c(5566L, expected_cells)), paste(dim(activity), collapse = "x"))
add_check("combined_activity_cell_order", identical(colnames(activity), eligible_cells), paste0("columns=", ncol(activity)))
add_check("combined_activity_unique_keys", uniqueN(rownames(activity)) == nrow(activity) && uniqueN(colnames(activity)) == ncol(activity), "regulator and cell keys")
add_check("combined_activity_finite", all(is.finite(activity)) && !anyNA(activity), paste0("values=", length(activity)))
reg_sd <- apply(activity, 1L, sd)
add_check("combined_activity_nonconstant", all(is.finite(reg_sd)) && all(reg_sd > 0), paste0("nonconstant=", sum(reg_sd > 0)))

meta_rds <- as.data.table(readRDS(file.path(mv_dir, "rds", "metaviper_cell_metadata_v1.rds")))
add_check("metadata_rds_cell_order", identical(meta_rds$exact_cell_id, eligible_cells), paste0("rows=", nrow(meta_rds)))
add_check("metadata_rds_patient_scope", setequal(unique(meta_rds$curated_patient_id), eligible_patients), paste0("patients=", uniqueN(meta_rds$curated_patient_id)))

run_qc <- fread(file.path(mv_dir, "tables", "metaviper_patient_run_qc_v1.tsv"), colClasses = c(curated_patient_id = "character"))
add_check("patient_run_qc_rows", nrow(run_qc) == expected_patients && identical(run_qc$curated_patient_id, eligible_patients), paste0("rows=", nrow(run_qc)))
add_check("patient_run_qc_totals", sum(run_qc$n_cells) == expected_cells && all(run_qc$n_regulators == nrow(activity)), paste0("cells=", sum(run_qc$n_cells)))

patient_results <- list()
for (i in seq_along(eligible_patients)) {
  patient <- eligible_patients[[i]]
  patient_cells <- meta[curated_patient_id == patient, exact_cell_id]
  path <- file.path(mv_dir, "rds", "patient_activity", paste0(gsub("[^A-Za-z0-9_.-]", "_", patient), ".rds"))
  checkpoint <- readRDS(path)
  dimension_ok <- identical(dim(checkpoint), c(nrow(activity), length(patient_cells)))
  row_unique <- uniqueN(rownames(checkpoint)) == nrow(checkpoint)
  col_unique <- uniqueN(colnames(checkpoint)) == ncol(checkpoint)
  row_set_equal <- setequal(rownames(checkpoint), rownames(activity))
  col_set_equal <- setequal(colnames(checkpoint), patient_cells)
  finite <- all(is.finite(checkpoint)) && !anyNA(checkpoint)
  row_order_identical <- identical(rownames(checkpoint), rownames(activity))
  col_order_identical <- identical(colnames(checkpoint), patient_cells)
  max_diff <- Inf
  if (dimension_ok && row_unique && col_unique && row_set_equal && col_set_equal && finite) {
    aligned <- checkpoint[rownames(activity), patient_cells, drop = FALSE]
    max_diff <- max(abs(aligned - activity[, patient_cells, drop = FALSE]))
  }
  passed <- dimension_ok && row_unique && col_unique && row_set_equal && col_set_equal && finite && is.finite(max_diff) && max_diff == 0
  add_check(paste0("patient_checkpoint_", i), passed, paste0(patient, ";row_order_identical=", row_order_identical, ";col_order_identical=", col_order_identical, ";max_abs_diff=", format(max_diff, scientific = TRUE)))
  hash <- strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1L]][1L]
  patient_results[[i]] <- data.table(
    curated_patient_id = patient,
    n_cells = length(patient_cells),
    n_regulators = nrow(checkpoint),
    row_order_identical = row_order_identical,
    row_set_equal = row_set_equal,
    col_order_identical = col_order_identical,
    col_set_equal = col_set_equal,
    max_abs_diff_after_key_alignment = max_diff,
    sha256 = hash,
    status = ifelse(passed, "PASS", "FAIL")
  )
  rm(checkpoint)
  gc(verbose = FALSE)
}
patient_table <- rbindlist(patient_results)
fwrite(patient_table, file.path(part, "patient_checkpoint_key_aligned_verification_v1.tsv"), sep = "\t")

mapping <- fread(file.path(mv_dir, "tables", "metaviper_regulator_id_map_v1.tsv"), colClasses = "character")
add_check("regulator_map_rows", nrow(mapping) == nrow(activity) && identical(mapping$entrez_id, rownames(activity)), paste0("rows=", nrow(mapping)))
add_check("regulator_symbol_coverage", sum(!is.na(mapping$symbol) & nzchar(mapping$symbol)) > 5000L, paste0("mapped=", sum(!is.na(mapping$symbol) & nzchar(mapping$symbol))))
network_qc <- fread(file.path(mv_dir, "tables", "kidney_interactome_qc_v1.tsv"))
add_check("kidney_network_qc", nrow(network_qc) == 2L && setequal(network_qc$network, c("KIRC", "KIRP")) && all(network_qc$n_regulators_minsize25 > 4000L), "KIRC and KIRP")

checks_dt <- rbindlist(checks)
fwrite(checks_dt, file.path(part, "verification_checks_v2.tsv"), sep = "\t")
failed <- checks_dt[status != "PASS", check]
status <- list(
  status = if (length(failed) == 0L) "PASS_INDEPENDENT_VALIDATION" else "FAIL_INDEPENDENT_VALIDATION",
  checks_passed = sum(checks_dt$status == "PASS"),
  checks_total = nrow(checks_dt),
  failed_checks = failed,
  n_regulators = nrow(activity),
  n_cells = ncol(activity),
  n_patients = length(eligible_patients),
  patient_checkpoints_exact_after_key_alignment = sum(patient_table$max_abs_diff_after_key_alignment == 0),
  patient_checkpoints_original_row_order_identical = sum(patient_table$row_order_identical),
  repair = "compare checkpoint matrices after unique regulator-key and cell-key alignment; no activity values recomputed"
)
write_json(status, file.path(part, "status.json"), auto_unbox = TRUE, pretty = TRUE)

files <- list.files(part, full.names = TRUE)
manifest <- rbindlist(lapply(files, function(path) data.table(
  file = basename(path),
  bytes = file.info(path)$size,
  sha256 = strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1L]][1L]
)))
fwrite(manifest, file.path(part, "VALIDATION_MANIFEST_v1.tsv"), sep = "\t")
if (!file.rename(part, output)) stop("atomic output rename failed")
if (length(failed) > 0L) quit(status = 1L)

