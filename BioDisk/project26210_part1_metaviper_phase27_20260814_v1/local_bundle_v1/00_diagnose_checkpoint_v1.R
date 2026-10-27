#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: 00_diagnose_checkpoint_v1.R COMBINED_ACTIVITY_RDS METADATA_RDS CHECKPOINT_RDS PATIENT_ID")
}

combined_path <- args[[1L]]
metadata_path <- args[[2L]]
checkpoint_path <- args[[3L]]
patient_id <- args[[4L]]
stopifnot(file.exists(combined_path), file.exists(metadata_path), file.exists(checkpoint_path))

activity <- readRDS(combined_path)
metadata <- as.data.frame(readRDS(metadata_path), stringsAsFactors = FALSE)
checkpoint <- readRDS(checkpoint_path)
patient_cells <- metadata$exact_cell_id[metadata$curated_patient_id == patient_id]

row_set_equal <- setequal(rownames(checkpoint), rownames(activity))
col_set_equal <- setequal(colnames(checkpoint), patient_cells)
aligned_diff <- Inf
if (row_set_equal && col_set_equal && all(is.finite(checkpoint))) {
  aligned <- checkpoint[rownames(activity), patient_cells, drop = FALSE]
  aligned_diff <- max(abs(aligned - activity[, patient_cells, drop = FALSE]))
}

out <- data.frame(
  patient = patient_id,
  checkpoint_rows = nrow(checkpoint),
  checkpoint_cols = ncol(checkpoint),
  combined_rows = nrow(activity),
  patient_cells = length(patient_cells),
  row_order_identical = identical(rownames(checkpoint), rownames(activity)),
  row_set_equal = row_set_equal,
  col_order_identical = identical(colnames(checkpoint), patient_cells),
  col_set_equal = col_set_equal,
  checkpoint_finite = all(is.finite(checkpoint)),
  aligned_max_abs_diff = aligned_diff,
  stringsAsFactors = FALSE
)
write.table(out, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
