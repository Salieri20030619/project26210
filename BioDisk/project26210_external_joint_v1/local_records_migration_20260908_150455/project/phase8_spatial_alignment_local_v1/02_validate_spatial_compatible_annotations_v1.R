#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing argument: ", flag)
  args[[idx + 1L]]
}

project_root <- normalizePath(arg_value("--project-root"), winslash = "/", mustWork = TRUE)
run_id <- arg_value("--run-id")
datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
phase7_root <- file.path(
  project_root, "phase7_annotation_local_v2", "runs", "20260810_234538_376",
  "objects_annotated"
)
run_root <- file.path(
  project_root, "phase8_spatial_alignment_local_v1", "annotation_spatial_crosswalk_v1",
  "runs", run_id
)
phase8_root <- file.path(run_root, "objects_annotated")
table_path <- file.path(run_root, "tables", "independent_output_validation.tsv")
status_path <- file.path(run_root, "PHASE8_SPATIAL_ANNOTATION_VALIDATION_STATUS")
if (file.exists(table_path) || file.exists(status_path)) stop("Refusing existing validation output")

expected_states <- c(
  "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
  "Epithelial-Normal", "F1-high Malignant", "F1-low Malignant", "Mast cell",
  "Monocyte", "NK cell", "Other Macrophage", "SMC/Pericyte",
  "SPP1-state Macrophage", "T cell"
)
phase7_columns <- c(
  "phase7_cluster_res0_8", "phase7_broad_annotation", "phase7_fine_annotation",
  "phase7_annotation_state", "phase7_annotation_confidence", "phase7_annotation_version"
)
phase8_columns <- c(
  "phase8_spatial_reference_state", "phase8_spatial_display_label",
  "phase8_spatial_mapping_status", "phase8_spatial_annotation_confidence",
  "phase8_spp1_program_score", "phase8_c1qc_program_score",
  "phase8_spp1_percentile_patient", "phase8_c1qc_percentile_patient",
  "phase8_f1_projection_score", "phase8_f1_percentile_patient",
  "phase8_f1_eligible_patient", "phase8_tymp_log_normalized",
  "phase8_malignancy_claim", "phase8_phase7_label_preserved",
  "phase8_annotation_version"
)

equal_with_na <- function(x, y) {
  identical(is.na(x), is.na(y)) && identical(as.character(x[!is.na(x)]), as.character(y[!is.na(y)]))
}

rows <- list()
for (dataset_current in datasets) {
  input_path <- file.path(phase7_root, paste0(dataset_current, "_annotated_local_v2.rds"))
  output_path <- file.path(phase8_root, paste0(dataset_current, "_spatial_annotated_v1.rds"))
  if (!file.exists(input_path) || !file.exists(output_path)) stop("Missing validation object")
  input_bytes_before <- file.info(input_path)$size
  input <- readRDS(input_path)
  output <- readRDS(output_path)
  input_meta <- input[[]]
  output_meta <- output[[]]
  missing_phase7 <- setdiff(phase7_columns, colnames(input_meta))
  missing_phase8 <- setdiff(phase8_columns, colnames(output_meta))
  phase7_equal <- !length(missing_phase7) && all(vapply(
    phase7_columns,
    function(column) equal_with_na(input_meta[[column]], output_meta[[column]]),
    logical(1)
  ))
  macrophage <- output_meta$phase7_fine_annotation %in% c("Macrophage_C1QC", "Macrophage_SPP1")
  tumor_like <- output_meta$phase7_fine_annotation == "Epithelial_tumor_like"
  valid_states <- all(
    is.na(output_meta$phase8_spatial_reference_state) |
      output_meta$phase8_spatial_reference_state %in% expected_states
  )
  checks <- c(
    inherits(input, "Seurat"),
    inherits(output, "Seurat"),
    identical(colnames(input), colnames(output)),
    ncol(input) == ncol(output),
    nrow(input) == nrow(output),
    !length(missing_phase8),
    phase7_equal,
    valid_states,
    !anyNA(output_meta$phase8_spatial_mapping_status),
    all(output_meta$phase8_phase7_label_preserved),
    !any(output_meta$phase8_malignancy_claim),
    !anyNA(output_meta$phase8_spp1_program_score[macrophage]),
    !anyNA(output_meta$phase8_c1qc_program_score[macrophage]),
    !anyNA(output_meta$phase8_spp1_percentile_patient[macrophage]),
    !anyNA(output_meta$phase8_c1qc_percentile_patient[macrophage]),
    !anyNA(output_meta$phase8_f1_projection_score[tumor_like]),
    file.info(input_path)$size == input_bytes_before
  )
  rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    status = if (all(checks)) "PASS" else "FAIL",
    cells_input = ncol(input),
    cells_output = ncol(output),
    features_input = nrow(input),
    features_output = nrow(output),
    cell_ids_identical = identical(colnames(input), colnames(output)),
    phase7_columns_identical = phase7_equal,
    phase8_required_columns_complete = !length(missing_phase8),
    phase8_mapping_status_complete = !anyNA(output_meta$phase8_spatial_mapping_status),
    valid_reference_states_only = valid_states,
    macrophage_program_complete = !anyNA(output_meta$phase8_spp1_program_score[macrophage]) &&
      !anyNA(output_meta$phase8_c1qc_program_score[macrophage]),
    tumor_like_f1_score_complete = !anyNA(output_meta$phase8_f1_projection_score[tumor_like]),
    malignancy_claims = sum(output_meta$phase8_malignancy_claim),
    input_bytes_before = input_bytes_before,
    input_bytes_after = file.info(input_path)$size,
    output_bytes = file.info(output_path)$size
  )
  rm(input, output, input_meta, output_meta)
  gc(verbose = FALSE)
}

validation <- rbindlist(rows)
if (any(validation$status != "PASS")) stop("Independent Phase 8 validation failed")
part <- paste0(table_path, ".part")
fwrite(validation, part, sep = "\t", quote = FALSE)
if (!file.rename(part, table_path)) stop("Validation table rename failed")

status_part <- paste0(status_path, ".part")
writeLines(c(
  "status=PASS",
  paste0("validated_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "datasets=5",
  paste0("cells=", sum(validation$cells_output)),
  "cell_ids_identical=5/5",
  "phase7_columns_identical=5/5",
  "phase8_required_columns_complete=5/5",
  "malignancy_claims=0"
), status_part, useBytes = TRUE)
if (!file.rename(status_part, status_path)) stop("Validation status rename failed")
