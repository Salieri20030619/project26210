#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_extract_part1_rds_report_data_v1.R <project_root> <output_dir>")
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (dir.exists(output_dir) || file.exists(output_dir)) stop("Refusing existing output: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
score_run <- file.path(
  project_root, "phase14_patient_program_validation_local_v1", "runs",
  "20260812_082000_provisional_v3"
)
final_run <- file.path(
  project_root, "phase14_patient_program_validation_local_v1", "runs",
  "20260812_final_cnv_gated_v6"
)
analysis_run <- file.path(
  project_root, "phase14_patient_program_validation_local_v1",
  "patient_analysis_final_cnv_v6"
)

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

clean_level <- function(x) {
  y <- as.character(x)
  y[is.na(y) | trimws(y) == ""] <- "<missing>"
  y
}

audit_rows <- list()
overview_rows <- list()
compartment_rows <- list()
fine_rows <- list()
state_rows <- list()
clinical_rows <- list()
patient_rows <- list()
score_rows <- list()

for (dataset_current in datasets) {
  rds_path <- file.path(
    score_run, dataset_current, "objects",
    paste0(dataset_current, "_cell_scores_compact_v1.rds")
  )
  if (!file.exists(rds_path)) stop("Missing RDS: ", rds_path)
  started <- Sys.time()
  cells <- as.data.table(readRDS(rds_path))
  finished <- Sys.time()

  required <- c(
    "dataset", "cell_id", "patient_key", "patient_id", "sample_id", "object_name",
    "timepoint", "treatment", "response", "phase7_broad_annotation",
    "phase7_fine_annotation", "phase8_spatial_reference_state",
    "phase8_spatial_display_label", "analysis_compartment"
  )
  missing_required <- setdiff(required, names(cells))
  if (length(missing_required)) stop(dataset_current, " missing fields: ", paste(missing_required, collapse = ";"))
  if (nrow(cells) == 0L || anyDuplicated(cells$cell_id) || !all(cells$dataset == dataset_current)) {
    stop("RDS identity failure: ", dataset_current)
  }
  score_cols <- grep("^(primary|sensitivity)__", names(cells), value = TRUE)
  if (length(score_cols) == 0L) stop("No program score columns: ", dataset_current)

  info <- file.info(rds_path)
  audit_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    rds_path = normalizePath(rds_path, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(info$size),
    md5 = unname(tools::md5sum(rds_path)),
    object_class = paste(class(cells), collapse = ";"),
    rows = nrow(cells),
    columns = ncol(cells),
    score_columns = length(score_cols),
    duplicate_cell_ids = anyDuplicated(cells$cell_id),
    read_started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    read_finished = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
    status = "PASS"
  )

  overview_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    cells_in_compact_rds = nrow(cells),
    patients = uniqueN(cells$patient_key),
    samples = uniqueN(cells$sample_id),
    objects = uniqueN(cells$object_name),
    broad_annotations = uniqueN(cells$phase7_broad_annotation),
    fine_annotations = uniqueN(cells$phase7_fine_annotation),
    spatial_reference_states = uniqueN(cells$phase8_spatial_reference_state),
    analysis_compartments = uniqueN(cells$analysis_compartment),
    missing_patient_key_cells = sum(is.na(cells$patient_key) | trimws(as.character(cells$patient_key)) == ""),
    missing_sample_id_cells = sum(is.na(cells$sample_id) | trimws(as.character(cells$sample_id)) == "")
  )

  compartment_rows[[dataset_current]] <- cells[, .(n_cells = .N), by = .(
    dataset, analysis_compartment = clean_level(analysis_compartment)
  )][, fraction_dataset_cells := n_cells / sum(n_cells), by = dataset]

  fine_rows[[dataset_current]] <- cells[, .(n_cells = .N), by = .(
    dataset,
    broad_annotation = clean_level(phase7_broad_annotation),
    fine_annotation = clean_level(phase7_fine_annotation)
  )][, fraction_dataset_cells := n_cells / sum(n_cells), by = dataset]

  state_rows[[dataset_current]] <- cells[, .(n_cells = .N), by = .(
    dataset,
    spatial_reference_state = clean_level(phase8_spatial_reference_state),
    spatial_display_label = clean_level(phase8_spatial_display_label)
  )][, fraction_dataset_cells := n_cells / sum(n_cells), by = dataset]

  for (field in c("timepoint", "treatment", "response", "tissue_confirmed")) {
    values <- clean_level(cells[[field]])
    clinical_rows[[paste(dataset_current, field, sep = "||")]] <- data.table(
      dataset = dataset_current, field = field, value = values
    )[, .(n_cells = .N), by = .(dataset, field, value)]
  }

  patient_rows[[dataset_current]] <- cells[, .(
    n_cells = .N,
    n_samples = uniqueN(sample_id),
    n_objects = uniqueN(object_name),
    n_compartments = uniqueN(analysis_compartment)
  ), by = .(dataset, patient_key, patient_id)]

  dataset_scores <- vector("list", length(score_cols))
  for (i in seq_along(score_cols)) {
    score_col <- score_cols[[i]]
    method <- sub("__.*$", "", score_col)
    score_id <- sub("^(primary|sensitivity)__", "", score_col)
    values <- suppressWarnings(as.numeric(cells[[score_col]]))
    finite <- is.finite(values)
    dataset_scores[[i]] <- data.table(
      dataset = dataset_current,
      method = method,
      score_id = score_id,
      score_column = score_col,
      n_cells_total = length(values),
      n_finite = sum(finite),
      fraction_finite = mean(finite),
      mean = if (any(finite)) mean(values[finite]) else NA_real_,
      sd = if (sum(finite) > 1L) sd(values[finite]) else NA_real_,
      median = if (any(finite)) median(values[finite]) else NA_real_,
      q25 = if (any(finite)) as.numeric(quantile(values[finite], 0.25, names = FALSE)) else NA_real_,
      q75 = if (any(finite)) as.numeric(quantile(values[finite], 0.75, names = FALSE)) else NA_real_,
      minimum = if (any(finite)) min(values[finite]) else NA_real_,
      maximum = if (any(finite)) max(values[finite]) else NA_real_
    )
  }
  score_rows[[dataset_current]] <- rbindlist(dataset_scores)

  rm(cells)
  invisible(gc())
}

atomic_fwrite(rbindlist(audit_rows), file.path(output_dir, "rds_object_audit_v1.tsv"))
atomic_fwrite(rbindlist(overview_rows), file.path(output_dir, "dataset_overview_from_rds_v1.tsv"))
atomic_fwrite(rbindlist(compartment_rows), file.path(output_dir, "analysis_compartment_counts_from_rds_v1.tsv"))
atomic_fwrite(rbindlist(fine_rows), file.path(output_dir, "fine_annotation_counts_from_rds_v1.tsv"))
atomic_fwrite(rbindlist(state_rows), file.path(output_dir, "spatial_state_counts_from_rds_v1.tsv"))
atomic_fwrite(rbindlist(clinical_rows), file.path(output_dir, "clinical_field_counts_from_rds_v1.tsv"))
atomic_fwrite(rbindlist(patient_rows), file.path(output_dir, "patient_coverage_from_rds_v1.tsv"))
atomic_fwrite(rbindlist(score_rows), file.path(output_dir, "program_score_distribution_from_rds_v1.tsv"))

final_status <- fread(file.path(final_run, "FINAL_CNV_GATED_INPUT_STATUS_v4.tsv"))
if (nrow(final_status) != 5L || !all(final_status$status == "PASS")) stop("Final CNV status gate failed")
atomic_fwrite(final_status, file.path(output_dir, "final_cnv_dataset_summary_v1.tsv"))

source_rows <- list()
for (dataset_current in datasets) {
  source_path <- file.path(final_run, dataset_current, "tables", "final_cnv_malignancy_counts_v4.tsv")
  source <- fread(source_path)
  source_rows[[dataset_current]] <- source[, .(
    tumor_like_cells = sum(n_cells),
    final_malignant_cells = sum(n_final_malignant)
  ), by = .(dataset, final_malignancy_source)]
}
atomic_fwrite(rbindlist(source_rows), file.path(output_dir, "final_cnv_source_summary_v1.tsv"))

associations <- fread(file.path(analysis_run, "within_dataset_patient_associations_v6.tsv"))
if (nrow(associations) != 120L) stop("Association grid is not 120 rows")
relation_summary <- associations[, .(
  grid_rows = .N,
  tested_rows = sum(status == "TESTED"),
  expected_direction_rows = sum(direction_matches_expected %in% TRUE, na.rm = TRUE),
  nominal_p_lt_0_05 = sum(p_value < 0.05, na.rm = TRUE),
  fdr_lt_0_05 = sum(p_fdr_bh_within_dataset_method < 0.05, na.rm = TRUE)
), by = relation]
atomic_fwrite(relation_summary, file.path(output_dir, "relation_result_summary_v1.tsv"))

longitudinal <- fread(file.path(analysis_run, "GSE264586_paired_longitudinal_v6.tsv"))
response <- fread(file.path(analysis_run, "GSE294109_response_stratified_v6.tsv"))
sensitivity <- fread(file.path(analysis_run, "GSE178481_full_vs_added_subset_sensitivity_v6.tsv"))
extensions <- data.table(
  analysis = c("GSE264586_paired_longitudinal", "GSE294109_response_stratified", "GSE178481_full_vs_added_subset"),
  rows = c(nrow(longitudinal), nrow(response), nrow(sensitivity)),
  source_file = c(
    "GSE264586_paired_longitudinal_v6.tsv",
    "GSE294109_response_stratified_v6.tsv",
    "GSE178481_full_vs_added_subset_sensitivity_v6.tsv"
  )
)
atomic_fwrite(extensions, file.path(output_dir, "extension_analysis_inventory_v1.tsv"))

software <- data.table(
  item = c("R", "platform", "data.table"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    R.version$platform,
    as.character(packageVersion("data.table"))
  )
)
atomic_fwrite(software, file.path(output_dir, "software_versions_v1.tsv"))

status <- data.table(
  status = "PASS_RDS_EXTRACTION",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = length(datasets),
  rds_rows = sum(rbindlist(overview_rows)$cells_in_compact_rds),
  report_figures_generated = FALSE,
  note = "Five compact score RDS files were read serially; final CNV and inference summaries were read from independently validated v6 tables."
)
atomic_fwrite(status, file.path(output_dir, "PART1_RDS_REPORT_EXTRACTION_STATUS_v1.tsv"))

cat("PASS_RDS_EXTRACTION\n")
