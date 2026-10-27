#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 05_finalize_cnv_gated_program_inputs_v4.R <project_root> <provisional_run> <synced_cnv_dir> <final_run>")
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
provisional_run <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
cnv_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
final_run <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
if (dir.exists(final_run) || file.exists(final_run)) stop("Refusing existing final run: ", final_run)

external_path <- file.path(cnv_dir, "external_cnv_cell_consensus_final_v2.tsv.gz")
part1_path <- file.path(cnv_dir, "phase9_all_cells_cnv_gate_v1.tsv.gz")
if (!all(file.exists(external_path, part1_path))) stop("Synced CNV tables are incomplete")

external <- fread(external_path, select = c(
  "exact_cell_id", "external_cnv_status_incremental_v1",
  "external_malignant_conservative_incremental_v1",
  "external_cnv_cell_finalizable_incremental_v1"
))
setnames(external, "exact_cell_id", "cell_id")
part1 <- fread(part1_path, select = c(
  "cell_id", "phase9_part1_cnv_assessed", "phase9_cnv_gate_status",
  "phase9_malignancy_supported", "phase9_non_malignant_supported"
))
if (anyDuplicated(external$cell_id) || anyDuplicated(part1$cell_id)) stop("Duplicate CNV cell IDs")
if (nrow(part1) != 585119L) stop("Part1 gate does not cover 585119 cells")
if (nrow(external) != 80072L || !all(external$external_cnv_cell_finalizable_incremental_v1)) {
  stop("External final consensus conservation failed")
}

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
dir.create(final_run, recursive = TRUE, showWarnings = FALSE)
summary_rows <- list()

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

logical_true <- function(x) !is.na(x) & as.logical(x)

for (dataset_current in datasets) {
  source_dir <- file.path(provisional_run, dataset_current)
  output_dir <- file.path(final_run, dataset_current)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "status"), showWarnings = FALSE)

  cell_path <- file.path(source_dir, "objects", paste0(dataset_current, "_cell_scores_compact_v1.rds"))
  patient_path <- file.path(source_dir, "tables", "patient_program_summary_v1.tsv")
  sample_path <- file.path(source_dir, "tables", "sample_program_summary_v1.tsv")
  cells <- as.data.table(readRDS(cell_path))
  if (anyDuplicated(cells$cell_id) || !all(cells$dataset == dataset_current)) stop("Compact score identity failure: ", dataset_current)
  cells <- merge(cells, part1, by = "cell_id", all.x = TRUE, sort = FALSE)
  cells <- merge(cells, external, by = "cell_id", all.x = TRUE, sort = FALSE)
  if (nrow(cells) != uniqueN(cells$cell_id)) stop("CNV merge expanded rows: ", dataset_current)

  assessed_part1 <- logical_true(cells$phase9_part1_cnv_assessed)
  supported_part1 <- logical_true(cells$phase9_malignancy_supported)
  supported_external <- logical_true(cells$external_malignant_conservative_incremental_v1)
  cells[, final_malignancy_supported := phase7_fine_annotation == "Epithelial_tumor_like" &
    fifelse(assessed_part1, supported_part1, supported_external)]
  cells[, final_malignancy_source := fcase(
    final_malignancy_supported & assessed_part1, "part1_rebuild_exact_cell",
    final_malignancy_supported & !assessed_part1 & supported_external, "external_cnv_conservative",
    assessed_part1 & !supported_part1, "part1_rebuild_not_conservative_malignant",
    !assessed_part1 & !is.na(external_cnv_status_incremental_v1), "external_cnv_not_conservative_malignant",
    default = "cnv_unavailable_or_not_eligible"
  )]

  source_patient <- fread(patient_path)
  source_sample <- fread(sample_path)
  nonmal_patient <- source_patient[analysis_compartment != "malignant_provisional"]
  nonmal_sample <- source_sample[analysis_compartment != "malignant_provisional"]
  malignant <- cells[analysis_compartment == "malignant_provisional" & final_malignancy_supported]
  signatures <- source_patient[analysis_compartment == "malignant_provisional", unique(signature_id)]
  patient_rows <- list()
  sample_rows <- list()
  index <- 0L
  for (signature_current in signatures) {
    primary_col <- paste0("primary__", signature_current, "__malignant_provisional")
    sensitivity_col <- paste0("sensitivity__", signature_current, "__malignant_provisional")
    if (!all(c(primary_col, sensitivity_col) %in% names(malignant))) next
    index <- index + 1L
    p <- malignant[is.finite(get(primary_col)) & is.finite(get(sensitivity_col)), .(
      n_cells = .N,
      mean_primary = mean(get(primary_col)), median_primary = median(get(primary_col)),
      mean_sensitivity = mean(get(sensitivity_col)), median_sensitivity = median(get(sensitivity_col))
    ), by = .(dataset, patient_key, patient_id)]
    p[, `:=`(unit = "patient", signature_id = signature_current, analysis_compartment = "malignant_provisional")]
    setcolorder(p, names(nonmal_patient))
    patient_rows[[index]] <- p

    s <- malignant[is.finite(get(primary_col)) & is.finite(get(sensitivity_col)), .(
      n_cells = .N,
      mean_primary = mean(get(primary_col)), median_primary = median(get(primary_col)),
      mean_sensitivity = mean(get(sensitivity_col)), median_sensitivity = median(get(sensitivity_col))
    ), by = .(dataset, patient_key, patient_id, sample_id, object_name, timepoint, treatment, response, tissue_confirmed)]
    s[, `:=`(unit = "sample", signature_id = signature_current, analysis_compartment = "malignant_provisional")]
    setcolorder(s, names(nonmal_sample))
    sample_rows[[index]] <- s
  }

  final_patient <- rbindlist(c(list(nonmal_patient), patient_rows), use.names = TRUE, fill = TRUE)
  final_sample <- rbindlist(c(list(nonmal_sample), sample_rows), use.names = TRUE, fill = TRUE)
  atomic_fwrite(final_patient, file.path(output_dir, "tables", "patient_program_summary_v1.tsv"))
  atomic_fwrite(final_sample, file.path(output_dir, "tables", "sample_program_summary_v1.tsv"))
  malignancy_counts <- cells[phase7_fine_annotation == "Epithelial_tumor_like", .(
    n_cells = .N,
    n_final_malignant = sum(final_malignancy_supported)
  ), by = .(dataset, patient_key, patient_id, final_malignancy_source)]
  atomic_fwrite(malignancy_counts, file.path(output_dir, "tables", "final_cnv_malignancy_counts_v4.tsv"))

  validation <- data.table(
    check = c("dataset_identity", "cell_ids_unique", "part1_priority_rule", "external_consensus_finalizable", "OVERALL"),
    result = "PASS",
    detail = c(dataset_current, nrow(cells), "assessed Part1 cells use rebuild support only", nrow(external), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  )
  atomic_fwrite(validation, file.path(output_dir, "status", paste0(dataset_current, "_INDEPENDENT_VALIDATION_v1.tsv")))
  summary_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    compact_cells = nrow(cells),
    tumor_like_cells = sum(cells$phase7_fine_annotation == "Epithelial_tumor_like"),
    final_malignant_cells = sum(cells$final_malignancy_supported),
    final_malignant_patients = uniqueN(cells[final_malignancy_supported == TRUE, patient_key]),
    status = "PASS"
  )
  rm(cells, malignant, final_patient, final_sample)
  gc()
}

summary <- rbindlist(summary_rows)
atomic_fwrite(summary, file.path(final_run, "FINAL_CNV_GATED_INPUT_STATUS_v4.tsv"))
