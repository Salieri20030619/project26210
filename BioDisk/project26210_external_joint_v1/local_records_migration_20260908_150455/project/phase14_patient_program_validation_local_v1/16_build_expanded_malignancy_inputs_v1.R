#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 16_build_expanded_malignancy_inputs_v1.R <project_root> <provisional_run> <synced_cnv_dir> <expanded_run>")
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
provisional_run <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
cnv_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
expanded_run <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
if (dir.exists(expanded_run) || file.exists(expanded_run)) stop("Refusing existing output: ", expanded_run)

external_path <- file.path(cnv_dir, "external_cnv_cell_consensus_final_v2.tsv.gz")
part1_path <- file.path(cnv_dir, "phase9_all_cells_cnv_gate_v1.tsv.gz")
if (!all(file.exists(external_path, part1_path))) stop("CNV input missing")

external <- fread(external_path, select = c(
  "exact_cell_id", "external_cnv_status_incremental_v1",
  "external_malignant_conservative_incremental_v1", "external_cnv_cell_finalizable_incremental_v1",
  "infercnv_call", "copykat_call"
))
setnames(external, "exact_cell_id", "cell_id")
part1 <- fread(part1_path, select = c(
  "cell_id", "phase9_part1_cnv_assessed", "phase9_cnv_gate_status",
  "phase9_malignancy_supported", "phase9_non_malignant_supported",
  "part1__infercnv_call_v1", "part1__copykat_call_v1"
))
if (anyDuplicated(external$cell_id) || anyDuplicated(part1$cell_id)) stop("Duplicate CNV cell ID")
if (nrow(part1) != 585119L || nrow(external) != 80072L) stop("CNV input conservation failure")
if (!all(external$external_cnv_cell_finalizable_incremental_v1)) stop("External CNV not finalizable")

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
dir.create(expanded_run, recursive = TRUE, showWarnings = FALSE)
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
  output_dir <- file.path(expanded_run, dataset_current)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "status"), showWarnings = FALSE)

  cell_path <- file.path(source_dir, "objects", paste0(dataset_current, "_cell_scores_compact_v1.rds"))
  cells <- as.data.table(readRDS(cell_path))
  if (anyDuplicated(cells$cell_id) || !all(cells$dataset == dataset_current)) stop("Cell identity failure: ", dataset_current)
  cells <- merge(cells, part1, by = "cell_id", all.x = TRUE, sort = FALSE)
  cells <- merge(cells, external, by = "cell_id", all.x = TRUE, sort = FALSE)
  if (nrow(cells) != uniqueN(cells$cell_id)) stop("CNV merge expanded rows: ", dataset_current)

  tumor_like <- cells$phase7_fine_annotation == "Epithelial_tumor_like"
  assessed_part1 <- logical_true(cells$phase9_part1_cnv_assessed)
  high_part1 <- logical_true(cells$phase9_malignancy_supported)
  high_external <- logical_true(cells$external_malignant_conservative_incremental_v1)
  part1_any <- cells$part1__infercnv_call_v1 == "high_cnv" | cells$part1__copykat_call_v1 == "aneuploid"
  external_any <- cells$infercnv_call == "high_cnv" | cells$copykat_call == "aneuploid"
  part1_any[is.na(part1_any)] <- FALSE
  external_any[is.na(external_any)] <- FALSE

  cells[, high_confidence_malignant := tumor_like & fifelse(assessed_part1, high_part1, high_external)]
  cells[, expanded_malignancy_supported := tumor_like & fifelse(assessed_part1, part1_any, external_any)]
  cells[, malignancy_evidence_tier := fcase(
    high_confidence_malignant & assessed_part1, "high_confidence_part1_rebuild",
    high_confidence_malignant & !assessed_part1, "high_confidence_external_consensus",
    expanded_malignancy_supported & !high_confidence_malignant & assessed_part1, "supportive_part1_single_method",
    expanded_malignancy_supported & !high_confidence_malignant & !assessed_part1, "supportive_external_single_method",
    tumor_like & assessed_part1, "not_supported_part1_rebuild",
    tumor_like & !assessed_part1 & !is.na(external_cnv_status_incremental_v1), "not_supported_external",
    tumor_like, "cnv_unavailable",
    default = "not_tumor_like"
  )]
  if (any(cells$high_confidence_malignant & !cells$expanded_malignancy_supported)) stop("Tier nesting failure: ", dataset_current)

  source_patient <- fread(file.path(source_dir, "tables", "patient_program_summary_v1.tsv"))
  source_sample <- fread(file.path(source_dir, "tables", "sample_program_summary_v1.tsv"))
  nonmal_patient <- source_patient[analysis_compartment != "malignant_provisional"]
  nonmal_sample <- source_sample[analysis_compartment != "malignant_provisional"]
  malignant <- cells[analysis_compartment == "malignant_provisional" & expanded_malignancy_supported]
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

  expanded_patient <- rbindlist(c(list(nonmal_patient), patient_rows), use.names = TRUE, fill = TRUE)
  expanded_sample <- rbindlist(c(list(nonmal_sample), sample_rows), use.names = TRUE, fill = TRUE)
  atomic_fwrite(expanded_patient, file.path(output_dir, "tables", "patient_program_summary_v1.tsv"))
  atomic_fwrite(expanded_sample, file.path(output_dir, "tables", "sample_program_summary_v1.tsv"))

  tier_counts <- cells[tumor_like, .(n_cells = .N), by = .(
    dataset, patient_key, patient_id, sample_id, object_name, timepoint, malignancy_evidence_tier
  )]
  atomic_fwrite(tier_counts, file.path(output_dir, "tables", "expanded_malignancy_evidence_counts_v1.tsv"))

  selected_cross_section <- malignant[
    !(dataset == "GSE178481" & grepl("Normal", sample_id, ignore.case = TRUE)) &
    !(dataset == "GSE264586" & timepoint != "Baseline")
  ]
  validation <- data.table(
    check = c("dataset_identity", "cell_ids_unique", "part1_rebuild_priority", "tier_nesting", "external_finalizable", "OVERALL"),
    result = "PASS",
    detail = c(
      dataset_current, nrow(cells),
      "Part1-assessed cells use only rebuild method calls; external calls do not override",
      paste0(sum(cells$high_confidence_malignant), " high within ", sum(cells$expanded_malignancy_supported), " expanded"),
      nrow(external), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
  )
  atomic_fwrite(validation, file.path(output_dir, "status", paste0(dataset_current, "_INDEPENDENT_VALIDATION_v1.tsv")))
  summary_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    compact_cells = nrow(cells),
    tumor_like_cells = sum(tumor_like),
    high_confidence_cells = sum(cells$high_confidence_malignant),
    supportive_cells = sum(cells$expanded_malignancy_supported & !cells$high_confidence_malignant),
    expanded_malignant_cells = sum(cells$expanded_malignancy_supported),
    high_confidence_patients_ge20_cross_section = selected_cross_section[high_confidence_malignant == TRUE, .N, by = patient_key][N >= 20, uniqueN(patient_key)],
    expanded_patients_ge20_cross_section = selected_cross_section[, .N, by = patient_key][N >= 20, uniqueN(patient_key)],
    status = "PASS_POSTHOC_EXPANDED_MALIGNANCY"
  )
  rm(cells, malignant, expanded_patient, expanded_sample, selected_cross_section)
  gc()
}

summary <- rbindlist(summary_rows)
if (sum(summary$high_confidence_cells) != 16645L) stop("High-confidence count does not reproduce strict v6")
atomic_fwrite(summary, file.path(expanded_run, "EXPANDED_MALIGNANCY_INPUT_STATUS_v1.tsv"))
