#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath("D:/project26210_external_joint_v1", winslash = "/", mustWork = TRUE)
stage <- file.path(project_root, "phase14_patient_program_validation_local_v1")
provisional <- file.path(stage, "runs", "20260812_082000_provisional_v3")
sync <- file.path(project_root, "phase11_external_cnv_consensus_local_v1", "final_server_sync_v5")
out <- file.path(stage, "posthoc_malignancy_evidence_tier_programs_v5")
if (dir.exists(out) || file.exists(out)) stop("Refusing existing output: ", out)

part1 <- fread(file.path(sync, "phase9_all_cells_cnv_gate_v1.tsv.gz"), select = c(
  "cell_id", "phase7_annotation_confidence", "phase9_part1_cnv_assessed", "phase9_malignancy_supported",
  "part1__infercnv_call_v1", "part1__copykat_call_v1"
))
external <- fread(file.path(sync, "external_cnv_cell_consensus_final_v2.tsv.gz"), select = c(
  "exact_cell_id", "external_malignant_conservative_incremental_v1", "infercnv_call", "copykat_call"
))
setnames(external, "exact_cell_id", "cell_id")
if (anyDuplicated(part1$cell_id) || anyDuplicated(external$cell_id)) stop("Duplicate CNV cell ID")

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
signatures <- c("f1_top100_weighted", "tymp_expression", "hypoxia_angiogenesis", "polyamine_mta_context")
dir.create(out, recursive = TRUE)

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}
logical_true <- function(x) !is.na(x) & as.logical(x)

audit_rows <- list()
patient_count_rows <- list()
program_rows <- list()

for (dataset_current in datasets) {
  cell_path <- file.path(provisional, dataset_current, "objects", paste0(dataset_current, "_cell_scores_compact_v1.rds"))
  cells <- as.data.table(readRDS(cell_path))
  if (anyDuplicated(cells$cell_id) || !all(cells$dataset == dataset_current)) stop("Cell identity failure: ", dataset_current)
  cells <- merge(cells, part1, by = "cell_id", all.x = TRUE, sort = FALSE)
  cells <- merge(cells, external, by = "cell_id", all.x = TRUE, sort = FALSE)
  if (nrow(cells) != uniqueN(cells$cell_id)) stop("CNV merge expansion: ", dataset_current)

  cells <- cells[phase7_fine_annotation == "Epithelial_tumor_like"]
  assessed <- logical_true(cells$phase9_part1_cnv_assessed)
  high_part1 <- logical_true(cells$phase9_malignancy_supported)
  high_external <- logical_true(cells$external_malignant_conservative_incremental_v1)
  p_infer <- cells$part1__infercnv_call_v1 == "high_cnv"; p_infer[is.na(p_infer)] <- FALSE
  p_copy <- cells$part1__copykat_call_v1 == "aneuploid"; p_copy[is.na(p_copy)] <- FALSE
  e_infer <- cells$infercnv_call == "high_cnv"; e_infer[is.na(e_infer)] <- FALSE
  e_copy <- cells$copykat_call == "aneuploid"; e_copy[is.na(e_copy)] <- FALSE

  cells[, evidence_detail := fcase(
    assessed & high_part1, "high_confidence_part1_rebuild",
    !assessed & high_external, "high_confidence_external_consensus",
    assessed & p_infer & !p_copy, "supportive_part1_infercnv_only",
    assessed & !p_infer & p_copy, "supportive_part1_copykat_only",
    !assessed & e_infer & !e_copy, "supportive_external_infercnv_only",
    !assessed & !e_infer & e_copy, "supportive_external_copykat_only",
    assessed, "not_supported_part1_rebuild",
    !assessed & (!is.na(infercnv_call) | !is.na(copykat_call)), "not_supported_external",
    default = "cnv_unavailable"
  )]
  cells[, evidence_tier := fcase(
    grepl("^high_confidence", evidence_detail), "high_confidence",
    grepl("^supportive", evidence_detail), "supportive",
    default = "not_in_expanded"
  )]
  cells[, expanded := evidence_tier %in% c("high_confidence", "supportive")]
  if (sum(cells$evidence_tier == "high_confidence") != sum(fifelse(assessed, high_part1, high_external))) {
    stop("High-confidence reproduction failure: ", dataset_current)
  }

  audit_rows[[dataset_current]] <- cells[, .(
    n_cells = .N,
    n_patients = uniqueN(patient_key),
    n_samples = uniqueN(object_name)
  ), by = .(dataset, evidence_tier, evidence_detail, phase7_annotation_confidence)]

  selected <- cells[
    !(dataset == "GSE178481" & grepl("Normal", sample_id, ignore.case = TRUE)) &
    !(dataset == "GSE264586" & timepoint != "Baseline")
  ]
  counts <- selected[evidence_tier %in% c("high_confidence", "supportive"), .N,
    by = .(dataset, patient_key, patient_id, evidence_tier)]
  expanded_counts <- counts[, .(n_cells = sum(N)), by = .(dataset, patient_key, patient_id)]
  expanded_counts[, evidence_tier := "expanded"]
  setnames(counts, "N", "n_cells")
  patient_count_rows[[dataset_current]] <- rbindlist(list(counts, expanded_counts), use.names = TRUE)

  program_groups <- c(
    "high_confidence", "supportive", "expanded",
    "expanded_exclude_low_annotation", "high_plus_infercnv_only_supportive",
    "supportive_part1_infercnv_only", "supportive_part1_copykat_only",
    "supportive_external_infercnv_only", "supportive_external_copykat_only"
  )
  for (tier_current in program_groups) {
    tier_cells <- if (tier_current == "expanded") {
      selected[expanded == TRUE]
    } else if (tier_current == "expanded_exclude_low_annotation") {
      selected[expanded == TRUE & phase7_annotation_confidence %in% c("high", "moderate")]
    } else if (tier_current == "high_plus_infercnv_only_supportive") {
      selected[evidence_tier == "high_confidence" | grepl("infercnv_only$", evidence_detail)]
    } else if (tier_current %in% c("high_confidence", "supportive")) {
      selected[evidence_tier == tier_current]
    } else {
      selected[evidence_detail == tier_current]
    }
    for (signature_current in signatures) {
      primary_col <- paste0("primary__", signature_current, "__malignant_provisional")
      sensitivity_col <- paste0("sensitivity__", signature_current, "__malignant_provisional")
      if (!all(c(primary_col, sensitivity_col) %in% names(tier_cells))) next
      value <- tier_cells[is.finite(get(primary_col)) & is.finite(get(sensitivity_col)), .(
        n_cells = .N,
        mean_primary = mean(get(primary_col)),
        mean_sensitivity = mean(get(sensitivity_col))
      ), by = .(dataset, patient_key, patient_id)]
      value[, `:=`(evidence_tier = tier_current, signature_id = signature_current)]
      program_rows[[length(program_rows) + 1L]] <- value
    }
  }
  rm(cells, selected)
  gc()
}

audit <- rbindlist(audit_rows, use.names = TRUE, fill = TRUE)
patient_counts <- rbindlist(patient_count_rows, use.names = TRUE, fill = TRUE)
programs <- rbindlist(program_rows, use.names = TRUE, fill = TRUE)
atomic_fwrite(audit, file.path(out, "cell_evidence_method_and_annotation_audit_v1.tsv"))
atomic_fwrite(patient_counts, file.path(out, "cross_sectional_patient_evidence_tier_counts_v1.tsv"))
atomic_fwrite(programs, file.path(out, "cross_sectional_patient_evidence_tier_programs_v1.tsv"))

status <- data.table(
  status = "PASS_POSTHOC_TIER_PROGRAMS",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  tumor_like_cells = sum(audit$n_cells),
  high_confidence_cells = sum(audit[evidence_tier == "high_confidence", n_cells]),
  supportive_cells = sum(audit[evidence_tier == "supportive", n_cells]),
  tier_program_rows = nrow(programs),
  rebuild_priority = TRUE
)
if (status$high_confidence_cells != 16645L || status$supportive_cells != 17775L) stop("Tier count failure")
atomic_fwrite(status, file.path(out, "MALIGNANCY_EVIDENCE_TIER_PROGRAM_STATUS_v1.tsv"))
