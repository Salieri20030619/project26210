#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

PROJECT <- "/mnt/BioDisk/ludingyu/project26210"
ROOT <- file.path(PROJECT, "results/part1_rebuild_v1")
INPUT <- file.path(ROOT, "05_cnmf/input_v1")
FULL <- file.path(ROOT, "05_cnmf/full_ranks4_8_v1")
EVIDENCE <- file.path(ROOT, "05_cnmf/tymp_rank_evidence_v1")
SENS <- file.path(ROOT, "05_cnmf/tymp_study_sensitivity_v1")
OUT <- file.path(ROOT, "05_cnmf/tymp_study_sensitivity_verification_v1")
MEMORY_GATE <- file.path(PROJECT, "check_memory_gate.sh")
RANKS <- 4:8

if (dir.exists(OUT) && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty verification directory: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(name, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, status = if (isTRUE(passed)) "PASS" else "FAIL",
    evidence = paste(evidence, collapse = ";")
  )
  fwrite(rbindlist(checks, fill = TRUE), file.path(OUT, "verification_checks_v1.tsv"), sep = "\t")
  if (!isTRUE(passed)) stop("Verification failed: ", name, ": ", paste(evidence, collapse = ";"))
}

gate_output <- system2(MEMORY_GATE, "cnmf_tymp_study_sensitivity_independent_verification_v1", stdout = TRUE, stderr = TRUE)
gate_status <- attr(gate_output, "status")
if (is.null(gate_status)) gate_status <- 0L
gate_text <- paste(gate_output, collapse = " ")
message(gate_text)
check("memory_gate", gate_status == 0L && grepl("PASS", gate_text, fixed = TRUE), gate_text)
check("sensitivity_complete", file.exists(file.path(SENS, "RUN_COMPLETE")), "terminal")
check("evidence_verification_complete", file.exists(file.path(ROOT, "05_cnmf/tymp_rank_evidence_verification_v1/RUN_COMPLETE")), "terminal")

metadata <- fromJSON(file.path(SENS, "run_metadata_v1.json"), simplifyVector = TRUE)
check("metadata_schema", identical(metadata$schema_version, "part1_cnmf_tymp_study_sensitivity_v1"), metadata$schema_version)
check("metadata_selection", metadata$selected_rank == 7L && identical(metadata$selected_TYMP_factor, "F1") && identical(metadata$selection_status, "selected_with_study_composition_warning"), c(metadata$selected_rank, metadata$selected_TYMP_factor, metadata$selection_status))
check("metadata_study_count", metadata$n_studies == 5L, metadata$n_studies)
check("metadata_patient_unit", identical(metadata$patient_statistical_unit, TRUE) && identical(metadata$cell_level_pseudoreplication_used, FALSE), c(metadata$patient_statistical_unit, metadata$cell_level_pseudoreplication_used))
check("metadata_warning", identical(metadata$study_composition_warning, TRUE), metadata$study_composition_warning)
check("metadata_no_dependency_changes", identical(metadata$dependency_changes, FALSE), metadata$dependency_changes)

gate_log <- fread(file.path(SENS, "memory_gate_log_v1.tsv"))
check("analysis_memory_gate", nrow(gate_log) == 1L && gate_log$exit_status == 0L && grepl("PASS", gate_log$gate_output, fixed = TRUE), gate_log$gate_output)

meta <- fread(file.path(INPUT, "cnmf_cell_metadata_v1.tsv.gz"),
              colClasses = c(exact_cell_id = "character", curated_patient_id = "character"))
candidate <- fread(file.path(EVIDENCE, "TYMP_candidate_by_rank_v1.tsv"))
patient_cor <- fread(file.path(EVIDENCE, "factor_patient_stratified_correlations_v1.tsv.gz"),
                     colClasses = c(curated_patient_id = "character"))
patient_study <- unique(meta[, .(curated_patient_id, study)])
studies <- sort(unique(meta$study))
candidate_patient <- merge(
  patient_cor[outcome == "TYMP_log1p_CP10K_v1"], candidate[, .(rank, factor)],
  by = c("rank", "factor"), all = FALSE
)
candidate_patient <- merge(candidate_patient, patient_study, by = "curated_patient_id", all.x = TRUE)
check("candidate_patient_rows", nrow(candidate_patient) == 100L && !anyNA(candidate_patient$study), nrow(candidate_patient))

study_summary <- fread(file.path(SENS, "candidate_study_stratified_patient_correlations_v1.tsv"))
check("study_summary_rows", nrow(study_summary) == 25L && !anyDuplicated(study_summary[, .(rank, factor, study)]), nrow(study_summary))
study_check <- candidate_patient[, .(
  n_patients_check = .N,
  median_check = median(rho[is.finite(rho)]),
  positive_check = mean(rho[is.finite(rho)] > 0)
), by = .(rank, factor, study)]
study_compare <- merge(study_summary, study_check, by = c("rank", "factor", "study"), all = TRUE)
check("study_summary_patient_counts", all(study_compare$n_patients == study_compare$n_patients_check), range(study_compare$n_patients - study_compare$n_patients_check))
check("study_summary_medians", max(abs(study_compare$median_patient_rho - study_compare$median_check)) < 1e-12, max(abs(study_compare$median_patient_rho - study_compare$median_check)))
check("study_summary_positive_fraction", max(abs(study_compare$fraction_patients_positive - study_compare$positive_check)) < 1e-12, max(abs(study_compare$fraction_patients_positive - study_compare$positive_check)))

loo <- fread(file.path(SENS, "candidate_leave_one_study_out_v1.tsv"))
check("loo_rows", nrow(loo) == 25L && !anyDuplicated(loo[, .(rank, factor, excluded_study)]), nrow(loo))
loo_check <- rbindlist(lapply(seq_len(nrow(candidate)), function(index) {
  source <- candidate_patient[rank == candidate$rank[index] & factor == candidate$factor[index]]
  rbindlist(lapply(studies, function(excluded_study) {
    values <- source[study != excluded_study & is.finite(rho), rho]
    data.table(
      rank = candidate$rank[index], factor = candidate$factor[index], excluded_study = excluded_study,
      n_check = length(values), median_check = median(values), positive_check = mean(values > 0)
    )
  }))
}))
loo_compare <- merge(loo, loo_check, by = c("rank", "factor", "excluded_study"), all = TRUE)
check("loo_patient_counts", all(loo_compare$n_patients_retained == loo_compare$n_check), range(loo_compare$n_patients_retained - loo_compare$n_check))
check("loo_medians", max(abs(loo_compare$median_patient_rho - loo_compare$median_check)) < 1e-12, max(abs(loo_compare$median_patient_rho - loo_compare$median_check)))
check("loo_positive_fraction", max(abs(loo_compare$fraction_patients_positive - loo_compare$positive_check)) < 1e-12, max(abs(loo_compare$fraction_patients_positive - loo_compare$positive_check)))
rank7_loo <- loo[rank == 7L & factor == "F1"]
check("rank7_loo_minimums", abs(min(rank7_loo$median_patient_rho) - 0.1018133) < 1e-6 && abs(min(rank7_loo$fraction_patients_positive) - 0.7) < 1e-12 && min(rank7_loo$n_patients_retained) == 10L, c(min(rank7_loo$median_patient_rho), min(rank7_loo$fraction_patients_positive), min(rank7_loo$n_patients_retained)))

composition <- fread(file.path(SENS, "candidate_dominant_study_composition_v1.tsv"))
check("composition_rows", nrow(composition) == 25L && !anyDuplicated(composition[, .(rank, factor, study)]), nrow(composition))
check("composition_baseline_counts", all(composition[, .(total = sum(n_all_cells)), by = .(rank, factor)]$total == 11950L), composition[, .(total = sum(n_all_cells)), by = .(rank, factor)]$total)
check("composition_dominant_fractions", all(abs(composition[, .(total = sum(fraction_of_candidate_dominant_cells)), by = .(rank, factor)]$total - 1) < 1e-12), composition[, .(total = sum(fraction_of_candidate_dominant_cells)), by = .(rank, factor)]$total)

for (index in seq_len(nrow(candidate))) {
  rank_value <- candidate$rank[index]
  factor_value <- candidate$factor[index]
  usage <- fread(file.path(FULL, paste0("consensus_usage_rank", rank_value, "_v1.tsv.gz")),
                 colClasses = c(exact_cell_id = "character"))
  factor_names <- paste0("F", seq_len(rank_value))
  dominant <- factor_names[max.col(as.matrix(usage[, ..factor_names]), ties.method = "first")]
  counts <- data.table(study = meta$study, dominant = dominant == factor_value)[, .(n_check = sum(dominant)), by = study]
  compare <- merge(composition[rank == rank_value & factor == factor_value], counts, by = "study", all = TRUE)
  check(paste0("rank", rank_value, "_dominant_counts"), all(compare$n_candidate_dominant_cells == compare$n_check), range(compare$n_candidate_dominant_cells - compare$n_check))
}

selection <- fread(file.path(SENS, "rank_selection_study_sensitivity_v1.tsv"))
selected <- fread(file.path(SENS, "selected_rank_v1.tsv"))
check("selection_rows", nrow(selection) == 5L && identical(selection$rank, RANKS), selection$rank)
loo_gate <- with(selection,
  min_leave_one_study_out_n_patients >= 10L &
  min_leave_one_study_out_median_rho >= 0.05 &
  min_leave_one_study_out_positive_fraction >= 0.60)
study_gate <- with(selection,
  n_studies_min3_patients >= 2L &
  n_studies_min3_patients_positive_median == n_studies_min3_patients)
final_gate <- selection$eligible_for_biological_rank_selection & loo_gate & study_gate
check("selection_loo_formula", identical(selection$leave_one_study_out_gate, loo_gate), selection$leave_one_study_out_gate)
check("selection_study_formula", identical(selection$multi_study_direction_gate, study_gate), selection$multi_study_direction_gate)
check("selection_final_formula", identical(selection$final_rank_eligible, final_gate), selection$final_rank_eligible)
check("unique_selected_rank", nrow(selected) == 1L && selected$selected_rank == 7L && selected$selected_TYMP_factor == "F1", c(selected$selected_rank, selected$selected_TYMP_factor))
check("selection_status_warning", selected$selection_status == "selected_with_study_composition_warning" && selected$study_composition_warning == TRUE && selected$candidate_dominant_study == "GSE207493", c(selected$selection_status, selected$candidate_dominant_study))
check("selection_interpretation_boundary", grepl("not genetic subclones", selected$interpretation_limit, fixed = TRUE), selected$interpretation_limit)
decision <- readLines(file.path(SENS, "selection_decision_v1.md"), warn = FALSE)
check("decision_record", any(grepl("Selected rank: 7", decision, fixed = TRUE)) && any(grepl("Warning", decision, fixed = TRUE)), decision)

manifest <- fread(file.path(SENS, "file_manifest_v1.tsv"), colClasses = "character")
check("manifest_unique", !anyDuplicated(manifest$file), uniqueN(manifest$file))
for (i in seq_len(nrow(manifest))) {
  path <- file.path(SENS, manifest$file[i])
  valid <- file.exists(path) && file.info(path)$size == as.numeric(manifest$size_bytes[i]) &&
    identical(digest(path, algo = "sha256", file = TRUE), manifest$sha256[i])
  check(paste0("hash_", manifest$file[i]), valid, manifest$sha256[i])
}

final_checks <- rbindlist(checks, fill = TRUE)
run_metadata <- list(
  schema_version = "part1_cnmf_tymp_study_sensitivity_verification_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_checks = nrow(final_checks), n_pass = sum(final_checks$status == "PASS"),
  n_fail = sum(final_checks$status == "FAIL"),
  selected_rank = 7L, selected_TYMP_factor = "F1",
  selection_status = "selected_with_study_composition_warning",
  selection_bundle = "PASS",
  dependency_changes = FALSE
)
write_json(run_metadata, file.path(OUT, "run_metadata_v1.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(OUT, "RUN_COMPLETE"))
print(run_metadata)
cat("PART1_CNMF_TYMP_STUDY_SENSITIVITY_VERIFICATION_V1_PASS\n")
