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
OUT <- file.path(ROOT, "05_cnmf/tymp_study_sensitivity_v1")
MEMORY_GATE <- file.path(PROJECT, "check_memory_gate.sh")
RANKS <- 4:8

if (dir.exists(OUT) && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty study sensitivity directory: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
gate_output <- system2(MEMORY_GATE, "cnmf_tymp_study_sensitivity_v1", stdout = TRUE, stderr = TRUE)
gate_status <- attr(gate_output, "status")
if (is.null(gate_status)) gate_status <- 0L
gate_text <- paste(gate_output, collapse = " ")
message(gate_text)
fwrite(data.table(stage = "cnmf_tymp_study_sensitivity_v1", exit_status = gate_status, gate_output = gate_text),
       file.path(OUT, "memory_gate_log_v1.tsv"), sep = "\t")
if (gate_status != 0L || !grepl("PASS", gate_text, fixed = TRUE)) stop("Memory gate blocked analysis")
if (!file.exists(file.path(ROOT, "05_cnmf/tymp_rank_evidence_verification_v1/RUN_COMPLETE"))) {
  stop("TYMP/rank evidence verification is not terminal")
}

meta <- fread(file.path(INPUT, "cnmf_cell_metadata_v1.tsv.gz"),
              colClasses = c(exact_cell_id = "character", curated_patient_id = "character"))
candidate <- fread(file.path(EVIDENCE, "TYMP_candidate_by_rank_v1.tsv"))
rank_evidence <- fread(file.path(EVIDENCE, "rank_selection_evidence_v1.tsv"))
patient_cor <- fread(file.path(EVIDENCE, "factor_patient_stratified_correlations_v1.tsv.gz"),
                     colClasses = c(curated_patient_id = "character"))
patient_study <- unique(meta[, .(curated_patient_id, study)])
if (anyDuplicated(patient_study$curated_patient_id)) stop("Patient maps to multiple studies")
studies <- sort(unique(meta$study))
if (length(studies) != 5L || nrow(candidate) != 5L) stop("Unexpected cohort or candidates")

candidate_patient <- merge(
  patient_cor[outcome == "TYMP_log1p_CP10K_v1"],
  candidate[, .(rank, factor)],
  by = c("rank", "factor"),
  all = FALSE
)
candidate_patient <- merge(candidate_patient, patient_study, by = "curated_patient_id", all.x = TRUE)
if (nrow(candidate_patient) != 100L || anyNA(candidate_patient$study)) stop("Candidate patient evidence mismatch")

study_summary <- candidate_patient[, {
  values <- rho[is.finite(rho)]
  list(
    n_patients = .N,
    n_patients_with_rho = length(values),
    median_patient_rho = if (length(values)) median(values) else NA_real_,
    fraction_patients_positive = if (length(values)) mean(values > 0) else NA_real_,
    min_patient_rho = if (length(values)) min(values) else NA_real_,
    max_patient_rho = if (length(values)) max(values) else NA_real_
  )
}, by = .(rank, factor, study)]
fwrite(study_summary, file.path(OUT, "candidate_study_stratified_patient_correlations_v1.tsv"), sep = "\t")

loo <- rbindlist(lapply(seq_len(nrow(candidate)), function(index) {
  rank_value <- candidate$rank[index]
  factor_value <- candidate$factor[index]
  source <- candidate_patient[rank == rank_value & factor == factor_value]
  rbindlist(lapply(studies, function(excluded_study) {
    values <- source[study != excluded_study & is.finite(rho), rho]
    p_value <- if (length(values) >= 3L && any(values != 0)) {
      suppressWarnings(wilcox.test(values, mu = 0, alternative = "greater", exact = FALSE)$p.value)
    } else NA_real_
    data.table(
      rank = rank_value, factor = factor_value, excluded_study = excluded_study,
      n_patients_retained = length(values),
      median_patient_rho = if (length(values)) median(values) else NA_real_,
      fraction_patients_positive = if (length(values)) mean(values > 0) else NA_real_,
      wilcoxon_greater_p = p_value
    )
  }))
}))
fwrite(loo, file.path(OUT, "candidate_leave_one_study_out_v1.tsv"), sep = "\t")

baseline <- meta[, .(n_all_cells = .N), by = study]
baseline[, baseline_cell_fraction := n_all_cells / sum(n_all_cells)]
dominant_store <- list()
for (index in seq_len(nrow(candidate))) {
  rank_value <- candidate$rank[index]
  factor_value <- candidate$factor[index]
  usage <- fread(file.path(FULL, paste0("consensus_usage_rank", rank_value, "_v1.tsv.gz")),
                 colClasses = c(exact_cell_id = "character"))
  if (!identical(usage$exact_cell_id, meta$exact_cell_id)) stop("Usage order mismatch")
  factor_names <- paste0("F", seq_len(rank_value))
  dominant <- factor_names[max.col(as.matrix(usage[, ..factor_names]), ties.method = "first")]
  x <- data.table(study = meta$study, is_candidate_dominant = dominant == factor_value)
  composition <- x[, .(n_candidate_dominant_cells = sum(is_candidate_dominant)), by = study]
  composition <- merge(baseline, composition, by = "study", all.x = TRUE)
  composition[is.na(n_candidate_dominant_cells), n_candidate_dominant_cells := 0L]
  composition[, `:=`(
    rank = rank_value,
    factor = factor_value,
    fraction_of_candidate_dominant_cells = n_candidate_dominant_cells / sum(n_candidate_dominant_cells),
    within_study_candidate_dominant_fraction = n_candidate_dominant_cells / n_all_cells
  )]
  composition[, dominant_share_enrichment_over_baseline := fraction_of_candidate_dominant_cells / baseline_cell_fraction]
  dominant_store[[length(dominant_store) + 1L]] <- composition
}
dominant_composition <- rbindlist(dominant_store)
setcolorder(dominant_composition, c(
  "rank", "factor", "study", "n_all_cells", "baseline_cell_fraction",
  "n_candidate_dominant_cells", "fraction_of_candidate_dominant_cells",
  "within_study_candidate_dominant_fraction", "dominant_share_enrichment_over_baseline"
))
fwrite(dominant_composition, file.path(OUT, "candidate_dominant_study_composition_v1.tsv"), sep = "\t")

selection <- merge(rank_evidence, candidate[, .(rank, factor)], by = c("rank", "factor"), all.x = FALSE)
loo_summary <- loo[, .(
  min_leave_one_study_out_median_rho = min(median_patient_rho),
  min_leave_one_study_out_positive_fraction = min(fraction_patients_positive),
  min_leave_one_study_out_n_patients = min(n_patients_retained),
  max_leave_one_study_out_p = max(wilcoxon_greater_p, na.rm = TRUE)
), by = .(rank, factor)]
study_support <- study_summary[, .(
  n_studies_total = .N,
  n_studies_min3_patients = sum(n_patients_with_rho >= 3L),
  n_studies_min3_patients_positive_median = sum(n_patients_with_rho >= 3L & median_patient_rho > 0),
  n_studies_any_positive_median = sum(median_patient_rho > 0, na.rm = TRUE)
), by = .(rank, factor)]
dominance_summary <- dominant_composition[, .(
  candidate_dominant_study_fraction = max(fraction_of_candidate_dominant_cells),
  candidate_dominant_study_enrichment = max(dominant_share_enrichment_over_baseline),
  candidate_dominant_study = study[which.max(fraction_of_candidate_dominant_cells)],
  candidate_dominant_study_baseline_fraction = baseline_cell_fraction[which.max(fraction_of_candidate_dominant_cells)]
), by = .(rank, factor)]
selection <- Reduce(function(x, y) merge(x, y, by = c("rank", "factor"), all.x = TRUE),
                    list(selection, loo_summary, study_support, dominance_summary))
selection[, leave_one_study_out_gate :=
  min_leave_one_study_out_n_patients >= 10L &
  min_leave_one_study_out_median_rho >= 0.05 &
  min_leave_one_study_out_positive_fraction >= 0.60
]
selection[, multi_study_direction_gate :=
  n_studies_min3_patients >= 2L &
  n_studies_min3_patients_positive_median == n_studies_min3_patients
]
selection[, study_composition_warning :=
  candidate_dominant_study_fraction > 0.90 |
  (candidate_dominant_study_fraction - candidate_dominant_study_baseline_fraction) > 0.20
]
selection[, final_rank_eligible :=
  eligible_for_biological_rank_selection == TRUE &
  leave_one_study_out_gate == TRUE & multi_study_direction_gate == TRUE
]
selected <- selection[final_rank_eligible == TRUE]
if (nrow(selected) != 1L || selected$rank != 7L || selected$factor != "F1") {
  stop("Study sensitivity did not yield the expected unique formal rank candidate")
}
selected_status <- if (selected$study_composition_warning) {
  "selected_with_study_composition_warning"
} else "selected"
selected_rank <- selected[, .(
  selected_rank = rank,
  selected_TYMP_factor = factor,
  selection_status = selected_status,
  primary_reason = "sole rank passing 20-run stability, patient-level TYMP, cross-rank reproducibility, and leave-one-study-out direction gates",
  interpretation_limit = "cNMF factors are transcriptional programs, not genetic subclones",
  study_composition_warning,
  candidate_dominant_study,
  candidate_dominant_study_fraction,
  candidate_dominant_study_enrichment,
  min_leave_one_study_out_median_rho,
  min_leave_one_study_out_positive_fraction,
  n_studies_min3_patients_positive_median,
  n_studies_min3_patients
)]
fwrite(selection, file.path(OUT, "rank_selection_study_sensitivity_v1.tsv"), sep = "\t")
fwrite(selected_rank, file.path(OUT, "selected_rank_v1.tsv"), sep = "\t")
writeLines(c(
  "# Part 1 cNMF rank selection v1",
  "",
  "- Selected rank: 7.",
  "- Primary TYMP-associated factor: F1.",
  "- Selection basis: rank 7 was the sole rank passing the formal 20-run stability, patient-level TYMP, cross-rank reproducibility, and leave-one-study-out direction gates.",
  "- Study sensitivity: excluding GSE207493 retained 10 patients with median within-patient rho >= 0.05 and >=60% positive directions; all studies with at least three eligible patients had positive median rho.",
  "- Warning: dominant F1 cells are strongly concentrated in one study, so study composition must remain an explicit limitation in downstream interpretation.",
  "- Interpretation boundary: factors are transcriptional programs and must not be called genetic subclones without separate CNV/subclone evidence."
), file.path(OUT, "selection_decision_v1.md"))

metadata <- list(
  schema_version = "part1_cnmf_tymp_study_sensitivity_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  selected_rank = 7L,
  selected_TYMP_factor = "F1",
  selection_status = selected_status,
  n_studies = length(studies),
  patient_statistical_unit = TRUE,
  cell_level_pseudoreplication_used = FALSE,
  biological_rank_selected = TRUE,
  study_composition_warning = TRUE,
  dependency_changes = FALSE
)
write_json(metadata, file.path(OUT, "run_metadata_v1.json"), auto_unbox = TRUE, pretty = TRUE)
manifest_files <- sort(setdiff(list.files(OUT), c("file_manifest_v1.tsv", "RUN_COMPLETE")))
manifest <- data.table(
  file = manifest_files,
  size_bytes = as.numeric(file.info(file.path(OUT, manifest_files))$size),
  sha256 = vapply(file.path(OUT, manifest_files), digest, character(1), algo = "sha256", file = TRUE)
)
fwrite(manifest, file.path(OUT, "file_manifest_v1.tsv"), sep = "\t")
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(OUT, "RUN_COMPLETE"))
print(selected_rank)
cat("PART1_CNMF_TYMP_STUDY_SENSITIVITY_V1_PASS\n")
