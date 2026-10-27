#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_run_phase60_f1_tymp_v1.R <project_root> <output_dir>")
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (dir.exists(output_dir) || file.exists(output_dir)) stop("Refusing existing output: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

phase14 <- file.path(project_root, "phase14_patient_program_validation_local_v1")
score_run <- file.path(phase14, "runs", "20260812_082000_provisional_v3")
cnv_dir <- file.path(project_root, "phase11_external_cnv_consensus_local_v1", "final_server_sync_v5")
crosswalk_path <- file.path(
  project_root, "phase6_metadata_confirmation_local_v1", "tables",
  "confirmed_object_patient_crosswalk.tsv"
)
datasets <- c("GSE264586", "GSE294109")

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

logical_true <- function(x) !is.na(x) & as.logical(x)

cor_values <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 3L || uniqueN(x) < 2L || uniqueN(y) < 2L) {
    return(list(n_cells = n, spearman_rho = NA_real_, p_value = NA_real_, status = "UNTESTABLE"))
  }
  result <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(
    n_cells = n,
    spearman_rho = unname(result$estimate),
    p_value = result$p.value,
    status = if (n >= 20L) "ELIGIBLE_GE20" else "DESCRIPTIVE_LT20"
  )
}

external <- fread(
  file.path(cnv_dir, "external_cnv_cell_consensus_final_v2.tsv.gz"),
  select = c(
    "exact_cell_id", "external_cnv_status_incremental_v1",
    "external_malignant_conservative_incremental_v1",
    "external_cnv_cell_finalizable_incremental_v1", "infercnv_call", "copykat_call"
  )
)
setnames(external, "exact_cell_id", "cell_id")
part1 <- fread(
  file.path(cnv_dir, "phase9_all_cells_cnv_gate_v1.tsv.gz"),
  select = c(
    "cell_id", "phase9_part1_cnv_assessed", "phase9_cnv_gate_status",
    "phase9_malignancy_supported", "phase9_non_malignant_supported",
    "part1__infercnv_call_v1", "part1__copykat_call_v1"
  )
)
if (anyDuplicated(external$cell_id) || anyDuplicated(part1$cell_id)) stop("Duplicate CNV cell IDs")
if (nrow(external) != 80072L || nrow(part1) != 585119L) stop("CNV input conservation failure")
if (!all(external$external_cnv_cell_finalizable_incremental_v1)) stop("External CNV is not finalizable")

crosswalk <- fread(crosswalk_path)
crosswalk <- crosswalk[dataset %chin% datasets]
crosswalk_audit <- crosswalk[, .(
  objects = .N,
  patients = uniqueN(patient_id_confirmed),
  samples = uniqueN(sample_id_confirmed),
  timepoints = uniqueN(timepoint_confirmed),
  responses = paste(sort(unique(response_confirmed)), collapse = ";"),
  mapping_pass = all(mapping_status == "PASS")
), by = dataset]
atomic_fwrite(crosswalk_audit, file.path(output_dir, "clinical_mapping_audit_v1.tsv"))

method_definitions <- data.table(
  method = c("primary_gene_z", "sensitivity_log_normalized"),
  f1_column = c(
    "primary__f1_top100_weighted__malignant_provisional",
    "sensitivity__f1_top100_weighted__malignant_provisional"
  ),
  tymp_column = c(
    "primary__tymp_expression__malignant_provisional",
    "sensitivity__tymp_expression__malignant_provisional"
  )
)
atomic_fwrite(method_definitions, file.path(output_dir, "score_method_definitions_v1.tsv"))

all_cells <- list()
correlation_rows <- list()
tier_count_rows <- list()
composition_rows <- list()
immune_state_rows <- list()
row_index <- 0L

for (dataset_current in datasets) {
  rds_path <- file.path(
    score_run, dataset_current, "objects",
    paste0(dataset_current, "_cell_scores_compact_v1.rds")
  )
  cells <- as.data.table(readRDS(rds_path))
  if (anyDuplicated(cells$cell_id) || !all(cells$dataset == dataset_current)) {
    stop("Compact RDS identity failure: ", dataset_current)
  }
  required_scores <- unique(c(method_definitions$f1_column, method_definitions$tymp_column))
  if (length(setdiff(required_scores, names(cells)))) stop("Missing F1/TYMP score columns: ", dataset_current)
  cells <- merge(cells, part1, by = "cell_id", all.x = TRUE, sort = FALSE)
  cells <- merge(cells, external, by = "cell_id", all.x = TRUE, sort = FALSE)
  if (nrow(cells) != uniqueN(cells$cell_id)) stop("CNV merge expansion: ", dataset_current)

  tumor_like <- cells$phase7_fine_annotation == "Epithelial_tumor_like"
  assessed_part1 <- logical_true(cells$phase9_part1_cnv_assessed)
  strict_part1 <- logical_true(cells$phase9_malignancy_supported)
  strict_external <- logical_true(cells$external_malignant_conservative_incremental_v1)
  part1_any <- cells$part1__infercnv_call_v1 == "high_cnv" | cells$part1__copykat_call_v1 == "aneuploid"
  external_any <- cells$infercnv_call == "high_cnv" | cells$copykat_call == "aneuploid"
  part1_any[is.na(part1_any)] <- FALSE
  external_any[is.na(external_any)] <- FALSE
  cells[, strict_malignant := tumor_like & fifelse(assessed_part1, strict_part1, strict_external)]
  cells[, expanded_malignant := tumor_like & fifelse(assessed_part1, part1_any, external_any)]
  if (any(cells$strict_malignant & !cells$expanded_malignant)) stop("Tier nesting failure: ", dataset_current)

  tier_count_rows[[dataset_current]] <- rbindlist(list(
    cells[strict_malignant == TRUE, .(n_cells = .N), by = .(
      dataset, patient_key, patient_id, sample_id, object_name, timepoint,
      treatment, response, tissue_confirmed
    )][, malignancy_tier := "strict_high_confidence"],
    cells[expanded_malignant == TRUE, .(n_cells = .N), by = .(
      dataset, patient_key, patient_id, sample_id, object_name, timepoint,
      treatment, response, tissue_confirmed
    )][, malignancy_tier := "expanded_any_method_positive"]
  ), use.names = TRUE)

  composition_sample <- cells[, .(n_cells = .N), by = .(
    dataset, patient_key, patient_id, sample_id, object_name, timepoint,
    treatment, response, broad_annotation = phase7_broad_annotation,
    fine_annotation = phase7_fine_annotation
  )]
  composition_sample[, fraction := n_cells / sum(n_cells), by = .(dataset, sample_id)]
  composition_sample[, unit := "sample_region"]
  composition_patient_time <- composition_sample[, .(n_cells = sum(n_cells)), by = .(
    dataset, patient_key, patient_id, timepoint, treatment, response,
    broad_annotation, fine_annotation
  )]
  composition_patient_time[, fraction := n_cells / sum(n_cells), by = .(dataset, patient_key, timepoint)]
  composition_patient_time[, `:=`(sample_id = "<pooled>", object_name = "<pooled>", unit = "patient_timepoint")]
  composition_rows[[dataset_current]] <- rbindlist(
    list(composition_sample, composition_patient_time), use.names = TRUE, fill = TRUE
  )

  immune_specs <- list(
    c("tcell_dysfunction", "t_nk"), c("cytotoxicity", "t_nk"),
    c("ifn_response", "t_nk"), c("spp1_state_no_spp1", "myeloid"),
    c("c1qc_state", "myeloid"), c("apoe_tam_expression_proxy", "myeloid"),
    c("monocyte_program", "myeloid"), c("myeloid_immunosuppression", "myeloid"),
    c("immunosuppression_core", "myeloid"), c("ifn_response", "myeloid"),
    c("hla_ii_supporting", "myeloid")
  )
  state_index <- 0L
  state_list <- list()
  for (spec in immune_specs) {
    score_id <- spec[[1L]]
    compartment <- spec[[2L]]
    score_column <- paste0("primary__", score_id, "__", compartment)
    if (!score_column %in% names(cells)) next
    state_index <- state_index + 1L
    state_list[[state_index]] <- cells[
      analysis_compartment == compartment & is.finite(get(score_column)),
      .(n_cells = .N, mean_score = mean(get(score_column)), median_score = median(get(score_column))),
      by = .(dataset, patient_key, patient_id, timepoint, treatment, response)
    ][, `:=`(score_id = score_id, analysis_compartment = compartment, score_column = score_column)]
  }
  immune_state_rows[[dataset_current]] <- rbindlist(state_list, use.names = TRUE, fill = TRUE)

  for (tier_current in c("strict_high_confidence", "expanded_any_method_positive")) {
    gate_column <- if (tier_current == "strict_high_confidence") "strict_malignant" else "expanded_malignant"
    malignant <- cells[get(gate_column) == TRUE]
    for (method_i in seq_len(nrow(method_definitions))) {
      method_current <- method_definitions$method[[method_i]]
      f1_column <- method_definitions$f1_column[[method_i]]
      tymp_column <- method_definitions$tymp_column[[method_i]]

      sample_result <- malignant[, cor_values(get(f1_column), get(tymp_column)), by = .(
        dataset, patient_key, patient_id, sample_id, object_name, timepoint,
        treatment, response, tissue_confirmed
      )]
      sample_result[, `:=`(
        unit = "sample_region", malignancy_tier = tier_current, method = method_current,
        f1_column = f1_column, tymp_column = tymp_column
      )]

      patient_time_result <- malignant[, cor_values(get(f1_column), get(tymp_column)), by = .(
        dataset, patient_key, patient_id, timepoint, treatment, response
      )]
      patient_time_result[, `:=`(
        sample_id = "<pooled>", object_name = "<pooled>", tissue_confirmed = "<pooled>",
        unit = "patient_timepoint", malignancy_tier = tier_current, method = method_current,
        f1_column = f1_column, tymp_column = tymp_column
      )]

      patient_result <- malignant[, cor_values(get(f1_column), get(tymp_column)), by = .(
        dataset, patient_key, patient_id
      )]
      patient_result[, `:=`(
        sample_id = "<pooled>", object_name = "<pooled>", timepoint = "<all>",
        treatment = "<mixed_or_pooled>", response = "<mixed_or_pooled>", tissue_confirmed = "<pooled>",
        unit = "patient_overall", malignancy_tier = tier_current, method = method_current,
        f1_column = f1_column, tymp_column = tymp_column
      )]

      row_index <- row_index + 1L
      correlation_rows[[row_index]] <- rbindlist(
        list(sample_result, patient_time_result, patient_result), use.names = TRUE, fill = TRUE
      )
    }
  }
  all_cells[[dataset_current]] <- cells
}

correlations <- rbindlist(correlation_rows, use.names = TRUE, fill = TRUE)
setcolorder(correlations, c(
  "dataset", "unit", "malignancy_tier", "method", "patient_key", "patient_id",
  "sample_id", "object_name", "timepoint", "treatment", "response", "tissue_confirmed",
  "n_cells", "spearman_rho", "p_value", "status", "f1_column", "tymp_column"
))
setorder(correlations, dataset, unit, malignancy_tier, method, patient_id, timepoint, sample_id)
atomic_fwrite(correlations, file.path(output_dir, "patient_f1_tymp_correlations_v1.tsv"))
atomic_fwrite(rbindlist(tier_count_rows, use.names = TRUE), file.path(output_dir, "malignancy_tier_cell_counts_v1.tsv"))

composition <- rbindlist(composition_rows, use.names = TRUE, fill = TRUE)
setorder(composition, dataset, unit, patient_id, timepoint, broad_annotation, fine_annotation)
atomic_fwrite(composition, file.path(output_dir, "immune_composition_by_unit_v1.tsv"))
immune_states <- rbindlist(immune_state_rows, use.names = TRUE, fill = TRUE)
setorder(immune_states, dataset, patient_id, timepoint, analysis_compartment, score_id)
atomic_fwrite(immune_states, file.path(output_dir, "immune_state_scores_patient_timepoint_v1.tsv"))

# GSE264586 paired changes relative to baseline. No response grouping is introduced.
long_source <- correlations[
  dataset == "GSE264586" & unit == "patient_timepoint" & status == "ELIGIBLE_GE20" &
    timepoint %chin% c("Baseline", "Cycle2", "End_of_treatment")
]
long_rows <- list()
long_index <- 0L
for (followup in c("Cycle2", "End_of_treatment")) {
  baseline <- long_source[timepoint == "Baseline", .(
    dataset, malignancy_tier, method, patient_key, patient_id,
    n_cells_baseline = n_cells, rho_baseline = spearman_rho
  )]
  after <- long_source[timepoint == followup, .(
    dataset, malignancy_tier, method, patient_key, patient_id,
    n_cells_followup = n_cells, rho_followup = spearman_rho
  )]
  paired <- merge(
    baseline, after,
    by = c("dataset", "malignancy_tier", "method", "patient_key", "patient_id"),
    all = FALSE
  )
  paired[, `:=`(comparison = paste0(followup, "_vs_Baseline"), delta_rho = rho_followup - rho_baseline)]
  long_index <- long_index + 1L
  long_rows[[long_index]] <- paired
}
longitudinal <- rbindlist(long_rows, use.names = TRUE, fill = TRUE)
atomic_fwrite(longitudinal, file.path(output_dir, "GSE264586_longitudinal_delta_rho_v1.tsv"))
longitudinal_summary <- longitudinal[, .(
  n_pairs = .N,
  median_rho_baseline = median(rho_baseline),
  median_rho_followup = median(rho_followup),
  median_delta_rho = median(delta_rho),
  wilcoxon_p_value = if (.N >= 2L) suppressWarnings(wilcox.test(rho_followup, rho_baseline, paired = TRUE, exact = FALSE)$p.value) else NA_real_
), by = .(malignancy_tier, method, comparison)]
atomic_fwrite(longitudinal_summary, file.path(output_dir, "GSE264586_longitudinal_summary_v1.tsv"))

# GSE294109: aggregate regions to patient before comparing ICB benefit groups.
response_patients <- correlations[
  dataset == "GSE294109" & unit == "patient_timepoint" & timepoint == "Surgery_after_treatment" &
    treatment != "Cabozantinib" & status == "ELIGIBLE_GE20" & is.finite(spearman_rho)
]
response_patients[, response_group := fifelse(response %chin% c("CR", "PR"), "benefit", fifelse(response == "No_effect", "no_effect", NA_character_))]
response_patients <- response_patients[!is.na(response_group)]
atomic_fwrite(response_patients, file.path(output_dir, "GSE294109_icb_patient_rho_v1.tsv"))
response_summary <- response_patients[, {
  benefit <- spearman_rho[response_group == "benefit"]
  no_effect <- spearman_rho[response_group == "no_effect"]
  list(
    n_benefit = length(benefit),
    n_no_effect = length(no_effect),
    median_benefit = if (length(benefit)) median(benefit) else NA_real_,
    median_no_effect = if (length(no_effect)) median(no_effect) else NA_real_,
    median_difference_benefit_minus_no_effect = if (length(benefit) && length(no_effect)) median(benefit) - median(no_effect) else NA_real_,
    wilcoxon_p_value = if (length(benefit) && length(no_effect)) suppressWarnings(wilcox.test(benefit, no_effect, exact = FALSE)$p.value) else NA_real_
  )
}, by = .(malignancy_tier, method)]
atomic_fwrite(response_summary, file.path(output_dir, "GSE294109_icb_response_comparison_v1.tsv"))

# Link malignant rho to immune-state means using one patient-level cross-section per cohort.
rho_cross <- correlations[
  unit == "patient_timepoint" & malignancy_tier == "strict_high_confidence" &
    method == "primary_gene_z" & status == "ELIGIBLE_GE20" &
    ((dataset == "GSE264586" & timepoint == "Baseline") |
       (dataset == "GSE294109" & timepoint == "Surgery_after_treatment" & treatment != "Cabozantinib")),
  .(dataset, patient_key, patient_id, timepoint, malignant_rho = spearman_rho)
]
state_join <- merge(
  rho_cross, immune_states,
  by = c("dataset", "patient_key", "patient_id", "timepoint"), all = FALSE
)
state_associations <- state_join[, {
  result <- cor_values(malignant_rho, mean_score)
  list(n_patients = result$n_cells, spearman_rho = result$spearman_rho, p_value = result$p_value, status = result$status)
}, by = .(dataset, analysis_compartment, score_id)]
atomic_fwrite(state_associations, file.path(output_dir, "malignant_rho_immune_state_associations_v1.tsv"))

broad_composition <- composition[
  unit == "patient_timepoint",
  .(fraction = sum(fraction), n_cells = sum(n_cells)),
  by = .(dataset, patient_key, patient_id, timepoint, treatment, response, broad_annotation)
]
composition_join <- merge(
  rho_cross, broad_composition,
  by = c("dataset", "patient_key", "patient_id", "timepoint"), all = FALSE
)
composition_associations <- composition_join[, {
  result <- cor_values(malignant_rho, fraction)
  list(n_patients = result$n_cells, spearman_rho = result$spearman_rho, p_value = result$p_value, status = result$status)
}, by = .(dataset, broad_annotation)]
atomic_fwrite(composition_associations, file.path(output_dir, "malignant_rho_immune_composition_associations_v1.tsv"))

status <- data.table(
  status = "PASS_PHASE60_F1_TYMP",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = length(datasets),
  correlation_rows = nrow(correlations),
  strict_eligible_patient_timepoints = correlations[
    unit == "patient_timepoint" & malignancy_tier == "strict_high_confidence" &
      method == "primary_gene_z" & status == "ELIGIBLE_GE20", .N
  ],
  expanded_eligible_patient_timepoints = correlations[
    unit == "patient_timepoint" & malignancy_tier == "expanded_any_method_positive" &
      method == "primary_gene_z" & status == "ELIGIBLE_GE20", .N
  ],
  gse264586_response_analysis = "NOT_RUN_RESPONSE_UNAVAILABLE",
  gse294109_icb_response_patients_expected = "9 benefit;3 no_effect",
  note = "Existing QC/annotation/CNV outputs reused; strict CNV is primary and any-method-positive is sensitivity."
)
atomic_fwrite(status, file.path(output_dir, "PHASE60_STATUS_v1.tsv"))

cat("PASS_PHASE60_F1_TYMP\n")
