#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath("D:/project26210_external_joint_v1", winslash = "/", mustWork = TRUE)
stage <- file.path(project_root, "phase14_patient_program_validation_local_v1")
provisional <- file.path(stage, "runs", "20260812_082000_provisional_v3")
sync <- file.path(project_root, "phase11_external_cnv_consensus_local_v1", "final_server_sync_v5")
out <- file.path(stage, "posthoc_fresh_external_consensus_inputs_v1")
if (dir.exists(out) || file.exists(out)) stop("Refusing existing output: ", out)

external <- fread(file.path(sync, "external_cnv_cell_consensus_final_v2.tsv.gz"), select = c(
  "exact_cell_id", "external_malignant_conservative_incremental_v1"
))
setnames(external, "exact_cell_id", "cell_id")
if (anyDuplicated(external$cell_id)) stop("Duplicate external CNV cell ID")

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
dir.create(out, recursive = TRUE)
summary_rows <- list()

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

for (dataset_current in datasets) {
  path <- file.path(provisional, dataset_current, "objects", paste0(dataset_current, "_cell_scores_compact_v1.rds"))
  cells <- as.data.table(readRDS(path))
  if (anyDuplicated(cells$cell_id) || !all(cells$dataset == dataset_current)) stop("Identity failure: ", dataset_current)
  cells <- merge(cells, external, by = "cell_id", all.x = TRUE, sort = FALSE)
  if (nrow(cells) != uniqueN(cells$cell_id)) stop("CNV merge expanded rows: ", dataset_current)
  supported <- !is.na(cells$external_malignant_conservative_incremental_v1) & cells$external_malignant_conservative_incremental_v1
  cells[, fresh_external_supported := phase7_fine_annotation == "Epithelial_tumor_like" & supported]

  selected <- cells[fresh_external_supported == TRUE]
  if (dataset_current == "GSE178481") selected <- selected[!grepl("Normal", sample_id, fixed = TRUE)]
  if (dataset_current == "GSE264586") selected <- selected[timepoint == "Baseline"]

  signature_ids <- c("f1_top100_weighted", "tymp_expression", "hypoxia_angiogenesis", "polyamine_mta_context")
  patient_rows <- list()
  for (signature_current in signature_ids) {
    primary_col <- paste0("primary__", signature_current, "__malignant_provisional")
    sensitivity_col <- paste0("sensitivity__", signature_current, "__malignant_provisional")
    value <- selected[is.finite(get(primary_col)) & is.finite(get(sensitivity_col)), .(
      n_cells = .N,
      mean_primary = mean(get(primary_col)),
      mean_sensitivity = mean(get(sensitivity_col))
    ), by = .(dataset, patient_key, patient_id)]
    value[, `:=`(signature_id = signature_current, analysis_compartment = "malignant_fresh_external_consensus")]
    patient_rows[[signature_current]] <- value
  }
  patient <- rbindlist(patient_rows, use.names = TRUE, fill = TRUE)
  dir.create(file.path(out, dataset_current), recursive = TRUE)
  atomic_fwrite(patient, file.path(out, dataset_current, "malignant_patient_program_summary_fresh_external_v1.tsv"))
  summary_rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    tumor_like_cells = sum(cells$phase7_fine_annotation == "Epithelial_tumor_like"),
    fresh_external_supported_cells = sum(cells$fresh_external_supported),
    selected_fresh_external_supported_cells = nrow(selected),
    patients_ge20 = uniqueN(patient[signature_id == "f1_top100_weighted" & n_cells >= 20, patient_key]),
    status = "PASS_POSTHOC_SENSITIVITY"
  )
  rm(cells, selected, patient)
  gc()
}

atomic_fwrite(rbindlist(summary_rows), file.path(out, "FRESH_EXTERNAL_INPUT_STATUS_v1.tsv"))

