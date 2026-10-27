#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 01_score_dataset_programs_v1.R <project_root> <dataset> <output_dir>")
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
dataset_current <- args[[2L]]
output_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)

allowed_datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
if (!dataset_current %chin% allowed_datasets) stop("Unexpected dataset: ", dataset_current)
if (dir.exists(output_dir)) stop("Refusing existing output directory: ", output_dir)

registry_dir <- file.path(project_root, "phase13_signature_lock_local_v1", "registry_v1")
registry_path <- file.path(registry_dir, "frozen_signature_registry_v1.tsv")
genes_path <- file.path(registry_dir, "frozen_signature_genes_v1.tsv")
rules_path <- file.path(registry_dir, "frozen_analysis_rules_v1.json")
phase8_run <- file.path(
  project_root, "phase8_spatial_alignment_local_v1", "annotation_spatial_crosswalk_v1",
  "runs", "20260811_152500_full_v1"
)
input_path <- file.path(phase8_run, "objects_annotated", paste0(dataset_current, "_spatial_annotated_v1.rds"))
if (!all(file.exists(registry_path, genes_path, rules_path, input_path))) {
  stop("Required Phase 13/Phase 8 input is absent")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("tables", "objects", "status", "logs")) {
  dir.create(file.path(output_dir, subdir), showWarnings = FALSE)
}

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

atomic_save_rds <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing output: ", path)
  saveRDS(x, part, compress = "gzip")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

sha256 <- function(path) {
  unname(tools::md5sum(path))
}

status_rows <- list()
add_status <- function(stage, result, detail) {
  status_rows[[length(status_rows) + 1L]] <<- data.table(
    time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    dataset = dataset_current,
    stage = stage,
    result = result,
    detail = as.character(detail)
  )
  message(sprintf("[%s] %s %s", stage, result, detail))
}

registry <- fread(registry_path)
signature_genes <- fread(genes_path)
if (anyDuplicated(registry$signature_id)) stop("Duplicate registry signature IDs")
if (anyDuplicated(signature_genes[, .(signature_id, gene)])) stop("Duplicate registry gene pairs")

add_status("load_input", "START", input_path)
input_bytes <- file.info(input_path)$size
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) stop("Input is not a Seurat object")
required_meta <- c(
  "dataset", "cell_id", "patient_key", "patient_id", "sample_id", "object_name",
  "timepoint", "treatment", "response", "tissue_confirmed",
  "phase7_broad_annotation", "phase7_fine_annotation",
  "phase8_spatial_reference_state", "phase8_spatial_display_label",
  "phase8_f1_projection_score", "phase8_malignancy_claim"
)
missing_meta <- setdiff(required_meta, colnames(object[[]]))
if (length(missing_meta)) stop("Missing metadata: ", paste(missing_meta, collapse = ","))
if (length(unique(object$dataset)) != 1L || unique(object$dataset) != dataset_current) {
  stop("Dataset identity mismatch")
}
if (anyNA(object$patient_key) || any(!nzchar(as.character(object$patient_key)))) {
  stop("Missing patient_key")
}

assay <- object[[DefaultAssay(object)]]
if (!"counts" %chin% Layers(assay)) stop("No counts layer")
counts <- LayerData(assay, layer = "counts")
if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
if (ncol(counts) != ncol(object) || !identical(colnames(counts), colnames(object))) {
  stop("Counts/object cell identity mismatch")
}
library_size <- Matrix::colSums(counts)
if (any(!is.finite(library_size)) || any(library_size <= 0)) stop("Invalid library size")
add_status("load_input", "PASS", sprintf("cells=%d features=%d bytes=%s", ncol(object), nrow(object), input_bytes))

meta <- as.data.table(object[[]], keep.rownames = "seurat_cell_name")
if (!identical(meta$seurat_cell_name, colnames(object))) stop("Metadata order mismatch")
meta[, analysis_compartment := fcase(
  phase7_fine_annotation == "Epithelial_tumor_like", "malignant_provisional",
  phase8_spatial_reference_state == "Epithelial-Normal", "normal_epithelial",
  phase7_broad_annotation == "Myeloid", "myeloid",
  phase7_broad_annotation == "Endothelial", "endothelial",
  phase7_broad_annotation %chin% c("T_cell", "NK_cell"), "t_nk",
  phase7_broad_annotation == "Fibroblast", "fibroblast",
  default = "other"
)]

requested_genes <- unique(signature_genes$gene)
present_genes <- intersect(requested_genes, rownames(counts))
target_counts <- counts[present_genes, , drop = FALSE]
target_norm <- target_counts %*% Matrix::Diagonal(x = 10000 / library_size)
target_norm@x <- log1p(target_norm@x)
dimnames(target_norm) <- dimnames(target_counts)
if (!identical(colnames(target_norm), colnames(object))) stop("Normalized target cell order mismatch")
rm(assay)
gc()

score_module <- function(cell_index, gene_table, weighted = FALSE) {
  requested <- gene_table$gene
  present <- requested[requested %chin% rownames(target_norm)]
  if (!length(present) || !length(cell_index)) {
    return(list(primary = rep(NA_real_, length(cell_index)), sensitivity = rep(NA_real_, length(cell_index)), present = present))
  }
  matrix_current <- target_norm[present, cell_index, drop = FALSE]
  sensitivity <- as.numeric(Matrix::colMeans(matrix_current))
  if (weighted) {
    weights <- as.numeric(gene_table$raw_weight[match(present, gene_table$gene)])
    weights <- weights / sum(weights)
    primary <- as.numeric(Matrix::crossprod(weights, matrix_current))
  } else if (length(present) == 1L) {
    primary <- as.numeric(matrix_current[1L, ])
  } else {
    means <- as.numeric(Matrix::rowMeans(matrix_current))
    means_sq <- as.numeric(Matrix::rowMeans(matrix_current ^ 2))
    sds <- sqrt(pmax(means_sq - means ^ 2, 0))
    inv_sd <- ifelse(is.finite(sds) & sds > 0, 1 / sds, 0)
    primary <- as.numeric(Matrix::crossprod(inv_sd, matrix_current)) / length(present) - sum(means * inv_sd) / length(present)
  }
  list(primary = primary, sensitivity = sensitivity, present = present)
}

compartment_map <- list(
  malignant = "malignant_provisional",
  myeloid = "myeloid",
  endothelial = "endothelial",
  t_nk = "t_nk",
  fibroblast = "fibroblast"
)

n_cells <- nrow(meta)
primary_scores <- list()
sensitivity_scores <- list()
summary_patient <- list()
summary_sample <- list()
availability <- list()
score_counter <- 0L

summarize_vector <- function(index, primary, sensitivity, signature_id, compartment, unit) {
  fields <- if (unit == "patient") c("dataset", "patient_key", "patient_id") else c(
    "dataset", "patient_key", "patient_id", "sample_id", "object_name", "timepoint",
    "treatment", "response", "tissue_confirmed"
  )
  frame <- meta[index, ..fields]
  frame[, `:=`(primary_score = primary, sensitivity_score = sensitivity)]
  answer <- frame[, .(
    n_cells = .N,
    mean_primary = mean(primary_score, na.rm = TRUE),
    median_primary = median(primary_score, na.rm = TRUE),
    mean_sensitivity = mean(sensitivity_score, na.rm = TRUE),
    median_sensitivity = median(sensitivity_score, na.rm = TRUE)
  ), by = fields]
  answer[, `:=`(signature_id = signature_id, analysis_compartment = compartment, unit = unit)]
  setcolorder(answer, c("dataset", "unit", setdiff(fields, "dataset"), "analysis_compartment", "signature_id", "n_cells", "mean_primary", "median_primary", "mean_sensitivity", "median_sensitivity"))
  answer
}

for (signature_current in registry$signature_id) {
  registry_row <- registry[signature_id == signature_current]
  gene_table <- signature_genes[signature_id == signature_current][order(as.integer(gene_order))]
  requested_compartments <- strsplit(registry_row$compartment, ";", fixed = TRUE)[[1L]]
  if (signature_current == "f1_top100_weighted") requested_compartments <- c(requested_compartments, "normal_epithelial")
  for (compartment_current in requested_compartments) {
    actual_compartment <- if (compartment_current == "normal_epithelial") "normal_epithelial" else compartment_map[[compartment_current]]
    if (is.null(actual_compartment)) stop("Unmapped compartment: ", compartment_current)
    index <- which(meta$analysis_compartment == actual_compartment)
    scored <- score_module(index, gene_table, weighted = signature_current == "f1_top100_weighted")
    gate <- length(scored$present) >= as.integer(registry_row$minimum_present_genes)
    availability[[length(availability) + 1L]] <- data.table(
      dataset = dataset_current,
      signature_id = signature_current,
      analysis_compartment = actual_compartment,
      requested_genes = nrow(gene_table),
      present_genes = length(scored$present),
      minimum_present_genes = as.integer(registry_row$minimum_present_genes),
      missing_genes = paste(setdiff(gene_table$gene, scored$present), collapse = ";"),
      n_cells = length(index),
      gate = ifelse(gate && length(index) > 0, "PASS", "FAIL")
    )
    if (!gate || !length(index)) next
    score_counter <- score_counter + 1L
    score_name <- paste(signature_current, actual_compartment, sep = "__")
    primary_full <- rep(NA_real_, n_cells)
    sensitivity_full <- rep(NA_real_, n_cells)
    primary_full[index] <- scored$primary
    sensitivity_full[index] <- scored$sensitivity
    primary_scores[[score_name]] <- primary_full
    sensitivity_scores[[score_name]] <- sensitivity_full
    summary_patient[[score_counter]] <- summarize_vector(index, scored$primary, scored$sensitivity, signature_current, actual_compartment, "patient")
    summary_sample[[score_counter]] <- summarize_vector(index, scored$primary, scored$sensitivity, signature_current, actual_compartment, "sample")
  }
}

availability_dt <- rbindlist(availability, fill = TRUE)
if (availability_dt[signature_id == "f1_top100_weighted" & analysis_compartment == "malignant_provisional", gate] != "PASS") {
  stop("F1 availability gate failed")
}
add_status("score_programs", "PASS", sprintf("score_contexts=%d requested_genes=%d present_genes=%d", score_counter, length(requested_genes), length(present_genes)))

patient_program <- rbindlist(summary_patient, fill = TRUE)
sample_program <- rbindlist(summary_sample, fill = TRUE)

count_fields_patient <- c("dataset", "patient_key", "patient_id", "analysis_compartment")
patient_compartment <- meta[, .(n_cells = .N), by = count_fields_patient]
patient_compartment[, fraction_all_cells := n_cells / sum(n_cells), by = .(dataset, patient_key)]
count_fields_sample <- c("dataset", "patient_key", "patient_id", "sample_id", "object_name", "timepoint", "treatment", "response", "tissue_confirmed", "analysis_compartment")
sample_compartment <- meta[, .(n_cells = .N), by = count_fields_sample]
sample_compartment[, fraction_all_cells := n_cells / sum(n_cells), by = .(dataset, patient_key, sample_id, object_name)]

state_fields_patient <- c("dataset", "patient_key", "patient_id", "phase8_spatial_reference_state", "phase8_spatial_display_label")
patient_state <- meta[, .(n_cells = .N), by = state_fields_patient]
patient_state[, fraction_patient_cells := n_cells / sum(n_cells), by = .(dataset, patient_key)]
state_fields_sample <- c("dataset", "patient_key", "patient_id", "sample_id", "object_name", "timepoint", "treatment", "response", "phase8_spatial_reference_state", "phase8_spatial_display_label")
sample_state <- meta[, .(n_cells = .N), by = state_fields_sample]
sample_state[, fraction_sample_cells := n_cells / sum(n_cells), by = .(dataset, patient_key, sample_id, object_name)]

f1_malignant_name <- "f1_top100_weighted__malignant_provisional"
f1_normal_name <- "f1_top100_weighted__normal_epithelial"
normal_values <- primary_scores[[f1_normal_name]][meta$analysis_compartment == "normal_epithelial"]
normal_values <- normal_values[is.finite(normal_values)]
thresholds <- data.table(
  dataset = dataset_current,
  normal_epithelial_cells = length(normal_values),
  threshold_type = c("normal_epithelial_q90", "normal_epithelial_q95"),
  threshold = if (length(normal_values) >= 100L) as.numeric(quantile(normal_values, c(0.90, 0.95), names = FALSE, type = 7)) else c(NA_real_, NA_real_),
  gate = ifelse(length(normal_values) >= 100L, "PASS", "FAIL_LT100_NORMAL_EPITHELIAL")
)

malignant_index <- which(meta$analysis_compartment == "malignant_provisional")
malignant_frame <- meta[malignant_index, .(
  dataset, patient_key, patient_id, sample_id, object_name, timepoint, treatment,
  response, tissue_confirmed
)]
malignant_frame[, f1_score := primary_scores[[f1_malignant_name]][malignant_index]]
prevalence_rows <- list()
for (i in seq_len(nrow(thresholds))) {
  threshold_row <- thresholds[i]
  patient_prev <- malignant_frame[, .(
    n_malignant_cells = .N,
    mean_f1 = mean(f1_score),
    median_f1 = median(f1_score),
    n_f1_above_threshold = if (is.finite(threshold_row$threshold)) sum(f1_score > threshold_row$threshold) else NA_integer_,
    fraction_f1_above_threshold = if (is.finite(threshold_row$threshold)) mean(f1_score > threshold_row$threshold) else NA_real_
  ), by = .(dataset, patient_key, patient_id)]
  patient_prev[, `:=`(unit = "patient", threshold_type = threshold_row$threshold_type, threshold = threshold_row$threshold)]
  sample_prev <- malignant_frame[, .(
    n_malignant_cells = .N,
    mean_f1 = mean(f1_score),
    median_f1 = median(f1_score),
    n_f1_above_threshold = if (is.finite(threshold_row$threshold)) sum(f1_score > threshold_row$threshold) else NA_integer_,
    fraction_f1_above_threshold = if (is.finite(threshold_row$threshold)) mean(f1_score > threshold_row$threshold) else NA_real_
  ), by = .(dataset, patient_key, patient_id, sample_id, object_name, timepoint, treatment, response, tissue_confirmed)]
  sample_prev[, `:=`(unit = "sample", threshold_type = threshold_row$threshold_type, threshold = threshold_row$threshold)]
  prevalence_rows[[length(prevalence_rows) + 1L]] <- patient_prev
  prevalence_rows[[length(prevalence_rows) + 1L]] <- sample_prev
}
f1_prevalence <- rbindlist(prevalence_rows, fill = TRUE)

phase8_f1 <- as.numeric(meta$phase8_f1_projection_score[malignant_index])
new_f1 <- as.numeric(primary_scores[[f1_malignant_name]][malignant_index])
finite_f1 <- is.finite(phase8_f1) & is.finite(new_f1)
f1_concordance <- data.table(
  dataset = dataset_current,
  n_compared = sum(finite_f1),
  max_abs_difference = if (any(finite_f1)) max(abs(phase8_f1[finite_f1] - new_f1[finite_f1])) else NA_real_,
  pearson_correlation = if (sum(finite_f1) >= 3L) cor(phase8_f1[finite_f1], new_f1[finite_f1]) else NA_real_,
  gate = ifelse(any(finite_f1) && max(abs(phase8_f1[finite_f1] - new_f1[finite_f1])) <= 1e-8, "PASS", "FAIL")
)
if (f1_concordance$gate != "PASS") stop("F1 score does not reproduce Phase 8")

cell_keep <- meta$analysis_compartment != "other"
cell_scores <- meta[cell_keep, .(
  dataset, cell_id, seurat_cell_name, patient_key, patient_id, sample_id,
  object_name, timepoint, treatment, response, tissue_confirmed,
  phase7_broad_annotation, phase7_fine_annotation,
  phase8_spatial_reference_state, phase8_spatial_display_label,
  phase8_malignancy_claim, analysis_compartment
)]
for (score_name in names(primary_scores)) {
  cell_scores[[paste0("primary__", score_name)]] <- primary_scores[[score_name]][cell_keep]
  cell_scores[[paste0("sensitivity__", score_name)]] <- sensitivity_scores[[score_name]][cell_keep]
}

tables <- file.path(output_dir, "tables")
atomic_fwrite(availability_dt, file.path(tables, "signature_gene_availability_v1.tsv"))
atomic_fwrite(patient_program, file.path(tables, "patient_program_summary_v1.tsv"))
atomic_fwrite(sample_program, file.path(tables, "sample_program_summary_v1.tsv"))
atomic_fwrite(patient_compartment, file.path(tables, "patient_compartment_counts_v1.tsv"))
atomic_fwrite(sample_compartment, file.path(tables, "sample_compartment_counts_v1.tsv"))
atomic_fwrite(patient_state, file.path(tables, "patient_phase8_state_counts_v1.tsv"))
atomic_fwrite(sample_state, file.path(tables, "sample_phase8_state_counts_v1.tsv"))
atomic_fwrite(thresholds, file.path(tables, "f1_normal_epithelial_thresholds_v1.tsv"))
atomic_fwrite(f1_prevalence, file.path(tables, "f1_malignant_prevalence_v1.tsv"))
atomic_fwrite(f1_concordance, file.path(tables, "f1_phase8_concordance_v1.tsv"))
atomic_save_rds(cell_scores, file.path(output_dir, "objects", paste0(dataset_current, "_cell_scores_compact_v1.rds")))

rm(object, counts, target_counts, target_norm, primary_scores, sensitivity_scores, cell_scores)
gc()

add_status("write_outputs", "PASS", sprintf("patient_program_rows=%d sample_program_rows=%d", nrow(patient_program), nrow(sample_program)))
add_status("OVERALL", "PASS_PROVISIONAL_MALIGNANCY", "Non-malignant compartments final; malignant summaries await external CNV consensus")
status <- rbindlist(status_rows)
atomic_fwrite(status, file.path(output_dir, "status", paste0(dataset_current, "_STATUS_v1.tsv")))

manifest_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[!grepl("FILE_MANIFEST", manifest_files, fixed = TRUE)]
manifest <- data.table(
  relative_path = substring(normalizePath(manifest_files, winslash = "/"), nchar(normalizePath(output_dir, winslash = "/")) + 2L),
  bytes = file.info(manifest_files)$size,
  md5 = unname(tools::md5sum(manifest_files))
)
atomic_fwrite(manifest, file.path(output_dir, "FILE_MANIFEST_MD5_v1.tsv"))
