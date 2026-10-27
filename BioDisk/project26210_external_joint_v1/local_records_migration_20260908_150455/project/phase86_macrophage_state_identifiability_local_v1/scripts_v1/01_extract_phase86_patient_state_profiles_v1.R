#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_extract_phase86_patient_state_profiles_v1.R <project_root> <phase86_root>")
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
phase_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- file.path(phase_root, "profiles_v1")
part_dir <- paste0(output_dir, ".part")
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) {
  stop("Refusing existing profiles output")
}
dir.create(part_dir, recursive = TRUE)

memory_free_gib <- function() {
  command <- "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  result <- suppressWarnings(system2(
    "powershell.exe",
    c("-NoProfile", "-Command", shQuote(command)),
    stdout = TRUE,
    stderr = TRUE
  ))
  values <- suppressWarnings(as.numeric(trimws(result)))
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  tail(values, 1L) / (1024^2)
}

free_start <- memory_free_gib()
if (!is.finite(free_start) || free_start < 15) {
  stop("Local free-memory gate failed: ", free_start, " GiB")
}

run_root <- file.path(
  project_root,
  "phase8_spatial_alignment_local_v1",
  "annotation_spatial_crosswalk_v1",
  "runs",
  "20260811_152500_full_v1"
)
objects_root <- file.path(run_root, "objects_annotated")
object_paths <- c(
  GSE207493 = file.path(objects_root, "GSE207493_spatial_annotated_v1.rds"),
  GSE264586 = file.path(objects_root, "GSE264586_spatial_annotated_v1.rds"),
  GSE294109 = file.path(objects_root, "GSE294109_spatial_annotated_v1.rds")
)
expected_bytes <- c(
  GSE207493 = 1923312451,
  GSE264586 = 941486138,
  GSE294109 = 1462703849
)

reference_root <- file.path(phase_root, "inputs_v1", "phase83_reference_v1")
coverage_path <- file.path(reference_root, "phase83_marker_coverage_v1.tsv")
reference_path <- file.path(reference_root, "phase83_reference_specificity_v1.tsv.gz")
if (!all(file.exists(c(object_paths, coverage_path, reference_path)))) {
  stop("Missing frozen input")
}

coverage <- fread(coverage_path)
coverage[, present_bool := tolower(as.character(present)) == "true"]
common <- coverage[, .(n_dataset = uniqueN(dataset), all_present = all(present_bool)), by = gene][
  n_dataset == 2L & all_present == TRUE,
  sort(unique(toupper(gene)))
]
if (length(common) != 478L) stop("Expected 478 common platform genes, found ", length(common))

reference <- fread(reference_path)
gene_column <- names(reference)[[1L]]
reference_genes <- toupper(reference[[gene_column]])
reference_states <- setdiff(names(reference), gene_column)
expected_states <- c(
  "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
  "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
  "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell"
)
if (!identical(reference_states, expected_states)) stop("Unexpected Phase83 state order")
if (!all(common %in% reference_genes)) stop("Common platform genes missing from Phase83 reference")

profile_rows <- list()
coverage_rows <- list()
status_rows <- list()
profile_i <- 0L

for (dataset in names(object_paths)) {
  free_before <- memory_free_gib()
  if (!is.finite(free_before) || free_before < 15) {
    stop("Free-memory gate failed before ", dataset, ": ", free_before, " GiB")
  }
  path <- object_paths[[dataset]]
  bytes <- file.info(path)$size
  if (!identical(as.numeric(bytes), as.numeric(expected_bytes[[dataset]]))) {
    stop("Input byte mismatch for ", dataset)
  }

  object <- readRDS(path)
  if (!inherits(object, "Seurat") || !"RNA" %in% Assays(object)) stop("Invalid Seurat object")
  metadata <- object[[]]
  required <- c("dataset", "patient_id", "phase8_spatial_reference_state")
  if (length(setdiff(required, names(metadata)))) stop("Missing Phase8 metadata in ", dataset)
  if (any(as.character(metadata$dataset) != dataset)) stop("Dataset metadata mismatch")
  if (anyNA(metadata$patient_id) || any(!nzchar(as.character(metadata$patient_id)))) stop("Missing patient IDs")

  counts <- LayerData(object[["RNA"]], layer = "counts")
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
  source_genes <- rownames(counts)
  upper <- toupper(source_genes)
  collision <- unique(upper[duplicated(upper) | duplicated(upper, fromLast = TRUE)])
  keep <- !upper %in% collision
  unique_upper <- upper[keep]
  unique_index <- which(keep)
  matched <- match(common, unique_upper)
  present <- !is.na(matched)
  if (mean(present) < 0.95) stop("Less than 95% common-gene coverage in ", dataset)

  selected <- Matrix(0, nrow = length(common), ncol = ncol(counts), sparse = TRUE)
  rownames(selected) <- common
  colnames(selected) <- colnames(counts)
  if (any(present)) {
    selected[present, ] <- counts[unique_index[matched[present]], , drop = FALSE]
  }

  state <- as.character(metadata$phase8_spatial_reference_state)
  state[state %in% c("F1-high Malignant", "F1-low Malignant")] <- "Malignant"
  patient <- as.character(metadata$patient_id)
  eligible <- state %in% expected_states
  if (!any(eligible)) stop("No eligible Phase17 states in ", dataset)

  group_table <- unique(data.table(patient = patient[eligible], state = state[eligible]))
  setorder(group_table, patient, state)
  for (j in seq_len(nrow(group_table))) {
    current_patient <- group_table$patient[[j]]
    current_state <- group_table$state[[j]]
    idx <- which(eligible & patient == current_patient & state == current_state)
    values <- as.numeric(Matrix::rowSums(selected[, idx, drop = FALSE]))
    profile_i <- profile_i + 1L
    row <- data.table(
      dataset = dataset,
      patient_id = current_patient,
      state = current_state,
      n_cells = length(idx)
    )
    row <- cbind(row, as.data.table(as.list(setNames(values, common))))
    profile_rows[[profile_i]] <- row
  }

  coverage_rows[[dataset]] <- data.table(
    dataset = dataset,
    panel_genes = length(common),
    present_genes = sum(present),
    missing_genes = paste(common[!present], collapse = ";"),
    uppercase_collision_keys = paste(sort(collision), collapse = ";")
  )
  status_rows[[dataset]] <- data.table(
    dataset = dataset,
    source_path = normalizePath(path, winslash = "/"),
    source_bytes = bytes,
    object_cells = ncol(counts),
    profile_rows = nrow(group_table),
    patients = uniqueN(patient[eligible]),
    free_memory_gib_before = free_before,
    status = "PASS_PROFILE_EXTRACTION"
  )

  rm(object, metadata, counts, selected)
  invisible(gc())
}

profiles <- rbindlist(profile_rows, use.names = TRUE, fill = TRUE)
if (anyDuplicated(profiles[, .(dataset, patient_id, state)])) stop("Duplicate patient-state profiles")
if (any(profiles$n_cells <= 0)) stop("Non-positive profile cell count")
gene_matrix <- as.matrix(profiles[, ..common])
if (any(!is.finite(gene_matrix)) || any(gene_matrix < 0) || any(abs(gene_matrix - round(gene_matrix)) > 1e-8)) {
  stop("Invalid aggregated counts")
}

fwrite(profiles, file.path(part_dir, "patient_state_profiles_v1.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(rbindlist(coverage_rows), file.path(part_dir, "profile_gene_coverage_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(status_rows), file.path(part_dir, "profile_extraction_status_v1.tsv"), sep = "\t", quote = FALSE)

state_counts <- profiles[, .(
  patients = uniqueN(patient_id),
  cells = sum(n_cells)
), by = .(dataset, state)]
setorder(state_counts, dataset, state)
fwrite(state_counts, file.path(part_dir, "patient_state_profile_counts_v1.tsv"), sep = "\t", quote = FALSE)

metadata_out <- data.table(
  record_key = c(
    "status", "started_free_memory_gib", "completed_free_memory_gib",
    "common_platform_genes", "datasets", "profile_rows", "total_cells"
  ),
  value = c(
    "PASS_PROFILE_EXTRACTION", sprintf("%.6f", free_start), sprintf("%.6f", memory_free_gib()),
    length(common), length(object_paths), nrow(profiles), sum(profiles$n_cells)
  )
)
fwrite(metadata_out, file.path(part_dir, "profile_run_metadata_v1.tsv"), sep = "\t", quote = FALSE)
writeLines("PASS_PROFILE_EXTRACTION", file.path(part_dir, "PHASE86_PROFILE_EXTRACTION_PASS_v1"))

if (!file.rename(part_dir, output_dir)) stop("Atomic profiles directory rename failed")
cat("PASS_PROFILE_EXTRACTION\n")
