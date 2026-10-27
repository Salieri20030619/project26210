#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 01_prepare_phase88_cellchat_inputs_v1.R <project_root> <phase88_root>")
project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
phase_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- file.path(phase_root, "prepared_inputs_v1")
part_dir <- paste0(output_dir, ".part")
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop("Refusing existing Phase88 prepared input")
dir.create(part_dir, recursive = TRUE)

memory_free_gib <- function() {
  command <- "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  result <- suppressWarnings(system2("powershell.exe", c("-NoProfile", "-Command", shQuote(command)), stdout = TRUE, stderr = TRUE))
  values <- suppressWarnings(as.numeric(trimws(result)))
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  tail(values, 1L) / (1024^2)
}

objects_root <- file.path(
  project_root, "phase8_spatial_alignment_local_v1", "annotation_spatial_crosswalk_v1",
  "runs", "20260811_152500_full_v1", "objects_annotated"
)
object_paths <- c(
  GSE264586 = file.path(objects_root, "GSE264586_spatial_annotated_v1.rds"),
  GSE294109 = file.path(objects_root, "GSE294109_spatial_annotated_v1.rds"),
  GSE207493 = file.path(objects_root, "GSE207493_spatial_annotated_v1.rds")
)
expected_bytes <- c(GSE207493 = 1923312451, GSE264586 = 941486138, GSE294109 = 1462703849)
db_path <- file.path(phase_root, "inputs_v1", "cellchat_resource_v1", "CellChatDB.human_v1.rds")
if (!all(file.exists(c(object_paths, db_path)))) stop("Missing frozen Phase88 input")

db <- readRDS(db_path)
if (!all(c("interaction", "complex", "cofactor") %in% names(db))) stop("Invalid CellChatDB")
interaction <- as.data.frame(db$interaction)
simple <- unique(c(as.character(interaction$ligand), as.character(interaction$receptor)))
complex_members <- unique(as.character(unlist(db$complex, use.names = FALSE)))
cofactor_members <- unique(as.character(unlist(db$cofactor, use.names = FALSE)))
required_tokens <- sort(unique(toupper(c(simple, complex_members, cofactor_members))))
required_tokens <- required_tokens[!is.na(required_tokens) & nzchar(required_tokens)]
fwrite(data.table(gene = required_tokens), file.path(part_dir, "cellchat_required_genes_regenerated_v1.tsv"), sep = "\t", quote = FALSE)

groups <- c(
  "Malignant", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage",
  "Monocyte", "T cell", "NK cell", "Endothelial"
)
required_unit_groups <- c("Malignant", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage")
cap <- 200L
seed <- 2621088L

deterministic_sample <- function(ids, cap_value) {
  ids <- as.character(ids)
  if (length(ids) <= cap_value) return(sort(ids))
  hashes <- vapply(ids, function(value) digest(paste0(seed, "|", value), algo = "sha256", serialize = FALSE), character(1))
  ids[order(hashes, ids)][seq_len(cap_value)]
}

status_rows <- list()
unit_rows <- list()
sampling_rows <- list()

for (dataset in names(object_paths)) {
  free_before <- memory_free_gib()
  if (!is.finite(free_before) || free_before < 15) stop("Memory gate failed before ", dataset, ": ", free_before)
  path <- object_paths[[dataset]]
  if (as.numeric(file.info(path)$size) != as.numeric(expected_bytes[[dataset]])) stop("Input byte mismatch: ", dataset)
  object <- readRDS(path)
  if (!inherits(object, "Seurat") || !"RNA" %in% Assays(object)) stop("Invalid Seurat object")
  metadata <- object[[]]
  required_meta <- c("patient_id", "sample_id", "phase8_spatial_reference_state", "phase8_tymp_log_normalized")
  if (length(setdiff(required_meta, names(metadata)))) stop("Missing Phase88 metadata: ", dataset)
  state <- as.character(metadata$phase8_spatial_reference_state)
  group <- state
  group[state %in% c("F1-high Malignant", "F1-low Malignant")] <- "Malignant"
  patient <- as.character(metadata$patient_id)
  sample_id <- as.character(metadata$sample_id)
  unit_id <- paste(dataset, patient, sample_id, sep = "__")
  eligible_cell <- group %in% groups

  count_table <- data.table(unit_id = unit_id[eligible_cell], patient_id = patient[eligible_cell], sample_id = sample_id[eligible_cell], communication_group = group[eligible_cell])[, .N, by = .(unit_id, patient_id, sample_id, communication_group)]
  wide <- dcast(count_table, unit_id + patient_id + sample_id ~ communication_group, value.var = "N", fill = 0)
  for (name in groups) if (!name %in% names(wide)) wide[[name]] <- 0L
  eligible_units <- wide[
    Malignant >= 20L &
      get("SPP1-state Macrophage") >= 10L &
      get("C1QC-state Macrophage") >= 10L &
      get("Other Macrophage") >= 10L,
    unit_id
  ]
  if (!length(eligible_units)) stop("No eligible Phase88 units: ", dataset)

  selected_ids <- character()
  for (current_unit in sort(eligible_units)) {
    full_unit <- which(unit_id == current_unit & eligible_cell)
    malignant_full <- full_unit[group[full_unit] == "Malignant"]
    malignant_tymp <- as.numeric(metadata$phase8_tymp_log_normalized[malignant_full])
    if (length(malignant_tymp) < 20L || any(!is.finite(malignant_tymp))) stop("Invalid malignant TYMP exposure: ", current_unit)
    unit_row <- wide[unit_id == current_unit]
    unit_rows[[length(unit_rows) + 1L]] <- data.table(
      dataset = dataset,
      unit_id = current_unit,
      patient_id = unit_row$patient_id[[1L]],
      sample_id = unit_row$sample_id[[1L]],
      malignant_tymp_median = median(malignant_tymp),
      malignant_cells_full = length(malignant_full),
      spp1_cells_full = unit_row[["SPP1-state Macrophage"]][[1L]],
      c1qc_cells_full = unit_row[["C1QC-state Macrophage"]][[1L]],
      other_macrophage_cells_full = unit_row[["Other Macrophage"]][[1L]]
    )
    for (current_group in groups) {
      ids <- rownames(metadata)[which(unit_id == current_unit & eligible_cell & !is.na(group) & group == current_group)]
      if (length(ids) < 10L) next
      chosen <- deterministic_sample(ids, cap)
      selected_ids <- c(selected_ids, chosen)
      sampling_rows[[length(sampling_rows) + 1L]] <- data.table(
        dataset = dataset,
        unit_id = current_unit,
        patient_id = unit_row$patient_id[[1L]],
        sample_id = unit_row$sample_id[[1L]],
        communication_group = current_group,
        available_cells = length(ids),
        sampled_cells = length(chosen),
        cap = cap
      )
    }
  }
  if (anyDuplicated(selected_ids)) stop("Duplicate sampled cell IDs: ", dataset)

  counts <- LayerData(object[["RNA"]], layer = "counts")
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
  source_genes <- rownames(counts)
  upper <- toupper(source_genes)
  collision <- unique(upper[duplicated(upper) | duplicated(upper, fromLast = TRUE)])
  keep <- !upper %in% collision
  unique_upper <- upper[keep]
  unique_index <- which(keep)
  present_tokens <- intersect(required_tokens, unique_upper)
  if (length(present_tokens) < 0.75 * length(required_tokens)) stop("Less than 75% CellChat token coverage: ", dataset)
  gene_index <- unique_index[match(present_tokens, unique_upper)]
  cell_index <- match(selected_ids, colnames(counts))
  if (anyNA(cell_index)) stop("Sampled cell join failure: ", dataset)
  libraries <- as.numeric(Matrix::colSums(counts[, cell_index, drop = FALSE]))
  if (any(!is.finite(libraries)) || any(libraries <= 0)) stop("Invalid library sizes: ", dataset)
  selected_counts <- counts[gene_index, cell_index, drop = FALSE]
  rownames(selected_counts) <- present_tokens
  colnames(selected_counts) <- selected_ids
  normalized <- selected_counts %*% Diagonal(x = 10000 / libraries)
  normalized@x <- log1p(normalized@x)
  if (any(!is.finite(normalized@x)) || any(normalized@x < 0)) stop("Invalid normalized CellChat data: ", dataset)

  selected_meta <- data.table(
    cell_id = selected_ids,
    dataset = dataset,
    unit_id = unit_id[match(selected_ids, rownames(metadata))],
    patient_id = patient[match(selected_ids, rownames(metadata))],
    sample_id = sample_id[match(selected_ids, rownames(metadata))],
    communication_group = group[match(selected_ids, rownames(metadata))]
  )
  dataset_dir <- file.path(part_dir, dataset)
  dir.create(dataset_dir)
  saveRDS(normalized, file.path(dataset_dir, "cellchat_log1p_cp10k_v1.rds"), version = 3)
  fwrite(selected_meta, file.path(dataset_dir, "cellchat_cell_metadata_v1.tsv"), sep = "\t", quote = FALSE)
  fwrite(data.table(gene = present_tokens), file.path(dataset_dir, "cellchat_genes_v1.tsv"), sep = "\t", quote = FALSE)
  status_rows[[dataset]] <- data.table(
    dataset = dataset,
    source_path = normalizePath(path, winslash = "/"),
    source_bytes = file.info(path)$size,
    source_cells = ncol(object),
    eligible_units = length(eligible_units),
    sampled_cells = length(selected_ids),
    cellchat_genes = length(present_tokens),
    required_tokens = length(required_tokens),
    free_memory_gib_before = free_before,
    status = "PASS_PHASE88_INPUT_PREPARATION"
  )
  rm(object, metadata, counts, selected_counts, normalized, selected_meta)
  invisible(gc())
}

units <- rbindlist(unit_rows)
sampling <- rbindlist(sampling_rows)
if (anyDuplicated(units$unit_id)) stop("Duplicate unit exposure rows")
fwrite(units, file.path(part_dir, "phase88_unit_exposure_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(sampling, file.path(part_dir, "phase88_sampling_qc_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(status_rows), file.path(part_dir, "phase88_input_status_v1.tsv"), sep = "\t", quote = FALSE)
writeLines("PASS_PHASE88_INPUT_PREPARATION", file.path(part_dir, "PHASE88_INPUT_PREPARATION_PASS_v1"))
if (!file.rename(part_dir, output_dir)) stop("Atomic Phase88 prepared-input rename failed")
cat("PASS_PHASE88_INPUT_PREPARATION\n")
