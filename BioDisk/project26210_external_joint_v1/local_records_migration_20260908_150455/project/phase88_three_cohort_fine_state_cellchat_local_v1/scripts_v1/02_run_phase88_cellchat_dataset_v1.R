#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 02_run_phase88_cellchat_dataset_v1.R <dataset> <prepared_root> <db_rds> <output_dir>")
dataset <- args[[1L]]
prepared_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
db_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
part_dir <- paste0(output_dir, ".part")
if (!dataset %in% c("GSE207493", "GSE264586", "GSE294109")) stop("Unexpected dataset")
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop("Refusing existing CellChat output")
dir.create(part_dir, recursive = TRUE)
dir.create(file.path(part_dir, "unit_rds"))

input_dir <- file.path(prepared_root, dataset)
expression <- readRDS(file.path(input_dir, "cellchat_log1p_cp10k_v1.rds"))
meta <- fread(file.path(input_dir, "cellchat_cell_metadata_v1.tsv"))
dataset_name <- dataset
exposure <- fread(file.path(prepared_root, "phase88_unit_exposure_v1.tsv"))[get("dataset") == dataset_name]
if (!inherits(expression, "sparseMatrix") || ncol(expression) != nrow(meta)) stop("Expression/metadata mismatch")
if (anyDuplicated(meta$cell_id) || anyDuplicated(exposure$unit_id)) stop("Duplicate Phase88 identifiers")
if (is.null(colnames(expression))) colnames(expression) <- as.character(meta$cell_id)
if (!identical(as.character(colnames(expression)), as.character(meta$cell_id))) stop("Cell order mismatch")
db <- readRDS(db_path)
future::plan("sequential")

macrophages <- c("SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage")
all_interactions <- list()
endpoint_rows <- list()
qc_rows <- list()

edge_sum <- function(tab, sources, targets) {
  if (!nrow(tab) || !length(sources) || !length(targets)) return(0)
  value <- tab$prob[tab$source %in% sources & tab$target %in% targets]
  if (!length(value)) return(0)
  sum(value, na.rm = TRUE)
}

specificity <- function(spp1, c1qc, other) log1p(spp1) - mean(c(log1p(c1qc), log1p(other)))

for (unit in sort(unique(meta$unit_id))) {
  message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", dataset, " ", unit)
  index <- which(meta$unit_id == unit)
  unit_meta <- as.data.frame(meta[index])
  rownames(unit_meta) <- unit_meta$cell_id
  present <- sort(unique(unit_meta$communication_group))
  required <- c("Malignant", macrophages)
  if (!all(required %in% present)) stop("Required fine-state group missing: ", unit)
  cellchat <- createCellChat(object = expression[, index, drop = FALSE], meta = unit_meta, group.by = "communication_group")
  cellchat@DB <- db
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = FALSE, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  tab <- subsetCommunication(cellchat, thresh = 1)
  if (is.null(tab) || !nrow(tab)) {
    tab <- data.frame(
      source = character(), target = character(), ligand = character(), receptor = character(),
      interaction_name = character(), pathway_name = character(), prob = numeric(), pval = numeric()
    )
  }
  tab$dataset <- dataset
  tab$unit_id <- unit
  tab$patient_id <- unit_meta$patient_id[[1L]]
  tab$sample_id <- unit_meta$sample_id[[1L]]
  all_interactions[[unit]] <- tab

  m_to_spp1 <- edge_sum(tab, "Malignant", "SPP1-state Macrophage")
  m_to_c1qc <- edge_sum(tab, "Malignant", "C1QC-state Macrophage")
  m_to_other <- edge_sum(tab, "Malignant", "Other Macrophage")
  spp1_to_m <- edge_sum(tab, "SPP1-state Macrophage", "Malignant")
  c1qc_to_m <- edge_sum(tab, "C1QC-state Macrophage", "Malignant")
  other_to_m <- edge_sum(tab, "Other Macrophage", "Malignant")
  spp1_to_t <- if ("T cell" %in% present) edge_sum(tab, "SPP1-state Macrophage", "T cell") else NA_real_
  c1qc_to_t <- if ("T cell" %in% present) edge_sum(tab, "C1QC-state Macrophage", "T cell") else NA_real_
  other_to_t <- if ("T cell" %in% present) edge_sum(tab, "Other Macrophage", "T cell") else NA_real_
  spp1_to_e <- if ("Endothelial" %in% present) edge_sum(tab, "SPP1-state Macrophage", "Endothelial") else NA_real_
  c1qc_to_e <- if ("Endothelial" %in% present) edge_sum(tab, "C1QC-state Macrophage", "Endothelial") else NA_real_
  other_to_e <- if ("Endothelial" %in% present) edge_sum(tab, "Other Macrophage", "Endothelial") else NA_real_
  spp1_out <- edge_sum(tab, "SPP1-state Macrophage", present)
  c1qc_out <- edge_sum(tab, "C1QC-state Macrophage", present)
  other_out <- edge_sum(tab, "Other Macrophage", present)
  spp1_in <- edge_sum(tab, present, "SPP1-state Macrophage")
  c1qc_in <- edge_sum(tab, present, "C1QC-state Macrophage")
  other_in <- edge_sum(tab, present, "Other Macrophage")

  values <- list(
    P01_Malignant_to_SPP1_specificity = specificity(m_to_spp1, m_to_c1qc, m_to_other),
    P02_SPP1_to_Malignant_specificity = specificity(spp1_to_m, c1qc_to_m, other_to_m),
    P03_SPP1_to_Tcell_specificity = if (all(is.finite(c(spp1_to_t, c1qc_to_t, other_to_t)))) specificity(spp1_to_t, c1qc_to_t, other_to_t) else NA_real_,
    P04_SPP1_to_Endothelial_specificity = if (all(is.finite(c(spp1_to_e, c1qc_to_e, other_to_e)))) specificity(spp1_to_e, c1qc_to_e, other_to_e) else NA_real_,
    P05_SPP1_total_outgoing_specificity = specificity(spp1_out, c1qc_out, other_out),
    P06_SPP1_total_incoming_specificity = specificity(spp1_in, c1qc_in, other_in),
    A01_Malignant_to_SPP1_absolute = m_to_spp1,
    A02_SPP1_to_Malignant_absolute = spp1_to_m,
    A03_SPP1_to_Tcell_absolute = spp1_to_t,
    A04_SPP1_to_Endothelial_absolute = spp1_to_e,
    A05_SPP1_total_outgoing_absolute = spp1_out,
    A06_SPP1_total_incoming_absolute = spp1_in
  )
  exposure_row <- exposure[unit_id == unit]
  if (nrow(exposure_row) != 1L) stop("Unit exposure mismatch: ", unit)
  for (endpoint in names(values)) {
    endpoint_rows[[length(endpoint_rows) + 1L]] <- data.table(
      dataset = dataset,
      unit_id = unit,
      patient_id = unit_meta$patient_id[[1L]],
      sample_id = unit_meta$sample_id[[1L]],
      endpoint_id = endpoint,
      tier = ifelse(startsWith(endpoint, "P"), "PRIMARY", "SECONDARY_ABSOLUTE"),
      malignant_tymp_median = exposure_row$malignant_tymp_median[[1L]],
      endpoint_value = as.numeric(values[[endpoint]])
    )
  }
  qc_rows[[unit]] <- data.table(
    dataset = dataset,
    unit_id = unit,
    patient_id = unit_meta$patient_id[[1L]],
    sample_id = unit_meta$sample_id[[1L]],
    sampled_cells = length(index),
    groups = length(present),
    group_names = paste(present, collapse = ";"),
    interaction_rows = nrow(tab),
    status = "PASS_CELLCHAT_UNIT"
  )
  saveRDS(cellchat, file.path(part_dir, "unit_rds", paste0(gsub("[^A-Za-z0-9]+", "_", unit), "_cellchat_v1.rds")), version = 3)
  rm(cellchat)
  invisible(gc())
}

interactions <- rbindlist(lapply(all_interactions, as.data.table), fill = TRUE)
endpoints <- rbindlist(endpoint_rows)
qc <- rbindlist(qc_rows)
fwrite(interactions, file.path(part_dir, paste0(dataset, "_cellchat_interactions_v1.tsv.gz")), sep = "\t", quote = FALSE)
fwrite(endpoints, file.path(part_dir, paste0(dataset, "_unit_endpoints_v1.tsv")), sep = "\t", quote = FALSE)
fwrite(qc, file.path(part_dir, paste0(dataset, "_unit_qc_v1.tsv")), sep = "\t", quote = FALSE)
status <- list(
  status = "PASS_PHASE88_CELLCHAT_DATASET",
  dataset = dataset,
  units = nrow(qc),
  sampled_cells = sum(qc$sampled_cells),
  interaction_rows = nrow(interactions),
  endpoint_rows = nrow(endpoints),
  cellchat_version = as.character(packageVersion("CellChat"))
)
write_json(status, file.path(part_dir, paste0(dataset, "_status_v1.json")), pretty = TRUE, auto_unbox = TRUE)
writeLines("PASS_PHASE88_CELLCHAT_DATASET", file.path(part_dir, paste0(dataset, "_PASS_v1")))
if (!file.rename(part_dir, output_dir)) stop("Atomic dataset CellChat rename failed")
cat("PASS_PHASE88_CELLCHAT_DATASET\n")
