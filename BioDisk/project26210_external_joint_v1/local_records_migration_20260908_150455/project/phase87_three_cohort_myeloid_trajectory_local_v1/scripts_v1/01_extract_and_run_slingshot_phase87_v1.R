#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(SingleCellExperiment)
  library(slingshot)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 01_extract_and_run_slingshot_phase87_v1.R <project_root> <phase87_root>")
project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
phase_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- file.path(phase_root, "trajectory_inputs_v1")
part_dir <- paste0(output_dir, ".part")
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop("Refusing existing output")
dir.create(part_dir, recursive = TRUE)

memory_free_gib <- function() {
  command <- "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  result <- suppressWarnings(system2("powershell.exe", c("-NoProfile", "-Command", shQuote(command)), stdout = TRUE, stderr = TRUE))
  values <- suppressWarnings(as.numeric(trimws(result)))
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  tail(values, 1L) / (1024^2)
}

free_start <- memory_free_gib()
if (!is.finite(free_start) || free_start < 15) stop("Free-memory gate failed: ", free_start)

objects_root <- file.path(
  project_root, "phase8_spatial_alignment_local_v1", "annotation_spatial_crosswalk_v1",
  "runs", "20260811_152500_full_v1", "objects_annotated"
)
object_paths <- c(
  GSE207493 = file.path(objects_root, "GSE207493_spatial_annotated_v1.rds"),
  GSE264586 = file.path(objects_root, "GSE264586_spatial_annotated_v1.rds"),
  GSE294109 = file.path(objects_root, "GSE294109_spatial_annotated_v1.rds")
)
expected_bytes <- c(GSE207493 = 1923312451, GSE264586 = 941486138, GSE294109 = 1462703849)
if (!all(file.exists(object_paths))) stop("Missing Phase8 object")

marker_path <- file.path(phase_root, "inputs_v1", "phase83_marker_members_v1.tsv")
markers <- fread(marker_path)
anchors <- c("SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL", "C1QA", "C1QB", "C1QC", "APOE", "TYMP")
marker_panel <- function(state_name) {
  values <- markers[state == state_name & !toupper(gene) %chin% anchors][order(rank), toupper(gene)]
  unique(values)[seq_len(min(20L, uniqueN(values)))]
}
programs <- list(
  SPP1_MARKER20_ANCHOR_EXCLUDED = marker_panel("SPP1-state Macrophage"),
  C1QC_MARKER20_ANCHOR_EXCLUDED = marker_panel("C1QC-state Macrophage"),
  OTHER_MARKER20_ANCHOR_EXCLUDED = marker_panel("Other Macrophage"),
  MONOCYTE_PROGRAM = c("FCN1", "VCAN", "S100A8", "S100A9", "CTSS", "LILRB1"),
  IMMUNOSUPPRESSION_NO_IL10 = c("TGFB1", "VSIG4", "MARCO", "CD274", "IDO1")
)
if (any(lengths(programs[1:3]) != 20L)) stop("Expected 20 anchor-excluded markers per macrophage state")

eligible_states <- c("Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage")
root_rows <- list()
status_rows <- list()
coverage_rows <- list()

for (dataset in names(object_paths)) {
  free_before <- memory_free_gib()
  if (!is.finite(free_before) || free_before < 15) stop("Memory gate failed before ", dataset, ": ", free_before)
  path <- object_paths[[dataset]]
  if (as.numeric(file.info(path)$size) != as.numeric(expected_bytes[[dataset]])) stop("Input byte mismatch: ", dataset)
  object <- readRDS(path)
  if (!inherits(object, "Seurat") || !"RNA" %in% Assays(object) || !"pca" %in% Reductions(object)) stop("Invalid Seurat trajectory input")
  metadata <- object[[]]
  required <- c("patient_id", "sample_id", "phase8_spatial_reference_state", "phase7_cluster_res0_8")
  if (length(setdiff(required, names(metadata)))) stop("Missing trajectory metadata: ", dataset)
  state_all <- as.character(metadata$phase8_spatial_reference_state)
  selected_cells <- rownames(metadata)[state_all %in% eligible_states]
  if (length(selected_cells) < 500L) stop("Too few trajectory cells: ", dataset)
  metadata <- metadata[selected_cells, , drop = FALSE]
  state <- as.character(metadata$phase8_spatial_reference_state)
  patient <- as.character(metadata$patient_id)
  sample_id <- as.character(metadata$sample_id)
  cluster <- as.character(metadata$phase7_cluster_res0_8)

  pca_all <- Embeddings(object[["pca"]])
  if (!all(selected_cells %in% rownames(pca_all)) || ncol(pca_all) < 30L) stop("PCA mismatch: ", dataset)
  pca <- pca_all[selected_cells, seq_len(30L), drop = FALSE]
  if (any(!is.finite(pca))) stop("Non-finite PCA: ", dataset)

  root_qc <- data.table(cluster = cluster, state = state)[, .(
    total_cells = .N,
    monocyte_cells = sum(state == "Monocyte")
  ), by = cluster]
  root_qc[, monocyte_fraction := monocyte_cells / total_cells]
  candidates <- root_qc[total_cells >= 50L & monocyte_cells > 0L]
  if (!nrow(candidates)) stop("No eligible root cluster: ", dataset)
  setorder(candidates, -monocyte_fraction, -monocyte_cells, -total_cells, cluster)
  root_cluster <- candidates$cluster[[1L]]
  root_indices <- which(cluster == root_cluster & state == "Monocyte")
  centroid <- colMeans(pca[root_indices, , drop = FALSE])
  root_distance <- rowSums((pca[root_indices, , drop = FALSE] - matrix(centroid, nrow = length(root_indices), ncol = ncol(pca), byrow = TRUE))^2)
  root_index <- root_indices[[which.min(root_distance)]]
  root_cell <- selected_cells[[root_index]]

  dummy <- sparseMatrix(i = integer(), j = integer(), dims = c(1L, length(selected_cells)))
  rownames(dummy) <- "DUMMY"
  colnames(dummy) <- selected_cells
  sce <- SingleCellExperiment(assays = list(counts = dummy), colData = DataFrame(cluster = cluster))
  reducedDims(sce) <- SimpleList(PCA = pca)
  sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "PCA", start.clus = root_cluster, stretch = 0)
  pseudo <- slingPseudotime(sce, na = FALSE)
  weights <- slingCurveWeights(sce)
  if (is.null(dim(pseudo))) pseudo <- matrix(pseudo, ncol = 1L)
  if (is.null(dim(weights))) weights <- matrix(weights, ncol = 1L)
  valid <- is.finite(pseudo) & is.finite(weights) & weights > 0
  numerator <- rowSums(ifelse(valid, pseudo * weights, 0))
  denominator <- rowSums(ifelse(valid, weights, 0))
  slingshot_pt <- ifelse(denominator > 0, numerator / denominator, NA_real_)
  finite_pt <- is.finite(slingshot_pt)
  if (sum(finite_pt) < 0.8 * length(slingshot_pt)) stop("Less than 80% finite Slingshot pseudotime: ", dataset)
  range_pt <- range(slingshot_pt[finite_pt])
  if (diff(range_pt) <= 0) stop("Degenerate Slingshot pseudotime: ", dataset)
  slingshot_pt[finite_pt] <- (slingshot_pt[finite_pt] - range_pt[[1L]]) / diff(range_pt)

  expression_all <- LayerData(object[["RNA"]], layer = "data")
  source_genes <- rownames(expression_all)
  upper <- toupper(source_genes)
  collision <- unique(upper[duplicated(upper) | duplicated(upper, fromLast = TRUE)])
  keep <- !upper %in% collision
  unique_upper <- upper[keep]
  unique_index <- which(keep)
  scores <- list()
  for (program_name in names(programs)) {
    frozen_genes <- unique(toupper(programs[[program_name]]))
    matched <- match(frozen_genes, unique_upper)
    present <- !is.na(matched)
    coverage <- mean(present)
    if (coverage < 0.70) stop("Program coverage below 70%: ", dataset, "/", program_name)
    values <- expression_all[unique_index[matched[present]], selected_cells, drop = FALSE]
    scores[[program_name]] <- as.numeric(Matrix::colMeans(values))
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      dataset = dataset,
      program = program_name,
      frozen_genes = length(frozen_genes),
      present_genes = sum(present),
      genes = paste(frozen_genes[present], collapse = ";"),
      coverage_fraction = coverage
    )
  }

  cell_table <- data.table(
    cell_id = selected_cells,
    dataset = dataset,
    patient_id = patient,
    sample_id = sample_id,
    phase8_state = state,
    phase7_cluster = cluster,
    root_cluster = root_cluster,
    root_cell = selected_cells == root_cell,
    slingshot_pseudotime = slingshot_pt
  )
  for (index in seq_len(ncol(pca))) cell_table[[paste0("PC_", index)]] <- pca[, index]
  for (program_name in names(scores)) cell_table[[program_name]] <- scores[[program_name]]
  fwrite(cell_table, file.path(part_dir, paste0(dataset, "_phase87_cells_v1.tsv.gz")), sep = "\t", quote = FALSE)

  root_rows[[dataset]] <- cbind(data.table(
    dataset = dataset,
    root_cluster = root_cluster,
    root_cell = root_cell,
    trajectory_cells = length(selected_cells),
    monocyte_cells = sum(state == "Monocyte"),
    slingshot_lineages = ncol(pseudo),
    finite_slingshot_cells = sum(finite_pt)
  ), candidates[1L, .(root_cluster_total = total_cells, root_cluster_monocytes = monocyte_cells, root_cluster_monocyte_fraction = monocyte_fraction)])
  status_rows[[dataset]] <- data.table(
    dataset = dataset,
    source_path = normalizePath(path, winslash = "/"),
    source_bytes = file.info(path)$size,
    source_cells = ncol(object),
    trajectory_cells = length(selected_cells),
    patients = uniqueN(patient),
    samples = uniqueN(sample_id),
    clusters = uniqueN(cluster),
    free_memory_gib_before = free_before,
    status = "PASS_SLINGSHOT_EXTRACTION"
  )

  rm(object, metadata, pca_all, pca, expression_all, sce, pseudo, weights, dummy, cell_table)
  invisible(gc())
}

fwrite(rbindlist(root_rows), file.path(part_dir, "phase87_root_qc_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(status_rows), file.path(part_dir, "phase87_dataset_status_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(coverage_rows), file.path(part_dir, "phase87_program_coverage_v1.tsv"), sep = "\t", quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(part_dir, "R_SESSION_INFO_v1.txt"))
writeLines("PASS_SLINGSHOT_EXTRACTION", file.path(part_dir, "PHASE87_SLINGSHOT_EXTRACTION_PASS_v1"))
if (!file.rename(part_dir, output_dir)) stop("Atomic trajectory-input directory rename failed")
cat("PASS_SLINGSHOT_EXTRACTION\n")
