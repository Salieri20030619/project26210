#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 01b_run_slingshot_phase87_v2.R <phase87_root>")
phase_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
profiles_root <- file.path(phase_root, "trajectory_profiles_v2")
output_dir <- file.path(phase_root, "trajectory_inputs_v1")
part_dir <- paste0(output_dir, ".part")
if (!dir.exists(profiles_root) || dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop("Invalid Slingshot v2 directory state")
dir.create(part_dir, recursive = TRUE)

datasets <- c("GSE207493", "GSE264586", "GSE294109")
profile_status <- fread(file.path(profiles_root, "phase87_dataset_status_v2.tsv"))
profile_root_qc <- fread(file.path(profiles_root, "phase87_root_qc_v2.tsv"))
if (!setequal(profile_status$dataset, datasets) || any(profile_status$status != "PASS_TRAJECTORY_PROFILE_EXTRACTION_V2")) stop("Profile extraction gate failed")

root_rows <- list()
status_rows <- list()

for (dataset in datasets) {
  input_path <- file.path(profiles_root, paste0(dataset, "_phase87_profiles_v2.tsv.gz"))
  frame <- fread(input_path)
  pc_names <- paste0("PC_", seq_len(30L))
  required <- c("cell_id", "dataset", "patient_id", "sample_id", "phase8_state", "phase7_cluster", "root_cluster", "root_cell", pc_names)
  if (length(setdiff(required, names(frame)))) stop("Missing compact profile columns: ", dataset)
  if (any(frame$dataset != dataset) || sum(frame$root_cell) != 1L) stop("Compact profile identity failure: ", dataset)
  pca <- as.matrix(frame[, ..pc_names])
  if (any(!is.finite(pca))) stop("Non-finite compact PCA: ", dataset)
  cluster <- as.character(frame$phase7_cluster)
  root_cluster <- unique(as.character(frame$root_cluster))
  if (length(root_cluster) != 1L) stop("Non-unique root cluster: ", dataset)

  dummy <- sparseMatrix(i = integer(), j = integer(), dims = c(1L, nrow(frame)))
  rownames(dummy) <- "DUMMY"
  colnames(dummy) <- frame$cell_id
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
  pseudotime <- ifelse(denominator > 0, numerator / denominator, NA_real_)
  finite <- is.finite(pseudotime)
  if (sum(finite) < 0.8 * length(pseudotime)) stop("Less than 80% finite Slingshot pseudotime: ", dataset)
  observed <- range(pseudotime[finite])
  if (diff(observed) <= 0) stop("Degenerate Slingshot pseudotime: ", dataset)
  pseudotime[finite] <- (pseudotime[finite] - observed[[1L]]) / diff(observed)
  frame[, slingshot_pseudotime := pseudotime]
  fwrite(frame, file.path(part_dir, paste0(dataset, "_phase87_cells_v1.tsv.gz")), sep = "\t", quote = FALSE)

  base_root <- profile_root_qc[profile_root_qc$dataset == dataset]
  root_rows[[dataset]] <- data.table(
    dataset = dataset,
    root_cluster = root_cluster,
    root_cell = frame[root_cell == TRUE, cell_id][[1L]],
    trajectory_cells = nrow(frame),
    monocyte_cells = sum(frame$phase8_state == "Monocyte"),
    slingshot_lineages = ncol(pseudo),
    finite_slingshot_cells = sum(finite),
    root_cluster_total = base_root$root_cluster_total[[1L]],
    root_cluster_monocytes = base_root$root_cluster_monocytes[[1L]],
    root_cluster_monocyte_fraction = base_root$root_cluster_monocyte_fraction[[1L]]
  )
  source_status <- profile_status[profile_status$dataset == dataset]
  status_rows[[dataset]] <- data.table(
    dataset = dataset,
    source_path = source_status$source_path[[1L]],
    source_bytes = source_status$source_bytes[[1L]],
    source_cells = source_status$source_cells[[1L]],
    trajectory_cells = nrow(frame),
    patients = uniqueN(frame$patient_id),
    samples = uniqueN(frame$sample_id),
    clusters = uniqueN(frame$phase7_cluster),
    free_memory_gib_before = source_status$free_memory_gib_before[[1L]],
    status = "PASS_SLINGSHOT_EXTRACTION"
  )
  rm(frame, pca, dummy, sce, pseudo, weights)
  invisible(gc())
}

fwrite(rbindlist(root_rows), file.path(part_dir, "phase87_root_qc_v1.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(status_rows), file.path(part_dir, "phase87_dataset_status_v1.tsv"), sep = "\t", quote = FALSE)
file.copy(file.path(profiles_root, "phase87_program_coverage_v2.tsv"), file.path(part_dir, "phase87_program_coverage_v1.tsv"), overwrite = FALSE)
writeLines(capture.output(sessionInfo()), file.path(part_dir, "R_SESSION_INFO_v1.txt"))
writeLines("PASS_SLINGSHOT_EXTRACTION", file.path(part_dir, "PHASE87_SLINGSHOT_EXTRACTION_PASS_v1"))
if (!file.rename(part_dir, output_dir)) stop("Atomic Slingshot directory rename failed")
cat("PASS_SLINGSHOT_EXTRACTION\n")
