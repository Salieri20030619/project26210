options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 32 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else getwd(), winslash = "\\", mustWork = TRUE)
run_id <- Sys.getenv("PHASE6_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
base_dir <- file.path(root, "phase6_processing_local_v1")
run_dir <- file.path(base_dir, "runs", run_id)
object_dir <- file.path(run_dir, "objects_processed")
table_dir <- file.path(run_dir, "tables")
figure_dir <- file.path(run_dir, "figures")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
lock_file <- file.path(base_dir, "RUN_LOCK")
on.exit(if (file.exists(lock_file)) unlink(lock_file), add = TRUE)

datasets <- c("GSE202374", "GSE178481", "GSE264586", "GSE294109", "GSE207493")
crosswalk_file <- file.path(root, "phase6_metadata_confirmation_local_v1", "tables", "confirmed_object_patient_crosswalk.tsv")
crosswalk <- fread(crosswalk_file, sep = "\t", encoding = "UTF-8")
if (nrow(crosswalk) != 95L || any(crosswalk$mapping_status != "PASS")) stop("Confirmed metadata crosswalk gate failed")

append_progress <- function(dataset, event, details = "") {
  p <- file.path(run_dir, "PHASE6_PROCESSING_PROGRESS.tsv")
  z <- data.table(timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), dataset = dataset, event = event, details = details)
  fwrite(z, p, sep = "\t", quote = TRUE, append = file.exists(p), col.names = !file.exists(p))
}
write_once <- function(x, path) {
  if (file.exists(path)) {
    prior <- fread(path, sep = "\t", encoding = "UTF-8")
    same <- identical(names(prior), names(x)) && nrow(prior) == nrow(x) &&
      isTRUE(all.equal(as.data.frame(prior), as.data.frame(x),
                       check.attributes = FALSE, tolerance = 1e-10))
    if (!same) stop("Pre-existing table differs; refusing to overwrite: ", path)
    return(invisible(FALSE))
  }
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic table rename failed: ", path)
  invisible(TRUE)
}
save_rds_atomic <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  part <- paste0(path, ".", run_id, ".", format(Sys.time(), "%Y%m%d_%H%M%S"), ".part")
  saveRDS(x, part, compress = FALSE)
  if (!file.rename(part, path)) stop("Atomic RDS rename failed: ", path)
}

seed <- 26210L
set.seed(seed)
processed_this_session <- 0L

for (ds in datasets) {
  started <- Sys.time()
  completed_summary <- file.path(table_dir, paste0(ds, "_processing_summary.tsv"))
  completed_object <- file.path(object_dir, paste0(ds, "_processed_clustered.rds"))
  if (file.exists(completed_summary) && file.exists(completed_object)) {
    prior <- fread(completed_summary)
    if (nrow(prior) == 1L && prior$status == "PASS" && isTRUE(prior$output_readRDS_pass)) {
      append_progress(ds, "RESUME_SKIP_COMPLETED", sprintf("cells=%d bytes=%.0f", prior$cells_after, prior$output_bytes))
      message("PHASE6_DATASET_RESUME_SKIP=", ds)
      next
    }
    stop(ds, ": pre-existing outputs are not a validated completed pair")
  }
  append_progress(ds, "START")
  message("PHASE6_DATASET_START=", ds)
  input_file <- file.path(root, "phase5_merge_by_dataset_local_v2", "objects_raw", paste0(ds, "_merged_raw.rds"))
  input_before <- file.info(input_file)
  obj0 <- readRDS(input_file)
  if (!inherits(obj0, "Seurat")) stop(ds, ": input is not Seurat")
  if (!identical(Layers(obj0[[DefaultAssay(obj0)]]), "counts")) stop(ds, ": input must have exactly one counts layer")
  md <- obj0[[]]
  map <- crosswalk[dataset == ds]
  idx <- match(as.character(md$object_name), map$object_name)
  if (anyNA(idx)) stop(ds, ": metadata crosswalk does not cover all cells")

  md$phase5_patient_id_original <- as.character(md$patient_id)
  md$phase5_sample_id_original <- as.character(md$sample_id)
  md$phase5_treatment_original <- as.character(md$treatment)
  md$phase5_timepoint_original <- as.character(md$timepoint)
  md$patient_id <- map$patient_id_confirmed[idx]
  md$patient_key <- map$patient_key[idx]
  md$sample_id <- map$sample_id_confirmed[idx]
  md$timepoint <- map$timepoint_confirmed[idx]
  md$treatment <- map$treatment_confirmed[idx]
  md$response <- map$response_confirmed[idx]
  md$tissue_confirmed <- map$tissue_confirmed[idx]
  md$sex_confirmed <- map$sex_confirmed[idx]
  md$age_confirmed <- map$age_confirmed[idx]
  md$metadata_mapping_confidence <- map$mapping_confidence[idx]

  counts <- LayerData(obj0, assay = DefaultAssay(obj0), layer = "counts")
  mt <- grep("^MT-", rownames(counts))
  pct_mt <- if (length(mt)) 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(1, Matrix::colSums(counts)) else rep(0, ncol(counts))
  md$percent.mt <- pct_mt
  pass <- md$nFeature_RNA >= 200 & md$nFeature_RNA <= 8000 & md$nCount_RNA >= 500 & md$nCount_RNA <= 100000 & md$percent.mt <= 25
  md$phase6_qc_pass <- pass
  qc_reason <- rep("PASS", length(pass))
  qc_reason[md$nFeature_RNA < 200] <- "nFeature_below_200"
  qc_reason[md$nFeature_RNA > 8000] <- "nFeature_above_8000"
  qc_reason[md$nCount_RNA < 500] <- "nCount_below_500"
  qc_reason[md$nCount_RNA > 100000] <- "nCount_above_100000"
  qc_reason[md$percent.mt > 25] <- "percent_mt_above_25"
  md$phase6_qc_primary_reason <- qc_reason

  qc_by_object <- data.table(object_name = as.character(md$object_name), patient_id = md$patient_id, pass = pass, reason = qc_reason)[,
    .(cells_before = .N, cells_after = sum(pass), cells_removed = sum(!pass), pass_fraction = mean(pass)),
    by = .(object_name, patient_id)]
  qc_by_object[, dataset := ds]
  write_once(qc_by_object, file.path(table_dir, paste0(ds, "_qc_reconciliation.tsv")))
  if (sum(pass) < 1000L) stop(ds, ": fewer than 1000 cells pass QC")

  counts_qc <- counts[, pass, drop = FALSE]
  md_qc <- md[pass, , drop = FALSE]
  rm(obj0, counts, md, pct_mt, pass, qc_reason); invisible(gc())
  obj <- CreateSeuratObject(counts = counts_qc, assay = "RNA", project = ds, meta.data = md_qc, min.cells = 0, min.features = 0)
  rm(counts_qc, md_qc); invisible(gc())

  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
  variable <- VariableFeatures(obj)
  obj <- ScaleData(obj, features = variable, block.size = 500, min.cells.to.block = 0, verbose = FALSE)
  obj <- RunPCA(obj, features = variable, npcs = 30, approx = TRUE, seed.use = seed, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, k.param = 20, nn.method = "annoy", annoy.metric = "euclidean", prune.SNN = 1/15, verbose = FALSE)
  obj <- FindClusters(obj, resolution = c(0.4, 0.8), algorithm = 1, random.seed = seed, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, umap.method = "uwot", metric = "cosine", n.neighbors = 30, min.dist = 0.3, seed.use = seed, verbose = FALSE)

  cluster_cols <- grep("_snn_res\\.(0\\.4|0\\.8)$", colnames(obj[[]]), value = TRUE)
  if (length(cluster_cols) != 2L) stop(ds, ": expected two cluster columns, got ", paste(cluster_cols, collapse = ","))
  primary_col <- grep("0\\.8$", cluster_cols, value = TRUE)
  Idents(obj) <- obj[[primary_col, drop = TRUE]]

  vf <- data.table(dataset = ds, rank = seq_along(variable), gene = variable)
  write_once(vf, file.path(table_dir, paste0(ds, "_variable_features.tsv")))
  pca_sd <- Stdev(obj, reduction = "pca")
  write_once(data.table(dataset = ds, PC = seq_along(pca_sd), stdev = pca_sd), file.path(table_dir, paste0(ds, "_pca_stdev.tsv")))
  cluster_counts <- rbindlist(lapply(cluster_cols, function(cc) {
    z <- as.character(obj[[cc, drop = TRUE]])
    data.table(dataset = ds, resolution = sub(".*res\\.", "", cc), cluster = z)[, .(cells = .N), by = .(dataset, resolution, cluster)]
  }))
  write_once(cluster_counts, file.path(table_dir, paste0(ds, "_cluster_cell_counts.tsv")))
  patient_counts <- data.table(patient_id = obj$patient_id, timepoint = obj$timepoint)[, .(cells = .N), by = .(patient_id, timepoint)][, dataset := ds]
  write_once(patient_counts, file.path(table_dir, paste0(ds, "_patient_cell_counts.tsv")))

  p1 <- DimPlot(obj, reduction = "umap", group.by = primary_col, label = TRUE, repel = TRUE, raster = TRUE) + ggtitle(paste(ds, "UMAP clusters, resolution 0.8"))
  ggsave(file.path(figure_dir, paste0(ds, "_umap_clusters_res0.8.png")), p1, width = 10, height = 8, dpi = 180)
  p2 <- DimPlot(obj, reduction = "umap", group.by = "patient_id", raster = TRUE) + ggtitle(paste(ds, "UMAP by confirmed patient"))
  ggsave(file.path(figure_dir, paste0(ds, "_umap_patient.png")), p2, width = 10, height = 8, dpi = 180)
  rm(p1, p2); invisible(gc())

  cells_after <- ncol(obj); features_after <- nrow(obj)
  clusters04 <- uniqueN(obj[[grep("0\\.4$", cluster_cols, value = TRUE), drop = TRUE]])
  clusters08 <- uniqueN(obj[[primary_col, drop = TRUE]])
  obj@misc$phase6_processing <- list(run_id = run_id, seed = seed, qc = c(nFeature_min = 200, nFeature_max = 8000, nCount_min = 500, nCount_max = 100000, percent_mt_max = 25), variable_features = 1500, pcs = 30, cluster_resolutions = c(0.4, 0.8))
  obj <- DietSeurat(obj, assays = "RNA", layers = c("counts", "data"), dimreducs = c("pca", "umap"), graphs = names(obj@graphs), misc = TRUE)

  output_file <- file.path(object_dir, paste0(ds, "_processed_clustered.rds"))
  save_rds_atomic(obj, output_file)
  output_bytes <- file.info(output_file)$size
  rm(obj); invisible(gc())

  check <- readRDS(output_file)
  check_ok <- inherits(check, "Seurat") && ncol(check) == cells_after && nrow(check) == features_after && all(c("pca", "umap") %in% names(check@reductions)) && all(cluster_cols %in% colnames(check[[]])) && all(c("counts", "data") %in% Layers(check[["RNA"]]))
  if (!check_ok) stop(ds, ": output reload validation failed")
  rm(check); invisible(gc())
  input_after <- file.info(input_file)
  input_unchanged <- identical(input_before$size, input_after$size) && identical(as.numeric(input_before$mtime), as.numeric(input_after$mtime))
  if (!input_unchanged) stop(ds, ": input RDS changed")

  summary <- data.table(
    dataset = ds,
    cells_before = sum(qc_by_object$cells_before),
    cells_after = cells_after,
    cells_removed = sum(qc_by_object$cells_removed),
    pass_fraction = cells_after / sum(qc_by_object$cells_before),
    features = features_after,
    variable_features = length(variable),
    pca_components = 30L,
    clusters_res0_4 = clusters04,
    clusters_res0_8 = clusters08,
    output_file = normalizePath(output_file, winslash = "\\", mustWork = TRUE),
    output_bytes = output_bytes,
    output_readRDS_pass = check_ok,
    source_rds_modified = !input_unchanged,
    elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
    status = "PASS"
  )
  write_once(summary, file.path(table_dir, paste0(ds, "_processing_summary.tsv")))
  append_progress(ds, "COMPLETE", sprintf("cells=%d clusters_res0.8=%d bytes=%.0f", cells_after, clusters08, output_bytes))
  message("PHASE6_DATASET_COMPLETE=", ds)
  invisible(gc())

  # A fresh R process per large dataset prevents unreleased/fragmented heap memory
  # from accumulating across datasets on this 32 GiB Windows workstation.
  processed_this_session <- processed_this_session + 1L
  if (processed_this_session >= 1L && ds != tail(datasets, 1L)) {
    append_progress(ds, "PROCESS_RECYCLE", "completed dataset; restart R before next dataset")
    if (file.exists(lock_file)) unlink(lock_file)
    message("PHASE6_PROCESS_RECYCLE_AFTER=", ds)
    quit(save = "no", status = 75L, runLast = FALSE)
  }
}

summary_files <- file.path(table_dir, paste0(datasets, "_processing_summary.tsv"))
all_summary <- rbindlist(lapply(summary_files, fread))
write_once(all_summary, file.path(table_dir, "processing_summary.tsv"))
all_clusters <- rbindlist(lapply(datasets, function(ds) fread(file.path(table_dir, paste0(ds, "_cluster_cell_counts.tsv")))))
write_once(all_clusters, file.path(table_dir, "cluster_cell_counts.tsv"))
status <- c(
  "status=PASS",
  paste0("run_id=", run_id),
  "datasets_expected=5",
  paste0("datasets_completed=", nrow(all_summary)),
  paste0("cells_input=", sum(all_summary$cells_before)),
  paste0("cells_after_qc=", sum(all_summary$cells_after)),
  paste0("all_outputs_readRDS=", all(all_summary$output_readRDS_pass)),
  paste0("source_rds_modified=", sum(all_summary$source_rds_modified)),
  "normalization_performed=1",
  "pca_performed=1",
  "umap_performed=1",
  "clustering_performed=1",
  "integration_performed=0",
  "annotation_performed=0",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(status, file.path(run_dir, "PHASE6_PROCESSING_STATUS"), useBytes = TRUE)
append_progress("ALL", "FINAL", "status=PASS")
message("PHASE6_PROCESSING=PASS")
