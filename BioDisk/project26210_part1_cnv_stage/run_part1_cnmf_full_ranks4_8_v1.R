#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(RcppML)
  library(data.table)
  library(cluster)
  library(jsonlite)
  library(digest)
})

PROJECT <- "/mnt/BioDisk/ludingyu/project26210"
ROOT <- file.path(PROJECT, "results/part1_rebuild_v1")
INPUT <- file.path(ROOT, "05_cnmf/input_v1")
OUT <- file.path(ROOT, "05_cnmf/full_ranks4_8_v1")
MEMORY_GATE <- file.path(PROJECT, "check_memory_gate.sh")
RANKS <- 4:8
N_RUNS <- 20L
MAXIT <- 200L
SEED_BASE <- 26210L
N_THREADS <- 4L

if (dir.exists(OUT) && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty full cNMF directory: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

gate_records <- list()
gate <- function(stage) {
  output <- system2(MEMORY_GATE, stage, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  text <- paste(output, collapse = " ")
  message(text)
  gate_records[[length(gate_records) + 1L]] <<- data.table(
    stage = stage, exit_status = status, gate_output = text
  )
  fwrite(rbindlist(gate_records), file.path(OUT, "memory_gate_log_v1.tsv"), sep = "\t")
  if (status != 0L || !grepl("PASS", text, fixed = TRUE)) {
    stop("Memory gate blocked stage ", stage, ": ", text)
  }
}

gate("full_load_input_v1")
for (terminal in c(
  file.path(INPUT, "RUN_COMPLETE"),
  file.path(ROOT, "05_cnmf/input_verification_v1/RUN_COMPLETE"),
  file.path(ROOT, "05_cnmf/pilot_verification_v1/RUN_COMPLETE")
)) if (!file.exists(terminal)) stop("Required terminal marker missing: ", terminal)

input_manifest <- fread(file.path(INPUT, "file_manifest_v1.tsv"), colClasses = "character")
matrix_name <- "cnmf_hvg_scaled_genes_by_cells_v1.mtx.gz"
matrix_row <- input_manifest[file == matrix_name]
if (nrow(matrix_row) != 1L) stop("Matrix manifest entry is absent or duplicated")
matrix_path <- file.path(INPUT, matrix_name)
matrix_sha256 <- digest(matrix_path, algo = "sha256", file = TRUE)
if (!identical(matrix_sha256, matrix_row$sha256)) stop("Input matrix hash mismatch")
cells <- fread(file.path(INPUT, "cnmf_cells_v1.tsv.gz"), colClasses = c(exact_cell_id = "character"))
genes <- fread(file.path(INPUT, "cnmf_hvg_genes_v1.tsv"), colClasses = c(gene = "character"))
if (nrow(cells) != 11950L || anyDuplicated(cells$exact_cell_id)) stop("Invalid cells")
if (nrow(genes) != 2500L || anyDuplicated(genes$gene)) stop("Invalid genes")

connection <- gzfile(matrix_path, open = "rb")
A <- readMM(connection)
close(connection)
A <- as(A, "CsparseMatrix")
if (!inherits(A, "dgCMatrix")) A <- as(A, "dgCMatrix")
dimnames(A) <- list(genes$gene, cells$exact_cell_id)
if (!identical(dim(A), c(2500L, 11950L)) || length(A@x) != 2416110L ||
    any(!is.finite(A@x)) || any(A@x < 0) || max(A@x) > 10 + 1e-12) {
  stop("Invalid cNMF matrix")
}
RcppML::setRcppMLthreads(N_THREADS)
if (RcppML::getRcppMLthreads() != N_THREADS) stop("Failed to set RcppML threads")

consensus_for_rank <- function(fits, rank_value) {
  loading_list <- list()
  usage_list <- list()
  meta_list <- list()
  counter <- 0L
  for (run_name in names(fits)) {
    fit <- fits[[run_name]]
    loading <- sweep(fit$w, 2L, fit$d, "*")
    usage <- sweep(fit$h, 1L, fit$d, "*")
    loading_sum <- colSums(loading)
    usage_sum <- colSums(usage)
    if (any(!is.finite(loading_sum)) || any(loading_sum <= 0) ||
        any(!is.finite(usage_sum)) || any(usage_sum <= 0)) stop("Invalid factor sums")
    loading <- sweep(loading, 2L, loading_sum, "/")
    usage <- sweep(usage, 2L, usage_sum, "/")
    for (component in seq_len(rank_value)) {
      counter <- counter + 1L
      component_id <- paste0("run", run_name, "_component", component)
      loading_list[[counter]] <- loading[, component]
      usage_list[[counter]] <- usage[component, ]
      meta_list[[counter]] <- data.table(
        component_id = component_id, run = as.integer(run_name), component = component
      )
    }
  }
  loading_matrix <- do.call(cbind, loading_list)
  usage_matrix <- do.call(cbind, usage_list)
  component_info <- rbindlist(meta_list)
  colnames(loading_matrix) <- component_info$component_id
  colnames(usage_matrix) <- component_info$component_id
  loading_l2 <- sqrt(colSums(loading_matrix^2))
  if (any(!is.finite(loading_l2)) || any(loading_l2 <= 0)) stop("Invalid loading norms")
  loading_unit <- sweep(loading_matrix, 2L, loading_l2, "/")
  set.seed(SEED_BASE + rank_value * 1000L)
  km <- kmeans(t(loading_unit), centers = rank_value, nstart = 100L, iter.max = 300L)
  component_info[, consensus_cluster := as.integer(km$cluster)]
  sil <- cluster::silhouette(km$cluster, stats::dist(t(loading_unit)))
  component_info[, silhouette := as.numeric(sil[, "sil_width"])]

  factor_names <- paste0("F", seq_len(rank_value))
  consensus_loading <- matrix(0, nrow = nrow(A), ncol = rank_value,
                              dimnames = list(rownames(A), factor_names))
  consensus_usage <- matrix(0, nrow = ncol(A), ncol = rank_value,
                            dimnames = list(colnames(A), factor_names))
  for (cl in seq_len(rank_value)) {
    idx <- which(component_info$consensus_cluster == cl)
    runs_present <- sort(unique(component_info$run[idx]))
    loading_by_run <- lapply(runs_present, function(run_id) {
      ii <- idx[component_info$run[idx] == run_id]
      rowMeans(loading_matrix[, ii, drop = FALSE])
    })
    usage_by_run <- lapply(runs_present, function(run_id) {
      ii <- idx[component_info$run[idx] == run_id]
      rowMeans(usage_matrix[, ii, drop = FALSE])
    })
    consensus_loading[, cl] <- rowMeans(do.call(cbind, loading_by_run))
    consensus_usage[, cl] <- rowMeans(do.call(cbind, usage_by_run))
  }
  usage_row_sum <- rowSums(consensus_usage)
  loading_col_sum <- colSums(consensus_loading)
  if (any(!is.finite(usage_row_sum)) || any(usage_row_sum <= 0) ||
      any(!is.finite(loading_col_sum)) || any(loading_col_sum <= 0)) stop("Invalid consensus")
  consensus_usage <- consensus_usage / usage_row_sum
  consensus_loading <- sweep(consensus_loading, 2L, loading_col_sum, "/")
  cluster_qc <- component_info[, .(
    n_components = .N,
    n_unique_runs = uniqueN(run),
    mean_silhouette = mean(silhouette),
    min_silhouette = min(silhouette)
  ), by = consensus_cluster]
  list(
    loading = consensus_loading,
    usage = consensus_usage,
    component_info = component_info,
    cluster_qc = cluster_qc,
    mean_silhouette = mean(component_info$silhouette),
    min_silhouette = min(component_info$silhouette),
    min_run_fraction = min(cluster_qc$n_unique_runs) / N_RUNS
  )
}

fit_store <- list()
run_qc <- list()
rank_results <- list()
for (rank_value in RANKS) {
  rank_fits <- list()
  for (run_id in seq_len(N_RUNS)) {
    gate(paste0("full_rank", rank_value, "_run", run_id))
    seed <- SEED_BASE + rank_value * 1000L + run_id
    start <- proc.time()[["elapsed"]]
    fit <- RcppML::nmf(A, k = rank_value, maxit = MAXIT, seed = seed, verbose = FALSE)
    elapsed <- proc.time()[["elapsed"]] - start
    valid <- identical(dim(fit$w), c(nrow(A), rank_value)) &&
      identical(dim(fit$h), c(rank_value, ncol(A))) && length(fit$d) == rank_value &&
      all(is.finite(fit$w)) && all(fit$w >= 0) && all(is.finite(fit$h)) &&
      all(fit$h >= 0) && all(is.finite(fit$d)) && all(fit$d >= 0)
    if (!valid) stop("Invalid fit at rank ", rank_value, " run ", run_id)
    mse <- RcppML::mse(A, w = fit$w, d = fit$d, h = fit$h)
    if (length(mse) != 1L || !is.finite(mse) || mse < 0) stop("Invalid fit MSE")
    rank_fits[[as.character(run_id)]] <- fit
    run_qc[[length(run_qc) + 1L]] <- data.table(
      rank = rank_value, run = run_id, seed = seed, iter = fit$iter,
      tolerance = fit$tol, mse = mse, elapsed_seconds = elapsed,
      fit_finite_nonnegative = valid
    )
    fwrite(rbindlist(run_qc), file.path(OUT, "run_qc_v1.tsv"), sep = "\t")
  }
  fit_store[[as.character(rank_value)]] <- rank_fits
  gate(paste0("full_rank", rank_value, "_consensus"))
  rank_results[[as.character(rank_value)]] <- consensus_for_rank(rank_fits, rank_value)
}

run_qc_table <- rbindlist(run_qc)
rank_qc <- rbindlist(lapply(names(rank_results), function(rank_name) {
  rank_value <- as.integer(rank_name)
  result <- rank_results[[rank_name]]
  rank_run <- run_qc_table[rank == rank_value]
  data.table(
    rank = rank_value,
    n_runs = N_RUNS,
    mean_component_silhouette = result$mean_silhouette,
    min_component_silhouette = result$min_silhouette,
    min_run_fraction = result$min_run_fraction,
    min_cluster_unique_runs = min(result$cluster_qc$n_unique_runs),
    max_cluster_unique_runs = max(result$cluster_qc$n_unique_runs),
    min_cluster_components = min(result$cluster_qc$n_components),
    max_cluster_components = max(result$cluster_qc$n_components),
    median_run_mse = median(rank_run$mse),
    iqr_run_mse = IQR(rank_run$mse),
    median_elapsed_seconds = median(rank_run$elapsed_seconds),
    stability_gate = result$mean_silhouette >= 0.70 && result$min_run_fraction >= 0.80
  )
}))
fwrite(rank_qc, file.path(OUT, "rank_qc_v1.tsv"), sep = "\t")
component_qc <- rbindlist(lapply(names(rank_results), function(rank_name) {
  x <- copy(rank_results[[rank_name]]$component_info)
  x[, rank := as.integer(rank_name)]
  setcolorder(x, c("rank", "component_id", "run", "component", "consensus_cluster", "silhouette"))
  x
}))
fwrite(component_qc, file.path(OUT, "component_assignments_v1.tsv"), sep = "\t")

for (rank_name in names(rank_results)) {
  rank_value <- as.integer(rank_name)
  result <- rank_results[[rank_name]]
  fwrite(as.data.table(result$usage, keep.rownames = "exact_cell_id"),
         file.path(OUT, paste0("consensus_usage_rank", rank_value, "_v1.tsv.gz")),
         sep = "\t", compress = "gzip")
  fwrite(as.data.table(result$loading, keep.rownames = "gene"),
         file.path(OUT, paste0("consensus_loading_rank", rank_value, "_v1.tsv.gz")),
         sep = "\t", compress = "gzip")
}

saveRDS(
  list(fits = fit_store, consensus = rank_results, ranks = RANKS, n_runs = N_RUNS,
       maxit = MAXIT, seed_base = SEED_BASE, input_matrix_sha256 = matrix_sha256),
  file.path(OUT, "full_models_v1.rds"), compress = "xz"
)
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info_v1.txt"))
metadata <- list(
  schema_version = "part1_cnmf_full_ranks4_8_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  purpose = "formal multi-start consensus for ranks 4-8; biological rank not yet selected",
  input_matrix = matrix_path,
  input_matrix_sha256 = matrix_sha256,
  matrix_orientation = "genes_by_cells",
  n_genes = nrow(A), n_cells = ncol(A), ranks = RANKS,
  n_runs_per_rank = N_RUNS, maxit = MAXIT, seed_base = SEED_BASE,
  rcppml_threads = N_THREADS, r_version = R.version.string,
  RcppML_version = as.character(packageVersion("RcppML")),
  consensus = "k-means on unit-L2 component loadings; per-run mean within cluster; equal run weighting",
  biological_rank_selected = FALSE,
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
print(rank_qc)
cat("PART1_CNMF_FULL_RANKS4_8_V1_PASS\n")
