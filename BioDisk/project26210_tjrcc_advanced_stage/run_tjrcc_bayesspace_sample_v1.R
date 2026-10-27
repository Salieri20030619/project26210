#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(BayesSpace)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(rhdf5)
  library(ggplot2)
  library(jsonlite)
})

project_root <- "/mnt/BioDisk/ludingyu/project26210"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: run_tjrcc_bayesspace_sample_v1.R SAMPLE INPUT_H5AD OUTPUT_DIR SEED")
}
sample_id <- args[[1]]
input_h5ad <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
seed <- as.integer(args[[4]])
if (!grepl("^[A-Za-z0-9_]+$", sample_id) || is.na(seed)) {
  stop("Invalid sample ID or seed")
}
gate_script <- file.path(project_root, "check_memory_gate.sh")
q_candidates <- 3:8
q_tune_nrep <- 2000L
q_tune_burn <- 500L
cluster_nrep <- 50000L
cluster_burn <- 5000L
cluster_thin <- 100L
completed <- FALSE

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_gate <- function(stage) {
  cat(sprintf("\n[%s] MEMORY_GATE stage=%s\n", format(Sys.time(), "%F %T %z"), stage))
  status <- system2(gate_script)
  if (!identical(status, 0L)) {
    stop("Memory gate failed at stage: ", stage)
  }
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

write_status <- function(status, note = NULL) {
  payload <- list(
    status = status,
    sample = sample_id,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    input_h5ad = input_h5ad,
    output_dir = output_dir,
    note = note
  )
  write_json(payload, file.path(output_dir, "status.json"), auto_unbox = TRUE, pretty = TRUE)
}

on.exit({
  if (!completed) {
    write_status("failed", "R BayesSpace stage exited before completion; inspect run.log")
  }
}, add = TRUE)

write_status("running", "BayesSpace preprocessing and spatial clustering")
set.seed(seed)
run_gate("read_h5ad_csr")

attrs <- h5readAttributes(input_h5ad, "/X")
if (!identical(attrs[["encoding-type"]], "csr_matrix")) {
  stop("Expected /X to be csr_matrix")
}
shape <- as.integer(attrs$shape)
if (length(shape) != 2L || shape[[1]] < 100L || shape[[2]] < 2000L) {
  stop("Unexpected H5AD shape: ", paste(shape, collapse = " x "))
}

x_data <- as.numeric(h5read(input_h5ad, "/X/data"))
x_indices <- as.integer(h5read(input_h5ad, "/X/indices"))
x_indptr <- as.integer(h5read(input_h5ad, "/X/indptr"))
barcodes <- as.character(h5read(input_h5ad, "/obs/_index"))
genes <- as.character(h5read(input_h5ad, "/var/_index"))
array_row <- as.integer(h5read(input_h5ad, "/obs/array_row"))
array_col <- as.integer(h5read(input_h5ad, "/obs/array_col"))
spatial_xy <- t(h5read(input_h5ad, "/obsm/spatial"))

if (length(x_data) < shape[[1]] || any(!is.finite(x_data)) || any(x_data < 0) ||
    any(abs(x_data - round(x_data)) > 1e-8)) {
  stop(sample_id, " count matrix failed integer/nonnegative/nnz validation")
}
if (length(unique(barcodes)) != shape[[1]] || length(unique(genes)) != shape[[2]]) {
  stop("Barcode or gene identifiers are not unique")
}
if (!identical(dim(spatial_xy), c(shape[[1]], 2L))) {
  stop("Unexpected spatial-coordinate dimensions")
}

counts_spot_gene <- new(
  "dgRMatrix",
  p = x_indptr,
  j = x_indices,
  x = x_data,
  Dim = shape
)
validObject(counts_spot_gene)
counts_gene_spot <- t(counts_spot_gene)
if (!inherits(counts_gene_spot, "dgCMatrix")) {
  counts_gene_spot <- as(counts_gene_spot, "CsparseMatrix")
}
rownames(counts_gene_spot) <- genes
colnames(counts_gene_spot) <- barcodes

col_data <- S4Vectors::DataFrame(
  array_row = array_row,
  array_col = array_col,
  pxl_row_in_fullres = spatial_xy[, 2],
  pxl_col_in_fullres = spatial_xy[, 1],
  row.names = barcodes
)
sce <- SingleCellExperiment(
  assays = list(counts = counts_gene_spot),
  colData = col_data
)
metadata(sce)$source_h5ad <- input_h5ad
metadata(sce)$source_h5ad_sha256 <- sha256_file(input_h5ad)
metadata(sce)$seed <- seed

observed_sums <- as.numeric(Matrix::colSums(counts(sce)))
expected_sums <- as.numeric(h5read(input_h5ad, "/obs/total_counts"))
if (length(expected_sums) != length(observed_sums) || any(!is.finite(expected_sums)) ||
    any(expected_sums < observed_sums)) {
  stop("H5AD obs.total_counts is invalid or smaller than shared-gene X library size")
}
total_counts_audit <- list(
  x_library_median = median(observed_sums),
  obs_total_counts_median = median(expected_sums),
  difference_min = min(expected_sums - observed_sums),
  difference_median = median(expected_sums - observed_sums),
  difference_max = max(expected_sums - observed_sums),
  interpretation = "X contains the 11702 genes shared with the cell2location reference; obs.total_counts retains the pre-intersection full-gene library size. BayesSpace uses X counts."
)

rm(counts_spot_gene, x_data, x_indices, x_indptr, observed_sums, expected_sums)
invisible(gc())

run_gate("spatial_preprocess")
sce <- spatialPreprocess(
  sce,
  platform = "Visium",
  n.PCs = 15,
  n.HVGs = 2000,
  skip.PCA = FALSE,
  log.normalize = TRUE,
  assay.type = "logcounts",
  BSPARAM = BiocSingular::ExactParam(),
  BPPARAM = BiocParallel::SerialParam()
)
if (!"PCA" %in% reducedDimNames(sce) || ncol(reducedDim(sce, "PCA")) != 15L) {
  stop("BayesSpace preprocessing did not generate 15 PCs")
}
saveRDS(sce, file.path(output_dir, paste0(sample_id, "_bayesspace_preprocessed.rds")), compress = FALSE)

run_gate("q_tune")
set.seed(seed)
sce <- qTune(
  sce,
  qs = q_candidates,
  burn.in = q_tune_burn,
  nrep = q_tune_nrep,
  cores = 1L,
  use.dimred = "PCA",
  d = 15,
  platform = "Visium",
  init.method = "mclust",
  model = "t",
  precision = "equal"
)
q_table <- attr(sce, "q.logliks")
if (!identical(as.integer(q_table$q), q_candidates) || any(!is.finite(q_table$loglik))) {
  stop("qTune returned incomplete or non-finite results")
}

break_candidates <- seq.int(2L, nrow(q_table) - 1L)
piecewise_sse <- vapply(break_candidates, function(index) {
  left <- q_table[seq_len(index), , drop = FALSE]
  right <- q_table[index:nrow(q_table), , drop = FALSE]
  sum(residuals(lm(loglik ~ q, data = left))^2) +
    sum(residuals(lm(loglik ~ q, data = right))^2)
}, numeric(1))
selected_index <- break_candidates[which.min(piecewise_sse)]
selected_q <- as.integer(q_table$q[[selected_index]])
q_table$selected <- q_table$q == selected_q
q_table$piecewise_sse_at_break <- NA_real_
q_table$piecewise_sse_at_break[break_candidates] <- piecewise_sse
write.table(
  q_table,
  file.path(output_dir, paste0(sample_id, "_q_tuning.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

q_plot <- ggplot(q_table, aes(x = q, y = -loglik)) +
  geom_line(linewidth = 0.7, color = "#355C7D") +
  geom_point(size = 2.5, color = "#355C7D") +
  geom_vline(xintercept = selected_q, linetype = 2, color = "#C0392B") +
  annotate("text", x = selected_q + 0.15, y = max(-q_table$loglik),
           label = paste0("Selected q=", selected_q), hjust = 0, vjust = 1,
           color = "#C0392B", size = 4) +
  labs(x = "Number of domains (q)", y = "Negative pseudo-log-likelihood",
       title = paste0(sample_id, " BayesSpace q tuning")) +
  theme_bw(base_size = 12)
ggsave(file.path(output_dir, paste0(sample_id, "_q_tuning.png")), q_plot, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, paste0(sample_id, "_q_tuning.pdf")), q_plot, width = 7, height = 5)

run_gate("spatial_cluster")
set.seed(seed)
sce <- spatialCluster(
  sce,
  q = selected_q,
  use.dimred = "PCA",
  d = 15,
  platform = "Visium",
  init.method = "mclust",
  model = "t",
  precision = "equal",
  nrep = cluster_nrep,
  burn.in = cluster_burn,
  thin = cluster_thin,
  save.chain = FALSE
)

cluster_labels <- as.integer(colData(sce)$spatial.cluster)
if (length(cluster_labels) != ncol(sce) || anyNA(cluster_labels) ||
    length(unique(cluster_labels)) != selected_q) {
  stop("Spatial clustering output failed label validation")
}

assignment <- data.frame(
  barcode = colnames(sce),
  array_row = as.integer(sce$array_row),
  array_col = as.integer(sce$array_col),
  spatial_x_fullres = as.numeric(sce$pxl_col_in_fullres),
  spatial_y_fullres = as.numeric(sce$pxl_row_in_fullres),
  bayesspace_cluster = cluster_labels
)
write.table(
  assignment,
  file.path(output_dir, paste0(sample_id, "_bayesspace_assignments.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cluster_sizes <- as.data.frame(table(bayesspace_cluster = cluster_labels))
cluster_sizes$proportion <- cluster_sizes$Freq / sum(cluster_sizes$Freq)
write.table(
  cluster_sizes,
  file.path(output_dir, paste0(sample_id, "_bayesspace_cluster_sizes.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

assignment$cluster_factor <- factor(assignment$bayesspace_cluster)
cluster_plot <- ggplot(assignment, aes(x = array_col, y = -array_row, color = cluster_factor)) +
  geom_point(size = 1.8) +
  coord_equal() +
  scale_color_brewer(palette = "Dark2", name = "Domain") +
  labs(x = "Array column", y = "Array row", title = paste0(sample_id, " BayesSpace domains")) +
  theme_void(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
ggsave(file.path(output_dir, paste0(sample_id, "_bayesspace_domains.png")), cluster_plot,
       width = 7, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, paste0(sample_id, "_bayesspace_domains.pdf")), cluster_plot,
       width = 7, height = 6.5, bg = "white")

metadata(sce)$BayesSpace.pilot <- list(
  selected_q = selected_q,
  q_candidates = q_candidates,
  q_selection_method = "minimum two-segment linear SSE on qTune pseudo-log-likelihood",
  q_tune_nrep = q_tune_nrep,
  q_tune_burn = q_tune_burn,
  cluster_nrep = cluster_nrep,
  cluster_burn = cluster_burn,
  cluster_thin = cluster_thin,
  seed = seed
)
saveRDS(sce, file.path(output_dir, paste0(sample_id, "_bayesspace_clustered.rds")), compress = FALSE)

run_metadata <- list(
  status = "bayesspace_complete",
  sample = sample_id,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_h5ad = input_h5ad,
  input_sha256 = metadata(sce)$source_h5ad_sha256,
  n_spots = ncol(sce),
  n_genes = nrow(sce),
  nnz = length(counts(sce)@x),
  total_counts_audit = total_counts_audit,
  selected_q = selected_q,
  cluster_sizes = as.list(setNames(cluster_sizes$Freq, cluster_sizes$bayesspace_cluster)),
  seed = seed,
  packages = list(
    R = R.version.string,
    BayesSpace = as.character(packageVersion("BayesSpace")),
    SingleCellExperiment = as.character(packageVersion("SingleCellExperiment")),
    rhdf5 = as.character(packageVersion("rhdf5"))
  ),
  evidence_boundary = "BayesSpace domains are computational domains, not pathologist-annotated tumor boundaries"
)
write_json(run_metadata, file.path(output_dir, "bayesspace_run_metadata.json"),
           auto_unbox = TRUE, pretty = TRUE)

completed <- TRUE
write_status("bayesspace_complete", paste0("selected_q=", selected_q, "; boundary postprocessing pending"))
cat(sprintf("\nBAYESSPACE_TJRCC_SAMPLE_PASS sample=%s selected_q=%d\n", sample_id, selected_q))
