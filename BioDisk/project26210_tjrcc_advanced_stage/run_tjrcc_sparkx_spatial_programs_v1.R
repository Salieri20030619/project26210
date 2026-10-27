#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
  library(SPARK)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: run_tjrcc_sparkx_spatial_programs_v1.R SAMPLE INPUT_H5AD OUTPUT_DIR SEED")
}

sample_id <- args[[1]]
input_h5ad <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
seed <- as.integer(args[[4]])
if (!grepl("^[A-Za-z0-9_]+$", sample_id) || is.na(seed)) {
  stop("Invalid sample ID or seed")
}
project_root <- "/mnt/BioDisk/ludingyu/project26210"
gate_script <- file.path(project_root, "check_memory_gate.sh")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completed <- FALSE

run_gate <- function(stage) {
  cat(sprintf("[%s] MEMORY_GATE stage=%s\n", format(Sys.time(), "%F %T %z"), stage))
  status <- system2(gate_script)
  if (!identical(status, 0L)) {
    stop("Memory gate failed at stage: ", stage)
  }
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

write_stage_status <- function(status, note = NULL) {
  payload <- list(
    status = status,
    sample = sample_id,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    input_h5ad = input_h5ad,
    output_dir = output_dir,
    note = note
  )
  write_json(payload, file.path(output_dir, "sparkx_status.json"), auto_unbox = TRUE, pretty = TRUE)
}

on.exit({
  if (!completed) {
    write_stage_status("failed", "SPARK-X stage exited before completion; inspect run.log")
  }
}, add = TRUE)

write_stage_status("running", "SPARK-X spatial variability testing")
set.seed(seed)
run_gate("sparkx_read_h5ad")

attrs <- h5readAttributes(input_h5ad, "/X")
if (!identical(attrs[["encoding-type"]], "csr_matrix")) {
  stop("Expected /X to be csr_matrix")
}
shape <- as.integer(attrs$shape)
if (length(shape) != 2L || shape[[1]] < 100L || shape[[2]] < 2000L) {
  stop("Unexpected H5AD shape: ", paste(shape, collapse = " x "))
}

x_data <- as.numeric(h5read(input_h5ad, "/X/data"))
nnz_input <- length(x_data)
x_indices <- as.integer(h5read(input_h5ad, "/X/indices"))
x_indptr <- as.integer(h5read(input_h5ad, "/X/indptr"))
barcodes <- as.character(h5read(input_h5ad, "/obs/_index"))
genes <- as.character(h5read(input_h5ad, "/var/_index"))
spatial_xy <- t(h5read(input_h5ad, "/obsm/spatial"))

if (length(x_data) < shape[[1]] || any(!is.finite(x_data)) || any(x_data < 0) ||
    any(abs(x_data - round(x_data)) > 1e-8)) {
  stop("Count matrix failed integer/nonnegative/nnz validation")
}
if (length(unique(barcodes)) != shape[[1]] || length(unique(genes)) != shape[[2]]) {
  stop("Barcode or gene identifiers are not unique")
}
if (!identical(dim(spatial_xy), c(shape[[1]], 2L)) || any(!is.finite(spatial_xy))) {
  stop("Spatial coordinates failed validation")
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

observed_sums <- as.numeric(Matrix::colSums(counts_gene_spot))
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
  interpretation = "X contains the 11702 genes shared with the cell2location reference; obs.total_counts retains the pre-intersection full-gene library size. SPARK-X and Hotspot use X counts and full X library size."
)

rm(counts_spot_gene, x_data, x_indices, x_indptr, observed_sums, expected_sums)
invisible(gc())
run_gate("sparkx_filter_genes")

detected_spots <- as.integer(Matrix::rowSums(counts_gene_spot > 0))
total_counts <- as.numeric(Matrix::rowSums(counts_gene_spot))
sum_squares <- as.numeric(Matrix::rowSums(counts_gene_spot ^ 2))
variance_numerator <- sum_squares - total_counts ^ 2 / shape[[1]]
minimum_detected_spots <- max(20L, as.integer(ceiling(0.01 * shape[[1]])))
keep <- detected_spots >= minimum_detected_spots & total_counts >= 50 & variance_numerator > 0
if (sum(keep) < 1000L) {
  stop("Gene filter retained too few genes: ", sum(keep))
}

gene_filter <- data.frame(
  gene = genes,
  detected_spots = detected_spots,
  detection_fraction = detected_spots / shape[[1]],
  total_counts = total_counts,
  nonzero_variance = variance_numerator > 0,
  tested = keep,
  stringsAsFactors = FALSE
)
write.table(
  gene_filter,
  file.path(output_dir, paste0(sample_id, "_gene_filter_qc.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

counts_tested <- counts_gene_spot[keep, , drop = FALSE]
if (any(Matrix::rowSums(counts_tested) <= 0) || any(Matrix::rowSums(counts_tested > 0) <= 1)) {
  stop("Filtered SPARK-X input contains an invalid gene")
}
rm(counts_gene_spot, sum_squares, variance_numerator)
invisible(gc())

run_gate("sparkx_mixture_kernels")
set.seed(seed)
spark_result <- SPARK::sparkx(
  counts_tested,
  spatial_xy,
  numCores = 4L,
  option = "mixture",
  verbose = FALSE
)
if (!all(c("stats", "res_stest", "res_mtest") %in% names(spark_result))) {
  stop("SPARK-X returned unexpected fields")
}
if (nrow(spark_result$res_mtest) != nrow(counts_tested)) {
  stop("SPARK-X result row count mismatch")
}

gene_names <- rownames(counts_tested)
res <- data.frame(
  gene = gene_names,
  sparkx_combined_p = as.numeric(spark_result$res_mtest$combinedPval),
  sparkx_by_fdr = as.numeric(spark_result$res_mtest$adjustedPval),
  sparkx_bh_fdr = p.adjust(as.numeric(spark_result$res_mtest$combinedPval), method = "BH"),
  sparkx_min_kernel_p = apply(spark_result$res_stest, 1, min),
  sparkx_max_kernel_stat = apply(spark_result$stats, 1, max),
  stringsAsFactors = FALSE
)
kernel_p <- as.data.frame(spark_result$res_stest)
kernel_stat <- as.data.frame(spark_result$stats)
names(kernel_p) <- paste0("sparkx_kernel_p__", names(kernel_p))
names(kernel_stat) <- paste0("sparkx_kernel_stat__", names(kernel_stat))
res <- cbind(res, kernel_p, kernel_stat)

if (any(!is.finite(as.matrix(res[, setdiff(names(res), "gene"), drop = FALSE]))) ||
    any(res$sparkx_combined_p < 0 | res$sparkx_combined_p > 1) ||
    any(res$sparkx_by_fdr < 0 | res$sparkx_by_fdr > 1)) {
  stop("SPARK-X result contains invalid values")
}
write.table(
  res,
  file.path(output_dir, paste0(sample_id, "_sparkx_gene_results.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

metadata <- list(
  status = "completed",
  sample = sample_id,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_h5ad = input_h5ad,
  input_sha256 = sha256_file(input_h5ad),
  n_spots = shape[[1]],
  n_genes_input = shape[[2]],
  nnz_input = nnz_input,
  total_counts_audit = total_counts_audit,
  minimum_detected_spots = minimum_detected_spots,
  minimum_total_counts = 50L,
  n_genes_tested = nrow(res),
  n_sparkx_by_fdr_0_05 = sum(res$sparkx_by_fdr < 0.05),
  n_sparkx_bh_fdr_0_05 = sum(res$sparkx_bh_fdr < 0.05),
  parameters = list(option = "mixture", kernels = 11L, num_cores = 4L, seed = seed),
  packages = list(R = R.version.string, SPARK = as.character(packageVersion("SPARK"))),
  inference_boundary = "Gene-level spatial variability within one section; spots are not biological replicates and significance is not a cohort-level claim."
)
write_json(metadata, file.path(output_dir, "sparkx_run_metadata.json"), auto_unbox = TRUE, pretty = TRUE)

completed <- TRUE
write_stage_status("completed", paste0("SPARK-X passed; genes_tested=", nrow(res)))
cat(sprintf("SPARKX_%s_PASS genes_tested=%d BY_FDR_0.05=%d\n", sample_id, nrow(res), sum(res$sparkx_by_fdr < 0.05)))
