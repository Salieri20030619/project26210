options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: validate_phase6_processing_local_v1.R <root> <run_id>")
root <- normalizePath(args[[1]], winslash = "\\", mustWork = TRUE)
run_id <- args[[2]]
run_dir <- file.path(root, "phase6_processing_local_v1", "runs", run_id)
table_dir <- file.path(run_dir, "tables")
object_dir <- file.path(run_dir, "objects_processed")
summary_file <- file.path(table_dir, "processing_summary.tsv")
phase5_summary_file <- file.path(root, "phase5_merge_by_dataset_local_v2", "tables", "dataset_merge_summary.tsv")
validation_file <- file.path(table_dir, "processing_output_validation.tsv")
software_file <- file.path(table_dir, "processing_validation_software.tsv")
status_file <- file.path(run_dir, "PHASE6_PROCESSING_VALIDATION_STATUS")

for (path in c(validation_file, software_file, status_file)) {
  if (file.exists(path)) stop("Refusing to overwrite existing validation output: ", path)
}

write_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

write_lines_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  writeLines(x, part, useBytes = TRUE)
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

expected <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
summary <- fread(summary_file, sep = "\t", encoding = "UTF-8")
phase5 <- fread(phase5_summary_file, sep = "\t", encoding = "UTF-8")
if (!setequal(summary$dataset, expected) || nrow(summary) != 5L) stop("Processing summary dataset gate failed")
if (!setequal(phase5$dataset, expected) || nrow(phase5) != 5L) stop("Phase 5 summary dataset gate failed")

rows <- vector("list", length(expected))
for (i in seq_along(expected)) {
  ds <- expected[[i]]
  message("PHASE6_VALIDATION_DATASET_START=", ds)
  s <- summary[dataset == ds]
  p5 <- phase5[dataset == ds]
  output <- file.path(object_dir, paste0(ds, "_processed_clustered.rds"))
  input <- file.path(root, "phase5_merge_by_dataset_local_v2", "objects_raw", paste0(ds, "_merged_raw.rds"))
  output_exists <- file.exists(output)
  input_exists <- file.exists(input)
  output_bytes_match <- output_exists && identical(as.numeric(file.info(output)$size), as.numeric(s$output_bytes))
  input_bytes_match <- input_exists && identical(as.numeric(file.info(input)$size), as.numeric(p5$output_bytes))
  if (!output_exists || !input_exists || !output_bytes_match || !input_bytes_match) stop(ds, ": file existence/size gate failed")

  obj <- readRDS(output)
  is_seurat <- inherits(obj, "Seurat")
  if (!is_seurat) stop(ds, ": output is not a Seurat object")
  assay <- DefaultAssay(obj)
  layers <- Layers(obj[[assay]])
  reductions <- names(obj@reductions)
  cluster_cols <- grep("_snn_res\\.(0\\.4|0\\.8)$", colnames(obj[[]]), value = TRUE)
  counts <- LayerData(obj, assay = assay, layer = "counts")
  count_values <- counts@x
  counts_finite <- all(is.finite(count_values))
  counts_nonnegative <- all(count_values >= 0)
  counts_integer <- all(count_values == floor(count_values))
  data_values <- LayerData(obj, assay = assay, layer = "data")@x
  data_finite <- all(is.finite(data_values))
  cells_match <- ncol(obj) == s$cells_after
  features_match <- nrow(obj) == s$features
  cells_unique <- !anyDuplicated(colnames(obj))
  qc_all_pass <- "phase6_qc_pass" %in% colnames(obj[[]]) && all(obj$phase6_qc_pass %in% TRUE)
  required_layers <- all(c("counts", "data") %in% layers)
  reductions_pass <- all(c("pca", "umap") %in% reductions)
  cluster_fields_pass <- length(cluster_cols) == 2L
  no_integrated_assay <- !any(tolower(names(obj@assays)) == "integrated")
  misc_run_id_match <- identical(as.character(obj@misc$phase6_processing$run_id), run_id)
  pass <- all(c(output_exists, input_exists, output_bytes_match, input_bytes_match, is_seurat,
                cells_match, features_match, cells_unique, qc_all_pass, required_layers,
                reductions_pass, cluster_fields_pass, no_integrated_assay, misc_run_id_match,
                counts_finite, counts_nonnegative, counts_integer, data_finite))

  rows[[i]] <- data.table(
    dataset = ds,
    output_file = normalizePath(output, winslash = "\\", mustWork = TRUE),
    output_bytes = as.numeric(file.info(output)$size),
    output_bytes_match = output_bytes_match,
    input_file = normalizePath(input, winslash = "\\", mustWork = TRUE),
    input_bytes = as.numeric(file.info(input)$size),
    input_bytes_match_phase5 = input_bytes_match,
    readRDS_pass = is_seurat,
    cells = ncol(obj),
    cells_match_summary = cells_match,
    features = nrow(obj),
    features_match_summary = features_match,
    cells_unique = cells_unique,
    qc_all_pass = qc_all_pass,
    assay = assay,
    layers = paste(layers, collapse = ","),
    required_layers_pass = required_layers,
    reductions = paste(reductions, collapse = ","),
    reductions_pass = reductions_pass,
    cluster_fields = paste(cluster_cols, collapse = ","),
    cluster_fields_pass = cluster_fields_pass,
    no_integrated_assay = no_integrated_assay,
    misc_run_id_match = misc_run_id_match,
    counts_nonzero_values = length(count_values),
    counts_finite = counts_finite,
    counts_nonnegative = counts_nonnegative,
    counts_integer = counts_integer,
    data_finite = data_finite,
    status = if (pass) "PASS" else "FAIL"
  )
  rm(obj, counts, count_values, data_values)
  invisible(gc())
  if (!pass) stop(ds, ": independent output validation failed")
  message("PHASE6_VALIDATION_DATASET_PASS=", ds)
}

validation <- rbindlist(rows)
write_atomic(validation, validation_file)
software <- data.table(
  component = c("R", "SeuratObject", "Matrix", "data.table"),
  version = c(as.character(getRversion()), as.character(packageVersion("SeuratObject")),
              as.character(packageVersion("Matrix")), as.character(packageVersion("data.table")))
)
write_atomic(software, software_file)
status <- c(
  "status=PASS",
  paste0("run_id=", run_id),
  "datasets_expected=5",
  paste0("datasets_validated=", nrow(validation)),
  paste0("cells_validated=", sum(validation$cells)),
  paste0("all_outputs_readRDS=", all(validation$readRDS_pass)),
  paste0("all_output_sizes_match=", all(validation$output_bytes_match)),
  paste0("all_phase5_input_sizes_match=", all(validation$input_bytes_match_phase5)),
  paste0("all_counts_safe=", all(validation$counts_finite & validation$counts_nonnegative & validation$counts_integer)),
  paste0("all_required_layers=", all(validation$required_layers_pass)),
  paste0("all_required_reductions=", all(validation$reductions_pass)),
  paste0("all_cluster_fields=", all(validation$cluster_fields_pass)),
  paste0("all_no_integrated_assay=", all(validation$no_integrated_assay)),
  paste0("validated_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
write_lines_atomic(status, status_file)
message("PHASE6_PROCESSING_VALIDATION=PASS")
