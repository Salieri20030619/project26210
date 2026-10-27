#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot determine this script path from --file")
script_file <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
expected_root <- normalizePath("D:/project26210_external_joint_v1", winslash = "/", mustWork = TRUE)
if (!identical(tolower(root), tolower(expected_root))) {
  stop("Project root mismatch: actual=", root, " expected=", expected_root)
}

run_id <- Sys.getenv("PHASE5_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
expected_datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
result_dir <- file.path(root, "phase5_merge_by_dataset_local_v2")
object_dir <- file.path(result_dir, "objects_raw")
table_dir <- file.path(result_dir, "tables")
work_dir <- file.path(result_dir, paste0("work_slots_", run_id))
progress_file <- file.path(result_dir, paste0("PHASE5_PROGRESS_", run_id, ".tsv"))
status_file <- file.path(result_dir, "PHASE5_MERGE_STATUS")
lock_dir <- file.path(result_dir, "RUN_LOCK")
manifest_file <- file.path(root, "local_audit_v1", "manifest_reconciliation.tsv")
clinical_file <- file.path(root, "phase4e_clinical_metadata_harmonization_v1", "tables", "object_clinical_manifest.tsv")
prior_audit_file <- file.path(root, "phase4d_95_object_compatibility_audit_v1", "qc", "object_level_audit.tsv")

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", paste0(..., collapse = ""), "\n", sep = "")
  flush.console()
}

write_text_once <- function(lines, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing status file: ", path)
  part <- paste0(path, ".", run_id, ".part")
  writeLines(enc2utf8(lines), con = part, useBytes = TRUE)
  if (!file.rename(part, path)) stop("Atomic status rename failed: ", part, " -> ", path)
}

write_tsv_atomic <- function(x, path) {
  part <- paste0(path, ".", run_id, ".part")
  data.table::fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (file.exists(path)) stop("Refusing to overwrite existing table: ", path)
  if (!file.rename(part, path)) stop("Atomic table rename failed: ", part, " -> ", path)
}

append_tsv <- function(x, path) {
  exists <- file.exists(path)
  data.table::fwrite(x, path, sep = "\t", quote = TRUE, na = "NA", append = exists, col.names = !exists)
}

append_progress <- function(event, dataset = NA_character_, object_name = NA_character_, detail = "") {
  append_tsv(data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    run_id = run_id,
    event = event,
    dataset = dataset,
    object_name = object_name,
    detail = detail
  ), progress_file)
}

required_packages <- c("Seurat", "SeuratObject", "Matrix", "data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

if (file.exists(status_file)) stop("Existing PHASE5_MERGE_STATUS found; use a higher version directory instead of overwriting")
existing_outputs <- file.path(object_dir, paste0(expected_datasets, "_merged_raw.rds"))
if (any(file.exists(existing_outputs))) stop("Existing merged output found; refusing overwrite: ", paste(existing_outputs[file.exists(existing_outputs)], collapse = ", "))
fixed_tables <- file.path(table_dir, c(
  "object_merge_audit.tsv", "metadata_crosswalk.tsv", "feature_reconciliation.tsv",
  "cell_count_reconciliation.tsv", "merged_metadata_schema.tsv", "source_file_integrity.tsv",
  "resource_preflight.tsv", "dataset_merge_summary.tsv"
))
if (any(file.exists(fixed_tables))) stop("Existing Phase 5 audit tables found; use a higher version directory instead of overwriting")
if (dir.exists(work_dir)) stop("Run-specific work directory already exists: ", work_dir)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- data.table::fread(manifest_file, sep = "\t", encoding = "UTF-8", data.table = FALSE)
clinical <- data.table::fread(clinical_file, sep = "\t", encoding = "UTF-8", data.table = FALSE)
prior_audit <- data.table::fread(prior_audit_file, sep = "\t", encoding = "UTF-8", data.table = FALSE)

required_manifest_fields <- c("dataset", "object_name", "local_absolute_path", "manifest_cells", "exists", "size_matches_manifest")
if (!all(required_manifest_fields %in% names(manifest))) stop("Manifest reconciliation is missing required fields")
if (nrow(manifest) != 95L) stop("Local reconciled manifest does not contain 95 rows: ", nrow(manifest))
if (!setequal(unique(manifest$dataset), expected_datasets)) stop("Local reconciled manifest does not cover the five expected datasets")
if (anyDuplicated(manifest$object_name)) stop("Duplicate object_name in local reconciled manifest")
if (anyDuplicated(manifest$local_absolute_path)) stop("Duplicate local source path in local reconciled manifest")
logical_true <- function(x) tolower(as.character(x)) == "true"
if (!all(logical_true(manifest$exists)) || !all(logical_true(manifest$size_matches_manifest))) stop("Local input gate is not complete")
if (!all(file.exists(manifest$local_absolute_path))) stop("One or more reconciled local source files are now missing")
if (nrow(clinical) != 95L || anyDuplicated(paste(clinical$dataset, clinical$object_name, sep = "\r"))) stop("Clinical manifest is not a one-to-one 95-row object manifest")

source_before <- data.frame(
  dataset = manifest$dataset,
  object_name = manifest$object_name,
  source_file = manifest$local_absolute_path,
  size_before = file.info(manifest$local_absolute_path)$size,
  mtime_before = format(file.info(manifest$local_absolute_path)$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
  mtime_before_epoch = as.numeric(file.info(manifest$local_absolute_path)$mtime),
  stringsAsFactors = FALSE
)
source_integrity_written <- FALSE
source_modified_count <- NA_integer_

record_source_integrity <- function() {
  if (source_integrity_written) return(source_modified_count)
  after <- file.info(source_before$source_file)
  out <- transform(source_before,
    exists_after = !is.na(after$size),
    size_after = after$size,
    mtime_after = format(after$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
    unchanged = !is.na(after$size) & size_before == after$size & mtime_before == format(after$mtime, "%Y-%m-%dT%H:%OS6%z")
  )
  # Compare POSIXct numerically as the authoritative mtime check; formatted text is retained as evidence.
  after_num <- as.numeric(after$mtime)
  out$unchanged <- out$exists_after & out$size_before == out$size_after & out$mtime_before_epoch == after_num
  write_tsv_atomic(out, file.path(table_dir, "source_file_integrity.tsv"))
  source_modified_count <<- sum(!out$unchanged)
  source_integrity_written <<- TRUE
  source_modified_count
}

# Conservative preflight based on previously measured sparse nonzero counts, then revalidated during this run.
prior_audit$nonzero_values <- as.numeric(prior_audit$nonzero_values)
prior_audit$cells <- as.numeric(prior_audit$cells)
prior_audit$file_bytes <- as.numeric(prior_audit$file_bytes)
resource_rows <- do.call(rbind, lapply(expected_datasets, function(ds) {
  z <- prior_audit[prior_audit$dataset == ds, , drop = FALSE]
  matrix_bytes <- 12 * sum(z$nonzero_values) + 4 * (sum(z$cells) + 1)
  max_source_matrix <- max(12 * z$nonzero_values + 4 * (z$cells + 1))
  estimated_peak <- ceiling(1.6 * matrix_bytes + max_source_matrix + 1024^3)
  data.frame(dataset = ds, prior_nnz = sum(z$nonzero_values), estimated_merged_dgC_bytes = matrix_bytes,
             max_source_dgC_bytes = max_source_matrix, estimated_peak_memory_bytes = estimated_peak)
}))
available_memory <- suppressWarnings(as.numeric(Sys.getenv("PHASE5_AVAILABLE_MEMORY_BYTES", unset = NA_character_)))
free_disk <- suppressWarnings(as.numeric(Sys.getenv("PHASE5_FREE_DISK_BYTES", unset = NA_character_)))
max_peak <- max(resource_rows$estimated_peak_memory_bytes)
estimated_output_bytes <- ceiling(sum(resource_rows$estimated_merged_dgC_bytes) * 1.15)
estimated_work_bytes <- ceiling(max(resource_rows$estimated_merged_dgC_bytes) * 1.25)
resource_rows$available_memory_bytes_at_launch <- available_memory
resource_rows$free_disk_bytes_at_launch <- free_disk
resource_rows$memory_gate <- !is.na(available_memory) & available_memory >= 1.25 * resource_rows$estimated_peak_memory_bytes
resource_rows$disk_gate <- !is.na(free_disk) & free_disk >= estimated_output_bytes + estimated_work_bytes + 10 * 1024^3
write_tsv_atomic(resource_rows, file.path(table_dir, "resource_preflight.tsv"))
if (is.na(available_memory) || available_memory < 1.25 * max_peak) {
  stop("Memory preflight failed: available=", available_memory, " required_with_margin=", ceiling(1.25 * max_peak))
}
if (is.na(free_disk) || free_disk < estimated_output_bytes + estimated_work_bytes + 10 * 1024^3) {
  stop("Disk preflight failed: free=", free_disk, " required=", estimated_output_bytes + estimated_work_bytes + 10 * 1024^3)
}

count_violations <- function(x, chunk_size = 5000000L) {
  n <- length(x)
  negative <- nonfinite <- noninteger <- 0
  if (!n) return(c(negative = 0, nonfinite = 0, noninteger = 0))
  starts <- seq.int(1, n, by = chunk_size)
  for (s in starts) {
    e <- min(n, s + chunk_size - 1)
    v <- x[s:e]
    finite <- is.finite(v)
    nonfinite <- nonfinite + sum(!finite)
    if (any(finite)) {
      vf <- v[finite]
      negative <- negative + sum(vf < 0)
      noninteger <- noninteger + sum(abs(vf - round(vf)) > 1e-8)
    }
  }
  c(negative = negative, nonfinite = nonfinite, noninteger = noninteger)
}

safe_value <- function(x, fallback = "Unknown") {
  if (length(x) == 0L || is.na(x[1]) || !nzchar(trimws(as.character(x[1])))) fallback else as.character(x[1])
}

source_meta_field <- function(meta, candidates, fallback) {
  key <- tolower(gsub("[^a-z0-9]+", "", names(meta)))
  candidate_key <- tolower(gsub("[^a-z0-9]+", "", candidates))
  hit <- match(candidate_key, key, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) return(list(values = rep(fallback, nrow(meta)), source = "fallback"))
  idx <- hit[1]
  values <- as.character(meta[[idx]])
  values[is.na(values) | !nzchar(trimws(values))] <- fallback
  list(values = values, source = names(meta)[idx])
}

make_metadata <- function(obj_meta, counts, ds, object_name, source_file, clinical_row, source_assay, count_location) {
  if (!identical(rownames(obj_meta), colnames(counts))) stop("Original metadata row names are not identical to count-matrix cell names")
  n <- ncol(counts)
  barcode <- colnames(counts)
  cell_id <- paste(ds, object_name, barcode, sep = "__")
  if (anyDuplicated(cell_id)) stop("Generated cell_id is duplicated within object: ", object_name)
  sex <- source_meta_field(obj_meta, c("sex", "gender"), "Unknown")
  age <- source_meta_field(obj_meta, c("age", "age_years", "ageyears"), "Unknown")
  batch <- source_meta_field(obj_meta, c("batch_id", "batch", "sequencing_batch"), object_name)
  tissue_status <- safe_value(clinical_row$tumor_normal)
  disease_status <- safe_value(clinical_row$tissue)
  timepoint <- safe_value(clinical_row$timepoint_standard)
  if (identical(timepoint, "Unknown")) timepoint <- safe_value(clinical_row$sampling_timepoint)
  core <- data.frame(
    cell_id = cell_id,
    barcode_original = barcode,
    dataset = rep(ds, n),
    object_name = rep(object_name, n),
    sample_id = rep(safe_value(clinical_row$sample_id, object_name), n),
    patient_id = rep(safe_value(clinical_row$patient_id), n),
    tissue_status = rep(tissue_status, n),
    disease_status = rep(disease_status, n),
    sex = sex$values,
    age = age$values,
    treatment = rep(safe_value(clinical_row$treatment), n),
    timepoint = rep(timepoint, n),
    batch_id = batch$values,
    source_file = rep(source_file, n),
    source_assay = rep(source_assay, n),
    source_count_location = rep(count_location, n),
    orig.ident = rep(object_name, n),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  src_names <- make.unique(paste0("src__", names(obj_meta)), sep = "__dup")
  source_copy <- obj_meta
  names(source_copy) <- src_names
  rownames(core) <- cell_id
  rownames(source_copy) <- cell_id
  combined <- cbind(core, source_copy)
  if (anyDuplicated(names(combined))) stop("Metadata name collision after src__ prefixing")

  core_sources <- c(
    cell_id = "constructed:dataset+object_name+original_barcode", barcode_original = "count_matrix:colnames",
    dataset = "transfer_manifest:dataset", object_name = "transfer_manifest:object_name",
    sample_id = "clinical_manifest:sample_id", patient_id = "clinical_manifest:patient_id",
    tissue_status = "clinical_manifest:tumor_normal", disease_status = "clinical_manifest:tissue",
    sex = paste0("source_metadata:", sex$source), age = paste0("source_metadata:", age$source),
    treatment = "clinical_manifest:treatment", timepoint = "clinical_manifest:timepoint_standard_or_sampling_timepoint",
    batch_id = paste0("source_metadata:", batch$source), source_file = "local_reconciled_manifest:local_absolute_path",
    source_assay = "object_audit:confirmed_assay", source_count_location = "object_layer:counts",
    orig.ident = "standardized:object_name"
  )
  core_crosswalk <- data.frame(
    dataset = ds, object_name = object_name, output_field = names(core),
    source_scope = sub(":.*$", "", core_sources[names(core)]),
    source_field = sub("^[^:]+:", "", core_sources[names(core)]),
    mapping_rule = c(
      "paste(dataset, object_name, barcode_original, sep='__')", "exact count-matrix column name",
      rep("exact manifest value repeated per cell", 6),
      "exact source field when present, otherwise Unknown", "exact source field when present, otherwise Unknown",
      "exact clinical manifest value", "prefer standardized timepoint, then sampling timepoint",
      "exact source field when present, otherwise object_name", "absolute local read path",
      "confirmed raw-count assay", "confirmed Seurat v5 counts layer", "object_name"
    ), stringsAsFactors = FALSE
  )
  source_crosswalk <- data.frame(
    dataset = ds, object_name = object_name, output_field = src_names,
    source_scope = "original_object_metadata", source_field = names(obj_meta),
    mapping_rule = "verbatim value copy; src__ prefix; make.unique collision suffix",
    stringsAsFactors = FALSE
  )
  list(metadata = combined, crosswalk = rbind(core_crosswalk, source_crosswalk))
}

object_audit_path <- file.path(table_dir, "object_merge_audit.tsv")
crosswalk_path <- file.path(table_dir, "metadata_crosswalk.tsv")
feature_path <- file.path(table_dir, "feature_reconciliation.tsv")
cell_path <- file.path(table_dir, "cell_count_reconciliation.tsv")
schema_path <- file.path(table_dir, "merged_metadata_schema.tsv")

completed_datasets <- character()
dataset_summaries <- list()
all_output_readable <- TRUE

process_dataset <- function(ds) {
  ds_manifest <- manifest[manifest$dataset == ds, , drop = FALSE]
  ds_manifest <- ds_manifest[order(ds_manifest$object_name), , drop = FALSE]
  ds_clinical <- clinical[clinical$dataset == ds, , drop = FALSE]
  ds_work <- file.path(work_dir, ds)
  dir.create(ds_work, recursive = TRUE, showWarnings = FALSE)
  log_line("DATASET_START dataset=", ds, " objects=", nrow(ds_manifest), " manifest_cells=", sum(as.numeric(ds_manifest$manifest_cells)))
  append_progress("DATASET_START", ds, detail = paste0("objects=", nrow(ds_manifest), ";manifest_cells=", sum(as.numeric(ds_manifest$manifest_cells))))

  feature_lists <- vector("list", nrow(ds_manifest))
  object_stats <- vector("list", nrow(ds_manifest))
  chunk_files <- character(nrow(ds_manifest))
  global_cell_ids <- character()

  for (j in seq_len(nrow(ds_manifest))) {
    row <- ds_manifest[j, , drop = FALSE]
    object_name <- row$object_name
    source_file <- row$local_absolute_path
    clinical_row <- ds_clinical[ds_clinical$object_name == object_name, , drop = FALSE]
    if (nrow(clinical_row) != 1L) stop("Clinical manifest match is not exactly one row for ", ds, "/", object_name)
    log_line("OBJECT_LOAD_START dataset=", ds, " index=", j, "/", nrow(ds_manifest), " object=", object_name)
    append_progress("OBJECT_LOAD_START", ds, object_name, paste0("index=", j, "/", nrow(ds_manifest)))
    audit_row <- data.frame(
      dataset = ds, object_name = object_name, source_file = source_file, source_assay = NA_character_,
      source_count_location = NA_character_, manifest_cells = as.numeric(row$manifest_cells), loaded_cells = NA_real_,
      features = NA_real_, nonzero_values = NA_real_, negative_count_values = NA_real_,
      non_finite_count_values = NA_real_, non_integer_count_values = NA_real_, metadata_field_count = NA_real_,
      staged_chunk_file = NA_character_, status = "FAIL", failure_reason = "unclassified error",
      stringsAsFactors = FALSE
    )
    ok <- FALSE
    tryCatch({
      obj <- readRDS(source_file)
      if (!inherits(obj, "Seurat")) stop("readRDS object is not a Seurat object; class=", paste(class(obj), collapse = ","))
      assays <- SeuratObject::Assays(obj)
      source_assay <- if ("RNA" %in% assays) "RNA" else SeuratObject::DefaultAssay(obj)
      if (!source_assay %in% assays) stop("No usable assay in object")
      layers <- SeuratObject::Layers(obj[[source_assay]])
      if (!"counts" %in% layers) stop("Exact raw counts layer is absent; layers=", paste(layers, collapse = ","))
      count_like_layers <- grep("^counts($|\\.)", layers, value = TRUE)
      if (!identical(count_like_layers, "counts")) stop("Ambiguous source counts layers: ", paste(count_like_layers, collapse = ","))
      count_location <- "LayerData:counts"
      counts <- SeuratObject::LayerData(obj, assay = source_assay, layer = "counts")
      if (!inherits(counts, "dgCMatrix")) stop("Raw counts are not dgCMatrix; class=", paste(class(counts), collapse = ","))
      if (is.null(rownames(counts)) || is.null(colnames(counts))) stop("Raw count matrix lacks feature or cell names")
      if (anyDuplicated(rownames(counts))) stop("Duplicate raw-count feature names")
      if (anyDuplicated(colnames(counts))) stop("Duplicate raw-count cell names within object")
      obj_meta <- obj[[]]
      if (nrow(obj_meta) != ncol(counts) || !identical(rownames(obj_meta), colnames(counts))) stop("Metadata/count cell alignment failure")
      violations <- count_violations(counts@x)
      if (any(violations > 0)) stop("Raw-count value violation: negative=", violations[["negative"]], ", nonfinite=", violations[["nonfinite"]], ", noninteger=", violations[["noninteger"]])
      if (ncol(counts) != as.numeric(row$manifest_cells)) stop("Cell count mismatch: manifest=", row$manifest_cells, " loaded=", ncol(counts))
      prior_row <- prior_audit[prior_audit$dataset == ds & prior_audit$object_name == object_name, , drop = FALSE]
      if (nrow(prior_row) != 1L) stop("Prior audit match is not exactly one row")
      if (as.numeric(prior_row$nonzero_values) != length(counts@x)) stop("Nonzero count differs from transferred object audit")
      if (as.numeric(prior_row$genes) != nrow(counts)) stop("Feature count differs from transferred object audit")

      meta_result <- make_metadata(obj_meta, counts, ds, object_name, source_file, clinical_row, source_assay, count_location)
      if (length(intersect(global_cell_ids, meta_result$metadata$cell_id))) stop("Generated cell_id overlaps a prior object in dataset")
      global_cell_ids <- c(global_cell_ids, meta_result$metadata$cell_id)
      feature_lists[[j]] <- rownames(counts)
      chunk <- list(
        dataset = ds, object_name = object_name, features = rownames(counts),
        i = counts@i, p = counts@p, x = counts@x,
        nrow = nrow(counts), ncol = ncol(counts), metadata = meta_result$metadata
      )
      chunk_file <- file.path(ds_work, sprintf("%03d_%s_slots.rds", j, gsub("[^A-Za-z0-9_.-]", "_", object_name)))
      chunk_part <- paste0(chunk_file, ".part")
      saveRDS(chunk, chunk_part, compress = FALSE, version = 3)
      if (!file.rename(chunk_part, chunk_file)) stop("Cannot atomically finalize staged slot file")
      chunk_files[j] <- chunk_file
      append_tsv(meta_result$crosswalk, crosswalk_path)

      audit_row$source_assay <- source_assay
      audit_row$source_count_location <- count_location
      audit_row$loaded_cells <- ncol(counts)
      audit_row$features <- nrow(counts)
      audit_row$nonzero_values <- length(counts@x)
      audit_row$negative_count_values <- violations[["negative"]]
      audit_row$non_finite_count_values <- violations[["nonfinite"]]
      audit_row$non_integer_count_values <- violations[["noninteger"]]
      audit_row$metadata_field_count <- ncol(obj_meta)
      audit_row$staged_chunk_file <- chunk_file
      audit_row$status <- "PASS"
      audit_row$failure_reason <- ""
      object_stats[[j]] <- audit_row
      ok <- TRUE
      rm(chunk, counts, obj, obj_meta, meta_result)
      gc(verbose = FALSE)
    }, error = function(e) {
      audit_row$failure_reason <<- conditionMessage(e)
    })
    append_tsv(audit_row, object_audit_path)
    if (!ok) {
      append_progress("OBJECT_FAIL", ds, object_name, audit_row$failure_reason)
      stop("Object audit failed for ", ds, "/", object_name, ": ", audit_row$failure_reason)
    }
    append_progress("OBJECT_PASS", ds, object_name, paste0("cells=", audit_row$loaded_cells, ";features=", audit_row$features, ";nnz=", audit_row$nonzero_values))
    log_line("OBJECT_PASS dataset=", ds, " object=", object_name, " cells=", audit_row$loaded_cells, " features=", audit_row$features, " nnz=", audit_row$nonzero_values)
  }

  feature_union <- unique(unlist(feature_lists, use.names = FALSE))
  feature_intersection <- Reduce(intersect, feature_lists)
  first_features <- feature_lists[[1]]
  merged_cells_expected <- sum(vapply(object_stats, function(x) x$loaded_cells, numeric(1)))
  manifest_cells_expected <- sum(as.numeric(ds_manifest$manifest_cells))
  total_nnz <- sum(vapply(object_stats, function(x) x$nonzero_values, numeric(1)))
  log_line("MATRIX_ASSEMBLY_START dataset=", ds, " cells=", merged_cells_expected, " union_features=", length(feature_union), " nnz=", total_nnz)
  append_progress("MATRIX_ASSEMBLY_START", ds, detail = paste0("cells=", merged_cells_expected, ";features=", length(feature_union), ";nnz=", total_nnz))

  x_out <- numeric(total_nnz)
  i_out <- integer(total_nnz)
  p_out <- integer(merged_cells_expected + 1)
  meta_list <- vector("list", length(chunk_files))
  nnz_offset <- 0L
  col_offset <- 0L
  for (j in seq_along(chunk_files)) {
    chunk <- readRDS(chunk_files[j])
    row_map <- match(chunk$features, feature_union)
    if (anyNA(row_map)) stop("Internal feature-union mapping failure")
    chunk_i <- chunk$i
    chunk_p <- chunk$p
    chunk_x <- chunk$x
    remapped_needed <- FALSE
    if (length(row_map) > 1L && any(diff(row_map) <= 0L)) {
      old_matrix <- methods::new("dgCMatrix", i = chunk$i, p = chunk$p, x = chunk$x,
        Dim = as.integer(c(chunk$nrow, chunk$ncol)), Dimnames = list(chunk$features, NULL))
      row_projection <- Matrix::sparseMatrix(i = row_map, j = seq_along(row_map), x = 1,
        dims = c(length(feature_union), length(row_map)))
      remapped <- methods::as(row_projection %*% old_matrix, "dgCMatrix")
      chunk_i <- remapped@i
      chunk_p <- remapped@p
      chunk_x <- remapped@x
      remapped_needed <- TRUE
      rm(old_matrix, row_projection, remapped)
    }
    n_chunk_nnz <- length(chunk_x)
    if (n_chunk_nnz) {
      idx <- seq.int(nnz_offset + 1, nnz_offset + n_chunk_nnz)
      x_out[idx] <- chunk_x
      if (!remapped_needed) {
        i_out[idx] <- row_map[chunk_i + 1L] - 1L
      } else {
        i_out[idx] <- chunk_i
      }
    }
    if (chunk$ncol) {
      p_idx <- seq.int(col_offset + 2L, col_offset + chunk$ncol + 1L)
      p_out[p_idx] <- as.integer(chunk_p[-1L] + nnz_offset)
    }
    meta_list[[j]] <- chunk$metadata
    nnz_offset <- as.integer(nnz_offset + n_chunk_nnz)
    col_offset <- col_offset + chunk$ncol
    rm(chunk, row_map)
    gc(verbose = FALSE)
  }
  if (nnz_offset != total_nnz || col_offset != merged_cells_expected || tail(p_out, 1) != total_nnz) stop("Sparse slot assembly conservation failure")
  merged_meta <- as.data.frame(data.table::rbindlist(meta_list, use.names = TRUE, fill = TRUE), check.names = FALSE)
  rownames(merged_meta) <- merged_meta$cell_id
  rm(meta_list)
  gc(verbose = FALSE)
  counts_merged <- methods::new("dgCMatrix", i = i_out, p = p_out, x = x_out,
    Dim = as.integer(c(length(feature_union), merged_cells_expected)),
    Dimnames = list(feature_union, merged_meta$cell_id))
  methods::validObject(counts_merged)
  rm(i_out, p_out, x_out)
  gc(verbose = FALSE)

  merged <- SeuratObject::CreateSeuratObject(counts = counts_merged, assay = "RNA", meta.data = merged_meta, project = ds)
  merged_layers <- SeuratObject::Layers(merged[["RNA"]])
  if (!identical(merged_layers, "counts")) stop("Merged object does not contain exactly one counts layer: ", paste(merged_layers, collapse = ","))
  if (length(SeuratObject::Reductions(merged)) != 0L || length(merged@graphs) != 0L || length(merged@neighbors) != 0L) stop("Merged raw object unexpectedly contains reductions, graphs, or neighbors")
  merged_counts <- SeuratObject::LayerData(merged, assay = "RNA", layer = "counts")
  if (ncol(merged_counts) != merged_cells_expected || ncol(merged_counts) != manifest_cells_expected) stop("Merged/manifest cell conservation failure")
  if (anyDuplicated(colnames(merged_counts)) || anyDuplicated(merged$cell_id)) stop("Merged cell_id is not unique")
  if (!identical(colnames(merged_counts), rownames(merged[[]])) || !identical(rownames(merged[[]]), unname(as.character(merged$cell_id)))) stop("Merged metadata row names are not identical to cell_id")
  if (!identical(rownames(merged_counts), feature_union)) stop("Merged feature set/order differs from constructed union")
  if (!identical(unique(as.character(merged$dataset)), ds)) stop("Merged dataset field is not restricted to current dataset")

  output_file <- file.path(object_dir, paste0(ds, "_merged_raw.rds"))
  output_part <- paste0(output_file, ".", run_id, ".part")
  if (file.exists(output_file) || file.exists(output_part)) stop("Output or part file already exists: ", output_file)
  saveRDS(merged, output_part, compress = FALSE, version = 3)
  if (!file.rename(output_part, output_file)) stop("Atomic output rename failed: ", output_part, " -> ", output_file)
  output_bytes <- file.info(output_file)$size
  rm(merged, merged_counts, counts_merged, merged_meta)
  gc(verbose = FALSE)

  reread <- readRDS(output_file)
  reread_counts <- SeuratObject::LayerData(reread, assay = "RNA", layer = "counts")
  reread_ok <- inherits(reread, "Seurat") && ncol(reread_counts) == merged_cells_expected &&
    nrow(reread_counts) == length(feature_union) && identical(colnames(reread_counts), rownames(reread[[]])) &&
    identical(rownames(reread[[]]), unname(as.character(reread$cell_id))) &&
    identical(rownames(reread_counts), feature_union) && identical(SeuratObject::Layers(reread[["RNA"]]), "counts")
  if (!reread_ok) stop("Final output failed readRDS/raw-count validation: ", output_file)

  feature_rows <- do.call(rbind, lapply(seq_along(feature_lists), function(j) data.frame(
    dataset = ds, object_name = ds_manifest$object_name[j], features = length(feature_lists[[j]]),
    feature_union = length(feature_union), feature_intersection = length(feature_intersection),
    same_order_as_dataset_first = identical(feature_lists[[j]], first_features),
    same_set_as_dataset_first = setequal(feature_lists[[j]], first_features),
    merged_features = nrow(reread_counts), stringsAsFactors = FALSE
  )))
  append_tsv(feature_rows, feature_path)
  cell_rows <- do.call(rbind, lapply(seq_along(object_stats), function(j) data.frame(
    record_type = "source_object", dataset = ds, object_name = object_stats[[j]]$object_name,
    manifest_cells = object_stats[[j]]$manifest_cells, loaded_cells = object_stats[[j]]$loaded_cells,
    merged_cells = NA_real_, conserved = object_stats[[j]]$manifest_cells == object_stats[[j]]$loaded_cells,
    stringsAsFactors = FALSE
  )))
  cell_rows <- rbind(cell_rows, data.frame(
    record_type = "dataset_total", dataset = ds, object_name = "__DATASET_TOTAL__",
    manifest_cells = manifest_cells_expected, loaded_cells = merged_cells_expected,
    merged_cells = ncol(reread_counts),
    conserved = manifest_cells_expected == merged_cells_expected && merged_cells_expected == ncol(reread_counts),
    stringsAsFactors = FALSE
  ))
  append_tsv(cell_rows, cell_path)
  schema <- data.frame(
    dataset = ds, field = names(reread[[]]), class = vapply(reread[[]], function(x) paste(class(x), collapse = ","), character(1)),
    non_missing_values = vapply(reread[[]], function(x) sum(!is.na(x)), integer(1)),
    unique_values = vapply(reread[[]], function(x) data.table::uniqueN(x, na.rm = FALSE), integer(1)),
    stringsAsFactors = FALSE
  )
  append_tsv(schema, schema_path)
  dataset_summary <- data.frame(
    dataset = ds, source_objects = nrow(ds_manifest), manifest_cells = manifest_cells_expected,
    loaded_cells = merged_cells_expected, merged_cells = ncol(reread_counts),
    feature_union = length(feature_union), feature_intersection = length(feature_intersection),
    merged_features = nrow(reread_counts), output_file = output_file, output_bytes = output_bytes,
    output_readRDS_pass = reread_ok, counts_layers = paste(SeuratObject::Layers(reread[["RNA"]]), collapse = ","),
    normalization_performed = 0L, clustering_performed = 0L, stringsAsFactors = FALSE
  )
  append_tsv(dataset_summary, file.path(table_dir, "dataset_merge_summary.tsv"))
  rm(reread, reread_counts)
  gc(verbose = FALSE)
  unlink(ds_work, recursive = TRUE, force = TRUE)
  append_progress("DATASET_PASS", ds, detail = paste0("cells=", merged_cells_expected, ";features=", length(feature_union), ";output_bytes=", output_bytes))
  log_line("DATASET_PASS dataset=", ds, " cells=", merged_cells_expected, " features=", length(feature_union), " output_bytes=", output_bytes)
  dataset_summary
}

if (!dir.create(lock_dir, showWarnings = FALSE)) stop("Phase 5 run lock already exists: ", lock_dir)
writeLines(c(paste0("run_id=", run_id), paste0("pid=", Sys.getpid()), paste0("started_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))), file.path(lock_dir, "owner.txt"))

main_result <- tryCatch({
  log_line("PHASE5_START run_id=", run_id, " root=", root, " source_objects=", nrow(manifest))
  append_progress("PHASE5_START", detail = paste0("source_objects=", nrow(manifest)))
  for (ds in expected_datasets) {
    dataset_summaries[[ds]] <- process_dataset(ds)
    completed_datasets <- c(completed_datasets, ds)
    gc(verbose = FALSE)
  }
  summary_table <- do.call(rbind, dataset_summaries)
  all_outputs_exist <- all(file.exists(file.path(object_dir, paste0(expected_datasets, "_merged_raw.rds"))))
  all_loaded_cells_conserved <- all(summary_table$loaded_cells == summary_table$merged_cells)
  all_manifest_cells_conserved <- all(summary_table$manifest_cells == summary_table$merged_cells)
  all_output_readable <- all(summary_table$output_readRDS_pass)
  source_modified_count <- record_source_integrity()
  pass <- length(completed_datasets) == 5L && nrow(manifest) == 95L && all_loaded_cells_conserved &&
    all_manifest_cells_conserved && all_outputs_exist && all_output_readable && source_modified_count == 0L
  if (!pass) stop("Final Phase 5 invariant check failed")
  write_text_once(c(
    "status=PASS",
    paste0("run_id=", run_id),
    "datasets_expected=5",
    paste0("datasets_completed=", length(completed_datasets)),
    paste0("source_objects=", nrow(manifest)),
    paste0("all_loaded_cells_conserved=", toupper(all_loaded_cells_conserved)),
    paste0("all_manifest_cells_conserved=", toupper(all_manifest_cells_conserved)),
    paste0("all_outputs_exist=", toupper(all_outputs_exist)),
    paste0("all_outputs_readRDS=", toupper(all_output_readable)),
    paste0("source_rds_modified=", source_modified_count),
    "normalization_performed=0",
    "clustering_performed=0",
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  ), status_file)
  append_progress("PHASE5_PASS", detail = "datasets_completed=5")
  log_line("PHASE5_PASS datasets_completed=5 source_objects=", nrow(manifest))
  0L
}, error = function(e) {
  msg <- conditionMessage(e)
  log_line("PHASE5_FAIL reason=", msg)
  try(append_progress("PHASE5_FAIL", detail = msg), silent = TRUE)
  if (!source_integrity_written) try(record_source_integrity(), silent = TRUE)
  if (!file.exists(status_file)) {
    try(write_text_once(c(
      "status=FAIL",
      paste0("run_id=", run_id),
      "datasets_expected=5",
      paste0("datasets_completed=", length(completed_datasets)),
      paste0("source_objects=", nrow(manifest)),
      paste0("failure_reason=", gsub("[\r\n]+", " ", msg)),
      paste0("source_rds_modified=", ifelse(is.na(source_modified_count), "UNKNOWN", source_modified_count)),
      "normalization_performed=0",
      "clustering_performed=0",
      paste0("failed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    ), status_file), silent = TRUE)
  }
  1L
})

unlink(lock_dir, recursive = TRUE, force = TRUE)
if (dir.exists(lock_dir)) {
  log_line("LOCK_CLEANUP_FAIL lock_dir=", lock_dir)
  main_result <- 1L
}
quit(save = "no", status = main_result, runLast = FALSE)
