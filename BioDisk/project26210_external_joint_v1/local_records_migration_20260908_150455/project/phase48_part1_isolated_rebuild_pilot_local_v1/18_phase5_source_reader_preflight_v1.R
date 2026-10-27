#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 18_phase5_source_reader_preflight_v1.R REBUILD_ROOT")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
scope_file <- file.path(root, "00_scope_and_provenance", "full_cohort_library_scope_v3.tsv")
inventory_file <- file.path(root, "01_input_inventory", "input_file_inventory_v2.tsv")
archive_dir <- file.path(root, "01_input_inventory", "public_geo_archives_v1")
h5ad_file <- file.path(root, "02_complete_cohort_build", "gse242299_h5ad_recovery_v3", "GSE242299_all_cells_50236_33538.h5ad")
final_dir <- file.path(root, "02_complete_cohort_build", "phase5_source_reader_preflight_v1")
part_dir <- paste0(final_dir, ".part")

if (dir.exists(final_dir) || dir.exists(part_dir)) stop("Refusing existing preflight output directory")
dir.create(part_dir, recursive = FALSE)
work_dir <- file.path(part_dir, "work")
dir.create(work_dir)

required <- c("Matrix", "data.table", "Seurat", "SeuratObject", "rhdf5", "R.utils", "jsonlite", "digest")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ","))

write_tsv <- function(x, path) {
  part <- paste0(path, ".part")
  data.table::fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (file.exists(path) || !file.rename(part, path)) stop("Cannot atomically write ", path)
}

count_violations <- function(x) {
  finite <- is.finite(x)
  c(
    nonfinite = sum(!finite),
    negative = if (any(finite)) sum(x[finite] < 0) else 0,
    noninteger = if (any(finite)) sum(x[finite] != floor(x[finite])) else 0
  )
}

matrix_row <- function(source_type, study, library_id, counts, original_metadata_fields, detail) {
  if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "dgCMatrix")
  methods::validObject(counts)
  violations <- count_violations(counts@x)
  data.frame(
    source_type = source_type,
    study = study,
    library_id = library_id,
    features = nrow(counts),
    cells = ncol(counts),
    nnz = length(counts@x),
    feature_names_unique = !anyDuplicated(rownames(counts)),
    cell_names_unique = !anyDuplicated(colnames(counts)),
    nonfinite_values = violations[["nonfinite"]],
    negative_values = violations[["negative"]],
    noninteger_values = violations[["noninteger"]],
    original_metadata_fields = original_metadata_fields,
    detail = detail,
    status = if (
      nrow(counts) > 0L && ncol(counts) > 0L &&
      !anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)) &&
      all(violations == 0)
    ) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}

scope <- data.table::fread(scope_file, sep = "\t", encoding = "UTF-8", data.table = FALSE)
scope <- scope[scope$inclusion_status == "include", , drop = FALSE]
inventory <- data.table::fread(inventory_file, sep = "\t", encoding = "UTF-8", data.table = FALSE)
if (nrow(scope) != 99L || nrow(inventory) != 45L) stop("Scope or RDS inventory count mismatch")
rows <- list()

# Existing RDS reader: one representative object from each RDS-backed dataset.
for (study in c("GSE178481", "GSE207493")) {
  item <- inventory[inventory$dataset == study, , drop = FALSE][1L, , drop = FALSE]
  obj <- readRDS(item$input_path)
  if (!inherits(obj, "Seurat")) stop("Representative RDS is not Seurat: ", item$input_path)
  assay <- if ("RNA" %in% SeuratObject::Assays(obj)) "RNA" else SeuratObject::DefaultAssay(obj)
  layers <- SeuratObject::Layers(obj[[assay]])
  if (!"counts" %in% layers) stop("Representative RDS lacks counts layer")
  counts <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  meta <- obj[[]]
  if (!identical(colnames(counts), rownames(meta))) stop("RDS count/metadata alignment failure")
  library_id <- sub("\\.rds$", "", basename(item$input_path))
  rows[[length(rows) + 1L]] <- matrix_row("Seurat_RDS", study, library_id, counts, ncol(meta), paste0("assay=", assay, ";layers=", paste(layers, collapse = "|")))
  rm(obj, counts, meta)
  gc(verbose = FALSE)
}

# GSE156632 dense CSV.gz reader from the formal archive.
csv_scope <- scope[scope$study == "GSE156632", , drop = FALSE][1L, , drop = FALSE]
csv_member <- csv_scope$archive_members
csv_archive <- file.path(archive_dir, "GSE156632_RAW.tar")
utils::untar(csv_archive, files = csv_member, exdir = work_dir)
csv_path <- file.path(work_dir, csv_member)
csv <- data.table::fread(csv_path, sep = ",", encoding = "UTF-8", data.table = FALSE)
if (ncol(csv) < 3L || !names(csv)[1] %in% c("Gene_ID", "Gene-ID") || names(csv)[2] != "Symbol") stop("Unexpected GSE156632 CSV schema")
feature_id <- as.character(csv[[1]])
if (anyNA(feature_id) || any(!nzchar(feature_id)) || anyDuplicated(feature_id)) stop("GSE156632 Gene_ID is empty or duplicated")
dense <- as.matrix(csv[, -(1:2), drop = FALSE])
storage.mode(dense) <- "double"
counts <- methods::as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix")
rownames(counts) <- feature_id
colnames(counts) <- names(csv)[-(1:2)]
rows[[length(rows) + 1L]] <- matrix_row("dense_CSV_gzip", "GSE156632", csv_scope$library_id, counts, 0L, paste0("feature_id=", names(csv)[1], ";feature_symbol=", names(csv)[2]))
rm(csv, dense, counts)
unlink(csv_path)
gc(verbose = FALSE)

# GSE159115 10x H5 reader from the formal archive.
h5_scope <- scope[scope$study == "GSE159115", , drop = FALSE][1L, , drop = FALSE]
h5_member <- h5_scope$archive_members
h5_archive <- file.path(archive_dir, "GSE159115_RAW.tar")
utils::untar(h5_archive, files = h5_member, exdir = work_dir)
h5_path <- file.path(work_dir, h5_member)
counts <- Seurat::Read10X_h5(h5_path, use.names = FALSE, unique.features = FALSE)
if (is.list(counts)) {
  if (!"Gene Expression" %in% names(counts)) stop("10x H5 contains multiple matrices but lacks Gene Expression")
  counts <- counts[["Gene Expression"]]
}
counts <- methods::as(counts, "dgCMatrix")
if (anyDuplicated(rownames(counts))) stop("GSE159115 feature IDs are duplicated")
rows[[length(rows) + 1L]] <- matrix_row("10x_H5", "GSE159115", h5_scope$library_id, counts, 0L, "Read10X_h5 use.names=FALSE unique.features=FALSE")
rm(counts)
unlink(h5_path)
gc(verbose = FALSE)

# GSE242299 AnnData CSR is directly a CSC matrix after transposition to genes x cells.
data <- rhdf5::h5read(h5ad_file, "/X/data")
indices <- rhdf5::h5read(h5ad_file, "/X/indices")
indptr <- rhdf5::h5read(h5ad_file, "/X/indptr")
features <- as.character(rhdf5::h5read(h5ad_file, "/var/_index"))
barcodes <- as.character(rhdf5::h5read(h5ad_file, "/obs/_index"))
if (length(indptr) != length(barcodes) + 1L) stop("H5AD CSR indptr/cell mismatch")
if (length(data) != length(indices) || tail(indptr, 1L) != length(data)) stop("H5AD CSR length mismatch")
counts <- methods::new(
  "dgCMatrix",
  i = as.integer(indices),
  p = as.integer(indptr),
  x = as.numeric(data),
  Dim = as.integer(c(length(features), length(barcodes))),
  Dimnames = list(features, barcodes)
)
methods::validObject(counts)
listing <- rhdf5::h5ls(h5ad_file, recursive = TRUE)
obs_fields <- listing[listing$group == "/obs" & listing$otype == "H5I_DATASET", , drop = FALSE]
rows[[length(rows) + 1L]] <- matrix_row("official_series_H5AD_X", "GSE242299", "ALL_31_LIBRARIES", counts, nrow(obs_fields), "AnnData CSR reinterpreted exactly as genes-by-cells CSC; no value transformation")
rm(data, indices, indptr, counts)
gc(verbose = FALSE)

result <- do.call(rbind, rows)
expected_cells <- c(GSE178481 = NA, GSE207493 = NA, GSE156632 = 1109, GSE159115 = 1704, GSE242299 = 50236)
for (idx in seq_len(nrow(result))) {
  expected <- expected_cells[[result$study[idx]]]
  if (!is.na(expected) && result$cells[idx] != expected) result$status[idx] <- "FAIL"
}
write_tsv(result, file.path(part_dir, "source_reader_preflight_v1.tsv"))

versions <- data.frame(
  component = c("R", required),
  version = c(R.version.string, vapply(required, function(package) as.character(utils::packageVersion(package)), character(1))),
  stringsAsFactors = FALSE
)
write_tsv(versions, file.path(part_dir, "source_reader_versions_v1.tsv"))
if (any(result$status != "PASS")) stop("One or more source readers failed")

writeLines(c(
  "status=PASS",
  paste0("reader_types=", nrow(result)),
  paste0("checks_pass=", sum(result$status == "PASS")),
  "normalization_performed=0",
  "integration_performed=0",
  "clustering_performed=0",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), file.path(part_dir, "PHASE5_SOURCE_READER_PREFLIGHT_PASS_v1"), useBytes = TRUE)
unlink(work_dir, recursive = TRUE, force = TRUE)
if (!file.rename(part_dir, final_dir)) stop("Cannot atomically finalize preflight directory")
cat("PASS final_dir=", final_dir, "\n", sep = "")
