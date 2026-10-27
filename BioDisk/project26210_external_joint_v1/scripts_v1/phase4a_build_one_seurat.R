#!/usr/bin/env Rscript

options(
    warn = 1,
    stringsAsFactors = FALSE,
    future.globals.maxSize = 8 * 1024^3
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10L) {
    stop(
        paste(
            "Expected 10 arguments:",
            "dataset sample_key gsm matrix features barcodes",
            "object_path summary_path done_path log_label"
        )
    )
}

dataset      <- args[[1]]
sample_key   <- args[[2]]
gsm          <- args[[3]]
matrix_path  <- args[[4]]
feature_path <- args[[5]]
barcode_path <- args[[6]]
object_path  <- args[[7]]
summary_path <- args[[8]]
done_path    <- args[[9]]
log_label    <- args[[10]]

required_packages <- c("Matrix", "Seurat")

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0L) {
    stop(
        "Missing R packages: ",
        paste(missing_packages, collapse = ", ")
    )
}

for (input_path in c(matrix_path, feature_path, barcode_path)) {
    if (!file.exists(input_path) || file.info(input_path)$size <= 0) {
        stop("Missing or empty input: ", input_path)
    }
}

dir.create(dirname(object_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(done_path), recursive = TRUE, showWarnings = FALSE)

message("START_SAMPLE=", log_label)
message("DATASET=", dataset)
message("SAMPLE_KEY=", sample_key)
message("GSM=", gsm)
message("MATRIX=", matrix_path)

counts <- Seurat::ReadMtx(
    mtx = matrix_path,
    cells = barcode_path,
    features = feature_path,
    cell.column = 1,
    feature.column = 2,
    unique.features = TRUE,
    strip.suffix = FALSE
)

counts <- methods::as(counts, "dgCMatrix")

feature_table <- utils::read.delim(
    gzfile(feature_path),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(feature_table) != nrow(counts)) {
    stop(
        "Feature row mismatch: feature_table=",
        nrow(feature_table),
        "; matrix=",
        nrow(counts)
    )
}

feature_type_filter_applied <- FALSE

if (ncol(feature_table) >= 3L) {
    feature_types <- as.character(feature_table[[3]])
    gene_expression_rows <- feature_types == "Gene Expression"

    if (any(gene_expression_rows) && !all(gene_expression_rows)) {
        counts <- counts[gene_expression_rows, , drop = FALSE]
        feature_type_filter_applied <- TRUE
    }
}

if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop("Expression matrix has zero rows or columns after feature filtering")
}

object <- Seurat::CreateSeuratObject(
    counts = counts,
    project = dataset,
    assay = "RNA",
    min.cells = 0,
    min.features = 0
)

object$dataset <- dataset
object$sample_key <- sample_key
object$gsm <- gsm

object$percent.mt <- Seurat::PercentageFeatureSet(
    object,
    pattern = "^MT-"
)

object$percent.ribo <- Seurat::PercentageFeatureSet(
    object,
    pattern = "^RP[SL]"
)

object$qc_pass_basic <- (
    object$nFeature_RNA >= 200 &
    object$nCount_RNA >= 500 &
    object$percent.mt <= 30
)

object@misc$phase4a_provenance <- list(
    dataset = dataset,
    sample_key = sample_key,
    gsm = gsm,
    matrix_path = matrix_path,
    feature_path = feature_path,
    barcode_path = barcode_path,
    feature_type_filter_applied = feature_type_filter_applied,
    qc_rule = paste(
        "nFeature_RNA >= 200;",
        "nCount_RNA >= 500;",
        "percent.mt <= 30;",
        "cells are flagged but not removed"
    ),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version = R.version.string,
    seurat_version = as.character(
        utils::packageVersion("Seurat")
    )
)

metadata <- object[[]]
n_cells <- ncol(object)
n_genes <- nrow(object)
n_pass <- sum(metadata$qc_pass_basic, na.rm = TRUE)

summary_table <- data.frame(
    dataset = dataset,
    sample_key = sample_key,
    gsm = gsm,
    genes = n_genes,
    cells = n_cells,
    nnz = Matrix::nnzero(counts),
    median_nCount = stats::median(metadata$nCount_RNA),
    median_nFeature = stats::median(metadata$nFeature_RNA),
    median_percent_mt = stats::median(metadata$percent.mt),
    median_percent_ribo = stats::median(metadata$percent.ribo),
    basic_qc_pass_cells = n_pass,
    basic_qc_pass_rate = n_pass / n_cells,
    feature_type_filter_applied = feature_type_filter_applied,
    stringsAsFactors = FALSE
)

temporary_object <- paste0(
    object_path,
    ".tmp_",
    Sys.getpid()
)

saveRDS(
    object,
    file = temporary_object,
    compress = "gzip"
)

if (!file.rename(temporary_object, object_path)) {
    stop("Failed to atomically move object to: ", object_path)
}

temporary_summary <- paste0(
    summary_path,
    ".tmp_",
    Sys.getpid()
)

utils::write.table(
    summary_table,
    file = temporary_summary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

if (!file.rename(temporary_summary, summary_path)) {
    stop("Failed to atomically move summary to: ", summary_path)
}

writeLines(
    c(
        paste0("status=PASS"),
        paste0("dataset=", dataset),
        paste0("sample_key=", sample_key),
        paste0("gsm=", gsm),
        paste0("cells=", n_cells),
        paste0("genes=", n_genes),
        paste0("basic_qc_pass_cells=", n_pass),
        paste0("completed_at=", format(
            Sys.time(),
            "%Y-%m-%dT%H:%M:%S%z"
        ))
    ),
    con = done_path
)

message("CELLS=", n_cells)
message("GENES=", n_genes)
message("BASIC_QC_PASS=", n_pass)
message("OBJECT=", object_path)
message("END_SAMPLE=", log_label)

rm(object, counts, metadata, feature_table)
invisible(gc())
