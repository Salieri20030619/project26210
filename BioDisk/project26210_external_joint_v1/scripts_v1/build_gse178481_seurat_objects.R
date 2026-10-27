#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    warn = 1
)

suppressPackageStartupMessages({
    library(Matrix)
    library(SeuratObject)
})

root_dir <- paste0(
    "/mnt/BioDisk/wangmingyuan/",
    "project26210_external_joint_v1"
)

phase_dir <- file.path(
    root_dir,
    "phase4c_gse178481_seurat_v1"
)

input_root <- file.path(
    phase_dir,
    "sparse_inputs"
)

conversion_manifest_file <- file.path(
    phase_dir,
    "sparse_conversion_manifest.tsv"
)

output_root <- file.path(
    phase_dir,
    "seurat_objects"
)

qc_dir <- file.path(
    phase_dir,
    "qc"
)

object_manifest_file <- file.path(
    qc_dir,
    "gse178481_seurat_object_manifest.tsv"
)

session_info_file <- file.path(
    qc_dir,
    "gse178481_seurat_session_info.txt"
)

completion_file <- file.path(
    phase_dir,
    ".seurat_build_complete"
)

dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    qc_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

stop_if_exists <- function(path, label) {
    if (file.exists(path)) {
        stop(
            label,
            " already exists; refusing to overwrite: ",
            path
        )
    }
}

stop_if_exists(
    object_manifest_file,
    "Object manifest"
)

stop_if_exists(
    completion_file,
    "Completion marker"
)

if (!file.exists(conversion_manifest_file)) {
    stop(
        "Conversion manifest does not exist: ",
        conversion_manifest_file
    )
}

conversion_manifest <- read.delim(
    conversion_manifest_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE
)

required_manifest_columns <- c(
    "sample_id",
    "input_file",
    "genes",
    "cells",
    "nonzero_values",
    "total_umi",
    "output_dir",
    "status"
)

missing_manifest_columns <- setdiff(
    required_manifest_columns,
    colnames(conversion_manifest)
)

if (length(missing_manifest_columns) > 0L) {
    stop(
        "Conversion manifest missing columns: ",
        paste(missing_manifest_columns, collapse = ",")
    )
}

if (nrow(conversion_manifest) != 26L) {
    stop(
        "Expected 26 conversion records, observed ",
        nrow(conversion_manifest)
    )
}

if (any(conversion_manifest$status != "PASS")) {
    failed_samples <- conversion_manifest$sample_id[
        conversion_manifest$status != "PASS"
    ]

    stop(
        "Non-PASS conversion records found: ",
        paste(failed_samples, collapse = ",")
    )
}

if (anyDuplicated(conversion_manifest$sample_id)) {
    stop("Duplicated sample_id values in conversion manifest")
}

read_gzip_lines <- function(path) {
    if (!file.exists(path)) {
        stop("Required file does not exist: ", path)
    }

    connection <- gzfile(
        path,
        open = "rt",
        encoding = "UTF-8"
    )

    on.exit(close(connection), add = TRUE)

    readLines(
        connection,
        warn = FALSE
    )
}

read_matrix_market_gzip <- function(path) {
    if (!file.exists(path)) {
        stop("Matrix file does not exist: ", path)
    }

    connection <- gzfile(
        path,
        open = "rt"
    )

    on.exit(close(connection), add = TRUE)

    matrix_object <- Matrix::readMM(connection)

    matrix_object <- methods::as(
        matrix_object,
        "CsparseMatrix"
    )

    if (!inherits(matrix_object, "dgCMatrix")) {
        matrix_object <- methods::as(
            matrix_object,
            "dgCMatrix"
        )
    }

    matrix_object
}

read_sample_metadata <- function(path) {
    if (!file.exists(path)) {
        stop("Sample metadata does not exist: ", path)
    }

    metadata <- read.delim(
        path,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE
    )

    required_columns <- c(
        "sample_id",
        "gsm",
        "patient_id",
        "tissue",
        "region",
        "source_dataset",
        "input_file"
    )

    missing_columns <- setdiff(
        required_columns,
        colnames(metadata)
    )

    if (length(missing_columns) > 0L) {
        stop(
            "Sample metadata missing columns: ",
            paste(missing_columns, collapse = ",")
        )
    }

    if (nrow(metadata) != 1L) {
        stop(
            "Sample metadata must contain exactly one row: ",
            path
        )
    }

    metadata
}

manifest_rows <- vector(
    mode = "list",
    length = nrow(conversion_manifest)
)

all_output_cell_names <- character(0)

build_started <- Sys.time()

for (index in seq_len(nrow(conversion_manifest))) {
    record <- conversion_manifest[index, , drop = FALSE]

    sample_id <- as.character(record$sample_id)
    sparse_dir <- as.character(record$output_dir)

    matrix_file <- file.path(
        sparse_dir,
        "matrix.mtx.gz"
    )

    features_file <- file.path(
        sparse_dir,
        "features.tsv.gz"
    )

    barcodes_file <- file.path(
        sparse_dir,
        "barcodes.tsv.gz"
    )

    sample_metadata_file <- file.path(
        sparse_dir,
        "sample_metadata.tsv"
    )

    conversion_complete_file <- file.path(
        sparse_dir,
        ".conversion_complete"
    )

    output_file <- file.path(
        output_root,
        paste0(sample_id, ".rds")
    )

    if (file.exists(output_file)) {
        stop(
            "Output RDS already exists; refusing to overwrite: ",
            output_file
        )
    }

    if (!file.exists(conversion_complete_file)) {
        stop(
            "Missing conversion completion marker: ",
            conversion_complete_file
        )
    }

    message(
        sprintf(
            "[%d/%d] START %s",
            index,
            nrow(conversion_manifest),
            sample_id
        )
    )

    sample_started <- Sys.time()

    features <- read_gzip_lines(features_file)
    barcodes <- read_gzip_lines(barcodes_file)
    sample_metadata <- read_sample_metadata(
        sample_metadata_file
    )

    if (anyNA(features) || any(features == "")) {
        stop(
            "Missing or empty feature names: ",
            sample_id
        )
    }

    if (anyDuplicated(features)) {
        duplicated_features <- unique(
            features[duplicated(features)]
        )

        stop(
            "Duplicated feature names in ",
            sample_id,
            ": ",
            paste(
                head(duplicated_features, 20L),
                collapse = ","
            )
        )
    }

    if (anyNA(barcodes) || any(barcodes == "")) {
        stop(
            "Missing or empty barcodes: ",
            sample_id
        )
    }

    if (anyDuplicated(barcodes)) {
        stop(
            "Duplicated barcodes within sample: ",
            sample_id
        )
    }

    expected_genes <- as.numeric(record$genes)
    expected_cells <- as.numeric(record$cells)
    expected_nnz <- as.numeric(record$nonzero_values)
    expected_umi <- as.numeric(record$total_umi)

    if (length(features) != expected_genes) {
        stop(
            "Feature count mismatch for ",
            sample_id,
            ": observed=",
            length(features),
            ", expected=",
            expected_genes
        )
    }

    if (length(barcodes) != expected_cells) {
        stop(
            "Barcode count mismatch for ",
            sample_id,
            ": observed=",
            length(barcodes),
            ", expected=",
            expected_cells
        )
    }

    counts <- read_matrix_market_gzip(matrix_file)

    rownames(counts) <- features
    colnames(counts) <- barcodes

    observed_genes <- nrow(counts)
    observed_cells <- ncol(counts)
    observed_nnz <- length(counts@x)
    observed_umi <- sum(counts@x)

    if (observed_genes != expected_genes) {
        stop(
            "Matrix gene count mismatch for ",
            sample_id
        )
    }

    if (observed_cells != expected_cells) {
        stop(
            "Matrix cell count mismatch for ",
            sample_id
        )
    }

    if (observed_nnz != expected_nnz) {
        stop(
            "Matrix nonzero count mismatch for ",
            sample_id,
            ": observed=",
            observed_nnz,
            ", expected=",
            expected_nnz
        )
    }

    if (observed_umi != expected_umi) {
        stop(
            "Matrix total UMI mismatch for ",
            sample_id,
            ": observed=",
            observed_umi,
            ", expected=",
            expected_umi
        )
    }

    if (any(counts@x < 0)) {
        stop(
            "Negative matrix value found in ",
            sample_id
        )
    }

    globally_unique_cell_names <- paste(
        sample_id,
        barcodes,
        sep = "__"
    )

    if (anyDuplicated(globally_unique_cell_names)) {
        stop(
            "Duplicated globally prefixed cell names within ",
            sample_id
        )
    }

    overlap <- intersect(
        globally_unique_cell_names,
        all_output_cell_names
    )

    if (length(overlap) > 0L) {
        stop(
            "Cell name collision across samples: ",
            sample_id
        )
    }

    colnames(counts) <- globally_unique_cell_names

    cell_metadata <- data.frame(
        sample_id = rep(
            sample_metadata$sample_id,
            observed_cells
        ),
        gsm = rep(
            sample_metadata$gsm,
            observed_cells
        ),
        patient_id = rep(
            sample_metadata$patient_id,
            observed_cells
        ),
        tissue = rep(
            sample_metadata$tissue,
            observed_cells
        ),
        region = rep(
            sample_metadata$region,
            observed_cells
        ),
        source_dataset = rep(
            sample_metadata$source_dataset,
            observed_cells
        ),
        input_file = rep(
            sample_metadata$input_file,
            observed_cells
        ),
        original_barcode = barcodes,
        stringsAsFactors = FALSE,
        row.names = globally_unique_cell_names
    )

    object <- CreateSeuratObject(
        counts = counts,
        assay = "RNA",
        project = "GSE178481",
        meta.data = cell_metadata,
        min.cells = 0,
        min.features = 0
    )

    DefaultAssay(object) <- "RNA"

    if (nrow(object) != expected_genes) {
        stop(
            "Seurat feature count mismatch for ",
            sample_id
        )
    }

    if (ncol(object) != expected_cells) {
        stop(
            "Seurat cell count mismatch for ",
            sample_id
        )
    }

    if (!identical(
        colnames(object),
        globally_unique_cell_names
    )) {
        stop(
            "Seurat cell names changed unexpectedly for ",
            sample_id
        )
    }

    required_metadata_columns <- c(
        "orig.ident",
        "nCount_RNA",
        "nFeature_RNA",
        "sample_id",
        "gsm",
        "patient_id",
        "tissue",
        "region",
        "source_dataset",
        "input_file",
        "original_barcode"
    )

    missing_object_metadata <- setdiff(
        required_metadata_columns,
        colnames(object[[]])
    )

    if (length(missing_object_metadata) > 0L) {
        stop(
            "Seurat metadata missing columns for ",
            sample_id,
            ": ",
            paste(
                missing_object_metadata,
                collapse = ","
            )
        )
    }

    object_umi <- sum(object$nCount_RNA)

    if (object_umi != expected_umi) {
        stop(
            "Seurat nCount_RNA sum mismatch for ",
            sample_id,
            ": observed=",
            object_umi,
            ", expected=",
            expected_umi
        )
    }

    saveRDS(
        object,
        file = output_file,
        compress = "gzip"
    )

    if (!file.exists(output_file)) {
        stop(
            "RDS was not created: ",
            output_file
        )
    }

    verification_object <- readRDS(output_file)

    verified_genes <- nrow(verification_object)
    verified_cells <- ncol(verification_object)
    verified_umi <- sum(
        verification_object$nCount_RNA
    )

    if (
        verified_genes != expected_genes ||
        verified_cells != expected_cells ||
        verified_umi != expected_umi
    ) {
        stop(
            "Post-save RDS verification failed for ",
            sample_id
        )
    }

    elapsed_seconds <- as.numeric(
        difftime(
            Sys.time(),
            sample_started,
            units = "secs"
        )
    )

    manifest_rows[[index]] <- data.frame(
        sample_id = sample_id,
        gsm = as.character(sample_metadata$gsm),
        patient_id = as.character(
            sample_metadata$patient_id
        ),
        tissue = as.character(sample_metadata$tissue),
        region = as.character(sample_metadata$region),
        genes = verified_genes,
        cells = verified_cells,
        nonzero_values = observed_nnz,
        total_umi = verified_umi,
        matrix_class = class(counts)[1],
        assay = DefaultAssay(verification_object),
        cell_names_unique_within_sample = (
            anyDuplicated(
                colnames(verification_object)
            ) == 0L
        ),
        output_file = output_file,
        output_bytes = file.info(output_file)$size,
        elapsed_seconds = sprintf(
            "%.3f",
            elapsed_seconds
        ),
        status = "PASS",
        reason = "OK",
        stringsAsFactors = FALSE
    )

    all_output_cell_names <- c(
        all_output_cell_names,
        colnames(verification_object)
    )

    rm(
        counts,
        object,
        verification_object,
        cell_metadata
    )

    invisible(gc())

    message(
        sprintf(
            paste0(
                "[%d/%d] PASS ",
                "cells=%d nnz=%.0f umi=%.0f seconds=%.1f"
            ),
            index,
            nrow(conversion_manifest),
            verified_cells,
            observed_nnz,
            verified_umi,
            elapsed_seconds
        )
    )
}

object_manifest <- do.call(
    rbind,
    manifest_rows
)

if (nrow(object_manifest) != 26L) {
    stop(
        "Final object manifest does not contain 26 records"
    )
}

if (any(object_manifest$status != "PASS")) {
    stop("Final object manifest contains non-PASS records")
}

if (anyDuplicated(all_output_cell_names)) {
    stop("Duplicated cell names across final Seurat objects")
}

manifest_temp <- paste0(
    object_manifest_file,
    ".tmp"
)

write.table(
    object_manifest,
    file = manifest_temp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

if (!file.rename(
    from = manifest_temp,
    to = object_manifest_file
)) {
    stop(
        "Failed to finalize object manifest: ",
        object_manifest_file
    )
}

session_connection <- file(
    session_info_file,
    open = "wt"
)

sink(session_connection)

cat("===== SESSION INFO =====\n")
print(sessionInfo())

sink()
close(session_connection)

total_elapsed_seconds <- as.numeric(
    difftime(
        Sys.time(),
        build_started,
        units = "secs"
    )
)

total_cells <- sum(object_manifest$cells)
total_nnz <- sum(object_manifest$nonzero_values)
total_umi <- sum(object_manifest$total_umi)

completion_temp <- paste0(
    completion_file,
    ".tmp"
)

writeLines(
    c(
        "status=PASS",
        paste0("objects=", nrow(object_manifest)),
        paste0("cells=", total_cells),
        paste0("nonzero_values=", total_nnz),
        paste0("total_umi=", total_umi),
        paste0(
            "unique_cell_names=",
            length(all_output_cell_names)
        ),
        paste0(
            "elapsed_seconds=",
            sprintf("%.3f", total_elapsed_seconds)
        ),
        paste0(
            "manifest=",
            object_manifest_file
        )
    ),
    con = completion_temp
)

if (!file.rename(
    from = completion_temp,
    to = completion_file
)) {
    stop(
        "Failed to finalize completion marker: ",
        completion_file
    )
}

cat("===== COMPLETION =====\n")
cat("status=PASS\n")
cat("objects=", nrow(object_manifest), "\n", sep = "")
cat("cells=", total_cells, "\n", sep = "")
cat("nonzero_values=", total_nnz, "\n", sep = "")
cat("total_umi=", total_umi, "\n", sep = "")
cat(
    "unique_cell_names=",
    length(all_output_cell_names),
    "\n",
    sep = ""
)
cat(
    "manifest=",
    object_manifest_file,
    "\n",
    sep = ""
)
