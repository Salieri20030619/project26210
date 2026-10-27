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

conversion_manifest_file <- file.path(
    phase_dir,
    "sparse_conversion_manifest.tsv"
)

object_manifest_file <- file.path(
    phase_dir,
    "qc",
    "gse178481_seurat_object_manifest.tsv"
)

audit_file <- file.path(
    phase_dir,
    "qc",
    "gse178481_independent_rds_audit.tsv"
)

rename_file <- file.path(
    phase_dir,
    "qc",
    "gse178481_feature_rename_map.tsv"
)

summary_file <- file.path(
    phase_dir,
    "qc",
    "gse178481_independent_rds_audit_summary.txt"
)

stop_if_exists <- function(path) {
    if (file.exists(path)) {
        stop("Refusing to overwrite existing output: ", path)
    }
}

stop_if_exists(audit_file)
stop_if_exists(rename_file)
stop_if_exists(summary_file)

conversion_manifest <- read.delim(
    conversion_manifest_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

object_manifest <- read.delim(
    object_manifest_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

if (nrow(conversion_manifest) != 26L) {
    stop(
        "Expected 26 conversion records, observed ",
        nrow(conversion_manifest)
    )
}

if (nrow(object_manifest) != 26L) {
    stop(
        "Expected 26 object records, observed ",
        nrow(object_manifest)
    )
}

if (anyDuplicated(conversion_manifest$sample_id)) {
    stop("Duplicated sample_id in conversion manifest")
}

if (anyDuplicated(object_manifest$sample_id)) {
    stop("Duplicated sample_id in object manifest")
}

audit_rows <- vector("list", nrow(object_manifest))
rename_rows <- list()
rename_index <- 0L
all_cell_names <- character(0)
reference_original_features <- NULL
reference_rds_features <- NULL

for (index in seq_len(nrow(object_manifest))) {
    record <- object_manifest[index, , drop = FALSE]
    sample_id <- as.character(record$sample_id)
    rds_file <- as.character(record$output_file)

    conversion_record <- conversion_manifest[
        conversion_manifest$sample_id == sample_id,
        ,
        drop = FALSE
    ]

    if (nrow(conversion_record) != 1L) {
        stop(
            "Cannot resolve one conversion record for ",
            sample_id
        )
    }

    sparse_dir <- as.character(
        conversion_record$output_dir
    )

    features_file <- file.path(
        sparse_dir,
        "features.tsv.gz"
    )

    if (!file.exists(rds_file)) {
        stop("Missing RDS: ", rds_file)
    }

    if (!file.exists(features_file)) {
        stop("Missing feature file: ", features_file)
    }

    feature_connection <- gzfile(
        features_file,
        open = "rt",
        encoding = "UTF-8"
    )

    original_features <- readLines(
        feature_connection,
        warn = FALSE
    )

    close(feature_connection)

    object <- readRDS(rds_file)

    if (!inherits(object, "Seurat")) {
        stop("Object is not a Seurat object: ", sample_id)
    }

    if (!"RNA" %in% names(object@assays)) {
        stop("RNA assay is absent: ", sample_id)
    }

    DefaultAssay(object) <- "RNA"

    counts <- LayerData(
        object = object,
        assay = "RNA",
        layer = "counts"
    )

    rds_features <- rownames(counts)
    rds_cells <- colnames(counts)

    if (length(original_features) != length(rds_features)) {
        stop(
            "Feature length mismatch for ",
            sample_id
        )
    }

    expected_seurat_features <- gsub(
        "_",
        "-",
        original_features,
        fixed = TRUE
    )

    renamed <- original_features != rds_features

    direct_replacement_match <- (
        expected_seurat_features == rds_features
    )

    unexpected_feature_changes <- sum(
        !direct_replacement_match
    )

    original_underscore_count <- sum(
        grepl("_", original_features, fixed = TRUE)
    )

    renamed_feature_count <- sum(renamed)

    transformed_duplicates <- sum(
        duplicated(expected_seurat_features) |
        duplicated(
            expected_seurat_features,
            fromLast = TRUE
        )
    )

    original_duplicate_count <- sum(
        duplicated(original_features)
    )

    rds_duplicate_count <- sum(
        duplicated(rds_features)
    )

    if (renamed_feature_count > 0L) {
        changed_indices <- which(renamed)

        rename_index <- rename_index + 1L

        rename_rows[[rename_index]] <- data.frame(
            sample_id = sample_id,
            feature_index = changed_indices,
            original_feature = original_features[
                changed_indices
            ],
            rds_feature = rds_features[
                changed_indices
            ],
            expected_underscore_replacement = (
                expected_seurat_features[
                    changed_indices
                ] == rds_features[
                    changed_indices
                ]
            ),
            stringsAsFactors = FALSE
        )
    }

    observed_genes <- nrow(counts)
    observed_cells <- ncol(counts)
    observed_nnz <- length(counts@x)
    observed_umi <- sum(counts@x)

    expected_genes <- as.numeric(
        conversion_record$genes
    )

    expected_cells <- as.numeric(
        conversion_record$cells
    )

    expected_nnz <- as.numeric(
        conversion_record$nonzero_values
    )

    expected_umi <- as.numeric(
        conversion_record$total_umi
    )

    metadata <- object[[]]

    required_metadata <- c(
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

    missing_metadata <- setdiff(
        required_metadata,
        colnames(metadata)
    )

    sample_metadata_consistent <- (
        length(unique(metadata$sample_id)) == 1L &&
        identical(
            unique(metadata$sample_id),
            sample_id
        )
    )

    cell_name_prefix_ok <- all(
        startsWith(
            rds_cells,
            paste0(sample_id, "__")
        )
    )

    cells_unique_within_object <- (
        anyDuplicated(rds_cells) == 0L
    )

    overlap_count <- length(
        intersect(rds_cells, all_cell_names)
    )

    if (is.null(reference_original_features)) {
        reference_original_features <- original_features
        reference_rds_features <- rds_features
        original_order_matches_reference <- TRUE
        rds_order_matches_reference <- TRUE
    } else {
        original_order_matches_reference <- identical(
            original_features,
            reference_original_features
        )

        rds_order_matches_reference <- identical(
            rds_features,
            reference_rds_features
        )
    }

    status <- if (
        observed_genes == expected_genes &&
        observed_cells == expected_cells &&
        observed_nnz == expected_nnz &&
        observed_umi == expected_umi &&
        length(missing_metadata) == 0L &&
        sample_metadata_consistent &&
        cell_name_prefix_ok &&
        cells_unique_within_object &&
        overlap_count == 0L &&
        original_duplicate_count == 0L &&
        rds_duplicate_count == 0L &&
        unexpected_feature_changes == 0L &&
        original_order_matches_reference &&
        rds_order_matches_reference
    ) {
        "PASS"
    } else {
        "FAIL"
    }

    audit_rows[[index]] <- data.frame(
        sample_id = sample_id,
        genes = observed_genes,
        cells = observed_cells,
        nonzero_values = observed_nnz,
        total_umi = observed_umi,
        matrix_class = class(counts)[1],
        original_underscore_features = original_underscore_count,
        renamed_features = renamed_feature_count,
        unexpected_feature_changes = unexpected_feature_changes,
        transformed_collision_features = transformed_duplicates,
        original_duplicate_features = original_duplicate_count,
        rds_duplicate_features = rds_duplicate_count,
        original_order_matches_reference = (
            original_order_matches_reference
        ),
        rds_order_matches_reference = (
            rds_order_matches_reference
        ),
        missing_metadata_columns = if (
            length(missing_metadata) == 0L
        ) {
            "NONE"
        } else {
            paste(missing_metadata, collapse = ",")
        },
        sample_metadata_consistent = (
            sample_metadata_consistent
        ),
        cell_name_prefix_ok = cell_name_prefix_ok,
        cells_unique_within_object = (
            cells_unique_within_object
        ),
        cross_object_cell_overlap = overlap_count,
        status = status,
        stringsAsFactors = FALSE
    )

    all_cell_names <- c(
        all_cell_names,
        rds_cells
    )

    message(
        sprintf(
            paste0(
                "[%d/26] %s ",
                "cells=%d umi=%.0f renamed=%d ",
                "unexpected=%d collisions=%d status=%s"
            ),
            index,
            sample_id,
            observed_cells,
            observed_umi,
            renamed_feature_count,
            unexpected_feature_changes,
            transformed_duplicates,
            status
        )
    )

    rm(object, counts, metadata)
    invisible(gc())
}

audit_table <- do.call(
    rbind,
    audit_rows
)

if (length(rename_rows) > 0L) {
    rename_table <- do.call(
        rbind,
        rename_rows
    )
} else {
    rename_table <- data.frame(
        sample_id = character(0),
        feature_index = integer(0),
        original_feature = character(0),
        rds_feature = character(0),
        expected_underscore_replacement = logical(0),
        stringsAsFactors = FALSE
    )
}

write.table(
    audit_table,
    file = audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    rename_table,
    file = rename_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

summary_lines <- c(
    "===== INDEPENDENT RDS AUDIT =====",
    paste0("objects=", nrow(audit_table)),
    paste0("cells=", sum(audit_table$cells)),
    paste0(
        "nonzero_values=",
        sum(audit_table$nonzero_values)
    ),
    paste0(
        "total_umi=",
        format(
            sum(audit_table$total_umi),
            scientific = FALSE,
            trim = TRUE
        )
    ),
    paste0(
        "unique_cell_names=",
        length(unique(all_cell_names))
    ),
    paste0(
        "reference_original_features=",
        length(reference_original_features)
    ),
    paste0(
        "reference_rds_features=",
        length(reference_rds_features)
    ),
    paste0(
        "renamed_features_per_object_min=",
        min(audit_table$renamed_features)
    ),
    paste0(
        "renamed_features_per_object_max=",
        max(audit_table$renamed_features)
    ),
    paste0(
        "unexpected_feature_changes_total=",
        sum(audit_table$unexpected_feature_changes)
    ),
    paste0(
        "transformed_collision_features_max=",
        max(audit_table$transformed_collision_features)
    ),
    paste0(
        "rds_duplicate_features_max=",
        max(audit_table$rds_duplicate_features)
    ),
    paste0(
        "cross_object_cell_overlap_total=",
        sum(audit_table$cross_object_cell_overlap)
    ),
    paste0(
        "status_PASS=",
        sum(audit_table$status == "PASS")
    ),
    paste0(
        "status_FAIL=",
        sum(audit_table$status == "FAIL")
    ),
    paste0("audit_file=", audit_file),
    paste0("rename_file=", rename_file)
)

writeLines(
    summary_lines,
    con = summary_file
)

cat(
    paste(summary_lines, collapse = "\n"),
    "\n",
    sep = ""
)
