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

phase4a_dir <- file.path(
    root_dir,
    "phase4a_10x_seurat_v1",
    "objects"
)

phase4c_dir <- file.path(
    root_dir,
    "phase4c_gse178481_seurat_v1",
    "seurat_objects"
)

output_dir <- file.path(
    root_dir,
    "phase4d_95_object_compatibility_audit_v1",
    "qc"
)

object_audit_file <- file.path(
    output_dir,
    "object_level_audit.tsv"
)

feature_summary_file <- file.path(
    output_dir,
    "feature_set_summary.tsv"
)

feature_pair_file <- file.path(
    output_dir,
    "feature_pairwise_overlap.tsv"
)

metadata_inventory_file <- file.path(
    output_dir,
    "metadata_column_inventory.tsv"
)

dataset_summary_file <- file.path(
    output_dir,
    "dataset_level_summary.tsv"
)

common_feature_file <- file.path(
    output_dir,
    "features_common_to_all_95_objects.txt"
)

canonical_common_feature_file <- file.path(
    output_dir,
    "canonical_features_common_to_all_95_objects.txt"
)

summary_file <- file.path(
    output_dir,
    "compatibility_audit_summary.txt"
)

session_info_file <- file.path(
    output_dir,
    "sessionInfo.txt"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

outputs <- c(
    object_audit_file,
    feature_summary_file,
    feature_pair_file,
    metadata_inventory_file,
    dataset_summary_file,
    common_feature_file,
    canonical_common_feature_file,
    summary_file,
    session_info_file
)

existing_outputs <- outputs[file.exists(outputs)]

if (length(existing_outputs) > 0L) {
    stop(
        "Refusing to overwrite existing outputs: ",
        paste(existing_outputs, collapse = ", ")
    )
}

phase4a_files <- sort(
    list.files(
        phase4a_dir,
        pattern = "\\.rds$",
        full.names = TRUE
    )
)

phase4c_files <- sort(
    list.files(
        phase4c_dir,
        pattern = "\\.rds$",
        full.names = TRUE
    )
)

if (length(phase4a_files) != 69L) {
    stop(
        "Expected 69 phase4a RDS files, observed ",
        length(phase4a_files)
    )
}

if (length(phase4c_files) != 26L) {
    stop(
        "Expected 26 phase4c RDS files, observed ",
        length(phase4c_files)
    )
}

input_table <- rbind(
    data.frame(
        input_group = "phase4a_existing_10x",
        file = phase4a_files,
        stringsAsFactors = FALSE
    ),
    data.frame(
        input_group = "phase4c_gse178481",
        file = phase4c_files,
        stringsAsFactors = FALSE
    )
)

input_table$object_name <- sub(
    "\\.rds$",
    "",
    basename(input_table$file),
    ignore.case = TRUE
)

input_table$dataset <- ifelse(
    grepl("^GSE[0-9]+", input_table$object_name),
    sub(
        "^(GSE[0-9]+).*",
        "\\1",
        input_table$object_name
    ),
    ifelse(
        input_table$input_group == "phase4c_gse178481",
        "GSE178481",
        "UNRESOLVED"
    )
)

if (nrow(input_table) != 95L) {
    stop("Final input table does not contain 95 records")
}

if (anyDuplicated(input_table$file)) {
    stop("Duplicated file paths in input table")
}

get_counts_layer <- function(object, assay_name) {
    assay_object <- object[[assay_name]]

    available_layers <- tryCatch(
        Layers(assay_object),
        error = function(error_condition) {
            character(0)
        }
    )

    if ("counts" %in% available_layers) {
        counts <- LayerData(
            object = object,
            assay = assay_name,
            layer = "counts"
        )

        return(
            list(
                matrix = counts,
                source = "LayerData:counts"
            )
        )
    }

    counts <- tryCatch(
        GetAssayData(
            object = object,
            assay = assay_name,
            slot = "counts"
        ),
        error = function(error_condition) {
            NULL
        }
    )

    if (!is.null(counts)) {
        return(
            list(
                matrix = counts,
                source = "GetAssayData:counts"
            )
        )
    }

    stop(
        "No readable counts layer in assay ",
        assay_name
    )
}

canonicalize_features <- function(features) {
    gsub(
        "_",
        "-",
        features,
        fixed = TRUE
    )
}

audit_rows <- vector(
    "list",
    nrow(input_table)
)

metadata_rows <- list()
metadata_row_index <- 0L

feature_list <- vector(
    "list",
    nrow(input_table)
)

canonical_feature_list <- vector(
    "list",
    nrow(input_table)
)

all_cell_names <- character(0)

for (index in seq_len(nrow(input_table))) {
    record <- input_table[index, , drop = FALSE]

    object_name <- as.character(record$object_name)
    object_file <- as.character(record$file)
    input_group <- as.character(record$input_group)
    dataset <- as.character(record$dataset)

    message(
        sprintf(
            "[%d/%d] START %s",
            index,
            nrow(input_table),
            object_name
        )
    )

    file_bytes <- file.info(object_file)$size

    audit_result <- tryCatch({
        object <- readRDS(object_file)

        if (!inherits(object, "Seurat")) {
            stop("RDS does not contain a Seurat object")
        }

        assay_names <- names(object@assays)

        if (length(assay_names) == 0L) {
            stop("No assay is present")
        }

        default_assay <- DefaultAssay(object)

        assay_name <- if ("RNA" %in% assay_names) {
            "RNA"
        } else {
            default_assay
        }

        counts_result <- get_counts_layer(
            object = object,
            assay_name = assay_name
        )

        counts <- counts_result$matrix
        counts_source <- counts_result$source

        features <- rownames(counts)
        cells <- colnames(counts)

        if (is.null(features) || anyNA(features)) {
            stop("Feature names are missing")
        }

        if (is.null(cells) || anyNA(cells)) {
            stop("Cell names are missing")
        }

        if (any(features == "")) {
            stop("Empty feature names are present")
        }

        if (any(cells == "")) {
            stop("Empty cell names are present")
        }

        canonical_features <- canonicalize_features(
            features
        )

        duplicate_features <- sum(
            duplicated(features)
        )

        duplicate_canonical_features <- sum(
            duplicated(canonical_features)
        )

        underscore_features <- sum(
            grepl("_", features, fixed = TRUE)
        )

        hyphen_features <- sum(
            grepl("-", features, fixed = TRUE)
        )

        all_upper_symbol_like <- sum(
            grepl(
                "^[A-Z0-9][A-Z0-9.-]*$",
                features
            )
        )

        ensembl_like <- sum(
            grepl(
                "^ENS[A-Z]*G[0-9]+",
                features
            )
        )

        mitochondrial_upper <- sum(
            grepl("^MT-", features)
        )

        mitochondrial_title <- sum(
            grepl("^mt-", features)
        )

        ribosomal_upper <- sum(
            grepl("^RP[SL][0-9]", features)
        )

        ribosomal_title <- sum(
            grepl("^Rp[sl][0-9]", features)
        )

        cells_unique <- (
            anyDuplicated(cells) == 0L
        )

        cross_object_overlap <- length(
            intersect(cells, all_cell_names)
        )

        metadata <- object[[]]

        metadata_cell_alignment <- identical(
            rownames(metadata),
            cells
        )

        ncount_available <- (
            "nCount_RNA" %in% colnames(metadata)
        )

        nfeature_available <- (
            "nFeature_RNA" %in% colnames(metadata)
        )

        observed_umi <- if (
            inherits(counts, "sparseMatrix")
        ) {
            sum(counts@x)
        } else {
            sum(counts)
        }

        observed_nnz <- if (
            inherits(counts, "sparseMatrix")
        ) {
            length(counts@x)
        } else {
            sum(counts != 0)
        }

        ncount_sum <- if (ncount_available) {
            sum(metadata$nCount_RNA)
        } else {
            NA_real_
        }

        ncount_matches_counts <- if (
            ncount_available
        ) {
            isTRUE(all.equal(
                as.numeric(ncount_sum),
                as.numeric(observed_umi),
                tolerance = 0
            ))
        } else {
            NA
        }

        feature_list[[index]] <- features
        canonical_feature_list[[index]] <- (
            canonical_features
        )

        if (ncol(metadata) > 0L) {
            for (
                metadata_column in colnames(metadata)
            ) {
                metadata_row_index <- (
                    metadata_row_index + 1L
                )

                metadata_value <- metadata[[metadata_column]]

                unique_values <- unique(
                    as.character(metadata_value)
                )

                unique_values <- unique_values[
                    !is.na(unique_values)
                ]

                metadata_rows[[metadata_row_index]] <- data.frame(
                    input_group = input_group,
                    dataset = dataset,
                    object_name = object_name,
                    metadata_column = metadata_column,
                    metadata_class = paste(
                        class(metadata_value),
                        collapse = ";"
                    ),
                    non_missing_values = sum(
                        !is.na(metadata_value)
                    ),
                    unique_values = length(
                        unique_values
                    ),
                    example_values = paste(
                        head(unique_values, 5L),
                        collapse = " | "
                    ),
                    stringsAsFactors = FALSE
                )
            }
        }

        all_cell_names <- c(
            all_cell_names,
            cells
        )

        status <- if (
            duplicate_features == 0L &&
            cells_unique &&
            cross_object_overlap == 0L &&
            metadata_cell_alignment &&
            duplicate_canonical_features == 0L
        ) {
            "PASS"
        } else {
            "FAIL"
        }

        result <- data.frame(
            input_group = input_group,
            dataset = dataset,
            object_name = object_name,
            file = object_file,
            file_bytes = file_bytes,
            object_class = class(object)[1],
            seurat_object_version = as.character(
                object@version
            ),
            assays = paste(
                assay_names,
                collapse = ";"
            ),
            default_assay = default_assay,
            audited_assay = assay_name,
            counts_source = counts_source,
            matrix_class = class(counts)[1],
            genes = nrow(counts),
            cells = ncol(counts),
            nonzero_values = observed_nnz,
            total_umi = observed_umi,
            metadata_columns = ncol(metadata),
            metadata_cell_alignment = (
                metadata_cell_alignment
            ),
            nCount_RNA_available = ncount_available,
            nCount_RNA_sum = ncount_sum,
            nCount_RNA_matches_counts = (
                ncount_matches_counts
            ),
            duplicate_features = duplicate_features,
            underscore_features = underscore_features,
            hyphen_features = hyphen_features,
            canonical_duplicate_features = (
                duplicate_canonical_features
            ),
            uppercase_symbol_like_features = (
                all_upper_symbol_like
            ),
            ensembl_like_features = ensembl_like,
            mitochondrial_upper_features = (
                mitochondrial_upper
            ),
            mitochondrial_title_features = (
                mitochondrial_title
            ),
            ribosomal_upper_features = (
                ribosomal_upper
            ),
            ribosomal_title_features = (
                ribosomal_title
            ),
            cells_unique_within_object = (
                cells_unique
            ),
            cross_object_cell_overlap = (
                cross_object_overlap
            ),
            status = status,
            reason = if (status == "PASS") {
                "OK"
            } else {
                "STRUCTURAL_OR_NAMING_FAILURE"
            },
            stringsAsFactors = FALSE
        )

        rm(
            object,
            counts,
            metadata
        )

        invisible(gc())

        result
    }, error = function(error_condition) {
        feature_list[[index]] <<- character(0)
        canonical_feature_list[[index]] <<- (
            character(0)
        )

        data.frame(
            input_group = input_group,
            dataset = dataset,
            object_name = object_name,
            file = object_file,
            file_bytes = file_bytes,
            object_class = NA_character_,
            seurat_object_version = NA_character_,
            assays = NA_character_,
            default_assay = NA_character_,
            audited_assay = NA_character_,
            counts_source = NA_character_,
            matrix_class = NA_character_,
            genes = NA_real_,
            cells = NA_real_,
            nonzero_values = NA_real_,
            total_umi = NA_real_,
            metadata_columns = NA_real_,
            metadata_cell_alignment = NA,
            nCount_RNA_available = NA,
            nCount_RNA_sum = NA_real_,
            nCount_RNA_matches_counts = NA,
            duplicate_features = NA_real_,
            underscore_features = NA_real_,
            hyphen_features = NA_real_,
            canonical_duplicate_features = NA_real_,
            uppercase_symbol_like_features = NA_real_,
            ensembl_like_features = NA_real_,
            mitochondrial_upper_features = NA_real_,
            mitochondrial_title_features = NA_real_,
            ribosomal_upper_features = NA_real_,
            ribosomal_title_features = NA_real_,
            cells_unique_within_object = NA,
            cross_object_cell_overlap = NA_real_,
            status = "ERROR",
            reason = conditionMessage(error_condition),
            stringsAsFactors = FALSE
        )
    })

    audit_rows[[index]] <- audit_result

    message(
        sprintf(
            paste0(
                "[%d/%d] %s genes=%s cells=%s ",
                "underscore=%s canonical_dup=%s status=%s"
            ),
            index,
            nrow(input_table),
            object_name,
            as.character(audit_result$genes),
            as.character(audit_result$cells),
            as.character(
                audit_result$underscore_features
            ),
            as.character(
                audit_result$canonical_duplicate_features
            ),
            audit_result$status
        )
    )
}

object_audit <- do.call(
    rbind,
    audit_rows
)

if (length(metadata_rows) > 0L) {
    metadata_inventory <- do.call(
        rbind,
        metadata_rows
    )
} else {
    metadata_inventory <- data.frame(
        input_group = character(0),
        dataset = character(0),
        object_name = character(0),
        metadata_column = character(0),
        metadata_class = character(0),
        non_missing_values = numeric(0),
        unique_values = numeric(0),
        example_values = character(0),
        stringsAsFactors = FALSE
    )
}

valid_feature_indices <- which(
    lengths(feature_list) > 0L
)

if (length(valid_feature_indices) != 95L) {
    common_features <- character(0)
    canonical_common_features <- character(0)
} else {
    common_features <- Reduce(
        intersect,
        feature_list
    )

    canonical_common_features <- Reduce(
        intersect,
        canonical_feature_list
    )
}

feature_keys <- vapply(
    feature_list,
    function(features) {
        paste(features, collapse = "\034")
    },
    character(1)
)

canonical_feature_keys <- vapply(
    canonical_feature_list,
    function(features) {
        paste(features, collapse = "\034")
    },
    character(1)
)

feature_set_id <- match(
    feature_keys,
    unique(feature_keys)
)

canonical_feature_set_id <- match(
    canonical_feature_keys,
    unique(canonical_feature_keys)
)

feature_summary <- data.frame(
    input_group = input_table$input_group,
    dataset = input_table$dataset,
    object_name = input_table$object_name,
    feature_count = lengths(feature_list),
    feature_set_id = feature_set_id,
    canonical_feature_set_id = (
        canonical_feature_set_id
    ),
    stringsAsFactors = FALSE
)

pair_rows <- list()
pair_index <- 0L

for (first_index in seq_len(nrow(input_table) - 1L)) {
    for (
        second_index in seq.int(
            first_index + 1L,
            nrow(input_table)
        )
    ) {
        first_features <- feature_list[[first_index]]

        second_features <- feature_list[[second_index]]

        first_canonical <- canonical_feature_list[[first_index]]

        second_canonical <- canonical_feature_list[[second_index]]

        raw_intersection <- length(
            intersect(
                first_features,
                second_features
            )
        )

        canonical_intersection <- length(
            intersect(
                first_canonical,
                second_canonical
            )
        )

        raw_union <- length(
            union(
                first_features,
                second_features
            )
        )

        canonical_union <- length(
            union(
                first_canonical,
                second_canonical
            )
        )

        pair_index <- pair_index + 1L

        pair_rows[[pair_index]] <- data.frame(
            object_a = input_table$object_name[
                first_index
            ],
            dataset_a = input_table$dataset[
                first_index
            ],
            group_a = input_table$input_group[
                first_index
            ],
            object_b = input_table$object_name[
                second_index
            ],
            dataset_b = input_table$dataset[
                second_index
            ],
            group_b = input_table$input_group[
                second_index
            ],
            features_a = length(first_features),
            features_b = length(second_features),
            raw_intersection = raw_intersection,
            raw_jaccard = if (raw_union > 0L) {
                raw_intersection / raw_union
            } else {
                NA_real_
            },
            canonical_intersection = (
                canonical_intersection
            ),
            canonical_jaccard = if (
                canonical_union > 0L
            ) {
                canonical_intersection /
                    canonical_union
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }
}

feature_pairwise <- do.call(
    rbind,
    pair_rows
)

dataset_levels <- sort(
    unique(object_audit$dataset)
)

dataset_rows <- vector(
    "list",
    length(dataset_levels)
)

for (index in seq_along(dataset_levels)) {
    dataset_name <- dataset_levels[index]

    subset_table <- object_audit[
        object_audit$dataset == dataset_name,
        ,
        drop = FALSE
    ]

    dataset_rows[[index]] <- data.frame(
        dataset = dataset_name,
        objects = nrow(subset_table),
        pass = sum(subset_table$status == "PASS"),
        fail = sum(subset_table$status == "FAIL"),
        error = sum(subset_table$status == "ERROR"),
        cells = sum(
            subset_table$cells,
            na.rm = TRUE
        ),
        genes_min = if (
            all(is.na(subset_table$genes))
        ) {
            NA_real_
        } else {
            min(
                subset_table$genes,
                na.rm = TRUE
            )
        },
        genes_max = if (
            all(is.na(subset_table$genes))
        ) {
            NA_real_
        } else {
            max(
                subset_table$genes,
                na.rm = TRUE
            )
        },
        feature_sets = length(unique(
            feature_summary$feature_set_id[
                feature_summary$dataset ==
                    dataset_name
            ]
        )),
        canonical_feature_sets = length(unique(
            feature_summary$canonical_feature_set_id[
                feature_summary$dataset ==
                    dataset_name
            ]
        )),
        stringsAsFactors = FALSE
    )
}

dataset_summary <- do.call(
    rbind,
    dataset_rows
)

write.table(
    object_audit,
    file = object_audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    feature_summary,
    file = feature_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    feature_pairwise,
    file = feature_pair_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    metadata_inventory,
    file = metadata_inventory_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    dataset_summary,
    file = dataset_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

writeLines(
    common_features,
    con = common_feature_file
)

writeLines(
    canonical_common_features,
    con = canonical_common_feature_file
)

phase4a_audit <- object_audit[
    object_audit$input_group ==
        "phase4a_existing_10x",
    ,
    drop = FALSE
]

phase4c_audit <- object_audit[
    object_audit$input_group ==
        "phase4c_gse178481",
    ,
    drop = FALSE
]

cross_group_pairs <- feature_pairwise[
    feature_pairwise$group_a !=
        feature_pairwise$group_b,
    ,
    drop = FALSE
]

summary_lines <- c(
    "===== 95 OBJECT COMPATIBILITY AUDIT =====",
    paste0("objects_total=", nrow(object_audit)),
    paste0("objects_phase4a=", nrow(phase4a_audit)),
    paste0("objects_phase4c=", nrow(phase4c_audit)),
    paste0(
        "status_PASS=",
        sum(object_audit$status == "PASS")
    ),
    paste0(
        "status_FAIL=",
        sum(object_audit$status == "FAIL")
    ),
    paste0(
        "status_ERROR=",
        sum(object_audit$status == "ERROR")
    ),
    paste0(
        "cells_total=",
        format(
            sum(object_audit$cells, na.rm = TRUE),
            scientific = FALSE,
            trim = TRUE
        )
    ),
    paste0(
        "cells_phase4a=",
        format(
            sum(phase4a_audit$cells, na.rm = TRUE),
            scientific = FALSE,
            trim = TRUE
        )
    ),
    paste0(
        "cells_phase4c=",
        format(
            sum(phase4c_audit$cells, na.rm = TRUE),
            scientific = FALSE,
            trim = TRUE
        )
    ),
    paste0(
        "unique_cell_names=",
        length(unique(all_cell_names))
    ),
    paste0(
        "cross_object_cell_overlap_total=",
        sum(
            object_audit$cross_object_cell_overlap,
            na.rm = TRUE
        )
    ),
    paste0(
        "raw_feature_sets=",
        length(unique(
            feature_summary$feature_set_id
        ))
    ),
    paste0(
        "canonical_feature_sets=",
        length(unique(
            feature_summary$canonical_feature_set_id
        ))
    ),
    paste0(
        "features_common_to_all_raw=",
        length(common_features)
    ),
    paste0(
        "features_common_to_all_canonical=",
        length(canonical_common_features)
    ),
    paste0(
        "canonical_collision_objects=",
        sum(
            object_audit$canonical_duplicate_features >
                0L,
            na.rm = TRUE
        )
    ),
    paste0(
        "duplicate_feature_objects=",
        sum(
            object_audit$duplicate_features > 0L,
            na.rm = TRUE
        )
    ),
    paste0(
        "cross_group_pair_count=",
        nrow(cross_group_pairs)
    ),
    paste0(
        "cross_group_raw_intersection_min=",
        min(
            cross_group_pairs$raw_intersection,
            na.rm = TRUE
        )
    ),
    paste0(
        "cross_group_raw_intersection_max=",
        max(
            cross_group_pairs$raw_intersection,
            na.rm = TRUE
        )
    ),
    paste0(
        "cross_group_canonical_intersection_min=",
        min(
            cross_group_pairs$canonical_intersection,
            na.rm = TRUE
        )
    ),
    paste0(
        "cross_group_canonical_intersection_max=",
        max(
            cross_group_pairs$canonical_intersection,
            na.rm = TRUE
        )
    ),
    paste0(
        "object_audit_file=",
        object_audit_file
    ),
    paste0(
        "dataset_summary_file=",
        dataset_summary_file
    ),
    paste0(
        "feature_pair_file=",
        feature_pair_file
    ),
    paste0(
        "metadata_inventory_file=",
        metadata_inventory_file
    )
)

writeLines(
    summary_lines,
    con = summary_file
)

session_connection <- file(
    session_info_file,
    open = "wt"
)

sink(session_connection)
print(sessionInfo())
sink()
close(session_connection)

cat(
    paste(summary_lines, collapse = "\n"),
    "\n",
    sep = ""
)
