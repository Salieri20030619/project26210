#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    warn = 1
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
    stop(
        "Usage: phase4b_qc_failure_audit.R ",
        "<phase4a_root> <output_root>"
    )
}

phase4a_root <- normalizePath(
    args[[1]],
    mustWork = TRUE
)

output_root <- args[[2]]

object_dir <- file.path(
    phase4a_root,
    "objects"
)

qc_summary_path <- file.path(
    phase4a_root,
    "sample_qc_summary.tsv"
)

if (!dir.exists(object_dir)) {
    stop("Object directory not found: ", object_dir)
}

if (!file.exists(qc_summary_path)) {
    stop("QC summary not found: ", qc_summary_path)
}

dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
)

object_files <- list.files(
    object_dir,
    pattern = "\\.rds$",
    full.names = TRUE
)

object_files <- sort(object_files)

if (length(object_files) != 69L) {
    stop(
        "Expected 69 RDS objects, found ",
        length(object_files)
    )
}

write_tsv <- function(data, path) {
    temporary_path <- paste0(
        path,
        ".tmp_",
        Sys.getpid()
    )

    utils::write.table(
        data,
        file = temporary_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "NA"
    )

    if (!file.rename(temporary_path, path)) {
        stop("Failed to create: ", path)
    }
}

get_provenance_value <- function(provenance, name) {
    value <- provenance[[name]]

    if (is.null(value) || length(value) == 0L) {
        return(NA_character_)
    }

    as.character(value[[1]])
}

metric_quantiles <- function(values) {
    probabilities <- c(
        0.01,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99
    )

    result <- stats::quantile(
        values,
        probs = probabilities,
        na.rm = TRUE,
        names = FALSE,
        type = 7
    )

    names(result) <- c(
        "p01",
        "p05",
        "p10",
        "p25",
        "p50",
        "p75",
        "p90",
        "p95",
        "p99"
    )

    result
}

sample_records <- vector(
    mode = "list",
    length = length(object_files)
)

quantile_records <- vector(
    mode = "list",
    length = length(object_files)
)

reason_records <- vector(
    mode = "list",
    length = length(object_files)
)

for (index in seq_along(object_files)) {
    object_path <- object_files[[index]]

    message(
        sprintf(
            "[%02d/%02d] %s",
            index,
            length(object_files),
            basename(object_path)
        )
    )

    object <- readRDS(object_path)

    if (!methods::is(object, "Seurat")) {
        stop("Not a Seurat object: ", object_path)
    }

    metadata <- object@meta.data

    required_columns <- c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "percent.ribo",
        "qc_pass_basic"
    )

    missing_columns <- setdiff(
        required_columns,
        colnames(metadata)
    )

    if (length(missing_columns) > 0L) {
        stop(
            "Missing metadata columns in ",
            object_path,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }

    provenance <- object@misc$phase4a_provenance

    if (is.null(provenance)) {
        stop(
            "Missing phase4a provenance: ",
            object_path
        )
    }

    dataset <- get_provenance_value(
        provenance,
        "dataset"
    )

    sample_key_raw <- get_provenance_value(
        provenance,
        "sample_key"
    )

    gsm <- get_provenance_value(
        provenance,
        "gsm"
    )

    if (is.na(dataset) || is.na(sample_key_raw)) {
        stop(
            "Incomplete provenance in: ",
            object_path
        )
    }

    missing_metric <- (
        is.na(metadata$nFeature_RNA) |
        is.na(metadata$nCount_RNA) |
        is.na(metadata$percent.mt)
    )

    low_feature <- (
        !is.na(metadata$nFeature_RNA) &
        metadata$nFeature_RNA < 200
    )

    low_count <- (
        !is.na(metadata$nCount_RNA) &
        metadata$nCount_RNA < 500
    )

    high_mt <- (
        !is.na(metadata$percent.mt) &
        metadata$percent.mt > 30
    )

    fail_any <- (
        missing_metric |
        low_feature |
        low_count |
        high_mt
    )

    pass_calculated <- !fail_any

    pass_recorded <- (
        !is.na(metadata$qc_pass_basic) &
        as.logical(metadata$qc_pass_basic)
    )

    if (!identical(
        as.logical(pass_calculated),
        as.logical(pass_recorded)
    )) {
        stop(
            "Recorded and recalculated QC flags differ: ",
            object_path
        )
    }

    failure_mask <- (
        as.integer(low_feature) +
        2L * as.integer(low_count) +
        4L * as.integer(high_mt)
    )

    reason_labels <- c(
        "Pass",
        "LowFeature",
        "LowCount",
        "LowFeature+LowCount",
        "HighMT",
        "LowFeature+HighMT",
        "LowCount+HighMT",
        "LowFeature+LowCount+HighMT"
    )

    failure_reason <- reason_labels[
        failure_mask + 1L
    ]

    failure_reason[missing_metric] <- "MissingMetric"

    reason_table <- as.data.frame(
        table(failure_reason),
        stringsAsFactors = FALSE
    )

    colnames(reason_table) <- c(
        "failure_reason",
        "cells"
    )

    reason_table$dataset <- dataset
    reason_table$sample_key_raw <- sample_key_raw
    reason_table$gsm <- gsm

    reason_table <- reason_table[
        ,
        c(
            "dataset",
            "sample_key_raw",
            "gsm",
            "failure_reason",
            "cells"
        )
    ]

    reason_records[[index]] <- reason_table

    cell_count <- nrow(metadata)
    pass_count <- sum(pass_calculated)
    fail_count <- sum(fail_any)

    sample_records[[index]] <- data.frame(
        dataset = dataset,
        sample_key_raw = sample_key_raw,
        gsm = gsm,
        object_file = object_path,
        cells = cell_count,
        qc_pass_cells = pass_count,
        qc_fail_cells = fail_count,
        qc_pass_rate = pass_count / cell_count,
        low_feature_cells = sum(low_feature),
        low_feature_rate = sum(low_feature) / cell_count,
        low_count_cells = sum(low_count),
        low_count_rate = sum(low_count) / cell_count,
        high_mt_cells = sum(high_mt),
        high_mt_rate = sum(high_mt) / cell_count,
        missing_metric_cells = sum(missing_metric),
        stringsAsFactors = FALSE
    )

    feature_q <- metric_quantiles(
        metadata$nFeature_RNA
    )

    count_q <- metric_quantiles(
        metadata$nCount_RNA
    )

    mt_q <- metric_quantiles(
        metadata$percent.mt
    )

    ribo_q <- metric_quantiles(
        metadata$percent.ribo
    )

    quantile_records[[index]] <- data.frame(
        dataset = dataset,
        sample_key_raw = sample_key_raw,
        gsm = gsm,
        nFeature_p01 = feature_q[["p01"]],
        nFeature_p05 = feature_q[["p05"]],
        nFeature_p10 = feature_q[["p10"]],
        nFeature_p25 = feature_q[["p25"]],
        nFeature_p50 = feature_q[["p50"]],
        nFeature_p75 = feature_q[["p75"]],
        nFeature_p90 = feature_q[["p90"]],
        nFeature_p95 = feature_q[["p95"]],
        nFeature_p99 = feature_q[["p99"]],
        nCount_p01 = count_q[["p01"]],
        nCount_p05 = count_q[["p05"]],
        nCount_p10 = count_q[["p10"]],
        nCount_p25 = count_q[["p25"]],
        nCount_p50 = count_q[["p50"]],
        nCount_p75 = count_q[["p75"]],
        nCount_p90 = count_q[["p90"]],
        nCount_p95 = count_q[["p95"]],
        nCount_p99 = count_q[["p99"]],
        percent_mt_p01 = mt_q[["p01"]],
        percent_mt_p05 = mt_q[["p05"]],
        percent_mt_p10 = mt_q[["p10"]],
        percent_mt_p25 = mt_q[["p25"]],
        percent_mt_p50 = mt_q[["p50"]],
        percent_mt_p75 = mt_q[["p75"]],
        percent_mt_p90 = mt_q[["p90"]],
        percent_mt_p95 = mt_q[["p95"]],
        percent_mt_p99 = mt_q[["p99"]],
        percent_ribo_p01 = ribo_q[["p01"]],
        percent_ribo_p05 = ribo_q[["p05"]],
        percent_ribo_p10 = ribo_q[["p10"]],
        percent_ribo_p25 = ribo_q[["p25"]],
        percent_ribo_p50 = ribo_q[["p50"]],
        percent_ribo_p75 = ribo_q[["p75"]],
        percent_ribo_p90 = ribo_q[["p90"]],
        percent_ribo_p95 = ribo_q[["p95"]],
        percent_ribo_p99 = ribo_q[["p99"]],
        stringsAsFactors = FALSE
    )

    rm(
        object,
        metadata,
        low_feature,
        low_count,
        high_mt,
        missing_metric,
        fail_any,
        pass_calculated,
        pass_recorded,
        failure_mask,
        failure_reason
    )

    invisible(gc())
}

sample_summary <- do.call(
    rbind,
    sample_records
)

quantile_summary <- do.call(
    rbind,
    quantile_records
)

reason_long <- do.call(
    rbind,
    reason_records
)

sample_summary$sample_key_clean <- gsub(
    "[._-]+$",
    "",
    sample_summary$sample_key_raw
)

sample_summary$canonical_sample_id <- paste(
    sample_summary$dataset,
    sample_summary$sample_key_clean,
    sep = "__"
)

if (anyDuplicated(
    sample_summary$canonical_sample_id
)) {
    duplicated_ids <- unique(
        sample_summary$canonical_sample_id[
            duplicated(
                sample_summary$canonical_sample_id
            ) |
            duplicated(
                sample_summary$canonical_sample_id,
                fromLast = TRUE
            )
        ]
    )

    stop(
        "Duplicate canonical sample IDs: ",
        paste(duplicated_ids, collapse = ", ")
    )
}

canonical_map <- sample_summary[
    ,
    c(
        "dataset",
        "sample_key_raw",
        "sample_key_clean",
        "canonical_sample_id",
        "gsm",
        "object_file"
    )
]

phase4a_summary <- utils::read.delim(
    qc_summary_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
)

phase4a_key <- paste(
    phase4a_summary$dataset,
    phase4a_summary$sample_key,
    sep = "\t"
)

audit_key <- paste(
    sample_summary$dataset,
    sample_summary$sample_key_raw,
    sep = "\t"
)

matched_index <- match(
    audit_key,
    phase4a_key
)

if (anyNA(matched_index)) {
    stop(
        "Audit samples do not match Phase 4A summary"
    )
}

if (any(
    sample_summary$cells !=
    phase4a_summary$cells[matched_index]
)) {
    stop(
        "Cell count mismatch against Phase 4A summary"
    )
}

if (any(
    sample_summary$qc_pass_cells !=
    phase4a_summary$basic_qc_pass_cells[
        matched_index
    ]
)) {
    stop(
        "QC pass count mismatch against Phase 4A summary"
    )
}

dataset_order <- c(
    "GSE207493",
    "GSE294109",
    "GSE264586",
    "GSE202374"
)

dataset_summary_list <- lapply(
    dataset_order,
    function(current_dataset) {
        subset_data <- sample_summary[
            sample_summary$dataset == current_dataset,
            ,
            drop = FALSE
        ]

        data.frame(
            dataset = current_dataset,
            samples = nrow(subset_data),
            cells = sum(subset_data$cells),
            qc_pass_cells = sum(
                subset_data$qc_pass_cells
            ),
            qc_fail_cells = sum(
                subset_data$qc_fail_cells
            ),
            qc_pass_rate = sum(
                subset_data$qc_pass_cells
            ) / sum(subset_data$cells),
            low_feature_cells = sum(
                subset_data$low_feature_cells
            ),
            low_count_cells = sum(
                subset_data$low_count_cells
            ),
            high_mt_cells = sum(
                subset_data$high_mt_cells
            ),
            missing_metric_cells = sum(
                subset_data$missing_metric_cells
            ),
            stringsAsFactors = FALSE
        )
    }
)

dataset_summary <- do.call(
    rbind,
    dataset_summary_list
)

overall_summary <- data.frame(
    dataset = "ALL",
    samples = nrow(sample_summary),
    cells = sum(sample_summary$cells),
    qc_pass_cells = sum(
        sample_summary$qc_pass_cells
    ),
    qc_fail_cells = sum(
        sample_summary$qc_fail_cells
    ),
    qc_pass_rate = sum(
        sample_summary$qc_pass_cells
    ) / sum(sample_summary$cells),
    low_feature_cells = sum(
        sample_summary$low_feature_cells
    ),
    low_count_cells = sum(
        sample_summary$low_count_cells
    ),
    high_mt_cells = sum(
        sample_summary$high_mt_cells
    ),
    missing_metric_cells = sum(
        sample_summary$missing_metric_cells
    ),
    stringsAsFactors = FALSE
)

dataset_summary <- rbind(
    dataset_summary,
    overall_summary
)

reason_by_dataset <- stats::aggregate(
    cells ~ dataset + failure_reason,
    data = reason_long,
    FUN = sum
)

reason_overall <- stats::aggregate(
    cells ~ failure_reason,
    data = reason_long,
    FUN = sum
)

reason_overall$rate <- (
    reason_overall$cells /
    sum(sample_summary$cells)
)

sample_summary <- sample_summary[
    order(
        sample_summary$qc_pass_rate,
        sample_summary$dataset,
        sample_summary$sample_key_raw
    ),
    ,
    drop = FALSE
]

canonical_map <- canonical_map[
    order(
        canonical_map$dataset,
        canonical_map$sample_key_clean
    ),
    ,
    drop = FALSE
]

write_tsv(
    sample_summary,
    file.path(
        output_root,
        "sample_failure_components.tsv"
    )
)

write_tsv(
    quantile_summary,
    file.path(
        output_root,
        "sample_qc_metric_quantiles.tsv"
    )
)

write_tsv(
    reason_long,
    file.path(
        output_root,
        "sample_failure_reason_combinations.tsv"
    )
)

write_tsv(
    dataset_summary,
    file.path(
        output_root,
        "dataset_failure_components.tsv"
    )
)

write_tsv(
    reason_by_dataset,
    file.path(
        output_root,
        "dataset_failure_reason_combinations.tsv"
    )
)

write_tsv(
    reason_overall,
    file.path(
        output_root,
        "overall_failure_reason_combinations.tsv"
    )
)

write_tsv(
    canonical_map,
    file.path(
        output_root,
        "canonical_sample_id_map.tsv"
    )
)

writeLines(
    capture.output(sessionInfo()),
    con = file.path(
        output_root,
        "sessionInfo.txt"
    )
)

completion_lines <- c(
    "status=PASS",
    paste0(
        "objects_audited=",
        nrow(sample_summary)
    ),
    paste0(
        "cells_audited=",
        sum(sample_summary$cells)
    ),
    paste0(
        "qc_pass_cells=",
        sum(sample_summary$qc_pass_cells)
    ),
    paste0(
        "qc_fail_cells=",
        sum(sample_summary$qc_fail_cells)
    ),
    paste0(
        "completed_at=",
        format(
            Sys.time(),
            "%Y-%m-%dT%H:%M:%S%z"
        )
    )
)

writeLines(
    completion_lines,
    con = file.path(
        output_root,
        "PHASE4B_AUDIT_COMPLETE"
    )
)

cat("\n===== DATASET SUMMARY =====\n")
print(
    dataset_summary,
    row.names = FALSE
)

cat("\n===== LOWEST 20 SAMPLE PASS RATES =====\n")
print(
    sample_summary[
        seq_len(min(20L, nrow(sample_summary))),
        c(
            "dataset",
            "sample_key_clean",
            "cells",
            "qc_fail_cells",
            "qc_pass_rate",
            "low_feature_cells",
            "low_count_cells",
            "high_mt_cells",
            "missing_metric_cells"
        )
    ],
    row.names = FALSE
)

cat("\n===== EXACT FAILURE COMBINATIONS =====\n")
print(
    reason_overall[
        order(
            -reason_overall$cells
        ),
        ,
        drop = FALSE
    ],
    row.names = FALSE
)

cat("\nPHASE4B_AUDIT_COMPLETE\n")
