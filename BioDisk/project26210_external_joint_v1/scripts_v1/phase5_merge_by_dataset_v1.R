options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 240
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

root <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"

manifest_file <- file.path(
  root,
  "phase4i_finalize_residual_treatment_v1",
  "tables",
  "object_clinical_manifest_phase4i_final.tsv"
)

out_dir <- file.path(
  root,
  "phase5_merge_by_dataset_v1"
)

object_dir <- file.path(out_dir, "objects_raw")
table_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

write_tsv <- function(x, path) {
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\r\n\t]+", " ", x)
  trimws(x)
}

empty_to_na <- function(x) {
  x <- clean_text(x)
  x[x == ""] <- NA_character_
  x
}

normalize_key <- function(x) {
  x <- iconv(
    as.character(x),
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  x[is.na(x)] <- ""
  tolower(gsub("[^a-z0-9]+", "", x))
}

safe_id <- function(x) {
  x <- iconv(
    clean_text(x),
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  x[is.na(x)] <- ""
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  ifelse(x == "", "unknown", x)
}

snake_name <- function(x) {
  x <- iconv(
    as.character(x),
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  empty_index <- x == ""

  if (any(empty_index)) {
    x[empty_index] <- paste0(
      "field_",
      which(empty_index)
    )
  }

  make.unique(x, sep = "__dup")
}

find_column <- function(column_names, candidates) {
  observed_key <- normalize_key(column_names)

  for (candidate in candidates) {
    hit <- which(
      observed_key == normalize_key(candidate)
    )

    if (length(hit)) {
      return(column_names[hit[1]])
    }
  }

  ""
}

standardize_source_column <- function(x) {
  if (
    inherits(x, "POSIXt") ||
    inherits(x, "Date") ||
    is.factor(x)
  ) {
    return(as.character(x))
  }

  if (is.list(x)) {
    return(
      vapply(
        x,
        function(value) {
          if (length(value) == 0L) {
            return(NA_character_)
          }

          paste(
            as.character(value),
            collapse = ";"
          )
        },
        character(1)
      )
    )
  }

  if (
    is.matrix(x) ||
    is.data.frame(x)
  ) {
    return(
      apply(
        as.matrix(x),
        1,
        function(value) {
          paste(value, collapse = ";")
        }
      )
    )
  }

  x
}

get_manifest_scalar <- function(
  manifest_row,
  candidates
) {
  for (candidate in candidates) {
    column_name <- find_column(
      names(manifest_row),
      candidate
    )

    if (column_name == "") {
      next
    }

    value <- clean_text(
      manifest_row[[column_name]][1]
    )

    if (value != "") {
      return(
        list(
          value = value,
          source = paste0(
            "manifest:",
            column_name
          )
        )
      )
    }
  }

  NULL
}

resolve_canonical_field <- function(
  manifest_row,
  source_meta,
  candidates,
  fallback_value = NULL,
  fallback_source = NULL
) {
  manifest_result <- get_manifest_scalar(
    manifest_row,
    candidates
  )

  if (!is.null(manifest_result)) {
    return(
      list(
        values = rep(
          manifest_result$value,
          nrow(source_meta)
        ),
        source = manifest_result$source
      )
    )
  }

  source_column <- find_column(
    names(source_meta),
    candidates
  )

  if (source_column != "") {
    values <- empty_to_na(
      source_meta[[source_column]]
    )

    if (any(!is.na(values))) {
      return(
        list(
          values = values,
          source = paste0(
            "object_metadata:",
            source_column
          )
        )
      )
    }
  }

  if (!is.null(fallback_value)) {
    return(
      list(
        values = rep(
          as.character(fallback_value),
          nrow(source_meta)
        ),
        source = fallback_source
      )
    )
  }

  list(
    values = rep(
      NA_character_,
      nrow(source_meta)
    ),
    source = "unavailable"
  )
}

get_raw_counts <- function(obj) {
  if (!inherits(obj, "Seurat")) {
    stop(
      "Input object is not a Seurat object; class=",
      paste(class(obj), collapse = ";")
    )
  }

  assay_names <- Assays(obj)

  assay_priority <- unique(
    c(
      "RNA",
      DefaultAssay(obj),
      assay_names
    )
  )

  assay_priority <- assay_priority[
    assay_priority %in% assay_names
  ]

  for (assay_name in assay_priority) {
    layer_names <- tryCatch(
      Layers(obj[[assay_name]]),
      error = function(e) character()
    )

    if ("counts" %in% layer_names) {
      counts <- tryCatch(
        LayerData(
          object = obj,
          assay = assay_name,
          layer = "counts"
        ),
        error = function(e) NULL
      )

      if (
        !is.null(counts) &&
        nrow(counts) > 0L &&
        ncol(counts) > 0L
      ) {
        return(
          list(
            counts = counts,
            assay = assay_name,
            source = "layer:counts"
          )
        )
      }
    }

    counts <- tryCatch(
      suppressWarnings(
        GetAssayData(
          object = obj,
          assay = assay_name,
          slot = "counts"
        )
      ),
      error = function(e) NULL
    )

    if (
      !is.null(counts) &&
      nrow(counts) > 0L &&
      ncol(counts) > 0L
    ) {
      return(
        list(
          counts = counts,
          assay = assay_name,
          source = "slot:counts"
        )
      )
    }
  }

  stop(
    "No non-empty raw counts found in assays: ",
    paste(assay_names, collapse = ";")
  )
}

audit_count_values <- function(counts) {
  values <- if (
    inherits(counts, "sparseMatrix")
  ) {
    counts@x
  } else {
    as.vector(counts)
  }

  values <- values[values != 0]

  if (!length(values)) {
    return(
      list(
        negative = 0L,
        non_finite = 0L,
        non_integer = 0L
      )
    )
  }

  list(
    negative = sum(values < 0, na.rm = TRUE),
    non_finite = sum(!is.finite(values)),
    non_integer = sum(
      abs(values - round(values)) > 1e-8,
      na.rm = TRUE
    )
  )
}

bind_rows_safe <- function(x) {
  x <- x[
    vapply(x, nrow, integer(1)) > 0L
  ]

  if (!length(x)) {
    return(data.frame())
  }

  do.call(rbind, x)
}

write_incremental_tables <- function(
  merge_rows,
  crosswalk_rows,
  feature_rows,
  summary_rows,
  schema_rows
) {
  write_tsv(
    bind_rows_safe(merge_rows),
    file.path(
      table_dir,
      "dataset_merge_manifest.tsv"
    )
  )

  write_tsv(
    bind_rows_safe(crosswalk_rows),
    file.path(
      table_dir,
      "metadata_field_crosswalk.tsv"
    )
  )

  write_tsv(
    bind_rows_safe(feature_rows),
    file.path(
      table_dir,
      "feature_set_audit.tsv"
    )
  )

  write_tsv(
    bind_rows_safe(summary_rows),
    file.path(
      table_dir,
      "cell_count_reconciliation.tsv"
    )
  )

  write_tsv(
    bind_rows_safe(schema_rows),
    file.path(
      table_dir,
      "dataset_metadata_schema.tsv"
    )
  )
}

if (!file.exists(manifest_file)) {
  stop(
    "Manifest not found: ",
    manifest_file
  )
}

manifest <- read_tsv(manifest_file)

required_manifest_fields <- c(
  "dataset",
  "object_name",
  "source_file",
  "cells",
  "rds_read_status"
)

missing_manifest_fields <- setdiff(
  required_manifest_fields,
  names(manifest)
)

if (length(missing_manifest_fields)) {
  stop(
    "Manifest missing required fields: ",
    paste(
      missing_manifest_fields,
      collapse = ", "
    )
  )
}

for (
  field_name in c(
    "dataset",
    "object_name",
    "source_file",
    "rds_read_status"
  )
) {
  manifest[[field_name]] <- clean_text(
    manifest[[field_name]]
  )
}

manifest$source_file <- normalizePath(
  manifest$source_file,
  winslash = "/",
  mustWork = FALSE
)

if (any(manifest$dataset == "")) {
  stop("Manifest contains empty dataset")
}

if (any(manifest$object_name == "")) {
  stop("Manifest contains empty object_name")
}

if (any(manifest$source_file == "")) {
  stop("Manifest contains empty source_file")
}

if (any(!file.exists(manifest$source_file))) {
  missing_paths <- manifest$source_file[
    !file.exists(manifest$source_file)
  ]

  writeLines(
    missing_paths,
    file.path(
      qc_dir,
      "missing_source_files.txt"
    )
  )

  stop(
    length(missing_paths),
    " source files are missing"
  )
}

object_key <- paste(
  manifest$dataset,
  manifest$object_name,
  sep = "\r"
)

if (anyDuplicated(object_key)) {
  duplicate_rows <- manifest[
    duplicated(object_key) |
    duplicated(object_key, fromLast = TRUE),
    ,
    drop = FALSE
  ]

  write_tsv(
    duplicate_rows,
    file.path(
      qc_dir,
      "duplicate_manifest_objects.tsv"
    )
  )

  stop(
    "Duplicate dataset/object_name keys found"
  )
}

existing_objects <- list.files(
  object_dir,
  pattern = "_merged_raw\\.rds$",
  full.names = TRUE
)

if (length(existing_objects)) {
  stop(
    "Output RDS already exists. No file was overwritten: ",
    paste(
      basename(existing_objects),
      collapse = ";"
    )
  )
}

canonical_specs <- list(
  sample_id = c(
    "sample_id",
    "sample",
    "sample_name",
    "sampleid",
    "geo_accession",
    "gsm"
  ),
  patient_id = c(
    "patient_id",
    "patient",
    "patientid",
    "donor_id",
    "donor",
    "subject_id",
    "subject",
    "case_id",
    "individual_id"
  ),
  tissue_status = c(
    "tissue_status",
    "tissue_type",
    "tissue",
    "sample_type",
    "tumor_normal",
    "tumour_normal",
    "site_status",
    "region",
    "location"
  ),
  disease_status = c(
    "disease_status",
    "disease",
    "diagnosis",
    "condition",
    "clinical_group"
  ),
  sex = c(
    "sex",
    "gender"
  ),
  age = c(
    "age",
    "age_years",
    "age_at_diagnosis"
  ),
  treatment = c(
    "treatment",
    "therapy",
    "treatment_status"
  ),
  timepoint = c(
    "timepoint",
    "time_point",
    "visit",
    "collection_time"
  ),
  batch_id = c(
    "batch_id",
    "batch",
    "library_id",
    "library",
    "sequencing_batch"
  )
)

datasets <- sort(unique(manifest$dataset))

merge_rows <- list()
crosswalk_rows <- list()
feature_rows <- list()
summary_rows <- list()
schema_rows <- list()

cat("===== PHASE 5 DATASET MERGE =====\n")
cat("started_at=", format(Sys.time()), "\n", sep = "")
cat("Seurat=", as.character(packageVersion("Seurat")), "\n", sep = "")
cat(
  "SeuratObject=",
  as.character(packageVersion("SeuratObject")),
  "\n",
  sep = ""
)
cat("manifest_objects=", nrow(manifest), "\n", sep = "")
cat("datasets=", length(datasets), "\n", sep = "")
cat(
  "dataset_names=",
  paste(datasets, collapse = ";"),
  "\n",
  sep = ""
)

for (
  dataset_index in seq_along(datasets)
) {
  dataset_name <- datasets[dataset_index]

  dataset_manifest <- manifest[
    manifest$dataset == dataset_name,
    ,
    drop = FALSE
  ]

  dataset_manifest <- dataset_manifest[
    order(dataset_manifest$object_name),
    ,
    drop = FALSE
  ]

  cat(
    "\n===== DATASET ",
    dataset_index,
    "/",
    length(datasets),
    ": ",
    dataset_name,
    " =====\n",
    sep = ""
  )

  sample_objects <- vector(
    "list",
    nrow(dataset_manifest)
  )

  feature_sets <- vector(
    "list",
    nrow(dataset_manifest)
  )

  dataset_loaded_cells <- 0

  for (
    object_index in seq_len(
      nrow(dataset_manifest)
    )
  ) {
    manifest_row <- dataset_manifest[
      object_index,
      ,
      drop = FALSE
    ]

    object_name <- manifest_row$object_name[1]
    source_file <- manifest_row$source_file[1]

    cat(
      "[",
      object_index,
      "/",
      nrow(dataset_manifest),
      "] ",
      object_name,
      "\n",
      sep = ""
    )

    obj <- readRDS(source_file)

    if (!inherits(obj, "Seurat")) {
      stop(
        "Unsupported object class for ",
        object_name,
        ": ",
        paste(class(obj), collapse = ";")
      )
    }

    source_meta <- as.data.frame(
      obj[[]],
      check.names = FALSE
    )

    raw_result <- get_raw_counts(obj)
    counts <- raw_result$counts

    if (is.null(rownames(counts))) {
      stop(
        "Counts have no feature names: ",
        object_name
      )
    }

    if (is.null(colnames(counts))) {
      stop(
        "Counts have no cell names: ",
        object_name
      )
    }

    if (anyDuplicated(rownames(counts))) {
      stop(
        "Duplicated feature names: ",
        object_name
      )
    }

    if (anyDuplicated(colnames(counts))) {
      stop(
        "Duplicated cell barcodes: ",
        object_name
      )
    }

    if (
      is.null(rownames(source_meta)) ||
      !setequal(
        rownames(source_meta),
        colnames(counts)
      )
    ) {
      stop(
        "Metadata/count cell names do not match: ",
        object_name
      )
    }

    source_meta <- source_meta[
      colnames(counts),
      ,
      drop = FALSE
    ]

    count_audit <- audit_count_values(counts)

    if (
      count_audit$negative > 0L ||
      count_audit$non_finite > 0L ||
      count_audit$non_integer > 0L
    ) {
      stop(
        "Counts are not valid non-negative integers for ",
        object_name,
        ": negative=",
        count_audit$negative,
        ", non_finite=",
        count_audit$non_finite,
        ", non_integer=",
        count_audit$non_integer
      )
    }

    original_barcode <- colnames(counts)

    dataset_safe <- safe_id(dataset_name)
    object_safe <- safe_id(object_name)

    cell_id <- paste(
      dataset_safe,
      object_safe,
      original_barcode,
      sep = "__"
    )

    if (anyDuplicated(cell_id)) {
      stop(
        "Cell ID collision inside object: ",
        object_name
      )
    }

    colnames(counts) <- cell_id

    canonical_meta <- data.frame(
      row.names = cell_id,
      cell_id = cell_id,
      barcode_original = original_barcode,
      dataset = rep(
        dataset_name,
        length(cell_id)
      ),
      object_name = rep(
        object_name,
        length(cell_id)
      ),
      source_file = rep(
        source_file,
        length(cell_id)
      ),
      source_assay = rep(
        raw_result$assay,
        length(cell_id)
      ),
      source_count_location = rep(
        raw_result$source,
        length(cell_id)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    resolved_fields <- list()

    for (
      canonical_field in names(
        canonical_specs
      )
    ) {
      if (canonical_field == "sample_id") {
        fallback_value <- object_name
        fallback_source <- "fallback:object_name"
      } else if (canonical_field == "batch_id") {
        sample_result <- resolved_fields$sample_id

        if (is.null(sample_result)) {
          fallback_value <- object_name
          fallback_source <- "fallback:object_name"
        } else {
          fallback_value <- sample_result$values
          fallback_source <- "fallback:sample_id"
        }
      } else {
        fallback_value <- NULL
        fallback_source <- NULL
      }

      result <- resolve_canonical_field(
        manifest_row = manifest_row,
        source_meta = source_meta,
        candidates = canonical_specs[[canonical_field]],
        fallback_value = fallback_value,
        fallback_source = fallback_source
      )

      if (
        canonical_field == "batch_id" &&
        length(result$values) == 1L
      ) {
        result$values <- rep(
          result$values,
          nrow(source_meta)
        )
      }

      resolved_fields[[canonical_field]] <- result

      canonical_meta[
        [canonical_field]
      ] <- result$values

      non_missing <- sum(
        !is.na(result$values) &
        clean_text(result$values) != ""
      )

      unique_values <- unique(
        clean_text(result$values)
      )

      unique_values <- unique_values[
        unique_values != ""
      ]

      crosswalk_rows[
        [length(crosswalk_rows) + 1L]
      ] <- data.frame(
        dataset = dataset_name,
        object_name = object_name,
        canonical_field = canonical_field,
        selected_source = result$source,
        non_missing_cells = non_missing,
        unique_values = length(
          unique_values
        ),
        stringsAsFactors = FALSE
      )
    }

    canonical_meta$orig.ident <- canonical_meta$sample_id

    source_meta_standardized <- source_meta

    if (ncol(source_meta_standardized)) {
      source_meta_standardized[] <- lapply(
        source_meta_standardized,
        standardize_source_column
      )

      names(source_meta_standardized) <- paste0(
        "src__",
        snake_name(
          names(source_meta_standardized)
        )
      )

      rownames(source_meta_standardized) <- cell_id
    }

    standardized_meta <- cbind(
      canonical_meta,
      source_meta_standardized
    )

    rownames(standardized_meta) <- cell_id

    sample_obj <- CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = dataset_name,
      meta.data = standardized_meta,
      min.cells = 0,
      min.features = 0
    )

    sample_obj$orig.ident <- sample_obj$sample_id

    if (ncol(sample_obj) != ncol(counts)) {
      stop(
        "Cell loss while rebuilding object: ",
        object_name
      )
    }

    expected_cells <- suppressWarnings(
      as.numeric(
        gsub(
          ",",
          "",
          manifest_row$cells[1],
          fixed = TRUE
        )
      )
    )

    observed_cells <- ncol(sample_obj)

    expected_match <- if (
      is.na(expected_cells)
    ) {
      NA
    } else {
      expected_cells == observed_cells
    }

    merge_rows[
      [length(merge_rows) + 1L]
    ] <- data.frame(
      dataset = dataset_name,
      object_name = object_name,
      source_file = source_file,
      source_assay = raw_result$assay,
      source_count_location = raw_result$source,
      manifest_cells = expected_cells,
      loaded_cells = observed_cells,
      cell_count_match = expected_match,
      features = nrow(sample_obj),
      original_metadata_fields = ncol(
        source_meta
      ),
      standardized_metadata_fields = ncol(
        standardized_meta
      ),
      negative_count_values = count_audit$negative,
      non_finite_count_values = count_audit$non_finite,
      non_integer_count_values = count_audit$non_integer,
      stringsAsFactors = FALSE
    )

    feature_sets[[object_index]] <- rownames(
      sample_obj
    )

    feature_rows[
      [length(feature_rows) + 1L]
    ] <- data.frame(
      dataset = dataset_name,
      object_name = object_name,
      features = nrow(sample_obj),
      same_order_as_dataset_first = if (
        object_index == 1L
      ) {
        TRUE
      } else {
        identical(
          feature_sets[[object_index]],
          feature_sets[[1]]
        )
      },
      same_set_as_dataset_first = if (
        object_index == 1L
      ) {
        TRUE
      } else {
        setequal(
          feature_sets[[object_index]],
          feature_sets[[1]]
        )
      },
      stringsAsFactors = FALSE
    )

    dataset_loaded_cells <- dataset_loaded_cells +
      observed_cells

    sample_objects[[object_index]] <- sample_obj

    rm(
      obj,
      counts,
      source_meta,
      source_meta_standardized,
      standardized_meta,
      canonical_meta,
      sample_obj
    )

    gc(verbose = FALSE)
  }

  all_cell_ids <- unlist(
    lapply(
      sample_objects,
      colnames
    ),
    use.names = FALSE
  )

  if (anyDuplicated(all_cell_ids)) {
    stop(
      "Cross-sample cell ID collision in dataset: ",
      dataset_name
    )
  }

  feature_union <- Reduce(
    union,
    feature_sets
  )

  feature_intersection <- Reduce(
    intersect,
    feature_sets
  )

  if (length(sample_objects) == 1L) {
    merged <- sample_objects[[1]]
  } else {
    merge_method <- getS3method(
      "merge",
      "Seurat"
    )

    merge_arguments <- list(
      x = sample_objects[[1]],
      y = sample_objects[-1],
      add.cell.ids = NULL,
      merge.data = FALSE,
      project = dataset_name
    )

    merge_formals <- names(
      formals(merge_method)
    )

    if ("collapse" %in% merge_formals) {
      merge_arguments$collapse <- TRUE
    }

    if ("merge.dr" %in% merge_formals) {
      merge_arguments$merge.dr <- FALSE
    }

    merged <- do.call(
      merge,
      merge_arguments
    )
  }

  DefaultAssay(merged) <- "RNA"

  count_layers <- tryCatch(
    Layers(
      merged[["RNA"]],
      search = "^counts"
    ),
    error = function(e) character()
  )

  if (length(count_layers) > 1L) {
    if (
      exists(
        "JoinLayers",
        envir = asNamespace("SeuratObject"),
        inherits = FALSE
      )
    ) {
      join_layers_function <- get(
        "JoinLayers",
        envir = asNamespace("SeuratObject")
      )

      merged <- join_layers_function(
        object = merged,
        assay = "RNA",
        layers = count_layers,
        new = "counts"
      )
    } else {
      stop(
        "Multiple count layers remain, but JoinLayers is unavailable: ",
        paste(count_layers, collapse = ";")
      )
    }
  }

  if (ncol(merged) != dataset_loaded_cells) {
    stop(
      "Dataset cell count is not conserved for ",
      dataset_name,
      ": before=",
      dataset_loaded_cells,
      ", after=",
      ncol(merged)
    )
  }

  if (!setequal(rownames(merged), feature_union)) {
    stop(
      "Merged feature set does not equal feature union: ",
      dataset_name
    )
  }

  merged_meta <- merged[[]]

  required_metadata_fields <- c(
    "cell_id",
    "barcode_original",
    "dataset",
    "object_name",
    "sample_id",
    "patient_id",
    "tissue_status",
    "disease_status",
    "sex",
    "age",
    "treatment",
    "timepoint",
    "batch_id",
    "source_file",
    "source_assay",
    "source_count_location",
    "orig.ident"
  )

  missing_metadata_fields <- setdiff(
    required_metadata_fields,
    names(merged_meta)
  )

  if (length(missing_metadata_fields)) {
    stop(
      "Merged object missing standardized fields: ",
      paste(
        missing_metadata_fields,
        collapse = ";"
      )
    )
  }

  if (!identical(
    rownames(merged_meta),
    merged_meta$cell_id
  )) {
    stop(
      "cell_id does not match metadata row names: ",
      dataset_name
    )
  }

  if (
    length(
      unique(
        clean_text(merged_meta$dataset)
      )
    ) != 1L ||
    unique(
      clean_text(merged_meta$dataset)
    ) != dataset_name
  ) {
    stop(
      "Dataset metadata is inconsistent: ",
      dataset_name
    )
  }

  merged@misc$phase5_merge_provenance <- list(
    dataset = dataset_name,
    created_at = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    ),
    manifest_file = manifest_file,
    source_objects = dataset_manifest$object_name,
    source_files = dataset_manifest$source_file,
    merge_unit = "dataset",
    count_assay = "RNA",
    count_policy = "raw_nonnegative_integer_counts_only",
    normalized_data_carried_forward = FALSE,
    reductions_carried_forward = FALSE
  )

  dataset_file_id <- safe_id(
    dataset_name
  )

  output_file <- file.path(
    object_dir,
    paste0(
      dataset_file_id,
      "_merged_raw.rds"
    )
  )

  temporary_file <- paste0(
    output_file,
    ".part.",
    Sys.getpid()
  )

  saveRDS(
    merged,
    temporary_file,
    compress = FALSE
  )

  if (!file.rename(
    temporary_file,
    output_file
  )) {
    stop(
      "Could not finalize output file: ",
      output_file
    )
  }

  output_size <- file.info(
    output_file
  )$size

  expected_dataset_cells <- suppressWarnings(
    as.numeric(
      gsub(
        ",",
        "",
        dataset_manifest$cells,
        fixed = TRUE
      )
    )
  )

  expected_cells_available <- all(
    !is.na(expected_dataset_cells)
  )

  expected_cells_total <- if (
    expected_cells_available
  ) {
    sum(expected_dataset_cells)
  } else {
    NA_real_
  }

  summary_rows[
    [length(summary_rows) + 1L]
  ] <- data.frame(
    dataset = dataset_name,
    source_objects = nrow(
      dataset_manifest
    ),
    manifest_cells = expected_cells_total,
    loaded_cells = dataset_loaded_cells,
    merged_cells = ncol(merged),
    loaded_vs_merged_match = (
      dataset_loaded_cells == ncol(merged)
    ),
    manifest_vs_merged_match = if (
      is.na(expected_cells_total)
    ) {
      NA
    } else {
      expected_cells_total == ncol(merged)
    },
    feature_union = length(
      feature_union
    ),
    feature_intersection = length(
      feature_intersection
    ),
    merged_features = nrow(merged),
    output_file = output_file,
    output_size_bytes = output_size,
    stringsAsFactors = FALSE
  )

  for (
    metadata_field in names(
      merged_meta
    )
  ) {
    field_values <- merged_meta[[metadata_field]]

    non_missing <- if (
      is.character(field_values) ||
      is.factor(field_values)
    ) {
      sum(
        !is.na(field_values) &
        clean_text(field_values) != ""
      )
    } else {
      sum(!is.na(field_values))
    }

    unique_values <- length(
      unique(
        field_values[
          !is.na(field_values)
        ]
      )
    )

    schema_rows[
      [length(schema_rows) + 1L]
    ] <- data.frame(
      dataset = dataset_name,
      metadata_field = metadata_field,
      field_class = paste(
        class(field_values),
        collapse = ";"
      ),
      cells = length(field_values),
      non_missing_cells = non_missing,
      missing_cells = (
        length(field_values) -
        non_missing
      ),
      unique_values = unique_values,
      standardized_core_field = (
        metadata_field %in%
          required_metadata_fields
      ),
      source_metadata_field = grepl(
        "^src__",
        metadata_field
      ),
      stringsAsFactors = FALSE
    )
  }

  write_incremental_tables(
    merge_rows = merge_rows,
    crosswalk_rows = crosswalk_rows,
    feature_rows = feature_rows,
    summary_rows = summary_rows,
    schema_rows = schema_rows
  )

  cat(
    "DATASET_COMPLETE=",
    dataset_name,
    "\n",
    sep = ""
  )
  cat(
    "SOURCE_OBJECTS=",
    nrow(dataset_manifest),
    "\n",
    sep = ""
  )
  cat(
    "CELLS=",
    ncol(merged),
    "\n",
    sep = ""
  )
  cat(
    "FEATURES=",
    nrow(merged),
    "\n",
    sep = ""
  )
  cat(
    "OUTPUT=",
    output_file,
    "\n",
    sep = ""
  )
  cat(
    "OUTPUT_SIZE_GB=",
    round(
      output_size / 1024^3,
      3
    ),
    "\n",
    sep = ""
  )

  rm(
    merged,
    merged_meta,
    sample_objects,
    feature_sets,
    feature_union,
    feature_intersection
  )

  gc(verbose = FALSE)
}

final_summary <- bind_rows_safe(
  summary_rows
)

all_cells_conserved <- all(
  final_summary$loaded_vs_merged_match
)

manifest_cells_conserved <- all(
  is.na(
    final_summary$manifest_vs_merged_match
  ) |
  final_summary$manifest_vs_merged_match
)

all_outputs_exist <- all(
  file.exists(
    final_summary$output_file
  )
)

status <- if (
  nrow(final_summary) == length(datasets) &&
  all_cells_conserved &&
  manifest_cells_conserved &&
  all_outputs_exist
) {
  "PASS"
} else {
  "FAIL"
}

status_lines <- c(
  paste0("status=", status),
  paste0(
    "datasets_expected=",
    length(datasets)
  ),
  paste0(
    "datasets_completed=",
    nrow(final_summary)
  ),
  paste0(
    "source_objects=",
    nrow(manifest)
  ),
  paste0(
    "merged_cells=",
    sum(final_summary$merged_cells)
  ),
  paste0(
    "all_loaded_cells_conserved=",
    all_cells_conserved
  ),
  paste0(
    "all_manifest_cells_conserved=",
    manifest_cells_conserved
  ),
  paste0(
    "all_outputs_exist=",
    all_outputs_exist
  ),
  paste0(
    "completed_at=",
    format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    )
  ),
  "source_rds_modified=0",
  "official_project_modified=0",
  "normalization_performed=0",
  "clustering_performed=0"
)

writeLines(
  status_lines,
  file.path(
    out_dir,
    "PHASE5_MERGE_STATUS"
  )
)

cat("\n===== FINAL DATASET SUMMARY =====\n")

print(
  final_summary[
    ,
    c(
      "dataset",
      "source_objects",
      "manifest_cells",
      "loaded_cells",
      "merged_cells",
      "feature_union",
      "feature_intersection",
      "merged_features",
      "output_size_bytes"
    ),
    drop = FALSE
  ],
  row.names = FALSE,
  right = FALSE
)

cat("\n===== FINAL STATUS =====\n")
cat(status_lines, sep = "\n")
cat("\n")

if (status != "PASS") {
  quit(
    save = "no",
    status = 2L
  )
}
