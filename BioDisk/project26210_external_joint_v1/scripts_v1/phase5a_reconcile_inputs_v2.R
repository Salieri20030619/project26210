options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 240
)

root <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"

manifest_file <- file.path(
  root,
  "phase4i_finalize_residual_treatment_v1",
  "tables",
  "object_clinical_manifest_phase4i_final.tsv"
)

phase5a_v1 <- file.path(
  root,
  "phase5a_audit_annotation_units_v1"
)

out_dir <- file.path(
  root,
  "phase5a_reconcile_inputs_v2"
)

table_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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

normalize_path <- function(x) {
  x <- clean_text(x)

  vapply(
    x,
    function(path) {
      if (path == "") {
        return("")
      }

      normalizePath(
        path,
        winslash = "/",
        mustWork = FALSE
      )
    },
    character(1)
  )
}

collapse_unique <- function(x) {
  x <- sort(unique(clean_text(x)))
  x <- x[x != ""]

  if (!length(x)) {
    return("")
  }

  paste(x, collapse = ";")
}

read_prefix <- function(path, n_bytes = 131072L) {
  connection <- NULL

  result <- tryCatch(
    {
      connection <- if (
        grepl("\\.gz$", path, ignore.case = TRUE)
      ) {
        gzfile(path, open = "rb")
      } else {
        file(path, open = "rb")
      }

      raw_value <- readBin(
        connection,
        what = "raw",
        n = n_bytes
      )

      if (!length(raw_value)) {
        return("")
      }

      raw_value[raw_value == as.raw(0)] <- as.raw(32)

      value <- rawToChar(raw_value)
      value <- enc2utf8(value)

      first_line <- strsplit(
        value,
        "\n",
        fixed = TRUE
      )[[1]][1]

      gsub("\r", "", first_line, fixed = TRUE)
    },
    error = function(e) {
      ""
    },
    finally = {
      if (!is.null(connection)) {
        try(close(connection), silent = TRUE)
      }
    }
  )

  result
}

manifest <- read_tsv(manifest_file)

required_columns <- c(
  "dataset",
  "object_name",
  "source_file",
  "cells",
  "rds_read_status"
)

missing_columns <- setdiff(
  required_columns,
  names(manifest)
)

if (length(missing_columns)) {
  stop(
    "Manifest missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

manifest$dataset <- clean_text(manifest$dataset)
manifest$object_name <- clean_text(manifest$object_name)
manifest$source_file <- clean_text(manifest$source_file)

source_path <- normalize_path(manifest$source_file)
source_exists <- file.exists(source_path)

source_size <- rep(NA_real_, nrow(manifest))

if (any(source_exists)) {
  source_size[source_exists] <- file.info(
    source_path[source_exists]
  )$size
}

expected_basename <- paste0(
  manifest$object_name,
  ".rds"
)

observed_basename <- basename(source_path)

path_validation <- data.frame(
  dataset = manifest$dataset,
  object_name = manifest$object_name,
  source_file = source_path,
  source_exists = source_exists,
  source_size_bytes = source_size,
  expected_basename = expected_basename,
  observed_basename = observed_basename,
  basename_match = (
    expected_basename == observed_basename
  ),
  upstream_rds_read_status = clean_text(
    manifest$rds_read_status
  ),
  stringsAsFactors = FALSE
)

path_validation$corrected_path_status <- ifelse(
  path_validation$source_exists,
  "SOURCE_FILE_CONFIRMED",
  "SOURCE_FILE_MISSING"
)

write_tsv(
  path_validation,
  file.path(
    table_dir,
    "object_source_file_validation.tsv"
  )
)

dataset_names <- sort(unique(manifest$dataset))

path_summary_rows <- lapply(
  dataset_names,
  function(dataset_name) {
    index <- manifest$dataset == dataset_name

    data.frame(
      dataset = dataset_name,
      manifest_objects = sum(index),
      source_files_confirmed = sum(
        path_validation$source_exists[index]
      ),
      source_files_missing = sum(
        !path_validation$source_exists[index]
      ),
      basename_matches = sum(
        path_validation$basename_match[index]
      ),
      total_size_gb = sum(
        path_validation$source_size_bytes[index],
        na.rm = TRUE
      ) / 1024^3,
      stringsAsFactors = FALSE
    )
  }
)

path_summary <- do.call(
  rbind,
  path_summary_rows
)

write_tsv(
  path_summary,
  file.path(
    table_dir,
    "dataset_source_file_summary.tsv"
  )
)

missing_source_files <- path_validation[
  !path_validation$source_exists,
  ,
  drop = FALSE
]

write_tsv(
  missing_source_files,
  file.path(
    qc_dir,
    "missing_source_files.tsv"
  )
)

object_audit_file <- file.path(
  phase5a_v1,
  "tables",
  "object_annotation_audit.tsv"
)

if (file.exists(object_audit_file)) {
  object_audit <- read_tsv(object_audit_file)

  annotation_evidence <- data.frame(
    metric = c(
      "Objects in manifest",
      "Objects read in Phase 5A v1",
      "Objects with selected annotation field",
      "Objects without selected annotation field"
    ),
    observed = c(
      nrow(manifest),
      sum(object_audit$read_status == "READ_OK"),
      sum(clean_text(object_audit$annotation_field) != ""),
      sum(
        object_audit$read_status == "READ_OK" &
        clean_text(object_audit$annotation_field) == ""
      )
    ),
    stringsAsFactors = FALSE
  )
} else {
  annotation_evidence <- data.frame(
    metric = c(
      "Objects in manifest",
      "Objects read in Phase 5A v1",
      "Objects with selected annotation field",
      "Objects without selected annotation field"
    ),
    observed = c(
      nrow(manifest),
      NA,
      NA,
      NA
    ),
    stringsAsFactors = FALSE
  )
}

write_tsv(
  annotation_evidence,
  file.path(
    table_dir,
    "annotation_evidence_summary.tsv"
  )
)

manifest_provenance_columns <- intersect(
  c(
    "dataset",
    "metadata_source",
    "source_reference",
    "evidence_level",
    "notes"
  ),
  names(manifest)
)

manifest_provenance <- do.call(
  rbind,
  lapply(
    dataset_names,
    function(dataset_name) {
      index <- manifest$dataset == dataset_name

      row <- data.frame(
        dataset = dataset_name,
        stringsAsFactors = FALSE
      )

      for (
        column_name in setdiff(
          manifest_provenance_columns,
          "dataset"
        )
      ) {
        row[[column_name]] <- collapse_unique(
          manifest[[column_name]][index]
        )
      }

      row
    }
  )
)

write_tsv(
  manifest_provenance,
  file.path(
    table_dir,
    "dataset_manifest_provenance.tsv"
  )
)

all_files <- list.files(
  root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)

all_files <- normalize_path(all_files)

exclude_pattern <- paste0(
  "^",
  gsub(
    "([][{}()+*^$|\\\\?.])",
    "\\\\\\1",
    normalize_path(out_dir)
  ),
  "/"
)

all_files <- all_files[
  !grepl(exclude_pattern, all_files)
]

source_file_set <- unique(
  source_path[source_path != ""]
)

is_source_rds <- all_files %in% source_file_set

file_info <- file.info(all_files)

file_size <- as.numeric(file_info$size)
file_name <- basename(all_files)
file_name_lower <- tolower(file_name)
path_lower <- tolower(all_files)

excluded_count_file <- grepl(
  paste(
    c(
      "^matrix\\.mtx",
      "^features\\.tsv",
      "^genes\\.tsv",
      "^barcodes\\.tsv",
      "raw_feature_bc_matrix",
      "filtered_feature_bc_matrix"
    ),
    collapse = "|"
  ),
  file_name_lower
)

supported_candidate <- grepl(
  paste0(
    "\\.(csv|tsv|txt|xlsx|xls|h5ad|h5mu|loom|",
    "rds|rda|rdata)(\\.gz)?$"
  ),
  path_lower
)

filename_keyword <- grepl(
  paste(
    c(
      "annot",
      "celltype",
      "cell_type",
      "cell-type",
      "celllabel",
      "cell_label",
      "cluster",
      "metadata",
      "meta_data",
      "meta-data",
      "lineage",
      "identity",
      "integrated",
      "merged"
    ),
    collapse = "|"
  ),
  path_lower
)

table_like <- grepl(
  "\\.(csv|tsv|txt)(\\.gz)?$",
  path_lower
)

header_keyword_pattern <- paste(
  c(
    "cell.?type",
    "cell.?annotation",
    "annotation",
    "cell.?label",
    "major.?type",
    "sub.?type",
    "lineage",
    "cluster",
    "seurat.?cluster",
    "identity",
    "cell.?identity"
  ),
  collapse = "|"
)

header_text <- rep("", length(all_files))
header_keyword <- rep(FALSE, length(all_files))

header_index <- which(
  supported_candidate &
  table_like &
  !excluded_count_file &
  file_size <= 1024^3
)

if (length(header_index)) {
  for (i in header_index) {
    header_text[i] <- read_prefix(all_files[i])

    header_keyword[i] <- grepl(
      header_keyword_pattern,
      header_text[i],
      ignore.case = TRUE,
      perl = TRUE
    )
  }
}

candidate_index <- which(
  supported_candidate &
  !excluded_count_file &
  !is_source_rds &
  (
    filename_keyword |
    header_keyword
  )
)

if (length(candidate_index)) {
  annotation_candidates <- data.frame(
    path = all_files[candidate_index],
    size_bytes = file_size[candidate_index],
    extension = sub(
      "^.*?(\\.[^.]+(?:\\.gz)?)$",
      "\\1",
      file_name_lower[candidate_index]
    ),
    filename_keyword = filename_keyword[candidate_index],
    header_keyword = header_keyword[candidate_index],
    header_excerpt = substr(
      clean_text(header_text[candidate_index]),
      1,
      500
    ),
    stringsAsFactors = FALSE
  )

  annotation_candidates$candidate_priority <- ifelse(
    annotation_candidates$header_keyword,
    1L,
    ifelse(
      annotation_candidates$filename_keyword,
      2L,
      3L
    )
  )

  annotation_candidates <- annotation_candidates[
    order(
      annotation_candidates$candidate_priority,
      annotation_candidates$path
    ),
    ,
    drop = FALSE
  ]
} else {
  annotation_candidates <- data.frame(
    path = character(),
    size_bytes = numeric(),
    extension = character(),
    filename_keyword = logical(),
    header_keyword = logical(),
    header_excerpt = character(),
    candidate_priority = integer(),
    stringsAsFactors = FALSE
  )
}

write_tsv(
  annotation_candidates,
  file.path(
    table_dir,
    "annotation_source_candidates.tsv"
  )
)

rds_index <- which(
  grepl("\\.rds$", path_lower) &
  !is_source_rds &
  !grepl(
    "/phase5a_audit_annotation_units_v1/",
    path_lower,
    fixed = TRUE
  )
)

if (length(rds_index)) {
  non_source_rds <- data.frame(
    path = all_files[rds_index],
    size_bytes = file_size[rds_index],
    size_gb = file_size[rds_index] / 1024^3,
    name_has_annotation_keyword = grepl(
      paste(
        c(
          "annot",
          "celltype",
          "cluster",
          "integrated",
          "merged",
          "final"
        ),
        collapse = "|"
      ),
      path_lower[rds_index]
    ),
    stringsAsFactors = FALSE
  )

  non_source_rds <- non_source_rds[
    order(
      -as.integer(
        non_source_rds$name_has_annotation_keyword
      ),
      -non_source_rds$size_bytes,
      non_source_rds$path
    ),
    ,
    drop = FALSE
  ]
} else {
  non_source_rds <- data.frame(
    path = character(),
    size_bytes = numeric(),
    size_gb = numeric(),
    name_has_annotation_keyword = logical(),
    stringsAsFactors = FALSE
  )
}

write_tsv(
  non_source_rds,
  file.path(
    table_dir,
    "non_source_rds_candidates.tsv"
  )
)

all_paths_confirmed <- all(
  path_validation$source_exists
)

path_status <- if (all_paths_confirmed) {
  "ALL_95_SOURCE_FILES_CONFIRMED"
} else {
  "SOURCE_FILES_MISSING"
}

annotation_status <- if (
  nrow(annotation_candidates) > 0L ||
  nrow(non_source_rds) > 0L
) {
  "CANDIDATE_SOURCES_REQUIRE_REVIEW"
} else {
  "NO_EXISTING_ANNOTATION_SOURCE_DETECTED"
}

status_lines <- c(
  paste0("path_status=", path_status),
  paste0(
    "manifest_objects=",
    nrow(manifest)
  ),
  paste0(
    "confirmed_source_files=",
    sum(path_validation$source_exists)
  ),
  paste0(
    "missing_source_files=",
    sum(!path_validation$source_exists)
  ),
  paste0(
    "basename_matches=",
    sum(path_validation$basename_match)
  ),
  paste0(
    "annotation_status=",
    annotation_status
  ),
  paste0(
    "annotation_file_candidates=",
    nrow(annotation_candidates)
  ),
  paste0(
    "non_source_rds_candidates=",
    nrow(non_source_rds)
  ),
  "rds_objects_read=0",
  "source_rds_modified=0",
  "phase4i_modified=0",
  "customer_report_modified=0"
)

writeLines(
  status_lines,
  file.path(
    out_dir,
    "PHASE5A_INPUT_RECONCILIATION_STATUS"
  )
)

report <- c(
  "# Phase 5A input reconciliation v2",
  "",
  paste0("- Path status: ", path_status),
  paste0(
    "- Confirmed source files: ",
    sum(path_validation$source_exists),
    "/",
    nrow(manifest)
  ),
  paste0(
    "- Existing annotation candidates: ",
    nrow(annotation_candidates)
  ),
  paste0(
    "- Non-source RDS candidates: ",
    nrow(non_source_rds)
  ),
  "",
  "## Interpretation",
  "",
  paste0(
    "Phase 5A v1 found no usable annotation field in ",
    "the successfully read source objects. The source_file ",
    "column is therefore treated as the authoritative object ",
    "path, but no cell-type label is inferred from file names."
  ),
  "",
  "## Safety",
  "",
  paste0(
    "No source RDS, Phase 4I output or customer report was ",
    "modified. No RDS object was loaded by this reconciliation."
  )
)

writeLines(
  report,
  file.path(
    out_dir,
    "PHASE5A_INPUT_RECONCILIATION_REPORT.md"
  ),
  useBytes = TRUE
)

cat("===== PATH VALIDATION SUMMARY =====\n")

print(
  path_summary,
  row.names = FALSE,
  right = FALSE
)

cat("\n===== MISSING SOURCE FILES =====\n")
cat(
  "objects=",
  nrow(missing_source_files),
  "\n",
  sep = ""
)

if (nrow(missing_source_files)) {
  print(
    missing_source_files[
      ,
      c(
        "dataset",
        "object_name",
        "source_file",
        "upstream_rds_read_status"
      ),
      drop = FALSE
    ],
    row.names = FALSE,
    right = FALSE
  )
} else {
  cat("NONE\n")
}

cat("\n===== ANNOTATION EVIDENCE =====\n")

print(
  annotation_evidence,
  row.names = FALSE,
  right = FALSE
)

cat("\n===== DATASET PROVENANCE =====\n")

print(
  manifest_provenance,
  row.names = FALSE,
  right = FALSE
)

cat("\n===== ANNOTATION SOURCE CANDIDATES =====\n")
cat(
  "candidates=",
  nrow(annotation_candidates),
  "\n",
  sep = ""
)

if (nrow(annotation_candidates)) {
  print(
    head(annotation_candidates, 80L),
    row.names = FALSE,
    right = FALSE
  )
} else {
  cat("NONE\n")
}

cat("\n===== NON-SOURCE RDS CANDIDATES =====\n")
cat(
  "candidates=",
  nrow(non_source_rds),
  "\n",
  sep = ""
)

if (nrow(non_source_rds)) {
  print(
    head(non_source_rds, 80L),
    row.names = FALSE,
    right = FALSE
  )
} else {
  cat("NONE\n")
}

cat("\n===== FINAL STATUS =====\n")
cat(status_lines, sep = "\n")
cat("\n")
