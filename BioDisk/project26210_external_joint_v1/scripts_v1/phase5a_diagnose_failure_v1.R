options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 240
)

root <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
phase4i <- file.path(
  root,
  "phase4i_finalize_residual_treatment_v1"
)
phase5a <- file.path(
  root,
  "phase5a_audit_annotation_units_v1"
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

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\r\n\t]+", " ", x)
  trimws(x)
}

manifest_file <- file.path(
  phase4i,
  "tables",
  "object_clinical_manifest_phase4i_final.tsv"
)

path_file <- file.path(
  phase5a,
  "tables",
  "object_rds_path_audit.tsv"
)

object_file <- file.path(
  phase5a,
  "tables",
  "object_annotation_audit.tsv"
)

structural_file <- file.path(
  phase5a,
  "qc",
  "phase5a_structural_qc.tsv"
)

manifest <- read_tsv(manifest_file)
path_audit <- read_tsv(path_file)
object_audit <- read_tsv(object_file)
structural_qc <- read_tsv(structural_file)

cat("===== PHASE5A STATUS FILE =====\n")

status_candidates <- c(
  file.path(phase5a, "PHASE5A_COMPLETE"),
  file.path(phase5a, "PHASE5A_INCOMPLETE")
)

existing_status <- status_candidates[
  file.exists(status_candidates)
]

if (length(existing_status)) {
  cat(readLines(existing_status[1]), sep = "\n")
  cat("\n")
} else {
  cat("STATUS_FILE_NOT_FOUND\n")
}

cat("\n===== STRUCTURAL QC =====\n")
print(
  structural_qc,
  row.names = FALSE,
  right = FALSE
)

cat("\n===== MANIFEST COLUMNS =====\n")

for (i in seq_along(names(manifest))) {
  cat(
    sprintf(
      "%03d\t%s\n",
      i,
      names(manifest)[i]
    )
  )
}

cat("\n===== PATH STATUS COUNTS =====\n")

path_status_count <- as.data.frame(
  table(path_audit$path_status),
  stringsAsFactors = FALSE
)

names(path_status_count) <- c(
  "path_status",
  "objects"
)

print(
  path_status_count,
  row.names = FALSE,
  right = FALSE
)

failed_path <- path_audit[
  path_audit$rds_path == "" |
  path_audit$path_status %in% c(
    "RDS_NOT_FOUND",
    "AMBIGUOUS_BASENAME"
  ),
  ,
  drop = FALSE
]

cat("\n===== UNRESOLVED PATH OBJECTS =====\n")
cat("objects=", nrow(failed_path), "\n", sep = "")

if (nrow(failed_path)) {
  print(
    failed_path,
    row.names = FALSE,
    right = FALSE
  )
} else {
  cat("NONE\n")
}

path_columns <- names(manifest)[
  grepl(
    "path|file|rds|object",
    names(manifest),
    ignore.case = TRUE
  )
]

path_columns <- unique(
  c(
    intersect(
      c(
        "dataset",
        "object_name"
      ),
      names(manifest)
    ),
    path_columns
  )
)

cat("\n===== MANIFEST PATH-LIKE FIELDS FOR FAILED OBJECTS =====\n")
cat(
  "candidate_columns=",
  paste(path_columns, collapse = ";"),
  "\n",
  sep = ""
)

if (nrow(failed_path)) {
  failed_key <- paste(
    failed_path$dataset,
    failed_path$object_name,
    sep = "\r"
  )

  manifest_key <- paste(
    clean_text(manifest$dataset),
    clean_text(manifest$object_name),
    sep = "\r"
  )

  failed_manifest <- manifest[
    manifest_key %in% failed_key,
    path_columns,
    drop = FALSE
  ]

  print(
    failed_manifest,
    row.names = FALSE,
    right = FALSE
  )
} else {
  cat("NONE\n")
}

cat("\n===== READ STATUS COUNTS =====\n")

read_status_count <- aggregate(
  x = list(
    objects = rep(1L, nrow(object_audit))
  ),
  by = list(
    dataset = object_audit$dataset,
    read_status = object_audit$read_status
  ),
  FUN = sum
)

print(
  read_status_count[
    order(
      read_status_count$dataset,
      read_status_count$read_status
    ),
    ,
    drop = FALSE
  ],
  row.names = FALSE,
  right = FALSE
)

cat("\n===== ANNOTATION FIELD COUNTS =====\n")

annotation_value <- clean_text(
  object_audit$annotation_field
)

annotation_value[
  annotation_value == ""
] <- "<EMPTY>"

annotation_count <- aggregate(
  x = list(
    objects = rep(1L, nrow(object_audit))
  ),
  by = list(
    dataset = object_audit$dataset,
    annotation_field = annotation_value
  ),
  FUN = sum
)

print(
  annotation_count[
    order(
      annotation_count$dataset,
      annotation_count$annotation_field
    ),
    ,
    drop = FALSE
  ],
  row.names = FALSE,
  right = FALSE
)

extract_metadata <- function(obj) {
  if (
    inherits(obj, "Seurat") ||
    inherits(obj, "SeuratObject")
  ) {
    return(obj[[]])
  }

  if (
    is.list(obj) &&
    "meta.data" %in% names(obj) &&
    is.data.frame(obj$meta.data)
  ) {
    return(obj$meta.data)
  }

  if (is.data.frame(obj)) {
    return(obj)
  }

  stop(
    "Unsupported class: ",
    paste(class(obj), collapse = ";")
  )
}

read_ok <- object_audit[
  object_audit$read_status == "READ_OK" &
  object_audit$rds_path != "",
  ,
  drop = FALSE
]

read_ok <- read_ok[
  !duplicated(read_ok$dataset),
  ,
  drop = FALSE
]

cat("\n===== ONE-OBJECT METADATA AUDIT PER READABLE DATASET =====\n")

if (!nrow(read_ok)) {
  cat("NO_READABLE_OBJECT\n")
} else {
  for (i in seq_len(nrow(read_ok))) {
    cat(
      "\n--- DATASET=",
      read_ok$dataset[i],
      " OBJECT=",
      read_ok$object_name[i],
      " ---\n",
      sep = ""
    )

    obj <- tryCatch(
      readRDS(read_ok$rds_path[i]),
      error = function(e) e
    )

    if (inherits(obj, "error")) {
      cat(
        "READ_ERROR=",
        conditionMessage(obj),
        "\n",
        sep = ""
      )
      next
    }

    meta <- tryCatch(
      extract_metadata(obj),
      error = function(e) e
    )

    if (inherits(meta, "error")) {
      cat(
        "METADATA_ERROR=",
        conditionMessage(meta),
        "\n",
        sep = ""
      )
      rm(obj)
      gc(verbose = FALSE)
      next
    }

    meta <- as.data.frame(
      meta,
      check.names = FALSE
    )

    cat(
      "class=",
      paste(class(obj), collapse = ";"),
      "\n",
      sep = ""
    )

    cat(
      "cells=",
      nrow(meta),
      "\n",
      sep = ""
    )

    cat("metadata_columns:\n")

    for (j in seq_along(names(meta))) {
      values <- clean_text(meta[[j]])
      values <- values[values != ""]

      unique_n <- length(unique(values))

      cat(
        sprintf(
          "%03d\t%s\tclass=%s\tunique=%d\n",
          j,
          names(meta)[j],
          paste(class(meta[[j]]), collapse = ";"),
          unique_n
        )
      )
    }

    candidate_columns <- names(meta)[
      grepl(
        paste(
          c(
            "cell",
            "type",
            "annot",
            "cluster",
            "ident",
            "subtype",
            "lineage",
            "seurat",
            "snn"
          ),
          collapse = "|"
        ),
        names(meta),
        ignore.case = TRUE
      )
    ]

    cat("candidate_categorical_fields:\n")

    if (!length(candidate_columns)) {
      cat("NONE\n")
    } else {
      for (column_name in candidate_columns) {
        values <- clean_text(meta[[column_name]])
        values <- values[values != ""]
        unique_values <- sort(unique(values))

        cat(
          column_name,
          "\tunique=",
          length(unique_values),
          "\tvalues=",
          paste(
            head(unique_values, 30L),
            collapse = ";"
          ),
          if (length(unique_values) > 30L) {
            ";<TRUNCATED>"
          } else {
            ""
          },
          "\n",
          sep = ""
        )
      }
    }

    rm(obj, meta)
    gc(verbose = FALSE)
  }
}

cat("\n===== DIAGNOSTIC SUMMARY =====\n")
cat(
  "manifest_objects=",
  nrow(manifest),
  "\n",
  sep = ""
)
cat(
  "resolved_paths=",
  sum(path_audit$rds_path != ""),
  "\n",
  sep = ""
)
cat(
  "unresolved_paths=",
  sum(path_audit$rds_path == ""),
  "\n",
  sep = ""
)
cat(
  "read_ok_objects=",
  sum(object_audit$read_status == "READ_OK"),
  "\n",
  sep = ""
)
cat(
  "objects_with_annotation_field=",
  sum(clean_text(object_audit$annotation_field) != ""),
  "\n",
  sep = ""
)
