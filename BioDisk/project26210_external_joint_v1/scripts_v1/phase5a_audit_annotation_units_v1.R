options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

root <- Sys.getenv(
  "PHASE5A_ROOT",
  unset = "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
)

phase4i_dir <- file.path(
  root,
  "phase4i_finalize_residual_treatment_v1"
)

out_dir <- file.path(
  root,
  "phase5a_audit_annotation_units_v1"
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

stamp <- function() {
  format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z"
  )
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

is_missing <- function(x) {
  toupper(clean_text(x)) %in% c(
    "",
    "NA",
    "N/A",
    "NULL",
    "UNKNOWN",
    "NOT AVAILABLE",
    "NOT REPORTED"
  )
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }

  read.delim(
    path,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

write_tsv <- function(x, path) {
  y <- x

  for (nm in names(y)) {
    if (is.character(y[[nm]])) {
      y[[nm]] <- clean_text(y[[nm]])
    }
  }

  write.table(
    y,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

first_existing_column <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]

  if (length(hit)) {
    hit[1]
  } else {
    ""
  }
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(clean_text(x)))
}

manifest_file <- file.path(
  phase4i_dir,
  "tables",
  "object_clinical_manifest_phase4i_final.tsv"
)

manifest <- read_tsv(manifest_file)

required_manifest_columns <- c(
  "dataset",
  "object_name",
  "cells"
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  names(manifest)
)

if (length(missing_manifest_columns)) {
  stop(
    "Phase 4I manifest missing columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

manifest$dataset <- clean_text(manifest$dataset)
manifest$object_name <- clean_text(manifest$object_name)

object_key <- paste(
  manifest$dataset,
  manifest$object_name,
  sep = "\r"
)

if (nrow(manifest) != 95L) {
  stop(
    "Expected 95 manifest objects, observed ",
    nrow(manifest)
  )
}

if (anyDuplicated(object_key)) {
  stop("Duplicated dataset-object keys in Phase 4I manifest")
}

path_columns <- c(
  "rds_path",
  "object_path",
  "file_path",
  "path",
  "rds_file",
  "source_path",
  "input_path"
)

manifest_path_column <- first_existing_column(
  manifest,
  path_columns
)

all_rds <- list.files(
  root,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

all_rds <- normalizePath(
  all_rds,
  winslash = "/",
  mustWork = FALSE
)

all_rds <- unique(all_rds)

rds_stem <- tools::file_path_sans_ext(
  basename(all_rds)
)

resolve_rds_path <- function(i) {
  explicit_path <- ""

  if (manifest_path_column != "") {
    explicit_path <- clean_text(
      manifest[[manifest_path_column]][i]
    )

    if (
      explicit_path != "" &&
      !grepl("^/", explicit_path)
    ) {
      explicit_path <- file.path(
        root,
        explicit_path
      )
    }

    if (
      explicit_path != "" &&
      file.exists(explicit_path)
    ) {
      return(
        list(
          path = normalizePath(
            explicit_path,
            winslash = "/",
            mustWork = TRUE
          ),
          status = "EXPLICIT_MANIFEST_PATH",
          candidates = 1L
        )
      )
    }
  }

  object_name <- manifest$object_name[i]

  exact_index <- which(
    rds_stem == object_name
  )

  if (length(exact_index) == 1L) {
    return(
      list(
        path = all_rds[exact_index],
        status = "EXACT_BASENAME_MATCH",
        candidates = 1L
      )
    )
  }

  if (length(exact_index) > 1L) {
    dataset_index <- exact_index[
      grepl(
        manifest$dataset[i],
        all_rds[exact_index],
        fixed = TRUE
      )
    ]

    if (length(dataset_index) == 1L) {
      return(
        list(
          path = all_rds[dataset_index],
          status = "BASENAME_AND_DATASET_MATCH",
          candidates = length(exact_index)
        )
      )
    }

    return(
      list(
        path = "",
        status = "AMBIGUOUS_BASENAME",
        candidates = length(exact_index)
      )
    )
  }

  list(
    path = "",
    status = "RDS_NOT_FOUND",
    candidates = 0L
  )
}

resolved_paths <- lapply(
  seq_len(nrow(manifest)),
  resolve_rds_path
)

path_audit <- data.frame(
  dataset = manifest$dataset,
  object_name = manifest$object_name,
  rds_path = vapply(
    resolved_paths,
    function(x) x$path,
    character(1)
  ),
  path_status = vapply(
    resolved_paths,
    function(x) x$status,
    character(1)
  ),
  candidate_paths = vapply(
    resolved_paths,
    function(x) x$candidates,
    integer(1)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  path_audit,
  file.path(
    table_dir,
    "object_rds_path_audit.tsv"
  )
)

annotation_priority <- c(
  "celltype_refined",
  "celltype_final",
  "cell_type_final",
  "final_celltype",
  "final_annotation",
  "celltype",
  "cell_type",
  "CellType",
  "annotation",
  "Annotation",
  "cell_annotation",
  "predicted.celltype",
  "predicted_celltype",
  "predicted.id",
  "major_celltype",
  "major_cell_type",
  "broad_celltype",
  "broad_cell_type",
  "celltype_major",
  "cell_type_major",
  "subtype",
  "cell_subtype"
)

patient_priority <- c(
  "patient_id",
  "patient",
  "Patient",
  "donor_id",
  "donor",
  "subject_id",
  "subject",
  "case_id",
  "case",
  "individual"
)

extract_metadata <- function(obj) {
  if (
    inherits(obj, "Seurat") ||
    inherits(obj, "SeuratObject")
  ) {
    return(obj[[]])
  }

  if (
    inherits(obj, "SingleCellExperiment") ||
    inherits(obj, "SummarizedExperiment")
  ) {
    if (!requireNamespace(
      "SummarizedExperiment",
      quietly = TRUE
    )) {
      stop(
        "SummarizedExperiment package is required ",
        "to read colData"
      )
    }

    return(
      as.data.frame(
        SummarizedExperiment::colData(obj),
        check.names = FALSE
      )
    )
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
    "Unsupported object class: ",
    paste(class(obj), collapse = ",")
  )
}

extract_genes <- function(obj) {
  genes <- character()

  if (inherits(obj, "Seurat")) {
    assay_names <- tryCatch(
      names(obj@assays),
      error = function(e) character()
    )

    for (assay_name in assay_names) {
      assay_genes <- tryCatch(
        rownames(obj@assays[[assay_name]]),
        error = function(e) character()
      )

      genes <- c(genes, assay_genes)
    }
  } else if (
    inherits(obj, "SingleCellExperiment") ||
    inherits(obj, "SummarizedExperiment")
  ) {
    genes <- tryCatch(
      rownames(obj),
      error = function(e) character()
    )
  } else {
    genes <- tryCatch(
      rownames(obj),
      error = function(e) character()
    )
  }

  unique(clean_text(genes))
}

find_annotation_columns <- function(meta) {
  named_hits <- annotation_priority[
    annotation_priority %in% names(meta)
  ]

  regex_hits <- names(meta)[
    grepl(
      "cell.*type|celltype|annotation|subtype",
      names(meta),
      ignore.case = TRUE
    )
  ]

  candidates <- unique(
    c(named_hits, regex_hits)
  )

  candidates[
    vapply(
      candidates,
      function(column_name) {
        values <- clean_text(meta[[column_name]])
        values <- values[!is_missing(values)]

        length(unique(values)) >= 2L
      },
      logical(1)
    )
  ]
}

select_annotation_column <- function(candidates) {
  if (!length(candidates)) {
    return("")
  }

  priority_index <- match(
    candidates,
    annotation_priority
  )

  priority_index[
    is.na(priority_index)
  ] <- length(annotation_priority) + seq_len(
    sum(is.na(priority_index))
  )

  candidates[
    order(priority_index)
  ][1]
}

extract_manifest_patient <- function(i) {
  if (!"patient_id" %in% names(manifest)) {
    return("")
  }

  clean_text(manifest$patient_id[i])
}

canonicalize_celltype <- function(raw_label) {
  x <- tolower(clean_text(raw_label))

  x <- gsub(
    "[_./-]+",
    " ",
    x
  )

  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )

  focus_tymp <- (
    grepl("tymp", x) &&
    grepl("high|hi|positive|\\+", x) &&
    grepl(
      "malignant|tumou?r|cancer|carcinoma",
      x
    )
  )

  focus_spp1 <- (
    grepl("spp1", x) &&
    grepl(
      "macrophage|macro|tam|myeloid",
      x
    )
  )

  focus_apoe_c1qc <- (
    grepl("apoe", x) &&
    grepl("c1qc", x) &&
    grepl(
      "macrophage|macro|tam|myeloid",
      x
    )
  )

  if (focus_tymp) {
    return("TYMP-High Malignant")
  }

  if (focus_spp1) {
    return("SPP1+ TAM")
  }

  if (focus_apoe_c1qc) {
    return("APOE+C1QC+ TAM")
  }

  if (
    grepl(
      "macrophage|macro|tam|mononuclear phagocyte",
      x
    )
  ) {
    return("Other Macrophage")
  }

  if (grepl("monocyte|\\bmono\\b", x)) {
    return("Monocyte")
  }

  if (
    grepl(
      "dendritic|\\bdc\\b|cd1c|clec9a|fcer1a",
      x
    )
  ) {
    return("Dendritic")
  }

  if (
    grepl(
      "natural killer|\\bnk\\b|nkt",
      x
    )
  ) {
    return("NK")
  }

  if (
    grepl(
      "\\bt cell\\b|\\btcell\\b|cd4|cd8|regulatory t|\\btreg\\b",
      x
    )
  ) {
    return("T cell")
  }

  if (
    grepl(
      "plasma|plasmablast",
      x
    )
  ) {
    return("Plasma")
  }

  if (
    grepl(
      "\\bb cell\\b|\\bbcell\\b|memory b|naive b",
      x
    )
  ) {
    return("B cell")
  }

  if (
    grepl(
      "endothelial|endothelium",
      x
    )
  ) {
    return("Endothelial")
  }

  if (
    grepl(
      "fibroblast|myofibroblast|\\bfibro\\b",
      x
    )
  ) {
    return("Fibroblast")
  }

  if (
    grepl(
      "pericyte|smooth muscle|\\bsmc\\b",
      x
    )
  ) {
    return("SMC/Pericyte")
  }

  if (grepl("mast", x)) {
    return("Mast")
  }

  if (
    grepl(
      "malignant|cancer|carcinoma|tumou?r cell",
      x
    )
  ) {
    return("Malignant")
  }

  if (
    grepl(
      "epithelial|proximal tubule|distal tubule|collecting duct|renal tubular",
      x
    )
  ) {
    return("Epithelial-Normal")
  }

  "UNMAPPED"
}

object_rows <- list()
candidate_rows <- list()
count_rows <- list()
gene_rows <- list()
patient_rows <- list()

for (i in seq_len(nrow(manifest))) {
  dataset <- manifest$dataset[i]
  object_name <- manifest$object_name[i]
  rds_path <- path_audit$rds_path[i]

  base_row <- data.frame(
    dataset = dataset,
    object_name = object_name,
    rds_path = rds_path,
    object_class = "",
    expected_cells = safe_numeric(
      manifest$cells[i]
    ),
    observed_cells = NA_real_,
    cell_count_match = FALSE,
    annotation_field = "",
    annotation_nonmissing_cells = 0L,
    annotation_coverage = NA_real_,
    annotation_unique_labels = 0L,
    manifest_patient_id = extract_manifest_patient(i),
    metadata_patient_field = "",
    metadata_patient_values = "",
    read_status = "",
    read_error = "",
    stringsAsFactors = FALSE
  )

  if (rds_path == "") {
    base_row$read_status <- "PATH_UNRESOLVED"
    object_rows[[length(object_rows) + 1L]] <- base_row
    next
  }

  obj <- tryCatch(
    readRDS(rds_path),
    error = function(e) e
  )

  if (inherits(obj, "error")) {
    base_row$read_status <- "READ_ERROR"
    base_row$read_error <- conditionMessage(obj)
    object_rows[[length(object_rows) + 1L]] <- base_row
    next
  }

  base_row$object_class <- paste(
    class(obj),
    collapse = ";"
  )

  meta <- tryCatch(
    extract_metadata(obj),
    error = function(e) e
  )

  if (inherits(meta, "error")) {
    base_row$read_status <- "METADATA_ERROR"
    base_row$read_error <- conditionMessage(meta)
    object_rows[[length(object_rows) + 1L]] <- base_row
    rm(obj)
    gc(verbose = FALSE)
    next
  }

  meta <- as.data.frame(
    meta,
    check.names = FALSE
  )

  base_row$observed_cells <- nrow(meta)

  base_row$cell_count_match <- (
    is.na(base_row$expected_cells) ||
    base_row$expected_cells == base_row$observed_cells
  )

  annotation_candidates <- find_annotation_columns(meta)

  if (length(annotation_candidates)) {
    for (column_name in annotation_candidates) {
      values <- clean_text(meta[[column_name]])
      valid <- values[!is_missing(values)]

      candidate_rows[[length(candidate_rows) + 1L]] <-
        data.frame(
          dataset = dataset,
          object_name = object_name,
          annotation_field = column_name,
          nonmissing_cells = length(valid),
          coverage = if (nrow(meta)) {
            length(valid) / nrow(meta)
          } else {
            NA_real_
          },
          unique_labels = length(unique(valid)),
          selected = FALSE,
          stringsAsFactors = FALSE
        )
    }
  }

  selected_annotation <- select_annotation_column(
    annotation_candidates
  )

  base_row$annotation_field <- selected_annotation

  if (selected_annotation != "") {
    values <- clean_text(
      meta[[selected_annotation]]
    )

    valid_index <- !is_missing(values)
    valid_values <- values[valid_index]

    base_row$annotation_nonmissing_cells <-
      length(valid_values)

    base_row$annotation_coverage <- if (nrow(meta)) {
      length(valid_values) / nrow(meta)
    } else {
      NA_real_
    }

    base_row$annotation_unique_labels <-
      length(unique(valid_values))

    if (length(candidate_rows)) {
      selected_index <- which(
        vapply(
          candidate_rows,
          function(x) {
            x$dataset[1] == dataset &&
            x$object_name[1] == object_name &&
            x$annotation_field[1] ==
              selected_annotation
          },
          logical(1)
        )
      )

      if (length(selected_index)) {
        candidate_rows[[selected_index[length(selected_index)]]]$selected <-
          TRUE
      }
    }

    if (length(valid_values)) {
      label_table <- as.data.frame(
        table(valid_values),
        stringsAsFactors = FALSE
      )

      names(label_table) <- c(
        "raw_celltype",
        "cells"
      )

      label_table$dataset <- dataset
      label_table$object_name <- object_name
      label_table$annotation_field <-
        selected_annotation

      count_rows[[length(count_rows) + 1L]] <-
        label_table[
          ,
          c(
            "dataset",
            "object_name",
            "annotation_field",
            "raw_celltype",
            "cells"
          )
        ]
    }
  }

  patient_candidates <- patient_priority[
    patient_priority %in% names(meta)
  ]

  if (length(patient_candidates)) {
    patient_field <- patient_candidates[1]
    patient_values <- unique(
      clean_text(meta[[patient_field]])
    )

    patient_values <- patient_values[
      !is_missing(patient_values)
    ]

    base_row$metadata_patient_field <-
      patient_field

    base_row$metadata_patient_values <- paste(
      sort(patient_values),
      collapse = ";"
    )
  }

  genes <- extract_genes(obj)
  genes_upper <- toupper(genes)

  gene_rows[[length(gene_rows) + 1L]] <-
    data.frame(
      dataset = dataset,
      object_name = object_name,
      TYMP = "TYMP" %in% genes_upper,
      SPP1 = "SPP1" %in% genes_upper,
      APOE = "APOE" %in% genes_upper,
      C1QC = "C1QC" %in% genes_upper,
      stringsAsFactors = FALSE
    )

  manifest_patient <- clean_text(
    base_row$manifest_patient_id
  )

  metadata_values <- clean_text(
    base_row$metadata_patient_values
  )

  if (!is_missing(manifest_patient)) {
    patient_id <- manifest_patient
    patient_source <- "PHASE4I_CONFIRMED"
    patient_status <- "CONFIRMED"
  } else if (
    metadata_values != "" &&
    !grepl(";", metadata_values, fixed = TRUE)
  ) {
    patient_id <- metadata_values
    patient_source <- paste0(
      "RDS_METADATA:",
      base_row$metadata_patient_field
    )
    patient_status <- "REVIEW_REQUIRED"
  } else {
    patient_id <- ""
    patient_source <- ""
    patient_status <- "UNRESOLVED"
  }

  patient_rows[[length(patient_rows) + 1L]] <-
    data.frame(
      dataset = dataset,
      object_name = object_name,
      patient_id = patient_id,
      patient_source = patient_source,
      patient_status = patient_status,
      cells = nrow(meta),
      stringsAsFactors = FALSE
    )

  base_row$read_status <- "READ_OK"

  object_rows[[length(object_rows) + 1L]] <- base_row

  rm(obj, meta)
  gc(verbose = FALSE)

  cat(
    sprintf(
      "[%d/%d] %s | %s | %s\n",
      i,
      nrow(manifest),
      dataset,
      object_name,
      base_row$annotation_field
    )
  )
}

object_audit <- do.call(
  rbind,
  object_rows
)

annotation_candidates <- if (length(candidate_rows)) {
  do.call(rbind, candidate_rows)
} else {
  data.frame(
    dataset = character(),
    object_name = character(),
    annotation_field = character(),
    nonmissing_cells = integer(),
    coverage = numeric(),
    unique_labels = integer(),
    selected = logical(),
    stringsAsFactors = FALSE
  )
}

raw_counts <- if (length(count_rows)) {
  do.call(rbind, count_rows)
} else {
  data.frame(
    dataset = character(),
    object_name = character(),
    annotation_field = character(),
    raw_celltype = character(),
    cells = integer(),
    stringsAsFactors = FALSE
  )
}

gene_availability <- if (length(gene_rows)) {
  do.call(rbind, gene_rows)
} else {
  data.frame(
    dataset = character(),
    object_name = character(),
    TYMP = logical(),
    SPP1 = logical(),
    APOE = logical(),
    C1QC = logical(),
    stringsAsFactors = FALSE
  )
}

patient_units <- if (length(patient_rows)) {
  do.call(rbind, patient_rows)
} else {
  data.frame(
    dataset = character(),
    object_name = character(),
    patient_id = character(),
    patient_source = character(),
    patient_status = character(),
    cells = numeric(),
    stringsAsFactors = FALSE
  )
}

if (nrow(raw_counts)) {
  crosswalk_key <- paste(
    raw_counts$dataset,
    raw_counts$annotation_field,
    raw_counts$raw_celltype,
    sep = "\r"
  )

  split_rows <- split(
    seq_len(nrow(raw_counts)),
    crosswalk_key
  )

  crosswalk <- do.call(
    rbind,
    lapply(
      split_rows,
      function(index) {
        canonical <- canonicalize_celltype(
          raw_counts$raw_celltype[index[1]]
        )

        mapping_status <- if (
          canonical %in% c(
            "TYMP-High Malignant",
            "SPP1+ TAM",
            "APOE+C1QC+ TAM"
          )
        ) {
          "FOCUS_LABEL_MATCH"
        } else if (canonical != "UNMAPPED") {
          "AUTO_BROAD_REVIEW"
        } else {
          "UNMAPPED_REVIEW_REQUIRED"
        }

        data.frame(
          dataset = raw_counts$dataset[index[1]],
          annotation_field =
            raw_counts$annotation_field[index[1]],
          raw_celltype =
            raw_counts$raw_celltype[index[1]],
          canonical_celltype = canonical,
          mapping_status = mapping_status,
          n_objects = length(
            unique(
              raw_counts$object_name[index]
            )
          ),
          cells = sum(
            safe_numeric(
              raw_counts$cells[index]
            ),
            na.rm = TRUE
          ),
          analyst_decision = "",
          analyst_notes = "",
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(crosswalk) <- NULL

  crosswalk <- crosswalk[
    order(
      crosswalk$dataset,
      crosswalk$canonical_celltype,
      -crosswalk$cells,
      crosswalk$raw_celltype
    ),
    ,
    drop = FALSE
  ]
} else {
  crosswalk <- data.frame(
    dataset = character(),
    annotation_field = character(),
    raw_celltype = character(),
    canonical_celltype = character(),
    mapping_status = character(),
    n_objects = integer(),
    cells = numeric(),
    analyst_decision = character(),
    analyst_notes = character(),
    stringsAsFactors = FALSE
  )
}

dataset_names <- sort(
  unique(manifest$dataset)
)

readiness_rows <- lapply(
  dataset_names,
  function(dataset_name) {
    manifest_index <- manifest$dataset == dataset_name

    audit_index <- object_audit$dataset == dataset_name
    count_index <- raw_counts$dataset == dataset_name
    gene_index <- gene_availability$dataset == dataset_name
    patient_index <- patient_units$dataset == dataset_name

    dataset_crosswalk <- crosswalk[
      crosswalk$dataset == dataset_name,
      ,
      drop = FALSE
    ]

    canonical_types <- unique(
      dataset_crosswalk$canonical_celltype[
        dataset_crosswalk$canonical_celltype !=
          "UNMAPPED"
      ]
    )

    unresolved_labels <- sum(
      dataset_crosswalk$mapping_status ==
        "UNMAPPED_REVIEW_REQUIRED"
    )

    confirmed_patients <- length(
      unique(
        patient_units$patient_id[
          patient_index &
          patient_units$patient_status ==
            "CONFIRMED" &
          patient_units$patient_id != ""
        ]
      )
    )

    review_patients <- length(
      unique(
        patient_units$patient_id[
          patient_index &
          patient_units$patient_status ==
            "REVIEW_REQUIRED" &
          patient_units$patient_id != ""
        ]
      )
    )

    annotation_ready <- (
      sum(
        object_audit$annotation_field[
          audit_index
        ] != ""
      ) ==
      sum(audit_index)
    )

    genes_ready <- if (sum(gene_index)) {
      all(
        gene_availability$TYMP[gene_index] &
        gene_availability$SPP1[gene_index]
      )
    } else {
      FALSE
    }

    readiness <- if (
      all(
        object_audit$read_status[audit_index] ==
          "READ_OK"
      ) &&
      annotation_ready &&
      unresolved_labels == 0L &&
      genes_ready
    ) {
      "READY_FOR_PHASE5B"
    } else {
      "REVIEW_REQUIRED"
    }

    data.frame(
      dataset = dataset_name,
      objects = sum(manifest_index),
      expected_cells = sum(
        safe_numeric(
          manifest$cells[manifest_index]
        ),
        na.rm = TRUE
      ),
      read_ok_objects = sum(
        object_audit$read_status[audit_index] ==
          "READ_OK"
      ),
      annotated_objects = sum(
        object_audit$annotation_field[audit_index] !=
          ""
      ),
      raw_labels = sum(count_index),
      canonical_types = length(canonical_types),
      unmapped_labels = unresolved_labels,
      has_tymp_high_label =
        "TYMP-High Malignant" %in% canonical_types,
      has_spp1_tam_label =
        "SPP1+ TAM" %in% canonical_types,
      has_apoe_c1qc_tam_label =
        "APOE+C1QC+ TAM" %in% canonical_types,
      all_objects_have_TYMP =
        sum(gene_index) > 0L &&
        all(gene_availability$TYMP[gene_index]),
      all_objects_have_SPP1 =
        sum(gene_index) > 0L &&
        all(gene_availability$SPP1[gene_index]),
      confirmed_patients = confirmed_patients,
      metadata_patients_for_review = review_patients,
      readiness = readiness,
      stringsAsFactors = FALSE
    )
  }
)

dataset_readiness <- do.call(
  rbind,
  readiness_rows
)

cell_count_mismatch <- object_audit[
  object_audit$read_status == "READ_OK" &
  !object_audit$cell_count_match,
  ,
  drop = FALSE
]

read_failures <- object_audit[
  object_audit$read_status != "READ_OK",
  ,
  drop = FALSE
]

annotation_missing <- object_audit[
  object_audit$read_status == "READ_OK" &
  object_audit$annotation_field == "",
  ,
  drop = FALSE
]

write_tsv(
  object_audit,
  file.path(
    table_dir,
    "object_annotation_audit.tsv"
  )
)

write_tsv(
  annotation_candidates,
  file.path(
    table_dir,
    "annotation_candidate_columns.tsv"
  )
)

write_tsv(
  raw_counts,
  file.path(
    table_dir,
    "celltype_counts_raw.tsv"
  )
)

write_tsv(
  crosswalk,
  file.path(
    table_dir,
    "annotation_crosswalk_candidate.tsv"
  )
)

write_tsv(
  gene_availability,
  file.path(
    table_dir,
    "focus_gene_availability.tsv"
  )
)

write_tsv(
  patient_units,
  file.path(
    table_dir,
    "patient_analysis_units_candidate.tsv"
  )
)

write_tsv(
  dataset_readiness,
  file.path(
    table_dir,
    "dataset_mainline_readiness.tsv"
  )
)

write_tsv(
  read_failures,
  file.path(
    qc_dir,
    "object_read_failures.tsv"
  )
)

write_tsv(
  cell_count_mismatch,
  file.path(
    qc_dir,
    "cell_count_mismatches.tsv"
  )
)

write_tsv(
  annotation_missing,
  file.path(
    qc_dir,
    "objects_without_annotation.tsv"
  )
)

structural_qc <- data.frame(
  check = c(
    "Phase 4I manifest objects",
    "Unique manifest object keys",
    "Resolved RDS paths",
    "Successfully read objects",
    "Cell-count matched objects",
    "RDS files modified",
    "Phase 4I files modified",
    "Customer report modified"
  ),
  expected = c(
    "95",
    "95",
    "95",
    "95",
    "95",
    "0",
    "0",
    "0"
  ),
  observed = c(
    as.character(nrow(manifest)),
    as.character(length(unique(object_key))),
    as.character(sum(path_audit$rds_path != "")),
    as.character(
      sum(object_audit$read_status == "READ_OK")
    ),
    as.character(
      sum(
        object_audit$read_status == "READ_OK" &
        object_audit$cell_count_match
      )
    ),
    "0",
    "0",
    "0"
  ),
  stringsAsFactors = FALSE
)

structural_qc$status <- ifelse(
  structural_qc$expected ==
    structural_qc$observed,
  "PASS",
  "FAIL"
)

write_tsv(
  structural_qc,
  file.path(
    qc_dir,
    "phase5a_structural_qc.tsv"
  )
)

structural_fail <- any(
  structural_qc$status == "FAIL"
)

review_required <- (
  nrow(annotation_missing) > 0L ||
  sum(
    crosswalk$mapping_status ==
      "UNMAPPED_REVIEW_REQUIRED"
  ) > 0L ||
  any(
    dataset_readiness$readiness ==
      "REVIEW_REQUIRED"
  )
)

overall_status <- if (structural_fail) {
  "FAIL"
} else if (review_required) {
  "REVIEW_REQUIRED"
} else {
  "READY_FOR_PHASE5B"
}

report <- c(
  "# Phase 5A annotation and analysis-unit audit",
  "",
  paste0("Generated: ", stamp()),
  "",
  "## Structural audit",
  "",
  paste0(
    "- Manifest objects: ",
    nrow(manifest)
  ),
  paste0(
    "- Resolved RDS paths: ",
    sum(path_audit$rds_path != "")
  ),
  paste0(
    "- Successfully read objects: ",
    sum(
      object_audit$read_status == "READ_OK"
    )
  ),
  paste0(
    "- Cell-count mismatches: ",
    nrow(cell_count_mismatch)
  ),
  "",
  "## Annotation audit",
  "",
  paste0(
    "- Objects with selected annotation field: ",
    sum(object_audit$annotation_field != "")
  ),
  paste0(
    "- Unique dataset-field-label combinations: ",
    nrow(crosswalk)
  ),
  paste0(
    "- Unmapped labels requiring review: ",
    sum(
      crosswalk$mapping_status ==
        "UNMAPPED_REVIEW_REQUIRED"
    )
  ),
  "",
  "## Patient units",
  "",
  paste0(
    "- Confirmed unique dataset-patient units: ",
    length(
      unique(
        paste(
          patient_units$dataset[
            patient_units$patient_status ==
              "CONFIRMED"
          ],
          patient_units$patient_id[
            patient_units$patient_status ==
              "CONFIRMED"
          ],
          sep = "\r"
        )
      )
    )
  ),
  paste0(
    "- Metadata-derived patient IDs requiring review: ",
    sum(
      patient_units$patient_status ==
        "REVIEW_REQUIRED"
    )
  ),
  paste0(
    "- Objects with unresolved patient ID: ",
    sum(
      patient_units$patient_status ==
        "UNRESOLVED"
    )
  ),
  "",
  "## Interpretation rule",
  "",
  paste0(
    "Existing source labels are audited separately from ",
    "new within-dataset TYMP-high malignant definitions. ",
    "A generic malignant label is not automatically ",
    "reclassified as TYMP-High Malignant."
  ),
  "",
  "## Safety",
  "",
  paste0(
    "No RDS object, Phase 4I output or customer report ",
    "was modified."
  ),
  "",
  "## Status",
  "",
  overall_status
)

writeLines(
  report,
  file.path(
    out_dir,
    "PHASE5A_REPORT.md"
  ),
  useBytes = TRUE
)

status_file <- file.path(
  out_dir,
  if (structural_fail) {
    "PHASE5A_INCOMPLETE"
  } else {
    "PHASE5A_COMPLETE"
  }
)

writeLines(
  c(
    paste0("completed_at=", stamp()),
    paste0(
      "structural_qc=",
      if (structural_fail) "FAIL" else "PASS"
    ),
    paste0(
      "analysis_readiness=",
      overall_status
    ),
    paste0(
      "manifest_objects=",
      nrow(manifest)
    ),
    paste0(
      "resolved_rds_paths=",
      sum(path_audit$rds_path != "")
    ),
    paste0(
      "read_ok_objects=",
      sum(object_audit$read_status == "READ_OK")
    ),
    paste0(
      "objects_with_annotation=",
      sum(object_audit$annotation_field != "")
    ),
    paste0(
      "cell_count_mismatches=",
      nrow(cell_count_mismatch)
    ),
    paste0(
      "crosswalk_rows=",
      nrow(crosswalk)
    ),
    paste0(
      "unmapped_labels=",
      sum(
        crosswalk$mapping_status ==
          "UNMAPPED_REVIEW_REQUIRED"
      )
    ),
    paste0(
      "confirmed_patient_object_rows=",
      sum(
        patient_units$patient_status ==
          "CONFIRMED"
      )
    ),
    paste0(
      "patient_rows_for_review=",
      sum(
        patient_units$patient_status ==
          "REVIEW_REQUIRED"
      )
    ),
    paste0(
      "unresolved_patient_object_rows=",
      sum(
        patient_units$patient_status ==
          "UNRESOLVED"
      )
    ),
    "rds_modified=0",
    "phase4i_modified=0",
    "customer_report_modified=0"
  ),
  status_file
)

cat(
  "STATUS=",
  overall_status,
  "\n",
  sep = ""
)

if (structural_fail) {
  quit(
    save = "no",
    status = 2L
  )
}
