options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

root <- Sys.getenv(
  "PHASE4H_ROOT",
  unset = "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
)

phase4e_dir <- file.path(
  root,
  "phase4e_clinical_metadata_harmonization_v1"
)

phase4g_dir <- file.path(
  root,
  "phase4g_resolve_modality_and_response_v1"
)

out_dir <- file.path(
  root,
  "phase4h_finalize_clinical_units_v1"
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

collapse_unique <- function(
  x,
  sep = " | ",
  default = ""
) {
  y <- clean_text(x)
  y <- y[!is_missing(y)]
  y <- unique(y)

  if (!length(y)) {
    return(default)
  }

  paste(y, collapse = sep)
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

require_columns <- function(x, columns, label) {
  missing_columns <- setdiff(columns, names(x))

  if (length(missing_columns)) {
    stop(
      label,
      " missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

extract_group <- function(x, pattern) {
  matches <- regexec(
    pattern,
    clean_text(x),
    ignore.case = TRUE,
    perl = TRUE
  )

  values <- regmatches(
    clean_text(x),
    matches
  )

  vapply(
    values,
    function(z) {
      if (length(z) >= 2L) {
        z[2]
      } else {
        ""
      }
    },
    character(1)
  )
}

make_updates <- function(
  dataset,
  object_name,
  field,
  value,
  metadata_source,
  source_reference,
  notes
) {
  data.frame(
    dataset = clean_text(dataset),
    object_name = clean_text(object_name),
    field = clean_text(field),
    value = clean_text(value),
    metadata_source = clean_text(metadata_source),
    source_reference = clean_text(source_reference),
    notes = clean_text(notes),
    stringsAsFactors = FALSE
  )
}

manifest_file <- file.path(
  phase4e_dir,
  "tables",
  "object_clinical_manifest.tsv"
)

resolved_file <- file.path(
  phase4g_dir,
  "tables",
  "object_geo_resolved_mapping.tsv"
)

gse207493_patient_file <- file.path(
  phase4g_dir,
  "tables",
  "gse207493_patient_id_evidence.tsv"
)

response_file <- file.path(
  phase4g_dir,
  "tables",
  "response_candidate_audit.tsv"
)

treatment_file <- file.path(
  phase4g_dir,
  "tables",
  "treatment_candidate_audit.tsv"
)

manifest <- read_tsv(manifest_file)
resolved <- read_tsv(resolved_file)
gse207493_patient <- read_tsv(
  gse207493_patient_file
)
response_audit <- read_tsv(response_file)
treatment_audit <- read_tsv(treatment_file)

require_columns(
  manifest,
  c(
    "dataset",
    "object_name",
    "source_file",
    "cells"
  ),
  "Phase 4E manifest"
)

require_columns(
  resolved,
  c(
    "dataset",
    "object_name",
    "resolved_match_status",
    "resolved_gsm",
    "resolved_geo_title"
  ),
  "Phase 4G resolved mapping"
)

require_columns(
  gse207493_patient,
  c(
    "dataset",
    "object_name",
    "patient_id",
    "source_reference"
  ),
  "GSE207493 patient evidence"
)

require_columns(
  response_audit,
  c(
    "dataset",
    "object_name",
    "raw_value",
    "source_reference"
  ),
  "Response audit"
)

require_columns(
  treatment_audit,
  c(
    "dataset",
    "object_name",
    "raw_key",
    "raw_value",
    "audit_decision",
    "source_reference"
  ),
  "Treatment audit"
)

manifest$key <- paste(
  manifest$dataset,
  manifest$object_name,
  sep = "\r"
)

resolved$key <- paste(
  resolved$dataset,
  resolved$object_name,
  sep = "\r"
)

if (nrow(manifest) != 95L) {
  stop(
    "Expected 95 manifest objects, observed ",
    nrow(manifest)
  )
}

if (anyDuplicated(manifest$key)) {
  stop("Duplicated dataset/object_name in Phase 4E manifest")
}

if (nrow(resolved) != 95L) {
  stop(
    "Expected 95 resolved mappings, observed ",
    nrow(resolved)
  )
}

if (anyDuplicated(resolved$key)) {
  stop("Duplicated dataset/object_name in resolved mapping")
}

if (
  sum(
    resolved$resolved_match_status == "MATCHED"
  ) != 95L
) {
  stop("Not all 95 objects have MATCHED status")
}

if (!setequal(manifest$key, resolved$key)) {
  stop(
    "Phase 4E manifest objects do not match ",
    "Phase 4G resolved objects"
  )
}

gse207493_patient$patient_id <- toupper(
  clean_text(gse207493_patient$patient_id)
)

if (nrow(gse207493_patient) != 19L) {
  stop(
    "Expected 19 GSE207493 patient rows, observed ",
    nrow(gse207493_patient)
  )
}

if (
  length(
    unique(gse207493_patient$patient_id)
  ) != 19L
) {
  stop(
    "Expected 19 unique GSE207493 patient IDs"
  )
}

if (
  !all(
    grepl(
      "^RCC[0-9]+$",
      gse207493_patient$patient_id
    )
  )
) {
  stop("Invalid GSE207493 patient ID format")
}

gse294109 <- resolved[
  resolved$dataset == "GSE294109",
  ,
  drop = FALSE
]

if (nrow(gse294109) != 25L) {
  stop(
    "Expected 25 GSE294109 objects, observed ",
    nrow(gse294109)
  )
}

gse294109$patient_id <- toupper(
  extract_group(
    gse294109$resolved_geo_title,
    "\\b(sc[0-9]+)\\b"
  )
)

if (
  any(
    gse294109$patient_id == ""
  )
) {
  print(
    gse294109[
      gse294109$patient_id == "",
      c(
        "object_name",
        "resolved_gsm",
        "resolved_geo_title"
      ),
      drop = FALSE
    ]
  )

  stop(
    "Failed to extract GSE294109 patient IDs"
  )
}

if (
  length(
    unique(gse294109$patient_id)
  ) != 15L
) {
  stop(
    "Expected 15 unique GSE294109 patients, observed ",
    length(
      unique(
        gse294109$patient_id
      )
    )
  )
}

title_lower <- tolower(
  clean_text(
    gse294109$resolved_geo_title
  )
)

gse294109$treatment <- ""

gse294109$treatment[
  grepl(
    "untreated",
    title_lower,
    fixed = TRUE
  )
] <- "Untreated"

gse294109$treatment[
  grepl(
    "ipilimumab/nivolumab treated",
    title_lower,
    fixed = TRUE
  )
] <- "Ipilimumab + Nivolumab"

gse294109$treatment[
  grepl(
    "cabozantinib treated",
    title_lower,
    fixed = TRUE
  )
] <- "Cabozantinib"

gse294109$treatment[
  grepl(
    "axitinib/pembrolizumab treated",
    title_lower,
    fixed = TRUE
  )
] <- "Axitinib + Pembrolizumab"

gse294109$treatment[
  grepl(
    "lenvatinib/pembrolizumab treated",
    title_lower,
    fixed = TRUE
  )
] <- "Lenvatinib + Pembrolizumab"

if (
  any(
    gse294109$treatment == ""
  )
) {
  print(
    gse294109[
      gse294109$treatment == "",
      c(
        "object_name",
        "resolved_gsm",
        "resolved_geo_title"
      ),
      drop = FALSE
    ]
  )

  stop(
    "Unresolved GSE294109 treatment titles"
  )
}

gse294109$treatment_status <- ifelse(
  gse294109$treatment == "Untreated",
  "Untreated",
  "Treated"
)

response_gse294109 <- response_audit[
  response_audit$dataset == "GSE294109",
  ,
  drop = FALSE
]

if (
  nrow(response_gse294109) != 25L
) {
  stop(
    "Expected 25 GSE294109 response audit rows, observed ",
    nrow(response_gse294109)
  )
}

if (
  anyDuplicated(
    response_gse294109$object_name
  )
) {
  stop(
    "Duplicated GSE294109 response audit objects"
  )
}

response_match <- match(
  gse294109$object_name,
  response_gse294109$object_name
)

if (
  anyNA(response_match)
) {
  stop(
    "GSE294109 response audit does not cover all objects"
  )
}

gse294109$response_raw <- clean_text(
  response_gse294109$raw_value[
    response_match
  ]
)

gse294109$response_reference <- clean_text(
  response_gse294109$source_reference[
    response_match
  ]
)

valid_response <- gse294109$response_raw %in% c(
  "",
  "CR",
  "PR",
  "No_effect"
)

if (!all(valid_response)) {
  print(
    unique(
      gse294109$response_raw[
        !valid_response
      ]
    )
  )

  stop(
    "Unsupported GSE294109 raw response value"
  )
}

if (
  any(
    gse294109$treatment_status == "Treated" &
      gse294109$response_raw == ""
  )
) {
  stop(
    "Treated GSE294109 object lacks pathological response"
  )
}

if (
  any(
    gse294109$treatment_status == "Untreated" &
      gse294109$response_raw != ""
  )
) {
  stop(
    "Untreated GSE294109 object has unexpected response"
  )
}

gse294109$response_group <- ifelse(
  gse294109$response_raw %in% c(
    "CR",
    "PR"
  ),
  "Responder",
  ifelse(
    gse294109$response_raw == "No_effect",
    "Non-responder",
    "Not applicable"
  )
)

gse294109$source_reference <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
  gse294109$resolved_gsm
)

updates <- rbind(
  make_updates(
    dataset = gse207493_patient$dataset,
    object_name =
      gse207493_patient$object_name,
    field = "patient_id",
    value = gse207493_patient$patient_id,
    metadata_source = paste0(
      "NCBI GEO GSE207493 overall design; ",
      "GEO Sample_title"
    ),
    source_reference =
      gse207493_patient$source_reference,
    notes = paste0(
      "Study-specific de-identified patient ID; ",
      "scRNA-seq GSM selected over paired scATAC-seq"
    )
  ),
  make_updates(
    dataset = gse294109$dataset,
    object_name = gse294109$object_name,
    field = "patient_id",
    value = gse294109$patient_id,
    metadata_source =
      "NCBI GEO GSE294109 Sample_title",
    source_reference =
      gse294109$source_reference,
    notes = paste0(
      "Study-specific patient ID extracted from ",
      "the exact matched GEO title"
    )
  ),
  make_updates(
    dataset = gse294109$dataset,
    object_name = gse294109$object_name,
    field = "treatment",
    value = gse294109$treatment,
    metadata_source =
      "NCBI GEO GSE294109 Sample_title",
    source_reference =
      gse294109$source_reference,
    notes = paste0(
      "Exact title-based treatment assignment; ",
      "not interpreted as a longitudinal timepoint"
    )
  ),
  make_updates(
    dataset = gse294109$dataset,
    object_name = gse294109$object_name,
    field = "treatment_status",
    value = gse294109$treatment_status,
    metadata_source =
      "Derived from exact GSE294109 treatment title",
    source_reference =
      gse294109$source_reference,
    notes = "Treated versus untreated exposure status"
  ),
  make_updates(
    dataset = gse294109$dataset,
    object_name = gse294109$object_name,
    field = "response_group",
    value = gse294109$response_group,
    metadata_source = paste0(
      "GEO Sample_characteristics_ch1; ",
      "Phase 4H controlled normalization"
    ),
    source_reference = ifelse(
      gse294109$response_reference == "",
      gse294109$source_reference,
      gse294109$response_reference
    ),
    notes = paste0(
      "CR/PR=Responder; No_effect=Non-responder; ",
      "untreated=Not applicable; No_effect is not ",
      "relabelled as PD"
    )
  )
)

response_present <- gse294109$response_raw != ""

updates <- rbind(
  updates,
  make_updates(
    dataset =
      gse294109$dataset[response_present],
    object_name =
      gse294109$object_name[response_present],
    field = "response_raw",
    value =
      gse294109$response_raw[response_present],
    metadata_source =
      "GEO SOFT: Sample_characteristics_ch1",
    source_reference =
      gse294109$response_reference[response_present],
    notes = paste0(
      "Exact raw pathological response retained; ",
      "No_effect is preserved verbatim"
    )
  )
)

updates$key <- paste(
  updates$dataset,
  updates$object_name,
  updates$field,
  sep = "\r"
)

if (anyDuplicated(updates$key)) {
  duplicated_updates <- updates[
    duplicated(updates$key) |
      duplicated(
        updates$key,
        fromLast = TRUE
      ),
    ,
    drop = FALSE
  ]

  write_tsv(
    duplicated_updates,
    file.path(
      qc_dir,
      "duplicate_updates.tsv"
    )
  )

  stop(
    "Duplicate Phase 4H update keys detected"
  )
}

updates$object_key <- paste(
  updates$dataset,
  updates$object_name,
  sep = "\r"
)

updates$manifest_row <- match(
  updates$object_key,
  manifest$key
)

if (
  anyNA(
    updates$manifest_row
  )
) {
  stop(
    "Phase 4H update contains object absent from manifest"
  )
}

enhanced <- manifest

new_fields <- setdiff(
  unique(updates$field),
  names(enhanced)
)

for (field in new_fields) {
  enhanced[[field]] <- ""
}

updates$current_value <- ""

for (i in seq_len(nrow(updates))) {
  field <- updates$field[i]
  row_index <- updates$manifest_row[i]

  updates$current_value[i] <- clean_text(
    enhanced[[field]][row_index]
  )
}

updates$conflict <- (
  !is_missing(updates$current_value) &
    updates$current_value != updates$value
)

updates$update_decision <- ifelse(
  updates$conflict,
  "CONFLICT_NOT_APPLIED",
  ifelse(
    updates$current_value == updates$value,
    "ALREADY_IDENTICAL",
    "APPLY_TO_ENHANCED_COPY"
  )
)

write_tsv(
  updates,
  file.path(
    table_dir,
    "phase4h_reviewed_updates.tsv"
  )
)

conflicts <- updates[
  updates$conflict,
  ,
  drop = FALSE
]

write_tsv(
  conflicts,
  file.path(
    qc_dir,
    "phase4h_conflicts.tsv"
  )
)

if (nrow(conflicts)) {
  writeLines(
    c(
      paste0("failed_at=", stamp()),
      "core_qc=FAIL",
      paste0(
        "conflicting_updates=",
        nrow(conflicts)
      ),
      "enhanced_manifest_written=0"
    ),
    file.path(
      out_dir,
      "PHASE4H_INCOMPLETE"
    )
  )

  print(conflicts)

  quit(
    save = "no",
    status = 2L
  )
}

for (i in seq_len(nrow(updates))) {
  field <- updates$field[i]
  row_index <- updates$manifest_row[i]

  enhanced[[field]][row_index] <- updates$value[i]
}

enhanced$key <- NULL

write_tsv(
  enhanced,
  file.path(
    table_dir,
    "object_clinical_manifest_phase4h_enhanced.tsv"
  )
)

override_columns <- c(
  "dataset",
  "object_name",
  "field",
  "value",
  "metadata_source",
  "source_reference",
  "notes"
)

write_tsv(
  updates[
    ,
    override_columns,
    drop = FALSE
  ],
  file.path(
    table_dir,
    paste0(
      "phase4e_manual_clinical_overrides_",
      "PHASE4H_REVIEWED.tsv"
    )
  )
)

patient_subset <- enhanced[
  enhanced$dataset %in% c(
    "GSE207493",
    "GSE294109"
  ) &
    !is_missing(enhanced$patient_id),
  ,
  drop = FALSE
]

patient_key <- paste(
  patient_subset$dataset,
  patient_subset$patient_id,
  sep = "\r"
)

patient_groups <- split(
  seq_len(nrow(patient_subset)),
  patient_key
)

patient_rows <- lapply(
  patient_groups,
  function(idx) {
    z <- patient_subset[
      idx,
      ,
      drop = FALSE
    ]

    data.frame(
      dataset = z$dataset[1],
      patient_id = z$patient_id[1],
      n_objects = nrow(z),
      n_cells = sum(
        suppressWarnings(
          as.numeric(z$cells)
        ),
        na.rm = TRUE
      ),
      treatment = collapse_unique(
        z$treatment
      ),
      treatment_status = collapse_unique(
        z$treatment_status
      ),
      response_raw = collapse_unique(
        z$response_raw
      ),
      response_group = collapse_unique(
        z$response_group
      ),
      object_names = paste(
        sort(
          unique(
            clean_text(
              z$object_name
            )
          )
        ),
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }
)

patient_table <- do.call(
  rbind,
  patient_rows
)

patient_table <- patient_table[
  order(
    patient_table$dataset,
    patient_table$patient_id
  ),
  ,
  drop = FALSE
]

patient_conflicts <- patient_table[
  grepl(
    " \\| ",
    patient_table$treatment,
    fixed = TRUE
  ) |
    grepl(
      " \\| ",
      patient_table$treatment_status,
      fixed = TRUE
    ) |
    grepl(
      " \\| ",
      patient_table$response_raw,
      fixed = TRUE
    ) |
    grepl(
      " \\| ",
      patient_table$response_group,
      fixed = TRUE
    ),
  ,
  drop = FALSE
]

write_tsv(
  patient_conflicts,
  file.path(
    qc_dir,
    "patient_level_conflicts.tsv"
  )
)

write_tsv(
  patient_table,
  file.path(
    table_dir,
    "patient_level_clinical_phase4h.tsv"
  )
)

gse294109_patient <- patient_table[
  patient_table$dataset == "GSE294109",
  ,
  drop = FALSE
]

residual_treatment <- treatment_audit[
  treatment_audit$audit_decision ==
    "REVIEW_REQUIRED" &
    treatment_audit$dataset !=
      "GSE294109",
  ,
  drop = FALSE
]

write_tsv(
  residual_treatment,
  file.path(
    table_dir,
    "residual_treatment_review.tsv"
  )
)

treatment_checks <- data.frame(
  treatment = c(
    "Untreated",
    "Ipilimumab + Nivolumab",
    "Cabozantinib",
    "Axitinib + Pembrolizumab",
    "Lenvatinib + Pembrolizumab"
  ),
  expected_objects = c(
    6L,
    12L,
    1L,
    4L,
    2L
  ),
  observed_objects = c(
    sum(
      gse294109$treatment == "Untreated"
    ),
    sum(
      gse294109$treatment ==
        "Ipilimumab + Nivolumab"
    ),
    sum(
      gse294109$treatment ==
        "Cabozantinib"
    ),
    sum(
      gse294109$treatment ==
        "Axitinib + Pembrolizumab"
    ),
    sum(
      gse294109$treatment ==
        "Lenvatinib + Pembrolizumab"
    )
  ),
  stringsAsFactors = FALSE
)

treatment_checks$status <- ifelse(
  treatment_checks$expected_objects ==
    treatment_checks$observed_objects,
  "PASS",
  "FAIL"
)

write_tsv(
  treatment_checks,
  file.path(
    qc_dir,
    "gse294109_treatment_counts.tsv"
  )
)

qc <- data.frame(
  check = c(
    "Input object count",
    "Final GEO matched objects",
    "Conflicting reviewed updates",
    "Patient-level conflicts",
    "GSE207493 patient objects",
    "GSE207493 unique patients",
    "GSE294109 object count",
    "GSE294109 unique patients",
    "GSE294109 treated objects",
    "GSE294109 untreated objects",
    "GSE294109 PR objects",
    "GSE294109 CR objects",
    "GSE294109 No_effect objects",
    "GSE294109 responder patients",
    "GSE294109 non-responder patients",
    "GSE294109 not-applicable patients",
    "Confirmed patients in patient table",
    "Residual treatment review rows",
    "Automatic Phase 4E modification"
  ),
  expected = c(
    "95",
    "95",
    "0",
    "0",
    "19",
    "19",
    "25",
    "15",
    "19",
    "6",
    "9",
    "5",
    "5",
    "10",
    "3",
    "2",
    "34",
    "5",
    "0"
  ),
  observed = c(
    as.character(nrow(manifest)),
    as.character(
      sum(
        resolved$resolved_match_status ==
          "MATCHED"
      )
    ),
    as.character(nrow(conflicts)),
    as.character(nrow(patient_conflicts)),
    as.character(nrow(gse207493_patient)),
    as.character(
      length(
        unique(
          gse207493_patient$patient_id
        )
      )
    ),
    as.character(nrow(gse294109)),
    as.character(
      length(
        unique(
          gse294109$patient_id
        )
      )
    ),
    as.character(
      sum(
        gse294109$treatment_status ==
          "Treated"
      )
    ),
    as.character(
      sum(
        gse294109$treatment_status ==
          "Untreated"
      )
    ),
    as.character(
      sum(
        gse294109$response_raw == "PR"
      )
    ),
    as.character(
      sum(
        gse294109$response_raw == "CR"
      )
    ),
    as.character(
      sum(
        gse294109$response_raw ==
          "No_effect"
      )
    ),
    as.character(
      sum(
        gse294109_patient$response_group ==
          "Responder"
      )
    ),
    as.character(
      sum(
        gse294109_patient$response_group ==
          "Non-responder"
      )
    ),
    as.character(
      sum(
        gse294109_patient$response_group ==
          "Not applicable"
      )
    ),
    as.character(nrow(patient_table)),
    as.character(nrow(residual_treatment)),
    "0"
  ),
  stringsAsFactors = FALSE
)

qc$status <- ifelse(
  qc$expected == qc$observed,
  "PASS",
  "FAIL"
)

qc$notes <- c(
  "",
  "",
  "",
  "Treatment and response must agree across objects from one patient",
  "",
  "",
  "",
  "Fifteen-patient human cohort",
  "",
  "",
  "",
  "",
  "Raw GEO value retained; not relabelled as PD",
  "CR or PR",
  "No_effect",
  "Untreated patients",
  "GSE207493 plus GSE294109",
  "Kept isolated and not written into enhanced manifest",
  "Only independent Phase 4H files were written"
)

write_tsv(
  qc,
  file.path(
    qc_dir,
    "phase4h_qc.tsv"
  )
)

core_fail <- (
  any(qc$status == "FAIL") ||
    any(treatment_checks$status == "FAIL")
)

report <- c(
  "# Phase 4H clinical unit finalization report",
  "",
  paste0("Generated: ", stamp()),
  "",
  "## Object and patient units",
  "",
  paste0(
    "- Input objects: ",
    nrow(manifest)
  ),
  paste0(
    "- Confirmed GSE207493 patients: ",
    length(
      unique(
        gse207493_patient$patient_id
      )
    )
  ),
  paste0(
    "- Confirmed GSE294109 patients: ",
    length(
      unique(
        gse294109$patient_id
      )
    )
  ),
  paste0(
    "- Patient-level table rows: ",
    nrow(patient_table)
  ),
  "",
  "## GSE294109 response",
  "",
  paste0(
    "- Treated patients: ",
    sum(
      gse294109_patient$treatment_status ==
        "Treated"
    )
  ),
  paste0(
    "- Responder patients: ",
    sum(
      gse294109_patient$response_group ==
        "Responder"
    )
  ),
  paste0(
    "- Non-responder patients: ",
    sum(
      gse294109_patient$response_group ==
        "Non-responder"
    )
  ),
  paste0(
    "- Untreated / response not applicable patients: ",
    sum(
      gse294109_patient$response_group ==
        "Not applicable"
    )
  ),
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "Multiple regional objects from the same patient are ",
    "not independent biological replicates. GSE294109 ",
    "response comparisons must use 13 treated patients, ",
    "not 19 treated objects or 25 total objects."
  ),
  "",
  paste0(
    "No_effect is retained as the source response and ",
    "normalized only to Non-responder. It is not assumed ",
    "to equal radiographic or pathological PD."
  ),
  "",
  paste0(
    "Treatment exposure is not a longitudinal timepoint. ",
    "No pre/post relationship is inferred."
  ),
  "",
  "## Safety",
  "",
  paste0(
    "Phase 4E and all RDS files were not modified. ",
    "An independent enhanced manifest and reviewed ",
    "override table were created."
  ),
  "",
  "## Status",
  "",
  if (core_fail) {
    "Phase 4H QC failed."
  } else {
    "Phase 4H QC passed."
  }
)

writeLines(
  report,
  file.path(
    out_dir,
    "PHASE4H_REPORT.md"
  ),
  useBytes = TRUE
)

status_file <- file.path(
  out_dir,
  if (core_fail) {
    "PHASE4H_INCOMPLETE"
  } else {
    "PHASE4H_COMPLETE"
  }
)

writeLines(
  c(
    paste0("completed_at=", stamp()),
    paste0(
      "core_qc=",
      if (core_fail) "FAIL" else "PASS"
    ),
    paste0(
      "object_rows=",
      nrow(manifest)
    ),
    paste0(
      "reviewed_update_rows=",
      nrow(updates)
    ),
    paste0(
      "confirmed_patient_rows=",
      nrow(patient_table)
    ),
    paste0(
      "gse294109_unique_patients=",
      length(
        unique(
          gse294109$patient_id
        )
      )
    ),
    paste0(
      "gse294109_treated_patients=",
      sum(
        gse294109_patient$treatment_status ==
          "Treated"
      )
    ),
    paste0(
      "gse294109_responder_patients=",
      sum(
        gse294109_patient$response_group ==
          "Responder"
      )
    ),
    paste0(
      "gse294109_nonresponder_patients=",
      sum(
        gse294109_patient$response_group ==
          "Non-responder"
      )
    ),
    paste0(
      "residual_treatment_review_rows=",
      nrow(residual_treatment)
    ),
    "phase4e_modified=0",
    "rds_modified=0"
  ),
  status_file
)

cat(
  "STATUS=",
  if (core_fail) "FAIL" else "PASS",
  "\n",
  sep = ""
)

if (core_fail) {
  quit(
    save = "no",
    status = 2L
  )
}
