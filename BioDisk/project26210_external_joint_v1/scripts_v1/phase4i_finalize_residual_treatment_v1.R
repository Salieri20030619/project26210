options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

root <- Sys.getenv(
  "PHASE4I_ROOT",
  unset = "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
)

phase4g_dir <- file.path(
  root,
  "phase4g_resolve_modality_and_response_v1"
)

phase4h_dir <- file.path(
  root,
  "phase4h_finalize_clinical_units_v1"
)

out_dir <- file.path(
  root,
  "phase4i_finalize_residual_treatment_v1"
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
  phase4h_dir,
  "tables",
  "object_clinical_manifest_phase4h_enhanced.tsv"
)

patient_file <- file.path(
  phase4h_dir,
  "tables",
  "patient_level_clinical_phase4h.tsv"
)

residual_file <- file.path(
  phase4h_dir,
  "tables",
  "residual_treatment_review.tsv"
)

resolved_file <- file.path(
  phase4g_dir,
  "tables",
  "object_geo_resolved_mapping.tsv"
)

manifest <- read_tsv(manifest_file)
patient_h <- read_tsv(patient_file)
residual <- read_tsv(residual_file)
resolved <- read_tsv(resolved_file)

require_columns(
  manifest,
  c(
    "dataset",
    "object_name",
    "cells"
  ),
  "Phase 4H enhanced manifest"
)

require_columns(
  patient_h,
  c(
    "dataset",
    "patient_id",
    "n_objects",
    "n_cells",
    "treatment",
    "treatment_status",
    "response_raw",
    "response_group",
    "object_names"
  ),
  "Phase 4H patient table"
)

require_columns(
  residual,
  c(
    "dataset",
    "resolved_gsm",
    "object_name",
    "raw_key",
    "raw_value",
    "source_reference",
    "audit_decision"
  ),
  "Phase 4H residual treatment review"
)

require_columns(
  resolved,
  c(
    "dataset",
    "object_name",
    "resolved_gsm",
    "resolved_geo_title",
    "resolved_match_status"
  ),
  "Phase 4G resolved mapping"
)

manifest_key <- paste(
  manifest$dataset,
  manifest$object_name,
  sep = "\r"
)

resolved_key <- paste(
  resolved$dataset,
  resolved$object_name,
  sep = "\r"
)

if (nrow(manifest) != 95L) {
  stop(
    "Expected 95 Phase 4H objects, observed ",
    nrow(manifest)
  )
}

if (nrow(resolved) != 95L) {
  stop(
    "Expected 95 resolved mappings, observed ",
    nrow(resolved)
  )
}

if (nrow(residual) != 5L) {
  stop(
    "Expected 5 residual treatment rows, observed ",
    nrow(residual)
  )
}

if (anyDuplicated(manifest_key)) {
  stop("Duplicated object keys in Phase 4H manifest")
}

if (anyDuplicated(resolved_key)) {
  stop("Duplicated object keys in resolved mapping")
}

if (!setequal(manifest_key, resolved_key)) {
  stop(
    "Phase 4H manifest and Phase 4G mapping objects differ"
  )
}

if (
  sum(
    resolved$resolved_match_status == "MATCHED"
  ) != 95L
) {
  stop("Not all resolved mappings have MATCHED status")
}

gse202374_map <- data.frame(
  resolved_gsm = c(
    "GSM6112169",
    "GSM6112170",
    "GSM6112171",
    "GSM6112172",
    "GSM6112173",
    "GSM6112174"
  ),
  patient_id = c(
    "C1",
    "C23",
    "P6",
    "P10",
    "P14",
    "P15"
  ),
  treatment = c(
    "None",
    "None",
    "15 Gy radiation dose",
    "15 Gy radiation dose",
    "15 Gy radiation dose",
    "15 Gy radiation dose"
  ),
  treatment_status = c(
    "Untreated",
    "Untreated",
    "Treated",
    "Treated",
    "Treated",
    "Treated"
  ),
  treatment_group = c(
    "No radiotherapy",
    "No radiotherapy",
    "Radiotherapy",
    "Radiotherapy",
    "Radiotherapy",
    "Radiotherapy"
  ),
  stringsAsFactors = FALSE
)

mapping_index <- match(
  gse202374_map$resolved_gsm,
  resolved$resolved_gsm
)

if (anyNA(mapping_index)) {
  stop("One or more expected GSE202374 GSMs are missing")
}

gse202374_objects <- resolved[
  mapping_index,
  ,
  drop = FALSE
]

if (
  !all(
    gse202374_objects$dataset == "GSE202374"
  )
) {
  stop("A GSE202374 GSM mapped to the wrong dataset")
}

gse202374_map$object_name <-
  gse202374_objects$object_name

gse202374_map$resolved_geo_title <-
  gse202374_objects$resolved_geo_title

gse202374_map$source_reference <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
  gse202374_map$resolved_gsm
)

expected_titles <- c(
  "Sample 1 (C1 - RS00123497)",
  "Sample 2 (C23 - RS01519713)",
  "Sample 3 (P6 - RS01813260)",
  "Sample 4 (P10 - RS01985670)",
  "Sample 5 (P14 - RS02065737)",
  "Sample 6 (P15 - RS02098409)"
)

if (
  !identical(
    clean_text(
      gse202374_map$resolved_geo_title
    ),
    expected_titles
  )
) {
  comparison <- data.frame(
    resolved_gsm = gse202374_map$resolved_gsm,
    expected_title = expected_titles,
    observed_title =
      gse202374_map$resolved_geo_title,
    stringsAsFactors = FALSE
  )

  write_tsv(
    comparison,
    file.path(
      qc_dir,
      "gse202374_title_mismatch.tsv"
    )
  )

  stop("Unexpected GSE202374 GEO titles")
}

gse178481_index <- which(
  resolved$dataset == "GSE178481" &
    resolved$resolved_gsm == "GSM5392403"
)

if (length(gse178481_index) != 1L) {
  stop(
    "Expected exactly one GSE178481 GSM5392403 object"
  )
}

gse178481_object <- resolved[
  gse178481_index,
  ,
  drop = FALSE
]

updates <- rbind(
  make_updates(
    dataset = "GSE202374",
    object_name =
      gse202374_map$object_name,
    field = "patient_id",
    value =
      gse202374_map$patient_id,
    metadata_source =
      "NCBI GEO GSE202374 Sample_title",
    source_reference =
      gse202374_map$source_reference,
    notes = paste0(
      "Study-specific de-identified patient ID ",
      "extracted from the exact matched GEO title"
    )
  ),
  make_updates(
    dataset = "GSE202374",
    object_name =
      gse202374_map$object_name,
    field = "treatment",
    value =
      gse202374_map$treatment,
    metadata_source = paste0(
      "NCBI GEO GSE202374 ",
      "Sample_characteristics_ch1"
    ),
    source_reference =
      gse202374_map$source_reference,
    notes = paste0(
      "Exact source treatment retained; ",
      "irradiated tumors were resected ",
      "four weeks after 15 Gy radiotherapy"
    )
  ),
  make_updates(
    dataset = "GSE202374",
    object_name =
      gse202374_map$object_name,
    field = "treatment_status",
    value =
      gse202374_map$treatment_status,
    metadata_source = paste0(
      "Derived from exact GSE202374 ",
      "treatment annotation"
    ),
    source_reference =
      gse202374_map$source_reference,
    notes =
      "Treated denotes preoperative 15 Gy radiotherapy"
  ),
  make_updates(
    dataset = "GSE202374",
    object_name =
      gse202374_map$object_name,
    field = "treatment_group",
    value =
      gse202374_map$treatment_group,
    metadata_source = paste0(
      "Controlled normalization of GSE202374 ",
      "treatment annotation"
    ),
    source_reference =
      gse202374_map$source_reference,
    notes = paste0(
      "Radiotherapy versus nephrectomy-only ",
      "comparison group"
    )
  ),
  make_updates(
    dataset = "GSE178481",
    object_name =
      gse178481_object$object_name,
    field = "source_treatment_annotation_raw",
    value =
      "chemotherapy for collicular lymphoma",
    metadata_source = paste0(
      "NCBI GEO GSM5392403 ",
      "Sample_characteristics_ch1"
    ),
    source_reference = paste0(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
      gse178481_object$resolved_gsm
    ),
    notes = paste0(
      "Exact source wording retained, including ",
      "the source spelling; not classified as ",
      "RCC-directed treatment"
    )
  ),
  make_updates(
    dataset = "GSE178481",
    object_name =
      gse178481_object$object_name,
    field = "treatment_context_note",
    value = paste0(
      "Excluded from RCC treatment harmonization: ",
      "source annotation describes chemotherapy ",
      "for another malignancy"
    ),
    metadata_source = paste0(
      "Phase 4I adjudication using GSE178481 ",
      "series design and GSM5392403 annotation"
    ),
    source_reference = paste0(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
      gse178481_object$resolved_gsm
    ),
    notes = paste0(
      "GSE178481 defines the ccRCC cohort as ",
      "treatment-naive; no RCC treatment is inferred"
    )
  )
)

update_key <- paste(
  updates$dataset,
  updates$object_name,
  updates$field,
  sep = "\r"
)

if (anyDuplicated(update_key)) {
  stop("Duplicated Phase 4I update keys")
}

updates$manifest_row <- match(
  paste(
    updates$dataset,
    updates$object_name,
    sep = "\r"
  ),
  manifest_key
)

if (anyNA(updates$manifest_row)) {
  stop("Phase 4I update refers to an unknown object")
}

new_fields <- setdiff(
  unique(updates$field),
  names(manifest)
)

for (field in new_fields) {
  manifest[[field]] <- ""
}

updates$current_value <- ""

for (i in seq_len(nrow(updates))) {
  updates$current_value[i] <- clean_text(
    manifest[
      updates$manifest_row[i],
      updates$field[i]
    ]
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
    "APPLY_TO_FINAL_COPY"
  )
)

write_tsv(
  updates,
  file.path(
    table_dir,
    "phase4i_reviewed_updates.tsv"
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
    "phase4i_conflicts.tsv"
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
      "phase4h_modified=0",
      "rds_modified=0"
    ),
    file.path(
      out_dir,
      "PHASE4I_INCOMPLETE"
    )
  )

  print(conflicts)

  quit(
    save = "no",
    status = 2L
  )
}

for (i in seq_len(nrow(updates))) {
  manifest[
    updates$manifest_row[i],
    updates$field[i]
  ] <- updates$value[i]
}

write_tsv(
  manifest,
  file.path(
    table_dir,
    "object_clinical_manifest_phase4i_final.tsv"
  )
)

adjudication <- residual

adjudication$final_decision <- ""
adjudication$final_field <- ""
adjudication$final_value <- ""
adjudication$rationale <- ""

accept_index <- (
  adjudication$dataset == "GSE202374" &
    adjudication$resolved_gsm %in%
      c(
        "GSM6112171",
        "GSM6112172",
        "GSM6112173",
        "GSM6112174"
      )
)

reject_index <- (
  adjudication$dataset == "GSE178481" &
    adjudication$resolved_gsm == "GSM5392403"
)

adjudication$final_decision[
  accept_index
] <- "ACCEPT_AS_RCC_RADIOTHERAPY"

adjudication$final_field[
  accept_index
] <- "treatment"

adjudication$final_value[
  accept_index
] <- "15 Gy radiation dose"

adjudication$rationale[
  accept_index
] <- paste0(
  "Official GEO treatment annotation and protocol ",
  "confirm preoperative 15 Gy radiotherapy with ",
  "tumor resection four weeks later"
)

adjudication$final_decision[
  reject_index
] <- "EXCLUDE_FROM_RCC_TREATMENT_RETAIN_CONTEXT"

adjudication$final_field[
  reject_index
] <- "source_treatment_annotation_raw"

adjudication$final_value[
  reject_index
] <- adjudication$raw_value[
  reject_index
]

adjudication$rationale[
  reject_index
] <- paste0(
  "The annotation refers to chemotherapy for another ",
  "malignancy, while GSE178481 defines the ccRCC ",
  "cohort as treatment-naive"
)

if (
  any(
    adjudication$final_decision == ""
  )
) {
  print(
    adjudication[
      adjudication$final_decision == "",
      ,
      drop = FALSE
    ]
  )

  stop("One or more residual rows remain unresolved")
}

write_tsv(
  adjudication,
  file.path(
    table_dir,
    "residual_treatment_final_adjudication.tsv"
  )
)

gse202_final_index <- match(
  gse202374_map$object_name,
  manifest$object_name
)

if (anyNA(gse202_final_index)) {
  stop("Final manifest lacks a GSE202374 object")
}

gse202_final <- manifest[
  gse202_final_index,
  ,
  drop = FALSE
]

gse202_patient <- data.frame(
  dataset = "GSE202374",
  patient_id = gse202_final$patient_id,
  n_objects = 1L,
  n_cells = suppressWarnings(
    as.numeric(
      gse202_final$cells
    )
  ),
  treatment = gse202_final$treatment,
  treatment_status =
    gse202_final$treatment_status,
  response_raw = "",
  response_group = "",
  object_names = gse202_final$object_name,
  stringsAsFactors = FALSE
)

patient_final <- rbind(
  patient_h[
    ,
    names(gse202_patient),
    drop = FALSE
  ],
  gse202_patient
)

patient_final <- patient_final[
  order(
    patient_final$dataset,
    patient_final$patient_id
  ),
  ,
  drop = FALSE
]

patient_key <- paste(
  patient_final$dataset,
  patient_final$patient_id,
  sep = "\r"
)

if (anyDuplicated(patient_key)) {
  duplicate_patients <- patient_final[
    duplicated(patient_key) |
      duplicated(
        patient_key,
        fromLast = TRUE
      ),
    ,
    drop = FALSE
  ]

  write_tsv(
    duplicate_patients,
    file.path(
      qc_dir,
      "duplicate_patient_rows.tsv"
    )
  )

  stop("Duplicated patient rows after Phase 4I")
}

write_tsv(
  patient_final,
  file.path(
    table_dir,
    "patient_level_clinical_phase4i.tsv"
  )
)

qc <- data.frame(
  check = c(
    "Phase 4H input objects",
    "Phase 4H patient rows",
    "Residual treatment rows",
    "Accepted residual rows",
    "Excluded residual rows",
    "Remaining unresolved residual rows",
    "GSE202374 objects",
    "GSE202374 unique patient IDs",
    "GSE202374 radiotherapy patients",
    "GSE202374 no-radiotherapy patients",
    "Phase 4I reviewed update rows",
    "Conflicting Phase 4I updates",
    "Final manifest objects",
    "Final confirmed patient rows",
    "Automatic Phase 4H modification",
    "RDS modification"
  ),
  expected = c(
    "95",
    "34",
    "5",
    "4",
    "1",
    "0",
    "6",
    "6",
    "4",
    "2",
    "26",
    "0",
    "95",
    "40",
    "0",
    "0"
  ),
  observed = c(
    as.character(nrow(manifest)),
    as.character(nrow(patient_h)),
    as.character(nrow(residual)),
    as.character(sum(accept_index)),
    as.character(sum(reject_index)),
    as.character(
      sum(
        adjudication$final_decision == ""
      )
    ),
    as.character(nrow(gse202_final)),
    as.character(
      length(
        unique(
          gse202_final$patient_id
        )
      )
    ),
    as.character(
      sum(
        gse202_final$treatment_group ==
          "Radiotherapy"
      )
    ),
    as.character(
      sum(
        gse202_final$treatment_group ==
          "No radiotherapy"
      )
    ),
    as.character(nrow(updates)),
    as.character(nrow(conflicts)),
    as.character(nrow(manifest)),
    as.character(nrow(patient_final)),
    "0",
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
  "GSE207493 plus GSE294109",
  "",
  "Four GSE202374 irradiated RCC tumors",
  paste0(
    "GSM5392403 excluded from RCC treatment; ",
    "raw context retained"
  ),
  "",
  "",
  "C1, C23, P6, P10, P14 and P15",
  "Preoperative 15 Gy radiotherapy",
  "Nephrectomy-only controls",
  "",
  "",
  "",
  paste0(
    "GSE207493, GSE294109 and GSE202374 ",
    "confirmed patient tables"
  ),
  "Only independent Phase 4I files were written",
  "No RDS file was opened or modified"
)

write_tsv(
  qc,
  file.path(
    qc_dir,
    "phase4i_qc.tsv"
  )
)

core_fail <- any(qc$status == "FAIL")

report <- c(
  "# Phase 4I residual treatment finalization report",
  "",
  paste0("Generated: ", stamp()),
  "",
  "## Residual adjudication",
  "",
  paste0(
    "- Residual rows reviewed: ",
    nrow(residual)
  ),
  paste0(
    "- Accepted as RCC radiotherapy: ",
    sum(accept_index)
  ),
  paste0(
    "- Excluded from RCC treatment: ",
    sum(reject_index)
  ),
  "- Remaining unresolved rows: 0",
  "",
  "## GSE202374",
  "",
  "- Confirmed patients: 6",
  "- Preoperative 15 Gy radiotherapy: 4",
  "- Nephrectomy-only controls: 2",
  paste0(
    "- Irradiated tumors were resected ",
    "four weeks after radiotherapy"
  ),
  "",
  "## GSE178481",
  "",
  paste0(
    "- GSM5392403 source wording was retained in ",
    "source_treatment_annotation_raw."
  ),
  paste0(
    "- It was not entered into the harmonized RCC ",
    "treatment field because it refers to another ",
    "malignancy and the ccRCC cohort is defined as ",
    "treatment-naive."
  ),
  "",
  "## Final units",
  "",
  paste0(
    "- Final object manifest rows: ",
    nrow(manifest)
  ),
  paste0(
    "- Confirmed patient table rows: ",
    nrow(patient_final)
  ),
  "",
  "## Safety",
  "",
  paste0(
    "Phase 4H files and all RDS objects were not modified. ",
    "Only independent Phase 4I outputs were created."
  ),
  "",
  "## Status",
  "",
  if (core_fail) {
    "Phase 4I QC failed."
  } else {
    "Phase 4I QC passed."
  }
)

writeLines(
  report,
  file.path(
    out_dir,
    "PHASE4I_REPORT.md"
  ),
  useBytes = TRUE
)

status_file <- file.path(
  out_dir,
  if (core_fail) {
    "PHASE4I_INCOMPLETE"
  } else {
    "PHASE4I_COMPLETE"
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
      "confirmed_patient_rows=",
      nrow(patient_final)
    ),
    paste0(
      "residual_rows_reviewed=",
      nrow(residual)
    ),
    paste0(
      "residual_rows_accepted=",
      sum(accept_index)
    ),
    paste0(
      "residual_rows_excluded=",
      sum(reject_index)
    ),
    "residual_rows_unresolved=0",
    "gse202374_patients=6",
    "gse202374_radiotherapy_patients=4",
    "gse202374_control_patients=2",
    "phase4h_modified=0",
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
