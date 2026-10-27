options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

root <- Sys.getenv(
  "PHASE4G_ROOT",
  unset = "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
)

phase4e_dir <- file.path(
  root,
  "phase4e_clinical_metadata_harmonization_v1"
)

phase4f_dir <- file.path(
  root,
  "phase4f_geo_clinical_evidence_v1"
)

out_dir <- file.path(
  root,
  "phase4g_resolve_modality_and_response_v1"
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
  y <- toupper(clean_text(x))

  y %in% c(
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
  missing_columns <- setdiff(
    columns,
    names(x)
  )

  if (length(missing_columns)) {
    stop(
      label,
      " missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

manifest_file <- file.path(
  phase4e_dir,
  "tables",
  "object_clinical_manifest.tsv"
)

mapping_file <- file.path(
  phase4f_dir,
  "tables",
  "object_geo_exact_match.tsv"
)

sample_metadata_file <- file.path(
  phase4f_dir,
  "tables",
  "geo_sample_metadata_long.tsv"
)

evidence_file <- file.path(
  phase4f_dir,
  "tables",
  "gsm_clinical_evidence.tsv"
)

manifest <- read_tsv(manifest_file)
mapping <- read_tsv(mapping_file)
sample_metadata <- read_tsv(sample_metadata_file)
gsm_evidence <- read_tsv(evidence_file)

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
  mapping,
  c(
    "dataset",
    "object_name",
    "match_status",
    "match_method",
    "gsm",
    "candidate_gsms"
  ),
  "Phase 4F object mapping"
)

require_columns(
  sample_metadata,
  c(
    "dataset",
    "gsm",
    "field",
    "value",
    "source_reference"
  ),
  "GEO sample metadata"
)

require_columns(
  gsm_evidence,
  c(
    "dataset",
    "gsm",
    "target_field",
    "raw_key",
    "raw_value",
    "standardized_value",
    "metadata_source",
    "source_reference"
  ),
  "GEO clinical evidence"
)

if (nrow(manifest) != 95L) {
  stop(
    "Expected 95 manifest objects, observed ",
    nrow(manifest)
  )
}

if (nrow(mapping) != 95L) {
  stop(
    "Expected 95 object mappings, observed ",
    nrow(mapping)
  )
}

sample_titles <- sample_metadata[
  sample_metadata$field == "Sample_title",
  c(
    "dataset",
    "gsm",
    "value",
    "source_reference"
  ),
  drop = FALSE
]

names(sample_titles)[
  names(sample_titles) == "value"
] <- "geo_title"

names(sample_titles)[
  names(sample_titles) == "source_reference"
] <- "title_source_reference"

sample_titles <- unique(sample_titles)

resolved_rows <- vector(
  "list",
  nrow(mapping)
)

modality_evidence_rows <- list()

for (i in seq_len(nrow(mapping))) {
  z <- mapping[i, , drop = FALSE]

  final_status <- clean_text(
    z$match_status
  )

  final_method <- clean_text(
    z$match_method
  )

  final_gsm <- clean_text(
    z$gsm
  )

  final_title <- ""

  if (
    z$dataset == "GSE207493" &&
    z$match_status == "AMBIGUOUS"
  ) {
    candidate_gsms <- unlist(
      strsplit(
        clean_text(z$candidate_gsms),
        ";",
        fixed = TRUE
      ),
      use.names = FALSE
    )

    candidate_gsms <- unique(
      clean_text(candidate_gsms)
    )

    candidate_gsms <- candidate_gsms[
      candidate_gsms != ""
    ]

    candidate_titles <- sample_titles[
      sample_titles$dataset == z$dataset &
        sample_titles$gsm %in% candidate_gsms,
      ,
      drop = FALSE
    ]

    candidate_titles$is_scrna <- grepl(
      "^scRNA[-_ ]?seq",
      candidate_titles$geo_title,
      ignore.case = TRUE,
      perl = TRUE
    )

    candidate_titles$is_scatac <- grepl(
      "^scATAC[-_ ]?seq",
      candidate_titles$geo_title,
      ignore.case = TRUE,
      perl = TRUE
    )

    modality_evidence_rows[[
      length(modality_evidence_rows) + 1L
    ]] <- data.frame(
      dataset = z$dataset,
      object_name = z$object_name,
      candidate_gsm = candidate_titles$gsm,
      geo_title = candidate_titles$geo_title,
      is_scrna = candidate_titles$is_scrna,
      is_scatac = candidate_titles$is_scatac,
      selection = ifelse(
        candidate_titles$is_scrna,
        "SELECT",
        "EXCLUDE_NON_EXPRESSION_MODALITY"
      ),
      source_reference =
        candidate_titles$title_source_reference,
      stringsAsFactors = FALSE
    )

    selected <- candidate_titles[
      candidate_titles$is_scrna,
      ,
      drop = FALSE
    ]

    if (nrow(selected) != 1L) {
      stop(
        "Expected exactly one scRNA-seq GSM for ",
        z$object_name,
        "; observed ",
        nrow(selected)
      )
    }

    final_status <- "MATCHED"
    final_method <- "exact_sample_id_plus_scrna_modality"
    final_gsm <- selected$gsm[1]
    final_title <- selected$geo_title[1]
  } else if (final_status == "MATCHED") {
    title_hit <- sample_titles[
      sample_titles$dataset == z$dataset &
        sample_titles$gsm == final_gsm,
      ,
      drop = FALSE
    ]

    final_title <- collapse_unique(
      title_hit$geo_title
    )
  }

  resolved_rows[[i]] <- data.frame(
    dataset = z$dataset,
    object_name = z$object_name,
    source_file = z$source_file,
    phase4f_match_status = z$match_status,
    phase4f_match_method = z$match_method,
    phase4f_gsm = z$gsm,
    candidate_gsms = z$candidate_gsms,
    resolved_match_status = final_status,
    resolved_match_method = final_method,
    resolved_gsm = final_gsm,
    resolved_geo_title = final_title,
    stringsAsFactors = FALSE
  )
}

resolved_mapping <- do.call(
  rbind,
  resolved_rows
)

modality_evidence <- if (
  length(modality_evidence_rows)
) {
  do.call(
    rbind,
    modality_evidence_rows
  )
} else {
  data.frame(
    dataset = character(),
    object_name = character(),
    candidate_gsm = character(),
    geo_title = character(),
    is_scrna = logical(),
    is_scatac = logical(),
    selection = character(),
    source_reference = character(),
    stringsAsFactors = FALSE
  )
}

initial_gse207493_ambiguous <- sum(
  mapping$dataset == "GSE207493" &
    mapping$match_status == "AMBIGUOUS"
)

resolved_gse207493 <- sum(
  resolved_mapping$dataset == "GSE207493" &
    resolved_mapping$phase4f_match_status ==
      "AMBIGUOUS" &
    resolved_mapping$resolved_match_status ==
      "MATCHED"
)

remaining_ambiguous <- sum(
  resolved_mapping$resolved_match_status ==
    "AMBIGUOUS"
)

remaining_unmatched <- sum(
  resolved_mapping$resolved_match_status ==
    "UNMATCHED"
)

write_tsv(
  modality_evidence,
  file.path(
    table_dir,
    "gse207493_modality_disambiguation_evidence.tsv"
  )
)

write_tsv(
  resolved_mapping,
  file.path(
    table_dir,
    "object_geo_resolved_mapping.tsv"
  )
)

gse207493_mapping <- resolved_mapping[
  resolved_mapping$dataset == "GSE207493" &
    resolved_mapping$resolved_match_status ==
      "MATCHED",
  ,
  drop = FALSE
]

gse207493_mapping$patient_id <- sub(
  "^scRNA[-_ ]?seq[_ -]*",
  "",
  gse207493_mapping$resolved_geo_title,
  ignore.case = TRUE,
  perl = TRUE
)

gse207493_mapping$patient_id <- clean_text(
  gse207493_mapping$patient_id
)

valid_patient_id <- grepl(
  "^RCC[0-9]+$",
  gse207493_mapping$patient_id,
  ignore.case = TRUE
)

if (!all(valid_patient_id)) {
  bad <- gse207493_mapping[
    !valid_patient_id,
    c(
      "object_name",
      "resolved_gsm",
      "resolved_geo_title",
      "patient_id"
    ),
    drop = FALSE
  ]

  print(bad)

  stop(
    "Invalid GSE207493 patient identifiers detected"
  )
}

gse207493_mapping$patient_id <- toupper(
  gse207493_mapping$patient_id
)

if (
  anyDuplicated(
    gse207493_mapping$patient_id
  )
) {
  stop(
    "GSE207493 patient identifiers are not unique"
  )
}

patient_id_evidence <- data.frame(
  dataset = gse207493_mapping$dataset,
  object_name = gse207493_mapping$object_name,
  gsm = gse207493_mapping$resolved_gsm,
  geo_title =
    gse207493_mapping$resolved_geo_title,
  patient_id =
    gse207493_mapping$patient_id,
  evidence_statement = paste0(
    "GSE207493 overall design states that 19 ccRCC ",
    "samples from 19 patients were each profiled by ",
    "scRNA-seq and scATAC-seq; the RCC identifier in ",
    "the selected scRNA-seq title is therefore retained ",
    "as a study-specific de-identified patient ID."
  ),
  source_reference = paste0(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
    gse207493_mapping$resolved_gsm
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  patient_id_evidence,
  file.path(
    table_dir,
    "gse207493_patient_id_evidence.tsv"
  )
)

resolved_evidence <- merge(
  resolved_mapping[
    resolved_mapping$resolved_match_status ==
      "MATCHED",
    c(
      "dataset",
      "object_name",
      "resolved_gsm",
      "resolved_geo_title",
      "resolved_match_method"
    ),
    drop = FALSE
  ],
  gsm_evidence,
  by.x = c(
    "dataset",
    "resolved_gsm"
  ),
  by.y = c(
    "dataset",
    "gsm"
  ),
  all = FALSE,
  sort = FALSE
)

response_evidence <- resolved_evidence[
  resolved_evidence$target_field ==
    "response_raw",
  ,
  drop = FALSE
]

response_audit_rows <- list()

if (nrow(response_evidence)) {
  response_groups <- split(
    seq_len(nrow(response_evidence)),
    interaction(
      response_evidence$dataset,
      response_evidence$object_name,
      drop = TRUE,
      lex.order = TRUE
    )
  )

  allowed_response <- c(
    "CR",
    "PR",
    "SD",
    "PD",
    "Responder",
    "Non-responder"
  )

  for (idx in response_groups) {
    z <- response_evidence[
      idx,
      ,
      drop = FALSE
    ]

    raw_values <- unique(
      clean_text(z$raw_value)
    )

    raw_values <- raw_values[
      !is_missing(raw_values)
    ]

    standardized_values <- unique(
      clean_text(z$standardized_value)
    )

    standardized_values <- standardized_values[
      !is_missing(standardized_values) &
        standardized_values != "Unknown"
    ]

    manifest_hit <- manifest[
      manifest$dataset == z$dataset[1] &
        manifest$object_name ==
          z$object_name[1],
      ,
      drop = FALSE
    ]

    existing_response <- if (
      "response_raw" %in% names(manifest_hit)
    ) {
      collapse_unique(
        manifest_hit$response_raw
      )
    } else {
      ""
    }

    eligibility <- if (
      length(raw_values) != 1L
    ) {
      "REJECT_MULTIPLE_OR_EMPTY_RAW_VALUES"
    } else if (
      length(standardized_values) != 1L
    ) {
      "REJECT_UNKNOWN_OR_CONFLICTING_STANDARDIZATION"
    } else if (
      !standardized_values %in%
        allowed_response
    ) {
      "REJECT_UNSUPPORTED_RESPONSE_CLASS"
    } else if (
      !is_missing(existing_response) &&
      existing_response != raw_values
    ) {
      "REJECT_CONFLICT_WITH_EXISTING_VALUE"
    } else if (
      !is_missing(existing_response) &&
      existing_response == raw_values
    ) {
      "ALREADY_PRESENT"
    } else {
      "ACCEPT_CANDIDATE"
    }

    response_audit_rows[[
      length(response_audit_rows) + 1L
    ]] <- data.frame(
      dataset = z$dataset[1],
      object_name = z$object_name[1],
      gsm = z$resolved_gsm[1],
      geo_title = z$resolved_geo_title[1],
      match_method = z$resolved_match_method[1],
      raw_key = collapse_unique(z$raw_key),
      raw_value = collapse_unique(raw_values),
      standardized_response =
        collapse_unique(
          standardized_values,
          default = "Unknown"
        ),
      existing_response = existing_response,
      eligibility = eligibility,
      metadata_source =
        collapse_unique(z$metadata_source),
      source_reference =
        collapse_unique(z$source_reference),
      stringsAsFactors = FALSE
    )
  }
}

response_audit <- if (
  length(response_audit_rows)
) {
  do.call(
    rbind,
    response_audit_rows
  )
} else {
  data.frame(
    dataset = character(),
    object_name = character(),
    gsm = character(),
    geo_title = character(),
    match_method = character(),
    raw_key = character(),
    raw_value = character(),
    standardized_response = character(),
    existing_response = character(),
    eligibility = character(),
    metadata_source = character(),
    source_reference = character(),
    stringsAsFactors = FALSE
  )
}

response_audit <- response_audit[
  order(
    response_audit$dataset,
    response_audit$object_name
  ),
  ,
  drop = FALSE
]

write_tsv(
  response_audit,
  file.path(
    table_dir,
    "response_candidate_audit.tsv"
  )
)

treatment_evidence <- resolved_evidence[
  resolved_evidence$target_field ==
    "treatment",
  ,
  drop = FALSE
]

if (nrow(treatment_evidence)) {
  treatment_value <- tolower(
    clean_text(
      treatment_evidence$raw_value
    )
  )

  invalid_exact <- treatment_value %in% c(
    "",
    "no",
    "none",
    "na",
    "n/a",
    "not applicable",
    "unknown",
    "yes"
  )

  specimen_pattern <- grepl(
    paste0(
      "tumou?r|thrombus|biopsy|tissue|",
      "primary rcc|metastatic site|",
      "normal kidney|adjacent"
    ),
    treatment_value,
    perl = TRUE
  )

  therapy_pattern <- grepl(
    paste0(
      "nivo|pembro|ipilimumab|avelumab|",
      "sunitinib|pazopanib|axitinib|",
      "cabozantinib|lenvatinib|",
      "everolimus|temsirolimus|",
      "bevacizumab|tivozanib|sorafenib|",
      "chemotherapy|immunotherapy|",
      "radiotherapy|anti.?pd|anti.?ctla|",
      "checkpoint|tki|treatment|therapy"
    ),
    treatment_value,
    perl = TRUE
  )

  treatment_evidence$audit_decision <- ifelse(
    invalid_exact,
    "REJECT_NON_TREATMENT_VALUE",
    ifelse(
      specimen_pattern & !therapy_pattern,
      "REJECT_SPECIMEN_DESCRIPTION",
      "REVIEW_REQUIRED"
    )
  )
} else {
  treatment_evidence$audit_decision <-
    character()
}

write_tsv(
  treatment_evidence,
  file.path(
    table_dir,
    "treatment_candidate_audit.tsv"
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

patient_overrides <- data.frame(
  dataset = patient_id_evidence$dataset,
  object_name = patient_id_evidence$object_name,
  field = "patient_id",
  value = patient_id_evidence$patient_id,
  metadata_source = paste0(
    "NCBI GEO GSE207493 overall design; ",
    "GEO Sample_title"
  ),
  source_reference =
    patient_id_evidence$source_reference,
  notes = paste0(
    "PHASE4G_ACCEPT_CANDIDATE; study-specific ",
    "de-identified patient ID; not a hospital MRN"
  ),
  stringsAsFactors = FALSE
)

accepted_response <- response_audit[
  response_audit$eligibility ==
    "ACCEPT_CANDIDATE",
  ,
  drop = FALSE
]

response_overrides <- if (
  nrow(accepted_response)
) {
  data.frame(
    dataset = accepted_response$dataset,
    object_name =
      accepted_response$object_name,
    field = "response_raw",
    value = accepted_response$raw_value,
    metadata_source =
      accepted_response$metadata_source,
    source_reference =
      accepted_response$source_reference,
    notes = paste0(
      "PHASE4G_ACCEPT_CANDIDATE; standardized_response=",
      accepted_response$standardized_response,
      "; exact object-to-GSM mapping; ",
      "requires Phase4E rerun for final harmonization"
    ),
    stringsAsFactors = FALSE
  )
} else {
  as.data.frame(
    setNames(
      replicate(
        length(override_columns),
        character(),
        simplify = FALSE
      ),
      override_columns
    ),
    stringsAsFactors = FALSE
  )
}

proposed_overrides <- rbind(
  patient_overrides[
    ,
    override_columns,
    drop = FALSE
  ],
  response_overrides[
    ,
    override_columns,
    drop = FALSE
  ]
)

proposed_overrides <- proposed_overrides[
  order(
    proposed_overrides$dataset,
    proposed_overrides$object_name,
    proposed_overrides$field
  ),
  ,
  drop = FALSE
]

write_tsv(
  proposed_overrides,
  file.path(
    table_dir,
    paste0(
      "phase4e_manual_clinical_overrides_",
      "PHASE4G_ACCEPT_CANDIDATES.tsv"
    )
  )
)

override_summary <- as.data.frame(
  table(
    proposed_overrides$dataset,
    proposed_overrides$field
  ),
  stringsAsFactors = FALSE
)

names(override_summary) <- c(
  "dataset",
  "field",
  "candidate_rows"
)

override_summary <- override_summary[
  override_summary$candidate_rows > 0,
  ,
  drop = FALSE
]

write_tsv(
  override_summary,
  file.path(
    table_dir,
    "phase4g_override_summary.tsv"
  )
)

qc <- data.frame(
  check = c(
    "Phase 4E object count",
    "Phase 4F mapping count",
    "Initial GSE207493 ambiguous objects",
    "GSE207493 ambiguity resolved by modality",
    "Remaining ambiguous objects",
    "Remaining unmatched objects",
    "Final matched objects",
    "GSE207493 unique patient IDs",
    "Automatic Phase 4E modification"
  ),
  expected = c(
    "95",
    "95",
    "19",
    "19",
    "0",
    "0",
    "95",
    "19",
    "0"
  ),
  observed = c(
    as.character(nrow(manifest)),
    as.character(nrow(mapping)),
    as.character(
      initial_gse207493_ambiguous
    ),
    as.character(resolved_gse207493),
    as.character(remaining_ambiguous),
    as.character(remaining_unmatched),
    as.character(
      sum(
        resolved_mapping$resolved_match_status ==
          "MATCHED"
      )
    ),
    as.character(
      length(
        unique(
          patient_id_evidence$patient_id
        )
      )
    ),
    "0"
  ),
  status = c(
    if (nrow(manifest) == 95L) "PASS" else "FAIL",
    if (nrow(mapping) == 95L) "PASS" else "FAIL",
    if (
      initial_gse207493_ambiguous == 19L
    ) "PASS" else "FAIL",
    if (
      resolved_gse207493 == 19L
    ) "PASS" else "FAIL",
    if (
      remaining_ambiguous == 0L
    ) "PASS" else "FAIL",
    if (
      remaining_unmatched == 0L
    ) "PASS" else "FAIL",
    if (
      sum(
        resolved_mapping$resolved_match_status ==
          "MATCHED"
      ) == 95L
    ) "PASS" else "FAIL",
    if (
      length(
        unique(
          patient_id_evidence$patient_id
        )
      ) == 19L
    ) "PASS" else "FAIL",
    "PASS"
  ),
  notes = c(
    "",
    "",
    paste0(
      "Each RCC sample has paired scRNA-seq ",
      "and scATAC-seq GEO records"
    ),
    paste0(
      "Expression RDS objects were mapped to ",
      "the unique scRNA-seq GEO record"
    ),
    "",
    "",
    "",
    paste0(
      "RCC identifiers are study-specific ",
      "de-identified patient IDs"
    ),
    paste0(
      "Only a separate candidate override ",
      "file was generated"
    )
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  qc,
  file.path(
    qc_dir,
    "phase4g_qc.tsv"
  )
)

core_fail <- any(
  qc$status == "FAIL"
)

response_status_counts <- if (
  nrow(response_audit)
) {
  as.data.frame(
    table(response_audit$eligibility),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    Var1 = character(),
    Freq = integer(),
    stringsAsFactors = FALSE
  )
}

names(response_status_counts) <- c(
  "eligibility",
  "objects"
)

treatment_status_counts <- if (
  nrow(treatment_evidence)
) {
  as.data.frame(
    table(
      treatment_evidence$audit_decision
    ),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    audit_decision = character(),
    rows = integer(),
    stringsAsFactors = FALSE
  )
}

if (nrow(treatment_status_counts)) {
  names(treatment_status_counts) <- c(
    "audit_decision",
    "rows"
  )
}

report <- c(
  "# Phase 4G modality and clinical evidence report",
  "",
  paste0("Generated: ", stamp()),
  "",
  "## Deterministic mapping result",
  "",
  paste0(
    "- Initial GSE207493 ambiguous objects: ",
    initial_gse207493_ambiguous
  ),
  paste0(
    "- Resolved to unique scRNA-seq GSM: ",
    resolved_gse207493
  ),
  paste0(
    "- Final matched objects: ",
    sum(
      resolved_mapping$resolved_match_status ==
        "MATCHED"
    ),
    "/",
    nrow(resolved_mapping)
  ),
  paste0(
    "- Remaining ambiguous objects: ",
    remaining_ambiguous
  ),
  paste0(
    "- Remaining unmatched objects: ",
    remaining_unmatched
  ),
  "",
  "## Clinical candidates",
  "",
  paste0(
    "- GSE207493 patient ID candidates: ",
    nrow(patient_overrides)
  ),
  paste0(
    "- Accepted response candidates: ",
    nrow(response_overrides)
  ),
  paste0(
    "- Total proposed Phase 4E override rows: ",
    nrow(proposed_overrides)
  ),
  "",
  "## Response audit",
  "",
  paste(
    capture.output(
      print(
        response_status_counts,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "## Treatment audit",
  "",
  paste(
    capture.output(
      print(
        treatment_status_counts,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "## Safety boundary",
  "",
  paste0(
    "Phase 4E inputs and outputs were not modified. ",
    "The generated override file is an independent ",
    "candidate file for review and controlled rerun."
  ),
  "",
  "## Status",
  "",
  if (core_fail) {
    "Core Phase 4G QC failed."
  } else {
    paste0(
      "Phase 4G structural QC passed. ",
      "GSE207493 assay-modality ambiguity is resolved."
    )
  }
)

writeLines(
  report,
  file.path(
    out_dir,
    "PHASE4G_REPORT.md"
  ),
  useBytes = TRUE
)

status_file <- file.path(
  out_dir,
  if (core_fail) {
    "PHASE4G_INCOMPLETE"
  } else {
    "PHASE4G_COMPLETE"
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
      "initial_gse207493_ambiguous=",
      initial_gse207493_ambiguous
    ),
    paste0(
      "gse207493_resolved=",
      resolved_gse207493
    ),
    paste0(
      "final_matched=",
      sum(
        resolved_mapping$resolved_match_status ==
          "MATCHED"
      )
    ),
    paste0(
      "remaining_ambiguous=",
      remaining_ambiguous
    ),
    paste0(
      "remaining_unmatched=",
      remaining_unmatched
    ),
    paste0(
      "patient_id_candidates=",
      nrow(patient_overrides)
    ),
    paste0(
      "response_candidates=",
      nrow(response_overrides)
    ),
    paste0(
      "proposed_override_rows=",
      nrow(proposed_overrides)
    )
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
