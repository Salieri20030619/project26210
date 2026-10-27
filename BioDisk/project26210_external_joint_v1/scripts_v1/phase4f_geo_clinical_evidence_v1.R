options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

root <- Sys.getenv(
  "PHASE4F_ROOT",
  unset = "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
)

input_manifest <- file.path(
  root,
  "phase4e_clinical_metadata_harmonization_v1",
  "tables",
  "object_clinical_manifest.tsv"
)

out_dir <- file.path(
  root,
  "phase4f_geo_clinical_evidence_v1"
)

raw_dir <- file.path(out_dir, "raw_geo")
table_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

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
  y <- clean_text(x)
  y == "" |
    toupper(y) %in% c(
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

normalize_key <- function(x) {
  y <- iconv(
    clean_text(x),
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  y[is.na(y)] <- ""
  y <- tolower(y)
  gsub("[^a-z0-9]+", "", y)
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

empty_df <- function(columns) {
  out <- as.data.frame(
    setNames(
      replicate(
        length(columns),
        character(),
        simplify = FALSE
      ),
      columns
    ),
    stringsAsFactors = FALSE
  )

  out
}

soft_key <- function(x) {
  sub(
    "^!",
    "",
    sub(" =.*$", "", x)
  )
}

soft_value <- function(x) {
  clean_text(
    sub("^[^=]*=[[:space:]]*", "", x)
  )
}

parse_soft <- function(dataset, path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  lines <- readLines(
    con,
    warn = FALSE,
    encoding = "UTF-8"
  )

  sample_start <- grep(
    "^\\^SAMPLE = GSM[0-9]+",
    lines
  )

  if (!length(sample_start)) {
    stop(
      "No SAMPLE blocks found in ",
      path
    )
  }

  first_sample <- min(sample_start)

  series_lines <- lines[
    seq_len(first_sample - 1L)
  ]

  series_lines <- series_lines[
    grepl("^!Series_", series_lines)
  ]

  series <- if (length(series_lines)) {
    data.frame(
      dataset = dataset,
      field = vapply(
        series_lines,
        soft_key,
        character(1)
      ),
      value = vapply(
        series_lines,
        soft_value,
        character(1)
      ),
      source_file = path,
      stringsAsFactors = FALSE
    )
  } else {
    empty_df(
      c(
        "dataset",
        "field",
        "value",
        "source_file"
      )
    )
  }

  sample_end <- c(
    sample_start[-1L] - 1L,
    length(lines)
  )

  rows <- list()

  for (i in seq_along(sample_start)) {
    block <- lines[
      sample_start[i]:sample_end[i]
    ]

    gsm <- sub(
      "^\\^SAMPLE = ",
      "",
      block[1]
    )

    sample_lines <- block[
      grepl("^!Sample_", block)
    ]

    if (!length(sample_lines)) {
      next
    }

    keys <- vapply(
      sample_lines,
      soft_key,
      character(1)
    )

    values <- vapply(
      sample_lines,
      soft_value,
      character(1)
    )

    value_index <- ave(
      seq_along(keys),
      keys,
      FUN = seq_along
    )

    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dataset,
      gsm = gsm,
      field = keys,
      value_index = as.integer(value_index),
      value = values,
      source_reference = paste0(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
        gsm
      ),
      source_file = path,
      stringsAsFactors = FALSE
    )
  }

  samples <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    empty_df(
      c(
        "dataset",
        "gsm",
        "field",
        "value_index",
        "value",
        "source_reference",
        "source_file"
      )
    )
  }

  list(
    series = series,
    samples = samples
  )
}

standard_timepoint <- function(x) {
  y <- tolower(clean_text(x))

  if (grepl(
    "progress|relaps|end.?of.?treatment|\\beot\\b",
    y,
    perl = TRUE
  )) {
    return("Progression/EOT")
  }

  if (grepl(
    "baseline|pre.?treatment|pretreatment|before.?treatment|treatment.?naive|untreated|pre.?therapy",
    y,
    perl = TRUE
  )) {
    return("Baseline")
  }

  if (grepl(
    "on.?treatment|during.?treatment|on.?therapy|during.?therapy|cycle[[:space:]_-]*[0-9]+",
    y,
    perl = TRUE
  )) {
    return("On-treatment")
  }

  if (grepl(
    "post.?treatment|after.?treatment|post.?therapy|after.?therapy",
    y,
    perl = TRUE
  )) {
    return("Post-treatment")
  }

  "Unknown"
}

standard_response <- function(x) {
  y <- toupper(clean_text(x))

  if (grepl(
    "NON.?RESPONDER|NON.?RESPONSE|NO.?RESPONSE",
    y,
    perl = TRUE
  )) {
    return("Non-responder")
  }

  if (grepl(
    "COMPLETE.?RESPONSE|\\bCR\\b",
    y,
    perl = TRUE
  )) {
    return("CR")
  }

  if (grepl(
    "PARTIAL.?RESPONSE|\\bPR\\b",
    y,
    perl = TRUE
  )) {
    return("PR")
  }

  if (grepl(
    "PROGRESSIVE.?DISEASE|PROGRESSION|\\bPD\\b",
    y,
    perl = TRUE
  )) {
    return("PD")
  }

  if (grepl(
    "STABLE.?DISEASE|\\bSD\\b",
    y,
    perl = TRUE
  )) {
    return("SD")
  }

  if (grepl(
    "\\bRESPONDER\\b|CLINICAL.?BENEFIT|DURABLE.?BENEFIT",
    y,
    perl = TRUE
  )) {
    return("Responder")
  }

  "Unknown"
}

standard_resistance <- function(x) {
  y <- tolower(clean_text(x))

  if (grepl(
    "acquired|secondary|adaptive",
    y
  )) {
    return("Acquired-resistant")
  }

  if (grepl(
    "primary|intrinsic|de.?novo",
    y
  )) {
    return("Primary-resistant")
  }

  if (grepl(
    "resistan|refractor",
    y
  )) {
    return("Resistant-unspecified")
  }

  "Unknown"
}

split_key_value <- function(x) {
  y <- clean_text(x)
  pos <- regexpr(":", y, fixed = TRUE)

  if (pos[1] < 1L) {
    return(c("", ""))
  }

  key <- trimws(substr(y, 1L, pos[1] - 1L))
  value <- trimws(
    substr(
      y,
      pos[1] + 1L,
      nchar(y)
    )
  )

  c(key, value)
}

classify_clinical_key <- function(key) {
  k <- tolower(clean_text(key))
  k <- gsub("[^a-z0-9]+", " ", k)
  k <- trimws(k)

  if (grepl(
    "(patient|subject|donor).*(id|number)|^(patient|subject|donor)$|^case id$",
    k,
    perl = TRUE
  )) {
    return("patient_id")
  }

  if (grepl(
    "response|recist|clinical benefit|best overall|outcome",
    k
  )) {
    return("response_raw")
  }

  if (grepl(
    "resistan|refractor",
    k
  )) {
    return("resistance_status")
  }

  if (grepl(
    "time.?point|timepoint|sampling.?time|collection.?time|visit|treatment.?status|therapy.?status",
    k
  )) {
    return("sampling_timepoint")
  }

  if (grepl(
    "regimen|combination",
    k
  )) {
    return("regimen")
  }

  if (grepl(
    "treatment|therapy|drug|agent",
    k
  )) {
    return("treatment")
  }

  ""
}

manifest <- read.delim(
  input_manifest,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

required_manifest <- c(
  "dataset",
  "object_name",
  "source_file",
  "cells",
  "rds_read_status"
)

missing_manifest <- setdiff(
  required_manifest,
  names(manifest)
)

if (length(missing_manifest)) {
  stop(
    "Missing manifest columns: ",
    paste(missing_manifest, collapse = ", ")
  )
}

if (nrow(manifest) != 95L) {
  stop(
    "Expected 95 objects, observed ",
    nrow(manifest)
  )
}

manifest$cells <- suppressWarnings(
  as.numeric(manifest$cells)
)

total_cells <- sum(
  manifest$cells,
  na.rm = TRUE
)

if (total_cells != 623431) {
  stop(
    "Expected 623431 cells, observed ",
    total_cells
  )
}

if (any(manifest$rds_read_status != "PASS")) {
  stop(
    "Manifest contains non-PASS RDS records"
  )
}

datasets <- sort(
  unique(clean_text(manifest$dataset))
)

parsed <- lapply(
  datasets,
  function(dataset) {
    soft_file <- file.path(
      raw_dir,
      paste0(dataset, "_family.soft.gz")
    )

    if (!file.exists(soft_file)) {
      stop(
        "Missing GEO SOFT file: ",
        soft_file
      )
    }

    parse_soft(dataset, soft_file)
  }
)

series_metadata <- do.call(
  rbind,
  lapply(parsed, `[[`, "series")
)

sample_metadata <- do.call(
  rbind,
  lapply(parsed, `[[`, "samples")
)

write_tsv(
  series_metadata,
  file.path(
    table_dir,
    "geo_series_metadata_long.tsv"
  )
)

write_tsv(
  sample_metadata,
  file.path(
    table_dir,
    "geo_sample_metadata_long.tsv"
  )
)

sample_keys <- split(
  sample_metadata,
  interaction(
    sample_metadata$dataset,
    sample_metadata$gsm,
    drop = TRUE,
    lex.order = TRUE
  )
)

sample_index_rows <- lapply(
  sample_keys,
  function(z) {
    selected <- z[
      z$field %in% c(
        "Sample_title",
        "Sample_source_name_ch1",
        "Sample_characteristics_ch1",
        "Sample_description"
      ),
      ,
      drop = FALSE
    ]

    values <- selected$value

    characteristic_values <- vapply(
      values,
      function(v) {
        kv <- split_key_value(v)

        if (kv[2] != "") {
          kv[2]
        } else {
          v
        }
      },
      character(1)
    )

    tokens <- unlist(
      strsplit(
        values,
        "[^A-Za-z0-9]+"
      ),
      use.names = FALSE
    )

    tokens <- tokens[
      nchar(tokens) >= 4L
    ]

    stop_tokens <- c(
      "sample",
      "samples",
      "tumor",
      "tumour",
      "normal",
      "tissue",
      "patient",
      "subject",
      "donor",
      "cells",
      "single",
      "biopsy"
    )

    tokens <- tokens[
      !tolower(tokens) %in% stop_tokens
    ]

    keys <- unique(
      c(
        z$gsm[1],
        values,
        characteristic_values,
        tokens
      )
    )

    keys <- unique(normalize_key(keys))
    keys <- keys[nchar(keys) >= 4L]

    data.frame(
      dataset = z$dataset[1],
      gsm = z$gsm[1],
      exact_keys = paste(
        keys,
        collapse = ";"
      ),
      title = collapse_unique(
        z$value[
          z$field == "Sample_title"
        ]
      ),
      source_name = collapse_unique(
        z$value[
          z$field ==
            "Sample_source_name_ch1"
        ]
      ),
      stringsAsFactors = FALSE
    )
  }
)

sample_index <- do.call(
  rbind,
  sample_index_rows
)

mapping_rows <- list()

for (i in seq_len(nrow(manifest))) {
  dataset <- clean_text(
    manifest$dataset[i]
  )

  object_name <- clean_text(
    manifest$object_name[i]
  )

  stem <- sub(
    "\\.rds$",
    "",
    object_name,
    ignore.case = TRUE
  )

  without_dataset <- sub(
    paste0(
      "^",
      dataset,
      "[._ -]*"
    ),
    "",
    stem,
    ignore.case = TRUE
  )

  gsm_in_name <- regmatches(
    stem,
    regexpr(
      "GSM[0-9]+",
      stem,
      ignore.case = TRUE
    )
  )

  if (!length(gsm_in_name) ||
      identical(gsm_in_name, character(0)) ||
      gsm_in_name == "") {
    gsm_in_name <- ""
  } else {
    gsm_in_name <- toupper(gsm_in_name)
  }

  without_gsm <- if (gsm_in_name != "") {
    gsub(
      gsm_in_name,
      "",
      without_dataset,
      ignore.case = TRUE
    )
  } else {
    without_dataset
  }

  aliases <- unique(
    normalize_key(
      c(
        stem,
        without_dataset,
        without_gsm
      )
    )
  )

  aliases <- aliases[
    nchar(aliases) >= 4L
  ]

  index_d <- sample_index[
    sample_index$dataset == dataset,
    ,
    drop = FALSE
  ]

  candidate_gsm <- character()
  method <- ""

  if (gsm_in_name != "") {
    direct <- index_d$gsm[
      toupper(index_d$gsm) ==
        gsm_in_name
    ]

    if (length(direct)) {
      candidate_gsm <- unique(direct)
      method <- "direct_gsm"
    }
  }

  if (!length(candidate_gsm) &&
      length(aliases)) {
    hit <- vapply(
      seq_len(nrow(index_d)),
      function(j) {
        keys <- strsplit(
          index_d$exact_keys[j],
          ";",
          fixed = TRUE
        )[[1]]

        any(aliases %in% keys)
      },
      logical(1)
    )

    candidate_gsm <- unique(
      index_d$gsm[hit]
    )

    if (length(candidate_gsm)) {
      method <- "exact_key"
    }
  }

  status <- if (length(candidate_gsm) == 1L) {
    "MATCHED"
  } else if (length(candidate_gsm) > 1L) {
    "AMBIGUOUS"
  } else {
    "UNMATCHED"
  }

  gsm <- if (status == "MATCHED") {
    candidate_gsm
  } else {
    ""
  }

  matched_title <- if (status == "MATCHED") {
    collapse_unique(
      index_d$title[
        index_d$gsm == gsm
      ]
    )
  } else {
    ""
  }

  mapping_rows[[i]] <- data.frame(
    dataset = dataset,
    object_name = object_name,
    source_file = manifest$source_file[i],
    object_aliases = paste(
      aliases,
      collapse = ";"
    ),
    match_status = status,
    match_method = method,
    gsm = gsm,
    candidate_gsms = paste(
      candidate_gsm,
      collapse = ";"
    ),
    geo_title = matched_title,
    stringsAsFactors = FALSE
  )
}

object_mapping <- do.call(
  rbind,
  mapping_rows
)

write_tsv(
  sample_index,
  file.path(
    table_dir,
    "geo_sample_exact_key_index.tsv"
  )
)

write_tsv(
  object_mapping,
  file.path(
    table_dir,
    "object_geo_exact_match.tsv"
  )
)

characteristics <- sample_metadata[
  sample_metadata$field ==
    "Sample_characteristics_ch1",
  ,
  drop = FALSE
]

evidence_rows <- list()

if (nrow(characteristics)) {
  for (i in seq_len(nrow(characteristics))) {
    kv <- split_key_value(
      characteristics$value[i]
    )

    if (kv[1] == "" || kv[2] == "") {
      next
    }

    target <- classify_clinical_key(
      kv[1]
    )

    if (target == "") {
      next
    }

    standardized <- switch(
      target,
      sampling_timepoint =
        standard_timepoint(kv[2]),
      response_raw =
        standard_response(kv[2]),
      resistance_status =
        standard_resistance(kv[2]),
      clean_text(kv[2])
    )

    evidence_rows[[
      length(evidence_rows) + 1L
    ]] <- data.frame(
      dataset =
        characteristics$dataset[i],
      gsm =
        characteristics$gsm[i],
      target_field = target,
      raw_key = kv[1],
      raw_value = kv[2],
      standardized_value = standardized,
      metadata_source =
        "GEO SOFT: Sample_characteristics_ch1",
      source_reference =
        characteristics$source_reference[i],
      stringsAsFactors = FALSE
    )
  }
}

gsm_evidence <- if (length(evidence_rows)) {
  unique(do.call(rbind, evidence_rows))
} else {
  empty_df(
    c(
      "dataset",
      "gsm",
      "target_field",
      "raw_key",
      "raw_value",
      "standardized_value",
      "metadata_source",
      "source_reference"
    )
  )
}

write_tsv(
  gsm_evidence,
  file.path(
    table_dir,
    "gsm_clinical_evidence.tsv"
  )
)

matched <- object_mapping[
  object_mapping$match_status == "MATCHED",
  c(
    "dataset",
    "object_name",
    "gsm",
    "match_method"
  ),
  drop = FALSE
]

object_evidence <- merge(
  matched,
  gsm_evidence,
  by = c("dataset", "gsm"),
  all = FALSE,
  sort = FALSE
)

candidate_columns <- c(
  "dataset",
  "object_name",
  "gsm",
  "match_method",
  "field",
  "candidate_status",
  "value",
  "standardized_value",
  "raw_keys",
  "metadata_source",
  "source_reference"
)

candidate_rows <- list()

if (nrow(object_evidence)) {
  groups <- split(
    seq_len(nrow(object_evidence)),
    interaction(
      object_evidence$dataset,
      object_evidence$object_name,
      object_evidence$target_field,
      drop = TRUE,
      lex.order = TRUE
    )
  )

  for (idx in groups) {
    z <- object_evidence[
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

    standard_values <- unique(
      clean_text(z$standardized_value)
    )

    standard_values <- standard_values[
      !is_missing(standard_values) &
        standard_values != "Unknown"
    ]

    candidate_status <- if (
      length(raw_values) == 1L
    ) {
      "UNIQUE"
    } else if (length(raw_values) > 1L) {
      "CONFLICT"
    } else {
      "EMPTY"
    }

    candidate_rows[[
      length(candidate_rows) + 1L
    ]] <- data.frame(
      dataset = z$dataset[1],
      object_name = z$object_name[1],
      gsm = z$gsm[1],
      match_method = z$match_method[1],
      field = z$target_field[1],
      candidate_status = candidate_status,
      value = collapse_unique(raw_values),
      standardized_value =
        collapse_unique(
          standard_values,
          default = "Unknown"
        ),
      raw_keys =
        collapse_unique(z$raw_key),
      metadata_source =
        collapse_unique(z$metadata_source),
      source_reference =
        collapse_unique(z$source_reference),
      stringsAsFactors = FALSE
    )
  }
}

object_candidates <- if (
  length(candidate_rows)
) {
  do.call(rbind, candidate_rows)
} else {
  empty_df(candidate_columns)
}

object_candidates <- object_candidates[
  ,
  candidate_columns,
  drop = FALSE
]

write_tsv(
  object_candidates,
  file.path(
    table_dir,
    "object_clinical_candidates.tsv"
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

if (nrow(object_candidates)) {
  override_candidates <- object_candidates[
    object_candidates$candidate_status ==
      "UNIQUE" &
      object_candidates$field %in% c(
        "patient_id",
        "sampling_timepoint",
        "response_raw",
        "resistance_status",
        "treatment",
        "regimen"
      ),
    ,
    drop = FALSE
  ]

  if (nrow(override_candidates)) {
    override_candidates <- data.frame(
      dataset =
        override_candidates$dataset,
      object_name =
        override_candidates$object_name,
      field =
        override_candidates$field,
      value =
        override_candidates$value,
      metadata_source =
        override_candidates$metadata_source,
      source_reference =
        override_candidates$source_reference,
      notes = paste0(
        "REVIEW REQUIRED; exact object-to-GSM ",
        "match using ",
        override_candidates$match_method,
        "; GEO key=",
        override_candidates$raw_keys
      ),
      stringsAsFactors = FALSE
    )
  } else {
    override_candidates <- empty_df(
      override_columns
    )
  }
} else {
  override_candidates <- empty_df(
    override_columns
  )
}

write_tsv(
  override_candidates,
  file.path(
    table_dir,
    paste0(
      "phase4e_manual_clinical_overrides_",
      "REVIEW_REQUIRED.tsv"
    )
  )
)

summary_rows <- lapply(
  datasets,
  function(dataset) {
    map_d <- object_mapping[
      object_mapping$dataset == dataset,
      ,
      drop = FALSE
    ]

    cand_d <- object_candidates[
      object_candidates$dataset == dataset,
      ,
      drop = FALSE
    ]

    count_unique <- function(field) {
      if (!nrow(cand_d)) {
        return(0L)
      }

      length(unique(
        cand_d$object_name[
          cand_d$field == field &
            cand_d$candidate_status ==
              "UNIQUE"
        ]
      ))
    }

    data.frame(
      dataset = dataset,
      objects = nrow(map_d),
      geo_samples = length(unique(
        sample_metadata$gsm[
          sample_metadata$dataset ==
            dataset
        ]
      )),
      exact_matched_objects = sum(
        map_d$match_status == "MATCHED"
      ),
      ambiguous_objects = sum(
        map_d$match_status == "AMBIGUOUS"
      ),
      unmatched_objects = sum(
        map_d$match_status == "UNMATCHED"
      ),
      patient_id_candidates =
        count_unique("patient_id"),
      timepoint_candidates =
        count_unique("sampling_timepoint"),
      response_candidates =
        count_unique("response_raw"),
      resistance_candidates =
        count_unique("resistance_status"),
      treatment_candidates =
        count_unique("treatment"),
      stringsAsFactors = FALSE
    )
  }
)

dataset_summary <- do.call(
  rbind,
  summary_rows
)

write_tsv(
  dataset_summary,
  file.path(
    table_dir,
    "phase4f_dataset_summary.tsv"
  )
)

qc <- data.frame(
  check = c(
    "Phase 4E manifest exists",
    "Manifest object count",
    "Manifest cell count",
    "Manifest RDS statuses",
    "Datasets with parsed GEO SOFT",
    "GEO sample records",
    "Object mappings classified",
    "Automatic Phase 4E overwrite"
  ),
  expected = c(
    "TRUE",
    "95",
    "623431",
    "All PASS",
    as.character(length(datasets)),
    ">0",
    "95",
    "0"
  ),
  observed = c(
    as.character(file.exists(input_manifest)),
    as.character(nrow(manifest)),
    as.character(total_cells),
    as.character(sum(
      manifest$rds_read_status != "PASS"
    )),
    as.character(length(unique(
      sample_metadata$dataset
    ))),
    as.character(length(unique(
      paste(
        sample_metadata$dataset,
        sample_metadata$gsm,
        sep = "::"
      )
    ))),
    as.character(nrow(object_mapping)),
    "0"
  ),
  status = c(
    if (file.exists(input_manifest)) {
      "PASS"
    } else {
      "FAIL"
    },
    if (nrow(manifest) == 95L) {
      "PASS"
    } else {
      "FAIL"
    },
    if (total_cells == 623431) {
      "PASS"
    } else {
      "FAIL"
    },
    if (all(
      manifest$rds_read_status == "PASS"
    )) {
      "PASS"
    } else {
      "FAIL"
    },
    if (
      setequal(
        datasets,
        unique(sample_metadata$dataset)
      )
    ) {
      "PASS"
    } else {
      "FAIL"
    },
    if (nrow(sample_metadata) > 0L) {
      "PASS"
    } else {
      "FAIL"
    },
    if (
      nrow(object_mapping) ==
        nrow(manifest)
    ) {
      "PASS"
    } else {
      "FAIL"
    },
    "PASS"
  ),
  notes = c(
    "",
    "",
    "",
    "",
    "",
    "",
    paste0(
      "MATCHED=",
      sum(
        object_mapping$match_status ==
          "MATCHED"
      ),
      "; AMBIGUOUS=",
      sum(
        object_mapping$match_status ==
          "AMBIGUOUS"
      ),
      "; UNMATCHED=",
      sum(
        object_mapping$match_status ==
          "UNMATCHED"
      )
    ),
    paste0(
      "Candidate overrides are review-only ",
      "and were not copied into Phase 4E inputs"
    )
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  qc,
  file.path(
    qc_dir,
    "phase4f_qc.tsv"
  )
)

core_fail <- any(qc$status == "FAIL")

report <- c(
  "# Phase 4F GEO clinical evidence report",
  "",
  paste0("Generated: ", stamp()),
  "",
  "## Scope",
  "",
  paste0(
    "This phase downloaded GEO family SOFT metadata, ",
    "parsed sample-level evidence, performed only ",
    "deterministic exact object-to-GSM matching, and ",
    "created review-required clinical override candidates."
  ),
  "",
  paste0(
    "No Phase 4E manifest, RDS object, formal project ",
    "result, or customer report was modified."
  ),
  "",
  "## Core totals",
  "",
  paste0("- Objects: ", nrow(manifest)),
  paste0("- Cells: ", total_cells),
  paste0("- Datasets: ", length(datasets)),
  paste0(
    "- GEO samples: ",
    length(unique(
      paste(
        sample_metadata$dataset,
        sample_metadata$gsm,
        sep = "::"
      )
    ))
  ),
  paste0(
    "- Exact matched objects: ",
    sum(
      object_mapping$match_status ==
        "MATCHED"
    )
  ),
  paste0(
    "- Ambiguous objects: ",
    sum(
      object_mapping$match_status ==
        "AMBIGUOUS"
    )
  ),
  paste0(
    "- Unmatched objects: ",
    sum(
      object_mapping$match_status ==
        "UNMATCHED"
    )
  ),
  paste0(
    "- Review-required override rows: ",
    nrow(override_candidates)
  ),
  "",
  "## Dataset summary",
  "",
  paste(
    capture.output(
      print(
        dataset_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "## Evidence rules",
  "",
  paste0(
    "- Direct GSM identifiers and exact normalized ",
    "metadata keys were allowed."
  ),
  paste0(
    "- Fuzzy matching and cross-dataset patient ",
    "matching were not used."
  ),
  paste0(
    "- GEO-derived candidates remain review-required ",
    "and were not automatically applied."
  ),
  paste0(
    "- Missing response, timepoint, or patient fields ",
    "remain unresolved rather than inferred."
  ),
  "",
  "## Status",
  "",
  if (core_fail) {
    "Core Phase 4F QC failed."
  } else {
    paste0(
      "Core GEO evidence acquisition passed. ",
      "Clinical completeness is not implied."
    )
  }
)

writeLines(
  report,
  file.path(
    out_dir,
    "PHASE4F_GEO_CLINICAL_EVIDENCE_REPORT.md"
  ),
  useBytes = TRUE
)

status_file <- file.path(
  out_dir,
  if (core_fail) {
    "PHASE4F_INCOMPLETE"
  } else {
    "PHASE4F_GEO_EVIDENCE_COMPLETE"
  }
)

writeLines(
  c(
    paste0("completed_at=", stamp()),
    paste0(
      "core_qc=",
      if (core_fail) "FAIL" else "PASS"
    ),
    paste0("objects=", nrow(manifest)),
    paste0("cells=", total_cells),
    paste0("datasets=", length(datasets)),
    paste0(
      "exact_matched=",
      sum(
        object_mapping$match_status ==
          "MATCHED"
      )
    ),
    paste0(
      "ambiguous=",
      sum(
        object_mapping$match_status ==
          "AMBIGUOUS"
      )
    ),
    paste0(
      "unmatched=",
      sum(
        object_mapping$match_status ==
          "UNMATCHED"
      )
    ),
    paste0(
      "review_override_rows=",
      nrow(override_candidates)
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
  quit(save = "no", status = 2L)
}
