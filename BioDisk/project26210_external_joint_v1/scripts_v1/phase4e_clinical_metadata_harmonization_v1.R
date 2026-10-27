#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, timeout = 1800)

WORK_ROOT <- normalizePath(
  Sys.getenv("PHASE4E_WORK_ROOT", "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"),
  winslash = "/", mustWork = FALSE
)
AUDIT_DIR <- file.path(WORK_ROOT, "phase4d_95_object_compatibility_audit_v1")
OUT_DIR <- file.path(WORK_ROOT, "phase4e_clinical_metadata_harmonization_v1")
TABLE_DIR <- file.path(OUT_DIR, "tables")
QC_DIR <- file.path(OUT_DIR, "qc")
CACHE_DIR <- file.path(OUT_DIR, "cache", "geo")
EXPECTED_OBJECTS <- 95L
EXPECTED_CELLS <- 623431
EXPECTED_DATASETS <- c(
  "GSE178481", "GSE202374", "GSE207493",
  "GSE264586", "GSE294109"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", stamp(), paste0(..., collapse = "")))
}
trim <- function(x) trimws(as.character(x))
is_missing <- function(x) {
  is.na(x) | trim(x) == "" |
    trim(x) %in% c("NA", "N/A", "Unknown", "unknown", "NULL")
}
first_nonmissing <- function(x, default = "Unknown") {
  x <- unique(trim(x[!is_missing(x)]))
  if (length(x) == 1L) x else default
}
collapse_unique <- function(x, sep = "; ", default = "Unknown") {
  x <- sort(unique(trim(x[!is_missing(x)])))
  if (length(x)) paste(x, collapse = sep) else default
}
norm_name <- function(x) {
  gsub("[^a-z0-9]+", "_", tolower(trim(x)))
}
norm_token <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(trim(x)))
}
safe_rate <- function(n, d) {
  if (d > 0) n / d else NA_real_
}
write_tsv <- function(x, path) {
  write.table(
    x, path,
    sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
  log_msg("WROTE ", path, " [", nrow(x), " rows]")
}

unresolved <- data.frame(
  severity = character(),
  dataset = character(),
  object_name = character(),
  field = character(),
  reason = character(),
  observed = character(),
  suggested_action = character(),
  stringsAsFactors = FALSE
)

add_unresolved <- function(
    severity, dataset, object_name, field, reason,
    observed = "", suggested_action = "") {
  unresolved <<- rbind(
    unresolved,
    data.frame(
      severity = severity,
      dataset = dataset,
      object_name = object_name,
      field = field,
      reason = reason,
      observed = observed,
      suggested_action = suggested_action,
      stringsAsFactors = FALSE
    )
  )
}

qc <- data.frame(
  test = character(),
  expected = character(),
  observed = character(),
  status = character(),
  details = character(),
  stringsAsFactors = FALSE
)

add_qc <- function(test, expected, observed, status, details = "") {
  qc <<- rbind(
    qc,
    data.frame(
      test = test,
      expected = as.character(expected),
      observed = as.character(observed),
      status = status,
      details = details,
      stringsAsFactors = FALSE
    )
  )
}

find_audit_file <- function() {
  preferred <- file.path(
    AUDIT_DIR, "qc", "object_level_audit.tsv"
  )
  if (file.exists(preferred)) return(preferred)

  hits <- list.files(
    WORK_ROOT,
    pattern = "^object_level_audit\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  hits <- hits[
    !grepl(
      "phase4e_clinical_metadata_harmonization",
      hits,
      fixed = TRUE
    )
  ]

  if (!length(hits)) {
    stop("Cannot find object_level_audit.tsv under ", WORK_ROOT)
  }

  info <- file.info(hits)
  hits[order(info$mtime, decreasing = TRUE)][1]
}

pick_col <- function(
    df, candidates, required = FALSE, label = "column") {
  nn <- norm_name(names(df))
  cc <- norm_name(candidates)
  hit <- which(nn %in% cc)

  if (length(hit)) return(names(df)[hit[1]])

  if (required) {
    stop(
      "Cannot identify ", label,
      ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  NA_character_
}

resolve_source_file <- function(
    path_value, object_name, audit_file, rds_index) {
  candidates <- character()

  if (!is_missing(path_value)) {
    p <- trim(path_value)
    candidates <- c(
      candidates,
      p,
      file.path(WORK_ROOT, p),
      file.path(dirname(audit_file), p)
    )
  }

  obj <- trim(object_name)
  if (!is_missing(obj)) {
    candidates <- c(
      candidates,
      file.path(WORK_ROOT, obj),
      file.path(WORK_ROOT, paste0(obj, ".rds"))
    )
  }

  candidates <- unique(candidates)
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 1L) {
    return(normalizePath(
      existing,
      winslash = "/",
      mustWork = TRUE
    ))
  }

  if (length(existing) > 1L) {
    same <- unique(normalizePath(
      existing,
      winslash = "/",
      mustWork = TRUE
    ))
    if (length(same) == 1L) return(same)
    return(NA_character_)
  }

  bases <- unique(c(
    basename(trim(path_value)),
    basename(obj),
    paste0(basename(obj), ".rds")
  ))
  bases <- bases[!is_missing(bases)]

  hit <- rds_index[basename(rds_index) %in% bases]
  if (length(hit) == 1L) return(hit)

  NA_character_
}

extract_meta <- function(obj) {
  if (isS4(obj) && "meta.data" %in% slotNames(obj)) {
    return(as.data.frame(slot(obj, "meta.data")))
  }
  if (is.list(obj) && !is.null(obj$meta.data)) {
    return(as.data.frame(obj$meta.data))
  }
  NULL
}

find_values <- function(df, patterns) {
  if (is.null(df) || !ncol(df)) {
    return(list(
      value = "Unknown",
      columns = "",
      conflict = FALSE
    ))
  }

  nn <- norm_name(names(df))
  nn_unwrapped <- sub("^characteristic_", "", nn)

  idx <- unique(unlist(lapply(
    patterns,
    function(p) {
      union(
        grep(p, nn, perl = TRUE),
        grep(p, nn_unwrapped, perl = TRUE)
      )
    }
  )))

  if (!length(idx)) {
    return(list(
      value = "Unknown",
      columns = "",
      conflict = FALSE
    ))
  }

  vals <- unique(unlist(
    lapply(
      df[idx],
      function(z) trim(z[!is_missing(z)])
    ),
    use.names = FALSE
  ))
  vals <- vals[!is_missing(vals)]

  list(
    value = if (length(vals) == 1L) vals else "Unknown",
    columns = paste(names(df)[idx], collapse = ";"),
    conflict = length(vals) > 1L,
    observed = paste(head(vals, 20), collapse = " | ")
  )
}

FIELD_PATTERNS <- list(
  gsm = c("^gsm$", "^geo_accession$"),
  sample_id = c(
    "^sample$", "^sample_id$", "^sampleid$",
    "^specimen$", "^specimen_id$", "^orig_ident$"
  ),
  patient_id = c(
    "^patient$", "^patient_id$", "^patientid$",
    "^donor$", "^donor_id$",
    "^subject$", "^subject_id$",
    "^case$", "^case_id$",
    "^individual$", "^individual_id$"
  ),
  tissue = c(
    "^tissue$", "^tissue_type$", "^source_name_ch1$"
  ),
  tumor_normal = c(
    "^tumor_normal$", "^tumour_normal$",
    "^sample_type$", "^tissue_status$"
  ),
  primary_metastatic = c(
    "^primary_metastatic$", "^lesion_type$",
    "^disease_site_type$"
  ),
  anatomical_site = c(
    "^anatomical_site$", "^anatomic_site$",
    "^biopsy_site$", "^organ$", "^site$"
  ),
  treatment = c(
    "^treatment$", "^therapy$",
    "^treatment_status$", "^therapy_status$"
  ),
  treatment_class = c(
    "^treatment_class$", "^therapy_class$"
  ),
  ici_exposure = c(
    "^ici$", "^ici_exposure$",
    "^immunotherapy$", "^checkpoint_inhibitor$"
  ),
  regimen = c(
    "^regimen$", "^treatment_regimen$",
    "^therapy_regimen$"
  ),
  sampling_timepoint = c(
    "^timepoint$", "^time_point$",
    "^sampling_timepoint$", "^visit$",
    "^collection_time$", "^biopsy_timepoint$"
  ),
  response_raw = c(
    "^response$", "^response_status$",
    "^best_response$", "^clinical_response$",
    "^outcome$"
  ),
  resistance_status = c(
    "^resistance$", "^resistance_status$",
    "^resistance_type$"
  ),
  resistance_definition = c(
    "^resistance_definition$",
    "^resistance_evidence$"
  ),
  title = c("^title$", "^sample_title$")
)

extract_fields <- function(df) {
  lapply(
    FIELD_PATTERNS,
    function(p) find_values(df, p)
  )
}

reconcile <- function(a, b, field, dataset, object_name) {
  av <- if (is.null(a$value)) "Unknown" else a$value
  bv <- if (is.null(b$value)) "Unknown" else b$value

  if (isTRUE(a$conflict)) {
    add_unresolved(
      "WARN", dataset, object_name, field,
      "Multiple values within RDS meta.data",
      a$observed
    )
  }

  if (isTRUE(b$conflict)) {
    add_unresolved(
      "WARN", dataset, object_name, field,
      "Multiple values within GEO/sample metadata",
      b$observed
    )
  }

  if (!is_missing(av) && !is_missing(bv)) {
    if (norm_token(av) == norm_token(bv)) {
      return(list(
        value = av,
        source = "RDS+GEO",
        evidence = "High"
      ))
    }

    add_unresolved(
      "WARN", dataset, object_name, field,
      "RDS and GEO values conflict",
      paste0("RDS=", av, " | GEO=", bv),
      "Resolve from primary paper/supplement"
    )

    return(list(
      value = "Unknown",
      source = "Conflict",
      evidence = "Unresolved"
    ))
  }

  if (!is_missing(av)) {
    return(list(
      value = av,
      source = paste0("RDS:", a$columns),
      evidence = "Medium"
    ))
  }

  if (!is_missing(bv)) {
    return(list(
      value = bv,
      source = paste0("GEO:", b$columns),
      evidence = "Medium"
    ))
  }

  list(
    value = "Unknown",
    source = "Not found",
    evidence = "Missing"
  )
}

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

rbind_fill <- function(...) {
  xs <- list(...)

  if (
    length(xs) == 1L &&
    is.list(xs[[1]]) &&
    !is.data.frame(xs[[1]])
  ) {
    xs <- xs[[1]]
  }

  xs <- xs[vapply(xs, is.data.frame, logical(1))]
  if (!length(xs)) return(data.frame())

  cols <- unique(unlist(lapply(xs, names)))

  xs <- lapply(xs, function(x) {
    miss <- setdiff(cols, names(x))
    for (m in miss) x[[m]] <- NA_character_
    x[cols]
  })

  do.call(rbind, xs)
}

geo_matrix_prefix <- function(gse) {
  sub("[0-9]{3}$", "nnn", gse)
}

download_geo_matrices <- function(gse) {
  local <- list.files(
    WORK_ROOT,
    pattern = paste0(
      gse,
      ".*series_matrix.*\\.txt(\\.gz)?$"
    ),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  local <- local[
    !grepl(
      "phase4e_clinical_metadata_harmonization",
      local,
      fixed = TRUE
    )
  ]

  cached <- list.files(
    CACHE_DIR,
    pattern = paste0(
      "^", gse,
      ".*series_matrix.*\\.txt(\\.gz)?$"
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )

  found <- unique(c(local, cached))
  if (length(found)) return(found)

  base_url <- sprintf(
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/",
      "%s/%s/matrix/"
    ),
    geo_matrix_prefix(gse),
    gse
  )
  index_file <- file.path(
    CACHE_DIR,
    paste0(gse, "_matrix_index.html")
  )

  urls <- character()

  ok <- tryCatch({
    suppressWarnings(download.file(
      base_url,
      index_file,
      mode = "wb",
      quiet = TRUE
    ))
    TRUE
  }, error = function(e) FALSE,
  warning = function(w) FALSE)

  if (ok && file.exists(index_file)) {
    txt <- paste(
      readLines(index_file, warn = FALSE),
      collapse = "\n"
    )
    names_hit <- unique(regmatches(
      txt,
      gregexpr(
        paste0(
          gse,
          "[^\"'<> ]*series_matrix",
          "[^\"'<> ]*\\.txt\\.gz"
        ),
        txt,
        perl = TRUE
      )
    )[[1]])
    names_hit <- names_hit[nzchar(names_hit)]
    urls <- paste0(base_url, names_hit)
  }

  if (!length(urls)) {
    urls <- paste0(
      base_url,
      gse,
      "_series_matrix.txt.gz"
    )
  }

  out <- character()

  for (u in urls) {
    dest <- file.path(CACHE_DIR, basename(u))

    ok2 <- tryCatch({
      suppressWarnings(download.file(
        u, dest,
        mode = "wb",
        quiet = TRUE
      ))
      file.exists(dest) &&
        file.info(dest)$size > 100
    }, error = function(e) FALSE,
    warning = function(w) FALSE)

    if (ok2) {
      out <- c(out, dest)
    } else if (file.exists(dest)) {
      unlink(dest)
    }
  }

  unique(out)
}

parse_geo_matrix_header <- function(path, dataset) {
  con <- if (
    grepl("\\.gz$", path, ignore.case = TRUE)
  ) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(con), add = TRUE)

  lines <- character()

  repeat {
    z <- readLines(con, n = 1L, warn = FALSE)
    if (
      !length(z) ||
      grepl("^!series_matrix_table_begin", z)
    ) {
      break
    }
    if (grepl("^!Sample_", z)) {
      lines <- c(lines, z)
    }
  }

  if (!length(lines)) return(data.frame())

  parsed <- lapply(lines, function(z) {
    x <- strsplit(z, "\\t", fixed = FALSE)[[1]]
    key <- sub("^!Sample_", "", x[1])
    vals <- gsub('^"|"$', "", x[-1])
    list(key = key, vals = vals)
  })

  n <- max(vapply(
    parsed,
    function(x) length(x$vals),
    integer(1)
  ))

  access <- parsed[vapply(
    parsed,
    function(x) x$key == "geo_accession",
    logical(1)
  )]

  if (!length(access)) return(data.frame())

  gsm <- access[[1]]$vals
  if (length(gsm) != n) return(data.frame())

  out <- data.frame(
    dataset = dataset,
    gsm = gsm,
    stringsAsFactors = FALSE
  )
  key_count <- list()

  for (x in parsed) {
    if (
      length(x$vals) != n ||
      x$key == "geo_accession"
    ) {
      next
    }

    key <- norm_name(x$key)

    if (key == "characteristics_ch1") {
      kv <- strsplit(x$vals, ":", fixed = TRUE)

      keys <- vapply(
        kv,
        function(y) {
          if (length(y) >= 2) {
            norm_name(y[1])
          } else {
            "characteristic"
          }
        },
        character(1)
      )

      vals <- vapply(
        kv,
        function(y) {
          if (length(y) >= 2) {
            trim(paste(y[-1], collapse = ":"))
          } else {
            trim(y[1])
          }
        },
        character(1)
      )

      for (k in unique(keys)) {
        col <- paste0("characteristic_", k)
        if (!col %in% names(out)) out[[col]] <- ""

        idx <- which(keys == k)
        out[[col]][idx] <- ifelse(
          out[[col]][idx] == "",
          vals[idx],
          paste(
            out[[col]][idx],
            vals[idx],
            sep = "; "
          )
        )
      }
    } else {
      key_count[[key]] <-
        (key_count[[key]] %||% 0L) + 1L

      col <- if (key_count[[key]] == 1L) {
        key
      } else {
        paste0(key, "_", key_count[[key]])
      }

      out[[col]] <- x$vals
    }
  }

  out
}

download_and_parse_miniml <- function(gse) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    return(data.frame())
  }

  local_xml <- list.files(
    WORK_ROOT,
    pattern = paste0(
      gse,
      ".*family.*\\.xml$"
    ),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  xml_file <- if (length(local_xml)) {
    local_xml[1]
  } else {
    NA_character_
  }

  if (is.na(xml_file)) {
    tgz <- file.path(
      CACHE_DIR,
      paste0(gse, "_family.xml.tgz")
    )

    if (
      !file.exists(tgz) ||
      file.info(tgz)$size < 100
    ) {
      url <- sprintf(
        paste0(
          "https://ftp.ncbi.nlm.nih.gov/geo/series/",
          "%s/%s/miniml/%s_family.xml.tgz"
        ),
        geo_matrix_prefix(gse),
        gse,
        gse
      )

      ok <- tryCatch({
        suppressWarnings(download.file(
          url, tgz,
          mode = "wb",
          quiet = TRUE
        ))
        file.exists(tgz) &&
          file.info(tgz)$size > 100
      }, error = function(e) FALSE,
      warning = function(w) FALSE)

      if (!ok) return(data.frame())
    }

    exdir <- file.path(
      CACHE_DIR,
      paste0(gse, "_miniml")
    )
    dir.create(
      exdir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    tryCatch(
      utils::untar(tgz, exdir = exdir),
      error = function(e) NULL
    )

    xmls <- list.files(
      exdir,
      pattern = "\\.xml$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (!length(xmls)) return(data.frame())
    xml_file <- xmls[1]
  }

  doc <- tryCatch(
    xml2::read_xml(xml_file),
    error = function(e) NULL
  )
  if (is.null(doc)) return(data.frame())

  samples <- xml2::xml_find_all(
    doc,
    ".//*[local-name()='Sample']"
  )
  if (!length(samples)) return(data.frame())

  rows <- lapply(samples, function(s) {
    gsm <- xml2::xml_attr(s, "iid")

    title <- xml2::xml_text(
      xml2::xml_find_first(
        s,
        ".//*[local-name()='Title']"
      )
    )

    source <- xml2::xml_text(
      xml2::xml_find_first(
        s,
        ".//*[local-name()='Source']"
      )
    )

    chars <- xml2::xml_find_all(
      s,
      ".//*[local-name()='Characteristics']"
    )

    row <- data.frame(
      dataset = gse,
      gsm = gsm,
      title = trim(title),
      source_name_ch1 = trim(source),
      stringsAsFactors = FALSE
    )

    if (length(chars)) {
      tags <- norm_name(
        xml2::xml_attr(chars, "tag")
      )
      vals <- trim(xml2::xml_text(chars))
      tags[is_missing(tags)] <- "characteristic"

      for (tag in unique(tags)) {
        col <- paste0("characteristic_", tag)
        row[[col]] <- collapse_unique(
          vals[tags == tag]
        )
      }
    }

    row
  })

  out <- do.call(rbind_fill, rows)
  out$geo_source_file <- xml_file
  out
}

load_geo_metadata <- function(datasets) {
  all <- list()

  for (gse in datasets) {
    log_msg("GEO metadata: ", gse)
    files <- download_geo_matrices(gse)

    if (!length(files)) {
      mini <- download_and_parse_miniml(gse)

      if (nrow(mini)) {
        all[[gse]] <- mini
        next
      }

      add_unresolved(
        "WARN", gse, "",
        "GEO metadata",
        "No series matrix or parseable MINiML metadata",
        "",
        paste0(
          "Install xml2 or provide primary ",
          "sample mapping/supplement"
        )
      )
      next
    }

    pieces <- lapply(files, function(f) {
      tryCatch(
        parse_geo_matrix_header(f, gse),
        error = function(e) data.frame()
      )
    })
    pieces <- pieces[
      vapply(pieces, nrow, integer(1)) > 0
    ]

    if (!length(pieces)) {
      mini <- download_and_parse_miniml(gse)

      if (nrow(mini)) {
        all[[gse]] <- mini
        next
      }

      add_unresolved(
        "WARN", gse, "",
        "GEO metadata",
        paste0(
          "Series matrix parsing failed and ",
          "MINiML was unavailable"
        ),
        paste(files, collapse = ";"),
        "Install xml2 or provide local sample mapping"
      )
      next
    }

    z <- do.call(rbind_fill, pieces)
    z$geo_source_file <- rep(
      paste(files, collapse = ";"),
      nrow(z)
    )
    all[[gse]] <- z
  }

  if (length(all)) {
    do.call(rbind_fill, all)
  } else {
    data.frame()
  }
}

standard_timepoint <- function(x) {
  y <- tolower(trim(x))
  if (is_missing(y)) return("Unknown")

  if (grepl(
    "progress|eot|end[ -]?of[ -]?treat|relapse",
    y
  )) {
    return("Progression/EOT")
  }

  if (grepl(
    "on[ -]?treat|cycle[ _-]?[0-9]+|^c[0-9]+$|during",
    y
  )) {
    return("On-treatment")
  }

  if (grepl(
    "post[ -]?treat|after[ -]?treat",
    y
  )) {
    return("Post-treatment")
  }

  if (grepl(
    paste0(
      "baseline|pre[ -]?treat|pretreat|",
      "treatment[ -]?naive|^pre$|^bl$"
    ),
    y
  )) {
    return("Baseline")
  }

  "Unknown"
}

standard_response <- function(x) {
  y <- tolower(trim(x))
  if (is_missing(y)) return("Unknown")

  if (grepl("complete response|^cr$", y)) {
    return("CR")
  }
  if (grepl("partial response|^pr$", y)) {
    return("PR")
  }
  if (grepl("stable disease|^sd$", y)) {
    return("SD")
  }
  if (grepl("progressive disease|^pd$", y)) {
    return("PD")
  }
  if (grepl(
    paste0(
      "non[- ]?responder|non[- ]?response|",
      "no response|lack of response"
    ),
    y
  )) {
    return("Non-responder")
  }
  if (grepl("responder|response", y)) {
    return("Responder")
  }

  "Unknown"
}

standard_tumor_normal <- function(x, tissue) {
  y <- tolower(paste(trim(x), trim(tissue)))

  if (grepl(
    "adjacent|normal|nat|non[- ]?tumou?r",
    y
  )) {
    return("Normal/NAT")
  }

  if (grepl(
    "tumou?r|ccrcc|carcinoma|metasta",
    y
  )) {
    return("Tumor")
  }

  "Unknown"
}

standard_primary_met <- function(x, site) {
  y <- tolower(paste(trim(x), trim(site)))

  if (grepl("metasta|secondary", y)) {
    return("Metastatic")
  }
  if (grepl("primary", y)) {
    return("Primary")
  }

  "Unknown"
}

standard_ici <- function(x, treatment, regimen) {
  y <- tolower(paste(
    trim(x),
    trim(treatment),
    trim(regimen)
  ))

  if (grepl(
    paste0(
      "nivolumab|ipilimumab|pembrolizumab|",
      "atezolizumab|avelumab|durvalumab|",
      "cemiplimab|anti[- ]?pd|anti[- ]?ctla|",
      "checkpoint|\\bici\\b|\\bicb\\b"
    ),
    y
  )) {
    return("Yes")
  }

  if (grepl(
    "no ici|ici[- ]?naive|without ici|no immunotherapy",
    y
  )) {
    return("No")
  }

  "Unknown"
}

explicit_resistance <- function(raw, definition) {
  y <- tolower(paste(
    trim(raw),
    trim(definition)
  ))

  if (grepl(
    "acquired|secondary resistance",
    y
  )) {
    return("Acquired-resistant")
  }

  if (grepl(
    "primary resistance|intrinsic resistance",
    y
  )) {
    return("Primary-resistant")
  }

  "Unknown"
}

log_msg("WORK_ROOT=", WORK_ROOT)

audit_file <- find_audit_file()
log_msg("AUDIT=", audit_file)

audit <- read.delim(
  audit_file,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

dataset_col <- pick_col(
  audit,
  c("dataset", "gse", "series"),
  FALSE
)
object_col <- pick_col(
  audit,
  c(
    "object_name", "object", "rds_object",
    "name", "file_name"
  ),
  FALSE
)
source_col <- pick_col(
  audit,
  c(
    "source_file", "rds_path", "file_path",
    "path", "rds_file", "file"
  ),
  FALSE
)
cells_col <- pick_col(
  audit,
  c("cells", "n_cells", "ncells", "cell_count"),
  FALSE
)

if (is.na(object_col) && is.na(source_col)) {
  stop(
    "Audit needs an object-name ",
    "or source-file column"
  )
}

rds_index <- list.files(
  WORK_ROOT,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
rds_index <- normalizePath(
  rds_index,
  winslash = "/",
  mustWork = TRUE
)
rds_index <- rds_index[
  !grepl(
    "phase4e_clinical_metadata_harmonization",
    rds_index,
    fixed = TRUE
  )
]

audit_std <- data.frame(
  dataset = if (!is.na(dataset_col)) {
    trim(audit[[dataset_col]])
  } else {
    ""
  },
  object_name = if (!is.na(object_col)) {
    trim(audit[[object_col]])
  } else {
    tools::file_path_sans_ext(
      basename(trim(audit[[source_col]]))
    )
  },
  source_raw = if (!is.na(source_col)) {
    trim(audit[[source_col]])
  } else {
    ""
  },
  audit_cells = if (!is.na(cells_col)) {
    suppressWarnings(
      as.numeric(audit[[cells_col]])
    )
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

if (any(is_missing(audit_std$dataset))) {
  src_text <- paste(
    audit_std$object_name,
    audit_std$source_raw
  )

  gse <- regmatches(
    src_text,
    regexpr(
      "GSE[0-9]+",
      src_text,
      ignore.case = TRUE
    )
  )

  gse[gse == ""] <- "Unknown"

  missing_dataset <- is_missing(
    audit_std$dataset
  )
  audit_std$dataset[missing_dataset] <-
    toupper(gse[missing_dataset])
}

audit_std$dataset <- toupper(
  audit_std$dataset
)

audit_std$source_file <- mapply(
  resolve_source_file,
  audit_std$source_raw,
  audit_std$object_name,
  MoreArgs = list(
    audit_file = audit_file,
    rds_index = rds_index
  ),
  USE.NAMES = FALSE
)

geo <- load_geo_metadata(
  sort(unique(
    audit_std$dataset[
      audit_std$dataset != "Unknown"
    ]
  ))
)

if (nrow(geo)) {
  write_tsv(
    geo,
    file.path(
      TABLE_DIR,
      "geo_sample_metadata_cache.tsv"
    )
  )
}

rows <- vector("list", nrow(audit_std))

for (i in seq_len(nrow(audit_std))) {
  a <- audit_std[i, ]

  log_msg(sprintf(
    "OBJECT %d/%d: %s",
    i,
    nrow(audit_std),
    a$object_name
  ))

  obj_meta <- NULL
  actual_cells <- NA_real_
  read_status <- "PASS"

  if (
    is.na(a$source_file) ||
    !file.exists(a$source_file)
  ) {
    read_status <- "SOURCE_NOT_FOUND"

    add_unresolved(
      "ERROR",
      a$dataset,
      a$object_name,
      "source_file",
      "RDS could not be resolved uniquely",
      a$source_raw,
      "Provide exact RDS path"
    )
  } else {
    obj <- tryCatch(
      readRDS(a$source_file),
      error = function(e) e
    )

    if (inherits(obj, "error")) {
      read_status <- "READ_ERROR"

      add_unresolved(
        "ERROR",
        a$dataset,
        a$object_name,
        "RDS",
        "readRDS failed",
        conditionMessage(obj)
      )
    } else {
      obj_meta <- extract_meta(obj)

      actual_cells <- if (!is.null(obj_meta)) {
        nrow(obj_meta)
      } else {
        tryCatch(
          ncol(obj),
          error = function(e) NA_real_
        )
      }

      rm(obj)
      invisible(gc(FALSE))
    }
  }

  rf <- extract_fields(obj_meta)

  gsm_candidates <- unique(c(
    regmatches(
      a$object_name,
      gregexpr(
        "GSM[0-9]+",
        a$object_name,
        ignore.case = TRUE
      )
    )[[1]],
    regmatches(
      a$source_raw,
      gregexpr(
        "GSM[0-9]+",
        a$source_raw,
        ignore.case = TRUE
      )
    )[[1]],
    if (!is_missing(rf$gsm$value)) {
      rf$gsm$value
    } else {
      character()
    }
  ))

  gsm_candidates <- toupper(
    gsm_candidates[
      grepl(
        "^GSM[0-9]+$",
        gsm_candidates,
        ignore.case = TRUE
      )
    ]
  )

  gsm <- if (length(gsm_candidates) == 1L) {
    gsm_candidates
  } else {
    "Unknown"
  }

  geo_row <- data.frame()
  geo_match_method <- "None"

  if (nrow(geo) && gsm != "Unknown") {
    hit <- which(
      geo$dataset == a$dataset &
        toupper(geo$gsm) == gsm
    )

    if (length(hit) == 1L) {
      geo_row <- geo[hit, , drop = FALSE]
      geo_match_method <- "Exact GSM"
    }
  }

  if (!nrow(geo_row) && nrow(geo)) {
    title_cols <- names(geo)[
      norm_name(names(geo)) %in%
        c("title", "sample_title")
    ]

    sample_try <- unique(c(
      a$object_name,
      rf$sample_id$value
    ))
    sample_try <- norm_token(
      sample_try[!is_missing(sample_try)]
    )

    hits <- integer()

    if (
      length(title_cols) &&
      length(sample_try)
    ) {
      for (tc in title_cols) {
        hits <- union(
          hits,
          which(
            geo$dataset == a$dataset &
              norm_token(geo[[tc]]) %in%
                sample_try
          )
        )
      }
    }

    if (length(hits) == 1L) {
      geo_row <- geo[hits, , drop = FALSE]
      geo_match_method <- "Exact normalized title"
      gsm <- first_nonmissing(geo_row$gsm)
    }
  }

  if (!nrow(geo_row)) {
    add_unresolved(
      "WARN",
      a$dataset,
      a$object_name,
      "GEO mapping",
      "No unique exact GEO sample match",
      gsm,
      paste0(
        "Resolve with primary sample mapping; ",
        "no fuzzy match was used"
      )
    )
  }

  gf <- extract_fields(
    if (nrow(geo_row)) geo_row else NULL
  )

  rec <- lapply(
    names(FIELD_PATTERNS),
    function(f) {
      reconcile(
        rf[[f]],
        gf[[f]],
        f,
        a$dataset,
        a$object_name
      )
    }
  )
  names(rec) <- names(FIELD_PATTERNS)

  val <- function(f) rec[[f]]$value
  src <- function(f) rec[[f]]$source
  ev <- function(f) rec[[f]]$evidence

  sample_id <- if (
    !is_missing(val("sample_id"))
  ) {
    val("sample_id")
  } else if (!is_missing(val("title"))) {
    val("title")
  } else {
    a$object_name
  }

  patient_id <- val("patient_id")
  tissue <- val("tissue")
  treatment <- val("treatment")
  regimen <- val("regimen")
  sampling_timepoint <- val("sampling_timepoint")
  response_raw <- val("response_raw")
  resistance_raw <- val("resistance_status")
  resistance_definition <- val(
    "resistance_definition"
  )

  timepoint_standard <- standard_timepoint(
    sampling_timepoint
  )
  response_standard <- standard_response(
    response_raw
  )
  resistance_status <- explicit_resistance(
    resistance_raw,
    resistance_definition
  )

  patient_key <- if (!is_missing(patient_id)) {
    paste(
      a$dataset,
      patient_id,
      sep = "::"
    )
  } else {
    "Unknown"
  }

  metadata_source <- collapse_unique(c(
    src("patient_id"),
    src("sampling_timepoint"),
    src("response_raw"),
    if (nrow(geo_row)) {
      "GEO series matrix/MINiML"
    } else {
      ""
    }
  ))

  source_reference <- if (
    nrow(geo_row) &&
    "geo_source_file" %in% names(geo_row)
  ) {
    first_nonmissing(
      geo_row$geo_source_file
    )
  } else {
    "Unknown"
  }

  evidence_level <- if (any(
    c(
      ev("patient_id"),
      ev("sampling_timepoint"),
      ev("response_raw")
    ) == "Unresolved"
  )) {
    "Unresolved"
  } else if (ev("patient_id") == "Missing") {
    "Low"
  } else {
    "Medium"
  }

  include_main <-
    patient_key != "Unknown" &&
    timepoint_standard != "Unknown" &&
    response_standard != "Unknown" &&
    evidence_level != "Unresolved"

  include_validation <-
    patient_key != "Unknown" &&
    read_status == "PASS"

  exclusion_reason <- collapse_unique(
    c(
      if (patient_key == "Unknown") {
        "Missing patient_id"
      } else {
        ""
      },
      if (timepoint_standard == "Unknown") {
        "Unknown timepoint"
      } else {
        ""
      },
      if (response_standard == "Unknown") {
        "Unknown response"
      } else {
        ""
      },
      if (read_status != "PASS") {
        read_status
      } else {
        ""
      }
    ),
    default = "None"
  )

  rows[[i]] <- data.frame(
    dataset = a$dataset,
    object_name = a$object_name,
    source_file = ifelse(
      is.na(a$source_file),
      "",
      a$source_file
    ),
    gsm = gsm,
    sample_id = sample_id,
    patient_id = patient_id,
    patient_key = patient_key,
    patient_id_source = src("patient_id"),
    patient_id_confidence = ev("patient_id"),
    tissue = tissue,
    tumor_normal = standard_tumor_normal(
      val("tumor_normal"),
      tissue
    ),
    primary_metastatic = standard_primary_met(
      val("primary_metastatic"),
      val("anatomical_site")
    ),
    anatomical_site = val("anatomical_site"),
    treatment = treatment,
    treatment_class = val("treatment_class"),
    ici_exposure = standard_ici(
      val("ici_exposure"),
      treatment,
      regimen
    ),
    regimen = regimen,
    sampling_timepoint = sampling_timepoint,
    timepoint_standard = timepoint_standard,
    response_raw = response_raw,
    response_standard = response_standard,
    resistance_status = resistance_status,
    resistance_definition = resistance_definition,
    cells = actual_cells,
    audit_cells = a$audit_cells,
    metadata_source = metadata_source,
    source_reference = source_reference,
    evidence_level = evidence_level,
    geo_match_method = geo_match_method,
    rds_read_status = read_status,
    include_main_analysis = include_main,
    include_validation = include_validation,
    exclusion_reason = exclusion_reason,
    notes = "",
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, rows)

override_candidates <- c(
  file.path(
    WORK_ROOT,
    "phase4e_manual_clinical_overrides.tsv"
  ),
  file.path(
    WORK_ROOT,
    "inputs",
    "phase4e_manual_clinical_overrides.tsv"
  )
)
override_file <- override_candidates[
  file.exists(override_candidates)
][1]

if (
  length(override_file) &&
  !is.na(override_file)
) {
  ov <- read.delim(
    override_file,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )

  needed <- c(
    "dataset", "object_name",
    "field", "value",
    "metadata_source", "source_reference"
  )

  if (!all(needed %in% names(ov))) {
    add_unresolved(
      "ERROR", "", "",
      "manual override",
      "Override file lacks required columns",
      paste(
        setdiff(needed, names(ov)),
        collapse = ","
      )
    )
  } else {
    allowed <- setdiff(
      names(manifest),
      c(
        "dataset", "object_name",
        "source_file", "cells",
        "audit_cells", "rds_read_status"
      )
    )

    for (j in seq_len(nrow(ov))) {
      idx <- which(
        manifest$dataset == ov$dataset[j] &
          manifest$object_name ==
            ov$object_name[j]
      )

      valid_source <-
        !is_missing(ov$metadata_source[j]) ||
        !is_missing(ov$source_reference[j])

      if (
        length(idx) != 1L ||
        !ov$field[j] %in% allowed ||
        !valid_source
      ) {
        add_unresolved(
          "ERROR",
          ov$dataset[j],
          ov$object_name[j],
          ov$field[j],
          paste0(
            "Invalid/unmatched manual override ",
            "or missing provenance"
          ),
          ov$value[j]
        )
        next
      }

      manifest[idx, ov$field[j]] <- ov$value[j]

      manifest$metadata_source[idx] <-
        collapse_unique(c(
          manifest$metadata_source[idx],
          paste0(
            "Manual:",
            ov$metadata_source[j]
          )
        ))

      manifest$source_reference[idx] <-
        collapse_unique(c(
          manifest$source_reference[idx],
          ov$source_reference[j]
        ))

      manifest$notes[idx] <-
        collapse_unique(c(
          manifest$notes[idx],
          if ("notes" %in% names(ov)) {
            ov$notes[j]
          } else {
            "manual override"
          }
        ))
    }
  }
}

for (i in seq_len(nrow(manifest))) {
  manifest$patient_key[i] <-
    if (!is_missing(manifest$patient_id[i])) {
      paste(
        manifest$dataset[i],
        manifest$patient_id[i],
        sep = "::"
      )
    } else {
      "Unknown"
    }

  tp <- standard_timepoint(
    manifest$sampling_timepoint[i]
  )
  if (
    tp != "Unknown" ||
    is_missing(manifest$timepoint_standard[i])
  ) {
    manifest$timepoint_standard[i] <- tp
  }

  rs <- standard_response(
    manifest$response_raw[i]
  )
  if (
    rs != "Unknown" ||
    is_missing(manifest$response_standard[i])
  ) {
    manifest$response_standard[i] <- rs
  }

  rr <- explicit_resistance(
    manifest$resistance_status[i],
    manifest$resistance_definition[i]
  )
  if (rr != "Unknown") {
    manifest$resistance_status[i] <- rr
  }

  manifest$ici_exposure[i] <- standard_ici(
    manifest$ici_exposure[i],
    manifest$treatment[i],
    manifest$regimen[i]
  )

  manifest$include_main_analysis[i] <-
    manifest$patient_key[i] != "Unknown" &&
    manifest$timepoint_standard[i] != "Unknown" &&
    manifest$response_standard[i] != "Unknown" &&
    manifest$evidence_level[i] != "Unresolved"

  manifest$include_validation[i] <-
    manifest$patient_key[i] != "Unknown" &&
    manifest$rds_read_status[i] == "PASS"

  manifest$exclusion_reason[i] <-
    collapse_unique(
      c(
        if (
          manifest$patient_key[i] == "Unknown"
        ) {
          "Missing patient_id"
        } else {
          ""
        },
        if (
          manifest$timepoint_standard[i] ==
            "Unknown"
        ) {
          "Unknown timepoint"
        } else {
          ""
        },
        if (
          manifest$response_standard[i] ==
            "Unknown"
        ) {
          "Unknown response"
        } else {
          ""
        },
        if (
          manifest$rds_read_status[i] != "PASS"
        ) {
          manifest$rds_read_status[i]
        } else {
          ""
        }
      ),
      default = "None"
    )
}

for (i in seq_len(nrow(manifest))) {
  for (
    f in c(
      "patient_id",
      "sampling_timepoint",
      "response_raw"
    )
  ) {
    val_f <- manifest[[f]][i]

    if (is_missing(val_f)) {
      add_unresolved(
        "WARN",
        manifest$dataset[i],
        manifest$object_name[i],
        f,
        "Required clinical field remains unresolved",
        "Unknown",
        paste0(
          "Resolve from GEO/paper/supplement ",
          "with source reference"
        )
      )
    }
  }

  if (
    !is.na(manifest$audit_cells[i]) &&
    !is.na(manifest$cells[i]) &&
    manifest$audit_cells[i] != manifest$cells[i]
  ) {
    add_unresolved(
      "ERROR",
      manifest$dataset[i],
      manifest$object_name[i],
      "cells",
      "Audit and RDS cell counts differ",
      paste0(
        "audit=", manifest$audit_cells[i],
        "; rds=", manifest$cells[i]
      )
    )
  }
}

known <- manifest[
  manifest$patient_key != "Unknown",
  ,
  drop = FALSE
]

patient_rows <- list()

if (nrow(known)) {
  keys <- unique(known$patient_key)

  for (k in keys) {
    z <- known[
      known$patient_key == k,
      ,
      drop = FALSE
    ]

    has_bl <- any(
      z$timepoint_standard == "Baseline"
    )
    has_on <- any(
      z$timepoint_standard == "On-treatment"
    )
    has_prog <- any(
      z$timepoint_standard == "Progression/EOT"
    )

    resp_bl <- unique(
      z$response_standard[
        z$timepoint_standard == "Baseline"
      ]
    )
    resp_bl <- resp_bl[resp_bl != "Unknown"]

    has_benefit_earlier <- any(
      z$response_standard %in%
        c("CR", "PR", "SD", "Responder") &
        z$timepoint_standard %in%
          c(
            "Baseline",
            "On-treatment",
            "Post-treatment"
          )
    )

    has_later_pd <- any(
      z$response_standard %in%
        c("PD", "Non-responder") &
        z$timepoint_standard ==
          "Progression/EOT"
    )

    explicit_acq <- any(
      z$resistance_status ==
        "Acquired-resistant"
    )

    acquired <-
      explicit_acq ||
      (has_benefit_earlier && has_later_pd)

    resistance <- if (acquired) {
      "Acquired-resistant"
    } else if (any(
      z$resistance_status ==
        "Primary-resistant"
    )) {
      "Primary-resistant"
    } else if (any(
      z$response_standard %in%
        c("CR", "PR", "SD", "Responder")
    )) {
      "Benefit"
    } else {
      "Unknown"
    }

    patient_rows[[k]] <- data.frame(
      dataset = z$dataset[1],
      patient_id = z$patient_id[1],
      patient_key = k,
      n_objects = nrow(z),
      n_cells = sum(z$cells, na.rm = TRUE),
      sample_ids = collapse_unique(
        z$sample_id
      ),
      timepoints = collapse_unique(
        z$timepoint_standard
      ),
      responses = collapse_unique(
        z$response_standard
      ),
      baseline_response = if (
        length(resp_bl) == 1L
      ) {
        resp_bl
      } else {
        "Unknown"
      },
      treatments = collapse_unique(
        z$treatment
      ),
      ici_exposure = collapse_unique(
        z$ici_exposure
      ),
      baseline_available = has_bl,
      on_treatment_available = has_on,
      progression_eot_available = has_prog,
      baseline_response_eligible =
        has_bl && length(resp_bl) == 1L,
      paired_change_eligible =
        has_bl && has_on,
      progression_evolution_eligible =
        has_prog && (has_bl || has_on),
      acquired_resistance_eligible =
        acquired,
      resistance_status = resistance,
      resistance_definition = if (
        explicit_acq
      ) {
        "Explicit source annotation"
      } else if (acquired) {
        paste0(
          "Derived within dataset+patient: ",
          "earlier benefit/stable control ",
          "followed by Progression/EOT PD"
        )
      } else {
        "Not established"
      },
      stringsAsFactors = FALSE
    )
  }
}

patient_manifest <- if (length(patient_rows)) {
  do.call(rbind, patient_rows)
} else {
  data.frame(
    dataset = character(),
    patient_id = character(),
    patient_key = character(),
    n_objects = integer(),
    n_cells = numeric(),
    sample_ids = character(),
    timepoints = character(),
    responses = character(),
    baseline_response = character(),
    treatments = character(),
    ici_exposure = character(),
    baseline_available = logical(),
    on_treatment_available = logical(),
    progression_eot_available = logical(),
    baseline_response_eligible = logical(),
    paired_change_eligible = logical(),
    progression_evolution_eligible = logical(),
    acquired_resistance_eligible = logical(),
    resistance_status = character(),
    resistance_definition = character()
  )
}

ds_rows <- lapply(
  sort(unique(manifest$dataset)),
  function(d) {
    z <- manifest[
      manifest$dataset == d,
      ,
      drop = FALSE
    ]
    p <- patient_manifest[
      patient_manifest$dataset == d,
      ,
      drop = FALSE
    ]

    data.frame(
      dataset = d,
      objects = nrow(z),
      cells = sum(z$cells, na.rm = TRUE),
      patients_resolved = nrow(p),
      patient_missing_rate = safe_rate(
        sum(z$patient_key == "Unknown"),
        nrow(z)
      ),
      timepoint_missing_rate = safe_rate(
        sum(
          z$timepoint_standard == "Unknown"
        ),
        nrow(z)
      ),
      response_missing_rate = safe_rate(
        sum(
          z$response_standard == "Unknown"
        ),
        nrow(z)
      ),
      baseline_response_patients = if (
        nrow(p)
      ) {
        sum(p$baseline_response_eligible)
      } else {
        0
      },
      paired_change_patients = if (
        nrow(p)
      ) {
        sum(p$paired_change_eligible)
      } else {
        0
      },
      progression_evolution_patients = if (
        nrow(p)
      ) {
        sum(
          p$progression_evolution_eligible
        )
      } else {
        0
      },
      strict_acquired_resistance_patients =
        if (nrow(p)) {
          sum(
            p$acquired_resistance_eligible
          )
        } else {
          0
        },
      main_analysis_objects = sum(
        z$include_main_analysis
      ),
      validation_objects = sum(
        z$include_validation
      ),
      stringsAsFactors = FALSE
    )
  }
)

dataset_summary <- do.call(rbind, ds_rows)

inventory_rows <- list()

for (d in sort(unique(manifest$dataset))) {
  p <- patient_manifest[
    patient_manifest$dataset == d,
    ,
    drop = FALSE
  ]

  bl <- p[
    p$baseline_response_eligible,
    ,
    drop = FALSE
  ]

  responders <- if (nrow(bl)) {
    sum(
      bl$baseline_response %in%
        c("CR", "PR", "Responder")
    )
  } else {
    0
  }

  nonresponders <- if (nrow(bl)) {
    sum(
      bl$baseline_response %in%
        c("SD", "PD", "Non-responder")
    )
  } else {
    0
  }

  specs <- list(
    c(
      "Baseline Responder vs Non-responder",
      nrow(bl),
      responders,
      nonresponders,
      if (
        responders >= 2 &&
        nonresponders >= 2
      ) {
        "Analysis-ready"
      } else if (nrow(bl)) {
        paste0(
          "Candidate; group sizes ",
          "insufficient/unresolved"
        )
      } else {
        "Not available"
      }
    ),
    c(
      "Baseline vs On-treatment paired",
      sum(p$paired_change_eligible),
      NA,
      NA,
      if (
        sum(p$paired_change_eligible) >= 2
      ) {
        "Analysis-ready"
      } else if (
        sum(p$paired_change_eligible) == 1
      ) {
        "Single pair only"
      } else {
        "Not available"
      }
    ),
    c(
      "Progression/EOT evolution",
      sum(
        p$progression_evolution_eligible
      ),
      NA,
      NA,
      if (
        sum(
          p$progression_evolution_eligible
        ) >= 2
      ) {
        "Analysis-ready"
      } else if (
        sum(
          p$progression_evolution_eligible
        ) == 1
      ) {
        "Single patient only"
      } else {
        "Not available"
      }
    ),
    c(
      "Strict acquired resistance",
      sum(
        p$acquired_resistance_eligible
      ),
      NA,
      NA,
      if (
        sum(
          p$acquired_resistance_eligible
        ) >= 2
      ) {
        "Analysis-ready"
      } else if (
        sum(
          p$acquired_resistance_eligible
        ) == 1
      ) {
        "Single patient only"
      } else {
        "Not available"
      }
    )
  )

  for (s in specs) {
    inventory_rows[[
      length(inventory_rows) + 1L
    ]] <- data.frame(
      dataset = d,
      analysis = s[1],
      eligible_patients =
        as.numeric(s[2]),
      responders = as.numeric(s[3]),
      nonresponders = as.numeric(s[4]),
      readiness = s[5],
      statistical_unit =
        "dataset + patient_id",
      notes = paste0(
        "Patient counts are not cell counts; ",
        "Unknown fields were not imputed"
      ),
      stringsAsFactors = FALSE
    )
  }
}

inventory <- do.call(
  rbind,
  inventory_rows
)

duplicate_objects <- sum(
  duplicated(manifest$object_name)
)

source_missing <- sum(
  manifest$source_file == "" |
    !file.exists(manifest$source_file)
)

read_fail <- sum(
  manifest$rds_read_status != "PASS"
)

cell_mismatch <- sum(
  !is.na(manifest$audit_cells) &
    !is.na(manifest$cells) &
    manifest$audit_cells != manifest$cells
)

total_cells <- sum(
  manifest$cells,
  na.rm = TRUE
)

reused_patient_ids <- if (nrow(known)) {
  patient_dataset_count <- aggregate(
    dataset ~ patient_id,
    known,
    function(x) length(unique(x))
  )
  sum(patient_dataset_count$dataset > 1)
} else {
  0
}

resistance_conflict <- if (
  nrow(patient_manifest)
) {
  sum(
    patient_manifest$resistance_status ==
      "Acquired-resistant" &
      !patient_manifest$progression_eot_available
  )
} else {
  0
}

manual_missing_source <- sum(
  grepl(
    "Manual:",
    manifest$metadata_source
  ) &
    is_missing(manifest$source_reference)
)

main_bad <- sum(
  manifest$include_main_analysis &
    (
      manifest$patient_key == "Unknown" |
        manifest$timepoint_standard ==
          "Unknown" |
        manifest$response_standard ==
          "Unknown"
    )
)

add_qc(
  "Manifest object count",
  EXPECTED_OBJECTS,
  nrow(manifest),
  if (
    nrow(manifest) == EXPECTED_OBJECTS
  ) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "Manifest total RDS cells",
  EXPECTED_CELLS,
  total_cells,
  if (total_cells == EXPECTED_CELLS) {
    "PASS"
  } else {
    "FAIL"
  },
  "Uses actual RDS/meta.data counts"
)

add_qc(
  "Expected datasets present",
  paste(EXPECTED_DATASETS, collapse = ","),
  paste(
    sort(unique(manifest$dataset)),
    collapse = ","
  ),
  if (setequal(
    EXPECTED_DATASETS,
    unique(manifest$dataset)
  )) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "object_name uniqueness",
  0,
  duplicate_objects,
  if (duplicate_objects == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "source_file exists",
  0,
  source_missing,
  if (source_missing == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "RDS read failures",
  0,
  read_fail,
  if (read_fail == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "Audit vs RDS cell mismatch",
  0,
  cell_mismatch,
  if (cell_mismatch == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "Patient missing rate",
  "Informational",
  sprintf(
    "%.3f",
    mean(manifest$patient_key == "Unknown")
  ),
  "INFO"
)

add_qc(
  "Timepoint missing rate",
  "Informational",
  sprintf(
    "%.3f",
    mean(
      manifest$timepoint_standard ==
        "Unknown"
    )
  ),
  "INFO"
)

add_qc(
  "Response missing rate",
  "Informational",
  sprintf(
    "%.3f",
    mean(
      manifest$response_standard ==
        "Unknown"
    )
  ),
  "INFO"
)

add_qc(
  "Raw patient IDs reused across datasets",
  "Composite key prevents merging",
  reused_patient_ids,
  "INFO",
  paste0(
    "All patient operations use ",
    "dataset::patient_id"
  )
)

add_qc(
  "Acquired resistance without Progression/EOT",
  0,
  resistance_conflict,
  if (resistance_conflict == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "Manual values without provenance",
  0,
  manual_missing_source,
  if (manual_missing_source == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

add_qc(
  "include_main_analysis lacking required evidence",
  0,
  main_bad,
  if (main_bad == 0) {
    "PASS"
  } else {
    "FAIL"
  }
)

write_tsv(
  manifest,
  file.path(
    TABLE_DIR,
    "object_clinical_manifest.tsv"
  )
)

write_tsv(
  patient_manifest,
  file.path(
    TABLE_DIR,
    "patient_longitudinal_manifest.tsv"
  )
)

write_tsv(
  dataset_summary,
  file.path(
    TABLE_DIR,
    "dataset_eligibility_summary.tsv"
  )
)

write_tsv(
  inventory,
  file.path(
    TABLE_DIR,
    "ici_resistance_analysis_inventory.tsv"
  )
)

write_tsv(
  unresolved,
  file.path(
    QC_DIR,
    "unresolved_metadata.tsv"
  )
)

write_tsv(
  qc,
  file.path(
    QC_DIR,
    "metadata_harmonization_qc.tsv"
  )
)

report <- c(
  "# Phase 4E clinical metadata harmonization report",
  "",
  paste0("Generated: ", stamp()),
  paste0("Audit input: `", audit_file, "`"),
  "",
  "## Scope",
  "",
  paste0(
    "This phase harmonizes object-, sample-, ",
    "and patient-level clinical metadata. ",
    "It does not merge Seurat objects, alter ",
    "RDS files, perform integration, ",
    "differential analysis, or claim external ",
    "biological validation."
  ),
  "",
  "## Core totals",
  "",
  paste0(
    "- Objects: ", nrow(manifest),
    " (expected ", EXPECTED_OBJECTS, ")"
  ),
  paste0(
    "- Actual RDS cells: ",
    format(
      total_cells,
      big.mark = ",",
      scientific = FALSE
    ),
    " (expected ",
    format(
      EXPECTED_CELLS,
      big.mark = ","
    ),
    ")"
  ),
  paste0(
    "- Resolved dataset-patient keys: ",
    nrow(patient_manifest)
  ),
  paste0(
    "- Unresolved records: ",
    nrow(unresolved)
  ),
  "",
  "## Dataset eligibility",
  "",
  paste(
    capture.output(print(
      dataset_summary,
      row.names = FALSE
    )),
    collapse = "\n"
  ),
  "",
  "## Interpretation rules",
  "",
  paste0(
    "- Patient matching uses ",
    "`dataset::patient_id`; raw patient IDs ",
    "are never merged across studies."
  ),
  paste0(
    "- PD alone is not acquired resistance. ",
    "Acquired resistance requires explicit ",
    "source annotation or within-dataset ",
    "longitudinal evidence of earlier benefit/",
    "stable control followed by ",
    "Progression/EOT PD."
  ),
  paste0(
    "- No fuzzy GEO-to-object matching was used. ",
    "Unmatched objects remain unresolved."
  ),
  paste0(
    "- Unknown clinical fields remain `Unknown`; ",
    "they are not inferred from GSE294109 ",
    "filename suffixes or other unverified ",
    "naming conventions."
  ),
  paste0(
    "- Cell and spot counts are not biological ",
    "replicates; all eligibility counts are ",
    "patient-level."
  ),
  "",
  "## Completion status",
  ""
)

core_fail <- any(qc$status == "FAIL")

report <- c(
  report,
  if (core_fail) {
    paste0(
      "Core QC failed. Phase 4E outputs were ",
      "written for diagnosis, but the phase ",
      "is not complete."
    )
  } else {
    paste0(
      "Core structural QC passed. Harmonization ",
      "completed with unresolved clinical fields ",
      "explicitly retained; this does not mean ",
      "that any downstream analysis is eligible."
    )
  }
)

writeLines(
  report,
  file.path(
    OUT_DIR,
    "PHASE4E_CLINICAL_METADATA_REPORT.md"
  ),
  useBytes = TRUE
)

status_file <- file.path(
  OUT_DIR,
  if (core_fail) {
    "PHASE4E_INCOMPLETE"
  } else {
    "PHASE4E_CLINICAL_METADATA_COMPLETE"
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
    paste0("cells=", total_cells)
  ),
  status_file
)

log_msg(
  "STATUS=",
  if (core_fail) "FAIL" else "PASS"
)

if (core_fail) {
  quit(save = "no", status = 2L)
}
