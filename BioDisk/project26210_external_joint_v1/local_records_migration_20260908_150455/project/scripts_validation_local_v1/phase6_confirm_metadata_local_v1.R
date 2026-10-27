options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "\\", mustWork = TRUE) else normalizePath(getwd(), winslash = "\\", mustWork = TRUE)
out_dir <- file.path(root, "phase6_metadata_confirmation_local_v1")
table_dir <- file.path(out_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

write_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  if (!file.rename(part, path)) stop("Rename failed: ", path)
}

manifest_file <- file.path(root, "phase4e_clinical_metadata_harmonization_v1", "tables", "object_clinical_manifest.tsv")
geo_file <- file.path(root, "phase4e_clinical_metadata_harmonization_v1", "tables", "geo_sample_metadata_cache.tsv")
manifest <- fread(manifest_file, sep = "\t", encoding = "UTF-8")
geo <- fread(geo_file, sep = "\t", encoding = "UTF-8")
geo <- unique(geo, by = c("dataset", "gsm"))
x <- merge(manifest, geo, by = c("dataset", "gsm"), all.x = TRUE, suffixes = c("_old", "_geo"))
if (nrow(x) != 95L || anyNA(x$title)) stop("Expected 95 exact GSM-title matches")

first_line <- function(z) sub("[\r\n].*$", "", z)
match_one <- function(pattern, text, group = 1L) {
  m <- regexec(pattern, text, perl = TRUE)
  z <- regmatches(text, m)
  vapply(z, function(a) if (length(a) > group) a[[group + 1L]] else NA_character_, character(1))
}
clean_unknown <- function(z) {
  z <- trimws(as.character(z))
  z[is.na(z) | !nzchar(z) | toupper(z) %in% c("NA", "N/A", "UNKNOWN", "NOT REPORTED")] <- "Unknown"
  z
}

x[, `:=`(
  patient_id_confirmed = "Unknown",
  sample_id_confirmed = "Unknown",
  timepoint_confirmed = "Unknown",
  treatment_confirmed = "Unknown",
  response_confirmed = "Unknown",
  tissue_confirmed = "Unknown",
  sex_confirmed = "Unknown",
  age_confirmed = "Unknown",
  mapping_rule = "",
  mapping_confidence = "",
  mapping_status = "FAIL"
)]

# GSE178481: the transferred source RDS contains an explicit patient_id, and GEO titles independently reproduce it.
i <- x$dataset == "GSE178481"
geo_patient <- match_one("^RCC-([A-Z]+[0-9]+)-", x$title[i])
x[i, `:=`(
  patient_id_confirmed = patient_id,
  sample_id_confirmed = sample_id,
  tissue_confirmed = clean_unknown(tissue),
  treatment_confirmed = clean_unknown(treatment),
  mapping_rule = "explicit RDS patient_id concordant with official GEO title RCC-<patient>-<specimen>",
  mapping_confidence = "High",
  mapping_status = ifelse(patient_id == geo_patient, "PASS", "FAIL")
)]

# GSE202374: six official GEO sample titles contain a unique public subject code (RS########) and sample code.
i <- x$dataset == "GSE202374"
sample_code <- match_one("^Sample [0-9]+ \\(([A-Z]+[0-9]+) - RS[0-9]+\\)$", x$title[i])
subject_code <- match_one("^Sample [0-9]+ \\([A-Z]+[0-9]+ - (RS[0-9]+)\\)$", x$title[i])
tx <- clean_unknown(x$characteristic_treatment[i])
x[i, `:=`(
  patient_id_confirmed = subject_code,
  sample_id_confirmed = sample_code,
  timepoint_confirmed = ifelse(tx == "15 Gy radiation dose", "Surgery_4wk_post_SBRT", "Nephrectomy_only"),
  treatment_confirmed = tx,
  tissue_confirmed = clean_unknown(x$source_name_ch1[i]),
  sex_confirmed = clean_unknown(x$characteristic_gender[i]),
  mapping_rule = "official GEO title Sample n (<sample_code> - <RS subject code>); one scRNA tumor per unique subject",
  mapping_confidence = "Medium",
  mapping_status = ifelse(!is.na(subject_code) & !is.na(sample_code), "PASS", "FAIL")
)]

# GSE207493: official series design states 19 patients; paired scRNA/scATAC titles share RCC patient IDs.
i <- x$dataset == "GSE207493"
rcc <- match_one("^scRNA-seq_(RCC[0-9]+)$", x$title[i])
x[i, `:=`(
  patient_id_confirmed = rcc,
  sample_id_confirmed = rcc,
  tissue_confirmed = clean_unknown(x$characteristic_tissue[i]),
  mapping_rule = "official GEO scRNA-seq_<RCC patient ID>; series design confirms 19 patients with paired scRNA/scATAC",
  mapping_confidence = "High",
  mapping_status = ifelse(!is.na(rcc), "PASS", "FAIL")
)]

# GSE264586: official series design and titles explicitly define 12 patients and Baseline/C2/EOT biopsies.
i <- x$dataset == "GSE264586"
rcc <- match_one("^(RCC[0-9]+), (?:Baseline|C2|EOT)$", x$title[i])
tp_raw <- match_one("^RCC[0-9]+, (Baseline|C2|EOT)$", x$title[i])
tp <- fifelse(tp_raw == "C2", "Cycle2", fifelse(tp_raw == "EOT", "End_of_treatment", tp_raw))
x[i, `:=`(
  patient_id_confirmed = rcc,
  sample_id_confirmed = paste(rcc, tp, sep = "__"),
  timepoint_confirmed = tp,
  treatment_confirmed = "sitravatinib+nivolumab+ipilimumab",
  tissue_confirmed = clean_unknown(x$characteristic_tissue[i]),
  mapping_rule = "official GEO title <RCC patient ID>, <Baseline/C2/EOT>; series design confirms 12 patients and 19 biopsies",
  mapping_confidence = "High",
  mapping_status = ifelse(!is.na(rcc) & !is.na(tp), "PASS", "FAIL")
)]

# GSE294109: official titles and library descriptions repeat scXX across multiple regions from the same case.
i <- x$dataset == "GSE294109"
sc <- match_one("(?:^|, )(sc[0-9]+)(?:[.,]|$)", x$title[i])
lib <- match_one("^Library name: ([^\r\n]+)", x$description[i])
tx <- clean_unknown(x$characteristic_pre_surgery_treatment[i])
resp <- clean_unknown(x$characteristic_pathological_response[i])
x[i, `:=`(
  patient_id_confirmed = sc,
  sample_id_confirmed = lib,
  timepoint_confirmed = ifelse(tx == "Untreated", "Surgery_untreated", "Surgery_after_treatment"),
  treatment_confirmed = tx,
  response_confirmed = resp,
  tissue_confirmed = clean_unknown(x$characteristic_tissue[i]),
  sex_confirmed = clean_unknown(x$characteristic_sex[i]),
  age_confirmed = clean_unknown(x$characteristic_age_at_time_of_surgery[i]),
  mapping_rule = "official GEO title contains scXX case ID; library description supplies region-specific scXXT[_region] ID",
  mapping_confidence = "High",
  mapping_status = ifelse(!is.na(sc) & !is.na(lib), "PASS", "FAIL")
)]

crosswalk <- x[, .(
  dataset, object_name, gsm, geo_title = title,
  patient_id_confirmed, patient_key = paste(dataset, patient_id_confirmed, sep = "::"),
  sample_id_confirmed, timepoint_confirmed, treatment_confirmed, response_confirmed,
  tissue_confirmed, sex_confirmed, age_confirmed,
  mapping_rule, mapping_confidence, mapping_status,
  source_file, geo_source_file
)]
setorder(crosswalk, dataset, gsm)

if (any(crosswalk$mapping_status != "PASS") || any(crosswalk$patient_id_confirmed == "Unknown")) stop("One or more mappings failed")
if (anyDuplicated(crosswalk[, .(dataset, object_name)]) || anyDuplicated(crosswalk[, .(dataset, gsm)])) stop("Duplicate object/GSM mappings")

expected <- data.table(dataset = c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"),
                       expected_objects = c(26L, 6L, 19L, 19L, 25L),
                       expected_patients = c(11L, 6L, 19L, 12L, 15L))
summary <- crosswalk[, .(objects = .N, patients = uniqueN(patient_id_confirmed),
                         samples = uniqueN(sample_id_confirmed), mapping_pass = all(mapping_status == "PASS")), by = dataset]
summary <- merge(expected, summary, by = "dataset", all = TRUE)
summary[, `:=`(objects_match = objects == expected_objects, patients_match = patients == expected_patients)]

longitudinal <- crosswalk[dataset == "GSE264586", .N, by = timepoint_confirmed]
gse294_consistency <- crosswalk[dataset == "GSE294109", .(
  libraries = .N,
  treatment_values = uniqueN(treatment_confirmed),
  response_values = uniqueN(response_confirmed),
  sex_values = uniqueN(sex_confirmed),
  age_values = uniqueN(age_confirmed),
  consistent = uniqueN(treatment_confirmed) == 1L && uniqueN(response_confirmed) == 1L && uniqueN(sex_confirmed) == 1L && uniqueN(age_confirmed) == 1L
), by = patient_id_confirmed]

if (!all(summary$objects_match & summary$patients_match & summary$mapping_pass)) stop("Dataset summary mismatch")
if (!identical(longitudinal[timepoint_confirmed == "Baseline", N], 12L) ||
    !identical(longitudinal[timepoint_confirmed == "Cycle2", N], 4L) ||
    !identical(longitudinal[timepoint_confirmed == "End_of_treatment", N], 3L)) stop("GSE264586 timepoint mismatch")
if (!all(gse294_consistency$consistent)) stop("GSE294109 within-patient clinical inconsistency")

sources <- data.table(
  dataset = c("GSE178481", "GSE202374", "GSE202374", "GSE207493", "GSE264586", "GSE294109"),
  source_type = c("NCBI_GEO", "NCBI_GEO", "Primary_publication", "NCBI_GEO", "NCBI_GEO", "NCBI_GEO_and_primary_publication"),
  url = c(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178481",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202374",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10124322/",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207493",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264586",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294109"
  ),
  checked_at = "2026-08-10",
  local_evidence = normalizePath(geo_file, winslash = "\\", mustWork = TRUE)
)

write_atomic(crosswalk, file.path(table_dir, "confirmed_object_patient_crosswalk.tsv"))
write_atomic(summary, file.path(table_dir, "metadata_confirmation_summary.tsv"))
write_atomic(longitudinal, file.path(table_dir, "GSE264586_timepoint_reconciliation.tsv"))
write_atomic(gse294_consistency, file.path(table_dir, "GSE294109_patient_consistency.tsv"))
write_atomic(sources, file.path(table_dir, "metadata_source_evidence.tsv"))

status <- c(
  "status=PASS",
  "objects_expected=95",
  paste0("objects_confirmed=", nrow(crosswalk)),
  paste0("patients_confirmed=", uniqueN(crosswalk$patient_key)),
  "datasets_confirmed=5",
  "GSE264586_baseline=12",
  "GSE264586_cycle2=4",
  "GSE264586_end_of_treatment=3",
  "source_rds_modified=0",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(status, file.path(out_dir, "PHASE6_METADATA_CONFIRMATION_STATUS"), useBytes = TRUE)
cat("PHASE6_METADATA_CONFIRMATION=PASS\n")
