#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 02_validate_nci60_preprocessing_v1.R <phase101_nci60_dir> <phase102_output_dir>")
input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
results_dir <- file.path(output_dir, "results_v1")
qc_dir <- file.path(output_dir, "qc_v1")

all_path <- file.path(results_dir, "nci60_all_lineages_preprocessed_v1.tsv")
renal_path <- file.path(results_dir, "nci60_renal_preprocessed_v1.tsv")
rep_path <- file.path(results_dir, "mta_technical_replicates_long_v1.tsv")
map_path <- file.path(qc_dir, "cell_line_name_mapping_v1.tsv")
manifest_path <- file.path(output_dir, "FILE_MANIFEST_SHA256_v1.tsv")
required <- c(all_path, renal_path, rep_path, map_path, manifest_path)
if (!all(file.exists(required))) stop("Missing Phase 102 output")

all_dt <- fread(all_path)
renal <- fread(renal_path)
rep_dt <- fread(rep_path)
mapping <- fread(map_path)
manifest <- fread(manifest_path)

source_raw <- fread(file.path(input_dir, "WEB_DATA_METABOLON_ALL.TXT"), showProgress = FALSE)
source_mta <- source_raw[tolower(TITLE) == "5'-s-methyl-5'-thioadenosine"]
source_mta[, cell_key := gsub("[^A-Z0-9]", "", sub("/ATCC$", "", toupper(trimws(cellname))))]
source_agg <- source_mta[, .(source_n = uniqueN(SAMPLENAME), source_mean = mean(as.numeric(Value)), source_sd = sd(as.numeric(Value))), by = cell_key]
joined <- merge(all_dt, source_agg, by = "cell_key", all.x = TRUE)

rna <- as.data.table(read_excel(file.path(input_dir, "RNA__RNA_seq_composite_expression.xls"), sheet = "Results", skip = 10L, .name_repair = "minimal"))
rna <- rna[get("Gene name d") %chin% c("TYMP", "MTAP")]
rna_cols <- names(rna)[-(1:6)]
independent_rna <- melt(rna, id.vars = "Gene name d", measure.vars = rna_cols, variable.name = "header", value.name = "value")
independent_rna[, cell_key := gsub("[^A-Z0-9]", "", sub("/ATCC$", "", toupper(sub("^[^:]+:", "", header))))]
independent_wide <- dcast(independent_rna, cell_key ~ `Gene name d`, value.var = "value")
rna_join <- merge(all_dt[, .(cell_key, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)], independent_wide, by = "cell_key", all.x = TRUE)

hash_rows <- manifest[, {
  p <- file.path(output_dir, relative_path)
  .(exists = file.exists(p), observed_sha256 = if (file.exists(p)) digest(p, algo = "sha256", file = TRUE) else NA_character_)
}, by = .(relative_path, sha256)]

expected_renal <- sort(c("7860", "A498", "ACHN", "CAKI1", "RXF393", "SN12C", "TK10", "UO31"))
checks <- data.table(
  check_id = c(
    "all_table_58_rows", "renal_table_8_rows", "renal_identity_exact", "all_values_complete",
    "replicate_source_row_count", "replicate_aggregation_n", "replicate_aggregation_mean",
    "replicate_aggregation_sd", "independent_rna_tymp", "independent_rna_mtap",
    "mapping_58_matched", "manifest_files_exist", "manifest_hashes_match",
    "no_inferential_columns", "analysis_flag_false"
  ),
  pass = c(
    nrow(all_dt) == 58L,
    nrow(renal) == 8L,
    identical(sort(renal$cell_key), expected_renal),
    all(complete.cases(all_dt[, .(mta_mean_raw, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)])),
    nrow(rep_dt) == nrow(source_mta),
    all(joined$mta_n_technical_repeats == joined$source_n),
    max(abs(joined$mta_mean_raw - joined$source_mean)) < 1e-12,
    max(abs(joined$mta_sd_raw - joined$source_sd)) < 1e-12,
    max(abs(rna_join$tymp_log2_fpkm_plus1 - as.numeric(rna_join$TYMP))) < 1e-12,
    max(abs(rna_join$mtap_log2_fpkm_plus1 - as.numeric(rna_join$MTAP))) < 1e-12,
    mapping[match_status == "MATCHED", .N] == 58L,
    all(hash_rows$exists),
    all(toupper(hash_rows$sha256) == toupper(hash_rows$observed_sha256)),
    !any(grepl("(^|_)(p|q|fdr|rho|correlation|coefficient|ci)(_|$)", names(all_dt), ignore.case = TRUE)),
    isFALSE(fromJSON(file.path(qc_dir, "preprocessing_status_v1.json"))$analysis_performed)
  )
)
checks[, detail := c(
  paste0("observed=", nrow(all_dt)), paste0("observed=", nrow(renal)), paste(renal$cell_key, collapse = ";"),
  paste0("complete=", sum(complete.cases(all_dt[, .(mta_mean_raw, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)]))),
  paste0("output=", nrow(rep_dt), ";source=", nrow(source_mta)),
  paste0("max_delta=", max(abs(joined$mta_n_technical_repeats - joined$source_n))),
  sprintf("max_delta=%.3g", max(abs(joined$mta_mean_raw - joined$source_mean))),
  sprintf("max_delta=%.3g", max(abs(joined$mta_sd_raw - joined$source_sd))),
  sprintf("max_delta=%.3g", max(abs(rna_join$tymp_log2_fpkm_plus1 - as.numeric(rna_join$TYMP)))),
  sprintf("max_delta=%.3g", max(abs(rna_join$mtap_log2_fpkm_plus1 - as.numeric(rna_join$MTAP)))),
  paste0("matched=", mapping[match_status == "MATCHED", .N]),
  paste0("present=", sum(hash_rows$exists), "/", nrow(hash_rows)),
  paste0("matching=", sum(toupper(hash_rows$sha256) == toupper(hash_rows$observed_sha256)), "/", nrow(hash_rows)),
  paste(names(all_dt), collapse = ";"),
  paste0("analysis_performed=", fromJSON(file.path(qc_dir, "preprocessing_status_v1.json"))$analysis_performed)
)]

fwrite(checks, file.path(qc_dir, "independent_validation_checks_v1.tsv"), sep = "\t", quote = FALSE)
validation <- list(
  status = if (all(checks$pass)) "PASS_INDEPENDENT_VALIDATION" else "FAIL_INDEPENDENT_VALIDATION",
  validated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks_passed = sum(checks$pass),
  checks_total = nrow(checks),
  failed_checks = checks[pass == FALSE, check_id]
)
write_json(validation, file.path(qc_dir, "independent_validation_status_v1.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(checks$pass)) stop("Independent validation failed: ", paste(validation$failed_checks, collapse = ", "))
cat("PASS_INDEPENDENT_VALIDATION\n")
cat("checks=", validation$checks_passed, "/", validation$checks_total, "\n", sep = "")
