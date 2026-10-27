#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_preprocess_nci60_tymp_mta_v1.R <phase101_nci60_dir> <output_dir>")
}

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
results_dir <- file.path(output_dir, "results_v1")
qc_dir <- file.path(output_dir, "qc_v1")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

raw_mta_path <- file.path(input_dir, "WEB_DATA_METABOLON_ALL.TXT")
avg_mta_path <- file.path(input_dir, "WEB_DATA_METABOLON_extracted", "WEB_DATA_METABOLON.TXT")
rna_path <- file.path(input_dir, "RNA__RNA_seq_composite_expression.xls")
required_paths <- c(raw_mta_path, avg_mta_path, rna_path)
if (!all(file.exists(required_paths))) {
  stop("Missing required input(s): ", paste(required_paths[!file.exists(required_paths)], collapse = "; "))
}

normalize_cell <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("/ATCC$", "", x)
  gsub("[^A-Z0-9]", "", x)
}

exact_mta <- "5'-s-methyl-5'-thioadenosine"
generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

raw <- fread(raw_mta_path, na.strings = c("", "NA"), showProgress = FALSE)
required_raw_cols <- c("SAMPLENAME", "TITLE", "panelnbr", "cellnbr", "pname", "cellname", "Value")
stopifnot(all(required_raw_cols %in% names(raw)))
raw_mta <- raw[tolower(TITLE) == tolower(exact_mta)]
if (nrow(raw_mta) == 0L) stop("Exact MTA identity was not found in replicate-level input")
raw_mta[, Value := as.numeric(Value)]
if (anyNA(raw_mta$Value)) stop("Non-numeric or missing exact-MTA replicate value")
raw_mta[, cell_key := normalize_cell(cellname)]
setorder(raw_mta, pname, cellname, SAMPLENAME)

global_floor <- min(raw_mta$Value)
mta_agg <- raw_mta[, .(
  tissue = unique(pname),
  panel_number = unique(panelnbr),
  cell_number = unique(cellnbr),
  mta_n_technical_repeats = uniqueN(SAMPLENAME),
  mta_mean_raw = mean(Value),
  mta_median_raw = median(Value),
  mta_sd_raw = if (.N > 1L) sd(Value) else NA_real_,
  mta_min_raw = min(Value),
  mta_max_raw = max(Value),
  mta_floor_repeat_count = sum(Value == global_floor),
  mta_any_repeat_at_global_floor = any(Value == global_floor),
  mta_all_repeats_at_global_floor = all(Value == global_floor)
), by = .(cell_key, metabolomics_cell_line = cellname)]
if (mta_agg[, any(lengths(tissue) != 1L)]) stop("Ambiguous tissue assignment")

avg <- fread(avg_mta_path, na.strings = c("", "NA"), showProgress = FALSE)
required_avg_cols <- c("TITLE", "pname", "cellname", "VALUE", "STD DEV")
stopifnot(all(required_avg_cols %in% names(avg)))
avg_mta <- avg[tolower(TITLE) == tolower(exact_mta)]
avg_mta[, cell_key := normalize_cell(cellname)]
avg_mta <- avg_mta[, .(
  official_cell_line = cellname,
  mta_mean_official = as.numeric(VALUE),
  mta_sd_official = as.numeric(`STD DEV`)
), by = cell_key]
if (anyDuplicated(avg_mta$cell_key)) stop("Duplicate exact-MTA rows in official averaged table")

rna <- as.data.table(read_excel(rna_path, sheet = "Results", skip = 10L, .name_repair = "minimal"))
gene_col <- "Gene name d"
if (!gene_col %in% names(rna)) stop("CellMiner gene-name column not found")
rna_target <- rna[get(gene_col) %chin% c("TYMP", "MTAP")]
if (!identical(sort(rna_target[[gene_col]]), c("MTAP", "TYMP"))) {
  stop("Expected one unambiguous MTAP row and one unambiguous TYMP row")
}
metadata_cols <- names(rna)[1:6]
rna_cols <- setdiff(names(rna), metadata_cols)
rna_long <- melt(
  rna_target,
  id.vars = gene_col,
  measure.vars = rna_cols,
  variable.name = "rna_header",
  value.name = "rna_log2_fpkm_plus1"
)
rna_long[, rna_panel_code := sub(":.*$", "", rna_header)]
rna_long[, rna_cell_line := sub("^[^:]+:", "", rna_header)]
rna_long[, cell_key := normalize_cell(rna_cell_line)]
rna_long[, rna_log2_fpkm_plus1 := as.numeric(rna_log2_fpkm_plus1)]
if (anyDuplicated(rna_long[, .(cell_key, gene = get(gene_col))])) stop("Duplicate normalized RNA cell-line/gene key")

rna_wide <- dcast(
  rna_long[, .(cell_key, rna_cell_line, rna_panel_code, gene = get(gene_col), rna_log2_fpkm_plus1)],
  cell_key + rna_cell_line + rna_panel_code ~ gene,
  value.var = "rna_log2_fpkm_plus1"
)
setnames(rna_wide, c("MTAP", "TYMP"), c("mtap_log2_fpkm_plus1", "tymp_log2_fpkm_plus1"))

paired <- merge(mta_agg, avg_mta, by = "cell_key", all.x = TRUE, sort = FALSE)
paired <- merge(paired, rna_wide, by = "cell_key", all.x = TRUE, sort = FALSE)
paired[, mta_mean_abs_delta_official := abs(mta_mean_raw - mta_mean_official)]
paired[, mta_sd_abs_delta_official := abs(mta_sd_raw - mta_sd_official)]
paired[, is_renal := tolower(tissue) == "renal"]
paired[, analysis_level_all_nci60 := TRUE]
paired[, analysis_level_renal := is_renal]
setcolorder(paired, c(
  "cell_key", "metabolomics_cell_line", "rna_cell_line", "tissue", "is_renal",
  "analysis_level_all_nci60", "analysis_level_renal", "rna_panel_code",
  "tymp_log2_fpkm_plus1", "mtap_log2_fpkm_plus1",
  "mta_n_technical_repeats", "mta_mean_raw", "mta_median_raw", "mta_sd_raw",
  "mta_min_raw", "mta_max_raw", "mta_mean_official", "mta_sd_official",
  "mta_mean_abs_delta_official", "mta_sd_abs_delta_official",
  "mta_floor_repeat_count", "mta_any_repeat_at_global_floor", "mta_all_repeats_at_global_floor",
  "panel_number", "cell_number", "official_cell_line"
))
setorder(paired, tissue, metabolomics_cell_line)
renal <- paired[is_renal == TRUE]

raw_out <- raw_mta[, .(
  sample_name = SAMPLENAME,
  cell_key,
  cell_line = cellname,
  tissue = pname,
  panel_number = panelnbr,
  cell_number = cellnbr,
  metabolite = TITLE,
  mta_value = Value,
  at_global_floor = Value == global_floor
)]

rna_union <- merge(
  unique(rna_wide[, .(cell_key, rna_cell_line, rna_panel_code)]),
  unique(mta_agg[, .(cell_key, metabolomics_cell_line, tissue)]),
  by = "cell_key", all = TRUE
)
rna_union[, match_status := fifelse(
  !is.na(rna_cell_line) & !is.na(metabolomics_cell_line), "MATCHED",
  fifelse(is.na(rna_cell_line), "METABOLOMICS_ONLY", "RNA_ONLY")
)]
setorder(rna_union, match_status, tissue, metabolomics_cell_line, rna_cell_line)

expected_renal_keys <- sort(normalize_cell(c("786-0", "A498", "ACHN", "CAKI-1", "RXF-393", "SN12C", "TK-10", "UO-31")))
checks <- data.table(
  check_id = c(
    "required_inputs_present", "exact_mta_identity_unique", "replicate_mta_58_lines",
    "official_mta_58_lines", "replicate_counts_at_least_two", "all_metabolomics_lines_match_rna",
    "paired_values_complete", "official_mean_crosscheck", "official_sd_crosscheck",
    "renal_subset_eight", "renal_subset_identity", "renal_values_complete",
    "technical_repeats_not_rows_in_analysis", "no_duplicate_analysis_cell_key"
  ),
  pass = c(
    all(file.exists(required_paths)),
    uniqueN(raw_mta$TITLE) == 1L && tolower(unique(raw_mta$TITLE)) == tolower(exact_mta),
    uniqueN(raw_mta$cell_key) == 58L,
    uniqueN(avg_mta$cell_key) == 58L,
    min(mta_agg$mta_n_technical_repeats) >= 2L,
    nrow(paired) == 58L && all(!is.na(paired$rna_cell_line)),
    all(complete.cases(paired[, .(mta_mean_raw, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)])),
    max(paired$mta_mean_abs_delta_official, na.rm = TRUE) <= 0.001,
    max(paired$mta_sd_abs_delta_official, na.rm = TRUE) <= 0.001,
    nrow(renal) == 8L,
    identical(sort(renal$cell_key), expected_renal_keys),
    all(complete.cases(renal[, .(mta_mean_raw, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)])),
    nrow(paired) == uniqueN(raw_mta$cell_key),
    anyDuplicated(paired$cell_key) == 0L
  )
)
checks[, detail := c(
  paste(basename(required_paths), collapse = "; "),
  paste(unique(raw_mta$TITLE), collapse = "; "),
  paste0("observed=", uniqueN(raw_mta$cell_key)),
  paste0("observed=", uniqueN(avg_mta$cell_key)),
  paste0("range=", min(mta_agg$mta_n_technical_repeats), "-", max(mta_agg$mta_n_technical_repeats)),
  paste0("matched=", sum(!is.na(paired$rna_cell_line)), "/", nrow(paired)),
  paste0("complete=", sum(complete.cases(paired[, .(mta_mean_raw, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)])), "/", nrow(paired)),
  sprintf("max_abs_delta=%.12g", max(paired$mta_mean_abs_delta_official, na.rm = TRUE)),
  sprintf("max_abs_delta=%.12g", max(paired$mta_sd_abs_delta_official, na.rm = TRUE)),
  paste0("observed=", nrow(renal)),
  paste(renal$metabolomics_cell_line, collapse = "; "),
  paste0("complete=", sum(complete.cases(renal[, .(mta_mean_raw, tymp_log2_fpkm_plus1, mtap_log2_fpkm_plus1)])), "/", nrow(renal)),
  paste0("analysis_rows=", nrow(paired), "; replicate_rows=", nrow(raw_mta)),
  paste0("duplicates=", anyDuplicated(paired$cell_key))
)]

if (!all(checks$pass)) {
  fwrite(checks, file.path(qc_dir, "preprocessing_checks_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")
  stop("Preprocessing checks failed: ", paste(checks[pass == FALSE, check_id], collapse = ", "))
}

fwrite(raw_out, file.path(results_dir, "mta_technical_replicates_long_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(paired, file.path(results_dir, "nci60_all_lineages_preprocessed_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(renal, file.path(results_dir, "nci60_renal_preprocessed_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(rna_union, file.path(qc_dir, "cell_line_name_mapping_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(checks, file.path(qc_dir, "preprocessing_checks_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")

dictionary <- data.table(
  field = names(paired),
  definition = c(
    "Punctuation-insensitive cell-line matching key", "Cell-line name in metabolomics source",
    "Cell-line name in original CellMiner workbook", "NCI60 tissue lineage", "Renal lineage flag",
    "Eligibility flag for future all-lineage analysis", "Eligibility flag for future renal-only analysis",
    "CellMiner RNA panel code", "Original CellMiner TYMP log2(FPKM+1)", "Original CellMiner MTAP log2(FPKM+1)",
    "Number of technical MTA repeat sample names", "Arithmetic mean of raw MTA repeats",
    "Median of raw MTA repeats", "Sample SD of raw MTA repeats", "Minimum raw MTA repeat",
    "Maximum raw MTA repeat", "MTA mean in recovered official averaged table",
    "MTA SD in recovered official averaged table", "Absolute raw-mean versus official-mean difference",
    "Absolute raw-SD versus official-SD difference", "Number of MTA repeats equal to global observed floor",
    "Whether any MTA repeat equals global observed floor", "Whether every MTA repeat equals global observed floor",
    "NCI60 panel number", "NCI60 cell number", "Cell-line name in official averaged table"
  )
)
fwrite(dictionary, file.path(results_dir, "data_dictionary_v1.tsv"), sep = "\t", quote = FALSE, na = "NA")

status <- list(
  status = "PASS_PREPROCESSING_COMPLETE",
  generated_at = generated_at,
  input_dir = input_dir,
  independent_unit = "cell_line",
  analysis_performed = FALSE,
  all_lineages_n = nrow(paired),
  renal_n = nrow(renal),
  raw_mta_repeat_rows = nrow(raw_mta),
  raw_mta_unique_sample_names = uniqueN(raw_mta$SAMPLENAME),
  global_observed_mta_floor = global_floor,
  lines_with_any_floor_repeat = paired[mta_any_repeat_at_global_floor == TRUE, .N],
  lines_with_all_repeats_at_floor = paired[mta_all_repeats_at_global_floor == TRUE, .N],
  renal_cell_lines = renal$metabolomics_cell_line,
  checks_passed = sum(checks$pass),
  checks_total = nrow(checks)
)
write_json(status, file.path(qc_dir, "preprocessing_status_v1.json"), auto_unbox = TRUE, pretty = TRUE)

summary_lines <- c(
  "# Phase 102 NCI60 preprocessing summary",
  "",
  paste0("- Status: **", status$status, "**; ", status$checks_passed, "/", status$checks_total, " internal checks passed."),
  paste0("- Independent unit: cell line. Technical MTA repeats retained separately and aggregated once per line."),
  paste0("- All-lineage analysis table: ", status$all_lineages_n, " cell lines with complete TYMP, MTAP, and exact-MTA values."),
  paste0("- Renal-restricted table: ", status$renal_n, " lines: ", paste(status$renal_cell_lines, collapse = ", "), "."),
  paste0("- Exact-MTA repeat rows: ", status$raw_mta_repeat_rows, " across ", status$raw_mta_unique_sample_names, " sample names."),
  paste0("- Global observed MTA floor: ", format(status$global_observed_mta_floor, scientific = FALSE), "; any-floor lines: ", status$lines_with_any_floor_repeat, "; all-repeats-at-floor lines: ", status$lines_with_all_repeats_at_floor, "."),
  "- RNA values come from the original CellMiner XLS, not the thresholded cleaned CSV.",
  "- No association, inferential test, model, P value, FDR, or causal claim was produced in this module.",
  "- Future analyses are frozen at two levels: all lineages with tissue-aware sensitivity, and renal-only exploratory analysis with independent n=8."
)
writeLines(summary_lines, file.path(output_dir, "PHASE102_PREPROCESSING_SUMMARY_v1.md"), useBytes = TRUE)

output_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[!grepl("FILE_MANIFEST_SHA256_v1\\.tsv$", output_files)]
manifest <- data.table(
  relative_path = gsub("\\\\", "/", substring(normalizePath(output_files, winslash = "/"), nchar(normalizePath(output_dir, winslash = "/")) + 2L)),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest, character(1), algo = "sha256", file = TRUE)
)
setorder(manifest, relative_path)
fwrite(manifest, file.path(output_dir, "FILE_MANIFEST_SHA256_v1.tsv"), sep = "\t", quote = FALSE)

cat("PASS_PREPROCESSING_COMPLETE\n")
cat("all_lineages_n=", nrow(paired), " renal_n=", nrow(renal), " checks=", sum(checks$pass), "/", nrow(checks), "\n", sep = "")
