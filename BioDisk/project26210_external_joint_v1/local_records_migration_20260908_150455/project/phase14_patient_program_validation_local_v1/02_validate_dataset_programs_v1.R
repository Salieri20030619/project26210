#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 02_validate_dataset_programs_v1.R <dataset> <output_dir>")
dataset_current <- args[[1L]]
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)

required <- c(
  "tables/signature_gene_availability_v1.tsv",
  "tables/patient_program_summary_v1.tsv",
  "tables/sample_program_summary_v1.tsv",
  "tables/patient_compartment_counts_v1.tsv",
  "tables/sample_compartment_counts_v1.tsv",
  "tables/f1_normal_epithelial_thresholds_v1.tsv",
  "tables/f1_malignant_prevalence_v1.tsv",
  "tables/f1_phase8_concordance_v1.tsv",
  paste0("objects/", dataset_current, "_cell_scores_compact_v1.rds"),
  paste0("status/", dataset_current, "_STATUS_v1.tsv")
)
checks <- list()
add_check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, result = ifelse(pass, "PASS", "FAIL"), detail = as.character(detail))
}

paths <- file.path(output_dir, required)
add_check("required_outputs_nonempty", all(file.exists(paths) & file.info(paths)$size > 0), paste(required, collapse = ";"))
availability <- fread(file.path(output_dir, required[[1L]]))
patient <- fread(file.path(output_dir, required[[2L]]))
sample <- fread(file.path(output_dir, required[[3L]]))
f1 <- fread(file.path(output_dir, "tables/f1_phase8_concordance_v1.tsv"))
scores <- readRDS(file.path(output_dir, paste0("objects/", dataset_current, "_cell_scores_compact_v1.rds")))
add_check("dataset_identity", all(c(availability$dataset, patient$dataset, sample$dataset, scores$dataset) == dataset_current), dataset_current)
add_check(
  "all_nonempty_signature_contexts_pass",
  all(availability[n_cells > 0, gate] == "PASS"),
  sprintf("contexts=%d empty_contexts=%d", availability[n_cells > 0, .N], availability[n_cells == 0, .N])
)
add_check("f1_exact_phase8_reproduction", nrow(f1) == 1L && f1$gate == "PASS" && f1$max_abs_difference <= 1e-8, paste(f1$max_abs_difference, collapse = ";"))
add_check("patient_keys_complete", !anyNA(patient$patient_key) && all(nzchar(patient$patient_key)), sprintf("patients=%d", uniqueN(patient$patient_key)))
add_check("sample_ids_complete", !anyNA(sample$sample_id) && all(nzchar(sample$sample_id)), sprintf("samples=%d", uniqueN(sample$sample_id)))
add_check("compact_scores_unique_cells", !anyDuplicated(scores$cell_id), sprintf("cells=%d", nrow(scores)))

overall <- all(vapply(checks, function(x) x$result == "PASS", logical(1L)))
checks[[length(checks) + 1L]] <- data.table(check = "OVERALL", result = ifelse(overall, "PASS", "FAIL"), detail = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
result <- rbindlist(checks)
path <- file.path(output_dir, "status", paste0(dataset_current, "_INDEPENDENT_VALIDATION_v1.tsv"))
part <- paste0(path, ".part")
if (file.exists(path) || file.exists(part)) stop("Refusing existing validator output")
fwrite(result, part, sep = "\t", quote = TRUE)
if (!file.rename(part, path)) stop("Atomic rename failed")
if (!overall) stop("Independent validation failed")
