#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(copykat)
  library(jsonlite)
})

options(scipen = 100)
set.seed(26210)

parse_args <- function(x) {
  if (length(x) %% 2 != 0) stop("Arguments must be --key value pairs")
  out <- list()
  for (i in seq(1, length(x), by = 2)) out[[sub("^--", "", x[[i]])]] <- x[[i + 1]]
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- args[["input-dir"]]
out_dir <- args[["output-dir"]]
n_cores <- as.integer(ifelse(is.null(args[["n-cores"]]), "8", args[["n-cores"]]))
memory_gate <- "/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh"

if (is.null(input_dir) || is.null(out_dir)) stop("Required: --input-dir PATH --output-dir PATH")
matrix_file <- file.path(input_dir, "counts_genes_by_cells_v1.mtx.gz")
mapping_file <- file.path(input_dir, "cell_id_mapping_v1.tsv.gz")
input_gene_order_file <- file.path(input_dir, "gene_order_v1.tsv.gz")
metadata_file <- file.path(input_dir, "unit_metadata_v1.json")
for (path in c(matrix_file, mapping_file, input_gene_order_file, metadata_file)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(out_dir, "status.json")
gate_file <- file.path(out_dir, "memory_gate_log.tsv")
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
gene_coverage_qc <- "not_assessed"
n_genes_after_filter <- NA_integer_

write_status <- function(status, message = "", extra = list()) {
  payload <- c(list(
    schema_version = "part1_copykat_unit_status_v2",
    status = status,
    message = message,
    started_at = started_at,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    input_dir = normalizePath(input_dir, mustWork = TRUE),
    output_dir = normalizePath(out_dir, mustWork = FALSE)
  ), extra)
  write_json(payload, status_file, pretty = TRUE, auto_unbox = TRUE)
}

check_gate <- function(stage) {
  output <- system2(memory_gate, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(output, "status")
  if (is.null(exit_status)) exit_status <- 0L
  line <- if (length(output)) paste(output, collapse = " ") else "NO_OUTPUT"
  record <- data.table(stage = stage, gate_output = line, exit_status = exit_status)
  fwrite(record, gate_file, sep = "\t", append = file.exists(gate_file), col.names = !file.exists(gate_file))
  cat(line, "\n")
  if (exit_status != 0L || !grepl("PASS", line, fixed = TRUE)) stop("Memory gate blocked stage: ", stage)
}

main <- function() {
  write_status("running", "initializing")
  check_gate("read_unit_input")
  attrs <- fromJSON(metadata_file, simplifyVector = TRUE)
  if (attrs$schema_version != "part1_cnv_unit_metadata_v1") stop("Unexpected input schema")
  if (attrs$unit_type != "copykat") stop("Input unit_type is not copykat")
  mapping <- as.data.table(read.delim(gzfile(mapping_file), check.names = FALSE))
  input_gene_order <- as.data.table(read.delim(gzfile(input_gene_order_file), check.names = FALSE))
  mat <- as(readMM(gzfile(matrix_file)), "dgCMatrix")
  genes <- as.character(input_gene_order$gene)
  tool_ids <- as.character(mapping$tool_id)
  exact_ids <- as.character(mapping$exact_cell_id)
  roles <- as.character(mapping$role)
  sample_ids <- as.character(mapping$sample_id)
  patient_ids <- as.character(mapping$curated_patient_id)
  studies <- as.character(mapping$study)
  dimnames(mat) <- list(genes, tool_ids)
  if (anyDuplicated(tool_ids) || anyDuplicated(exact_ids)) stop("Duplicate cell IDs")
  if (any(!is.finite(mat@x)) || any(mat@x < 0) || any(abs(mat@x - round(mat@x)) > 1e-8)) stop("Invalid counts")

  fwrite(mapping, file.path(out_dir, "cell_id_mapping_v2.tsv"), sep = "\t")

  keep_cells <- Matrix::colSums(mat > 0) >= 200L
  mat <- mat[, keep_cells, drop = FALSE]
  kept_mapping <- mapping[match(colnames(mat), tool_id)]
  query_ids <- kept_mapping[role == "query", tool_id]
  reference_ids <- kept_mapping[role == "reference", tool_id]
  if (length(query_ids) < 50 || length(reference_ids) < 50) stop("Too few query/reference cells after filter")
  keep_genes <- Matrix::rowSums(mat > 0) / ncol(mat) > 0.05
  mat <- mat[keep_genes, , drop = FALSE]
  n_genes_after_filter <<- nrow(mat)
  gene_coverage_qc <<- if (nrow(mat) < 5000) "warning_below_v1_guard" else "pass_v1_guard"
  if (gene_coverage_qc != "pass_v1_guard") {
    warning(
      "Gene coverage below the retired v1 guard; continuing with unchanged CopyKAT parameters: ",
      nrow(mat)
    )
  }

  check_gate("densify_copykat_matrix")
  raw_matrix <- as.matrix(mat)
  storage.mode(raw_matrix) <- "numeric"
  check_gate("run_copykat")
  analysis_dir <- file.path(out_dir, "copykat_run")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(analysis_dir)
  result <- tryCatch(
    copykat::copykat(
      rawmat = raw_matrix,
      id.type = "S",
      cell.line = "no",
      ngene.chr = 5,
      min.gene.per.cell = 200,
      LOW.DR = 0.05,
      UP.DR = 0.1,
      win.size = 25,
      norm.cell.names = reference_ids,
      KS.cut = 0.1,
      sam.name = as.character(attrs$unit_id),
      distance = "euclidean",
      test.emd = "FALSE",
      output.seg = FALSE,
      plot.genes = FALSE,
      genome = "hg20",
      n.cores = n_cores
    ),
    finally = setwd(old_wd)
  )

  check_gate("save_copykat_result")
  saveRDS(result, file.path(out_dir, "copykat_result_v2.rds"), compress = "gzip")
  prediction <- as.data.table(result$prediction)
  if (!all(c("cell.names", "copykat.pred") %in% names(prediction))) {
    stop("Unexpected CopyKAT prediction columns: ", paste(names(prediction), collapse = ","))
  }
  prediction[, tool_id := as.character(cell.names)]
  prediction <- merge(prediction, mapping, by = "tool_id", all.x = TRUE, sort = FALSE)
  if (anyNA(prediction$exact_cell_id)) stop("CopyKAT returned unmapped tool IDs")
  if (anyDuplicated(prediction$tool_id)) stop("Duplicate CopyKAT predictions")
  prediction[, unit_id := as.character(attrs$unit_id)]
  prediction[, reference_mode := as.character(attrs$reference_mode)]
  fwrite(prediction, file.path(out_dir, "copykat_predictions_mapped_v2.tsv"), sep = "\t")

  summary <- data.table(
    parameter = c(
      "unit_id", "curated_patient_id", "study", "reference_mode",
      "n_query_input", "n_reference_input", "n_query_after_filter",
      "n_reference_after_filter", "n_genes_used", "n_predictions",
      "n_cores", "copykat_version", "seed", "exact_id_mapping_complete",
      "gene_coverage_qc", "retired_v1_guard_genes"
    ),
    value = as.character(c(
      attrs$unit_id, attrs$curated_patient_id, attrs$study, attrs$reference_mode,
      sum(roles == "query"), sum(roles == "reference"), length(query_ids),
      length(reference_ids), nrow(mat), nrow(prediction), n_cores,
      packageVersion("copykat"), 26210, TRUE,
      gene_coverage_qc, 5000
    ))
  )
  fwrite(summary, file.path(out_dir, "run_summary_v2.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_v2.txt"))
  writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(out_dir, "RUN_COMPLETE"))
  write_status("completed", "CopyKAT unit completed", list(
    unit_id = attrs$unit_id,
    n_query_cells = length(query_ids),
    n_reference_cells = length(reference_ids),
    n_genes_used = nrow(mat),
    n_predictions = nrow(prediction),
    gene_coverage_qc = gene_coverage_qc,
    retired_v1_guard_genes = 5000
  ))
  cat("PART1_COPYKAT_UNIT_V2_PASS unit=", attrs$unit_id, "\n", sep = "")
}

tryCatch(
  main(),
  error = function(e) {
    message_text <- conditionMessage(e)
    writeLines(message_text, file.path(out_dir, "ERROR.txt"))
    method_insufficient <- gene_coverage_qc == "warning_below_v1_guard" &&
      grepl("dim\\(X\\)|positive|正数", message_text)
    if (method_insufficient) {
      write_status("not_estimable", message_text, list(
        n_genes_used = n_genes_after_filter,
        gene_coverage_qc = gene_coverage_qc,
        retired_v1_guard_genes = 5000,
        reason_class = "copykat_internal_dimension_failure_low_coverage"
      ))
      writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(out_dir, "RUN_NOT_ESTIMABLE"))
      cat("PART1_COPYKAT_UNIT_V2_NOT_ESTIMABLE genes=", n_genes_after_filter,
          " message=", message_text, "\n", sep = "")
    } else {
      write_status("failed", message_text)
      stop(e)
    }
  }
)
