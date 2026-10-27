#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(infercnv)
  library(jsonlite)
})

options(scipen = 100)
set.seed(26210)

parse_args <- function(x) {
  if (length(x) %% 2 != 0) stop("Arguments must be --key value pairs")
  out <- list()
  for (i in seq(1, length(x), by = 2)) {
    key <- sub("^--", "", x[[i]])
    out[[key]] <- x[[i + 1]]
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- args[["input-dir"]]
out_dir <- args[["output-dir"]]
threads <- as.integer(ifelse(is.null(args[["threads"]]), "8", args[["threads"]]))
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

write_status <- function(status, message = "", extra = list()) {
  payload <- c(list(
    schema_version = "part1_infercnv_unit_status_v2",
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
  if (attrs$unit_type != "infercnv") stop("Input unit_type is not infercnv")
  mapping <- as.data.table(read.delim(gzfile(mapping_file), check.names = FALSE))
  input_gene_order <- as.data.table(read.delim(gzfile(input_gene_order_file), check.names = FALSE))
  mat <- as(readMM(gzfile(matrix_file)), "dgCMatrix")
  genes <- as.character(input_gene_order$gene)
  chromosomes <- as.character(input_gene_order$chromosome)
  starts <- as.numeric(input_gene_order$start)
  ends <- as.numeric(input_gene_order$end)
  tool_ids <- as.character(mapping$tool_id)
  exact_ids <- as.character(mapping$exact_cell_id)
  roles <- as.character(mapping$role)
  sample_ids <- as.character(mapping$sample_id)
  patient_ids <- as.character(mapping$curated_patient_id)
  studies <- as.character(mapping$study)
  dimnames(mat) <- list(genes, tool_ids)
  if (nrow(mat) != length(genes) || ncol(mat) != length(tool_ids)) stop("Matrix dimensions do not match names")
  if (anyDuplicated(tool_ids) || anyDuplicated(exact_ids)) stop("Duplicate cell IDs")
  if (any(!is.finite(mat@x)) || any(mat@x < 0) || any(abs(mat@x - round(mat@x)) > 1e-8)) stop("Invalid counts")
  if (sum(roles == "query") < 50 || sum(roles == "reference") < 50) stop("Too few query/reference cells")

  fwrite(mapping, file.path(out_dir, "cell_id_mapping_v2.tsv"), sep = "\t")

  expressed <- Matrix::rowSums(mat > 0) >= 5L
  mat <- mat[expressed, , drop = FALSE]
  gene_order <- data.table(
    gene = genes,
    chromosome = chromosomes,
    start = starts,
    end = ends
  )[expressed]
  gene_coverage_qc <- if (nrow(mat) < 10000) "warning_below_v1_guard" else "pass_v1_guard"
  if (gene_coverage_qc != "pass_v1_guard") {
    warning(
      "Gene coverage below the retired v1 guard; continuing with unchanged inferCNV parameters: ",
      nrow(mat)
    )
  }

  annotation <- data.table(
    tool_id = tool_ids,
    infercnv_group = ifelse(roles == "reference", "ref_normal", "obs_tumor")
  )
  annotation_file <- file.path(out_dir, "cell_annotations_tool_ids_v2.txt")
  gene_order_file <- file.path(out_dir, "gene_order_used_v2.txt")
  fwrite(annotation, annotation_file, sep = "\t", col.names = FALSE)
  fwrite(gene_order, gene_order_file, sep = "\t", col.names = FALSE)

  check_gate("create_infercnv_object")
  obj <- CreateInfercnvObject(
    raw_counts_matrix = mat,
    annotations_file = annotation_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = "ref_normal"
  )
  saveRDS(obj, file.path(out_dir, "infercnv_object_pre_run_v2.rds"), compress = "gzip")

  check_gate("run_infercnv")
  analysis_dir <- file.path(out_dir, "infercnv_run")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- infercnv::run(
    obj,
    cutoff = 0.1,
    out_dir = analysis_dir,
    cluster_by_groups = FALSE,
    cluster_references = FALSE,
    denoise = TRUE,
    HMM = TRUE,
    HMM_type = "i6",
    analysis_mode = "subclusters",
    BayesMaxPNormal = 0.5,
    num_threads = threads,
    per_chr_hmm_subclusters_references = FALSE,
    resume_mode = FALSE,
    no_plot = TRUE,
    no_prelim_plot = TRUE,
    up_to_step = 17
  )

  check_gate("save_infercnv_object")
  saveRDS(obj, file.path(out_dir, "infercnv_object_post_run_v2.rds"), compress = "gzip")
  summary <- data.table(
    parameter = c(
      "unit_id", "curated_patient_id", "study", "reference_mode",
      "n_query_cells", "n_reference_cells", "n_input_cells", "n_genes_used",
      "threads", "infercnv_version", "seed", "exact_id_mapping_complete",
      "gene_coverage_qc", "retired_v1_guard_genes"
    ),
    value = as.character(c(
      attrs$unit_id, attrs$curated_patient_id, attrs$study, attrs$reference_mode,
      sum(roles == "query"), sum(roles == "reference"), length(tool_ids), nrow(mat),
      threads, packageVersion("infercnv"), 26210, TRUE,
      gene_coverage_qc, 10000
    ))
  )
  fwrite(summary, file.path(out_dir, "run_summary_v2.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_v2.txt"))
  writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(out_dir, "RUN_COMPLETE"))
  write_status("completed", "inferCNV unit completed", list(
    unit_id = attrs$unit_id,
    n_query_cells = sum(roles == "query"),
    n_reference_cells = sum(roles == "reference"),
    n_genes_used = nrow(mat),
    gene_coverage_qc = gene_coverage_qc,
    retired_v1_guard_genes = 10000
  ))
  cat("PART1_INFERCNV_UNIT_V2_PASS unit=", attrs$unit_id, "\n", sep = "")
}

tryCatch(
  main(),
  error = function(e) {
    writeLines(conditionMessage(e), file.path(out_dir, "ERROR.txt"))
    write_status("failed", conditionMessage(e))
    stop(e)
  }
)
