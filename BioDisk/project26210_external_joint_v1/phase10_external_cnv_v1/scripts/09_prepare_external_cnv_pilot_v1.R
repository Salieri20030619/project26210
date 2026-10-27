#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
})

parse_args <- function(x) {
  if (length(x) %% 2L != 0L) stop("Arguments must be --key value pairs")
  out <- list()
  for (i in seq(1L, length(x), by = 2L)) out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("crosswalk", "source-rds", "gene-order", "output-root")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

crosswalk_path <- normalizePath(args[["crosswalk"]], mustWork = TRUE)
source_rds <- normalizePath(args[["source-rds"]], mustWork = TRUE)
gene_order_input <- normalizePath(args[["gene-order"]], mustWork = TRUE)
output_root <- args[["output-root"]]
dataset_target <- "GSE202374"
patient_target <- "RS01985670"
sample_target <- "P10"
unit_ids <- c(infercnv = paste(dataset_target, patient_target, sep = "__"), copykat = paste(dataset_target, sample_target, sep = "__"))

core_barcode <- function(x) {
  hit <- regmatches(x, regexpr("[ACGTN]+-[0-9]+$", x, perl = TRUE))
  hit[hit == ""] <- NA_character_
  toupper(hit)
}

write_gzip_from_file <- function(source, target_part) {
  input <- file(source, open = "rb")
  output <- gzfile(target_part, open = "wb", compression = 6)
  on.exit({
    try(close(input), silent = TRUE)
    try(close(output), silent = TRUE)
  }, add = TRUE)
  repeat {
    block <- readBin(input, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
  close(input)
  close(output)
}

write_unit <- function(unit_type, query, reference, counts, gene_order) {
  unit_id <- unit_ids[[unit_type]]
  safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", unit_id)
  unit_dir <- file.path(output_root, unit_type, safe_id)
  if (dir.exists(unit_dir)) stop("No-overwrite gate: unit directory exists: ", unit_dir)
  dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)

  selected <- rbind(
    copy(query)[, role := "query"],
    copy(reference)[, role := "reference"]
  )
  selected[, tool_id := c(
    sprintf("Q%07d", seq_len(sum(role == "query"))),
    sprintf("R%07d", seq_len(sum(role == "reference")))
  )]
  setkey(selected, core_10x_barcode)
  source_lookup <- data.table(source_cell_name = colnames(counts), core_10x_barcode = core_barcode(colnames(counts)))
  if (anyNA(source_lookup$core_10x_barcode) || anyDuplicated(source_lookup$core_10x_barcode)) {
    stop("Source RDS terminal barcode keys are missing or duplicated")
  }
  selected <- source_lookup[selected, on = "core_10x_barcode"]
  if (anyNA(selected$source_cell_name)) stop("One or more selected cells are absent from source RDS")
  if (anyDuplicated(selected$source_cell_name) || anyDuplicated(selected$cell_id)) stop("Selected cell mapping is not unique")

  cell_indices <- match(selected$source_cell_name, colnames(counts))
  gene_indices <- match(gene_order$gene, rownames(counts))
  use_gene <- !is.na(gene_indices)
  gene_order_used <- gene_order[use_gene]
  gene_indices <- gene_indices[use_gene]
  matrix <- as(counts[gene_indices, cell_indices, drop = FALSE], "dgCMatrix")
  colnames(matrix) <- selected$tool_id
  rownames(matrix) <- gene_order_used$gene
  if (nrow(matrix) < 10000L) stop("Too few ordered genes: ", nrow(matrix))
  if (any(!is.finite(matrix@x)) || any(matrix@x < 0) || any(abs(matrix@x - round(matrix@x)) > 1e-8)) {
    stop("Invalid raw counts")
  }

  matrix_final <- file.path(unit_dir, "counts_genes_by_cells_v1.mtx.gz")
  matrix_plain_part <- file.path(unit_dir, "counts_genes_by_cells_v1.part.mtx")
  matrix_gzip_part <- file.path(unit_dir, "counts_genes_by_cells_v1.part.mtx.gz")
  Matrix::writeMM(matrix, matrix_plain_part)
  write_gzip_from_file(matrix_plain_part, matrix_gzip_part)
  unlink(matrix_plain_part)
  if (!file.rename(matrix_gzip_part, matrix_final)) stop("Atomic matrix rename failed")

  mapping <- selected[, .(
    tool_id,
    exact_cell_id = cell_id,
    role,
    sample_id,
    curated_patient_id = paste(dataset_target, patient_target, sep = "__"),
    study = dataset_target,
    source_tissue = "Tumor",
    cell_type_broad_v1 = ifelse(role == "query", "Epithelial cell", "Lymphoid reference"),
    source_cell_name,
    source_rds = source_rds
  )]
  mapping_part <- file.path(unit_dir, "cell_id_mapping_v1.part.tsv.gz")
  mapping_final <- file.path(unit_dir, "cell_id_mapping_v1.tsv.gz")
  fwrite(mapping, mapping_part, sep = "\t", quote = FALSE, compress = "gzip")
  if (!file.rename(mapping_part, mapping_final)) stop("Atomic mapping rename failed")

  gene_part <- file.path(unit_dir, "gene_order_v1.part.tsv.gz")
  gene_final <- file.path(unit_dir, "gene_order_v1.tsv.gz")
  fwrite(gene_order_used, gene_part, sep = "\t", quote = FALSE, compress = "gzip")
  if (!file.rename(gene_part, gene_final)) stop("Atomic gene-order rename failed")

  metadata <- list(
    schema_version = "part1_cnv_unit_metadata_v1",
    status = "completed",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    unit_type = unit_type,
    unit_id = unit_id,
    safe_unit_id = safe_id,
    curated_patient_id = paste(dataset_target, patient_target, sep = "__"),
    study = dataset_target,
    reference_mode = "matched_patient_tumor_lymphoid",
    seed = 26210,
    n_query_cells = nrow(query),
    n_reference_cells = nrow(reference),
    n_input_cells = nrow(selected),
    n_genes = nrow(matrix),
    matrix_nnz = length(matrix@x),
    source_rds = source_rds,
    source_rds_sha256 = digest(file = source_rds, algo = "sha256", serialize = FALSE),
    crosswalk_sha256 = digest(file = crosswalk_path, algo = "sha256", serialize = FALSE),
    matrix_mtx = matrix_final,
    matrix_mtx_size_bytes = file.info(matrix_final)$size,
    matrix_mtx_sha256 = digest(file = matrix_final, algo = "sha256", serialize = FALSE),
    mapping_tsv = mapping_final,
    mapping_tsv_size_bytes = file.info(mapping_final)$size,
    mapping_tsv_sha256 = digest(file = mapping_final, algo = "sha256", serialize = FALSE),
    gene_order_tsv = gene_final,
    gene_order_tsv_size_bytes = file.info(gene_final)$size,
    gene_order_tsv_sha256 = digest(file = gene_final, algo = "sha256", serialize = FALSE),
    exact_cell_id_policy = "Phase 8 global ID retained; GSM/sample plus terminal 10x barcode matched uniquely to source RDS"
  )
  metadata_part <- file.path(unit_dir, "unit_metadata_v1.part.json")
  metadata_final <- file.path(unit_dir, "unit_metadata_v1.json")
  write_json(metadata, metadata_part, pretty = TRUE, auto_unbox = TRUE)
  if (!file.rename(metadata_part, metadata_final)) stop("Atomic metadata rename failed")
  metadata
}

cells <- fread(crosswalk_path, select = c(
  "cell_id", "dataset", "patient_id", "sample_id", "barcode_original",
  "phase7_broad_annotation"
))
pilot <- cells[
  dataset == dataset_target & patient_id == patient_target & sample_id == sample_target
]
if (nrow(pilot) == 0L) stop("Pilot sample absent from crosswalk")
pilot[, core_10x_barcode := core_barcode(barcode_original)]
if (anyNA(pilot$core_10x_barcode) || anyDuplicated(pilot$core_10x_barcode)) stop("Pilot keys invalid")
query <- pilot[phase7_broad_annotation == "Epithelial"]
reference_candidates <- pilot[phase7_broad_annotation %chin% c("T_cell", "NK_cell", "B_cell")]
if (nrow(query) != 55L) stop("Expected 55 query cells, found ", nrow(query))
if (nrow(reference_candidates) != 2285L) stop("Expected 2285 reference candidates, found ", nrow(reference_candidates))

set.seed(26210)
reference_infercnv <- reference_candidates[sample(seq_len(.N), 200L)]
set.seed(26210)
reference_copykat <- reference_candidates[sample(seq_len(.N), 1000L)]

obj <- readRDS(source_rds)
if (!inherits(obj, "Seurat")) stop("Source object is not Seurat")
counts <- LayerData(obj, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
gene_order <- fread(gene_order_input, header = FALSE, col.names = c("gene", "chromosome", "start", "end"))
gene_order <- gene_order[!duplicated(gene)]

records <- list(
  write_unit("infercnv", query, reference_infercnv, counts, gene_order),
  write_unit("copykat", query, reference_copykat, counts, gene_order)
)
summary_part <- file.path(output_root, "pilot_input_summary_v1.part.json")
summary_final <- file.path(output_root, "pilot_input_summary_v1.json")
write_json(records, summary_part, pretty = TRUE, auto_unbox = TRUE)
if (!file.rename(summary_part, summary_final)) stop("Atomic pilot summary rename failed")
cat("EXTERNAL_CNV_PILOT_INPUT_PASS\n")
