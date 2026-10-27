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
required <- c("crosswalk", "preflight-dir", "source-root", "gene-order", "output-root")
missing <- required[!required %in% names(args)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

crosswalk_path <- normalizePath(args[["crosswalk"]], mustWork = TRUE)
preflight_dir <- normalizePath(args[["preflight-dir"]], mustWork = TRUE)
source_root <- normalizePath(args[["source-root"]], mustWork = TRUE)
gene_order_input <- normalizePath(args[["gene-order"]], mustWork = TRUE)
output_root <- args[["output-root"]]
unit_type_requested <- if (!is.null(args[["unit-type"]])) args[["unit-type"]] else "all"
unit_filter <- if (!is.null(args[["unit-filter"]]) && nzchar(args[["unit-filter"]])) {
  trimws(strsplit(args[["unit-filter"]], ",", fixed = TRUE)[[1]])
} else character()
if (!unit_type_requested %chin% c("all", "infercnv", "copykat")) stop("Invalid --unit-type")
if (dir.exists(output_root)) stop("No-overwrite gate: output root already exists: ", output_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

core_barcode <- function(x) {
  hit <- regmatches(x, regexpr("[ACGTN]+-[0-9]+$", x, perl = TRUE))
  hit[hit == ""] <- NA_character_
  toupper(hit)
}

write_gzip_from_file <- function(source, target_part) {
  input <- file(source, open = "rb")
  output <- gzfile(target_part, open = "wb", compression = 6)
  repeat {
    block <- readBin(input, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
  close(input)
  close(output)
}

stable_seed <- function(unit_type, unit_id) {
  as.integer(26210L + (sum(utf8ToInt(paste(unit_type, unit_id, sep = "|"))) %% 1000000L))
}

balanced_reference <- function(candidates, target, seed) {
  candidates <- unique(copy(candidates), by = "cell_id")
  setorder(candidates, sample_id, object_name, cell_id)
  if (nrow(candidates) <= target) return(candidates)
  set.seed(seed)
  samples <- sort(unique(candidates$sample_id))
  quota <- max(1L, target %/% length(samples))
  selected <- integer()
  for (sample_value in samples) {
    positions <- which(candidates$sample_id == sample_value)
    take <- min(quota, length(positions))
    selected <- c(selected, sample(positions, take))
  }
  selected <- unique(selected)
  if (length(selected) < target) {
    remaining <- setdiff(seq_len(nrow(candidates)), selected)
    selected <- c(selected, sample(remaining, min(target - length(selected), length(remaining))))
  }
  candidates[sort(selected[seq_len(min(target, length(selected)))])]
}

source_files <- list.files(source_root, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (!length(source_files)) stop("No source RDS files found")
source_map <- data.table(
  object_name = sub("[.]rds$", "", basename(source_files), ignore.case = TRUE),
  source_rds = normalizePath(source_files, mustWork = TRUE)
)
if (anyDuplicated(source_map$object_name)) {
  duplicates <- source_map[duplicated(object_name) | duplicated(object_name, fromLast = TRUE)]
  stop("Duplicate source object basenames: ", paste(unique(duplicates$object_name), collapse = ", "))
}
setkey(source_map, object_name)

cells <- fread(crosswalk_path, select = c(
  "cell_id", "dataset", "patient_id", "sample_id", "object_name", "barcode_original",
  "phase7_broad_annotation"
))
if (nrow(cells) != 585119L || anyDuplicated(cells$cell_id)) stop("Crosswalk conservation gate failed")
cells[, crosswalk_row := .I]
cells[, core_10x_barcode := core_barcode(barcode_original)]
if (anyNA(cells$core_10x_barcode)) stop("Missing terminal barcode in crosswalk")
cells[, cnv_tissue_class := "Tumor"]
cells[
  dataset == "GSE178481" &
    (grepl("Normal", sample_id, ignore.case = TRUE) | grepl("Normal", object_name, ignore.case = TRUE)),
  cnv_tissue_class := "Normal"
]
cells[, current_epithelial := phase7_broad_annotation == "Epithelial"]
cells[, current_lymphoid := phase7_broad_annotation %chin% c("T_cell", "NK_cell", "B_cell")]
required_objects <- unique(cells$object_name)
missing_objects <- setdiff(required_objects, source_map$object_name)
if (length(missing_objects)) stop("Missing source RDS objects: ", paste(head(missing_objects, 20L), collapse = ", "))

gene_order_all <- fread(gene_order_input, header = FALSE, col.names = c("gene", "chromosome", "start", "end"))
gene_order_all <- gene_order_all[!duplicated(gene)]

reference_candidates <- function(dataset_value, patient_value, mode) {
  if (mode == "matched_patient_normal_epithelial") {
    cells[dataset == dataset_value & patient_id == patient_value & current_epithelial & cnv_tissue_class == "Normal"]
  } else if (mode == "same_study_normal_epithelial_pool") {
    cells[dataset == dataset_value & current_epithelial & cnv_tissue_class == "Normal"]
  } else if (mode == "matched_patient_tumor_lymphoid") {
    cells[dataset == dataset_value & patient_id == patient_value & current_lymphoid & cnv_tissue_class == "Tumor"]
  } else {
    stop("Unsupported reference mode: ", mode)
  }
}

build_selected_matrix <- function(selected) {
  selected <- copy(selected)
  selected[, selected_order := .I]
  object_groups <- split(selected, selected$object_name)
  pieces_raw <- vector("list", length(object_groups))
  feature_lists <- vector("list", length(object_groups))
  source_records <- vector("list", length(object_groups))
  group_names <- names(object_groups)
  for (i in seq_along(object_groups)) {
    group <- as.data.table(object_groups[[i]])
    source_path <- source_map[group_names[[i]], source_rds]
    if (length(source_path) != 1L || is.na(source_path)) stop("Source lookup failed for ", group_names[[i]])
    obj <- readRDS(source_path)
    if (!inherits(obj, "Seurat")) stop("Source is not Seurat: ", source_path)
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
    if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
    source_keys <- core_barcode(colnames(counts))
    if (anyNA(source_keys) || anyDuplicated(source_keys)) stop("Invalid source barcode keys: ", source_path)
    indices <- match(group$core_10x_barcode, source_keys)
    if (anyNA(indices)) stop("Selected barcode absent from source: ", source_path)
    piece <- as(counts[, indices, drop = FALSE], "dgCMatrix")
    colnames(piece) <- group$cell_id
    pieces_raw[[i]] <- piece
    feature_lists[[i]] <- rownames(piece)
    source_records[[i]] <- data.table(
      object_name = group_names[[i]],
      source_rds = source_path,
      source_rds_bytes = file.info(source_path)$size,
      source_rds_sha256 = digest(file = source_path, algo = "sha256", serialize = FALSE)
    )
    rm(obj, counts, piece)
    invisible(gc())
  }
  shared_genes <- Reduce(intersect, feature_lists)
  gene_order_used <- gene_order_all[gene %chin% shared_genes]
  if (nrow(gene_order_used) < 10000L) stop("Too few shared ordered genes: ", nrow(gene_order_used))
  pieces <- vector("list", length(pieces_raw))
  for (i in seq_along(pieces_raw)) {
    indices <- match(gene_order_used$gene, rownames(pieces_raw[[i]]))
    if (anyNA(indices)) stop("Shared gene alignment failed")
    pieces[[i]] <- pieces_raw[[i]][indices, , drop = FALSE]
  }
  matrix <- do.call(cbind, pieces)
  matrix <- as(matrix[, match(selected$cell_id, colnames(matrix)), drop = FALSE], "dgCMatrix")
  colnames(matrix) <- selected$tool_id
  rownames(matrix) <- gene_order_used$gene
  if (any(!is.finite(matrix@x)) || any(matrix@x < 0) || any(abs(matrix@x - round(matrix@x)) > 1e-8)) {
    stop("Invalid raw counts")
  }
  list(matrix = matrix, gene_order = gene_order_used, sources = rbindlist(source_records))
}

write_unit <- function(unit_type, unit_row) {
  unit_id <- as.character(unit_row$unit_id)
  dataset_value <- as.character(unit_row$dataset)
  patient_value <- as.character(unit_row$patient_id)
  sample_value <- if (unit_type == "copykat") as.character(unit_row$sample_id) else NA_character_
  safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", unit_id)
  unit_dir <- file.path(output_root, unit_type, safe_id)
  if (dir.exists(unit_dir)) stop("Unit directory unexpectedly exists: ", unit_dir)
  dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)

  if (unit_type == "infercnv") {
    query <- cells[
      dataset == dataset_value & patient_id == patient_value & current_epithelial & cnv_tissue_class == "Tumor"
    ]
    target_reference <- 200L
  } else {
    query <- cells[
      dataset == dataset_value & patient_id == patient_value & sample_id == sample_value &
        current_epithelial & cnv_tissue_class == "Tumor"
    ]
    target_reference <- 1000L
  }
  if (nrow(query) != as.integer(unit_row$n_tumor_epithelial)) {
    stop("Query count mismatch for ", unit_id, ": ", nrow(query), " != ", unit_row$n_tumor_epithelial)
  }
  candidates <- reference_candidates(dataset_value, patient_value, as.character(unit_row$reference_mode))
  if (nrow(candidates) != as.integer(unit_row$n_reference_candidates)) {
    stop("Reference candidate count mismatch for ", unit_id)
  }
  reference <- balanced_reference(candidates, min(target_reference, nrow(candidates)), stable_seed(unit_type, unit_id))
  if (nrow(query) < 50L || nrow(reference) < 50L) stop("Eligibility count gate failed for ", unit_id)
  if (length(intersect(query$cell_id, reference$cell_id))) stop("Query/reference overlap for ", unit_id)

  selected <- rbind(
    copy(query)[, role := "query"],
    copy(reference)[, role := "reference"]
  )
  selected[, tool_id := ifelse(
    role == "query",
    sprintf("Q%07d", seq_len(.N)),
    sprintf("R%07d", seq_len(.N))
  ), by = role]
  built <- build_selected_matrix(selected)
  matrix <- built$matrix

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
    curated_patient_id = paste(dataset, patient_id, sep = "__"),
    study = dataset,
    source_tissue = cnv_tissue_class,
    cell_type_broad_v1 = ifelse(role == "query", "Epithelial cell", "Reference cell"),
    object_name,
    original_barcode = barcode_original
  )]
  mapping_part <- file.path(unit_dir, "cell_id_mapping_v1.part.tsv.gz")
  mapping_final <- file.path(unit_dir, "cell_id_mapping_v1.tsv.gz")
  fwrite(mapping, mapping_part, sep = "\t", quote = FALSE, compress = "gzip")
  if (!file.rename(mapping_part, mapping_final)) stop("Atomic mapping rename failed")

  gene_part <- file.path(unit_dir, "gene_order_v1.part.tsv.gz")
  gene_final <- file.path(unit_dir, "gene_order_v1.tsv.gz")
  fwrite(built$gene_order, gene_part, sep = "\t", quote = FALSE, compress = "gzip")
  if (!file.rename(gene_part, gene_final)) stop("Atomic gene-order rename failed")

  metadata <- list(
    schema_version = "part1_cnv_unit_metadata_v1",
    status = "completed",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    unit_type = unit_type,
    unit_id = unit_id,
    safe_unit_id = safe_id,
    curated_patient_id = paste(dataset_value, patient_value, sep = "__"),
    study = dataset_value,
    reference_mode = as.character(unit_row$reference_mode),
    seed = stable_seed(unit_type, unit_id),
    n_query_cells = nrow(query),
    n_reference_cells = nrow(reference),
    n_input_cells = nrow(selected),
    n_genes = nrow(matrix),
    matrix_nnz = length(matrix@x),
    source_objects = built$sources,
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
    exact_cell_id_policy = "Phase 8 global ID retained; object name plus terminal 10x barcode matched uniquely to source RDS"
  )
  metadata_part <- file.path(unit_dir, "unit_metadata_v1.part.json")
  metadata_final <- file.path(unit_dir, "unit_metadata_v1.json")
  write_json(metadata, metadata_part, pretty = TRUE, auto_unbox = TRUE)
  if (!file.rename(metadata_part, metadata_final)) stop("Atomic metadata rename failed")
  data.table(
    unit_type = unit_type,
    unit_id = unit_id,
    dataset = dataset_value,
    patient_id = patient_value,
    sample_id = sample_value,
    reference_mode = as.character(unit_row$reference_mode),
    n_query_cells = nrow(query),
    n_reference_cells = nrow(reference),
    n_genes = nrow(matrix),
    matrix_nnz = length(matrix@x),
    matrix_bytes = file.info(matrix_final)$size,
    status = "PASS"
  )
}

patient_manifest <- fread(file.path(preflight_dir, "external_infercnv_patient_preflight_v1.tsv"))
sample_manifest <- fread(file.path(preflight_dir, "external_copykat_sample_preflight_v1.tsv"))
patient_manifest <- patient_manifest[eligible_infercnv == TRUE]
sample_manifest <- sample_manifest[eligible_copykat == TRUE]
if (unit_type_requested == "infercnv") sample_manifest <- sample_manifest[0]
if (unit_type_requested == "copykat") patient_manifest <- patient_manifest[0]
if (length(unit_filter)) {
  patient_manifest <- patient_manifest[unit_id %chin% unit_filter]
  sample_manifest <- sample_manifest[unit_id %chin% unit_filter]
}
if (nrow(patient_manifest) + nrow(sample_manifest) == 0L) stop("No eligible units selected")

records <- list()
progress_file <- file.path(output_root, "input_build_progress_v1.tsv")
for (unit_type in c("infercnv", "copykat")) {
  manifest <- if (unit_type == "infercnv") patient_manifest else sample_manifest
  if (!nrow(manifest)) next
  for (i in seq_len(nrow(manifest))) {
    unit_row <- manifest[i]
    message(sprintf("[%s] build %s %s (%d/%d)", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), unit_type, unit_row$unit_id, i, nrow(manifest)))
    record <- write_unit(unit_type, unit_row)
    records[[length(records) + 1L]] <- record
    progress <- rbindlist(records, use.names = TRUE, fill = TRUE)
    progress_part <- paste0(progress_file, ".part")
    fwrite(progress, progress_part, sep = "\t", quote = FALSE)
    if (file.exists(progress_file)) unlink(progress_file)
    if (!file.rename(progress_part, progress_file)) stop("Progress atomic rename failed")
    invisible(gc())
  }
}

summary <- rbindlist(records, use.names = TRUE, fill = TRUE)
status <- list(
  status = "PASS",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  crosswalk = crosswalk_path,
  crosswalk_sha256 = digest(file = crosswalk_path, algo = "sha256", serialize = FALSE),
  source_root = source_root,
  gene_order_input = gene_order_input,
  unit_type_requested = unit_type_requested,
  unit_filter = unit_filter,
  n_units = nrow(summary),
  n_infercnv_units = sum(summary$unit_type == "infercnv"),
  n_copykat_units = sum(summary$unit_type == "copykat"),
  total_query_cells_across_units = sum(summary$n_query_cells),
  total_matrix_bytes = sum(summary$matrix_bytes),
  software = list(
    R = R.version.string,
    SeuratObject = as.character(packageVersion("SeuratObject")),
    Matrix = as.character(packageVersion("Matrix")),
    data.table = as.character(packageVersion("data.table")),
    jsonlite = as.character(packageVersion("jsonlite")),
    digest = as.character(packageVersion("digest"))
  )
)
status_part <- file.path(output_root, "FULL_INPUT_BUILD_STATUS.part.json")
status_final <- file.path(output_root, "FULL_INPUT_BUILD_STATUS.json")
write_json(status, status_part, pretty = TRUE, auto_unbox = TRUE)
if (!file.rename(status_part, status_final)) stop("Status atomic rename failed")
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(output_root, "RUN_COMPLETE.part"))
if (!file.rename(file.path(output_root, "RUN_COMPLETE.part"), file.path(output_root, "RUN_COMPLETE"))) stop("Completion rename failed")
cat("EXTERNAL_CNV_FULL_INPUT_BUILD_PASS units=", nrow(summary), "\n", sep = "")
