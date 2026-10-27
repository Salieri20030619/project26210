#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, required = TRUE) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    if (required) stop("Missing argument: ", flag)
    return(NULL)
  }
  args[[idx + 1L]]
}

project_root <- normalizePath(arg_value("--project-root"), winslash = "/", mustWork = TRUE)
run_id <- arg_value("--run-id")
only_dataset <- arg_value("--dataset", required = FALSE)

phase8_root <- file.path(project_root, "phase8_spatial_alignment_local_v1")
reference_root <- file.path(phase8_root, "frozen_reference_inputs_v1")
input_root <- file.path(
  project_root,
  "phase7_annotation_local_v2", "runs", "20260810_234538_376", "objects_annotated"
)
run_root <- file.path(phase8_root, "annotation_spatial_crosswalk_v1", "runs", run_id)
objects_out <- file.path(run_root, "objects_annotated")
tables_out <- file.path(run_root, "tables")
status_out <- file.path(run_root, "status")
dir.create(objects_out, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_out, recursive = TRUE, showWarnings = FALSE)
dir.create(status_out, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
if (!is.null(only_dataset)) {
  if (!only_dataset %in% datasets) stop("Unsupported dataset: ", only_dataset)
  datasets <- only_dataset
}

input_paths <- setNames(
  file.path(input_root, paste0(datasets, "_annotated_local_v2.rds")),
  datasets
)
if (any(!file.exists(input_paths))) {
  stop("Missing input objects: ", paste(input_paths[!file.exists(input_paths)], collapse = "; "))
}

marker_path <- file.path(reference_root, "reference_v2", "reference_v2_marker_ranking.tsv")
f1_path <- file.path(reference_root, "cnmf_rank7", "rank7_program_top100_genes_v1.tsv")
state_qc_path <- file.path(
  reference_root, "cell2location_reference_v3", "reference_v3_state_qc.tsv"
)
for (path in c(marker_path, f1_path, state_qc_path)) {
  if (!file.exists(path)) stop("Missing frozen reference input: ", path)
}

reference_states <- fread(state_qc_path)$state_v2
expected_states <- c(
  "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
  "Epithelial-Normal", "F1-high Malignant", "F1-low Malignant", "Mast cell",
  "Monocyte", "NK cell", "Other Macrophage", "SMC/Pericyte",
  "SPP1-state Macrophage", "T cell"
)
if (!setequal(reference_states, expected_states)) {
  stop("Frozen reference states do not match the expected 14-state ontology")
}

reference_markers <- fread(marker_path)
reference_markers <- reference_markers[rank <= 25L & specificity_log1p_cpm > 0]
if (!setequal(unique(reference_markers$state_v2), expected_states)) {
  stop("Top-25 marker table does not cover all 14 states")
}

f1_top <- fread(f1_path)[factor == "F1" & gene_rank <= 100L]
if (nrow(f1_top) != 100L || anyDuplicated(f1_top$gene)) {
  stop("F1 top-100 loading table is incomplete or duplicated")
}
if (any(!is.finite(f1_top$loading)) || any(f1_top$loading < 0)) {
  stop("Invalid F1 loadings")
}

spp1_program <- c("SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL")
c1qc_program <- c("C1QA", "C1QB", "C1QC", "APOE")
normal_epithelial <- c(
  "Proximal_tubule", "Collecting_duct", "Loop_of_Henle", "Distal_tubule", "Podocyte_like"
)
t_cell_labels <- c("CD4_naive_memory", "CD8_cytotoxic", "Treg", "T_cell_unresolved")
dendritic_labels <- c("cDC1", "cDC2", "Dendritic_pDC")
macrophage_labels <- c("Macrophage_C1QC", "Macrophage_SPP1")

atomic_fwrite <- function(x, path) {
  part <- paste0(path, ".part")
  if (file.exists(path) || file.exists(part)) stop("Refusing existing table output: ", path)
  fwrite(x, part, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic table rename failed: ", path)
}

write_status <- function(dataset, status, fields = list()) {
  path <- file.path(status_out, paste0(dataset, "_status.tsv"))
  part <- paste0(path, ".", Sys.getpid(), ".part")
  row <- as.data.table(c(list(
    dataset = dataset,
    status = status,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), fields))
  fwrite(row, part, sep = "\t", quote = FALSE, na = "NA")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace mutable status: ", path)
  if (!file.rename(part, path)) stop("Status rename failed: ", path)
}

percent_rank_by <- function(values, groups, eligible) {
  result <- rep(NA_real_, length(values))
  positions <- which(eligible & is.finite(values) & !is.na(groups) & nzchar(groups))
  if (!length(positions)) return(result)
  split_positions <- split(positions, groups[positions], drop = TRUE)
  for (pos in split_positions) {
    result[pos] <- rank(values[pos], ties.method = "average") / length(pos)
  }
  result
}

safe_colmeans <- function(expression_matrix, genes) {
  present <- intersect(genes, rownames(expression_matrix))
  if (!length(present)) return(rep(NA_real_, ncol(expression_matrix)))
  as.numeric(Matrix::colMeans(expression_matrix[present, , drop = FALSE]))
}

weighted_score <- function(expression_matrix, gene_weights) {
  present <- intersect(gene_weights$gene, rownames(expression_matrix))
  if (!length(present)) return(rep(NA_real_, ncol(expression_matrix)))
  weights <- gene_weights$loading[match(present, gene_weights$gene)]
  weights <- weights / sum(weights)
  as.numeric(Matrix::crossprod(weights, expression_matrix[present, , drop = FALSE]))
}

direct_state <- function(fine) {
  result <- rep(NA_character_, length(fine))
  result[fine == "B_cell"] <- "B cell"
  result[fine %in% dendritic_labels] <- "Dendritic cell"
  result[fine == "Endothelial"] <- "Endothelial"
  result[fine %in% normal_epithelial] <- "Epithelial-Normal"
  result[fine == "Mast_cell"] <- "Mast cell"
  result[fine == "Monocyte_FCN1"] <- "Monocyte"
  result[fine == "NK_cell"] <- "NK cell"
  result[fine == "Pericyte_SMC"] <- "SMC/Pericyte"
  result[fine %in% t_cell_labels] <- "T cell"
  result
}

display_label <- function(state, status) {
  result <- state
  result[state == "C1QC-state Macrophage"] <- "APOE+C1QC+ TAM"
  result[state == "SPP1-state Macrophage"] <- "SPP1+ TAM"
  result[state == "F1-high Malignant"] <- "F1-high malignant candidate"
  result[state == "F1-low Malignant"] <- "F1-low malignant candidate"
  result[status == "auxiliary_not_in_locked14"] <- "Auxiliary/non-14-state"
  result[status == "unresolved"] <- "Unresolved"
  result[status == "f1_unresolved_lt20_tumor_like_cells_per_patient"] <-
    "Tumor-like epithelial; F1 unresolved (<20 cells/patient)"
  result
}

all_cluster_scores <- list()
all_marker_overlap <- list()
all_summaries <- list()

for (dataset_current in datasets) {
  input_path <- input_paths[[dataset_current]]
  output_path <- file.path(objects_out, paste0(dataset_current, "_spatial_annotated_v1.rds"))
  output_part <- paste0(output_path, ".part")
  if (file.exists(output_path) || file.exists(output_part)) {
    stop("Refusing existing object output: ", output_path)
  }
  write_status(dataset_current, "loading", list(input_bytes = file.info(input_path)$size))
  input_bytes_before <- file.info(input_path)$size
  object <- readRDS(input_path)
  if (!inherits(object, "Seurat")) stop("Input is not a Seurat object: ", dataset_current)

  metadata <- object[[]]
  required_metadata <- c(
    "dataset", "patient_id", "sample_id", "phase7_cluster_res0_8",
    "phase7_broad_annotation", "phase7_fine_annotation", "phase7_annotation_confidence"
  )
  missing_metadata <- setdiff(required_metadata, colnames(metadata))
  if (length(missing_metadata)) {
    stop("Missing metadata in ", dataset_current, ": ", paste(missing_metadata, collapse = ", "))
  }
  if (any(metadata$dataset != dataset_current)) stop("Dataset metadata mismatch")
  if (anyNA(metadata$patient_id) || any(!nzchar(as.character(metadata$patient_id)))) {
    stop("Missing patient_id in ", dataset_current)
  }

  expression_matrix <- SeuratObject::LayerData(object[["RNA"]], layer = "data")
  if (!inherits(expression_matrix, "sparseMatrix")) expression_matrix <- as(expression_matrix, "dgCMatrix")
  if (ncol(expression_matrix) != nrow(metadata)) stop("Expression/metadata cell mismatch")

  fine <- as.character(metadata$phase7_fine_annotation)
  patient <- as.character(metadata$patient_id)
  sample <- as.character(metadata$sample_id)
  cluster <- as.character(metadata$phase7_cluster_res0_8)
  phase7_confidence <- as.character(metadata$phase7_annotation_confidence)
  n_cells <- nrow(metadata)

  state <- direct_state(fine)
  status <- rep("unresolved", n_cells)
  status[!is.na(state)] <- "direct_locked14"

  macrophage <- fine %in% macrophage_labels
  spp1_score <- rep(NA_real_, n_cells)
  c1qc_score <- rep(NA_real_, n_cells)
  spp1_pct_patient <- rep(NA_real_, n_cells)
  c1qc_pct_patient <- rep(NA_real_, n_cells)
  spp1_pct_sample <- rep(NA_real_, n_cells)
  c1qc_pct_sample <- rep(NA_real_, n_cells)
  if (any(macrophage)) {
    spp1_score[macrophage] <- safe_colmeans(expression_matrix[, macrophage, drop = FALSE], spp1_program)
    c1qc_score[macrophage] <- safe_colmeans(expression_matrix[, macrophage, drop = FALSE], c1qc_program)
    spp1_pct_patient <- percent_rank_by(spp1_score, patient, macrophage)
    c1qc_pct_patient <- percent_rank_by(c1qc_score, patient, macrophage)
    spp1_pct_sample <- percent_rank_by(spp1_score, sample, macrophage)
    c1qc_pct_sample <- percent_rank_by(c1qc_score, sample, macrophage)
    spp1_high <- macrophage & spp1_pct_patient >= 0.75 & spp1_pct_patient >= c1qc_pct_patient
    c1qc_high <- macrophage & c1qc_pct_patient >= 0.75 & !spp1_high
    state[macrophage] <- "Other Macrophage"
    state[spp1_high] <- "SPP1-state Macrophage"
    state[c1qc_high] <- "C1QC-state Macrophage"
    status[macrophage] <- "macrophage_program_projected"
  }

  tumor_like <- fine == "Epithelial_tumor_like"
  f1_score <- rep(NA_real_, n_cells)
  f1_pct_patient <- rep(NA_real_, n_cells)
  f1_pct_sample <- rep(NA_real_, n_cells)
  f1_eligible <- rep(FALSE, n_cells)
  f1_overlap <- intersect(f1_top$gene, rownames(expression_matrix))
  if (length(f1_overlap) < 70L) {
    stop("F1 top-100 overlap below 70 genes in ", dataset_current, ": ", length(f1_overlap))
  }
  if (any(tumor_like)) {
    f1_score[tumor_like] <- weighted_score(
      expression_matrix[, tumor_like, drop = FALSE], f1_top
    )
    tumor_counts <- table(patient[tumor_like])
    eligible_patients <- names(tumor_counts)[tumor_counts >= 20L]
    f1_eligible <- tumor_like & patient %in% eligible_patients
    f1_pct_patient <- percent_rank_by(f1_score, patient, f1_eligible)
    f1_pct_sample <- percent_rank_by(f1_score, sample, f1_eligible)
    high <- f1_eligible & f1_pct_patient >= 0.75
    low <- f1_eligible & f1_pct_patient < 0.75
    state[high] <- "F1-high Malignant"
    state[low] <- "F1-low Malignant"
    status[high | low] <- "f1_program_projected_no_external_cnv"
    status[tumor_like & !f1_eligible] <- "f1_unresolved_lt20_tumor_like_cells_per_patient"
  }

  auxiliary <- fine %in% c(
    "Fibroblast", "Plasma_cell", "Neutrophil", "Myeloid_unresolved",
    "Epithelial_unresolved", "Epithelial_alveolar_like", "Epithelial_ciliated_like",
    "Epithelial_hepatocyte_like", "Unresolved", "Unresolved_low_cell_count"
  )
  state[auxiliary] <- NA_character_
  status[auxiliary] <- "auxiliary_not_in_locked14"

  if (any(!is.na(state) & !state %in% expected_states)) stop("Unexpected Phase 8 state")
  if (anyNA(status)) stop("Missing Phase 8 status")

  tymp <- rep(NA_real_, n_cells)
  if ("TYMP" %in% rownames(expression_matrix)) {
    tymp <- as.numeric(expression_matrix["TYMP", ])
  }

  confidence <- rep("low", n_cells)
  direct <- status == "direct_locked14"
  confidence[direct & phase7_confidence == "high"] <- "high"
  confidence[direct & phase7_confidence != "high"] <- "moderate"
  confidence[status == "macrophage_program_projected"] <- "moderate"
  confidence[status == "f1_program_projected_no_external_cnv"] <- "low"
  confidence[status == "auxiliary_not_in_locked14" & phase7_confidence == "high"] <- "moderate"

  marker_score_rows <- list()
  marker_overlap_rows <- list()
  for (state_current in expected_states) {
    genes_requested <- reference_markers[state_v2 == state_current, gene]
    genes_present <- intersect(genes_requested, rownames(expression_matrix))
    marker_overlap_rows[[state_current]] <- data.table(
      dataset = dataset_current,
      reference_state = state_current,
      markers_requested = length(genes_requested),
      markers_present = length(genes_present),
      overlap_fraction = length(genes_present) / length(genes_requested),
      marker_genes_present = paste(genes_present, collapse = ";")
    )
    score <- safe_colmeans(expression_matrix, genes_present)
    score_sum <- rowsum(score, group = cluster, reorder = FALSE)
    cluster_n <- table(factor(cluster, levels = rownames(score_sum)))
    marker_score_rows[[state_current]] <- data.table(
      dataset = dataset_current,
      cluster = rownames(score_sum),
      reference_state = state_current,
      mean_top25_marker_score = as.numeric(score_sum[, 1L]) / as.numeric(cluster_n),
      markers_present = length(genes_present)
    )
  }
  dataset_cluster_scores <- rbindlist(marker_score_rows)
  setorder(dataset_cluster_scores, cluster, -mean_top25_marker_score, reference_state)
  dataset_cluster_scores[, reference_score_rank := seq_len(.N), by = cluster]
  top_by_cluster <- dataset_cluster_scores[reference_score_rank == 1L]
  top_match <- match(cluster, top_by_cluster$cluster)

  new_metadata <- data.frame(
    phase8_spatial_reference_state = state,
    phase8_spatial_display_label = display_label(state, status),
    phase8_spatial_mapping_status = status,
    phase8_spatial_annotation_confidence = confidence,
    phase8_spp1_program_score = spp1_score,
    phase8_c1qc_program_score = c1qc_score,
    phase8_spp1_percentile_patient = spp1_pct_patient,
    phase8_c1qc_percentile_patient = c1qc_pct_patient,
    phase8_spp1_percentile_sample = spp1_pct_sample,
    phase8_c1qc_percentile_sample = c1qc_pct_sample,
    phase8_f1_projection_score = f1_score,
    phase8_f1_percentile_patient = f1_pct_patient,
    phase8_f1_percentile_sample = f1_pct_sample,
    phase8_f1_eligible_patient = f1_eligible,
    phase8_tymp_log_normalized = tymp,
    phase8_reference_marker_top1 = top_by_cluster$reference_state[top_match],
    phase8_reference_marker_top1_score = top_by_cluster$mean_top25_marker_score[top_match],
    phase8_malignancy_claim = FALSE,
    phase8_phase7_label_preserved = TRUE,
    phase8_annotation_version = "spatial_crosswalk_v1",
    row.names = rownames(metadata),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  object <- SeuratObject::AddMetaData(object, metadata = new_metadata)

  write_status(dataset_current, "writing", list(
    cells = n_cells,
    f1_top100_genes_present = length(f1_overlap)
  ))
  saveRDS(object, output_part, compress = "gzip")
  if (!file.rename(output_part, output_path)) stop("Atomic RDS rename failed: ", output_path)
  rm(object, expression_matrix)
  gc(verbose = FALSE)

  check <- readRDS(output_path)
  check_metadata <- check[[]]
  required_phase8 <- colnames(new_metadata)
  if (!all(required_phase8 %in% colnames(check_metadata))) stop("Output metadata validation failed")
  if (nrow(check_metadata) != n_cells) stop("Output cell-count validation failed")
  if (anyNA(check_metadata$phase8_spatial_mapping_status)) stop("Output status contains NA")
  if (!all(check_metadata$phase8_phase7_label_preserved)) stop("Phase 7 preservation flag failed")
  rm(check, check_metadata)
  gc(verbose = FALSE)

  if (file.info(input_path)$size != input_bytes_before) stop("Input object size changed")
  summary_dt <- as.data.table(new_metadata)[, .N, by = .(
    phase8_spatial_reference_state,
    phase8_spatial_display_label,
    phase8_spatial_mapping_status,
    phase8_spatial_annotation_confidence
  )]
  summary_dt[, dataset := dataset_current]
  setcolorder(summary_dt, c("dataset", setdiff(names(summary_dt), "dataset")))
  all_summaries[[dataset_current]] <- summary_dt
  all_cluster_scores[[dataset_current]] <- dataset_cluster_scores
  all_marker_overlap[[dataset_current]] <- rbindlist(marker_overlap_rows)
  write_status(dataset_current, "PASS", list(
    cells = n_cells,
    output_bytes = file.info(output_path)$size,
    input_bytes_unchanged = TRUE,
    phase7_labels_preserved = TRUE,
    malignancy_claims = 0L
  ))
}

atomic_fwrite(rbindlist(all_summaries, fill = TRUE), file.path(tables_out, "cell_state_summary.tsv"))
atomic_fwrite(
  rbindlist(all_cluster_scores, fill = TRUE),
  file.path(tables_out, "cluster_reference_marker_scores.tsv")
)
atomic_fwrite(
  rbindlist(all_marker_overlap, fill = TRUE),
  file.path(tables_out, "reference_marker_overlap.tsv")
)

software <- data.table(
  component = c("R", "SeuratObject", "Matrix", "data.table"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("SeuratObject")),
    as.character(packageVersion("Matrix")),
    as.character(packageVersion("data.table"))
  )
)
atomic_fwrite(software, file.path(tables_out, "software_versions.tsv"))

final_status <- file.path(run_root, "PHASE8_SPATIAL_ANNOTATION_STATUS")
if (file.exists(final_status)) stop("Refusing existing final status")
writeLines(c(
  "status=PASS",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("datasets=", length(datasets)),
  paste0("run_id=", run_id),
  "phase7_labels_preserved=TRUE",
  "malignancy_claims=0",
  "reference_states=14"
), paste0(final_status, ".part"), useBytes = TRUE)
if (!file.rename(paste0(final_status, ".part"), final_status)) stop("Final status rename failed")
