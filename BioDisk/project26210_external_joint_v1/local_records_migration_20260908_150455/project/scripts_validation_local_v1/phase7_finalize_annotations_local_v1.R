options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 32 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: phase7_finalize_annotations_local_v1.R <root> <processing_run_id> <output_version>")
dry_run <- identical(Sys.getenv("PHASE7_DRY_RUN", unset = "0"), "1")
root <- normalizePath(args[[1]], winslash = "\\", mustWork = TRUE)
run_id <- args[[2]]
output_version <- args[[3]]
if (!grepl("^local_v[0-9]+$", output_version)) stop("Invalid output_version: ", output_version)
annotation_version <- paste0("phase7_annotation_", output_version)
datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
input_run_dir <- file.path(root, "phase7_annotation_local_v1", "runs", run_id)
run_dir <- file.path(root, annotation_version, "runs", run_id)
table_dir <- file.path(run_dir, "tables")
figure_dir <- file.path(run_dir, "figures_annotation")
object_dir <- file.path(run_dir, "objects_annotated")
source_dir <- file.path(root, "phase6_processing_local_v1", "runs", run_id, "objects_processed")
review_file <- file.path(input_run_dir, "tables", "cluster_annotation_review_packet.tsv")
override_file <- file.path(root, "phase7_annotation_local_v1", "tables", "manual_cluster_annotation_overrides.tsv")
mapping_file <- file.path(table_dir, "final_cluster_annotations.tsv")
summary_file <- file.path(table_dir, "final_annotation_summary.tsv")
state_summary_file <- file.path(table_dir, "final_annotation_state_summary.tsv")
confidence_summary_file <- file.path(table_dir, "final_annotation_confidence_summary.tsv")
software_file <- file.path(table_dir, "annotation_software.tsv")
status_file <- file.path(run_dir, "PHASE7_ANNOTATION_STATUS")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

object_paths <- file.path(object_dir, paste0(datasets, "_annotated_", output_version, ".rds"))
dataset_status_paths <- file.path(run_dir, paste0(datasets, "_ANNOTATION_STATUS"))
figure_paths <- file.path(figure_dir, paste0(datasets, "_annotation_umap.png"))
all_outputs <- c(mapping_file, summary_file, state_summary_file, confidence_summary_file,
                 software_file, status_file, object_paths, dataset_status_paths, figure_paths)
for (p in all_outputs) if (file.exists(p)) stop("Refusing to overwrite: ", p)
if (!all(file.exists(c(review_file, override_file)))) stop("Review packet or override table missing")

write_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}
write_lines_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  writeLines(x, part, useBytes = TRUE)
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

review <- fread(review_file, sep = "\t", encoding = "UTF-8")
overrides <- fread(override_file, sep = "\t", encoding = "UTF-8")
review[, cluster := as.character(cluster)]
overrides[, cluster := as.character(cluster)]
if (nrow(review) != 178L || uniqueN(paste(review$dataset, review$cluster, sep = "::")) != 178L) stop("Review packet key failure")
if (anyDuplicated(overrides, by = c("dataset", "cluster"))) stop("Duplicate override keys")
if (nrow(overrides[!review, on = .(dataset, cluster)])) stop("Override key absent from review packet")

ann <- copy(review)
ann[, broad_annotation := fifelse(candidate_label %chin% c("Mixed_lineage_candidate", "Unresolved"),
                                  "Unresolved", candidate_label)]
ann[, fine_annotation := broad_annotation]
ann[broad_annotation == "Epithelial", fine_annotation := fifelse(
  as.numeric(top_score__epithelial) >= 0.5 & as.numeric(score_margin__epithelial) >= 0.2 &
    as.integer(top_markers_detected_ge10pct__epithelial) >= 1L,
  top_label__epithelial, "Epithelial_unresolved")]
ann[broad_annotation == "Myeloid", fine_annotation := fifelse(
  as.numeric(top_score__myeloid) >= 0.5 & as.numeric(score_margin__myeloid) >= 0.15 &
    as.integer(top_markers_detected_ge10pct__myeloid) >= 2L,
  top_label__myeloid, "Myeloid_unresolved")]
ann[broad_annotation == "T_cell", fine_annotation := fifelse(
  as.numeric(top_score__t_nk) >= 0.5 & as.numeric(score_margin__t_nk) >= 0.15 &
    as.integer(top_markers_detected_ge10pct__t_nk) >= 2L,
  top_label__t_nk, "T_cell_unresolved")]
ann[, annotation_state := fifelse(as.numeric(top_score__state) >= 2 &
                                    as.integer(top_markers_detected_ge10pct__state) == 5L,
                                  "Cycling", "None")]
ann[, confidence := "moderate"]
ann[cluster_cells < 50L | max_patient_fraction >= 0.98 |
      broad_annotation == "Unresolved" | grepl("_unresolved$", fine_annotation, ignore.case = TRUE) |
      as.numeric(score_margin__broad) < 0.15, confidence := "low"]
ann[confidence == "moderate" & cluster_cells >= 100L & max_patient_fraction <= 0.90 &
      as.numeric(score_margin__broad) >= 0.50 &
      as.integer(top_markers_detected_ge10pct__broad) >= 3L, confidence := "high"]
ann[, manual_override := FALSE]
ann[, review_note := "Hierarchical canonical-panel rule with independent cluster-marker review packet."]

setkey(ann, dataset, cluster)
setkey(overrides, dataset, cluster)
idx <- ann[overrides, which = TRUE]
ann[idx, `:=`(
  broad_annotation = overrides$broad_annotation,
  fine_annotation = overrides$fine_annotation,
  annotation_state = overrides$annotation_state,
  confidence = overrides$confidence,
  manual_override = TRUE,
  review_note = overrides$review_note
)]
ann[, annotation_method := "independent_cluster_markers_plus_canonical_panel_with_manual_audit"]
ann[, malignancy_claim := FALSE]
ann[, original_author_labels_copied := FALSE]
ann[, evidence_summary := paste0(
  "top10=", top10_markers,
  " | broad=", top_label__broad, ":", sprintf("%.3f", as.numeric(top_score__broad)),
  ";second=", second_label__broad, ":", sprintf("%.3f", as.numeric(second_score__broad)),
  " | patients=", n_patients, ";max_patient_fraction=", sprintf("%.3f", as.numeric(max_patient_fraction))
)]
ann[, cluster_order__ := suppressWarnings(as.integer(cluster))]
setorder(ann, dataset, cluster_order__, cluster)
ann[, cluster_order__ := NULL]

mapping_cols <- c(
  "dataset", "cluster", "cluster_cells", "broad_annotation", "fine_annotation", "annotation_state",
  "confidence", "manual_override", "review_note", "annotation_method", "malignancy_claim",
  "original_author_labels_copied", "evidence_summary", "top20_markers", "candidate_label",
  "top_label__broad", "top_score__broad", "second_label__broad", "second_score__broad",
  "score_margin__broad", "n_patients", "max_patient_fraction", "n_samples", "max_sample_fraction",
  "tissue_composition", "disease_composition"
)
mapping <- ann[, ..mapping_cols]
if (anyNA(mapping[, .(dataset, cluster, broad_annotation, fine_annotation, annotation_state, confidence)])) {
  stop("Final mapping contains missing required values")
}
if (!all(mapping$confidence %chin% c("high", "moderate", "low"))) stop("Invalid confidence value")
if (dry_run) {
  message("PHASE7_ANNOTATION_DRY_RUN_PASS clusters=", nrow(mapping),
          " cells=", sum(mapping$cluster_cells),
          " overrides=", sum(mapping$manual_override))
  quit(save = "no", status = 0L, runLast = FALSE)
}
write_atomic(mapping, mapping_file)
write_atomic(mapping[, .(clusters = .N, cells = sum(cluster_cells)),
                     by = .(dataset, broad_annotation, fine_annotation)], summary_file)
write_atomic(mapping[, .(clusters = .N, cells = sum(cluster_cells)),
                     by = .(dataset, annotation_state)], state_summary_file)
write_atomic(mapping[, .(clusters = .N, cells = sum(cluster_cells)),
                     by = .(dataset, confidence)], confidence_summary_file)
write_atomic(data.table(
  component = c("R", "Seurat", "SeuratObject", "data.table", "ggplot2"),
  version = c(as.character(getRversion()), as.character(packageVersion("Seurat")),
              as.character(packageVersion("SeuratObject")), as.character(packageVersion("data.table")),
              as.character(packageVersion("ggplot2")))
), software_file)

validated <- logical(length(datasets))
names(validated) <- datasets
for (dataset_current in datasets) {
  dataset <- dataset_current
  source_file <- file.path(source_dir, paste0(dataset, "_processed_clustered.rds"))
  output_file <- file.path(object_dir, paste0(dataset, "_annotated_", output_version, ".rds"))
  dataset_status <- file.path(run_dir, paste0(dataset, "_ANNOTATION_STATUS"))
  figure_file <- file.path(figure_dir, paste0(dataset, "_annotation_umap.png"))
  source_before <- file.info(source_file)
  obj <- readRDS(source_file)
  cluster_col <- grep("_snn_res\\.0\\.8$", colnames(obj[[]]), value = TRUE)
  if (length(cluster_col) != 1L) stop(dataset, ": expected one res0.8 cluster field")
  new_cols <- c("phase7_cluster_res0_8", "phase7_broad_annotation", "phase7_fine_annotation",
                "phase7_annotation_state", "phase7_annotation_confidence", "phase7_annotation_version")
  if (any(new_cols %in% colnames(obj[[]]))) stop(dataset, ": refusing metadata collision")
  cluster_value <- as.character(obj[[cluster_col, drop = TRUE]])
  map <- mapping[mapping$dataset == dataset_current]
  match_idx <- match(cluster_value, map$cluster)
  if (anyNA(match_idx)) stop(dataset, ": unmapped cells")
  obj$phase7_cluster_res0_8 <- cluster_value
  obj$phase7_broad_annotation <- factor(map$broad_annotation[match_idx], levels = unique(map$broad_annotation))
  obj$phase7_fine_annotation <- factor(map$fine_annotation[match_idx], levels = unique(map$fine_annotation))
  obj$phase7_annotation_state <- factor(map$annotation_state[match_idx], levels = unique(map$annotation_state))
  obj$phase7_annotation_confidence <- factor(map$confidence[match_idx], levels = c("high", "moderate", "low"))
  obj$phase7_annotation_version <- annotation_version
  if (anyNA(obj$phase7_broad_annotation) || anyNA(obj$phase7_fine_annotation)) stop(dataset, ": NA cell labels")

  p1 <- DimPlot(obj, reduction = "umap", group.by = "phase7_broad_annotation", label = TRUE,
                repel = TRUE, raster = TRUE) + ggtitle(paste(dataset, "broad annotation")) + NoLegend()
  p2 <- DimPlot(obj, reduction = "umap", group.by = "phase7_fine_annotation", label = TRUE,
                repel = TRUE, raster = TRUE) + ggtitle(paste(dataset, "fine annotation")) + NoLegend()
  combined <- patchwork::wrap_plots(p1, p2, ncol = 2)
  figure_part <- paste0(figure_file, ".part.png")
  ggsave(figure_part, combined, device = "png", width = 22, height = 10, dpi = 160, limitsize = FALSE)
  if (!file.rename(figure_part, figure_file)) stop(dataset, ": atomic figure rename failed")

  output_part <- paste0(output_file, ".part")
  saveRDS(obj, output_part, compress = TRUE)
  if (!file.rename(output_part, output_file)) stop(dataset, ": atomic RDS rename failed")
  source_after <- file.info(source_file)
  unchanged <- identical(source_before$size, source_after$size) &&
    identical(as.numeric(source_before$mtime), as.numeric(source_after$mtime))
  if (!unchanged) stop(dataset, ": source RDS changed")
  rm(obj); gc(verbose = FALSE)
  check <- readRDS(output_file)
  validated[[dataset]] <- inherits(check, "Seurat") && ncol(check) == sum(map$cluster_cells) &&
    all(new_cols %in% colnames(check[[]])) && !anyNA(check$phase7_broad_annotation) &&
    !anyNA(check$phase7_fine_annotation)
  if (!validated[[dataset]]) stop(dataset, ": annotated RDS validation failed")
  rm(check); gc(verbose = FALSE)
  write_lines_atomic(c(
    "status=PASS",
    paste0("dataset=", dataset),
    paste0("processing_run_id=", run_id),
    paste0("cells=", sum(map$cluster_cells)),
    paste0("clusters=", nrow(map)),
    paste0("broad_labels=", uniqueN(map$broad_annotation)),
    paste0("fine_labels=", uniqueN(map$fine_annotation)),
    paste0("low_confidence_clusters=", sum(map$confidence == "low")),
    paste0("manual_override_clusters=", sum(map$manual_override)),
    "source_rds_modified=FALSE",
    "annotated_rds_readable=TRUE",
    "original_author_labels_copied=FALSE",
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  ), dataset_status)
}

if (!all(validated)) stop("Not all annotated RDS outputs validated")
write_lines_atomic(c(
  "status=PASS",
  paste0("processing_run_id=", run_id),
  paste0("datasets=", length(datasets)),
  paste0("clusters_annotated=", nrow(mapping)),
  paste0("cells_annotated=", sum(mapping$cluster_cells)),
  paste0("broad_labels=", uniqueN(mapping$broad_annotation)),
  paste0("fine_labels=", uniqueN(mapping$fine_annotation)),
  paste0("low_confidence_clusters=", sum(mapping$confidence == "low")),
  paste0("unresolved_clusters=", sum(grepl("Unresolved", mapping$fine_annotation, ignore.case = TRUE))),
  paste0("manual_override_clusters=", sum(mapping$manual_override)),
  "all_outputs_readRDS=TRUE",
  "source_rds_modified=0",
  "integration_performed=0",
  "malignancy_claims_made=0",
  "original_author_labels_copied=0",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), status_file)
message("PHASE7_ANNOTATION_PASS")
