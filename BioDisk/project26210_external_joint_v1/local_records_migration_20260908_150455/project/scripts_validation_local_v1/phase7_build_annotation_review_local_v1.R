options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: phase7_build_annotation_review_local_v1.R <root> <processing_run_id>")
root <- normalizePath(args[[1]], winslash = "\\", mustWork = TRUE)
run_id <- args[[2]]
datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
run_dir <- file.path(root, "phase7_annotation_local_v1", "runs", run_id)
table_dir <- file.path(run_dir, "tables")
object_dir <- file.path(root, "phase6_processing_local_v1", "runs", run_id, "objects_processed")
out_review <- file.path(table_dir, "cluster_annotation_review_packet.tsv")
out_scores <- file.path(table_dir, "hierarchical_annotation_scores.tsv")
out_meta <- file.path(table_dir, "cluster_metadata_coverage.tsv")
out_original <- file.path(table_dir, "original_annotation_metadata_columns.tsv")
out_validation <- file.path(table_dir, "marker_discovery_validation.tsv")
out_status <- file.path(run_dir, "PHASE7_ANNOTATION_REVIEW_STATUS")
for (p in c(out_review, out_scores, out_meta, out_original, out_validation, out_status)) {
  if (file.exists(p)) stop("Refusing to overwrite: ", p)
}

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

reviews <- list()
score_list <- list()
meta_list <- list()
original_list <- list()
validation_list <- list()

for (dataset in datasets) {
  status_file <- file.path(run_dir, paste0(dataset, "_MARKER_DISCOVERY_STATUS"))
  canonical_file <- file.path(table_dir, paste0(dataset, "_canonical_marker_expression.tsv"))
  candidate_file <- file.path(table_dir, paste0(dataset, "_broad_lineage_candidates.tsv"))
  sampling_file <- file.path(table_dir, paste0(dataset, "_marker_sampling.tsv"))
  marker_file <- file.path(table_dir, paste0(dataset, "_cluster_markers.tsv"))
  top_file <- file.path(table_dir, paste0(dataset, "_cluster_top_markers.tsv"))
  rds_file <- file.path(object_dir, paste0(dataset, "_processed_clustered.rds"))
  required <- c(status_file, canonical_file, candidate_file, sampling_file, marker_file, top_file, rds_file)
  if (!all(file.exists(required))) stop(dataset, ": missing required input")
  status_lines <- readLines(status_file, warn = FALSE)
  if (!"status=PASS" %in% status_lines) stop(dataset, ": marker status is not PASS")

  canonical <- fread(canonical_file, sep = "\t", encoding = "UTF-8")
  candidates <- fread(candidate_file, sep = "\t", encoding = "UTF-8")
  sampling <- fread(sampling_file, sep = "\t", encoding = "UTF-8")
  top <- fread(top_file, sep = "\t", encoding = "UTF-8")
  candidates[, cluster := as.character(cluster)]
  sampling[, cluster := as.character(cluster)]
  candidates[, cluster_cells := NULL]
  candidates <- merge(candidates, sampling[, .(cluster, cluster_cells = cells_total)], by = "cluster", all.x = TRUE)

  canonical[, cluster := as.character(cluster)]
  canonical[, gene_z := {
    s <- sd(avg_log_normalized)
    if (!is.finite(s) || s == 0) rep(0, .N) else (avg_log_normalized - mean(avg_log_normalized)) / s
  }, by = .(hierarchy, gene)]
  canonical[, gene_evidence := gene_z * sqrt(pct_detected)]
  scores <- canonical[, .(
    score = mean(gene_evidence),
    mean_expression = mean(avg_log_normalized),
    mean_pct_detected = mean(pct_detected),
    markers_available = uniqueN(gene),
    markers_detected_ge10pct = uniqueN(gene[pct_detected >= 0.10])
  ), by = .(dataset, cluster, hierarchy, label)]
  setorder(scores, dataset, cluster, hierarchy, -score, label)
  scores[, rank := seq_len(.N), by = .(dataset, cluster, hierarchy)]
  score_list[[dataset]] <- scores

  top[, cluster := as.character(cluster)]
  setorder(top, dataset, cluster, marker_rank)
  marker_summary <- top[marker_rank <= 20L, .(
    top20_markers = paste(gene, collapse = ";"),
    top10_markers = paste(head(gene, 10L), collapse = ";")
  ), by = .(dataset, cluster)]

  hierarchy_summary <- scores[rank <= 2L, {
    first <- .SD[1L]
    second <- .SD[min(2L, .N)]
    .(top_label = as.character(first$label), top_score = first$score,
      second_label = as.character(second$label), second_score = second$score,
      score_margin = first$score - second$score,
      top_markers_detected_ge10pct = first$markers_detected_ge10pct)
  }, by = .(dataset, cluster, hierarchy)]
  hierarchy_wide <- dcast(
    hierarchy_summary, dataset + cluster ~ hierarchy,
    value.var = c("top_label", "top_score", "second_label", "second_score", "score_margin",
                  "top_markers_detected_ge10pct"), sep = "__"
  )

  input_before <- file.info(rds_file)
  obj <- readRDS(rds_file)
  md <- as.data.table(obj[[]], keep.rownames = "cell_id__")
  cluster_col <- grep("_snn_res\\.0\\.8$", names(md), value = TRUE)
  if (length(cluster_col) != 1L) stop(dataset, ": expected one res0.8 cluster column")
  md[, cluster := as.character(get(cluster_col))]
  patient_col <- if ("patient_id" %in% names(md)) "patient_id" else NA_character_
  sample_col <- if ("sample_id" %in% names(md)) "sample_id" else NA_character_
  tissue_col <- if ("tissue_status" %in% names(md)) "tissue_status" else NA_character_
  disease_col <- if ("disease_status" %in% names(md)) "disease_status" else NA_character_
  normalize_text <- function(x) {
    x <- as.character(x); x[is.na(x) | !nzchar(trimws(x))] <- "Unknown"; x
  }
  md[, patient__ := if (is.na(patient_col)) "Unknown" else normalize_text(get(patient_col))]
  md[, sample__ := if (is.na(sample_col)) "Unknown" else normalize_text(get(sample_col))]
  md[, tissue__ := if (is.na(tissue_col)) "Unknown" else normalize_text(get(tissue_col))]
  md[, disease__ := if (is.na(disease_col)) "Unknown" else normalize_text(get(disease_col))]
  pat <- md[, .N, by = .(cluster, patient__)][order(cluster, -N, patient__)]
  pat_sum <- pat[, .(n_patients = uniqueN(patient__), top_patient = patient__[1L],
                     max_patient_fraction = N[1L] / sum(N)), by = cluster]
  samp <- md[, .N, by = .(cluster, sample__)][order(cluster, -N, sample__)]
  samp_sum <- samp[, .(n_samples = uniqueN(sample__), top_sample = sample__[1L],
                       max_sample_fraction = N[1L] / sum(N)), by = cluster]
  tissue <- md[, .N, by = .(cluster, tissue__)][order(cluster, -N, tissue__)]
  tissue_sum <- tissue[, .(tissue_composition = paste0(tissue__, ":", N, collapse = ";")), by = cluster]
  disease <- md[, .N, by = .(cluster, disease__)][order(cluster, -N, disease__)]
  disease_sum <- disease[, .(disease_composition = paste0(disease__, ":", N, collapse = ";")), by = cluster]
  meta <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE),
                 list(md[, .(cluster_cells_metadata = .N), by = cluster], pat_sum, samp_sum, tissue_sum, disease_sum))
  meta[, dataset := dataset]
  setcolorder(meta, c("dataset", "cluster", setdiff(names(meta), c("dataset", "cluster"))))
  meta_list[[dataset]] <- meta

  excluded <- c("seurat_clusters", cluster_col)
  annotation_cols <- setdiff(grep("cell.?type|annot|lineage|identity|cell.?class", names(md),
                                  value = TRUE, ignore.case = TRUE), excluded)
  if (length(annotation_cols)) {
    original_list[[dataset]] <- rbindlist(lapply(annotation_cols, function(nm) data.table(
      dataset = dataset, metadata_column = nm,
      nonmissing_cells = sum(!is.na(md[[nm]]) & nzchar(trimws(as.character(md[[nm]])))),
      unique_nonmissing_values = uniqueN(as.character(md[[nm]])[!is.na(md[[nm]]) & nzchar(trimws(as.character(md[[nm]])))])
    )), fill = TRUE)
  } else {
    original_list[[dataset]] <- data.table(dataset = dataset, metadata_column = NA_character_,
                                           nonmissing_cells = 0L, unique_nonmissing_values = 0L)
  }
  input_after <- file.info(rds_file)
  unchanged <- identical(input_before$size, input_after$size) &&
    identical(as.numeric(input_before$mtime), as.numeric(input_after$mtime))
  if (!unchanged) stop(dataset, ": source RDS changed")
  rm(obj, md); gc(verbose = FALSE)

  review <- Reduce(function(x, y) merge(x, y, by = c("dataset", "cluster"), all = TRUE),
                   list(candidates, marker_summary, hierarchy_wide, meta))
  review[, cluster_order__ := suppressWarnings(as.integer(cluster))]
  setorder(review, dataset, cluster_order__, cluster)
  review[, cluster_order__ := NULL]
  reviews[[dataset]] <- review

  marker_rows <- nrow(fread(marker_file, sep = "\t", select = "dataset"))
  validation_list[[dataset]] <- data.table(
    dataset = dataset, status_pass = TRUE, clusters = uniqueN(sampling$cluster),
    marker_rows = marker_rows, top_marker_rows = nrow(top),
    canonical_rows = nrow(canonical), candidate_rows = nrow(candidates),
    cluster_cell_counts_complete = all(is.finite(candidates$cluster_cells)),
    source_rds_unchanged = unchanged
  )
  rm(canonical, candidates, sampling, top, scores, meta, review); gc(verbose = FALSE)
}

review_all <- rbindlist(reviews, fill = TRUE)
scores_all <- rbindlist(score_list, fill = TRUE)
meta_all <- rbindlist(meta_list, fill = TRUE)
original_all <- rbindlist(original_list, fill = TRUE)
validation_all <- rbindlist(validation_list, fill = TRUE)
if (nrow(review_all) != 178L || uniqueN(review_all[, paste(dataset, cluster)]) != 178L) {
  stop("Expected exactly 178 unique dataset-clusters")
}
if (!all(validation_all$status_pass & validation_all$cluster_cell_counts_complete & validation_all$source_rds_unchanged)) {
  stop("Marker validation invariant failed")
}
write_atomic(review_all, out_review)
write_atomic(scores_all, out_scores)
write_atomic(meta_all, out_meta)
write_atomic(original_all, out_original)
write_atomic(validation_all, out_validation)
write_lines_atomic(c(
  "status=PASS",
  paste0("processing_run_id=", run_id),
  paste0("datasets=", length(datasets)),
  paste0("clusters_reviewed=", nrow(review_all)),
  paste0("marker_rows=", sum(validation_all$marker_rows)),
  paste0("top_marker_rows=", sum(validation_all$top_marker_rows)),
  paste0("source_rds_modified=", sum(!validation_all$source_rds_unchanged)),
  "automatic_candidates_are_final_annotations=FALSE",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), out_status)
message("PHASE7_ANNOTATION_REVIEW_PASS")
