options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 32 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: phase7_discover_markers_local_v1.R <root> <processing_run_id> <dataset>")
root <- normalizePath(args[[1]], winslash = "\\", mustWork = TRUE)
processing_run_id <- args[[2]]
dataset <- args[[3]]
expected <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
if (!dataset %in% expected) stop("Unexpected dataset: ", dataset)

input_file <- file.path(root, "phase6_processing_local_v1", "runs", processing_run_id,
                        "objects_processed", paste0(dataset, "_processed_clustered.rds"))
panel_file <- file.path(root, "phase7_annotation_local_v1", "tables", "canonical_marker_panel.tsv")
run_dir <- file.path(root, "phase7_annotation_local_v1", "runs", processing_run_id)
table_dir <- file.path(run_dir, "tables")
figure_dir <- file.path(run_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  markers = file.path(table_dir, paste0(dataset, "_cluster_markers.tsv")),
  top = file.path(table_dir, paste0(dataset, "_cluster_top_markers.tsv")),
  canonical = file.path(table_dir, paste0(dataset, "_canonical_marker_expression.tsv")),
  scores = file.path(table_dir, paste0(dataset, "_broad_lineage_scores.tsv")),
  candidates = file.path(table_dir, paste0(dataset, "_broad_lineage_candidates.tsv")),
  sampling = file.path(table_dir, paste0(dataset, "_marker_sampling.tsv")),
  software = file.path(table_dir, paste0(dataset, "_marker_software.tsv")),
  figure = file.path(figure_dir, paste0(dataset, "_canonical_marker_dotplot.png")),
  status = file.path(run_dir, paste0(dataset, "_MARKER_DISCOVERY_STATUS"))
)
for (p in unlist(paths, use.names = FALSE)) if (file.exists(p)) stop("Refusing to overwrite: ", p)

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

set.seed(26210L)
input_before <- file.info(input_file)
obj <- readRDS(input_file)
if (!inherits(obj, "Seurat")) stop(dataset, ": input is not Seurat")
cluster_col <- grep("_snn_res\\.0\\.8$", colnames(obj[[]]), value = TRUE)
if (length(cluster_col) != 1L) stop(dataset, ": expected one resolution-0.8 cluster field")
cluster <- as.character(obj[[cluster_col, drop = TRUE]])
clusters <- sort(unique(cluster), method = "radix")
Idents(obj) <- factor(cluster, levels = clusters)

panel <- fread(panel_file, sep = "\t", encoding = "UTF-8")
panel <- panel[marker_role == "positive"]
panel_genes <- unique(panel$gene)
present <- intersect(panel_genes, rownames(obj))
missing <- setdiff(panel_genes, present)
if (length(present) < 80L) stop(dataset, ": too few marker panel genes present")

membership <- sparse.model.matrix(~ 0 + factor(cluster, levels = clusters))
colnames(membership) <- clusters
cluster_n <- as.numeric(colSums(membership))
names(cluster_n) <- clusters
data_mat <- LayerData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]
counts_mat <- LayerData(obj, assay = "RNA", layer = "counts")[present, , drop = FALSE]
det_mat <- counts_mat
det_mat@x[] <- 1
avg <- as.matrix(data_mat %*% membership)
avg <- sweep(avg, 2, cluster_n, "/")
pct <- as.matrix(det_mat %*% membership)
pct <- sweep(pct, 2, cluster_n, "/")

avg_long <- as.data.table(as.table(avg))
setnames(avg_long, c("gene", "cluster", "avg_log_normalized"))
pct_long <- as.data.table(as.table(pct))
setnames(pct_long, c("gene", "cluster", "pct_detected"))
canonical <- merge(avg_long, pct_long, by = c("gene", "cluster"), all = TRUE)
canonical <- merge(canonical, panel, by = "gene", allow.cartesian = TRUE)
canonical[, dataset := dataset]
canonical[, cluster_cells := cluster_n[cluster]]
setcolorder(canonical, c("dataset", "cluster", "cluster_cells", "hierarchy", "label", "gene",
                         "avg_log_normalized", "pct_detected", "marker_role", "notes"))
canonical[, cluster_order__ := suppressWarnings(as.integer(as.character(cluster)))]
setorder(canonical, hierarchy, label, gene, cluster_order__, cluster)
canonical[, cluster_order__ := NULL]
write_atomic(canonical, paths$canonical)

broad <- canonical[hierarchy == "broad"]
broad[, gene_z := {
  s <- sd(avg_log_normalized)
  if (!is.finite(s) || s == 0) rep(0, .N) else (avg_log_normalized - mean(avg_log_normalized)) / s
}, by = gene]
broad[, gene_evidence := gene_z * sqrt(pct_detected)]
scores <- broad[, .(
  score = mean(gene_evidence),
  mean_expression = mean(avg_log_normalized),
  mean_pct_detected = mean(pct_detected),
  markers_available = uniqueN(gene),
  markers_detected_ge10pct = uniqueN(gene[pct_detected >= 0.10])
), by = .(dataset, cluster, label)]
setorder(scores, cluster, -score, label)
scores[, rank := seq_len(.N), by = .(dataset, cluster)]
write_atomic(scores, paths$scores)

candidates <- scores[, {
  first <- .SD[1L]
  second <- .SD[min(2L, .N)]
  margin <- first$score - second$score
  decision <- if (first$markers_detected_ge10pct < 2L || !is.finite(first$score) || first$score <= 0) {
    "Unresolved"
  } else if (margin < 0.15) {
    "Mixed_lineage_candidate"
  } else {
    as.character(first$label)
  }
  .(candidate_label = decision,
    top_label = as.character(first$label), top_score = first$score,
    second_label = as.character(second$label), second_score = second$score,
    score_margin = margin,
    top_markers_detected_ge10pct = first$markers_detected_ge10pct,
    cluster_cells = cluster_n[as.character(.BY$cluster)])
}, by = .(dataset, cluster)]
candidates[, cluster_order__ := suppressWarnings(as.integer(as.character(cluster)))]
setorder(candidates, cluster_order__, cluster)
candidates[, cluster_order__ := NULL]
write_atomic(candidates, paths$candidates)

cells_by_cluster <- split(colnames(obj), cluster)
sampled <- unlist(lapply(names(cells_by_cluster), function(cl) {
  x <- cells_by_cluster[[cl]]
  if (length(x) <= 1500L) x else sample(x, 1500L, replace = FALSE)
}), use.names = FALSE)
sampling <- data.table(dataset = dataset, cluster = names(cells_by_cluster),
                       cells_total = lengths(cells_by_cluster),
                       cells_sampled = pmin(lengths(cells_by_cluster), 1500L), seed = 26210L)
sampling[, cluster_order__ := suppressWarnings(as.integer(as.character(cluster)))]
setorder(sampling, cluster_order__, cluster)
sampling[, cluster_order__ := NULL]
write_atomic(sampling, paths$sampling)

obj_small <- subset(obj, cells = sampled)
Idents(obj_small) <- factor(as.character(obj_small[[cluster_col, drop = TRUE]]), levels = clusters)
markers <- FindAllMarkers(
  obj_small, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
  min.pct = 0.10, logfc.threshold = 0.25, min.diff.pct = 0.05,
  return.thresh = 0.05, random.seed = 26210L, verbose = TRUE
)
markers <- as.data.table(markers, keep.rownames = FALSE)
if (!"gene" %in% names(markers)) stop(dataset, ": FindAllMarkers returned no gene column")
fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))
if (length(fc_col) != 1L) stop(dataset, ": could not identify fold-change column")
markers[, dataset := dataset]
markers[, marker_test_scope := "cell_level_exploratory_ranking_not_patient_inference"]
setcolorder(markers, c("dataset", "cluster", "gene", setdiff(names(markers), c("dataset", "cluster", "gene"))))
markers[, cluster_order__ := suppressWarnings(as.integer(as.character(cluster)))]
setorderv(markers, c("cluster_order__", "cluster", "p_val_adj", fc_col), c(1L, 1L, 1L, -1L))
markers[, cluster_order__ := NULL]
write_atomic(markers, paths$markers)

top <- markers[!grepl("^(MT-|RPL|RPS)", gene)]
setorderv(top, c("dataset", "cluster", "p_val_adj", fc_col), c(1L, 1L, 1L, -1L))
top <- top[, head(.SD, 50L), by = .(dataset, cluster)]
top[, marker_rank := seq_len(.N), by = .(dataset, cluster)]
setcolorder(top, c("dataset", "cluster", "marker_rank", "gene", setdiff(names(top), c("dataset", "cluster", "marker_rank", "gene"))))
write_atomic(top, paths$top)

plot_dt <- copy(broad)
plot_dt[, cluster := factor(cluster, levels = rev(clusters))]
plot_dt[, gene := factor(gene, levels = unique(panel[hierarchy == "broad", gene]))]
p <- ggplot(plot_dt, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_detected, color = gene_z)) +
  scale_size(range = c(0, 5), limits = c(0, 1)) +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
  facet_grid(. ~ label, scales = "free_x", space = "free_x") +
  labs(title = paste(dataset, "canonical marker evidence by cluster"),
       subtitle = "Color: within-gene z score of mean log-normalized expression; size: detection fraction",
       x = NULL, y = "resolution 0.8 cluster", color = "gene z", size = "detected") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = 8), panel.grid.minor = element_blank())
plot_part <- paste0(paths$figure, ".part.png")
ggsave(plot_part, p, device = "png", width = 24, height = max(9, length(clusters) * 0.32), dpi = 160, limitsize = FALSE)
if (!file.rename(plot_part, paths$figure)) stop("Atomic figure rename failed")

write_atomic(data.table(
  component = c("R", "Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2"),
  version = c(as.character(getRversion()), as.character(packageVersion("Seurat")),
              as.character(packageVersion("SeuratObject")), as.character(packageVersion("Matrix")),
              as.character(packageVersion("data.table")), as.character(packageVersion("ggplot2")))
), paths$software)

input_after <- file.info(input_file)
input_unchanged <- identical(input_before$size, input_after$size) &&
  identical(as.numeric(input_before$mtime), as.numeric(input_after$mtime))
if (!input_unchanged) stop(dataset, ": Phase 6 input changed")
write_lines_atomic(c(
  "status=PASS",
  paste0("dataset=", dataset),
  paste0("processing_run_id=", processing_run_id),
  paste0("clusters=", length(clusters)),
  paste0("cells_total=", ncol(obj)),
  paste0("cells_marker_sample=", length(sampled)),
  paste0("panel_genes_present=", length(present)),
  paste0("panel_genes_missing=", length(missing)),
  paste0("marker_rows=", nrow(markers)),
  paste0("top_marker_rows=", nrow(top)),
  paste0("source_rds_modified=", !input_unchanged),
  "final_annotation_performed=0",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), paths$status)
message("PHASE7_MARKER_DISCOVERY_PASS=", dataset)
