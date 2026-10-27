suppressPackageStartupMessages({
  library(data.table)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 00_extract_atlas_umap_v2.R <project_root> <output_dir>")
project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
part_dir <- paste0(output_dir, ".part")
if (dir.exists(output_dir) || dir.exists(part_dir)) stop("Output or .part already exists")
dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)

objects_root <- file.path(
  project_root, "phase8_spatial_alignment_local_v1", "annotation_spatial_crosswalk_v1",
  "runs", "20260811_152500_full_v1", "objects_annotated"
)
object_paths <- c(
  GSE207493 = file.path(objects_root, "GSE207493_spatial_annotated_v1.rds"),
  GSE264586 = file.path(objects_root, "GSE264586_spatial_annotated_v1.rds"),
  GSE294109 = file.path(objects_root, "GSE294109_spatial_annotated_v1.rds")
)
expected_states <- c("Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage")

choose_numeric_column <- function(meta, pattern, preferred = character()) {
  candidates <- grep(pattern, names(meta), value = TRUE, ignore.case = TRUE)
  candidates <- unique(c(preferred[preferred %in% candidates], candidates))
  candidates <- candidates[vapply(meta[, candidates, drop = FALSE], is.numeric, logical(1))]
  if (!length(candidates)) return(NA_character_)
  candidates[[1]]
}

set.seed(26210)
status_rows <- list()
for (dataset in names(object_paths)) {
  path <- object_paths[[dataset]]
  if (!file.exists(path)) stop("Missing object: ", path)
  object <- readRDS(path)
  if (!inherits(object, "Seurat") || !"umap" %in% Reductions(object)) stop("Invalid Seurat/UMAP: ", dataset)
  meta <- object[[]]
  umap <- Embeddings(object, reduction = "umap")[rownames(meta), , drop = FALSE]
  broad_candidates <- grep("phase7.*broad", names(meta), value = TRUE, ignore.case = TRUE)
  fine_candidates <- grep("phase7.*fine", names(meta), value = TRUE, ignore.case = TRUE)
  state_candidates <- grep("phase8.*state", names(meta), value = TRUE, ignore.case = TRUE)
  broad_col <- if ("phase7_broad_annotation" %in% names(meta)) "phase7_broad_annotation" else broad_candidates[[1]]
  fine_col <- if ("phase7_fine_annotation" %in% names(meta)) "phase7_fine_annotation" else fine_candidates[[1]]
  state_col <- if ("phase8_spatial_reference_state" %in% names(meta)) "phase8_spatial_reference_state" else state_candidates[[1]]
  f1_col <- choose_numeric_column(meta, "f1.*score|score.*f1", c("phase8_f1_score", "phase8_f1_loading_score"))

  tymp <- rep(NA_real_, nrow(meta))
  if ("TYMP" %in% rownames(object[["RNA"]])) {
    tymp_layer <- LayerData(object[["RNA"]], layer = "data", features = "TYMP")
    tymp <- as.numeric(tymp_layer[1, rownames(meta), drop = TRUE])
  }
  frame <- data.table(
    cell_id = rownames(meta), dataset = dataset,
    patient_id = as.character(meta[["patient_id"]]), sample_id = as.character(meta[["sample_id"]]),
    broad = as.character(meta[[broad_col]]), fine = as.character(meta[[fine_col]]),
    phase8_state = as.character(meta[[state_col]]), UMAP_1 = as.numeric(umap[, 1]), UMAP_2 = as.numeric(umap[, 2]),
    TYMP_log1p = tymp, F1_score = if (is.na(f1_col)) NA_real_ else as.numeric(meta[[f1_col]])
  )
  frame[is.na(broad) | broad == "", broad := "Unresolved"]
  frame[, plot_group := fifelse(phase8_state %chin% expected_states, phase8_state, broad)]
  sampled <- frame[, {
    cap <- if (.BY$plot_group %chin% expected_states) 6000L else 2200L
    if (.N <= cap) .SD else .SD[sample.int(.N, cap)]
  }, by = plot_group]
  sampled[, plot_group := NULL]
  setorder(sampled, broad, phase8_state, cell_id)
  fwrite(sampled, file.path(part_dir, paste0(dataset, "_atlas_umap_sample_v1.tsv.gz")), sep = "\t", quote = FALSE)
  fwrite(data.table(column = names(meta)), file.path(part_dir, paste0(dataset, "_metadata_columns_v1.tsv")), sep = "\t", quote = FALSE)
  status_rows[[dataset]] <- data.table(
    dataset = dataset, source_path = path, source_bytes = file.info(path)$size,
    total_cells = nrow(frame), sampled_cells = nrow(sampled), patients = uniqueN(frame$patient_id), samples = uniqueN(frame$sample_id),
    broad_column = broad_col, fine_column = fine_col, state_column = state_col,
    f1_column = ifelse(is.na(f1_col), "UNAVAILABLE", f1_col), tymp_finite = sum(is.finite(frame$TYMP_log1p)), status = "PASS"
  )
  rm(object, meta, umap, frame, sampled)
  invisible(gc(verbose = FALSE))
}
fwrite(rbindlist(status_rows), file.path(part_dir, "atlas_extraction_status_v1.tsv"), sep = "\t", quote = FALSE)
if (!file.rename(part_dir, output_dir)) stop("Atomic rename failed")
