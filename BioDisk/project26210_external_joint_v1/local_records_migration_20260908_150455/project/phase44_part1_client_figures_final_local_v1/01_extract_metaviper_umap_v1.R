options(stringsAsFactors = FALSE, width = 10000)

root <- "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/08_metaviper_v1/rds"
metadata_path <- file.path(root, "metaviper_cell_metadata_v1.rds")
activity_path <- file.path(root, "metaviper_activity_entrez_by_cell_v1.rds")

meta <- readRDS(metadata_path)
activity <- readRDS(activity_path)

if (is.null(colnames(activity))) stop("activity matrix has no cell names")

pick_col <- function(candidates, label) {
  hit <- candidates[candidates %in% names(meta)]
  if (length(hit) == 0L) stop(sprintf("missing metadata column for %s; columns=%s", label, paste(names(meta), collapse=",")))
  hit[[1L]]
}

patient_col <- pick_col(c("patient", "curated_patient_id", "patient_id"), "patient")
study_col <- pick_col(c("study", "study_v1", "dataset"), "study")

target_ids <- c(CASP9 = "842", CHMP1A = "5119", PRSS8 = "5652", PXN = "5829", HIF1A = "3091")
missing_rows <- setdiff(unname(target_ids), rownames(activity))
if (length(missing_rows) > 0L) stop(sprintf("missing Entrez rows: %s", paste(missing_rows, collapse=",")))

if (!"exact_cell_id" %in% names(meta)) stop("metadata has no exact_cell_id")
meta_index <- match(colnames(activity), meta$exact_cell_id)
common <- colnames(activity)[!is.na(meta_index)]
if (length(common) == 0L) stop("no common cell IDs")
meta <- meta[match(common, meta$exact_cell_id), , drop = FALSE]
activity <- activity[unname(target_ids), common, drop = FALSE]

out <- data.frame(
  exact_cell_id = common,
  study = meta[[study_col]],
  patient = meta[[patient_col]],
  t(activity),
  check.names = FALSE
)
names(out)[(ncol(out) - length(target_ids) + 1L):ncol(out)] <- paste0(names(target_ids), "_activity")

write.table(out, file = "", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
