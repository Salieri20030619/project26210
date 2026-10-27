options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(SeuratObject)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: review_phase6_clustering_local_v1.R <root> <run_id>")
root <- normalizePath(args[[1]], winslash = "\\", mustWork = TRUE)
run_id <- args[[2]]
run_dir <- file.path(root, "phase6_processing_local_v1", "runs", run_id)
review_dir <- file.path(root, "phase6_clustering_review_local_v1")
table_dir <- file.path(review_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

composition_file <- file.path(table_dir, "cluster_patient_cell_counts.tsv")
summary_file <- file.path(table_dir, "cluster_patient_review_summary.tsv")
software_file <- file.path(table_dir, "clustering_review_software.tsv")
for (path in c(composition_file, summary_file, software_file)) {
  if (file.exists(path)) stop("Refusing to overwrite existing review output: ", path)
}

write_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}

datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
composition <- vector("list", length(datasets))
review <- vector("list", length(datasets))

for (i in seq_along(datasets)) {
  ds <- datasets[[i]]
  message("CLUSTER_REVIEW_DATASET_START=", ds)
  object_file <- file.path(run_dir, "objects_processed", paste0(ds, "_processed_clustered.rds"))
  obj <- readRDS(object_file)
  if (!inherits(obj, "Seurat")) stop(ds, ": processed object is not Seurat")
  md <- obj[[]]
  if (!"patient_id" %in% names(md) || anyNA(md$patient_id) || any(!nzchar(as.character(md$patient_id)))) {
    stop(ds, ": patient_id completeness gate failed")
  }
  cluster_cols <- grep("_snn_res\\.(0\\.4|0\\.8)$", names(md), value = TRUE)
  if (length(cluster_cols) != 2L) stop(ds, ": expected two cluster fields")
  per_resolution <- lapply(cluster_cols, function(cc) {
    resolution <- sub(".*res\\.", "", cc)
    z <- data.table(
      dataset = ds,
      resolution = resolution,
      cluster = as.character(md[[cc]]),
      patient_id = as.character(md$patient_id)
    )[, .(cells = .N), by = .(dataset, resolution, cluster, patient_id)]
    setorder(z, dataset, resolution, cluster, -cells, patient_id)
    z
  })
  z <- rbindlist(per_resolution)
  composition[[i]] <- z
  s <- z[, {
    total <- sum(cells)
    top <- which.max(cells)
    probs <- cells / total
    .(
      cells = total,
      patients = .N,
      top_patient = patient_id[[top]],
      top_patient_cells = cells[[top]],
      top_patient_fraction = cells[[top]] / total,
      patient_entropy = -sum(probs * log(probs))
    )
  }, by = .(dataset, resolution, cluster)]
  s[, tiny_cluster_lt100 := cells < 100L]
  s[, very_tiny_cluster_lt50 := cells < 50L]
  s[, patient_dominant_ge90pct := top_patient_fraction >= 0.90]
  s[, review_flag := fifelse(very_tiny_cluster_lt50 & patient_dominant_ge90pct,
                             "VERY_TINY_AND_PATIENT_DOMINANT",
                      fifelse(very_tiny_cluster_lt50, "VERY_TINY_LT50",
                      fifelse(tiny_cluster_lt100 & patient_dominant_ge90pct,
                              "TINY_AND_PATIENT_DOMINANT",
                      fifelse(tiny_cluster_lt100, "TINY_LT100",
                      fifelse(patient_dominant_ge90pct, "PATIENT_DOMINANT_GE90PCT", "NONE")))))]
  s[, cluster_order__ := suppressWarnings(as.integer(cluster))]
  setorder(s, dataset, resolution, cluster_order__, cluster)
  s[, cluster_order__ := NULL]
  review[[i]] <- s
  rm(obj, md, z, s, per_resolution)
  invisible(gc())
  message("CLUSTER_REVIEW_DATASET_PASS=", ds)
}

composition <- rbindlist(composition)
review <- rbindlist(review)
write_atomic(composition, composition_file)
write_atomic(review, summary_file)
write_atomic(data.table(
  component = c("R", "SeuratObject", "data.table"),
  version = c(as.character(getRversion()), as.character(packageVersion("SeuratObject")),
              as.character(packageVersion("data.table")))
), software_file)
message("CLUSTER_REVIEW_QUANTITATIVE=PASS")
