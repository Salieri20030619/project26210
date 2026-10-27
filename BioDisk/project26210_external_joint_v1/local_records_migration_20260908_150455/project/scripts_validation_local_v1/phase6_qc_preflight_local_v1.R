options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(SeuratObject); library(Matrix); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else getwd(), winslash = "\\", mustWork = TRUE)
out <- file.path(root, "phase6_processing_local_v1", "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
rows <- list(); object_rows <- list()
qprobs <- c(0, .01, .05, .5, .95, .99, 1)
for (ds in datasets) {
  message("QC preflight: ", ds)
  f <- file.path(root, "phase5_merge_by_dataset_local_v2", "objects_raw", paste0(ds, "_merged_raw.rds"))
  obj <- readRDS(f); md <- obj[[]]
  mat <- LayerData(obj, assay = DefaultAssay(obj), layer = "counts")
  mt <- grep("^MT-", rownames(mat))
  pct <- if (length(mt)) 100 * Matrix::colSums(mat[mt, , drop = FALSE]) / pmax(1, Matrix::colSums(mat)) else rep(0, ncol(mat))
  nfeature <- as.numeric(md$nFeature_RNA); ncount <- as.numeric(md$nCount_RNA)
  pass <- nfeature >= 200 & nfeature <= 8000 & ncount >= 500 & ncount <= 100000 & pct <= 25
  metrics <- list(nFeature_RNA = nfeature, nCount_RNA = ncount, percent_mt = pct)
  rows[[ds]] <- rbindlist(lapply(names(metrics), function(nm) {
    q <- quantile(metrics[[nm]], probs = qprobs, na.rm = TRUE, names = FALSE)
    data.table(dataset = ds, metric = nm, statistic = c("min", "p01", "p05", "median", "p95", "p99", "max"), value = q)
  }))
  tmp <- data.table(object_name = as.character(md$object_name), pass = pass)
  object_rows[[ds]] <- tmp[, .(cells = .N, cells_pass = sum(pass), cells_fail = sum(!pass), pass_fraction = mean(pass)), by = object_name][, dataset := ds]
  rm(obj, md, mat, pct, tmp); invisible(gc())
}
summary <- rbindlist(rows); by_object <- rbindlist(object_rows, fill = TRUE)
by_dataset <- by_object[, .(objects = .N, cells = sum(cells), cells_pass = sum(cells_pass), cells_fail = sum(cells_fail), pass_fraction = sum(cells_pass)/sum(cells)), by = dataset]
for (z in list(list(summary, "qc_distribution_preflight.tsv"), list(by_object, "qc_gate_by_object_preflight.tsv"), list(by_dataset, "qc_gate_by_dataset_preflight.tsv"))) {
  p <- file.path(out, z[[2]]); part <- paste0(p, ".part"); fwrite(z[[1]], part, sep = "\t", quote = TRUE); if (file.exists(p)) stop("Refusing overwrite: ", p); if (!file.rename(part, p)) stop("rename failed")
}
cat("PHASE6_QC_PREFLIGHT=PASS\n")
