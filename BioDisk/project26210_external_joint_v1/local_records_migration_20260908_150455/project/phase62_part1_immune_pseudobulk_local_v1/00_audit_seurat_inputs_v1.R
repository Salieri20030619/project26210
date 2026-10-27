#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: 00_audit_seurat_inputs_v1.R <project_root> <output_tsv>')
root <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = '/', mustWork = FALSE)
run <- file.path(root, 'phase6_processing_local_v1', 'runs', '20260810_234538_376', 'objects_processed')
datasets <- c('GSE264586', 'GSE294109')
rows <- list()
for (dataset_current in datasets) {
  path <- file.path(run, paste0(dataset_current, '_processed_clustered.rds'))
  object <- readRDS(path)
  assay <- DefaultAssay(object)
  layers <- Layers(object[[assay]])
  metadata <- colnames(object@meta.data)
  rows[[dataset_current]] <- data.table(
    dataset = dataset_current,
    path = normalizePath(path, winslash = '/', mustWork = TRUE),
    bytes = file.info(path)$size,
    object_class = paste(class(object), collapse = ';'),
    genes = nrow(object),
    cells = ncol(object),
    assays = paste(Assays(object), collapse = ';'),
    default_assay = assay,
    layers = paste(layers, collapse = ';'),
    has_counts = 'counts' %in% layers,
    has_data = 'data' %in% layers,
    metadata_columns = paste(metadata, collapse = ';'),
    has_patient_id = 'patient_id' %in% metadata,
    has_sample_id = 'sample_id' %in% metadata,
    has_phase7_fine = 'phase7_fine_annotation' %in% metadata,
    has_response = 'response' %in% metadata,
    has_timepoint = 'timepoint' %in% metadata
  )
  rm(object)
  invisible(gc())
}
result <- rbindlist(rows)
part <- paste0(output, '.part')
fwrite(result, part, sep = '\t', quote = TRUE, na = 'NA')
if (!file.rename(part, output)) stop('Atomic rename failed')
print(result[, .(dataset, genes, cells, default_assay, layers, has_counts, has_patient_id, has_phase7_fine)])
