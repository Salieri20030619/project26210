#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: 00_export_phase67_metadata_v1.R <compact_run_root> <output_dir>')
compact_root <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = '/', mustWork = FALSE)
part_dir <- paste0(output_dir, '.part')
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop('Refusing existing output')
dir.create(part_dir, recursive = TRUE)

sha256 <- function(path) unname(tools::md5sum(path))
atomic_fwrite <- function(x, path) {
  part <- if (grepl('\\.gz$', path)) sub('\\.gz$', '.part.gz', path) else paste0(path, '.part')
  fwrite(x, part, sep = '\t', quote = TRUE, na = 'NA')
  if (!file.rename(part, path)) stop('Atomic rename failed: ', path)
}

datasets <- c('GSE264586', 'GSE294109')
columns <- c(
  'dataset', 'cell_id', 'seurat_cell_name', 'patient_key', 'patient_id',
  'sample_id', 'object_name', 'timepoint', 'treatment', 'response',
  'phase7_broad_annotation', 'phase7_fine_annotation'
)
audits <- list()
for (dataset in datasets) {
  path <- file.path(compact_root, dataset, 'objects', paste0(dataset, '_cell_scores_compact_v1.rds'))
  if (!file.exists(path)) stop('Missing compact RDS: ', path)
  value <- as.data.table(readRDS(path))
  if (!all(columns %in% names(value))) stop('Missing metadata columns: ', dataset)
  selected <- copy(value[, ..columns])
  for (column in columns) set(selected, j = column, value = as.character(selected[[column]]))
  if (anyDuplicated(selected$cell_id) || anyNA(selected$cell_id)) stop('Invalid cell IDs: ', dataset)
  selected[, lineage_v1 := fifelse(
    phase7_broad_annotation == 'Myeloid', 'Myeloid_DC',
    fifelse(phase7_fine_annotation == 'Treg', 'Treg', 'Other')
  )]
  atomic_fwrite(selected, file.path(part_dir, paste0(dataset, '_phase67_metadata_v1.tsv.gz')))
  audits[[dataset]] <- selected[, .(
    n_cells = .N,
    n_patients = uniqueN(patient_key),
    n_patient_timepoints = uniqueN(paste(patient_key, timepoint))
  ), by = .(dataset, lineage_v1, phase7_fine_annotation)]
  rm(value, selected)
  invisible(gc())
}

audit <- rbindlist(audits)
setorder(audit, dataset, lineage_v1, phase7_fine_annotation)
atomic_fwrite(audit, file.path(part_dir, 'phase67_metadata_audit_v1.tsv'))
status <- data.table(
  status = 'PASS_PHASE67_METADATA_EXPORT',
  datasets = length(datasets),
  total_cells = sum(audit$n_cells),
  myeloid_dc_cells = audit[lineage_v1 == 'Myeloid_DC', sum(n_cells)],
  treg_cells = audit[lineage_v1 == 'Treg', sum(n_cells)]
)
atomic_fwrite(status, file.path(part_dir, 'PHASE67_METADATA_STATUS_v1.tsv'))
if (!file.rename(part_dir, output_dir)) stop('Final rename failed')
print(status)
