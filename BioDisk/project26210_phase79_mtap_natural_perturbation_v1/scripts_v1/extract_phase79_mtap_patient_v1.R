#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5L) {
  stop("Usage: extract_phase79_mtap_patient_v1.R dataset seurat_rds metadata patient_summary output")
}
dataset <- args[[1L]]
rds_path <- normalizePath(args[[2L]], winslash="/", mustWork=TRUE)
metadata_path <- normalizePath(args[[3L]], winslash="/", mustWork=TRUE)
patient_path <- normalizePath(args[[4L]], winslash="/", mustWork=TRUE)
output_path <- args[[5L]]
if (file.exists(output_path) || file.exists(paste0(output_path, ".part"))) stop("Refusing existing output")

meta <- fread(metadata_path)
if (dataset == "GSE207493") {
  if (anyDuplicated(meta$cell_barcode)) stop("Duplicate metadata ID")
  meta[, `:=`(
    cell_id=cell_barcode,
    seurat_cell_name=paste0("GSE207493__GSE207493__", sub("^([^_]+)_", "\\1__", cell_barcode)),
    patient_key=paste0("GSE207493::", sample),
    timepoint="Single",
    is_tumor=cell_type == "Tumor"
  )]
} else {
  if (anyDuplicated(meta$cell_id) || anyDuplicated(meta$seurat_cell_name)) stop("Duplicate metadata ID")
  meta[, is_tumor := phase7_fine_annotation == "Epithelial_tumor_like"]
}
meta[, unit_id := paste(patient_key, timepoint, sep="::")]

object <- readRDS(rds_path)
if (!inherits(object, "Seurat") || !("RNA" %chin% Assays(object))) stop("Invalid Seurat object")
eall <- LayerData(object[["RNA"]], layer="data")
missing <- setdiff(meta$seurat_cell_name, colnames(eall))
if (length(missing)) stop("Metadata/Seurat mismatch: ", length(missing))
source_names <- rownames(eall)
upper <- toupper(source_names)
collisions <- unique(upper[duplicated(upper) | duplicated(upper, fromLast=TRUE)])
if ("MTAP" %chin% collisions || !("MTAP" %chin% upper)) stop("MTAP unavailable or ambiguous")
mtap_row <- which(upper == "MTAP")[[1L]]
tumor_meta <- meta[is_tumor == TRUE]
mtap <- as.numeric(eall[mtap_row, tumor_meta$seurat_cell_name, drop=TRUE])
rm(object, eall)
invisible(gc())

cells <- data.table(unit_id=tumor_meta$unit_id, patient_key=tumor_meta$patient_key, mtap=mtap)
units <- cells[, .(
  tumor_mtap_mean=mean(mtap),
  tumor_mtap_detect_fraction=mean(mtap > 0),
  n_tumor_cells_mtap=.N
), by=.(patient_key, unit_id)]
patients_mtap <- units[, .(
  tumor_mtap_mean=median(tumor_mtap_mean),
  tumor_mtap_detect_fraction=median(tumor_mtap_detect_fraction),
  n_tumor_cells_mtap=sum(n_tumor_cells_mtap)
), by=patient_key]

patient <- fread(patient_path)
merged <- merge(patient, patients_mtap, by="patient_key", all.x=TRUE, sort=FALSE)
if (nrow(merged) != nrow(patient) || anyNA(merged$tumor_mtap_mean)) stop("Patient MTAP merge failed")
merged[, mtap_source_dataset := dataset]
temp <- paste0(output_path, ".part")
fwrite(merged, temp, sep="\t", quote=TRUE, na="NA")
if (!file.rename(temp, output_path)) stop("Atomic rename failed")
