options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], winslash = "\\", mustWork = TRUE) else normalizePath(getwd(), winslash = "\\", mustWork = TRUE)
run_id <- Sys.getenv("PHASE6A_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

expected_datasets <- c("GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109")
phase5_dir <- file.path(root, "phase5_merge_by_dataset_local_v2")
object_dir <- file.path(phase5_dir, "objects_raw")
clinical_file <- file.path(root, "phase4e_clinical_metadata_harmonization_v1", "tables", "object_clinical_manifest.tsv")
run_dir <- file.path(root, "singlecell_validation_local_v1", "runs", run_id)
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv_atomic <- function(x, path) {
  part <- paste0(path, ".part")
  data.table::fwrite(x, part, sep = "\t", quote = TRUE, na = "NA")
  if (file.exists(path)) stop("Refusing to overwrite existing output: ", path)
  if (!file.rename(part, path)) stop("Atomic rename failed: ", part, " -> ", path)
}

append_progress <- function(dataset, event, details = "") {
  path <- file.path(run_dir, "PHASE6A_PROGRESS.tsv")
  row <- data.table(timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), dataset = dataset, event = event, details = details)
  data.table::fwrite(row, path, sep = "\t", append = file.exists(path), col.names = !file.exists(path), quote = TRUE)
}

is_informative <- function(x) {
  z <- trimws(as.character(x))
  !is.na(z) & nzchar(z) & !toupper(z) %in% c("NA", "N/A", "UNKNOWN", "UNRESOLVED", "NOT FOUND", "MISSING", "NONE PROVIDED")
}

phase5_status <- file.path(phase5_dir, "PHASE5_MERGE_STATUS")
if (!file.exists(phase5_status) || !any(grepl("^status=PASS$", readLines(phase5_status, warn = FALSE)))) {
  stop("Phase 5 PASS status is absent: ", phase5_status)
}
if (!file.exists(clinical_file)) stop("Clinical manifest is absent: ", clinical_file)

marker_registry <- data.table::rbindlist(list(
  data.table(program = "malignant_TYMP_hypoxia_EMT_anchors", gene = c("TYMP", "EPCAM", "KRT8", "KRT18", "CA9", "VEGFA", "HIF1A", "VIM")),
  data.table(program = "myeloid_SPP1_APOE_C1QC_HLAII_anchors", gene = c("SPP1", "APOE", "C1QA", "C1QB", "C1QC", "TREM2", "GPNMB", "FOLR2", "CXCL10", "HLA-DRA", "CD74", "LST1")),
  data.table(program = "vascular_endothelial_anchors", gene = c("PECAM1", "VWF", "EMCN", "KDR", "FLT1", "ENG", "RAMP2")),
  data.table(program = "T_NK_function_anchors", gene = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "PDCD1", "LAG3", "HAVCR2")),
  data.table(program = "fibroblast_ECM_anchors", gene = c("POSTN", "COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1")),
  data.table(program = "signaling_metabolic_proxy_anchors", gene = c("MIF", "CXCR4", "TGFB1", "CSF1", "CSF1R", "MTAP", "MAT2A", "AMD1", "SRM", "SMS"))
))
marker_registry[, interpretation := "Technical feature/detectability anchor only; not a frozen signature and not a confirmatory score"]
write_tsv_atomic(marker_registry, file.path(table_dir, "anchor_gene_registry.tsv"))

signature_candidates <- list.files(root, recursive = TRUE, full.names = TRUE,
                                   pattern = "(signature|geneset|gene_set|module).*(tsv|csv|txt|gmt|rds)$",
                                   ignore.case = TRUE)
signature_candidates <- signature_candidates[!grepl("singlecell_validation_local_v1|local_reference_review_v1", signature_candidates, ignore.case = TRUE)]
signature_inventory <- data.table(path = normalizePath(signature_candidates, winslash = "\\", mustWork = FALSE),
                                  bytes = if (length(signature_candidates)) file.info(signature_candidates)$size else numeric())
if (!nrow(signature_inventory)) signature_inventory <- data.table(path = character(), bytes = numeric())
write_tsv_atomic(signature_inventory, file.path(table_dir, "frozen_signature_file_inventory.tsv"))

clinical <- data.table::fread(clinical_file, sep = "\t", encoding = "UTF-8")
object_audit <- list()
metadata_audit <- list()
unit_counts <- list()
label_inventory <- list()
marker_audit <- list()
input_integrity <- list()

for (ds in expected_datasets) {
  append_progress(ds, "START")
  object_file <- file.path(object_dir, paste0(ds, "_merged_raw.rds"))
  if (!file.exists(object_file)) stop("Missing Phase 5 object: ", object_file)
  fi_before <- file.info(object_file)
  message(sprintf("[%s] reading %s", format(Sys.time(), "%H:%M:%S"), ds))
  obj <- readRDS(object_file)
  if (!inherits(obj, "Seurat")) stop(ds, " is not a Seurat object")

  assays <- SeuratObject::Assays(obj)
  assay <- SeuratObject::DefaultAssay(obj)
  layers <- SeuratObject::Layers(obj[[assay]])
  counts_exact <- identical(layers, "counts")
  if (!"counts" %in% layers) stop(ds, " has no counts layer")
  counts <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  md <- obj[[]]
  ids <- colnames(obj)
  core_fields <- c("cell_id", "barcode_original", "dataset", "object_name", "sample_id", "patient_id", "tissue_status", "disease_status", "treatment", "timepoint", "source_file", "source_assay", "source_count_location")
  missing_core <- setdiff(core_fields, colnames(md))

  candidate_fields <- grep("cell.?type|annotation|lineage|major.?type|minor.?type|celltypist|predicted.?label|cluster", colnames(md), value = TRUE, ignore.case = TRUE)
  candidate_fields <- setdiff(candidate_fields, c("nFeature_RNA", "nCount_RNA"))
  informative_label_fields <- candidate_fields[vapply(candidate_fields, function(f) length(unique(as.character(md[[f]])[is_informative(md[[f]])])) > 1L, logical(1))]

  ident_values <- as.character(SeuratObject::Idents(obj))
  ident_unique <- length(unique(ident_values))
  ident_matches_object <- "object_name" %in% colnames(md) && identical(ident_values, as.character(md$object_name))
  ident_matches_orig <- "orig.ident" %in% colnames(md) && identical(ident_values, as.character(md$orig.ident))

  object_audit[[ds]] <- data.table(
    dataset = ds,
    object_file = normalizePath(object_file, winslash = "\\", mustWork = TRUE),
    file_bytes = fi_before$size,
    object_class = paste(class(obj), collapse = ";"),
    assays = paste(assays, collapse = ";"),
    default_assay = assay,
    layers = paste(layers, collapse = ";"),
    counts_layer_only = counts_exact,
    counts_class = paste(class(counts), collapse = ";"),
    features = nrow(counts),
    cells = ncol(counts),
    nonzero_values = Matrix::nnzero(counts),
    core_fields_missing = paste(missing_core, collapse = ";"),
    cell_id_unique = !anyDuplicated(ids),
    metadata_rownames_match = identical(rownames(md), ids),
    dataset_field_valid = "dataset" %in% colnames(md) && identical(unique(as.character(md$dataset)), ds),
    identity_unique_values = ident_unique,
    identities_equal_object_name = ident_matches_object,
    identities_equal_orig_ident = ident_matches_orig,
    reductions = paste(names(obj@reductions), collapse = ";"),
    graphs = paste(names(obj@graphs), collapse = ";"),
    commands = paste(names(obj@commands), collapse = ";"),
    technical_pass = counts_exact && !length(missing_core) && !anyDuplicated(ids) && identical(rownames(md), ids)
  )

  fields_to_audit <- unique(c(core_fields, candidate_fields))
  metadata_audit[[ds]] <- data.table::rbindlist(lapply(fields_to_audit, function(f) {
    if (!f %in% colnames(md)) return(data.table(dataset = ds, field = f, present = FALSE, class = NA_character_, nonmissing_values = 0L, informative_values = 0L, unique_informative_values = 0L, examples = ""))
    z <- as.character(md[[f]])
    keep <- is_informative(z)
    examples <- head(sort(unique(z[keep])), 8L)
    data.table(dataset = ds, field = f, present = TRUE, class = paste(class(md[[f]]), collapse = ";"),
               nonmissing_values = sum(!is.na(z) & nzchar(trimws(z))), informative_values = sum(keep),
               unique_informative_values = length(unique(z[keep])), examples = paste(examples, collapse = ";"))
  }), fill = TRUE)

  if (length(candidate_fields)) {
    label_inventory[[ds]] <- data.table(dataset = ds, field = candidate_fields,
                                       unique_informative_values = vapply(candidate_fields, function(f) length(unique(as.character(md[[f]])[is_informative(md[[f]])])), integer(1)),
                                       examples = vapply(candidate_fields, function(f) paste(head(sort(unique(as.character(md[[f]])[is_informative(md[[f]])])), 12L), collapse = ";"), character(1)))
  } else {
    label_inventory[[ds]] <- data.table(dataset = ds, field = NA_character_, unique_informative_values = 0L, examples = "")
  }

  uc <- data.table(dataset = ds,
                   object_name = if ("object_name" %in% colnames(md)) as.character(md$object_name) else NA_character_,
                   sample_id = if ("sample_id" %in% colnames(md)) as.character(md$sample_id) else NA_character_,
                   patient_id = if ("patient_id" %in% colnames(md)) as.character(md$patient_id) else NA_character_,
                   tissue_status = if ("tissue_status" %in% colnames(md)) as.character(md$tissue_status) else NA_character_,
                   treatment = if ("treatment" %in% colnames(md)) as.character(md$treatment) else NA_character_,
                   timepoint = if ("timepoint" %in% colnames(md)) as.character(md$timepoint) else NA_character_)
  unit_counts[[ds]] <- uc[, .(cells = .N), by = .(dataset, object_name, sample_id, patient_id, tissue_status, treatment, timepoint)]

  present <- marker_registry$gene %in% rownames(counts)
  m <- marker_registry[, .(program, gene, interpretation)]
  m[, `:=`(dataset = ds, exact_feature_present = present, positive_cells = 0, total_raw_counts = 0, detection_fraction = 0)]
  if (any(present)) {
    panel <- counts[marker_registry$gene[present], , drop = FALSE]
    pos <- as.numeric(Matrix::rowSums(panel > 0))
    totals <- as.numeric(Matrix::rowSums(panel))
    m[present, `:=`(positive_cells = pos, total_raw_counts = totals, detection_fraction = pos / ncol(counts))]
    rm(panel)
  }
  marker_audit[[ds]] <- m[, .(dataset, program, gene, exact_feature_present, positive_cells, cells = ncol(counts), detection_fraction, total_raw_counts, interpretation)]

  fi_after <- file.info(object_file)
  input_integrity[[ds]] <- data.table(dataset = ds, object_file = normalizePath(object_file, winslash = "\\", mustWork = TRUE),
                                     size_before = fi_before$size, size_after = fi_after$size,
                                     mtime_before = format(fi_before$mtime, "%Y-%m-%dT%H:%M:%S%z"),
                                     mtime_after = format(fi_after$mtime, "%Y-%m-%dT%H:%M:%S%z"),
                                     unchanged = identical(fi_before$size, fi_after$size) && identical(as.numeric(fi_before$mtime), as.numeric(fi_after$mtime)))
  append_progress(ds, "COMPLETE", sprintf("cells=%d features=%d label_fields=%d", ncol(counts), nrow(counts), length(informative_label_fields)))
  rm(obj, counts, md, uc, ident_values)
  invisible(gc())
}

object_audit_dt <- data.table::rbindlist(object_audit, fill = TRUE)
metadata_audit_dt <- data.table::rbindlist(metadata_audit, fill = TRUE)
unit_counts_dt <- data.table::rbindlist(unit_counts, fill = TRUE)
label_inventory_dt <- data.table::rbindlist(label_inventory, fill = TRUE)
marker_audit_dt <- data.table::rbindlist(marker_audit, fill = TRUE)
input_integrity_dt <- data.table::rbindlist(input_integrity, fill = TRUE)

clinical_reconciliation <- clinical[, {
  patient_ok <- is_informative(patient_id) & is_informative(patient_id_source) & toupper(patient_id_source) != "NOT FOUND"
  timepoint_ok <- is_informative(sampling_timepoint)
  response_ok <- is_informative(response_standard)
  .(manifest_objects = .N,
    manifest_cells = sum(as.numeric(cells), na.rm = TRUE),
    verified_patient_objects = sum(patient_ok),
    verified_patients = uniqueN(patient_key[patient_ok & is_informative(patient_key)]),
    verified_timepoint_objects = sum(timepoint_ok),
    verified_response_objects = sum(response_ok),
    validation_included_objects = sum(toupper(as.character(include_validation)) == "TRUE", na.rm = TRUE),
    evidence_levels = paste(sort(unique(as.character(evidence_level))), collapse = ";"))
}, by = dataset]
clinical_reconciliation <- merge(data.table(dataset = expected_datasets), clinical_reconciliation, by = "dataset", all.x = TRUE)

metadata_gate <- metadata_audit_dt[field %in% c("sample_id", "patient_id", "treatment", "timepoint"),
                                   .(unique_informative_values = max(unique_informative_values)), by = .(dataset, field)]
metadata_gate <- data.table::dcast(metadata_gate, dataset ~ field, value.var = "unique_informative_values", fill = 0)
label_gate <- label_inventory_dt[, .(informative_label_fields = sum(!is.na(field) & unique_informative_values > 1L)), by = dataset]
marker_gate <- marker_audit_dt[, .(anchor_genes_present = sum(exact_feature_present), anchor_genes_total = .N), by = dataset]
gate <- Reduce(function(x, y) merge(x, y, by = "dataset", all = TRUE), list(
  object_audit_dt[, .(dataset, technical_pass)], metadata_gate, label_gate, marker_gate, clinical_reconciliation
))
gate[, frozen_signature_files := nrow(signature_inventory)]
gate[, `:=`(
  patient_unit_ready = verified_patients >= 2L,
  compartment_labels_ready = informative_label_fields >= 1L,
  frozen_signatures_ready = frozen_signature_files >= 1L,
  confirmatory_validation_ready = technical_pass & verified_patients >= 2L & informative_label_fields >= 1L & frozen_signature_files >= 1L
)]
gate[, blocker := vapply(seq_len(.N), function(i) {
  z <- character()
  if (!technical_pass[i]) z <- c(z, "STRUCTURAL_OBJECT_FAILURE")
  if (verified_patients[i] < 2L) z <- c(z, "PATIENT_MAPPING_UNRESOLVED")
  if (informative_label_fields[i] < 1L) z <- c(z, "COMPARTMENT_LABELS_ABSENT")
  if (frozen_signature_files[i] < 1L) z <- c(z, "FROZEN_SIGNATURES_ABSENT")
  if (!length(z)) "NONE" else paste(z, collapse = ";")
}, character(1))]

write_tsv_atomic(object_audit_dt, file.path(table_dir, "dataset_object_readiness.tsv"))
write_tsv_atomic(metadata_audit_dt, file.path(table_dir, "metadata_field_readiness.tsv"))
write_tsv_atomic(unit_counts_dt, file.path(table_dir, "independent_unit_cell_counts.tsv"))
write_tsv_atomic(label_inventory_dt, file.path(table_dir, "candidate_label_field_inventory.tsv"))
write_tsv_atomic(marker_audit_dt, file.path(table_dir, "anchor_gene_detectability.tsv"))
write_tsv_atomic(clinical_reconciliation, file.path(table_dir, "clinical_manifest_reconciliation.tsv"))
write_tsv_atomic(input_integrity_dt, file.path(table_dir, "input_integrity_after_read.tsv"))
write_tsv_atomic(gate, file.path(table_dir, "validation_gate_summary.tsv"))

technical_all <- all(object_audit_dt$technical_pass)
source_unchanged <- all(input_integrity_dt$unchanged)
ready_n <- sum(gate$confirmatory_validation_ready)
status <- if (technical_all && source_unchanged && ready_n == length(expected_datasets)) "PASS" else if (technical_all && source_unchanged) "BLOCKED_VALIDATION_INPUTS" else "FAIL_TECHNICAL_READINESS"
status_lines <- c(
  paste0("status=", status),
  paste0("run_id=", run_id),
  paste0("datasets_expected=", length(expected_datasets)),
  paste0("datasets_technical_pass=", sum(object_audit_dt$technical_pass)),
  paste0("datasets_confirmatory_ready=", ready_n),
  paste0("source_rds_modified=", if (source_unchanged) 0 else sum(!input_integrity_dt$unchanged)),
  paste0("frozen_signature_files=", nrow(signature_inventory)),
  "normalization_performed=0",
  "integration_performed=0",
  "clustering_performed=0",
  "annotation_performed=0",
  "confirmatory_hypothesis_tests_performed=0",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(status_lines, file.path(run_dir, "PHASE6A_VALIDATION_READINESS_STATUS"), useBytes = TRUE)

summary_lines <- c(
  "# Phase 6A validation-readiness summary",
  "",
  paste0("Run ID: `", run_id, "`"),
  paste0("Status: **", status, "**"),
  "",
  "This run is a technical and analytical-readiness audit. It is not biological validation.",
  "",
  paste0("- Five raw Seurat objects technical PASS: ", sum(object_audit_dt$technical_pass), "/5."),
  paste0("- Source RDS unchanged after read: ", sum(input_integrity_dt$unchanged), "/5."),
  paste0("- Datasets ready for confirmatory validation: ", ready_n, "/5."),
  paste0("- Frozen signature files found: ", nrow(signature_inventory), "."),
  "- Normalization, integration, clustering, annotation and confirmatory testing performed: 0.",
  "",
  "## Gate evidence",
  "",
  paste(capture.output(print(gate[, .(dataset, verified_patients, informative_label_fields, anchor_genes_present, anchor_genes_total, confirmatory_validation_ready, blocker)])), collapse = "\n"),
  "",
  "Confirmatory work remains blocked until traceable patient mappings, compartment labels/approved annotation and project-frozen signatures are available."
)
writeLines(summary_lines, file.path(run_dir, "PHASE6A_VALIDATION_READINESS_SUMMARY.md"), useBytes = TRUE)
append_progress("ALL", "FINAL", paste0("status=", status))
message("PHASE6A_FINAL_STATUS=", status)

