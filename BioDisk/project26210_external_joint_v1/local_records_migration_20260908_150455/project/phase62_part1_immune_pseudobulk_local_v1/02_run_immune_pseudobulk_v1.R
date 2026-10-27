#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: 02_run_immune_pseudobulk_v1.R <project_root> <output_dir>')
project_root <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = '/', mustWork = FALSE)
if (dir.exists(output_dir) || file.exists(output_dir)) stop('Refusing existing output: ', output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c('GSE264586', 'GSE294109')
object_root <- file.path(
  project_root, 'phase6_processing_local_v1', 'runs', '20260810_234538_376', 'objects_processed'
)
compact_root <- file.path(
  project_root, 'phase14_patient_program_validation_local_v1', 'runs',
  '20260812_082000_provisional_v3'
)
minimum_cells <- 50L

atomic_fwrite <- function(x, path) {
  part <- if (grepl('\\.gz$', path)) sub('\\.gz$', '.part.gz', path) else paste0(path, '.part')
  fwrite(x, part, sep = '\t', quote = TRUE, na = 'NA')
  if (!file.rename(part, path)) stop('Atomic rename failed: ', path)
}

hallmark_table <- as.data.table(msigdbr(species = 'Homo sapiens', category = 'H'))
hallmark_pathways <- split(hallmark_table$gene_symbol, hallmark_table$gs_name)

aggregate_subtype <- function(counts, annotation, subtype) {
  selected <- annotation[phase7_fine_annotation == subtype]
  selected[, group_id := paste(patient_key, timepoint, sep = '||')]
  group_counts <- selected[, .(n_cells = .N), by = .(
    group_id, patient_key, patient_id, timepoint, treatment, response
  )]
  eligible <- group_counts[n_cells >= minimum_cells]
  selected <- selected[group_id %chin% eligible$group_id]
  if (!nrow(selected)) return(NULL)
  group_factor <- factor(selected$group_id, levels = eligible$group_id)
  if (nlevels(group_factor) == 1L) {
    aggregation <- Matrix(1, nrow = nrow(selected), ncol = 1L, sparse = TRUE)
  } else {
    aggregation <- sparse.model.matrix(~ 0 + group_factor)
  }
  pseudobulk <- counts[, selected$matrix_column, drop = FALSE] %*% aggregation
  colnames(pseudobulk) <- levels(group_factor)
  eligible <- eligible[match(colnames(pseudobulk), group_id)]
  list(counts = pseudobulk, metadata = eligible)
}

run_model <- function(pseudobulk, metadata, design, coefficient, model_info) {
  if (ncol(pseudobulk) != nrow(metadata) || nrow(design) != nrow(metadata)) stop('Model dimension mismatch')
  dge <- DGEList(counts = pseudobulk)
  keep <- filterByExpr(dge, design = design)
  if (sum(keep) < 100L) {
    return(list(
      genes = NULL, pathways = NULL,
      summary = c(model_info, list(status = 'SKIP_LT100_GENES', genes_tested = sum(keep), fdr_lt_0_05 = 0L))
    ))
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  fit_result <- tryCatch({
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmQLFit(dge, design, robust = TRUE)
    test <- glmQLFTest(fit, coef = coefficient)
    table <- as.data.table(topTags(test, n = Inf, sort.by = 'none')$table, keep.rownames = 'gene')
    table[, signed_rank := sign(logFC) * sqrt(pmax(F, 0))]
    table[, `:=`(
      dataset = model_info$dataset,
      contrast = model_info$contrast,
      fine_annotation = model_info$fine_annotation,
      effect_definition = model_info$effect_definition
    )]
    ranks <- table$signed_rank
    names(ranks) <- table$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    pathway <- as.data.table(suppressWarnings(fgseaMultilevel(
      pathways = hallmark_pathways, stats = ranks, minSize = 10L, maxSize = 500L, nproc = 1L
    )))
    if (nrow(pathway)) {
      pathway[, leadingEdge := vapply(leadingEdge, paste, collapse = ';', character(1))]
      pathway[, `:=`(
        dataset = model_info$dataset,
        contrast = model_info$contrast,
        fine_annotation = model_info$fine_annotation,
        effect_definition = model_info$effect_definition
      )]
    }
    list(table = table, pathway = pathway)
  }, error = function(e) e)
  if (inherits(fit_result, 'error')) {
    return(list(
      genes = NULL, pathways = NULL,
      summary = c(model_info, list(status = paste0('ERROR:', conditionMessage(fit_result)), genes_tested = sum(keep), fdr_lt_0_05 = NA_integer_))
    ))
  }
  list(
    genes = fit_result$table,
    pathways = fit_result$pathway,
    summary = c(model_info, list(
      status = 'TESTED', genes_tested = nrow(fit_result$table),
      fdr_lt_0_05 = sum(fit_result$table$FDR < 0.05, na.rm = TRUE),
      hallmark_fdr_lt_0_05 = sum(fit_result$pathway$padj < 0.05, na.rm = TRUE)
    ))
  )
}

gene_results <- list()
pathway_results <- list()
model_summaries <- list()
unit_audits <- list()
input_audits <- list()
result_index <- 0L

for (dataset_current in datasets) {
  object_path <- file.path(object_root, paste0(dataset_current, '_processed_clustered.rds'))
  compact_path <- file.path(
    compact_root, dataset_current, 'objects', paste0(dataset_current, '_cell_scores_compact_v1.rds')
  )
  object <- readRDS(object_path)
  compact <- as.data.table(readRDS(compact_path))
  metadata <- as.data.table(object[[]], keep.rownames = 'seurat_cell')
  metadata[, matrix_column := .I]
  annotation <- compact[, .(
    cell_id, patient_key, patient_id, sample_id, object_name, timepoint,
    treatment, response, phase7_broad_annotation, phase7_fine_annotation
  )]
  if (anyDuplicated(metadata$cell_id) || anyDuplicated(annotation$cell_id)) stop('Duplicate cell IDs: ', dataset_current)
  annotation <- merge(
    metadata[, .(cell_id, matrix_column)], annotation,
    by = 'cell_id', all = FALSE, sort = FALSE
  )
  setorder(annotation, matrix_column)
  if (nrow(annotation) != nrow(compact)) stop('Compact annotation did not map completely: ', dataset_current)
  counts <- LayerData(object[[DefaultAssay(object)]], layer = 'counts')
  if (ncol(counts) != nrow(metadata)) stop('Counts/metadata mismatch: ', dataset_current)
  immune <- annotation[phase7_broad_annotation %chin% c('T_cell', 'NK_cell', 'Myeloid')]
  subtypes <- sort(unique(immune$phase7_fine_annotation))
  input_audits[[dataset_current]] <- data.table(
    dataset = dataset_current,
    object_cells = ncol(object),
    compact_cells = nrow(compact),
    mapped_cells = nrow(annotation),
    immune_cells = nrow(immune),
    immune_subtypes = length(subtypes),
    counts_layer = TRUE,
    mapping_status = 'PASS'
  )

  for (subtype in subtypes) {
    aggregated <- aggregate_subtype(counts, immune, subtype)
    if (is.null(aggregated)) next
    pb <- aggregated$counts
    sample_meta <- aggregated$metadata

    if (dataset_current == 'GSE294109') {
      eligible <- sample_meta[
        timepoint == 'Surgery_after_treatment' & treatment != 'Cabozantinib' &
          response %chin% c('CR', 'PR', 'No_effect')
      ]
      eligible[, response_group := factor(
        fifelse(response %chin% c('CR', 'PR'), 'benefit', 'no_effect'),
        levels = c('no_effect', 'benefit')
      )]
      if (eligible[response_group == 'benefit', uniqueN(patient_key)] >= 3L &&
          eligible[response_group == 'no_effect', uniqueN(patient_key)] >= 2L) {
        eligible <- eligible[order(patient_key)]
        model_counts <- pb[, match(eligible$group_id, colnames(pb)), drop = FALSE]
        design <- model.matrix(~ response_group, data = eligible)
        info <- list(
          dataset = dataset_current,
          contrast = 'ICB_benefit_vs_no_effect',
          fine_annotation = subtype,
          effect_definition = 'logFC benefit minus no_effect',
          n_units = nrow(eligible),
          n_benefit = eligible[response_group == 'benefit', uniqueN(patient_key)],
          n_no_effect = eligible[response_group == 'no_effect', uniqueN(patient_key)],
          n_pairs = NA_integer_
        )
        unit_audits[[paste(dataset_current, subtype, 'response', sep = '||')]] <- eligible[, .(
          dataset = dataset_current, contrast = 'ICB_benefit_vs_no_effect',
          fine_annotation = subtype, patient_key, patient_id, timepoint, response,
          response_group = as.character(response_group), n_cells
        )]
        result <- run_model(model_counts, eligible, design, 2L, info)
        result_index <- result_index + 1L
        model_summaries[[result_index]] <- as.data.table(result$summary)
        gene_results[[result_index]] <- result$genes
        pathway_results[[result_index]] <- result$pathways
      }
    } else {
      for (followup in c('Cycle2', 'End_of_treatment')) {
        candidate <- sample_meta[timepoint %chin% c('Baseline', followup)]
        paired_patients <- candidate[, .(timepoints = uniqueN(timepoint)), by = patient_key][timepoints == 2L, patient_key]
        eligible <- candidate[patient_key %chin% paired_patients]
        if (uniqueN(eligible$patient_key) < 3L) next
        eligible <- eligible[order(patient_key, factor(timepoint, levels = c('Baseline', followup)))]
        eligible[, patient_factor := factor(patient_key)]
        eligible[, timepoint_factor := factor(timepoint, levels = c('Baseline', followup))]
        model_counts <- pb[, match(eligible$group_id, colnames(pb)), drop = FALSE]
        design <- model.matrix(~ patient_factor + timepoint_factor, data = eligible)
        coefficient <- ncol(design)
        contrast_name <- paste0(followup, '_vs_Baseline')
        info <- list(
          dataset = dataset_current,
          contrast = contrast_name,
          fine_annotation = subtype,
          effect_definition = paste0('logFC ', followup, ' minus Baseline'),
          n_units = nrow(eligible),
          n_benefit = NA_integer_,
          n_no_effect = NA_integer_,
          n_pairs = uniqueN(eligible$patient_key)
        )
        unit_audits[[paste(dataset_current, subtype, followup, sep = '||')]] <- eligible[, .(
          dataset = dataset_current, contrast = contrast_name,
          fine_annotation = subtype, patient_key, patient_id, timepoint, response,
          response_group = NA_character_, n_cells
        )]
        result <- run_model(model_counts, eligible, design, coefficient, info)
        result_index <- result_index + 1L
        model_summaries[[result_index]] <- as.data.table(result$summary)
        gene_results[[result_index]] <- result$genes
        pathway_results[[result_index]] <- result$pathways
      }
    }
  }
  rm(object, compact, metadata, annotation, counts, immune)
  invisible(gc())
}

genes <- rbindlist(gene_results, use.names = TRUE, fill = TRUE)
pathways <- rbindlist(pathway_results, use.names = TRUE, fill = TRUE)
models <- rbindlist(model_summaries, use.names = TRUE, fill = TRUE)
units <- rbindlist(unit_audits, use.names = TRUE, fill = TRUE)
inputs <- rbindlist(input_audits, use.names = TRUE, fill = TRUE)

setorder(genes, dataset, contrast, fine_annotation, PValue, gene)
setorder(pathways, dataset, contrast, fine_annotation, padj, pval, pathway)
setorder(models, dataset, contrast, fine_annotation)
setorder(units, dataset, contrast, fine_annotation, patient_id, timepoint)
atomic_fwrite(inputs, file.path(output_dir, 'pseudobulk_input_audit_v1.tsv'))
atomic_fwrite(units, file.path(output_dir, 'pseudobulk_unit_audit_v1.tsv'))
atomic_fwrite(models, file.path(output_dir, 'pseudobulk_model_summary_v1.tsv'))
atomic_fwrite(genes, file.path(output_dir, 'pseudobulk_gene_differential_expression_v1.tsv.gz'))
atomic_fwrite(pathways, file.path(output_dir, 'pseudobulk_hallmark_fgsea_v1.tsv'))

status <- data.table(
  status = 'PASS_PHASE62_IMMUNE_PSEUDOBULK',
  completed_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'),
  input_datasets = nrow(inputs),
  models_total = nrow(models),
  models_tested = sum(models$status == 'TESTED'),
  gene_rows = nrow(genes),
  pathway_rows = nrow(pathways),
  minimum_cells_per_patient_timepoint_subtype = minimum_cells,
  response_unit = 'patient pooled across regions',
  longitudinal_unit = 'paired patient',
  de_method = 'edgeR quasi-likelihood with TMM',
  pathway_method = 'fgseaMultilevel Hallmark; signed sqrt QL-F ranking'
)
atomic_fwrite(status, file.path(output_dir, 'PHASE62_STATUS_v1.tsv'))
print(status)
