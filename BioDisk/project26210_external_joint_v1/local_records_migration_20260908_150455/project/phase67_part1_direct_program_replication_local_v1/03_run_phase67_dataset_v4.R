#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop('Usage: 03_run_phase67_dataset_v4.R <dataset> <seurat_rds> <metadata_tsv_gz> <input_dir> <output_dir>')
dataset <- args[[1L]]
seurat_path <- normalizePath(args[[2L]], winslash = '/', mustWork = TRUE)
metadata_path <- normalizePath(args[[3L]], winslash = '/', mustWork = TRUE)
input_dir <- normalizePath(args[[4L]], winslash = '/', mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], winslash = '/', mustWork = FALSE)
part_dir <- paste0(output_dir, '.part')
if (!dataset %chin% c('GSE264586', 'GSE294109')) stop('Unexpected dataset')
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop('Refusing existing output')
dir.create(part_dir, recursive = TRUE)

SEED <- 26210L
MIN_CELLS <- 50L

atomic_fwrite <- function(x, path) {
  part <- if (grepl('\\.gz$', path)) sub('\\.gz$', '.part.gz', path) else paste0(path, '.part')
  fwrite(x, part, sep = '\t', quote = TRUE, na = 'NA')
  if (!file.rename(part, path)) stop('Atomic rename failed: ', path)
}

score_matched_control <- function(expression, genes, ctrl_size, n_bins = 25L, seed = SEED) {
  genes <- unique(intersect(toupper(genes), rownames(expression)))
  if (length(genes) < 10L) return(NULL)
  averages <- Matrix::rowMeans(expression)
  finite <- is.finite(averages)
  pool_names <- rownames(expression)[finite]
  pool_average <- averages[finite]
  n_items <- as.integer(round(length(pool_average) / (n_bins - 1L)))
  if (n_items < 1L) stop('Invalid expression bins')
  gene_rank <- data.table::frank(pool_average, ties.method = 'min')
  bins <- floor(gene_rank / n_items)
  names(bins) <- pool_names
  missing_bins <- setdiff(genes, names(bins))
  if (length(missing_bins)) stop('Signature genes missing expression bins')
  set.seed(seed)
  controls <- character()
  for (bin in unique(bins[genes])) {
    candidates <- names(bins)[bins == bin & !names(bins) %chin% genes]
    if (!length(candidates)) next
    if (length(candidates) > ctrl_size) candidates <- sample(candidates, ctrl_size, replace = FALSE)
    controls <- union(controls, candidates)
  }
  if (!length(controls)) stop('No matched control genes')
  score <- Matrix::colMeans(expression[genes, , drop = FALSE]) - Matrix::colMeans(expression[controls, , drop = FALSE])
  list(score = as.numeric(score), genes = genes, controls = controls)
}

residualize <- function(values, s_score, g2m_score) {
  design <- cbind(Intercept = 1, S_score = s_score, G2M_score = g2m_score)
  fit <- lm.fit(design, values)
  as.numeric(fit$residuals)
}

signed_rank <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values) || all(abs(values) < 1e-15)) return(c(statistic = 0, p_value = 1))
  result <- suppressWarnings(wilcox.test(values, mu = 0, alternative = 'two.sided', exact = length(values) < 50L))
  c(statistic = unname(result$statistic), p_value = result$p.value)
}

parse_reactome <- function(path) {
  lines <- readLines(unz(path, 'ReactomePathways.gmt'))
  parts <- strsplit(lines, '\t', fixed = TRUE)
  setNames(lapply(parts, function(x) unique(toupper(x[-c(1L, 2L)]))), vapply(parts, `[[`, character(1), 2L))
}

build_signatures <- function() {
  targets <- fread(file.path(input_dir, 'phase64_frozen_regulon_targets_v1.tsv'))
  cross <- fread(file.path(input_dir, 'de_novo_tf_cross_method_FDR_v1.tsv'))
  direction <- cross[canonical_metric == 'continuous', .(
    expected_phase64_direction = fifelse(effect_SCENIC > 0 & effect_VIPER > 0, 'positive', 'negative')
  ), by = .(endpoint, tf = TF)]
  setorder(targets, endpoint, tf_order, tf, -weight, target_gene)
  regulon <- targets[, .(genes = list(setdiff(unique(toupper(target_gene))[seq_len(min(.N, 200L))], 'TYMP'))), by = .(endpoint, tf_order, tf)]
  regulon <- merge(regulon, direction, by = c('endpoint', 'tf'), all.x = TRUE, sort = FALSE)
  if (anyNA(regulon$expected_phase64_direction)) stop('Missing regulon direction')
  regulon[, `:=`(
    family = 'regulon',
    signature_id = paste(endpoint, 'REGULON', sprintf('%02d', tf_order), tf, sep = '__'),
    label = paste0(tf, ' positive regulon top-weighted targets'),
    source = 'Phase64 de novo SCENIC frozen cross-method TF list'
  )]

  reactome <- parse_reactome(file.path(input_dir, 'ReactomePathways.gmt.zip'))
  metabolic <- fread(file.path(input_dir, 'gse207493_metabolic_contrast_summary_retained_v1.tsv'))
  metabolic <- metabolic[adjusted_q_value < 0.05 & contrast %chin% c('TYMP_high_vs_low_MyeloidDC', 'TYMP_high_vs_low_Treg')]
  metabolic[, endpoint := fifelse(contrast == 'TYMP_high_vs_low_MyeloidDC', 'myeloid_dc_tymp', 'treg_tymp')]
  metabolic[, genes := lapply(stable_id, function(id) setdiff(reactome[[id]], 'TYMP'))]
  metabolic[, `:=`(
    family = 'metabolic',
    signature_id = paste(endpoint, 'REACTOME', stable_id, sep = '__'),
    label = pathway,
    expected_phase64_direction = fifelse(direction == 'higher_in_feature_high', 'positive', 'negative'),
    source = 'Phase64 cycle-adjusted TYMP-excluded Reactome contrast FDR<0.05'
  )]
  signatures <- rbindlist(list(
    regulon[, .(endpoint, family, signature_id, label, expected_phase64_direction, source, genes)],
    metabolic[, .(endpoint, family, signature_id, label, expected_phase64_direction, source, genes)]
  ))
  if (nrow(signatures) != 50L || anyDuplicated(signatures$signature_id)) stop('Invalid signature registry')
  signatures
}

metadata <- fread(metadata_path)
if (anyDuplicated(metadata$cell_id) || anyDuplicated(metadata$seurat_cell_name)) stop('Duplicate metadata IDs')
signatures <- build_signatures()
cycle <- fread(file.path(input_dir, 'seurat_cell_cycle_genes_v1.tsv'))
s_genes <- toupper(cycle[phase == 'S', gene])
g2m_genes <- toupper(cycle[phase == 'G2M', gene])
if (length(s_genes) != 43L || length(g2m_genes) != 54L) stop('Invalid cycle genes')

object <- readRDS(seurat_path)
if (!inherits(object, 'Seurat')) stop('Input is not Seurat')
if (!'RNA' %chin% Assays(object)) stop('RNA assay missing')
source_object_cells <- ncol(object)
expression_all <- LayerData(object[['RNA']], layer = 'data')
if (ncol(expression_all) != ncol(object)) stop('Expression/object mismatch')
missing_cells <- setdiff(metadata$seurat_cell_name, colnames(expression_all))
if (length(missing_cells)) stop('Metadata/Seurat mismatch: ', length(missing_cells))

metadata[, lineage_v1 := fifelse(
  phase7_broad_annotation == 'Myeloid', 'Myeloid_DC',
  fifelse(phase7_fine_annotation == 'Treg', 'Treg', 'Other')
)]
selected_meta <- metadata[lineage_v1 %chin% c('Myeloid_DC', 'Treg')]
selected_columns <- match(selected_meta$seurat_cell_name, colnames(expression_all))
expression <- expression_all[, selected_columns, drop = FALSE]
colnames(expression) <- selected_meta$cell_id
rm(object, expression_all)
invisible(gc())

source_gene_names <- rownames(expression)
upper_gene_names <- toupper(source_gene_names)
collision_keys <- unique(upper_gene_names[duplicated(upper_gene_names) | duplicated(upper_gene_names, fromLast = TRUE)])
keep_genes <- !upper_gene_names %chin% collision_keys
collision_audit <- data.table(
  uppercase_key = collision_keys,
  source_names = vapply(collision_keys, function(key) paste(source_gene_names[upper_gene_names == key], collapse = ';'), character(1))
)
expression <- expression[keep_genes, , drop = FALSE]
rownames(expression) <- upper_gene_names[keep_genes]
if (anyDuplicated(rownames(expression)) || !'TYMP' %chin% rownames(expression)) stop('Invalid post-collision gene names')

correlation_rows <- list()
coverage_rows <- list()
lineage_audit <- list()
result_index <- 0L

for (endpoint in c('myeloid_dc_tymp', 'treg_tymp')) {
  endpoint_value <- endpoint
  lineage <- if (endpoint == 'myeloid_dc_tymp') 'Myeloid_DC' else 'Treg'
  lineage_meta <- selected_meta[lineage_v1 == lineage]
  if (!nrow(lineage_meta)) {
    lineage_audit[[endpoint]] <- data.table(dataset = dataset, endpoint = endpoint, status = 'UNEVALUABLE_NO_LINEAGE', n_cells = 0L, n_patients = 0L, eligible_patient_timepoints = 0L)
    next
  }
  lineage_expression <- expression[, match(lineage_meta$cell_id, colnames(expression)), drop = FALSE]
  cycle_ctrl <- min(length(intersect(s_genes, rownames(lineage_expression))), length(intersect(g2m_genes, rownames(lineage_expression))))
  s_score <- score_matched_control(lineage_expression, s_genes, ctrl_size = cycle_ctrl, seed = SEED)
  g2m_score <- score_matched_control(lineage_expression, g2m_genes, ctrl_size = cycle_ctrl, seed = SEED)
  if (is.null(s_score) || is.null(g2m_score)) stop('Cycle score failure: ', endpoint)
  tymp <- as.numeric(lineage_expression['TYMP', ])

  endpoint_signatures <- signatures[signatures[['endpoint']] == endpoint_value]
  score_values <- list()
  for (i in seq_len(nrow(endpoint_signatures))) {
    signature <- endpoint_signatures[i]
    result <- score_matched_control(
      lineage_expression, signature$genes[[1L]],
      ctrl_size = max(20L, min(50L, length(signature$genes[[1L]]))),
      seed = SEED
    )
    eligible <- !is.null(result)
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      dataset = dataset, endpoint = endpoint, family = signature$family,
      signature_id = signature$signature_id,
      n_frozen_genes = length(signature$genes[[1L]]),
      n_genes_present = if (eligible) length(result$genes) else length(intersect(signature$genes[[1L]], rownames(lineage_expression))),
      n_control_genes = if (eligible) length(result$controls) else 0L,
      eligible = eligible
    )
    if (eligible) score_values[[signature$signature_id]] <- result$score
  }

  eligible_groups <- 0L
  group_keys <- unique(lineage_meta[, .(patient_key, timepoint)])
  for (g in seq_len(nrow(group_keys))) {
    key <- group_keys[g]
    indices <- which(lineage_meta$patient_key == key$patient_key & lineage_meta$timepoint == key$timepoint)
    if (length(indices) < MIN_CELLS || uniqueN(tymp[indices]) < 3L) next
    tymp_residual <- residualize(tymp[indices], s_score$score[indices], g2m_score$score[indices])
    eligible_groups <- eligible_groups + 1L
    for (signature_id in names(score_values)) {
      signature_id_value <- signature_id
      values <- score_values[[signature_id]][indices]
      if (uniqueN(values) < 3L) next
      values_residual <- residualize(values, s_score$score[indices], g2m_score$score[indices])
      info <- endpoint_signatures[endpoint_signatures[['signature_id']] == signature_id_value]
      result_index <- result_index + 1L
      correlation_rows[[result_index]] <- data.table(
        dataset = dataset, endpoint = endpoint, family = info$family,
        signature_id = signature_id, patient_key = key$patient_key,
        patient_id = lineage_meta$patient_id[indices[1L]], timepoint = key$timepoint,
        treatment = lineage_meta$treatment[indices[1L]], response = lineage_meta$response[indices[1L]],
        n_cells = length(indices),
        raw_spearman_rho = suppressWarnings(cor(tymp[indices], values, method = 'spearman')),
        cycle_adjusted_spearman_rho = suppressWarnings(cor(tymp_residual, values_residual, method = 'spearman'))
      )
    }
  }
  lineage_audit[[endpoint]] <- data.table(
    dataset = dataset, endpoint = endpoint, status = 'PASS', n_cells = nrow(lineage_meta),
    n_patients = uniqueN(lineage_meta$patient_key), eligible_patient_timepoints = eligible_groups
  )
  rm(lineage_expression, score_values)
  invisible(gc())
}

correlations <- rbindlist(correlation_rows, use.names = TRUE, fill = TRUE)
if (!nrow(correlations)) stop('No correlations produced')
patient <- correlations[, .(
  patient_median_raw_rho = median(raw_spearman_rho, na.rm = TRUE),
  patient_median_adjusted_rho = median(cycle_adjusted_spearman_rho, na.rm = TRUE),
  n_timepoints = uniqueN(timepoint), total_cells = sum(n_cells)
), by = .(dataset, endpoint, family, signature_id, patient_key)]

summary <- patient[, {
  raw <- signed_rank(patient_median_raw_rho)
  adjusted <- signed_rank(patient_median_adjusted_rho)
  list(
    n_patients = uniqueN(patient_key), n_patient_timepoints = sum(n_timepoints),
    total_cells = sum(total_cells), median_raw_rho = median(patient_median_raw_rho),
    raw_wilcoxon_statistic = raw[['statistic']], raw_p_value = raw[['p_value']],
    median_cycle_adjusted_rho = median(patient_median_adjusted_rho),
    adjusted_wilcoxon_statistic = adjusted[['statistic']], adjusted_p_value = adjusted[['p_value']]
  )
}, by = .(dataset, endpoint, family, signature_id)]
summary[, raw_q_value := p.adjust(raw_p_value, method = 'BH'), by = .(dataset, endpoint, family)]
summary[, adjusted_q_value := p.adjust(adjusted_p_value, method = 'BH'), by = .(dataset, endpoint, family)]
registry <- signatures[, .(endpoint, family, signature_id, label, expected_phase64_direction, source, n_frozen_genes = lengths(genes))]
summary <- merge(summary, registry, by = c('endpoint', 'family', 'signature_id'), all.x = TRUE, sort = FALSE)
summary[, observed_direction := fifelse(median_cycle_adjusted_rho > 0, 'positive', 'negative')]
summary[, direction_matches_phase64 := observed_direction == expected_phase64_direction]
summary[, external_replication := direction_matches_phase64 & adjusted_q_value < 0.05]

setorder(correlations, endpoint, family, signature_id, patient_key, timepoint)
setorder(summary, endpoint, family, adjusted_q_value, signature_id)
coverage <- rbindlist(coverage_rows)
audit <- rbindlist(lineage_audit)
atomic_fwrite(registry, file.path(part_dir, paste0(dataset, '_signature_registry_v1.tsv')))
atomic_fwrite(collision_audit, file.path(part_dir, paste0(dataset, '_uppercase_collision_audit_v1.tsv')))
atomic_fwrite(coverage, file.path(part_dir, paste0(dataset, '_signature_coverage_v1.tsv')))
atomic_fwrite(audit, file.path(part_dir, paste0(dataset, '_lineage_audit_v1.tsv')))
atomic_fwrite(correlations, file.path(part_dir, paste0(dataset, '_correlations_by_patient_timepoint_v1.tsv')))
atomic_fwrite(summary, file.path(part_dir, paste0(dataset, '_external_replication_summary_v1.tsv')))
status <- data.table(
  status = 'PASS_PHASE67_DATASET_V4', dataset = dataset,
  source_object_cells = source_object_cells,
  selected_cells = ncol(expression), genes_after_collision_exclusion = nrow(expression),
  uppercase_collision_keys = length(collision_keys), correlation_rows = nrow(correlations),
  summary_rows = nrow(summary), external_fdr_replications = sum(summary$external_replication)
)
atomic_fwrite(status, file.path(part_dir, paste0(dataset, '_status_v1.tsv')))
if (!file.rename(part_dir, output_dir)) stop('Final directory rename failed')
print(status)
