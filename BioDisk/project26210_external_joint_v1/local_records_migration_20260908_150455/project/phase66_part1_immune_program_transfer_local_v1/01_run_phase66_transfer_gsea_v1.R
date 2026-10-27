#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: 01_run_phase66_transfer_gsea_v1.R <input_dir> <output_dir>')
input_dir <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = '/', mustWork = FALSE)
part_dir <- paste0(output_dir, '.part')
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) {
  stop('Refusing existing output')
}
dir.create(part_dir, recursive = TRUE)

atomic_fwrite <- function(x, path) {
  part <- paste0(path, '.part')
  fwrite(x, part, sep = '\t', quote = TRUE, na = 'NA')
  if (!file.rename(part, path)) stop('Atomic rename failed: ', path)
}

phase62_status <- fread(file.path(input_dir, 'PHASE62_STATUS_v1.tsv'))
if (nrow(phase62_status) != 1L || phase62_status$status != 'PASS_PHASE62_IMMUNE_PSEUDOBULK') {
  stop('Phase62 status is not PASS')
}

model_summary <- fread(file.path(input_dir, 'pseudobulk_model_summary_v1.tsv'))
gene_results <- fread(file.path(input_dir, 'pseudobulk_gene_differential_expression_v1.tsv.gz'))
reg_targets <- fread(file.path(input_dir, 'phase64_frozen_regulon_targets_v1.tsv'))
reg_summary <- fread(file.path(input_dir, 'phase64_frozen_regulon_summary_v1.tsv'))
tf_cross <- fread(file.path(input_dir, 'de_novo_tf_cross_method_FDR_v1.tsv'))
met_results <- fread(file.path(input_dir, 'gse207493_metabolic_contrast_summary_retained_v1.tsv'))
met_panel <- fread(file.path(input_dir, 'gse207493_metabolic_pathway_panel_v1.tsv'))

gmt_lines <- readLines(unz(file.path(input_dir, 'ReactomePathways.gmt.zip'), 'ReactomePathways.gmt'))
gmt_parts <- strsplit(gmt_lines, '\t', fixed = TRUE)
reactome <- setNames(lapply(gmt_parts, function(x) unique(toupper(x[-c(1L, 2L)]))), vapply(gmt_parts, `[[`, character(1), 2L))

tf_direction <- tf_cross[canonical_metric == 'continuous', .(
  expected_phase64_direction = ifelse(effect_SCENIC > 0 & effect_VIPER > 0, 'positive', 'negative')
), by = .(endpoint, tf = TF)]
if (tf_direction[, anyDuplicated(paste(endpoint, tf))]) stop('Duplicate continuous TF direction rows')

setorder(reg_targets, endpoint, tf_order, tf, -weight, target_gene)
reg_signatures <- reg_targets[, .(genes = list(unique(target_gene)[seq_len(min(.N, 200L))])), by = .(endpoint, tf_order, tf)]
reg_signatures <- merge(reg_signatures, tf_direction, by = c('endpoint', 'tf'), all.x = TRUE, sort = FALSE)
if (anyNA(reg_signatures$expected_phase64_direction)) stop('Missing TF directions')
reg_keys <- paste(reg_signatures$endpoint, reg_signatures$tf)
reg_summary_keys <- paste(reg_summary$endpoint, reg_summary$tf)
reg_source_sizes <- reg_summary$n_targets[match(reg_keys, reg_summary_keys)]
if (anyNA(reg_source_sizes)) stop('Missing regulon source sizes')
reg_signatures[, `:=`(
  family = 'regulon',
  signature_id = paste(endpoint, 'REGULON', sprintf('%02d', tf_order), tf, sep = '__'),
  label = paste0(tf, ' positive regulon top-weighted targets'),
  source = 'Phase64 de novo SCENIC; frozen cross-method TF list',
  n_source_genes = reg_source_sizes,
  n_frozen_genes = lengths(genes)
)]

contrast_endpoint <- c(
  F1_high_vs_low_Tumor = 'tumor_f1',
  TYMP_high_vs_low_MyeloidDC = 'myeloid_dc_tymp',
  TYMP_high_vs_low_Treg = 'treg_tymp'
)
met_primary <- met_results[adjusted_q_value < 0.05 & contrast %chin% names(contrast_endpoint)]
met_primary[, endpoint := unname(contrast_endpoint[contrast])]
met_primary <- met_primary[endpoint %chin% c('myeloid_dc_tymp', 'treg_tymp')]
met_signatures <- met_primary[, {
  genes <- reactome[[stable_id]]
  genes <- setdiff(genes, 'TYMP')
  list(genes = list(genes))
}, by = .(endpoint, contrast, stable_id, pathway, direction)]
met_source_sizes <- met_panel$genes_total[match(met_signatures$stable_id, met_panel$stable_id)]
if (anyNA(met_source_sizes)) stop('Missing Reactome source sizes')
met_signatures[, `:=`(
  family = 'metabolic',
  signature_id = paste(endpoint, 'REACTOME', stable_id, sep = '__'),
  label = pathway,
  expected_phase64_direction = direction,
  source = 'Phase64 cycle-adjusted TYMP-excluded Reactome contrast FDR<0.05',
  n_source_genes = met_source_sizes,
  n_frozen_genes = lengths(genes),
  tf_order = NA_integer_,
  tf = NA_character_
)]

signature_list <- rbindlist(list(
  reg_signatures[, .(endpoint, family, signature_id, label, expected_phase64_direction, source, n_source_genes, n_frozen_genes, genes)],
  met_signatures[, .(endpoint, family, signature_id, label, expected_phase64_direction, source, n_source_genes, n_frozen_genes, genes)]
), use.names = TRUE)
if (anyDuplicated(signature_list$signature_id)) stop('Duplicate signature IDs')

registry <- signature_list[, .(endpoint, family, signature_id, label, expected_phase64_direction, source, n_source_genes, n_frozen_genes)]
setorder(registry, endpoint, family, signature_id)
atomic_fwrite(registry, file.path(part_dir, 'TRANSFER_SIGNATURE_REGISTRY_v1.tsv'))

models <- model_summary[
  status == 'TESTED' & (
    fine_annotation %chin% c('Macrophage_C1QC', 'Monocyte_FCN1', 'Treg')
  )
]
models <- models[
  (fine_annotation %chin% c('Macrophage_C1QC', 'Monocyte_FCN1')) |
    (dataset == 'GSE294109' & contrast == 'ICB_benefit_vs_no_effect' & fine_annotation == 'Treg')
]
models[, endpoint := fifelse(fine_annotation == 'Treg', 'treg_tymp', 'myeloid_dc_tymp')]
models[, model_id := paste(dataset, contrast, fine_annotation, sep = '__')]
if (nrow(models) != 5L || anyDuplicated(models$model_id)) stop('Expected exactly five transfer models')

all_rows <- list()
coverage_rows <- list()
row_index <- 0L
for (i in seq_len(nrow(models))) {
  model <- models[i]
  de <- gene_results[dataset == model$dataset & contrast == model$contrast & fine_annotation == model$fine_annotation]
  de <- de[is.finite(signed_rank) & !is.na(gene)]
  de[, `:=`(gene_upper = toupper(gene), abs_signed_rank = abs(signed_rank))]
  setorder(de, gene_upper, -abs_signed_rank, PValue)
  de <- de[!duplicated(gene_upper)]
  ranks <- de$signed_rank
  names(ranks) <- de$gene_upper
  ranks <- sort(ranks, decreasing = TRUE)

  selected <- signature_list[endpoint == model$endpoint]
  pathways <- list()
  for (j in seq_len(nrow(selected))) {
    genes <- intersect(selected$genes[[j]], names(ranks))
    eligible <- length(genes) >= 10L && length(genes) <= 500L
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      model_id = model$model_id,
      dataset = model$dataset,
      contrast = model$contrast,
      fine_annotation = model$fine_annotation,
      endpoint = model$endpoint,
      family = selected$family[j],
      signature_id = selected$signature_id[j],
      n_frozen_genes = selected$n_frozen_genes[j],
      n_genes_in_model = length(genes),
      eligible = eligible
    )
    if (eligible) pathways[[selected$signature_id[j]]] <- genes
  }

  set.seed(26210L + i)
  result <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 10L, maxSize = 500L,
    eps = 1e-50, nproc = 1L
  )))
  if (nrow(result)) result[, leadingEdge := vapply(leadingEdge, paste, collapse = ';', character(1))]
  tested <- data.table(signature_id = names(pathways))
  tested <- merge(tested, result, by.x = 'signature_id', by.y = 'pathway', all.x = TRUE, sort = FALSE)
  tested <- merge(tested, registry, by = 'signature_id', all.x = TRUE, sort = FALSE)
  tested[, `:=`(
    model_id = model$model_id,
    dataset = model$dataset,
    contrast = model$contrast,
    fine_annotation = model$fine_annotation,
    effect_definition = model$effect_definition,
    n_units = model$n_units,
    n_benefit = model$n_benefit,
    n_no_effect = model$n_no_effect,
    n_pairs = model$n_pairs,
    status = 'TESTED'
  )]
  row_index <- row_index + 1L
  all_rows[[row_index]] <- tested
}

coverage <- rbindlist(coverage_rows)
results <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
results[, transfer_q_value := p.adjust(pval, method = 'BH'), by = .(model_id, family)]
results[, transfer_grade := fifelse(
  transfer_q_value < 0.05, 'T1_FDR_CONTEXT_TRANSFER',
  fifelse(pval < 0.05, 'T2_NOMINAL_CONTEXT_TRANSFER', 'T0_NOT_TRANSFERRED')
)]
setorder(results, dataset, contrast, fine_annotation, family, transfer_q_value, pval, signature_id)
setorder(coverage, dataset, contrast, fine_annotation, family, signature_id)

summary <- results[, .(
  signatures_tested = .N,
  FDR_lt_0_05 = sum(transfer_q_value < 0.05, na.rm = TRUE),
  nominal_p_lt_0_05 = sum(pval < 0.05, na.rm = TRUE),
  positive_NES = sum(NES > 0, na.rm = TRUE),
  negative_NES = sum(NES < 0, na.rm = TRUE)
), by = .(model_id, dataset, contrast, fine_annotation, endpoint, family, effect_definition, n_units, n_benefit, n_no_effect, n_pairs)]
setorder(summary, dataset, contrast, fine_annotation, family)

retained <- results[transfer_grade != 'T0_NOT_TRANSFERRED']
atomic_fwrite(coverage, file.path(part_dir, 'TRANSFER_SIGNATURE_COVERAGE_v1.tsv'))
atomic_fwrite(results, file.path(part_dir, 'TRANSFER_GSEA_FULL_v1.tsv'))
atomic_fwrite(retained, file.path(part_dir, 'TRANSFER_GSEA_RETAINED_v1.tsv'))
atomic_fwrite(summary, file.path(part_dir, 'TRANSFER_MODEL_SUMMARY_v1.tsv'))

fdr_total <- results[transfer_q_value < 0.05, .N]
nominal_only <- results[pval < 0.05 & transfer_q_value >= 0.05, .N]
report <- c(
  '# Phase66：免疫程序独立临床情境转移筛查', '',
  '## 设计', '',
  '- 输入为 Phase62 已完成的 12 个患者级 edgeR quasi-likelihood 模型逐基因统计；不重新读取单细胞对象。',
  '- 转移签名来自 Phase64 GSE207493：Myeloid/DC-TYMP 与 Treg-TYMP 的冻结 de novo regulon，以及循环校正、TYMP 排除后的 Reactome FDR 通路。',
  '- regulon 使用冻结 TF 的正向 regulon 中权重最高的 200 个靶基因；代谢签名保持完整 Reactome 定义并再次排除 TYMP。',
  '- 统计使用 Phase62 同一 signed sqrt(QL-F) 排名与 fgseaMultilevel；BH 在每个模型、每个签名家族内计算。', '',
  '## 范围', '',
  sprintf('- 共检验 %d 个模型：4 个髓系模型和 1 个 GSE294109 Treg 疗效模型。', nrow(models)),
  sprintf('- 签名登记 %d 个：Myeloid/DC %d 个，Treg %d 个。', nrow(registry), registry[endpoint == 'myeloid_dc_tymp', .N], registry[endpoint == 'treg_tymp', .N]),
  sprintf('- 总检验 %d 条；FDR 转移 %d 条；仅名义转移 %d 条。', nrow(results), fdr_total, nominal_only), '',
  '## 解释边界', '',
  '本阶段检验 GSE207493 程序能否在独立免疫治疗队列的临床/时间对比中形成一致的基因集富集。它不是 TYMP 连续表达相关性的直接复算，也不能单独把程序升级为疗效预测标志物。不同临床对比之间的 NES 方向不能机械视为重复。', '',
  '完整结果见 `TRANSFER_GSEA_FULL_v1.tsv`，保留结果见 `TRANSFER_GSEA_RETAINED_v1.tsv`。'
)
writeLines(report, file.path(part_dir, 'PHASE66_TRANSFER_SUMMARY_v1.md'), useBytes = TRUE)

status <- data.table(
  status = 'PASS_PHASE66_TRANSFER_GSEA',
  completed_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'),
  models = nrow(models),
  signatures = nrow(registry),
  tests = nrow(results),
  fdr_transfers = fdr_total,
  nominal_only_transfers = nominal_only,
  method = 'fgseaMultilevel on Phase62 signed sqrt QL-F ranks; BH within model and signature family'
)
atomic_fwrite(status, file.path(part_dir, 'PHASE66_STATUS_v1.tsv'))
if (!file.rename(part_dir, output_dir)) stop('Final directory rename failed')
print(status)
