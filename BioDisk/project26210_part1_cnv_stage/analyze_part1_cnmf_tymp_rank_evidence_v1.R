#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

PROJECT <- "/mnt/BioDisk/ludingyu/project26210"
ROOT <- file.path(PROJECT, "results/part1_rebuild_v1")
INPUT <- file.path(ROOT, "05_cnmf/input_v1")
FULL <- file.path(ROOT, "05_cnmf/full_ranks4_8_v1")
AUDIT <- file.path(ROOT, "04_subclones_tymp/input_audit_v1")
OUT <- file.path(ROOT, "05_cnmf/tymp_rank_evidence_v1")
MEMORY_GATE <- file.path(PROJECT, "check_memory_gate.sh")
RANKS <- 4:8

if (dir.exists(OUT) && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty TYMP/rank evidence directory: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

gate_records <- list()
gate <- function(stage) {
  output <- system2(MEMORY_GATE, stage, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  text <- paste(output, collapse = " ")
  message(text)
  gate_records[[length(gate_records) + 1L]] <<- data.table(stage = stage, exit_status = status, gate_output = text)
  fwrite(rbindlist(gate_records), file.path(OUT, "memory_gate_log_v1.tsv"), sep = "\t")
  if (status != 0L || !grepl("PASS", text, fixed = TRUE)) stop("Memory gate blocked ", stage)
}

gate("tymp_rank_evidence_load_v1")
for (terminal in c(
  file.path(ROOT, "05_cnmf/input_verification_v1/RUN_COMPLETE"),
  file.path(ROOT, "05_cnmf/full_ranks4_8_verification_v1/RUN_COMPLETE")
)) if (!file.exists(terminal)) stop("Missing prerequisite: ", terminal)

meta <- fread(file.path(INPUT, "cnmf_cell_metadata_v1.tsv.gz"),
              colClasses = c(exact_cell_id = "character", curated_patient_id = "character"))
patients <- fread(file.path(AUDIT, "subclone_cnmf_patient_inventory_v1.tsv"),
                  colClasses = c(curated_patient_id = "character"))
eligible_patients <- patients[eligible_for_patient_level_test_min20 == TRUE, curated_patient_id]
if (length(eligible_patients) != 20L || nrow(meta) != 11950L || anyDuplicated(meta$exact_cell_id)) {
  stop("Unexpected formal metadata or patient eligibility")
}
if (any(!is.finite(meta$TYMP_log1p_CP10K_v1)) || any(meta$TYMP_log1p_CP10K_v1 < 0) ||
    any(!is.finite(meta$infercnv_burden_over_threshold_v1))) stop("Invalid outcomes")

gene_sets <- list(
  TYMP_axis = c("TYMP", "RRM2B", "TK1", "SLC29A1", "SLC29A2", "VEGFA", "SPP1"),
  hypoxia_glycolysis = c("CA9", "NDUFA4L2", "EGLN3", "VEGFA", "BNIP3", "SLC2A1", "LDHA", "ENO1", "PGK1", "HK2"),
  cell_cycle = c("MKI67", "TOP2A", "UBE2C", "CDC20", "CCNB1", "CCNB2", "CDK1", "TK1"),
  interferon = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "STAT1"),
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT17"),
  mesenchymal_invasion = c("VIM", "FN1", "TGFBI", "MMP2", "MMP9", "COL1A1", "COL1A2")
)

patient_correlations <- function(score, outcome, outcome_name, rank_value, factor_name) {
  dt <- data.table(
    curated_patient_id = meta$curated_patient_id,
    score = score,
    outcome = outcome
  )[curated_patient_id %in% eligible_patients]
  result <- dt[, {
    rho <- if (.N >= 20L && uniqueN(score) >= 3L && uniqueN(outcome) >= 3L) {
      suppressWarnings(cor(score, outcome, method = "spearman"))
    } else NA_real_
    list(n_cells = .N, n_unique_score = uniqueN(score), n_unique_outcome = uniqueN(outcome), rho = rho)
  }, by = curated_patient_id]
  result[, `:=`(rank = rank_value, factor = factor_name, outcome = outcome_name)]
  setcolorder(result, c("rank", "factor", "outcome", "curated_patient_id", "n_cells", "n_unique_score", "n_unique_outcome", "rho"))
  result
}

usage_store <- list()
top_gene_store <- list()
patient_cor_store <- list()
factor_meta_store <- list()

for (rank_value in RANKS) {
  gate(paste0("tymp_rank_evidence_rank", rank_value))
  usage_dt <- fread(file.path(FULL, paste0("consensus_usage_rank", rank_value, "_v1.tsv.gz")),
                    colClasses = c(exact_cell_id = "character"))
  loading_dt <- fread(file.path(FULL, paste0("consensus_loading_rank", rank_value, "_v1.tsv.gz")),
                      colClasses = c(gene = "character"))
  if (!identical(usage_dt$exact_cell_id, meta$exact_cell_id)) stop("Usage order mismatch rank ", rank_value)
  factor_names <- paste0("F", seq_len(rank_value))
  if (!identical(names(usage_dt)[-1L], factor_names) || !identical(names(loading_dt)[-1L], factor_names)) {
    stop("Factor columns mismatch rank ", rank_value)
  }
  usage <- as.matrix(usage_dt[, ..factor_names])
  loading <- as.matrix(loading_dt[, ..factor_names])
  rownames(usage) <- usage_dt$exact_cell_id
  rownames(loading) <- loading_dt$gene
  if (any(!is.finite(usage)) || any(usage < 0) || max(abs(rowSums(usage) - 1)) > 1e-12) stop("Invalid usage")
  if (any(!is.finite(loading)) || any(loading < 0) || max(abs(colSums(loading) - 1)) > 1e-12) stop("Invalid loading")
  usage_store[[as.character(rank_value)]] <- usage
  dominant <- factor_names[max.col(usage, ties.method = "first")]

  for (factor_name in factor_names) {
    patient_cor_store[[length(patient_cor_store) + 1L]] <- patient_correlations(
      usage[, factor_name], meta$TYMP_log1p_CP10K_v1, "TYMP_log1p_CP10K_v1", rank_value, factor_name
    )
    patient_cor_store[[length(patient_cor_store) + 1L]] <- patient_correlations(
      usage[, factor_name], meta$infercnv_burden_over_threshold_v1,
      "infercnv_burden_over_threshold_v1", rank_value, factor_name
    )
    values <- sort(loading[, factor_name], decreasing = TRUE)
    n_export <- min(100L, length(values))
    top_genes <- names(values)[seq_len(n_export)]
    top_gene_store[[paste(rank_value, factor_name, sep = "_")]] <- top_genes[seq_len(min(50L, n_export))]
    top_gene_table <- data.table(
      rank = rank_value, factor = factor_name, gene_rank = seq_len(n_export),
      gene = top_genes, normalized_hvg_loading = unname(values[seq_len(n_export)])
    )
    top_gene_store[[paste0("table_", rank_value, "_", factor_name)]] <- top_gene_table
    dominant_index <- dominant == factor_name
    dominant_studies <- meta$study[dominant_index]
    dominant_patients <- meta$curated_patient_id[dominant_index]
    factor_meta_store[[length(factor_meta_store) + 1L]] <- data.table(
      rank = rank_value,
      factor = factor_name,
      n_dominant_cells = sum(dominant_index),
      n_dominant_studies = uniqueN(dominant_studies),
      n_dominant_patients = uniqueN(dominant_patients),
      max_dominant_study_fraction = if (sum(dominant_index) > 0L) max(table(dominant_studies)) / sum(dominant_index) else NA_real_,
      max_dominant_patient_fraction = if (sum(dominant_index) > 0L) max(table(dominant_patients)) / sum(dominant_index) else NA_real_,
      top15_genes = paste(head(top_genes, 15L), collapse = ","),
      TYMP_axis_top50 = sum(head(top_genes, 50L) %in% gene_sets$TYMP_axis),
      hypoxia_glycolysis_top50 = sum(head(top_genes, 50L) %in% gene_sets$hypoxia_glycolysis),
      cell_cycle_top50 = sum(head(top_genes, 50L) %in% gene_sets$cell_cycle),
      interferon_top50 = sum(head(top_genes, 50L) %in% gene_sets$interferon),
      epithelial_top50 = sum(head(top_genes, 50L) %in% gene_sets$epithelial),
      mesenchymal_invasion_top50 = sum(head(top_genes, 50L) %in% gene_sets$mesenchymal_invasion)
    )
  }
}

patient_cor <- rbindlist(patient_cor_store)
fwrite(patient_cor, file.path(OUT, "factor_patient_stratified_correlations_v1.tsv.gz"), sep = "\t", compress = "gzip")

factor_summary <- patient_cor[, {
  values <- rho[is.finite(rho)]
  p_value <- if (length(values) >= 3L && any(values != 0)) {
    suppressWarnings(wilcox.test(values, mu = 0, alternative = "greater", exact = FALSE)$p.value)
  } else NA_real_
  list(
    n_eligible_patients = .N,
    n_patients_with_rho = length(values),
    median_patient_rho = if (length(values)) median(values) else NA_real_,
    mean_patient_rho = if (length(values)) mean(values) else NA_real_,
    q25_patient_rho = if (length(values)) unname(quantile(values, 0.25)) else NA_real_,
    q75_patient_rho = if (length(values)) unname(quantile(values, 0.75)) else NA_real_,
    fraction_patients_positive = if (length(values)) mean(values > 0) else NA_real_,
    wilcoxon_greater_p = p_value
  )
}, by = .(rank, factor, outcome)]
factor_summary[, wilcoxon_greater_fdr := p.adjust(wilcoxon_greater_p, method = "BH"), by = outcome]

factor_meta <- rbindlist(factor_meta_store)
tymp_summary <- factor_summary[outcome == "TYMP_log1p_CP10K_v1"]
setnames(
  tymp_summary,
  c("n_patients_with_rho", "median_patient_rho", "mean_patient_rho", "q25_patient_rho", "q75_patient_rho", "fraction_patients_positive", "wilcoxon_greater_p", "wilcoxon_greater_fdr"),
  paste0(c("n_patients_with_rho", "median_patient_rho", "mean_patient_rho", "q25_patient_rho", "q75_patient_rho", "fraction_patients_positive", "wilcoxon_greater_p", "wilcoxon_greater_fdr"), "_TYMP")
)
tymp_summary[, outcome := NULL]
cnv_summary <- factor_summary[outcome == "infercnv_burden_over_threshold_v1"]
setnames(
  cnv_summary,
  c("n_patients_with_rho", "median_patient_rho", "mean_patient_rho", "q25_patient_rho", "q75_patient_rho", "fraction_patients_positive", "wilcoxon_greater_p", "wilcoxon_greater_fdr"),
  paste0(c("n_patients_with_rho", "median_patient_rho", "mean_patient_rho", "q25_patient_rho", "q75_patient_rho", "fraction_patients_positive", "wilcoxon_greater_p", "wilcoxon_greater_fdr"), "_CNV")
)
cnv_summary[, c("outcome", "n_eligible_patients") := NULL]
factor_evidence <- merge(merge(factor_meta, tymp_summary, by = c("rank", "factor")), cnv_summary, by = c("rank", "factor"), all.x = TRUE)
factor_evidence[, robust_TYMP_candidate :=
  n_patients_with_rho_TYMP >= 10L & median_patient_rho_TYMP >= 0.10 &
  fraction_patients_positive_TYMP >= 0.60 & wilcoxon_greater_fdr_TYMP <= 0.10 &
  n_dominant_cells >= 50L & n_dominant_studies >= 2L & n_dominant_patients >= 3L
]
setorder(factor_evidence, rank, -median_patient_rho_TYMP, -fraction_patients_positive_TYMP, factor)
candidate <- factor_evidence[, .SD[1L], by = rank]
candidate[, selected_TYMP_candidate_within_rank := TRUE]
factor_evidence <- merge(factor_evidence, candidate[, .(rank, selected_factor = factor)], by = "rank", all.x = TRUE)
factor_evidence[, selected_TYMP_candidate_within_rank := factor == selected_factor]
factor_evidence[, selected_factor := NULL]

top_gene_tables <- top_gene_store[grepl("^table_", names(top_gene_store))]
top_gene_table <- rbindlist(top_gene_tables)
fwrite(top_gene_table, file.path(OUT, "factor_top100_hvg_genes_v1.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(factor_summary, file.path(OUT, "factor_outcome_summary_v1.tsv"), sep = "\t")
fwrite(factor_evidence, file.path(OUT, "factor_evidence_v1.tsv"), sep = "\t")

pair_store <- list()
for (i in seq_len(length(RANKS) - 1L)) {
  for (j in (i + 1L):length(RANKS)) {
    rank_a <- RANKS[i]
    rank_b <- RANKS[j]
    usage_a <- usage_store[[as.character(rank_a)]]
    usage_b <- usage_store[[as.character(rank_b)]]
    for (factor_a in colnames(usage_a)) {
      for (factor_b in colnames(usage_b)) {
        genes_a <- top_gene_store[[paste(rank_a, factor_a, sep = "_")]]
        genes_b <- top_gene_store[[paste(rank_b, factor_b, sep = "_")]]
        pair_store[[length(pair_store) + 1L]] <- data.table(
          rank_a = rank_a, factor_a = factor_a, rank_b = rank_b, factor_b = factor_b,
          cell_usage_spearman = suppressWarnings(cor(usage_a[, factor_a], usage_b[, factor_b], method = "spearman")),
          top50_gene_jaccard = length(intersect(genes_a, genes_b)) / length(union(genes_a, genes_b))
        )
      }
    }
  }
}
cross_rank <- rbindlist(pair_store)
candidate_keys <- candidate[, .(rank, factor)]
candidate_cross <- cross_rank[
  candidate_keys, on = .(rank_a = rank, factor_a = factor), nomatch = 0L
][candidate_keys, on = .(rank_b = rank, factor_b = factor), nomatch = 0L]
fwrite(cross_rank, file.path(OUT, "factor_cross_rank_similarity_v1.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(candidate_cross, file.path(OUT, "TYMP_candidate_cross_rank_similarity_v1.tsv"), sep = "\t")

rank_qc <- fread(file.path(FULL, "rank_qc_v1.tsv"))
rank_evidence <- merge(rank_qc, candidate, by = "rank", all.x = TRUE)
stable_ranks <- rank_qc[stability_gate == TRUE, rank]
rank_evidence[, max_candidate_similarity_to_other_stable_rank := vapply(seq_len(.N), function(index) {
  rank_value <- rank[index]
  factor_value <- factor[index]
  values <- cross_rank[
    (rank_a == rank_value & factor_a == factor_value & rank_b %in% stable_ranks) |
    (rank_b == rank_value & factor_b == factor_value & rank_a %in% stable_ranks),
    cell_usage_spearman
  ]
  if (length(values)) max(values, na.rm = TRUE) else NA_real_
}, numeric(1))]
rank_evidence[, eligible_for_biological_rank_selection :=
  stability_gate == TRUE & robust_TYMP_candidate == TRUE &
  max_candidate_similarity_to_other_stable_rank >= 0.50
]
fwrite(candidate, file.path(OUT, "TYMP_candidate_by_rank_v1.tsv"), sep = "\t")
fwrite(rank_evidence, file.path(OUT, "rank_selection_evidence_v1.tsv"), sep = "\t")

metadata <- list(
  schema_version = "part1_cnmf_tymp_rank_evidence_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_cells = nrow(meta),
  n_total_patients = uniqueN(meta$curated_patient_id),
  n_patient_level_test_patients = length(eligible_patients),
  patient_minimum_cells = 20L,
  statistical_unit = "patient-specific Spearman rho; one-sample Wilcoxon across patient rho values",
  cell_level_pseudoreplication_used = FALSE,
  ranks = RANKS,
  biological_rank_selected = FALSE,
  dependency_changes = FALSE
)
write_json(metadata, file.path(OUT, "run_metadata_v1.json"), auto_unbox = TRUE, pretty = TRUE)
manifest_files <- sort(setdiff(list.files(OUT), c("file_manifest_v1.tsv", "RUN_COMPLETE")))
manifest <- data.table(
  file = manifest_files,
  size_bytes = as.numeric(file.info(file.path(OUT, manifest_files))$size),
  sha256 = vapply(file.path(OUT, manifest_files), digest, character(1), algo = "sha256", file = TRUE)
)
fwrite(manifest, file.path(OUT, "file_manifest_v1.tsv"), sep = "\t")
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(OUT, "RUN_COMPLETE"))
print(rank_evidence[, .(
  rank, stability_gate, factor, median_patient_rho_TYMP,
  fraction_patients_positive_TYMP, wilcoxon_greater_fdr_TYMP,
  robust_TYMP_candidate, max_candidate_similarity_to_other_stable_rank,
  eligible_for_biological_rank_selection
)])
cat("PART1_CNMF_TYMP_RANK_EVIDENCE_V1_PASS\n")
