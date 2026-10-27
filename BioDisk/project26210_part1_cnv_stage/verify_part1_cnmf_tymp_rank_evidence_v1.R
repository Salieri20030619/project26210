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
EVIDENCE <- file.path(ROOT, "05_cnmf/tymp_rank_evidence_v1")
OUT <- file.path(ROOT, "05_cnmf/tymp_rank_evidence_verification_v1")
MEMORY_GATE <- file.path(PROJECT, "check_memory_gate.sh")
RANKS <- 4:8

if (dir.exists(OUT) && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty verification directory: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(name, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, status = if (isTRUE(passed)) "PASS" else "FAIL",
    evidence = paste(evidence, collapse = ";")
  )
  fwrite(rbindlist(checks, fill = TRUE), file.path(OUT, "verification_checks_v1.tsv"), sep = "\t")
  if (!isTRUE(passed)) stop("Verification failed: ", name, ": ", paste(evidence, collapse = ";"))
}

gate_output <- system2(MEMORY_GATE, "cnmf_tymp_rank_evidence_independent_verification_v1", stdout = TRUE, stderr = TRUE)
gate_status <- attr(gate_output, "status")
if (is.null(gate_status)) gate_status <- 0L
gate_text <- paste(gate_output, collapse = " ")
message(gate_text)
check("memory_gate", gate_status == 0L && grepl("PASS", gate_text, fixed = TRUE), gate_text)
check("evidence_complete", file.exists(file.path(EVIDENCE, "RUN_COMPLETE")), "terminal")
check("full_verification_complete", file.exists(file.path(ROOT, "05_cnmf/full_ranks4_8_verification_v1/RUN_COMPLETE")), "terminal")

metadata <- fromJSON(file.path(EVIDENCE, "run_metadata_v1.json"), simplifyVector = TRUE)
check("metadata_schema", identical(metadata$schema_version, "part1_cnmf_tymp_rank_evidence_v1"), metadata$schema_version)
check("metadata_counts", metadata$n_cells == 11950L && metadata$n_total_patients == 26L && metadata$n_patient_level_test_patients == 20L, c(metadata$n_cells, metadata$n_total_patients, metadata$n_patient_level_test_patients))
check("metadata_patient_unit", metadata$patient_minimum_cells == 20L && grepl("patient-specific", metadata$statistical_unit, fixed = TRUE), metadata$statistical_unit)
check("metadata_no_pseudoreplication", identical(metadata$cell_level_pseudoreplication_used, FALSE), metadata$cell_level_pseudoreplication_used)
check("metadata_no_rank_selection", identical(metadata$biological_rank_selected, FALSE), metadata$biological_rank_selected)
check("metadata_no_dependency_changes", identical(metadata$dependency_changes, FALSE), metadata$dependency_changes)

gates <- fread(file.path(EVIDENCE, "memory_gate_log_v1.tsv"))
expected_stages <- c("tymp_rank_evidence_load_v1", paste0("tymp_rank_evidence_rank", RANKS))
check("memory_gate_count", nrow(gates) == 6L, nrow(gates))
check("memory_gate_stages", identical(gates$stage, expected_stages), gates$stage)
check("memory_gates_all_pass", all(gates$exit_status == 0L) && all(grepl("PASS", gates$gate_output, fixed = TRUE)), sum(grepl("PASS", gates$gate_output, fixed = TRUE)))

meta <- fread(file.path(INPUT, "cnmf_cell_metadata_v1.tsv.gz"),
              colClasses = c(exact_cell_id = "character", curated_patient_id = "character"))
patients <- fread(file.path(AUDIT, "subclone_cnmf_patient_inventory_v1.tsv"),
                  colClasses = c(curated_patient_id = "character"))
eligible_patients <- patients[eligible_for_patient_level_test_min20 == TRUE, curated_patient_id]
check("eligible_patient_count", length(eligible_patients) == 20L && !anyDuplicated(eligible_patients), length(eligible_patients))

patient_cor <- fread(file.path(EVIDENCE, "factor_patient_stratified_correlations_v1.tsv.gz"),
                     colClasses = c(curated_patient_id = "character"))
check("patient_correlation_rows", nrow(patient_cor) == 1200L && !anyDuplicated(patient_cor[, .(rank, factor, outcome, curated_patient_id)]), nrow(patient_cor))
combo_counts <- patient_cor[, .N, by = .(rank, factor, outcome)]
check("patient_correlation_combos", nrow(combo_counts) == 60L && all(combo_counts$N == 20L), c(nrow(combo_counts), range(combo_counts$N)))
check("patient_correlation_patient_set", setequal(unique(patient_cor$curated_patient_id), eligible_patients), uniqueN(patient_cor$curated_patient_id))
check("patient_correlation_cell_floor", all(patient_cor$n_cells >= 20L), range(patient_cor$n_cells))

recomputed <- list()
usage_store <- list()
top50_store <- list()
for (rank_value in RANKS) {
  usage_dt <- fread(file.path(FULL, paste0("consensus_usage_rank", rank_value, "_v1.tsv.gz")),
                    colClasses = c(exact_cell_id = "character"))
  loading_dt <- fread(file.path(FULL, paste0("consensus_loading_rank", rank_value, "_v1.tsv.gz")),
                      colClasses = c(gene = "character"))
  check(paste0("rank", rank_value, "_usage_order"), identical(usage_dt$exact_cell_id, meta$exact_cell_id), "exact")
  factor_names <- paste0("F", seq_len(rank_value))
  usage <- as.matrix(usage_dt[, ..factor_names])
  rownames(usage) <- usage_dt$exact_cell_id
  usage_store[[as.character(rank_value)]] <- usage
  loading <- as.matrix(loading_dt[, ..factor_names])
  rownames(loading) <- loading_dt$gene
  for (factor_name in factor_names) {
    top50_store[[paste(rank_value, factor_name, sep = "_")]] <- names(sort(loading[, factor_name], decreasing = TRUE))[1:50]
    for (outcome_name in c("TYMP_log1p_CP10K_v1", "infercnv_burden_over_threshold_v1")) {
      dt <- data.table(
        curated_patient_id = meta$curated_patient_id,
        score = usage[, factor_name], outcome_value = meta[[outcome_name]]
      )[curated_patient_id %in% eligible_patients]
      x <- dt[, .(
        n_cells_check = .N,
        rho_check = if (.N >= 20L && uniqueN(score) >= 3L && uniqueN(outcome_value) >= 3L) suppressWarnings(cor(score, outcome_value, method = "spearman")) else NA_real_
      ), by = curated_patient_id]
      x[, `:=`(rank = rank_value, factor = factor_name, outcome = outcome_name)]
      recomputed[[length(recomputed) + 1L]] <- x
    }
  }
}
recomputed <- rbindlist(recomputed)
cor_compare <- merge(patient_cor, recomputed, by = c("rank", "factor", "outcome", "curated_patient_id"), all = TRUE)
check("patient_correlation_merge", nrow(cor_compare) == 1200L && !anyNA(cor_compare$n_cells) && !anyNA(cor_compare$n_cells_check), nrow(cor_compare))
rho_equal <- (is.na(cor_compare$rho) & is.na(cor_compare$rho_check)) |
  (!is.na(cor_compare$rho) & !is.na(cor_compare$rho_check) & abs(cor_compare$rho - cor_compare$rho_check) < 1e-12)
check("patient_correlations_recomputed", all(rho_equal) && all(cor_compare$n_cells == cor_compare$n_cells_check), max(abs(cor_compare$rho - cor_compare$rho_check), na.rm = TRUE))

summary_table <- fread(file.path(EVIDENCE, "factor_outcome_summary_v1.tsv"))
check("factor_summary_rows", nrow(summary_table) == 60L && !anyDuplicated(summary_table[, .(rank, factor, outcome)]), nrow(summary_table))
summary_check <- recomputed[, {
  values <- rho_check[is.finite(rho_check)]
  list(
    n_patients_check = length(values),
    median_check = if (length(values)) median(values) else NA_real_,
    fraction_positive_check = if (length(values)) mean(values > 0) else NA_real_
  )
}, by = .(rank, factor, outcome)]
summary_compare <- merge(summary_table, summary_check, by = c("rank", "factor", "outcome"), all = TRUE)
check("factor_summary_patient_counts", all(summary_compare$n_patients_with_rho == summary_compare$n_patients_check), range(summary_compare$n_patients_with_rho - summary_compare$n_patients_check))
check("factor_summary_medians", max(abs(summary_compare$median_patient_rho - summary_compare$median_check), na.rm = TRUE) < 1e-12, max(abs(summary_compare$median_patient_rho - summary_compare$median_check), na.rm = TRUE))
check("factor_summary_positive_fraction", max(abs(summary_compare$fraction_patients_positive - summary_compare$fraction_positive_check), na.rm = TRUE) < 1e-12, max(abs(summary_compare$fraction_patients_positive - summary_compare$fraction_positive_check), na.rm = TRUE))
check("factor_summary_fdr_range", all(summary_table$wilcoxon_greater_fdr[is.finite(summary_table$wilcoxon_greater_fdr)] >= 0 & summary_table$wilcoxon_greater_fdr[is.finite(summary_table$wilcoxon_greater_fdr)] <= 1), range(summary_table$wilcoxon_greater_fdr, na.rm = TRUE))

factor_evidence <- fread(file.path(EVIDENCE, "factor_evidence_v1.tsv"))
candidate <- fread(file.path(EVIDENCE, "TYMP_candidate_by_rank_v1.tsv"))
check("factor_evidence_rows", nrow(factor_evidence) == sum(RANKS) && !anyDuplicated(factor_evidence[, .(rank, factor)]), nrow(factor_evidence))
check("candidate_rows", nrow(candidate) == length(RANKS) && identical(candidate$rank, RANKS), candidate$rank)
check("candidate_factors", identical(candidate$factor, c("F4", "F4", "F3", "F1", "F3")), candidate$factor)
candidate_rule <- factor_evidence[order(rank, -median_patient_rho_TYMP, -fraction_patients_positive_TYMP, factor), .SD[1L], by = rank]
check("candidate_selection_rule", identical(candidate_rule$factor, candidate$factor), candidate_rule$factor)
robust_rule <- with(factor_evidence,
  n_patients_with_rho_TYMP >= 10L & median_patient_rho_TYMP >= 0.10 &
  fraction_patients_positive_TYMP >= 0.60 & wilcoxon_greater_fdr_TYMP <= 0.10 &
  n_dominant_cells >= 50L & n_dominant_studies >= 2L & n_dominant_patients >= 3L)
check("robust_candidate_formula", identical(factor_evidence$robust_TYMP_candidate, robust_rule), sum(robust_rule))
check("candidate_robust_ranks", identical(candidate[robust_TYMP_candidate == TRUE]$rank, c(7L, 8L)), candidate[robust_TYMP_candidate == TRUE]$rank)

top_genes <- fread(file.path(EVIDENCE, "factor_top100_hvg_genes_v1.tsv.gz"))
check("top_gene_rows", nrow(top_genes) == sum(RANKS) * 100L && !anyDuplicated(top_genes[, .(rank, factor, gene_rank)]), nrow(top_genes))
check("top_gene_rank_range", all(top_genes$gene_rank >= 1L & top_genes$gene_rank <= 100L), range(top_genes$gene_rank))

cross <- fread(file.path(EVIDENCE, "factor_cross_rank_similarity_v1.tsv.gz"))
check("cross_rank_rows", nrow(cross) == 355L && !anyDuplicated(cross[, .(rank_a, factor_a, rank_b, factor_b)]), nrow(cross))
check("cross_rank_ranges", all(is.finite(cross$cell_usage_spearman)) && all(cross$cell_usage_spearman >= -1 & cross$cell_usage_spearman <= 1) && all(cross$top50_gene_jaccard >= 0 & cross$top50_gene_jaccard <= 1), c(range(cross$cell_usage_spearman), range(cross$top50_gene_jaccard)))
cross_check <- cross[, {
  usage_a <- usage_store[[as.character(rank_a[1L])]][, factor_a[1L]]
  usage_b <- usage_store[[as.character(rank_b[1L])]][, factor_b[1L]]
  genes_a <- top50_store[[paste(rank_a[1L], factor_a[1L], sep = "_")]]
  genes_b <- top50_store[[paste(rank_b[1L], factor_b[1L], sep = "_")]]
  list(
    rho_check = suppressWarnings(cor(usage_a, usage_b, method = "spearman")),
    jaccard_check = length(intersect(genes_a, genes_b)) / length(union(genes_a, genes_b))
  )
}, by = .(rank_a, factor_a, rank_b, factor_b)]
check("cross_rank_rho_recomputed", max(abs(cross$cell_usage_spearman - cross_check$rho_check)) < 1e-12, max(abs(cross$cell_usage_spearman - cross_check$rho_check)))
check("cross_rank_jaccard_recomputed", max(abs(cross$top50_gene_jaccard - cross_check$jaccard_check)) < 1e-12, max(abs(cross$top50_gene_jaccard - cross_check$jaccard_check)))

rank_evidence <- fread(file.path(EVIDENCE, "rank_selection_evidence_v1.tsv"))
check("rank_evidence_rows", nrow(rank_evidence) == 5L && identical(rank_evidence$rank, RANKS), rank_evidence$rank)
check("rank_stability_states", identical(rank_evidence[stability_gate == TRUE]$rank, c(4L, 6L, 7L)), rank_evidence[stability_gate == TRUE]$rank)
check("rank_selection_eligibility", identical(rank_evidence[eligible_for_biological_rank_selection == TRUE]$rank, 7L), rank_evidence[eligible_for_biological_rank_selection == TRUE]$rank)
check("rank_not_yet_selected", all(!grepl("selected_rank", names(rank_evidence), fixed = TRUE)), names(rank_evidence))

manifest <- fread(file.path(EVIDENCE, "file_manifest_v1.tsv"), colClasses = "character")
check("manifest_unique", !anyDuplicated(manifest$file), uniqueN(manifest$file))
for (i in seq_len(nrow(manifest))) {
  path <- file.path(EVIDENCE, manifest$file[i])
  valid <- file.exists(path) && file.info(path)$size == as.numeric(manifest$size_bytes[i]) &&
    identical(digest(path, algo = "sha256", file = TRUE), manifest$sha256[i])
  check(paste0("hash_", manifest$file[i]), valid, manifest$sha256[i])
}

final_checks <- rbindlist(checks, fill = TRUE)
run_metadata <- list(
  schema_version = "part1_cnmf_tymp_rank_evidence_verification_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_checks = nrow(final_checks), n_pass = sum(final_checks$status == "PASS"),
  n_fail = sum(final_checks$status == "FAIL"),
  evidence_bundle = "PASS",
  sole_pre_study_sensitivity_eligible_rank = 7L,
  biological_rank_selected = FALSE,
  dependency_changes = FALSE
)
write_json(run_metadata, file.path(OUT, "run_metadata_v1.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(OUT, "RUN_COMPLETE"))
print(run_metadata)
cat("PART1_CNMF_TYMP_RANK_EVIDENCE_VERIFICATION_V1_PASS\n")
