#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: 05_analyze_hif1a_sensitivity_v1.R METADATA SCORE_DIR METAVIPER_DIR CTX_DIR OUTPUT")
}
metadata_path <- args[[1L]]
score_dir <- args[[2L]]
mv_dir <- args[[3L]]
ctx_dir <- args[[4L]]
output <- args[[5L]]
part <- paste0(output, ".part")
stopifnot(file.exists(metadata_path), dir.exists(score_dir), dir.exists(mv_dir), dir.exists(ctx_dir))
if (file.exists(output) || dir.exists(output) || file.exists(part) || dir.exists(part)) stop("refusing existing output")
dir.create(part, recursive = TRUE, showWarnings = FALSE)

meta <- fread(metadata_path, colClasses = c(exact_cell_id = "character", curated_patient_id = "character"))
required <- c("exact_cell_id", "study", "curated_patient_id", "F1", "dominant_factor_rank7_v1")
if (!all(required %chin% names(meta)) || nrow(meta) != 11950L || uniqueN(meta$exact_cell_id) != 11950L) stop("metadata gate failed")

score_path <- file.path(score_dir, "cell_hif1a_downstream_score_v1.tsv.gz")
sc <- fread(score_path, colClasses = c(exact_cell_id = "character"))
if (nrow(sc) != 11950L || !identical(sc$exact_cell_id, meta$exact_cell_id)) stop("signature score cell gate failed")

activity <- readRDS(file.path(mv_dir, "rds", "metaviper_activity_entrez_by_cell_v1.rds"))
map <- fread(file.path(mv_dir, "tables", "metaviper_regulator_id_map_v1.tsv"), colClasses = "character")
hif_ids <- map[symbol == "HIF1A", entrez_id]
if (length(hif_ids) != 1L || !(hif_ids %chin% rownames(activity))) stop("unique metaVIPER HIF1A mapping gate failed")
common_cells <- colnames(activity)
if (length(common_cells) != 11921L || uniqueN(common_cells) != 11921L) stop("metaVIPER cell gate failed")
idx <- match(common_cells, meta$exact_cell_id)
if (anyNA(idx)) stop("common cell join failed")
common_meta <- copy(meta[idx])
if (!identical(common_meta$exact_cell_id, common_cells) || uniqueN(common_meta$curated_patient_id) != 20L) stop("common metadata gate failed")

score_list <- list(
  data.table(exact_cell_id = common_cells, modality = "metaVIPER_HIF1A", evidence_family = "independent", score = as.numeric(activity[hif_ids, ])),
  data.table(exact_cell_id = common_cells, modality = "frozen_downstream_5gene", evidence_family = "independent", score = sc$hif1a_downstream_score[idx]),
  data.table(exact_cell_id = common_cells, modality = "HIF1A_RNA", evidence_family = "descriptive_rna", score = sc$HIF1A_log1p_cp10k[idx])
)

auc_files <- sort(list.files(file.path(ctx_dir, "auc"), pattern = "_auc[.]csv$", full.names = TRUE))
auc_registry <- list()
if (length(auc_files)) {
  for (path in auc_files) {
    x <- fread(path, check.names = FALSE)
    if (ncol(x) < 2L) stop("invalid targeted AUCell file: ", path)
    ids <- as.character(x[[1L]])
    score_cols <- setdiff(names(x), names(x)[1L])
    hif_col <- score_cols[grepl("^HIF1A\\(\\+\\)$", score_cols)]
    if (length(hif_col) != 1L || uniqueN(ids) != 11950L) stop("targeted AUCell schema gate failed: ", path)
    stem <- sub("_auc[.]csv$", "", basename(path))
    pos <- match(common_cells, ids)
    if (anyNA(pos)) stop("targeted AUCell cell join failed: ", path)
    score_list[[length(score_list) + 1L]] <- data.table(
      exact_cell_id = common_cells,
      modality = paste0("pySCENIC_", stem),
      evidence_family = "targeted_pyscenic",
      score = as.numeric(x[[hif_col]][pos])
    )
    auc_registry[[length(auc_registry) + 1L]] <- data.table(file = path, stem = stem, regulon = hif_col, n_cells = length(ids))
  }
}
scores <- rbindlist(score_list)
if (scores[, any(!is.finite(score))] || scores[, uniqueN(exact_cell_id), by = modality][, any(V1 != 11921L)]) stop("score numeric/key gate failed")
fwrite(rbindlist(auc_registry, fill = TRUE), file.path(part, "targeted_auc_registry_v1.tsv"), sep = "\t")

patient_effects <- list()
for (mod in unique(scores$modality)) {
  values <- scores[modality == mod]
  values[, row_id := match(exact_cell_id, common_meta$exact_cell_id)]
  values[, `:=`(
    study = common_meta$study[row_id],
    curated_patient_id = common_meta$curated_patient_id[row_id],
    F1 = as.numeric(common_meta$F1[row_id]),
    dominant = common_meta$dominant_factor_rank7_v1[row_id]
  )]
  family <- unique(values$evidence_family)
  for (patient in unique(common_meta$curated_patient_id)) {
    d <- values[curated_patient_id == patient]
    study <- unique(d$study)
    eligible_cont <- nrow(d) >= 20L && length(unique(d$F1)) > 1L
    rho <- if (eligible_cont) suppressWarnings(cor(d$score, d$F1, method = "spearman")) else NA_real_
    is_f1 <- d$dominant == "F1"
    eligible_group <- sum(is_f1) >= 10L && sum(!is_f1) >= 10L
    diff <- if (eligible_group) median(d$score[is_f1]) - median(d$score[!is_f1]) else NA_real_
    patient_effects[[length(patient_effects) + 1L]] <- data.table(
      modality = mod, evidence_family = family, study = study,
      curated_patient_id = patient, endpoint = "continuous_F1_spearman",
      n_cells = nrow(d), n_f1_cells = sum(is_f1), n_other_cells = sum(!is_f1), effect = rho
    )
    patient_effects[[length(patient_effects) + 1L]] <- data.table(
      modality = mod, evidence_family = family, study = study,
      curated_patient_id = patient, endpoint = "F1_dominant_median_difference",
      n_cells = nrow(d), n_f1_cells = sum(is_f1), n_other_cells = sum(!is_f1), effect = diff
    )
  }
}
patient_effects <- rbindlist(patient_effects)
fwrite(patient_effects, file.path(part, "patient_hif1a_sensitivity_effects_v1.tsv"), sep = "\t", na = "NA")

summarize_one <- function(v) {
  v <- v[is.finite(v)]
  p <- if (length(v) >= 3L) {
    if (all(v == 0)) 1 else suppressWarnings(wilcox.test(v, mu = 0, alternative = "two.sided", exact = NULL)$p.value)
  } else NA_real_
  list(n_eligible_patients = length(v), n_positive = sum(v > 0), n_negative = sum(v < 0),
       n_zero = sum(v == 0), median_patient_effect = if (length(v)) median(v) else NA_real_, wilcoxon_p_two_sided = p)
}
summary <- patient_effects[, summarize_one(effect), by = .(modality, evidence_family, endpoint)]
summary[, BH_FDR := p.adjust(wilcoxon_p_two_sided, method = "BH"), by = .(evidence_family, endpoint)]
summary[, reporting_tier := fifelse(n_eligible_patients >= 10L, "full_statistics_fdr",
                             fifelse(n_eligible_patients >= 3L & wilcoxon_p_two_sided < 0.05,
                                     "small_n_nominal_significant", "direction_only"))]
summary[, endpoint_pass := median_patient_effect > 0 &
          ((n_eligible_patients >= 10L & BH_FDR < 0.05) |
           (n_eligible_patients >= 3L & n_eligible_patients < 10L & wilcoxon_p_two_sided < 0.05))]
setorder(summary, evidence_family, modality, endpoint)
fwrite(summary, file.path(part, "hif1a_sensitivity_patient_summary_v1.tsv"), sep = "\t", na = "NA")

meta_scores <- scores[modality == "metaVIPER_HIF1A", .(exact_cell_id, meta_score = score)]
sig_scores <- scores[modality == "frozen_downstream_5gene", .(exact_cell_id, signature_score = score)]
concordance <- merge(meta_scores, sig_scores, by = "exact_cell_id")
concordance[, curated_patient_id := common_meta$curated_patient_id[match(exact_cell_id, common_meta$exact_cell_id)]]
concordance_patient <- concordance[, .(n_cells = .N, spearman_rho = suppressWarnings(cor(meta_score, signature_score, method = "spearman"))), by = curated_patient_id]
concordance_summary <- as.data.table(summarize_one(concordance_patient$spearman_rho))
fwrite(concordance_patient, file.path(part, "metaviper_downstream_patient_concordance_v1.tsv"), sep = "\t")
fwrite(concordance_summary, file.path(part, "metaviper_downstream_concordance_summary_v1.tsv"), sep = "\t")

target_files <- sort(list.files(file.path(ctx_dir, "regulons"), pattern = "_targets[.]tsv$", full.names = TRUE))
target_sets <- list()
for (path in target_files) {
  x <- fread(path)
  stem <- sub("_targets[.]tsv$", "", basename(path))
  target_sets[[stem]] <- unique(x[regulon == "HIF1A(+)", target])
}
overlap <- list()
for (variant in c("default", "min10_only", "nes2_5_only", "combined_light")) {
  a <- target_sets[[paste0("GENIE3_", variant)]]
  b <- target_sets[[paste0("GRNBoost2_", variant)]]
  if (is.null(a)) a <- character()
  if (is.null(b)) b <- character()
  overlap[[length(overlap) + 1L]] <- data.table(
    variant = variant, n_genie3_targets = length(a), n_grnboost2_targets = length(b),
    n_intersection = length(intersect(a, b)), n_union = length(union(a, b)),
    jaccard = if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA_real_,
    shared_targets = paste(sort(intersect(a, b)), collapse = ";")
  )
}
fwrite(rbindlist(overlap), file.path(part, "cross_method_target_overlap_v1.tsv"), sep = "\t", na = "NA")

status <- list(
  status = "PASS_ANALYSIS_COMPLETE",
  n_common_cells = length(common_cells), n_patients = uniqueN(common_meta$curated_patient_id),
  n_modalities = uniqueN(scores$modality), n_targeted_auc_profiles = length(auc_files),
  independent_modalities = c("metaVIPER_HIF1A", "frozen_downstream_5gene"),
  no_existing_pyscenic_output_modified = TRUE
)
write(toJSON(status, auto_unbox = TRUE, pretty = TRUE), file.path(part, "status.json"))
file.rename(part, output)
