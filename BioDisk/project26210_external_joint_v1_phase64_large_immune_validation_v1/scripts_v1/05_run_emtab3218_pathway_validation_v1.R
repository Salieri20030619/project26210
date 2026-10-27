args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript 05_run_emtab3218_pathway_validation_v1.R <rma_matrix> <geneset_annotation> <sample_scores> <output_dir> <minimum_genes>")
}

rma_path <- args[[1L]]
geneset_path <- args[[2L]]
sample_scores_path <- args[[3L]]
output_dir <- args[[4L]]
minimum_genes <- as.integer(args[[5L]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(frame, path) {
  write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
z_vector <- function(value) {
  value <- as.numeric(value)
  standard_deviation <- sd(value, na.rm = TRUE)
  if (!is.finite(standard_deviation) || standard_deviation == 0) return(rep(NA_real_, length(value)))
  (value - mean(value, na.rm = TRUE)) / standard_deviation
}
auc_higher <- function(score, response) {
  keep <- is.finite(score) & response %in% c(0, 1)
  score <- score[keep]
  response <- response[keep]
  if (length(unique(response)) != 2L) return(NA_real_)
  n_positive <- sum(response == 1)
  n_negative <- sum(response == 0)
  rank_sum <- sum(rank(score, ties.method = "average")[response == 1])
  (rank_sum - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
}

empty_results <- function() data.frame(
  cohort = character(), model_id = character(), model_type = character(), timepoint = character(),
  outcome = character(), pathway_id = character(), pathway_collection = character(),
  n_genes = integer(), n = integer(), n_event = integer(), effect_type = character(),
  effect = numeric(), secondary_metric_type = character(), secondary_metric = numeric(),
  p_value = numeric(), q_value = numeric(), direction = character(), theme_related = logical(),
  nominal_rank = integer(), retention_reason = character(), evidence_class = character(),
  stringsAsFactors = FALSE
)

make_row <- function(model_id, model_type, timepoint, outcome, pathway_id, collection,
                     n_genes, n, n_event, effect_type, effect, secondary_type = "",
                     secondary = NA_real_, p_value = NA_real_, direction = "") {
  data.frame(
    cohort = "E-MTAB-3218", model_id = model_id, model_type = model_type,
    timepoint = timepoint, outcome = outcome, pathway_id = pathway_id,
    pathway_collection = collection, n_genes = as.integer(n_genes), n = as.integer(n),
    n_event = as.integer(n_event), effect_type = effect_type, effect = as.numeric(effect),
    secondary_metric_type = secondary_type, secondary_metric = as.numeric(secondary),
    p_value = as.numeric(p_value), q_value = NA_real_, direction = direction,
    theme_related = FALSE, nominal_rank = NA_integer_, retention_reason = "",
    evidence_class = "NOT_RETAINED", stringsAsFactors = FALSE
  )
}

expr <- as.matrix(read.delim(
  rma_path, header = TRUE, row.names = 1L, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = ""
))
storage.mode(expr) <- "double"
genesets <- read.delim(
  geneset_path, header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE
)
samples <- read.delim(sample_scores_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(samples$matrix_column %in% colnames(expr))) stop("Sample score mapping is incomplete")
expr <- expr[, samples$matrix_column, drop = FALSE]

geneset_columns <- grep("^Geneset_", names(genesets), value = TRUE)
if (length(geneset_columns) < 5L) stop("Too few pathway annotation columns")
memberships <- list()
for (column_name in geneset_columns) {
  collection <- sub("^Geneset_", "", column_name)
  for (row_index in seq_len(nrow(genesets))) {
    probeset <- genesets$Probeset[[row_index]]
    if (!probeset %in% rownames(expr)) next
    raw_value <- trimws(genesets[[column_name]][[row_index]])
    if (is.na(raw_value) || !nzchar(raw_value) || raw_value == "NA") next
    labels <- unique(strsplit(raw_value, "[[:space:];,]+")[[1L]])
    labels <- labels[nzchar(labels) & labels != "NA"]
    for (label in labels) {
      pathway_id <- paste(collection, label, sep = "::")
      memberships[[pathway_id]] <- unique(c(memberships[[pathway_id]], probeset))
    }
  }
}
memberships <- memberships[vapply(memberships, length, integer(1L)) >= minimum_genes]
if (length(memberships) < 10L) stop("Too few pathways passed coverage")

expression_z <- t(apply(expr, 1L, z_vector))
pathway_scores <- vapply(memberships, function(probesets) {
  apply(expression_z[probesets, , drop = FALSE], 2L, median, na.rm = TRUE)
}, numeric(ncol(expr)))
if (is.vector(pathway_scores)) pathway_scores <- matrix(pathway_scores, ncol = 1L)
rownames(pathway_scores) <- samples$sample_id
colnames(pathway_scores) <- names(memberships)

score_frame <- cbind(samples, as.data.frame(pathway_scores, check.names = FALSE))
pathway_meta <- data.frame(
  pathway_id = names(memberships),
  pathway_collection = sub("::.*$", "", names(memberships)),
  n_genes = vapply(memberships, length, integer(1L)),
  genes = vapply(memberships, function(probesets) {
    symbols <- unique(genesets$Symbol[match(probesets, genesets$Probeset)])
    paste(symbols[!is.na(symbols) & nzchar(symbols)], collapse = ";")
  }, character(1L)),
  stringsAsFactors = FALSE
)

results <- empty_results()
binary_specs <- data.frame(
  timepoint = rep(c("baseline", "week4"), each = 2L),
  response_column = rep(c("orr_primary", "response20_sensitivity"), 2L),
  stringsAsFactors = FALSE
)
for (spec_index in seq_len(nrow(binary_specs))) {
  timepoint_name <- binary_specs$timepoint[[spec_index]]
  response_column <- binary_specs$response_column[[spec_index]]
  current <- score_frame[score_frame$timepoint == timepoint_name, , drop = FALSE]
  response <- suppressWarnings(as.numeric(current[[response_column]]))
  for (pathway_index in seq_len(nrow(pathway_meta))) {
    pathway_id <- pathway_meta$pathway_id[[pathway_index]]
    score <- suppressWarnings(as.numeric(current[[pathway_id]]))
    keep <- is.finite(score) & response %in% c(0, 1)
    if (sum(keep) < 10L || length(unique(response[keep])) != 2L || min(table(response[keep])) < 3L) next
    test <- suppressWarnings(wilcox.test(score[keep & response == 1], score[keep & response == 0], exact = FALSE))
    effect <- median(score[keep & response == 1]) - median(score[keep & response == 0])
    results <- rbind(results, make_row(
      paste(timepoint_name, response_column, sep = "__"), "binary_wilcoxon",
      timepoint_name, response_column, pathway_id,
      pathway_meta$pathway_collection[[pathway_index]], pathway_meta$n_genes[[pathway_index]],
      sum(keep), sum(response[keep] == 1), "median_difference_responder_minus_nonresponder",
      effect, "auc_higher_is_responder", auc_higher(score[keep], response[keep]), test$p.value,
      ifelse(effect > 0, "higher_in_response", "lower_in_response")
    ))
  }
}

paired <- merge(
  score_frame[score_frame$timepoint == "baseline", , drop = FALSE],
  score_frame[score_frame$timepoint == "week4", , drop = FALSE],
  by = "patient_id", suffixes = c("_baseline", "_week4"), all = FALSE
)
for (pathway_index in seq_len(nrow(pathway_meta))) {
  pathway_id <- pathway_meta$pathway_id[[pathway_index]]
  delta <- suppressWarnings(as.numeric(paired[[paste0(pathway_id, "_week4")]])) -
    suppressWarnings(as.numeric(paired[[paste0(pathway_id, "_baseline")]]))
  keep <- is.finite(delta)
  if (sum(keep) >= 10L) {
    test <- suppressWarnings(wilcox.test(delta[keep], mu = 0, exact = FALSE))
    effect <- median(delta[keep])
    results <- rbind(results, make_row(
      "paired__overall_change", "paired_wilcoxon", "paired", "week4_minus_baseline",
      pathway_id, pathway_meta$pathway_collection[[pathway_index]],
      pathway_meta$n_genes[[pathway_index]], sum(keep), 0L, "median_paired_delta",
      effect, p_value = test$p.value,
      direction = ifelse(effect > 0, "increased_week4", "decreased_week4")
    ))
  }
  for (response_column in c("orr_primary", "response20_sensitivity")) {
    response <- suppressWarnings(as.numeric(paired[[paste0(response_column, "_baseline")]]))
    response_keep <- keep & response %in% c(0, 1)
    if (sum(response_keep) < 10L || length(unique(response[response_keep])) != 2L || min(table(response[response_keep])) < 3L) next
    test <- suppressWarnings(wilcox.test(
      delta[response_keep & response == 1], delta[response_keep & response == 0], exact = FALSE
    ))
    effect <- median(delta[response_keep & response == 1]) - median(delta[response_keep & response == 0])
    results <- rbind(results, make_row(
      paste("paired_delta", response_column, sep = "__"), "delta_wilcoxon", "paired",
      response_column, pathway_id, pathway_meta$pathway_collection[[pathway_index]],
      pathway_meta$n_genes[[pathway_index]], sum(response_keep), sum(response[response_keep] == 1),
      "median_delta_difference_responder_minus_nonresponder", effect,
      "auc_higher_delta_is_responder", auc_higher(delta[response_keep], response[response_keep]),
      test$p.value, ifelse(effect > 0, "greater_in_response", "lower_in_response")
    ))
  }
}

for (timepoint_name in c("baseline", "week4")) {
  current <- score_frame[score_frame$timepoint == timepoint_name, , drop = FALSE]
  for (axis_name in c("f1_score_gene_z_weighted", "tymp_expression_z")) {
    axis <- suppressWarnings(as.numeric(current[[axis_name]]))
    for (pathway_index in seq_len(nrow(pathway_meta))) {
      pathway_id <- pathway_meta$pathway_id[[pathway_index]]
      score <- suppressWarnings(as.numeric(current[[pathway_id]]))
      keep <- is.finite(score) & is.finite(axis)
      if (sum(keep) < 10L) next
      test <- suppressWarnings(cor.test(score[keep], axis[keep], method = "spearman", exact = FALSE))
      effect <- unname(test$estimate)
      results <- rbind(results, make_row(
        paste(timepoint_name, "correlation", axis_name, sep = "__"), "spearman",
        timepoint_name, axis_name, pathway_id, pathway_meta$pathway_collection[[pathway_index]],
        pathway_meta$n_genes[[pathway_index]], sum(keep), 0L, "spearman_rho", effect,
        p_value = test$p.value, direction = ifelse(effect > 0, "positive", "negative")
      ))
    }
  }
}

theme_pattern <- paste(c(
  "CD3TCR", "Tcell", "IFN", "Javelin", "Merck", "RIG1", "CM009", "Adenosine",
  "IMmotion150", "Ascierto", "EMT", "HIF1A", "INFLAM", "IMMUNE", "ANGIO",
  "MYELOID", "TCELL", "T_CELL", "HYPOX", "GLYCOL", "TNFA", "IL6", "COMPLEMENT"
), collapse = "|")
results$theme_related <- grepl(theme_pattern, results$pathway_id, ignore.case = TRUE)
for (model_name in unique(results$model_id)) {
  model_rows <- which(results$model_id == model_name & is.finite(results$p_value))
  results$q_value[model_rows] <- p.adjust(results$p_value[model_rows], method = "BH")
  nominal_rows <- model_rows[results$p_value[model_rows] < 0.05]
  nominal_rows <- nominal_rows[order(results$p_value[nominal_rows], results$pathway_id[nominal_rows])]
  if (length(nominal_rows)) results$nominal_rank[nominal_rows] <- seq_along(nominal_rows)
}
is_fdr <- is.finite(results$q_value) & results$q_value < 0.05
is_theme_nominal <- results$theme_related & is.finite(results$p_value) & results$p_value < 0.05
is_top15_nominal <- is.finite(results$nominal_rank) & results$nominal_rank <= 15L
results$retention_reason <- ifelse(
  is_fdr, "FDR_LT_0.05",
  ifelse(is_theme_nominal, "THEME_NOMINAL_P_LT_0.05",
         ifelse(is_top15_nominal, "TOP15_NOMINAL_P_LT_0.05", ""))
)
results$evidence_class <- ifelse(is_fdr, "C1_FDR", ifelse(nzchar(results$retention_reason), "C2_NOMINAL", "NOT_RETAINED"))
results <- results[order(results$model_id, results$p_value, results$pathway_id), , drop = FALSE]
retained <- results[nzchar(results$retention_reason), , drop = FALSE]

write_tsv(pathway_meta, file.path(output_dir, "emtab3218_pathway_coverage_v1.tsv"))
write_tsv(score_frame, file.path(output_dir, "emtab3218_pathway_scores_v1.tsv"))
write_tsv(results, file.path(output_dir, "emtab3218_pathway_models_all_v1.tsv"))
write_tsv(retained, file.path(output_dir, "emtab3218_pathway_models_retained_v1.tsv"))

model_gate <- do.call(rbind, lapply(unique(results$model_id), function(model_name) {
  current <- results[results$model_id == model_name, , drop = FALSE]
  data.frame(
    model_id = model_name, n_pathways_tested = nrow(current),
    n_nominal_p_lt_0.05 = sum(current$p_value < 0.05, na.rm = TRUE),
    n_fdr_q_lt_0.05 = sum(current$q_value < 0.05, na.rm = TRUE),
    n_retained = sum(nzchar(current$retention_reason)),
    top15_complete = sum(current$retention_reason == "TOP15_NOMINAL_P_LT_0.05") <= 15L,
    stringsAsFactors = FALSE
  )
}))
write_tsv(model_gate, file.path(output_dir, "emtab3218_pathway_model_gate_v1.tsv"))

status <- data.frame(
  gate = c("pathways_scored", "model_rows", "model_count", "retained_rows", "top15_rule_valid"),
  observed = c(length(memberships), nrow(results), length(unique(results$model_id)), nrow(retained), all(model_gate$top15_complete)),
  expected_or_threshold = c(10, 100, 5, 0, TRUE),
  status = c(
    ifelse(length(memberships) >= 10L, "PASS", "FAIL"),
    ifelse(nrow(results) >= 100L, "PASS", "FAIL"),
    ifelse(length(unique(results$model_id)) >= 5L, "PASS", "FAIL"),
    "INFO",
    ifelse(all(model_gate$top15_complete), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_dir, "emtab3218_pathway_validation_status_v1.tsv"))
print(status, row.names = FALSE)
print(model_gate, row.names = FALSE)
if (any(status$status == "FAIL")) quit(status = 2L)
cat("PASS_E_MTAB_3218_PATHWAY_VALIDATION\n")
