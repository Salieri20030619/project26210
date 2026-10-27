#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: 01_analyze_m11_mr_effect_v1.R FILTERED_MR_XLSX METAVIPER_DIR OUTPUT")
mr_path <- args[[1L]]; mv_dir <- args[[2L]]; output <- args[[3L]]
part <- paste0(output, ".part")
stopifnot(file.exists(mr_path), dir.exists(mv_dir))
if (file.exists(output) || dir.exists(output) || file.exists(part) || dir.exists(part)) stop("refusing existing output")
dir.create(part, recursive = TRUE, showWarnings = FALSE)

raw <- as.data.table(read_excel(mr_path, sheet = 1, col_names = FALSE))
if (nrow(raw) < 2L || ncol(raw) < 7L) stop("unexpected MR workbook structure")
setnames(raw, names(raw)[1:7], c("protein", "mr_tier", "or_ci", "mr_p", "mr_fdr", "q_statistic", "q_p"))
mr <- raw[-1L]
mr[, `:=`(protein = trimws(as.character(protein)), mr_tier = trimws(as.character(mr_tier)), or_ci = trimws(as.character(or_ci)))]
mr[, odds_ratio := suppressWarnings(as.numeric(sub("^\\s*([0-9]+(?:\\.[0-9]+)?).*$", "\\1", or_ci, perl = TRUE)))]
for (column in c("mr_p", "mr_fdr", "q_statistic", "q_p")) mr[, (column) := suppressWarnings(as.numeric(get(column)))]
mr[, mr_beta_log_or := log(odds_ratio)]
if (nrow(mr) != 8L || uniqueN(mr$protein) != 8L || any(!is.finite(mr$mr_beta_log_or))) stop("final MR target gate failed")

activity <- readRDS(file.path(mv_dir, "rds", "metaviper_activity_entrez_by_cell_v1.rds"))
meta <- as.data.table(readRDS(file.path(mv_dir, "rds", "metaviper_cell_metadata_v1.rds")))
map <- fread(file.path(mv_dir, "tables", "metaviper_regulator_id_map_v1.tsv"), colClasses = "character")
if (!identical(meta$exact_cell_id, colnames(activity)) || ncol(activity) != 11921L || uniqueN(meta$curated_patient_id) != 20L) stop("validated metaVIPER scope gate failed")
alias <- c(KDELC2 = "POGLUT3")
mr[, lookup_symbol := fifelse(protein %chin% names(alias), unname(alias[protein]), protein)]
mr <- merge(mr, map[, .(lookup_symbol = symbol, entrez_id)], by = "lookup_symbol", all.x = TRUE, sort = FALSE)
mr[, activity_row_present := !is.na(entrez_id) & entrez_id %chin% rownames(activity)]
mr[, inference_status := fifelse(activity_row_present, "inferable_validated_metaVIPER", "unavailable_no_eligible_kidney_regulon")]
mr <- mr[match(trimws(as.character(raw[-1L][[1L]])), protein)]
fwrite(mr, file.path(part, "m11_final_mr_target_inference_status_v1.tsv"), sep = "\t", na = "NA")

inferable <- mr[activity_row_present == TRUE]
if (nrow(inferable) < 3L) stop("fewer than three inferable final MR proteins")
patient_means <- rbindlist(lapply(seq_len(nrow(inferable)), function(i) {
  values <- as.numeric(activity[inferable$entrez_id[[i]], ])
  data.table(curated_patient_id = meta$curated_patient_id, activity = values)[,
    .(mean_activity = mean(activity), median_activity = median(activity), n_cells = .N),
    by = curated_patient_id][, `:=`(protein = inferable$protein[[i]], entrez_id = inferable$entrez_id[[i]])]
}))
setcolorder(patient_means, c("protein", "entrez_id", "curated_patient_id", "n_cells", "mean_activity", "median_activity"))
fwrite(patient_means, file.path(part, "m11_patient_mean_activity_v1.tsv"), sep = "\t")

points <- patient_means[, .(patient_equal_mean_activity = mean(mean_activity), patient_sd = sd(mean_activity), n_patients = .N), by = .(protein, entrez_id)]
points <- merge(points, inferable[, .(protein, entrez_id, mr_tier, odds_ratio, mr_beta_log_or, mr_p, mr_fdr)], by = c("protein", "entrez_id"), sort = FALSE)
setorder(points, mr_beta_log_or)
fwrite(points, file.path(part, "m11_mr_effect_activity_points_v1.tsv"), sep = "\t")

safe_cor <- function(x, y, method) {
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) return(c(estimate = NA_real_, p = NA_real_))
  z <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  c(estimate = unname(z$estimate), p = z$p.value)
}
pearson <- safe_cor(points$mr_beta_log_or, points$patient_equal_mean_activity, "pearson")
spearman <- safe_cor(points$mr_beta_log_or, points$patient_equal_mean_activity, "spearman")
global <- data.table(
  n_proteins = nrow(points), pearson_r = pearson[["estimate"]], pearson_p_descriptive = pearson[["p"]],
  spearman_rho = spearman[["estimate"]], spearman_p_descriptive = spearman[["p"]],
  inference_tier = ifelse(nrow(points) < 10L, "direction_only_lt10_proteins", "full_cross_protein_inference")
)
fwrite(global, file.path(part, "m11_global_correlation_v1.tsv"), sep = "\t", na = "NA")

patient_corr <- patient_means[, {
  d <- merge(.SD, inferable[, .(protein, mr_beta_log_or)], by = "protein")
  p <- safe_cor(d$mr_beta_log_or, d$mean_activity, "pearson")
  s <- safe_cor(d$mr_beta_log_or, d$mean_activity, "spearman")
  .(n_proteins = nrow(d), pearson_r = p[["estimate"]], spearman_rho = s[["estimate"]])
}, by = curated_patient_id]
fwrite(patient_corr, file.path(part, "m11_patient_cross_protein_correlations_v1.tsv"), sep = "\t")

summarize <- function(v, metric) {
  v <- v[is.finite(v)]
  p <- if (length(v) >= 3L && !all(v == 0)) suppressWarnings(wilcox.test(v, mu = 0, exact = NULL)$p.value) else NA_real_
  data.table(metric = metric, n_patients = length(v), n_positive = sum(v > 0), n_negative = sum(v < 0),
             median_patient_correlation = median(v), wilcoxon_p_two_sided = p)
}
patient_summary <- rbindlist(list(summarize(patient_corr$pearson_r, "pearson_r"), summarize(patient_corr$spearman_rho, "spearman_rho")))
patient_summary[, BH_FDR := p.adjust(wilcoxon_p_two_sided, method = "BH")]
fwrite(patient_summary, file.path(part, "m11_patient_correlation_summary_v1.tsv"), sep = "\t")

status <- list(status = "PASS_M11_ANALYSIS_COMPLETE", n_final_mr_targets = nrow(mr), n_inferable_targets = nrow(inferable),
               inferable_targets = inferable$protein, unavailable_targets = mr[activity_row_present == FALSE, protein],
               n_cells = ncol(activity), n_patients = uniqueN(meta$curated_patient_id), figures_created = FALSE)
write(toJSON(status, auto_unbox = TRUE, pretty = TRUE), file.path(part, "status.json"))
file.rename(part, output)
