suppressPackageStartupMessages({library(Matrix); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript script contrast.tsv.gz ligand_target.rds output_dir")
contrast_path <- args[[1]]; matrix_path <- args[[2]]; out <- args[[3]]
if (dir.exists(out)) stop("Refusing to overwrite output directory")
dir.create(out, recursive = TRUE)
contrast <- read.delim(gzfile(contrast_path), check.names = FALSE, stringsAsFactors = FALSE)
lt <- readRDS(matrix_path)
ligands <- c("MST1", "ICAM5")
if (!all(ligands %in% colnames(lt))) stop("Candidate ligands absent from NicheNet matrix")
hif_axis <- c("HIF1A", "CA9", "VEGFA", "ADM", "ANGPTL4", "LDHA")
direct_targets <- c("TYMP", hif_axis)

patient_rows <- list(); direct_rows <- list()
for (ligand in ligands) {
  weights <- as.numeric(lt[, ligand]); names(weights) <- rownames(lt)
  direct_rows[[ligand]] <- data.frame(ligand = ligand, target = direct_targets,
                                      nichenet_weight = weights[direct_targets],
                                      target_weight_percentile = vapply(weights[direct_targets], function(x) mean(weights <= x), numeric(1)))
  for (patient in unique(contrast$patient)) {
    z <- contrast[contrast$patient == patient & (contrast$detected_fraction_high + contrast$detected_fraction_low) / 2 >= 0.05, ]
    common <- intersect(z$gene, names(weights))
    z <- z[match(common, z$gene), ]
    w <- weights[common]
    keep <- is.finite(w) & is.finite(z$high_minus_low)
    pearson <- if (sum(keep) >= 100 && sd(w[keep]) > 0 && sd(z$high_minus_low[keep]) > 0) cor(w[keep], z$high_minus_low[keep], method = "pearson") else NA_real_
    spearman <- if (sum(keep) >= 100 && sd(w[keep]) > 0 && sd(z$high_minus_low[keep]) > 0) cor(w[keep], z$high_minus_low[keep], method = "spearman") else NA_real_
    module <- z[z$gene %in% hif_axis, ]
    module_w <- weights[module$gene]
    module_score <- if (nrow(module) == length(hif_axis) && sum(abs(module_w)) > 0) sum(module$high_minus_low * module_w) / sum(abs(module_w)) else NA_real_
    patient_rows[[paste(ligand, patient)]] <- data.frame(study = z$study[[1]], patient = patient, ligand = ligand,
      n_background_genes = sum(keep), complete_target_pearson = pearson,
      complete_target_spearman = spearman, hif1a_axis_weighted_response = module_score)
  }
}
patients <- do.call(rbind, patient_rows)
summarize_endpoint <- function(endpoint, family) {
  rows <- lapply(split(patients, patients$ligand), function(z) {
    v <- z[[endpoint]]; v <- v[is.finite(v)]
    p <- if (length(v) >= 3L && any(v != 0)) wilcox.test(v, mu = 0, alternative = "two.sided", exact = FALSE)$p.value else 1
    data.frame(family = family, endpoint = endpoint, ligand = z$ligand[[1]], n_eligible_patients = length(v),
      n_positive = sum(v > 0), n_negative = sum(v < 0), n_zero = sum(v == 0), median_patient_effect = median(v), wilcoxon_p_two_sided = p)
  })
  result <- do.call(rbind, rows)
  result$BH_FDR_2_ligands <- p.adjust(result$wilcoxon_p_two_sided, method = "BH")
  result$statistical_support <- ifelse(result$n_eligible_patients >= 10, result$BH_FDR_2_ligands < 0.05,
    ifelse(result$n_eligible_patients >= 3, result$wilcoxon_p_two_sided < 0.05, FALSE))
  result$direction <- ifelse(result$median_patient_effect > 0, "positive", ifelse(result$median_patient_effect < 0, "negative", "zero"))
  result
}
summary <- rbind(summarize_endpoint("complete_target_pearson", "PRIMARY_COMPLETE_TARGET"),
                 summarize_endpoint("complete_target_spearman", "SENSITIVITY_COMPLETE_TARGET"),
                 summarize_endpoint("hif1a_axis_weighted_response", "SECONDARY_HIF1A_AXIS"))
write.table(patients, file.path(out, "c03_nichenet_patient_ligand_activity_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out, "c03_nichenet_ligand_activity_summary_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, direct_rows), file.path(out, "c03_nichenet_direct_target_weights_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
payload <- list(status = "PASS_C03_NICHENET_ANALYSIS_COMPLETE_PRESENTATION_PENDING",
                completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                n_patients = length(unique(patients$patient)), n_ligands = 2L,
                two_sided_testing = TRUE, direction_prespecified = FALSE, figures_created = FALSE,
                interpretation = "NicheNet prior concordance; not measured ligand exposure, perturbation response or causality")
part <- file.path(out, "status.json.part"); write_json(payload, part, pretty = TRUE, auto_unbox = TRUE)
if (!file.rename(part, file.path(out, "status.json"))) stop("Atomic status rename failed")
