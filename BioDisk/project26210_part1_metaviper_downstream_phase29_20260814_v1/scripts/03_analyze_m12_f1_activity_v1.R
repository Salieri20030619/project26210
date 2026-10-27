#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: 03_analyze_m12_f1_activity_v1.R METADATA METAVIPER_DIR M11_TARGET_STATUS OUTPUT")
metadata_path <- args[[1L]]; mv_dir <- args[[2L]]; target_path <- args[[3L]]; output <- args[[4L]]
part <- paste0(output, ".part")
stopifnot(file.exists(metadata_path), dir.exists(mv_dir), file.exists(target_path))
if (file.exists(output) || dir.exists(output) || file.exists(part) || dir.exists(part)) stop("refusing existing output")
dir.create(part, recursive = TRUE, showWarnings = FALSE)

all_meta <- fread(metadata_path, colClasses = c(exact_cell_id="character", curated_patient_id="character"))
required <- c("exact_cell_id", "study", "curated_patient_id", "F1", "dominant_factor_rank7_v1")
if (!all(required %chin% names(all_meta)) || nrow(all_meta) != 11950L || uniqueN(all_meta$exact_cell_id) != 11950L) stop("metadata gate failed")
activity <- readRDS(file.path(mv_dir, "rds", "metaviper_activity_entrez_by_cell_v1.rds"))
map <- fread(file.path(mv_dir, "tables", "metaviper_regulator_id_map_v1.tsv"), colClasses="character")
idx <- match(colnames(activity), all_meta$exact_cell_id)
if (ncol(activity) != 11921L || anyNA(idx)) stop("metaVIPER join gate failed")
meta <- copy(all_meta[idx])
if (!identical(meta$exact_cell_id, colnames(activity)) || uniqueN(meta$curated_patient_id) != 20L) stop("common scope gate failed")

targets <- fread(target_path, colClasses=c(entrez_id="character"))[activity_row_present == TRUE, .(protein, entrez_id)]
hif <- map[symbol == "HIF1A", .(protein="HIF1A", entrez_id)]
if (nrow(targets) != 4L || nrow(hif) != 1L || !(hif$entrez_id %chin% rownames(activity))) stop("target gate failed")
features <- rbindlist(list(targets[, evidence_family := "final_MR_inferable"], hif[, evidence_family := "contextual_HIF1A"]))
fwrite(features, file.path(part, "m12_feature_registry_v1.tsv"), sep="\t")

effects <- list()
for (i in seq_len(nrow(features))) {
  values <- as.numeric(activity[features$entrez_id[[i]], ])
  d <- data.table(exact_cell_id=meta$exact_cell_id, study=meta$study, curated_patient_id=meta$curated_patient_id,
                  F1=as.numeric(meta$F1), dominant=meta$dominant_factor_rank7_v1, score=values)
  for (patient in unique(meta$curated_patient_id)) {
    x <- d[curated_patient_id == patient]
    cont <- if (nrow(x) >= 20L && uniqueN(x$F1) > 1L) suppressWarnings(cor(x$score, x$F1, method="spearman")) else NA_real_
    is_f1 <- x$dominant == "F1"
    dominant_diff <- if (sum(is_f1) >= 10L && sum(!is_f1) >= 10L) median(x$score[is_f1]) - median(x$score[!is_f1]) else NA_real_
    q25 <- as.numeric(quantile(x$F1, 0.25, names=FALSE, type=7)); q75 <- as.numeric(quantile(x$F1, 0.75, names=FALSE, type=7))
    low <- x$F1 <= q25; high <- x$F1 >= q75
    quartile_diff <- if (q75 > q25 && sum(low) >= 10L && sum(high) >= 10L) median(x$score[high]) - median(x$score[low]) else NA_real_
    common <- list(protein=features$protein[[i]], entrez_id=features$entrez_id[[i]], evidence_family=features$evidence_family[[i]],
                   study=unique(x$study), curated_patient_id=patient, n_cells=nrow(x), n_f1_cells=sum(is_f1), n_other_cells=sum(!is_f1),
                   n_low_quartile=sum(low), n_high_quartile=sum(high))
    effects[[length(effects)+1L]] <- as.data.table(c(common, list(endpoint="continuous_F1_spearman", effect=cont)))
    effects[[length(effects)+1L]] <- as.data.table(c(common, list(endpoint="F1_dominant_median_difference", effect=dominant_diff)))
    effects[[length(effects)+1L]] <- as.data.table(c(common, list(endpoint="F1_Q4_minus_Q1_median_difference", effect=quartile_diff)))
  }
}
effects <- rbindlist(effects, fill=TRUE)
fwrite(effects, file.path(part, "m12_patient_activity_effects_v1.tsv"), sep="\t", na="NA")

summarize <- function(v) {
  v <- v[is.finite(v)]
  p <- if (length(v) >= 3L) { if (all(v == 0)) 1 else suppressWarnings(wilcox.test(v, mu=0, exact=NULL)$p.value) } else NA_real_
  list(n_eligible_patients=length(v), n_positive=sum(v>0), n_negative=sum(v<0), n_zero=sum(v==0),
       median_patient_effect=if(length(v)) median(v) else NA_real_, wilcoxon_p_two_sided=p)
}
summary <- effects[, summarize(effect), by=.(protein, entrez_id, evidence_family, endpoint)]
summary[, BH_FDR := p.adjust(wilcoxon_p_two_sided, method="BH"), by=.(evidence_family, endpoint)]
summary[, reporting_tier := fifelse(n_eligible_patients >= 10L, "full_statistics_fdr",
                             fifelse(n_eligible_patients >= 3L & wilcoxon_p_two_sided < 0.05, "small_n_nominal_significant", "direction_only"))]
summary[, positive_support := median_patient_effect > 0 & ((n_eligible_patients >= 10L & BH_FDR < 0.05) |
                                  (n_eligible_patients >= 3L & n_eligible_patients < 10L & wilcoxon_p_two_sided < 0.05))]
setorder(summary, evidence_family, protein, endpoint)
fwrite(summary, file.path(part, "m12_f1_activity_patient_summary_v1.tsv"), sep="\t", na="NA")
status <- list(status="PASS_M12_ANALYSIS_COMPLETE", n_cells=ncol(activity), n_patients=uniqueN(meta$curated_patient_id),
               n_primary_mr_features=nrow(targets), contextual_features="HIF1A", endpoints=3, figures_created=FALSE)
write(toJSON(status, auto_unbox=TRUE, pretty=TRUE), file.path(part, "status.json"))
file.rename(part, output)

