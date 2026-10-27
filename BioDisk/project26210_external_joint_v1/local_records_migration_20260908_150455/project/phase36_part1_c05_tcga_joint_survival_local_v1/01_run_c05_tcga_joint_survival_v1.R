suppressPackageStartupMessages({library(data.table); library(survival); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript script expression.tsv clinical.tsv output_dir")
expr_path <- args[[1]]; clinical_path <- args[[2]]; out <- args[[3]]
if (dir.exists(out)) stop("Refusing to overwrite output directory")
dir.create(out, recursive = TRUE)

genes <- c("TYMP", "MST1R", "ITGAL", "ITGB2", "ITGAM", "CD209", "PRTG", "IGDCC4")
expr <- fread(expr_path, sep = "\t", data.table = FALSE, check.names = FALSE)
expr <- expr[!is.na(expr$Hugo_Symbol) & expr$Hugo_Symbol %in% genes, , drop = FALSE]
if (!all(genes %in% expr$Hugo_Symbol)) stop("Required genes absent from TCGA expression")
sample_cols <- setdiff(names(expr), c("Hugo_Symbol", "Entrez_Gene_Id"))
primary <- sample_cols[substr(sample_cols, 14, 15) == "01"]
if (!length(primary)) stop("No TCGA primary-tumor samples")
mat <- as.matrix(expr[, primary, drop = FALSE]); storage.mode(mat) <- "double"
mat <- log2(mat + 1)
rownames(mat) <- expr$Hugo_Symbol
if (anyDuplicated(rownames(mat))) {
  mat <- rowsum(mat, group = rownames(mat), reorder = FALSE) / as.numeric(table(factor(rownames(mat), levels = unique(rownames(mat)))))
}
patients <- substr(colnames(mat), 1, 12)
if (anyDuplicated(patients)) {
  mat <- t(rowsum(t(mat), group = patients, reorder = FALSE) / as.numeric(table(factor(patients, levels = unique(patients)))))
  patients <- colnames(mat)
}

clinical <- fread(clinical_path, sep = "\t", skip = 4, data.table = FALSE, check.names = FALSE)
clinical$OS_event <- suppressWarnings(as.numeric(sub(":.*", "", clinical$OS_STATUS)))
clinical$OS_time <- suppressWarnings(as.numeric(clinical$OS_MONTHS))
clinical$PFS_event <- suppressWarnings(as.numeric(sub(":.*", "", clinical$PFS_STATUS)))
clinical$PFS_time <- suppressWarnings(as.numeric(clinical$PFS_MONTHS))
clinical$age <- suppressWarnings(as.numeric(clinical$AGE))
clinical$sex <- factor(clinical$SEX)
stage <- toupper(as.character(clinical$AJCC_PATHOLOGIC_TUMOR_STAGE))
stage <- gsub("STAGE |[ABC]$", "", stage)
stage[!stage %in% c("I", "II", "III", "IV")] <- NA
clinical$stage_simple <- factor(stage, levels = c("I", "II", "III", "IV"))

z <- function(x) as.numeric(scale(x))
sample_df <- data.frame(patient = patients, t(mat[genes, , drop = FALSE]), check.names = FALSE)
sample_df$ICAM5_ITGAL_ITGB2_primary <- (z(sample_df$ITGAL) + z(sample_df$ITGB2)) / 2
sample_df$ICAM5_ITGAM_ITGB2_sensitivity <- (z(sample_df$ITGAM) + z(sample_df$ITGB2)) / 2
axes <- c("MST1R", "ICAM5_ITGAL_ITGB2_primary", "CD209", "PRTG", "IGDCC4", "ICAM5_ITGAM_ITGB2_sensitivity")
df <- merge(sample_df, clinical[, c("PATIENT_ID", "OS_event", "OS_time", "PFS_event", "PFS_time", "age", "sex", "stage_simple")],
            by.x = "patient", by.y = "PATIENT_ID", all = FALSE)

group_rows <- list(); logrank_rows <- list(); cox_rows <- list(); interaction_rows <- list()
tidy_cox <- function(fit, receptor, endpoint, model) {
  s <- summary(fit)
  data.frame(receptor_axis = receptor, endpoint = endpoint, model = model,
             term = rownames(s$coefficients), beta = s$coefficients[, "coef"],
             hazard_ratio = s$coefficients[, "exp(coef)"], standard_error = s$coefficients[, "se(coef)"],
             z = s$coefficients[, "z"], p_two_sided = s$coefficients[, "Pr(>|z|)"],
             ci95_low = s$conf.int[, "lower .95"], ci95_high = s$conf.int[, "upper .95"], row.names = NULL)
}

for (receptor in axes) {
  for (endpoint in c("OS", "PFS")) {
    time_col <- paste0(endpoint, "_time"); event_col <- paste0(endpoint, "_event")
    d <- df[is.finite(df[[time_col]]) & df[[time_col]] > 0 & is.finite(df[[event_col]]) & is.finite(df$TYMP) & is.finite(df[[receptor]]), , drop = FALSE]
    tymp_high <- d$TYMP > median(d$TYMP)
    receptor_high <- d[[receptor]] > median(d[[receptor]])
    d$joint_group <- factor(ifelse(!tymp_high & !receptor_high, "Both_low", ifelse(tymp_high & receptor_high, "Both_high", "Single_high")),
                            levels = c("Both_low", "Single_high", "Both_high"))
    d$time <- d[[time_col]]; d$event <- d[[event_col]]
    fit <- survfit(Surv(time, event) ~ joint_group, data = d)
    fit_table <- summary(fit)$table
    for (group in levels(d$joint_group)) {
      g <- d[d$joint_group == group, ]
      median_survival <- unname(fit_table[paste0("joint_group=", group), "median"])
      group_rows[[paste(receptor, endpoint, group)]] <- data.frame(receptor_axis = receptor, endpoint = endpoint, joint_group = group,
        n_patients = nrow(g), n_events = sum(g$event), median_followup_or_survival_months = median_survival,
        median_TYMP = median(g$TYMP), median_receptor = median(g[[receptor]]))
    }
    sd <- survdiff(Surv(time, event) ~ joint_group, data = d)
    logrank_rows[[paste(receptor, endpoint)]] <- data.frame(receptor_axis = receptor, endpoint = endpoint, n_patients = nrow(d),
      n_events = sum(d$event), chisq = unname(sd$chisq), degrees_freedom = length(sd$n) - 1L,
      logrank_p_two_sided = pchisq(sd$chisq, df = length(sd$n) - 1L, lower.tail = FALSE))
    cox_rows[[paste(receptor, endpoint, "unadjusted")]] <- tidy_cox(coxph(Surv(time, event) ~ joint_group, data = d), receptor, endpoint, "UNADJUSTED_GROUP")
    adjusted <- d[complete.cases(d[, c("age", "sex", "stage_simple")]), ]
    if (nrow(adjusted) >= 50 && sum(adjusted$event) >= 20) {
      cox_rows[[paste(receptor, endpoint, "adjusted")]] <- tidy_cox(coxph(Surv(time, event) ~ joint_group + age + sex + stage_simple, data = adjusted), receptor, endpoint, "ADJUSTED_GROUP_AGE_SEX_STAGE")
      adjusted$TYMP_z <- z(adjusted$TYMP); adjusted$receptor_z <- z(adjusted[[receptor]])
      interaction_rows[[paste(receptor, endpoint)]] <- tidy_cox(coxph(Surv(time, event) ~ TYMP_z * receptor_z + age + sex + stage_simple, data = adjusted), receptor, endpoint, "ADJUSTED_CONTINUOUS_INTERACTION")
    }
  }
}
groups <- do.call(rbind, group_rows); logrank <- do.call(rbind, logrank_rows); cox <- do.call(rbind, cox_rows); interactions <- do.call(rbind, interaction_rows)
logrank$BH_FDR_within_endpoint6axes <- ave(logrank$logrank_p_two_sided, logrank$endpoint, FUN = function(p) p.adjust(p, method = "BH"))
logrank$statistical_support <- logrank$BH_FDR_within_endpoint6axes < 0.05
write.table(df, file.path(out, "c05_tcga_analysis_cohort_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(groups, file.path(out, "c05_joint_group_survival_summary_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(logrank, file.path(out, "c05_joint_group_logrank_summary_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cox, file.path(out, "c05_joint_group_cox_models_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(interactions, file.path(out, "c05_continuous_interaction_cox_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
payload <- list(status = "PASS_C05_TCGA_JOINT_SURVIVAL_ANALYSIS_COMPLETE_PRESENTATION_PENDING",
                completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), n_expression_patients = ncol(mat),
                n_merged_patients = nrow(df), receptor_axes = axes, endpoints = c("OS", "PFS"),
                median_cutoffs_only = TRUE, direction_prespecified = FALSE, figures_created = FALSE,
                interpretation = "TCGA bulk-RNA prognosis association; not cell-specific receptor activity, circulating ligand exposure, treatment resistance or causality")
part <- file.path(out, "status.json.part"); write_json(payload, part, pretty = TRUE, auto_unbox = TRUE)
if (!file.rename(part, file.path(out, "status.json"))) stop("Atomic status rename failed")
