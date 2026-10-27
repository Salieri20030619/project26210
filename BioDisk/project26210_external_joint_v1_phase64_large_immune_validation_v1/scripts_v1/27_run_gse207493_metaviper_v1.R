#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS="8", OPENBLAS_NUM_THREADS="8", MKL_NUM_THREADS="8",
  VECLIB_MAXIMUM_THREADS="8", NUMEXPR_NUM_THREADS="8"
)
.libPaths(c(
  "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/08_metaviper_v1/r_library",
  .libPaths()
))
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(aracne.networks)
  library(data.table)
  library(Matrix)
  library(org.Hs.eg.db)
  library(viper)
})

PROJECT <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1"
INPUT <- file.path(PROJECT, "inputs_v1/GSE207493/tf_regulon_v1")
OUT <- file.path(PROJECT, "results_v1/GSE207493/tf_regulon_v1/metaviper_v1")
TF_FILE <- "/mnt/BioDisk/ludingyu/project26210/references/pyscenic/hg38/allTFs_hg38.txt"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
CHECKPOINTS <- file.path(OUT, "patient_checkpoints_v1")
dir.create(CHECKPOINTS, recursive=TRUE, showWarnings=FALSE)
MIN_REGULON_SIZE <- 25L
SEED <- 26210L

if (file.exists(file.path(OUT, "RUN_COMPLETE"))) {
  cat("GSE207493_METAVIPER_V1_REUSED\n")
  quit(status=0)
}
required <- c("RUN_COMPLETE", "raw_counts_genes_by_cells_v1.mtx.gz", "genes_v1.tsv",
              "cells_v1.tsv", "cell_metadata_v1.tsv.gz")
if (any(!file.exists(file.path(INPUT, required)))) stop("Missing frozen input")

genes <- fread(file.path(INPUT, "genes_v1.tsv"), colClasses="character")$gene
cells <- fread(file.path(INPUT, "cells_v1.tsv"), colClasses="character")$cell_barcode
meta <- fread(file.path(INPUT, "cell_metadata_v1.tsv.gz"), colClasses=c(cell_barcode="character", sample="character"))
if (!identical(meta$cell_barcode, cells) || nrow(meta) != 39912L) stop("Metadata order/count mismatch")
con <- gzfile(file.path(INPUT, "raw_counts_genes_by_cells_v1.mtx.gz"), open="rb")
counts <- as(readMM(con), "dgCMatrix"); close(con)
dimnames(counts) <- list(genes, cells)
if (!identical(dim(counts), c(length(genes), length(cells)))) stop("Matrix mismatch")
if (any(!is.finite(counts@x)) || any(counts@x < 0)) stop("Invalid counts")

mapped <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=genes, keytype="SYMBOL", column="ENTREZID", multiVals="first")
map_dt <- data.table(gene=genes, entrez_id=as.character(mapped[genes]))
idx <- which(!is.na(map_dt$entrez_id))
entrez <- unique(map_dt$entrez_id[idx])
aggregation <- sparseMatrix(i=match(map_dt$entrez_id[idx], entrez), j=seq_along(idx), x=1,
                            dims=c(length(entrez), length(idx)))
counts_entrez <- aggregation %*% counts[idx, , drop=FALSE]
dimnames(counts_entrez) <- list(entrez, cells)
lib_size <- Matrix::colSums(counts)
if (any(lib_size <= 0)) stop("Zero library")
log_cp10k <- counts_entrez %*% Diagonal(x=1e4/lib_size)
log_cp10k@x <- log1p(log_cp10k@x)
dimnames(log_cp10k) <- list(entrez, cells)
rm(counts, counts_entrez); gc(verbose=FALSE)

data(regulonkirc); data(regulonkirp)
filter_regulon <- function(regulon, observed) {
  result <- lapply(regulon, function(entry) {
    keep <- names(entry$tfmode) %in% observed
    likelihood <- entry$likelihood
    if (length(likelihood) == length(keep)) likelihood <- likelihood[keep]
    list(tfmode=entry$tfmode[keep], likelihood=likelihood)
  })
  class(result) <- class(regulon)
  result
}
kirc <- filter_regulon(regulonkirc, rownames(log_cp10k))
kirp <- filter_regulon(regulonkirp, rownames(log_cp10k))
network_qc <- data.table(
  network=c("KIRC", "KIRP"), n_regulators=c(length(kirc), length(kirp)),
  n_regulators_minsize25=c(
    sum(vapply(kirc, function(x) length(x$tfmode) >= MIN_REGULON_SIZE, logical(1))),
    sum(vapply(kirp, function(x) length(x$tfmode) >= MIN_REGULON_SIZE, logical(1)))
  )
)
fwrite(network_qc, file.path(OUT, "metaviper_network_qc_v1.tsv"), sep="\t")

endpoint_specs <- list(
  tumor_f1=list(types="Tumor", exposure="f1_score"),
  myeloid_dc_tymp=list(types=c("Myeloid", "Dendritic"), exposure="TYMP_log1p"),
  treg_tymp=list(types="Treg", exposure="TYMP_log1p")
)

residualize <- function(y, s, g) {
  design <- cbind(1, s, g)
  coefficients <- solve(crossprod(design), crossprod(design, y))
  y - design %*% coefficients
}

effect_rows <- list(); run_rows <- list(); row_counter <- 0L
samples <- sort(unique(meta$sample))
for (i in seq_along(samples)) {
  sample_id <- samples[i]
  patient_cells <- meta[sample == sample_id, cell_barcode]
  checkpoint <- file.path(CHECKPOINTS, paste0(sample_id, ".rds"))
  started <- proc.time()[["elapsed"]]
  if (file.exists(checkpoint)) {
    activity <- readRDS(checkpoint); reused <- TRUE
  } else {
    expression_dense <- as.matrix(log_cp10k[, patient_cells, drop=FALSE])
    set.seed(SEED+i)
    activity <- viper::viper(
      expression_dense, list(KIRC=kirc, KIRP=kirp), method="rank",
      minsize=MIN_REGULON_SIZE, pleiotropy=FALSE, cores=2, verbose=FALSE
    )
    rm(expression_dense); gc(verbose=FALSE)
    if (!identical(colnames(activity), patient_cells) || any(!is.finite(activity))) {
      stop("Invalid patient activity: ", sample_id)
    }
    saveRDS(activity, checkpoint, compress="gzip"); reused <- FALSE
  }
  patient_meta <- meta[match(colnames(activity), cell_barcode)]
  for (endpoint in names(endpoint_specs)) {
    spec <- endpoint_specs[[endpoint]]
    use <- which(patient_meta$cell_type %chin% spec$types)
    if (length(use) < 40L) next
    exposure <- as.numeric(patient_meta[[spec$exposure]][use])
    s_score <- as.numeric(patient_meta$S_score[use])
    g_score <- as.numeric(patient_meta$G2M_score[use])
    if (sd(exposure) == 0) next
    exposure_resid <- as.numeric(residualize(matrix(exposure, ncol=1), s_score, g_score))
    activity_resid <- residualize(t(activity[, use, drop=FALSE]), s_score, g_score)
    continuous <- as.numeric(cor(activity_resid, exposure_resid, method="pearson"))
    names(continuous) <- colnames(activity_resid)
    lower <- as.numeric(quantile(exposure_resid, 0.25, na.rm=TRUE, type=7))
    upper <- as.numeric(quantile(exposure_resid, 0.75, na.rm=TRUE, type=7))
    low_idx <- which(exposure_resid <= lower); high_idx <- which(exposure_resid >= upper)
    contrast <- colMeans(activity_resid[high_idx, , drop=FALSE]) - colMeans(activity_resid[low_idx, , drop=FALSE])
    for (metric in c("cycle_adjusted_continuous", "cycle_adjusted_high_minus_low")) {
      effect <- if (metric == "cycle_adjusted_continuous") continuous else contrast
      row_counter <- row_counter + 1L
      effect_rows[[row_counter]] <- data.table(
        sample=sample_id, endpoint=endpoint, metric=metric,
        entrez_id=names(effect), effect=as.numeric(effect), n_cells=length(use),
        n_low=length(low_idx), n_high=length(high_idx)
      )
    }
  }
  run_rows[[i]] <- data.table(
    sample=sample_id, n_cells=length(patient_cells), n_regulators=nrow(activity),
    elapsed_seconds=proc.time()[["elapsed"]]-started, checkpoint_reused=reused
  )
  fwrite(rbindlist(run_rows), file.path(OUT, "metaviper_run_qc_v1.tsv"), sep="\t")
  cat(sprintf("[metaVIPER] %d/%d %s cells=%d regulators=%d\n",
              i, length(samples), sample_id, length(patient_cells), nrow(activity)))
  rm(activity); gc(verbose=FALSE)
}

effects <- rbindlist(effect_rows)
if (!nrow(effects) || any(!is.finite(effects$effect))) stop("No valid endpoint effects")
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=unique(effects$entrez_id), keytype="ENTREZID",
                                 column="SYMBOL", multiVals="first")
effects[, symbol := as.character(symbols[entrez_id])]
setcolorder(effects, c("sample", "endpoint", "metric", "entrez_id", "symbol", "effect",
                      "n_cells", "n_low", "n_high"))
fwrite(effects, file.path(OUT, "metaviper_effects_by_sample_v1.tsv.gz"), sep="\t")

safe_wilcox <- function(x) {
  if (length(x) < 3L || all(x == 0)) return(NA_real_)
  tryCatch(wilcox.test(x, mu=0, alternative="two.sided", exact=FALSE)$p.value,
           error=function(e) NA_real_)
}
summary <- effects[, .(
  n_eligible_samples=.N, median_effect=median(effect), mean_effect=mean(effect),
  positive_fraction=mean(effect > 0), wilcoxon_p_two_sided=safe_wilcox(effect)
), by=.(endpoint, metric, entrez_id, symbol)]
summary[, BH_FDR := NA_real_]
summary[n_eligible_samples >= 10L, BH_FDR := p.adjust(wilcoxon_p_two_sided, method="BH"),
        by=.(endpoint, metric)]
tf_symbols <- unique(trimws(readLines(TF_FILE, warn=FALSE)))
summary[, is_TF := !is.na(symbol) & symbol %chin% tf_symbols]
topic_symbols <- c(
  "ILF2","HOXD3","HNF4A","HIF1A","EPAS1","NFE2L2","SREBF1","SREBF2",
  "MYC","PPARA","PPARG","STAT1","STAT3","IRF1","IRF5","RELA","RELB",
  "NFKB1","NFKB2","FOXO1","FOXO3","TFEB","MITF","ATF4","XBP1","CREB1",
  "SP1","ETS1","BATF","NFATC1","NFATC2","TOX","FOXP3","IKZF2","E2F1","E2F8"
)
summary[, topic_related := !is.na(symbol) & symbol %chin% topic_symbols]
summary[, nominal_rank_all := frank(wilcoxon_p_two_sided, ties.method="min", na.last="keep"),
        by=.(endpoint, metric)]
summary[, nominal_rank_TF := frank(fifelse(is_TF, wilcoxon_p_two_sided, NA_real_),
                                     ties.method="min", na.last="keep"),
        by=.(endpoint, metric)]
summary[, retained := is_TF & (
  (!is.na(BH_FDR) & BH_FDR < 0.05) |
  (topic_related & !is.na(wilcoxon_p_two_sided) & wilcoxon_p_two_sided < 0.05) |
  (!is.na(nominal_rank_TF) & nominal_rank_TF <= 15)
)]
summary[, evidence_tier := fifelse(!is.na(BH_FDR) & BH_FDR < 0.05, "FDR_significant",
  fifelse(topic_related & !is.na(wilcoxon_p_two_sided) & wilcoxon_p_two_sided < 0.05,
          "topic_nominal_p_lt_0.05",
    fifelse(is_TF & nominal_rank_TF <= 15, "top15_nominal_TF", "not_retained")))]
setorder(summary, endpoint, metric, wilcoxon_p_two_sided, -abs(median_effect), na.last=TRUE)
fwrite(summary, file.path(OUT, "metaviper_regulator_summary_full_v1.tsv"), sep="\t")
fwrite(summary[retained == TRUE], file.path(OUT, "metaviper_tf_summary_retained_v1.tsv"), sep="\t")
fwrite(summary[symbol %chin% topic_symbols], file.path(OUT, "metaviper_topic_regulators_v1.tsv"), sep="\t")

old_drivers <- c("ILF2", "HOXD3", "HNF4A", "HIF1A")
replication <- summary[symbol %chin% old_drivers]
replication[, prespecified_role := fifelse(symbol == "HIF1A", "contextual", "old_part1_SCENIC_driver")]
fwrite(replication, file.path(OUT, "metaviper_old_part1_driver_replication_v1.tsv"), sep="\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info_v1.txt"))
writeLines(format(Sys.time(), tz="", usetz=TRUE), file.path(OUT, "RUN_COMPLETE"))
cat(sprintf("GSE207493_METAVIPER_V1_PASS rows=%d retained=%d\n", nrow(summary), sum(summary$retained)))
