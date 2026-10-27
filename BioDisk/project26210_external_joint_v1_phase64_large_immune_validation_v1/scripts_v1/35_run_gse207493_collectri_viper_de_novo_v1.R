#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS="8", OPENBLAS_NUM_THREADS="8", MKL_NUM_THREADS="8", NUMEXPR_NUM_THREADS="8")
suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(Matrix)
})

PROJECT <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1"
INPUT <- file.path(PROJECT, "inputs_v1/GSE207493/tf_regulon_v1")
NETWORK_FILE <- file.path(PROJECT, "inputs_v1/GSE207493/tf_regulon_de_novo_v1/public_regulons_v1/collectri_human_frozen_v1.tsv.gz")
OUT <- file.path(PROJECT, "results_v1/GSE207493/tf_regulon_de_novo_v1/collectri_viper_v1")
CHECKPOINTS <- file.path(OUT, "endpoint_checkpoints_v1")
dir.create(CHECKPOINTS, recursive=TRUE, showWarnings=FALSE)
SEED <- 26210L

if (file.exists(file.path(OUT, "ACTIVITY_COMPLETE"))) {
  cat("COLLECTRI_VIPER_ACTIVITY_V1_REUSED\n")
  quit(status=0)
}
genes <- fread(file.path(INPUT, "genes_v1.tsv"), colClasses="character")$gene
cells <- fread(file.path(INPUT, "cells_v1.tsv"), colClasses="character")$cell_barcode
meta <- fread(file.path(INPUT, "cell_metadata_v1.tsv.gz"), colClasses=c(cell_barcode="character", sample="character"))
if (!identical(meta$cell_barcode, cells) || nrow(meta) != 39912L) stop("Input metadata mismatch")
con <- gzfile(file.path(INPUT, "raw_counts_genes_by_cells_v1.mtx.gz"), open="rb")
counts <- as(readMM(con), "dgCMatrix"); close(con)
dimnames(counts) <- list(genes, cells)
if (!identical(dim(counts), c(length(genes), length(cells)))) stop("Input matrix mismatch")
lib_size <- Matrix::colSums(counts)
log_cp10k <- counts %*% Diagonal(x=1e4/lib_size)
log_cp10k@x <- log1p(log_cp10k@x)
dimnames(log_cp10k) <- list(genes, cells)
rm(counts); gc(verbose=FALSE)

network <- fread(NETWORK_FILE)
network <- network[target %chin% genes]
network_qc <- data.table(
  resource="CollecTRI", n_edges_input=nrow(fread(NETWORK_FILE)), n_edges_observed=nrow(network),
  n_sources_observed=uniqueN(network$source), n_targets_observed=uniqueN(network$target),
  minsize=5L, old_part1_network_used=FALSE
)
if (network_qc$n_sources_observed < 500L) stop("Insufficient CollecTRI coverage")
fwrite(network_qc, file.path(OUT, "collectri_network_coverage_v1.tsv"), sep="\t")

endpoint_specs <- list(
  tumor_f1=list(types="Tumor", exposure="f1_score"),
  myeloid_dc_tymp=list(types=c("Myeloid", "Dendritic"), exposure="TYMP_log1p"),
  treg_tymp=list(types="Treg", exposure="TYMP_log1p")
)
residualize <- function(y, s, g) {
  design <- cbind(1, s, g)
  y - design %*% solve(crossprod(design), crossprod(design, y))
}

effect_files <- character(); run_rows <- list(); run_counter <- 0L
samples <- sort(unique(meta$sample))
for (sample_id in samples) {
  for (endpoint in names(endpoint_specs)) {
    spec <- endpoint_specs[[endpoint]]
    one_meta <- meta[sample == sample_id & cell_type %chin% spec$types]
    if (nrow(one_meta) < 40L || sd(one_meta[[spec$exposure]]) == 0) next
    checkpoint <- file.path(CHECKPOINTS, paste0(sample_id, "__", endpoint, ".tsv.gz"))
    effect_files <- c(effect_files, checkpoint)
    started <- proc.time()[["elapsed"]]
    if (file.exists(checkpoint)) {
      result_effect <- fread(checkpoint); reused <- TRUE
    } else {
      expression_dense <- as.matrix(log_cp10k[, one_meta$cell_barcode, drop=FALSE])
      set.seed(SEED + match(sample_id, samples))
      activity_long <- as.data.table(decoupleR::run_viper(
        mat=expression_dense, network=network, .source="source", .target="target",
        .mor="mor", .likelihood="likelihood", verbose=FALSE, minsize=5,
        pleiotropy=TRUE, eset.filter=FALSE, cores=2
      ))
      rm(expression_dense); gc(verbose=FALSE)
      if (!all(c("source", "condition", "score") %in% names(activity_long))) stop("Unexpected VIPER output")
      activity_wide <- dcast(activity_long, condition ~ source, value.var="score")
      order_idx <- match(one_meta$cell_barcode, activity_wide$condition)
      if (anyNA(order_idx)) stop("VIPER cell order mismatch: ", sample_id, " ", endpoint)
      activity <- as.matrix(activity_wide[order_idx, -1, with=FALSE])
      colnames(activity) <- names(activity_wide)[-1]
      if (any(!is.finite(activity))) stop("Nonfinite VIPER activity")
      exposure <- as.numeric(one_meta[[spec$exposure]])
      exposure_resid <- as.numeric(residualize(matrix(exposure, ncol=1), one_meta$S_score, one_meta$G2M_score))
      activity_resid <- residualize(activity, one_meta$S_score, one_meta$G2M_score)
      continuous <- as.numeric(cor(activity_resid, exposure_resid, method="pearson"))
      names(continuous) <- colnames(activity_resid)
      lower <- quantile(exposure_resid, 0.25, type=7)
      upper <- quantile(exposure_resid, 0.75, type=7)
      low <- exposure_resid <= lower; high <- exposure_resid >= upper
      contrast <- colMeans(activity_resid[high, , drop=FALSE]) - colMeans(activity_resid[low, , drop=FALSE])
      result_effect <- rbind(
        data.table(sample=sample_id, endpoint=endpoint, metric="cycle_adjusted_continuous",
                   source=names(continuous), effect=continuous, n_cells=nrow(one_meta),
                   n_low=sum(low), n_high=sum(high)),
        data.table(sample=sample_id, endpoint=endpoint, metric="cycle_adjusted_high_minus_low",
                   source=names(contrast), effect=as.numeric(contrast), n_cells=nrow(one_meta),
                   n_low=sum(low), n_high=sum(high))
      )
      fwrite(result_effect, checkpoint, sep="\t"); reused <- FALSE
      rm(activity_long, activity_wide, activity, activity_resid); gc(verbose=FALSE)
    }
    run_counter <- run_counter + 1L
    run_rows[[run_counter]] <- data.table(
      sample=sample_id, endpoint=endpoint, n_cells=nrow(one_meta),
      n_sources=uniqueN(result_effect$source), elapsed_seconds=proc.time()[["elapsed"]]-started,
      checkpoint_reused=reused
    )
    fwrite(rbindlist(run_rows), file.path(OUT, "collectri_viper_run_qc_v1.tsv"), sep="\t")
    cat(sprintf("[CollecTRI-VIPER] %s %s cells=%d sources=%d reused=%s\n",
                sample_id, endpoint, nrow(one_meta), uniqueN(result_effect$source), reused))
  }
}
effects <- rbindlist(lapply(unique(effect_files), fread), use.names=TRUE)
if (uniqueN(effects$sample) != 19L || uniqueN(effects$source) < 500L) stop("Incomplete CollecTRI-VIPER effects")
fwrite(effects, file.path(OUT, "collectri_viper_effects_by_sample_v1.tsv.gz"), sep="\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info_v1.txt"))
writeLines(format(Sys.time(), tz="", usetz=TRUE), file.path(OUT, "ACTIVITY_COMPLETE"))
cat(sprintf("COLLECTRI_VIPER_ACTIVITY_V1_PASS rows=%d sources=%d\n", nrow(effects), uniqueN(effects$source)))
