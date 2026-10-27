#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(monocle3)
  library(SingleCellExperiment)
  library(jsonlite)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) stop("Missing argument: ", flag)
  args[[pos + 1L]]
}
input_root <- normalizePath(get_arg("--input"), mustWork=TRUE)
output_root <- get_arg("--output")
threads <- as.integer(get_arg("--threads"))
if (dir.exists(output_root) || file.exists(output_root)) stop("Refusing to overwrite output: ", output_root)
dir.create(output_root, recursive=TRUE)
set.seed(20260814)
options(mc.cores=threads)

scale01 <- function(x) {
  out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (!any(ok)) return(out)
  r <- range(x[ok]); out[ok] <- if (diff(r) == 0) 0.5 else (x[ok] - r[1]) / diff(r)
  out
}

study_dirs <- sort(list.dirs(input_root, recursive=FALSE, full.names=TRUE))
all_cells <- list(); roots <- list(); qc <- list()
for (study_dir in study_dirs) {
  study <- basename(study_dir)
  message("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] ", study)
  genes <- fread(file.path(study_dir, "genes_v1.tsv"), header=FALSE)$V1
  md <- fread(file.path(study_dir, "cell_metadata_v1.tsv.gz"))
  counts <- as(readMM(gzfile(file.path(study_dir, "raw_counts_genes_by_cells_v1.mtx.gz"))), "dgCMatrix")
  latent_dt <- fread(file.path(study_dir, "scvi_latent_v1.tsv.gz"), check.names=FALSE)
  cells <- as.character(md$cell_id)
  stopifnot(nrow(counts) == length(genes), ncol(counts) == length(cells), identical(as.character(latent_dt[[1]]), cells))
  rownames(counts) <- genes; colnames(counts) <- cells
  latent <- as.matrix(latent_dt[, -1]); rownames(latent) <- cells
  cell_md <- as.data.frame(md); rownames(cell_md) <- cells
  gene_md <- data.frame(gene_short_name=genes, row.names=genes)
  cds <- new_cell_data_set(counts, cell_metadata=cell_md, gene_metadata=gene_md)
  reducedDims(cds)$PCA <- latent[, seq_len(min(30L, ncol(latent))), drop=FALSE]
  cds <- reduce_dimension(cds, reduction_method="UMAP", preprocess_method="PCA", umap.metric="cosine", umap.min_dist=0.1, umap.n_neighbors=30L, cores=threads, verbose=TRUE)
  cds <- cluster_cells(cds, reduction_method="UMAP", k=30L, cluster_method="leiden", num_iter=2L, random_seed=20260814L, verbose=TRUE)
  cds <- learn_graph(cds, use_partition=FALSE, close_loop=TRUE, learn_graph_control=list(minimal_branch_len=10L, euclidean_distance_ratio=1), verbose=TRUE)
  mono <- which(md$state_v2 == "Monocyte")
  if (length(mono) < 50L) stop(study, ": fewer than 50 Monocytes")
  center <- apply(latent[mono, , drop=FALSE], 2L, median)
  root_idx <- mono[which.min(rowSums((latent[mono, , drop=FALSE] - matrix(center, nrow=length(mono), ncol=length(center), byrow=TRUE))^2))]
  root_cell <- cells[root_idx]
  cds <- order_cells(cds, reduction_method="UMAP", root_cells=root_cell)
  pt_raw <- pseudotime(cds, reduction_method="UMAP")[cells]
  pt <- scale01(pt_raw)
  finite_fraction <- mean(is.finite(pt))
  if (finite_fraction < 0.95) stop(study, ": finite pseudotime fraction below 0.95: ", finite_fraction)
  cell_out <- data.table(
    cell_id=cells,
    study=study,
    patient=md$patient,
    state_v2=md$state_v2,
    monocle3_pseudotime_v1=pt,
    monocle3_cluster_v1=as.character(clusters(cds, reduction_method="UMAP")[cells]),
    monocle3_partition_v1=as.character(partitions(cds, reduction_method="UMAP")[cells])
  )
  all_cells[[study]] <- cell_out
  roots[[study]] <- data.table(study=study, root_cell=root_cell, root_state=md$state_v2[root_idx], n_cells=length(cells), finite_fraction=finite_fraction)
  qc[[study]] <- data.table(study=study, n_cells=length(cells), n_patients=uniqueN(md$patient), n_clusters=uniqueN(cell_out$monocle3_cluster_v1), n_partitions=uniqueN(cell_out$monocle3_partition_v1), n_finite=sum(is.finite(pt)))
  saveRDS(cds, file.path(output_root, paste0(study, "_monocle3_cds_v1.rds")), compress="xz")
  rm(cds, counts, latent); gc()
}

fwrite(rbindlist(all_cells), file.path(output_root, "monocle3_cell_pseudotime_v1.tsv.gz"), sep="\t")
fwrite(rbindlist(roots), file.path(output_root, "monocle3_root_qc_v1.tsv"), sep="\t")
fwrite(rbindlist(qc), file.path(output_root, "monocle3_study_qc_v1.tsv"), sep="\t")
writeLines(capture.output(sessionInfo()), file.path(output_root, "R_session_info_v1.txt"))
status <- list(
  status="PASS_Y02_MONOCLE3_COMPUTE",
  completed_at=format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_cells=sum(vapply(all_cells, nrow, integer(1))),
  n_studies=length(all_cells),
  n_patients=uniqueN(rbindlist(all_cells)$patient),
  study_separate=TRUE,
  root_rule="Monocyte cell nearest the within-study Monocyte median in X_scVI",
  graph_method="Monocle3 UMAP plus Leiden plus learn_graph; use_partition=FALSE",
  interpretation="cross-sectional state ordering; not lineage tracing, RNA velocity, longitudinal time or causality",
  figures_created=FALSE
)
write_json(status, file.path(output_root, "status.json"), auto_unbox=TRUE, pretty=TRUE)
cat("Y02_MONOCLE3_COMPUTE_PASS\n")
