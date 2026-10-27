#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('Usage: 00_preflight_phase73_v1.R <output_tsv>')
output_tsv <- args[[1L]]

base71 <- '/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1'
base64 <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1'
base67 <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1'

cfg <- list(
  GSE207493 = list(
    rds=file.path(base71,'inputs_v1/seurat/GSE207493_annotated_local_v2.rds'),
    meta=file.path(base64,'results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz')),
  GSE264586 = list(
    rds=file.path(base71,'inputs_v1/seurat/GSE264586_annotated_local_v2.rds'),
    meta=file.path(base67,'inputs_v1/metadata_v2/GSE264586_phase67_metadata_v1.tsv.gz')),
  GSE294109 = list(
    rds=file.path(base71,'inputs_v1/seurat/GSE294109_annotated_local_v2.rds'),
    meta=file.path(base67,'inputs_v1/metadata_v2/GSE294109_phase67_metadata_v1.tsv.gz'))
)

targets <- c('TYMP','STAT3','SPP1','MTAP')
rows <- list()
for (dataset in names(cfg)) {
  message('PREFLIGHT ',dataset)
  x <- cfg[[dataset]]
  stopifnot(file.exists(x$rds),file.exists(x$meta))
  metadata <- fread(x$meta)
  object <- readRDS(x$rds)
  stopifnot(inherits(object,'Seurat'),'RNA' %in% Assays(object))
  layers <- Layers(object[['RNA']])
  count_layers <- layers[grepl('^counts($|\\.)',layers)]
  if (length(count_layers) != 1L) stop(dataset,': expected exactly one counts layer, found ',paste(count_layers,collapse=';'))
  if (dataset == 'GSE207493') {
    metadata <- metadata[cell_type %in% c('Myeloid','Dendritic')]
    metadata[,seurat_cell_name:=paste0('GSE207493__GSE207493__',sub('^([^_]+)_','\\1__',cell_barcode))]
  } else {
    metadata <- metadata[phase7_broad_annotation=='Myeloid']
  }
  if (anyDuplicated(metadata$seurat_cell_name)) stop(dataset,': duplicate selected Seurat IDs')
  missing <- setdiff(metadata$seurat_cell_name,colnames(object))
  if (length(missing)) stop(dataset,': missing selected cells=',length(missing))
  counts <- LayerData(object[['RNA']],layer=count_layers[[1L]])[,metadata$seurat_cell_name,drop=FALSE]
  upper <- toupper(rownames(counts))
  collision <- upper %in% unique(upper[duplicated(upper)|duplicated(upper,fromLast=TRUE)])
  counts <- counts[!collision,,drop=FALSE]
  rownames(counts) <- upper[!collision]
  target_present <- targets %in% rownames(counts)
  target_detected <- vapply(targets,function(g) if (g %in% rownames(counts)) sum(counts[g,]!=0) else 0L,integer(1L))
  rows[[length(rows)+1L]] <- data.table(
    dataset=dataset,
    source_cells=ncol(object),
    selected_myeloid_dc_cells=ncol(counts),
    genes=nrow(counts),
    layers=paste(layers,collapse=';'),
    counts_layer=count_layers[[1L]],
    target=targets,
    target_present=target_present,
    target_detected_cells=target_detected,
    min_library_size=min(Matrix::colSums(counts)),
    median_library_size=median(Matrix::colSums(counts)),
    max_library_size=max(Matrix::colSums(counts)))
  rm(object,counts,metadata); invisible(gc())
}

out <- rbindlist(rows)
dir.create(dirname(output_tsv),recursive=TRUE,showWarnings=FALSE)
fwrite(out,output_tsv,sep='\t',quote=TRUE,na='NA')
cat(sprintf('PHASE73_PREFLIGHT_PASS rows=%d datasets=%d\n',nrow(out),uniqueN(out$dataset)))
