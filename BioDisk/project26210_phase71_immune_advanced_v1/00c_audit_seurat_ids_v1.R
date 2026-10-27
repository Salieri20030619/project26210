#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat); library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4L) stop('Usage: dataset rds metadata output')
dataset <- args[[1L]]; obj <- readRDS(args[[2L]]); meta <- fread(args[[3L]]); output <- args[[4L]]
object_ids <- colnames(obj)
if (dataset=='GSE207493') {
  source_ids <- meta$cell_barcode
} else {
  source_ids <- meta$seurat_cell_name
}
candidates <- list(
  original=source_ids,
  add_dataset=paste0(dataset,'__',source_ids),
  remove_dataset_once=sub(paste0('^',dataset,'__'),'',source_ids),
  remove_dataset_twice=sub(paste0('^',dataset,'__',dataset,'__'),paste0(dataset,'__'),source_ids),
  gsm_barcode=paste0(sub('.*(GSM[0-9]+).*','\\1',source_ids),'_',sub('.*__','',source_ids)),
  sample_barcode=paste0(sub('.*__([^_]+)__[^_]+$','\\1',source_ids),'_',sub('.*__','',source_ids))
)
rows <- rbindlist(lapply(names(candidates),function(name){x<-candidates[[name]]; data.table(dataset=dataset,candidate=name,n_source=length(x),n_unique=uniqueN(x),n_matched=sum(x %chin% object_ids),n_missing=sum(!x %chin% object_ids),example_source=source_ids[1],example_candidate=x[1],example_object=object_ids[1])}))
fwrite(rows,output,sep='\t',quote=TRUE)
cat(sprintf('{"dataset":"%s","object_cells":%d,"best_match":%d}\n',dataset,length(object_ids),max(rows$n_matched)))
