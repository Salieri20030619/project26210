#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
obj <- readRDS(args[[1L]]); meta <- fread(args[[2L]])
mapped <- paste0('GSE207493__GSE207493__',sub('^([^_]+)_','\\1__',meta$cell_barcode))
rows <- data.table(n_metadata=nrow(meta),n_unique=uniqueN(mapped),n_matched=sum(mapped %chin% colnames(obj)),n_missing=sum(!mapped %chin% colnames(obj)),example_source=meta$cell_barcode[1],example_mapped=mapped[1],example_object=colnames(obj)[1])
fwrite(rows,args[[3L]],sep='\t',quote=TRUE)
cat(sprintf('{"matched":%d,"missing":%d,"unique":%d}\n',rows$n_matched,rows$n_missing,rows$n_unique))
