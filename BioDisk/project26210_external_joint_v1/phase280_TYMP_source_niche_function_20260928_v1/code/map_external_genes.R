p <- commandArgs(trailingOnly=TRUE)[1]
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Hs.eg.db))
ids <- read.delim(file.path(p,"resources/GSE171542/ENSEMBL_IDS.tsv"),stringsAsFactors=FALSE)$gene_id
valid <- intersect(ids, keys(org.Hs.eg.db,keytype="ENSEMBL"))
mapped <- AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype="ENSEMBL",columns="SYMBOL")
write.table(as.data.frame(mapped),file.path(p,"resources/GSE171542/ENSEMBL_SYMBOL_MAPPING.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
writeLines(capture.output(sessionInfo()),file.path(p,"qc/ANNOTATION_R_SESSION.txt"))
cat("ANNOTATION_MAPPING_COMPLETE",nrow(mapped),"\n")
