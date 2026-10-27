options(stringsAsFactors=FALSE,width=160)
suppressPackageStartupMessages(library(SeuratObject))
root <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1'
out <- file.path(root,'legacy_object');dir.create(out)
p <- '/mnt/BioDisk/ludingyu/project26210/object/02_macrophage_final_subtypes.rds'
cat('START',format(Sys.time()),'\n');flush.console()
obj <- readRDS(p)
cat('LOADED',format(Sys.time()),'\n');flush.console()
meta <- obj@meta.data
write.table(data.frame(cell_id=rownames(meta),meta,check.names=FALSE),gzfile(file.path(out,'metadata.tsv.gz')) ,sep='\t',quote=FALSE,row.names=FALSE)
summary <- c(paste('class',paste(class(obj),collapse=',')),paste('object_version',as.character(obj@version)),paste('dimensions',paste(dim(obj),collapse='x')),paste('object_size_GB',as.numeric(object.size(obj))/1024^3),paste('assays',paste(names(obj@assays),collapse=',')),paste('reductions',paste(names(obj@reductions),collapse=',')),paste('graphs',paste(names(obj@graphs),collapse=',')),paste('neighbors',paste(names(obj@neighbors),collapse=',')),paste('metadata_columns',paste(colnames(meta),collapse=',')),paste('misc_keys',paste(names(obj@misc),collapse=',')),paste('tools_keys',paste(names(obj@tools),collapse=',')))
writeLines(summary,file.path(out,'OBJECT_SUMMARY.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_SESSION_INFO.txt'))
writeLines(capture.output(str(obj@misc,max.level=3,list.len=40,vec.len=8)),file.path(out,'MISC_STRUCTURE.txt'))
writeLines(capture.output(str(obj@tools,max.level=3,list.len=40,vec.len=8)),file.path(out,'TOOLS_STRUCTURE.txt'))
cmds <- obj@commands
writeLines(capture.output(lapply(cmds,function(z) {list(name=z@name,time=as.character(z@time.stamp),assay=z@assay.used,call=z@call.string,params=z@params)})),file.path(out,'COMMANDS.txt'))
writeLines(capture.output(lapply(obj@reductions,function(z) list(dim=dim(z@cell.embeddings),assay=z@assay.used,key=z@key))),file.path(out,'REDUCTIONS.txt'))
for (nm in names(obj@assays)) {
 a <- obj@assays[[nm]]
 writeLines(capture.output(str(a,max.level=1)),file.path(out,paste0('ASSAY_',nm,'_STRUCTURE.txt')))
}
rows <- list()
for (col in names(meta)) {
 x <- meta[[col]];n <- length(unique(x))
 if (n<=200) {
  t <- as.data.frame(table(x,useNA='ifany'));names(t)<-c('value','n');t$column<-col;rows[[col]]<-t[,c('column','value','n')]
 }
}
write.table(do.call(rbind,rows),file.path(out,'CATEGORICAL_COUNTS.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(cell_id=names(obj@active.ident),active_ident=as.character(obj@active.ident)),gzfile(file.path(out,'active_ident.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
# Small expression summaries verify marker labeling without fitting models or testing cells as replicates.
assay <- if ('RNA' %in% names(obj@assays)) 'RNA' else names(obj@assays)[1]
features <- intersect(c('SPP1','C1QA','C1QB','C1QC','APOE','GPNMB','LGALS3','FABP5','CD9','LPL','CX3CR1','IL1B','CXCL8','HLA-DRA','CD74','CCL3','MKI67','TOP2A','LST1','TYROBP','FCN1','S100A8','S100A9','CD68','CSF1R','TREM2'),rownames(obj[[assay]]))
mat <- tryCatch(GetAssayData(obj,assay=assay,layer='data')[features,,drop=FALSE],error=function(e)NULL)
if (!is.null(mat) && 'macro_subtype' %in% names(meta)) {
 res<-list()
 for (g in unique(as.character(meta$macro_subtype))) {
  ids<-which(as.character(meta$macro_subtype)==g)
  if (!length(ids)) next
  res[[g]]<-data.frame(subtype=g,gene=features,n_cells=length(ids),mean_expression=Matrix::rowMeans(mat[,ids,drop=FALSE]),fraction_positive=Matrix::rowMeans(mat[,ids,drop=FALSE]>0))
 }
 write.table(do.call(rbind,res),file.path(out,'MARKER_DESCRIPTIVE_MEANS.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
cat('COMPLETE',format(Sys.time()),'\n')
writeLines(format(Sys.time()),file.path(out,'COMPLETE.txt'))
