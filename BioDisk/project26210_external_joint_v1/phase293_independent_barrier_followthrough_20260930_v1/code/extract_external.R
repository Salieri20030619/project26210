suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite)})
H<-'/mnt/BioDisk/wangmingyuan';R<-file.path(H,'project26210_external_joint_v1')
P<-file.path(R,'phase293_independent_barrier_followthrough_20260930_v1')
gate<-function(stage){
 out<-system2('bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh',stdout=TRUE,stderr=TRUE)
 cat(stage,'\n',paste(out,collapse='\n'),'\n',file=file.path(P,'logs/extraction_memory.log'),append=TRUE)
 if(!is.null(attr(out,'status')) && attr(out,'status')!=0) stop('MEMORY_GATE')
}
id<-fromJSON(file.path(R,'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json'))
prog<-fromJSON(file.path(R,'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json'))
genes<-sort(unique(c(unlist(id),unlist(prog),'TYMP','SPP1','CD68','C1QA','C1QB','C1QC','CSF1R','FCN1','S100A8','S100A9','CD8A','CD8B','TRAC')))
for(ds in c('GSE264586','GSE294109')){
 gate(paste0('load_',ds))
 input<-file.path(H,'project26210_phase71_immune_advanced_v1/inputs_v1/seurat',paste0(ds,'_annotated_local_v2.rds'))
 before<-file.info(input)[,c('size','mtime')]
 obj<-readRDS(input)
 m<-fread(file.path(H,'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2',paste0(ds,'_phase67_metadata_v1.tsv.gz')))
 stopifnot(!anyDuplicated(m$cell_id),!anyDuplicated(m$seurat_cell_name),all(m$seurat_cell_name %in% rownames(obj@meta.data)))
 idx<-match(m$seurat_cell_name,rownames(obj@meta.data))
 stopifnot(all(m$phase7_fine_annotation==obj@meta.data$phase7_fine_annotation[idx]))
 x<-GetAssayData(obj,assay='RNA',layer='counts')
 stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
 x<-x[,m$seurat_cell_name,drop=FALSE]
 stopifnot(all(is.finite(x@x)),all(x@x>=0),all(x@x==floor(x@x)))
 total<-Matrix::colSums(x);ng<-Matrix::colSums(x>0)
 m[,total_counts:=as.numeric(total)];m[,n_genes:=as.integer(ng)]
 for(field in names(obj@meta.data)[grepl('tissue|organ|site|hist|tumor|normal',names(obj@meta.data),ignore.case=TRUE)]){
   m[[paste0('original_',field)]]<-as.character(obj@meta.data[[field]][idx])
 }
 g<-intersect(genes,rownames(x))
 y<-as.data.table(as.matrix(t(x[g,,drop=FALSE])));y[,cell_id:=m$cell_id];setcolorder(y,c('cell_id',g))
 fwrite(m,file.path(P,'inputs',paste0(ds,'_METADATA.tsv.gz')),sep='\t')
 fwrite(y,file.path(P,'inputs',paste0(ds,'_PANEL_COUNTS.tsv.gz')),sep='\t')
 fwrite(data.table(gene=genes,present=genes %in% g),file.path(P,'tables',paste0(ds,'_GENE_COVERAGE.tsv')),sep='\t')
 stopifnot(identical(before,file.info(input)[,c('size','mtime')]))
 cat('EXTRACTED',ds,'cells',nrow(m),'genes',length(g),'UMI',sum(total),'\n')
 rm(obj,x,m,y,total,ng);gc()
}
writeLines(capture.output(sessionInfo()),file.path(P,'resources/EXTRACTION_R_SESSION.txt'))
writeLines(format(Sys.time(),tz='Asia/Shanghai'),file.path(P,'qc/EXTRACTION_COMPLETE'))
