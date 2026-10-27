suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite)})
H<-'/mnt/BioDisk/wangmingyuan';R<-file.path(H,'project26210_external_joint_v1');P<-file.path(R,'phase293_independent_barrier_followthrough_20260930_v1')
checks<-list();sources<-list()
check<-function(name,value){checks[[length(checks)+1]]<<-data.table(check=name,passed=isTRUE(value));if(!isTRUE(value))stop(name)}
id<-fromJSON(file.path(R,'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json'))
for(ds in c('GSE264586','GSE294109')){
 out<-system2('bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh',stdout=TRUE,stderr=TRUE)
 cat(ds,paste(out,collapse='\n'),'\n',file=file.path(P,'logs/raw_verify_memory.log'),append=TRUE);stopifnot(is.null(attr(out,'status'))||attr(out,'status')==0)
 f<-file.path(H,'project26210_phase71_immune_advanced_v1/inputs_v1/seurat',paste0(ds,'_annotated_local_v2.rds'))
 before<-system2('sha256sum',shQuote(f),stdout=TRUE)
 obj<-readRDS(f);m<-fread(file.path(P,'inputs',paste0(ds,'_METADATA.tsv.gz')))
 y<-fread(file.path(P,'inputs',paste0(ds,'_PANEL_COUNTS.tsv.gz')))
 lab<-fread(file.path(P,'inputs',paste0(ds,'_IDENTITY_FREEZE.tsv.gz')))
 check(paste0(ds,'_id'),identical(m$cell_id,y$cell_id)&&identical(m$cell_id,lab$cell_id))
 x<-GetAssayData(obj,assay='RNA',layer='counts')[,m$seurat_cell_name,drop=FALSE]
 gs<-setdiff(names(y),'cell_id');got<-as.matrix(t(x[gs,,drop=FALSE]))
 check(paste0(ds,'_all_panel_counts'),isTRUE(all.equal(unname(got),unname(as.matrix(y[,..gs])),tolerance=0)))
 check(paste0(ds,'_totals'),all(as.numeric(Matrix::colSums(x))==m$total_counts))
 check(paste0(ds,'_ngenes'),all(as.numeric(Matrix::colSums(x>0))==m$n_genes))
 d<-lapply(id,function(g)rowSums(got[,intersect(g,colnames(got)),drop=FALSE]>0))
 core<-d$fibro_core>=2 & d$fibro_support>=1
 comp<-d$immune>=2|d$endothelial>=2|d$epithelial>=2
 mural<-d$pericyte>=2|(rowSums(got[,c('MYH11','CNN1'),drop=FALSE]>0)>=1 & rowSums(got[,c('ACTA2','TAGLN'),drop=FALSE]>0)>=1)
 expected<-m$phase7_broad_annotation=='Fibroblast'&core&!comp&!mural
 check(paste0(ds,'_strict_members'),all(expected==lab$strict_fibro))
 check(paste0(ds,'_source_unchanged'),identical(before,system2('sha256sum',shQuote(f),stdout=TRUE)))
 sources[[length(sources)+1]]<-list(dataset=ds,path=f,sha256=strsplit(before,' ')[[1]][1],cells=nrow(m),genes=length(gs),values=length(got),strict_cells=sum(expected))
 cat('VERIFIED',ds,'cells',nrow(m),'panel_values',length(got),'strict_fibro',sum(expected),'\n')
 rm(obj,x,m,y,lab,got);gc()
}
fwrite(rbindlist(checks),file.path(P,'qc/RAW_EXTRACTION_CHECKS.tsv'),sep='\t')
write_json(list(status='PASS',checks=length(checks),sources=sources),file.path(P,'qc/RAW_EXTRACTION_VALIDATION.json'),pretty=TRUE,auto_unbox=TRUE)
