suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table)})
H<-'/mnt/BioDisk/wangmingyuan';R<-file.path(H,'project26210_external_joint_v1');O<-file.path(R,'phase213_immune_restart_execution_20260919_v1');D<-file.path(O,'I5_external_description');if(dir.exists(D))stop('refuse overwrite');dir.create(D)
F<-file.path(R,'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2');P<-file.path(R,'phase211_immune_receiver_refreeze_20260919_v1/external')
rc<-system('/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh');if(rc!=0)stop('memory gate')
ds<-'GSE294109';path<-file.path(H,'project26210_phase71_immune_advanced_v1/inputs_v1/seurat/GSE294109_annotated_local_v2.rds');before<-file.info(path)[,c('size','mtime')];obj<-readRDS(path);x<-GetAssayData(obj,assay='RNA',layer='counts');total<-Matrix::colSums(x);stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
src<-as.data.frame(fread(file.path(P,'GSE294109_CELL_MAPPING.tsv.gz')));rec<-as.data.frame(fread(file.path(F,'GSE294109_RECEIVER_CELL_FREEZE.tsv.gz')));elig<-as.data.frame(fread(file.path(F,'EXTERNAL_ENDPOINT_SAMPLE_FREEZE.tsv')));elig<-elig[elig$dataset==ds & elig$eligible==TRUE,,drop=FALSE]
stopifnot(!anyDuplicated(src$cell_id),!anyDuplicated(rec$cell_id),all(src$seurat_cell_name %in% colnames(x)),all(rec$cell_id %in% src$cell_id));stopifnot(all(rec$seurat_cell_name==src$seurat_cell_name[match(rec$cell_id,src$cell_id)]))
programs<-list(CD8_cytotoxicity=c('NKG7','GNLY','GZMB','PRF1','CTSW'),CD8_dysfunction=c('PDCD1','CTLA4','LAG3','HAVCR2','TIGIT','TOX','CXCL13','ENTPD1'));stopifnot(all(unique(unlist(programs)) %in% rownames(x)),'SPP1' %in% rownames(x))
rows<-list();coverage<-list()
for(i in seq_len(nrow(elig))){
 e<-elig[i,,drop=FALSE];a<-src[src$patient_key==e$patient & src$sample_id==e$sample_id & src$timepoint==e$timepoint & src$treatment==e$treatment,,drop=FALSE];m<-a[a$phase7_fine_annotation %in% c('Macrophage_C1QC','Macrophage_SPP1'),,drop=FALSE]
 labs<-if(e$endpoint=='Treg_CD4')c('CD4_Tcon','CD4_Treg') else 'CD8_T';b<-rec[rec$patient_key==e$patient & rec$sample_id==e$sample_id & rec$timepoint==e$timepoint & rec$treatment==e$treatment & rec$receiver_label_v1 %in% labs,,drop=FALSE];stopifnot(nrow(m)==e$n_mature,nrow(b)==e$n_receiver,nrow(m)>=50,all(total[m$seurat_cell_name]>0),all(total[b$seurat_cell_name]>0))
 sp<-sum(x['SPP1',m$seurat_cell_name]);den<-sum(total[m$seurat_cell_name]);pred<-log2(1+1e6*sp/den)
 if(e$endpoint=='Treg_CD4'){out<-sum(b$receiver_label_v1=='CD4_Treg')/nrow(b)}else{
  gs<-programs[[e$endpoint]];mat<-as.matrix(x[gs,b$seurat_cell_name,drop=FALSE]);norm<-log1p(sweep(mat,2,total[b$seurat_cell_name],'/')*10000);out<-mean(norm)
  for(g in gs)coverage[[length(coverage)+1L]]<-data.frame(dataset=ds,patient=e$patient,sample_id=e$sample_id,endpoint=e$endpoint,gene=g,n_receiver=nrow(b),detected_cells=sum(mat[g,]>0),mean_log1pCP10K=mean(norm[g,]))
 }
 rows[[length(rows)+1L]]<-data.frame(dataset=ds,patient=e$patient,sample_id=e$sample_id,timepoint=e$timepoint,treatment=e$treatment,endpoint=e$endpoint,n_mature=nrow(m),n_receiver=nrow(b),mature_SPP1_counts=sp,mature_total_counts=den,x_mature_SPP1=pred,y_outcome=out,n_Treg=if(e$endpoint=='Treg_CD4')sum(b$receiver_label_v1=='CD4_Treg') else NA_integer_,C4_used=FALSE)
}
out<-do.call(rbind,rows);fwrite(out,file.path(D,'EXTERNAL_SAMPLE_METRICS.tsv'),sep='\t');fwrite(do.call(rbind,coverage),file.path(D,'EXTERNAL_GENE_COVERAGE.tsv'),sep='\t')
stopifnot(identical(before,file.info(path)[,c('size','mtime')]))
writeLines(capture.output(sessionInfo()),file.path(D,'R_SESSION_INFO.txt'));writeLines('PASS_SAMPLE_METRICS_NO_INFERENCE',file.path(D,'R_EXTRACTION_STATUS.txt'));cat('EXTERNAL_SAMPLE_METRICS_COMPLETE',nrow(out),'rows\n')
