suppressPackageStartupMessages({library(Seurat);library(data.table);library(Matrix)})
H<-'/mnt/BioDisk/wangmingyuan';O<-file.path(H,'project26210_external_joint_v1/phase211_immune_receiver_refreeze_20260919_v1');dir.create(file.path(O,'external'),showWarnings=FALSE)
genes<-c('CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','CD4','IL7R','CD8A','CD8B','FOXP3','IKZF2','IL2RA','CTLA4','KLRD1','NCR1','KLRF1','NCAM1','FCGR3A','TRDC','TRGC1','TRGC2','NKG7','GNLY','GZMB','PRF1','CTSW','PDCD1','LAG3','HAVCR2','TIGIT','TOX','CXCL13','ENTPD1','SPP1')
for(ds in c('GSE264586','GSE294109')){
 rc<-system('/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh');if(rc!=0)stop('memory gate')
 p<-file.path(H,'project26210_phase71_immune_advanced_v1/inputs_v1/seurat',paste0(ds,'_annotated_local_v2.rds'));before<-file.info(p)[,c('size','mtime')];obj<-readRDS(p);md<-obj@meta.data
 cat('DATASET',ds,'METADATA_FIELDS',paste(names(md),collapse='|'),'\n')
 frozen<-fread(file.path(H,'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2',paste0(ds,'_phase67_metadata_v1.tsv.gz')))
 stopifnot(!anyDuplicated(frozen$cell_id),!anyDuplicated(frozen$seurat_cell_name),all(frozen$seurat_cell_name %in% rownames(md)))
 ix<-match(frozen$seurat_cell_name,rownames(md));needed<-c('phase7_cluster_res0_8','phase7_broad_annotation','phase7_fine_annotation','phase7_annotation_state','phase7_annotation_confidence')
 stopifnot(all(needed %in% names(md)),all(as.character(md$phase7_fine_annotation[ix])==frozen$phase7_fine_annotation))
 for(n in needed)frozen[[paste0('rds_',n)]]<-as.character(md[[n]][ix])
 tissuefields<-names(md)[grepl('tissue|organ|site|tumor|normal|disease',names(md),ignore.case=TRUE)]
 for(n in tissuefields)frozen[[paste0('rds_source_',n)]]<-as.character(md[[n]][ix])
 fwrite(frozen,file.path(O,'external',paste0(ds,'_CELL_MAPPING.tsv.gz')),sep='\t')
 x<-GetAssayData(obj,assay='RNA',layer='counts');stopifnot(all(frozen$seurat_cell_name %in% colnames(x)),!anyDuplicated(rownames(x)))
 g<-intersect(genes,rownames(x));sel<-frozen$phase7_broad_annotation %in% c('T_cell','NK_cell');ids<-frozen$seurat_cell_name[sel]
 panel<-as.data.table(as.matrix(t(x[g,ids,drop=FALSE])));panel[,cell_id:=frozen$cell_id[sel]];setcolorder(panel,c('cell_id',g));fwrite(panel,file.path(O,'external',paste0(ds,'_TNK_PANEL_RAW_COUNTS.tsv.gz')),sep='\t')
 writeLines(paste(setdiff(genes,g),collapse=','),file.path(O,'external',paste0(ds,'_MISSING_GENES.txt')))
 stopifnot(identical(before,file.info(p)[,c('size','mtime')]))
 cat('EXTERNAL_COMPLETE',ds,'mapped',nrow(frozen),'TNK',sum(sel),'genes',length(g),'\n');rm(obj,x,md,frozen,panel);gc()
}
writeLines(capture.output(sessionInfo()),file.path(O,'external/R_SESSION_INFO.txt'))
writeLines('PASS',file.path(O,'external/EXTRACTION_COMPLETE'))
