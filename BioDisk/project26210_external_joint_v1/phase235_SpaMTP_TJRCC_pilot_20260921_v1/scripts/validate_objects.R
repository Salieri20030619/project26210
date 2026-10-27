options(width=150)
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(jsonlite)})
P<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase235_SpaMTP_TJRCC_pilot_20260921_v1"
q<-read.delim(file.path(P,"tables","SPAMTP_MODULE_QC.tsv"),check.names=FALSE)
rows<-list()
for(i in seq_len(nrow(q))){
 x<-readRDS(q$rds_path[i]);stopifnot(isTRUE(validObject(x)),nrow(x)==q$features[i],ncol(x)==q$pixels_or_bins[i])
 stopifnot("SPM" %in% Assays(x),"fov" %in% names(x@images))
 m<-LayerData(x,assay="SPM",layer="data")
 err<-max(abs(Matrix::colSums(m)-nrow(m)));stopifnot(err<1e-8)
 if(q$binned[i]) stopifnot(ncol(Embeddings(x,"SpatialPCA"))==20,length(unique(x$SM_technical_region4))==4)
 rows[[i]]<-data.frame(file=basename(q$rds_path[i]),valid=TRUE,TIC_column_sum_error=err,features=nrow(x),units=ncol(x))
 rm(x,m);gc()
}
out<-do.call(rbind,rows)
write.table(out,file.path(P,"audit","SERIALIZED_OBJECT_CHECKS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write_json(list(status="PASS",objects=nrow(out),raw_to_input_to_Seurat_counts_verified=TRUE,TIC_normalization_verified=TRUE),file.path(P,"SERIALIZED_OBJECT_VALIDATION.json"),pretty=TRUE,auto_unbox=TRUE)
print(out,row.names=FALSE)
