p <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase316_same_paper_repository_recovery_20261004_v1"
stopifnot(system("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")==0)
suppressPackageStartupMessages(library(SeuratObject))
x <- readRDS(file.path(p,"downloads/Example_Data_Human_CRC_ST_data.rds"))
im <- x@images[[1]];a<-x@assays[[1]]
cat("SCALE\n");print(im@scale.factors);cat("ASSAY_CLASS",class(a),"SLOTS",slotNames(a),"\n")
co<-as.data.frame(im@coordinates);m<-as.data.frame(x@meta.data)
stopifnot(identical(rownames(co),colnames(x)),identical(rownames(m),colnames(x)))
low<-im@scale.factors$lowres;hi<-im@scale.factors$hires
co$source_imagerow<-co$imagerow*low;co$source_imagecol<-co$imagecol*low
write.table(cbind(barcode=rownames(co),co),file.path(p,"tables/CANDIDATE_SCALED_COORDINATES.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
genes<-c('CD3D','CD8A','CD4','CD19','MS4A1','GZMA','PRF1','GNLY','FAP','SPP1','COL1A1','COL1A2','EPCAM')
stopifnot(all(genes %in% rownames(a@counts)),all(genes %in% rownames(a@data)))
counts<-a@counts;data<-a@data
cat("COUNTS",dim(counts),"DATA",dim(data),"NONZERO",length(counts@x),"\n")
cm<-Matrix::colSums(counts);nf<-Matrix::colSums(counts>0)
qc<-data.frame(barcode=colnames(x),counts_library_size=cm,counts_detected_features=nf,meta_library_size=m$nCount_Spatial,meta_detected_features=m$nFeature_Spatial,original_label=m$DefineTypes)
write.table(qc,file.path(p,"tables/CANDIDATE_COUNT_QC.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for (kind in c('counts','data')) {
 z<-if(kind=='counts')counts else data
 vals<-as.data.frame(t(as.matrix(z[genes,,drop=FALSE])))
 write.table(cbind(barcode=rownames(vals),vals),file.path(p,paste0("tables/CANDIDATE_",toupper(kind),"_TARGETS.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
}
writeLines(c(paste('lowres',low),paste('hires',hi),paste('R',R.version.string),paste('SeuratObject',as.character(packageVersion('SeuratObject')))),file.path(p,"resources/CANDIDATE_OBJECT_SCALE_AND_ENV.txt"))
cat("EXTRACT_COMPLETE\n")
