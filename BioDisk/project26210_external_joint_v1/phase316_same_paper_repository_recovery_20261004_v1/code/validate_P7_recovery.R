p <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase316_same_paper_repository_recovery_20261004_v1"
b <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1"
stopifnot(system("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")==0)
suppressPackageStartupMessages(library(SeuratObject));suppressPackageStartupMessages(library(Matrix))
x<-readRDS(file.path(p,"downloads/Example_Data_Human_CRC_ST_data.rds"));co<-x@images[[1]]@coordinates;scale<-x@images[[1]]@scale.factors$lowres
labels<-read.delim(file.path(b,"tables/SOURCE_Figure_6e.tsv"),check.names=FALSE)
key<-function(a,b)paste(round(a,8),round(b,8),sep="_")
j<-match(key(co$imagerow*scale,co$imagecol*scale),key(labels$imagerow,labels$imagecol))
checks<-data.frame(check=character(),pass=logical(),detail=character())
ck<-function(n,v,d="") {checks<<-rbind(checks,data.frame(check=n,pass=isTRUE(v),detail=as.character(d)))}
ck("all_coordinates_match",length(j)==4248&&!anyNA(j)&&!anyDuplicated(j));ck("raw_RDS_original_labels",identical(as.character(x@meta.data$DefineTypes),as.character(labels$DefineTypes[j])))
oldqc<-readxl::read_excel(file.path(b,"downloads/source.xlsx"),sheet="Figure 8d-f");qj<-match(key(co$imagerow*scale,co$imagecol*scale),key(oldqc$imageRow,oldqc$imageCol))
cts<-x@assays[[1]]@counts;total<-Matrix::colSums(cts);nf<-Matrix::colSums(cts>0)
ck("raw_UMI_original_QC",all(total==oldqc$nCount_Spatial[qj]));ck("raw_features_original_QC",all(nf==oldqc$nFeature_Spatial[qj]))
stopifnot(system("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")==0)
res<-read.delim(file.path(p,"tables/P7_RECOVERED_IMMUNE_RESULTS.tsv"),check.names=FALSE);lab<-as.character(labels$DefineTypes[j]);new<-list()
for(i in seq_len(nrow(res))) {
 row<-res[i,];values<-log1p(as.numeric(cts[row$gene,])/total*1e4)
 a<-lab=="Malignant epithelial cells";z<-!a
 if(row$comparison=="exclude_unknown")z<-z&!tolower(lab)%in%c("unknown","unkown")
 computed<-c(mean_malignant=mean(values[a]),mean_other=mean(values[z]),delta_mean=mean(values[a])-mean(values[z]),detection_malignant=mean(values[a]>0),detection_other=mean(values[z]>0),delta_detection=mean(values[a]>0)-mean(values[z]>0))
 for(nm in names(computed))ck(paste(row$comparison,row$gene,nm,sep="_"),abs(computed[nm]-as.numeric(row[[nm]]))<1e-12,abs(computed[nm]-as.numeric(row[[nm]])))
 new[[i]]<-cbind(data.frame(comparison=row$comparison,gene=row$gene),as.data.frame(t(computed)))
}
write.table(checks,file.path(p,"qc/R_INDEPENDENT_CHECKS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,new),file.path(p,"qc/R_INDEPENDENT_P7_RESULTS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(p,"resources/VALIDATION_R_SESSION.txt"))
cat("INDEPENDENT_CHECKS",nrow(checks),"FAIL",sum(!checks$pass),"\n");stopifnot(all(checks$pass))
