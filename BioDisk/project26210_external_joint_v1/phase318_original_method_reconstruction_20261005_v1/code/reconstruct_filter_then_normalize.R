suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(jsonlite))
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase318_original_method_reconstruction_20261005_v1";r<-dirname(p)
stopifnot(system("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")==0)
env<-new.env(parent=asNamespace("Seurat"));sys.source(file.path(p,"resources/AddModuleScore_3.2.1.R"),envir=env);score<-env$AddModuleScore
x<-readRDS(file.path(r,"phase316_same_paper_repository_recovery_20261004_v1/downloads/Example_Data_Human_CRC_ST_data.rds"))
f<-fromJSON(file.path(r,"phase317_RCC_external_paper_relation_validation_20261005_v1/resources/FEATURES_FROM_PAPER.json"))[c("FAP_CAF25","SPP1_MAC25")]
ct<-GetAssayData(x,slot="counts");keep<-Matrix::rowSums(ct)>=10 & Matrix::rowSums(ct>0)>=3
co<-x@images[[1]]@coordinates;sc<-x@images[[1]]@scale.factors$lowres
key<-function(a,b)paste(round(a,8),round(b,8),sep="_")
base<-file.path(r,"phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1/tables")
a<-read.delim(file.path(base,"SOURCE_Figure_6g.tsv"),check.names=FALSE);b<-read.delim(file.path(base,"SOURCE_Figure_6h.tsv"),check.names=FALSE)
j<-match(key(co$imagerow*sc,co$imagecol*sc),key(a$imagerow,a$imagecol));k<-match(key(co$imagerow*sc,co$imagecol*sc),key(b$imagerow,b$imagecol));stopifnot(!anyNA(j),!anyNA(k))
res<-list()
for(v in c("paper_gene_filter_then_LogNormalize")) {
 y<-if(v=="existing_universe")x else subset(x,features=rownames(x)[keep])
 cat("POST_FILTER_SPOTS_BELOW_200",sum(Matrix::colSums(GetAssayData(y,slot="counts")>0)<200),"\n");y<-subset(y,cells=colnames(y)[Matrix::colSums(GetAssayData(y,slot="counts")>0)>=200]);y<-NormalizeData(y,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
 y<-score(y,features=f,assay="Spatial",name="P318")
 out<-data.frame(spot_id=colnames(y),FAP=y$P3181,SPP1=y$P3182)
 write.table(out,file.path(p,"tables",paste0("P7_",v,"_scores.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 for(g in 1:2){target<-if(g==1)a[j[match(out$spot_id,colnames(x))],3] else b[k[match(out$spot_id,colnames(x))],3];z<-out[,g+1];res[[length(res)+1]]<-data.frame(variant=v,program=names(f)[g],spots=ncol(y),genes=nrow(y),pearson=cor(z,target),spearman=cor(z,target,method="spearman"),max_abs_error=max(abs(z-target)))}
}
res<-do.call(rbind,res);write.table(res,file.path(p,"tables/P7_PAPER_QC_NORMALIZATION_CALIBRATION.tsv"),sep="\t",quote=FALSE,row.names=FALSE);print(as.data.frame(res))
write_json(list(Seurat=as.character(packageVersion("Seurat")),R=R.version.string,scoring_function="verbatim Seurat v3.2.1",min_features=min(Matrix::colSums(ct>0)),genes_total=nrow(x),paper_filter_genes=sum(keep),original_commands=lapply(x@commands,function(z)z@params)),file.path(p,"qc/P7_PAPER_QC_NORMALIZATION_PROVENANCE.json"),pretty=TRUE,auto_unbox=TRUE)
