p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase317_RCC_external_paper_relation_validation_20261005_v1"
r<-dirname(p)
stopifnot(system("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")==0)
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(jsonlite))
f<-fromJSON(file.path(p,"resources/FEATURES_FROM_PAPER.json"))
mods<-c('FAP_CAF25','SPP1_MAC25','FAP_CAF25_disjoint','SPP1_MAC25_disjoint','ECM_org_heldout','Collagen_org_heldout','TGFb_response_heldout')
featurelist<-f[mods]
score_fun<-Seurat::AddModuleScore
original<-deparse(body(score_fun));loc<-which(trimws(original)=='return(object)');stopifnot(length(loc)==1)
original[loc]<-"attr(object, 'phase317_controls') <- ctrl.use; attr(object, 'phase317_features') <- features; return(object)"
body(score_fun)<-parse(text=paste(original,collapse='\n'))[[1]]
writeLines(deparse(score_fun),file.path(p,'code/AddModuleScore_with_audit.R'))
x<-readRDS(file.path(r,'phase316_same_paper_repository_recovery_20261004_v1/downloads/Example_Data_Human_CRC_ST_data.rds'))
x<-score_fun(x,features=featurelist,nbin=24,ctrl=100,seed=1,name='P317',assay='Spatial',search=FALSE)
m<-as.data.frame(x@meta.data[,paste0('P317',seq_along(mods)),drop=FALSE]);names(m)<-mods
co<-x@images[[1]]@coordinates;sc<-x@images[[1]]@scale.factors$lowres
out<-cbind(data.frame(spot_id=rownames(m),imagerow=co$imagerow*sc,imagecol=co$imagecol*sc),m)
write.table(out,file.path(p,'tables/P7_REFERENCE_CALIBRATION_SCORES.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
ctl<-attr(x,'phase317_controls');names(ctl)<-mods;write_json(ctl,file.path(p,'resources/P7_REFERENCE_CONTROLS.json'),pretty=TRUE,auto_unbox=TRUE)
base<-file.path(r,'phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1/tables')
a<-read.delim(file.path(base,'SOURCE_Figure_6g.tsv'),check.names=FALSE);b<-read.delim(file.path(base,'SOURCE_Figure_6h.tsv'),check.names=FALSE)
key<-function(a,b)paste(round(a,8),round(b,8),sep='_');j<-match(key(out$imagerow,out$imagecol),key(a$imagerow,a$imagecol));k<-match(key(out$imagerow,out$imagecol),key(b$imagerow,b$imagecol));stopifnot(!anyNA(j),!anyNA(k))
vals<-data.frame(program=mods[1:2],pearson=c(cor(out$FAP_CAF25,a[j,3]),cor(out$SPP1_MAC25,b[k,3])),spearman=c(cor(out$FAP_CAF25,a[j,3],method='spearman'),cor(out$SPP1_MAC25,b[k,3],method='spearman')))
vals$pass<-vals$pearson>=.90 & vals$spearman>=.90
write.table(vals,file.path(p,'qc/REFERENCE_CALIBRATION.tsv'),sep='\t',quote=FALSE,row.names=FALSE);print(as.data.frame(vals));cat('ALL_PASS',all(vals$pass),'\n')
