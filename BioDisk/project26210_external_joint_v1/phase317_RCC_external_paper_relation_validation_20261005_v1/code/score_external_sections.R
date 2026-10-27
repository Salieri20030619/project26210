p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase317_RCC_external_paper_relation_validation_20261005_v1"
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(Matrix))
reg<-read.delim(file.path(p,'resources/EXTERNAL_SECTION_REGISTER.tsv'),check.names=FALSE)
f<-fromJSON(file.path(p,'resources/FEATURES_FROM_PAPER.json'));mods<-c('FAP_CAF25','SPP1_MAC25','FAP_CAF25_disjoint','SPP1_MAC25_disjoint','ECM_org_heldout','Collagen_org_heldout','TGFb_response_heldout')
score_fun<-Seurat::AddModuleScore;src<-deparse(body(score_fun));ii<-which(trimws(src)=='return(object)');stopifnot(length(ii)==1);src[ii]<-"attr(object,'controls')<-ctrl.use; attr(object,'used_features')<-features; return(object)";body(score_fun)<-parse(text=paste(src,collapse='\n'))[[1]]
readbin<-function(file,what,n,size) {con<-file(file,'rb');on.exit(close(con));readBin(con,what,n=n,size=size,endian='little')}
for (sample in reg$sample) {
 stopifnot(system('/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh')==0)
 d<-file.path(p,'inputs',sample);inf<-fromJSON(file.path(d,'INPUT_QC.json'));genes<-readLines(file.path(d,'genes.txt'));spots<-readLines(file.path(d,'spots.txt'))
 xx<-readbin(file.path(d,'values.bin'),'double',inf$nnz,8);ii<-readbin(file.path(d,'indices.bin'),'integer',inf$nnz,4);pp<-readbin(file.path(d,'indptr.bin'),'integer',inf$n_spots+1,4)
 mat<-new('dgCMatrix',x=xx,i=ii,p=pp,Dim=as.integer(c(inf$n_genes,inf$n_spots)),Dimnames=list(genes,spots))
 obj<-CreateSeuratObject(counts=mat,min.cells=0,min.features=0);obj<-NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
 featurelist<-lapply(f[mods],function(gs)intersect(gs,rownames(obj)));stopifnot(all(mapply(function(a,b)length(a)>=max(5,ceiling(.8*length(b))),featurelist,f[mods])))
 stopifnot(system('/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh')==0)
 obj<-score_fun(obj,features=featurelist,nbin=24,ctrl=100,seed=1,name='Ext',search=FALSE)
 result<-as.data.frame(obj@meta.data[,paste0('Ext',seq_along(mods)),drop=FALSE]);names(result)<-mods
 write.table(cbind(spot_id=rownames(result),result),file.path(d,'MODULE_SCORES.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 ctl<-attr(obj,'controls');names(ctl)<-mods;used<-attr(obj,'used_features');names(used)<-mods
 write_json(list(controls=ctl,features=used),file.path(d,'SCORING_MEMBERSHIP.json'),pretty=TRUE,auto_unbox=FALSE)
 write.table(data.frame(input_gene=genes,scoring_gene=rownames(obj)),file.path(d,'GENE_NAME_MAPPING.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 cat('SCORED',sample,nrow(result),'\n');flush.console();rm(obj,mat,xx,ii,pp,result);gc()
}
writeLines(capture.output(sessionInfo()),file.path(p,'resources/EXTERNAL_SCORING_R_SESSION.txt'));cat('ALL_SCORED',nrow(reg),'\n')
