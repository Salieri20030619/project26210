#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(data.table);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop('Usage: dataset rds metadata cycle output')
dataset<-args[[1]]; rds<-args[[2]]; metadata_path<-args[[3]]; cycle_path<-args[[4]]; out<-args[[5]]; part<-paste0(out,'.part')
if(dir.exists(out)||dir.exists(part))stop('Refusing existing output');dir.create(part,recursive=TRUE)
SEED<-2621072L; MIN_CELLS<-50L; MIN_STRATUM<-30L
GENES<-list(RESPONSE=c('ARG1','VEGFA','IL10','TIMP1'),SPP1=c('GPNMB','LGALS3','FABP5','CD9','LPL'),C1QC=c('C1QA','C1QB','C1QC'),IMMUNE=c('TGFB1','VSIG4','MARCO','CD274','IDO1'),MTA_PROD=c('MAT2A','MAT2B','AMD1','SRM','SMS'),MONO=c('FCN1','VCAN','S100A8','S100A9','CTSS','LILRB1'))
REL<-data.table(relation_id=paste0('A0',1:6),predictor=c('TYMP','RESPONSE','TYMP','MTA_PROD','RESPONSE','RESPONSE'),outcome=c('RESPONSE','SPP1','SPP1','RESPONSE','IMMUNE','C1QC'),expected=c('positive','positive','positive','positive','positive','two_sided'))
score_ctrl<-function(e,genes,ctrl_size,n_bins=25L){genes<-unique(intersect(toupper(genes),rownames(e)));av<-Matrix::rowMeans(e);pool<-rownames(e)[is.finite(av)];v<-av[is.finite(av)];n<-as.integer(round(length(v)/(n_bins-1L)));b<-floor(frank(v,ties.method='min')/n);names(b)<-pool;set.seed(SEED);ct<-character();for(x in unique(b[genes])){z<-names(b)[b==x & !names(b)%chin%genes];if(length(z)>ctrl_size)z<-sample(z,ctrl_size);ct<-union(ct,z)};as.numeric(Matrix::colMeans(e[genes,,drop=FALSE])-Matrix::colMeans(e[ct,,drop=FALSE]))}
resid<-function(v,s,g,state=NULL){if(is.null(state))return(as.numeric(lm.fit(cbind(1,s,g),v)$residuals));d<-data.frame(v=v,s=s,g=g,state=factor(state));as.numeric(residuals(lm(v~s+g+state,data=d)))}
rho<-function(x,y)if(uniqueN(x)<3||uniqueN(y)<3)NA_real_ else suppressWarnings(cor(x,y,method='spearman'))
sr<-function(v){v<-v[is.finite(v)];if(!length(v)||all(abs(v)<1e-15))return(1);suppressWarnings(wilcox.test(v,mu=0,exact=length(v)<50)$p.value)}
meta<-fread(metadata_path);obj<-readRDS(rds);eall<-LayerData(obj[['RNA']],layer='data')
if(dataset=='GSE207493'){
 meta<-meta[cell_type%chin%c('Myeloid','Dendritic')];meta[,`:=`(cell_id=cell_barcode,seurat_cell_name=paste0('GSE207493__GSE207493__',sub('^([^_]+)_','\\1__',cell_barcode)),patient_key=paste0('GSE207493::',sample),patient_id=sample,timepoint='Single',treatment='Unknown',response='Unknown',fine_state=paste(cell_type,leiden,sep='_'),state_class=fifelse(cell_type=='Dendritic','DC','Myeloid_unspecified'))]
}else{
 meta<-meta[phase7_broad_annotation=='Myeloid'];meta[,fine_state:=phase7_fine_annotation];meta[,state_class:=fifelse(grepl('^Macrophage',fine_state),'Macrophage',fifelse(grepl('^Monocyte',fine_state),'Monocyte',fifelse(grepl('DC|Dendritic',fine_state),'DC','Other')))]
}
if(length(setdiff(meta$seurat_cell_name,colnames(eall))))stop('ID mismatch')
e<-eall[,match(meta$seurat_cell_name,colnames(eall)),drop=FALSE];colnames(e)<-meta$cell_id;rm(obj,eall);gc()
u<-toupper(rownames(e));dup<-unique(u[duplicated(u)|duplicated(u,fromLast=TRUE)]);e<-e[!u%chin%dup,,drop=FALSE];rownames(e)<-u[!u%chin%dup]
cy<-fread(cycle_path);ss<-score_ctrl(e,cy[phase=='S',gene],43L);gg<-score_ctrl(e,cy[phase=='G2M',gene],54L)
need<-unique(c(unlist(GENES),'TYMP'));m<-as.matrix(e[intersect(need,rownames(e)),,drop=FALSE]);sdv<-apply(m,1,sd);z<-sweep(sweep(m,1,rowMeans(m),'-'),1,sdv,'/');z[!is.finite(z)]<-0
sc<-data.table(cell_id=colnames(e),TYMP=as.numeric(m['TYMP',]),S=ss,G=gg)
for(nm in names(GENES))sc[[nm]]<-colMeans(z[intersect(GENES[[nm]],rownames(z)),,drop=FALSE])
for(gene in GENES$RESPONSE)sc[[paste0('RESPONSE_NO_',gene)]]<-colMeans(z[setdiff(GENES$RESPONSE,gene),,drop=FALSE])
sc<-cbind(meta[,.(cell_id,patient_key,patient_id,timepoint,fine_state,state_class)],sc[,!"cell_id",with=FALSE])
unit_rows<-list();strata_rows<-list();loo_rows<-list();i<-0L;j<-0L;k<-0L
for(key_i in seq_len(nrow(unique(sc[,.(patient_key,timepoint)])))){
 key<-unique(sc[,.(patient_key,timepoint)])[key_i];d<-sc[patient_key==key$patient_key & timepoint==key$timepoint];if(nrow(d)<MIN_CELLS)next
 for(r in seq_len(nrow(REL))){rr<-REL[r];x<-d[[rr$predictor]];y<-d[[rr$outcome]];ok<-is.finite(x)&is.finite(y);if(sum(ok)<MIN_CELLS)next;xr<-resid(x[ok],d$S[ok],d$G[ok],d$fine_state[ok]);yr<-resid(y[ok],d$S[ok],d$G[ok],d$fine_state[ok]);i<-i+1L;unit_rows[[i]]<-data.table(dataset=dataset,relation_id=rr$relation_id,expected=rr$expected,patient_key=key$patient_key,timepoint=key$timepoint,n_cells=sum(ok),state_adjusted_rho=rho(xr,yr))}
 for(cls in unique(d$state_class)){dd<-d[state_class==cls];if(nrow(dd)<MIN_STRATUM)next;for(rid in c('A01','A02','A03','A05')){rr<-REL[relation_id==rid];x<-dd[[rr$predictor]];y<-dd[[rr$outcome]];ok<-is.finite(x)&is.finite(y);if(sum(ok)<MIN_STRATUM)next;j<-j+1L;strata_rows[[j]]<-data.table(dataset=dataset,state_class=cls,relation_id=rid,expected=rr$expected,patient_key=key$patient_key,timepoint=key$timepoint,n_cells=sum(ok),rho=rho(resid(x[ok],dd$S[ok],dd$G[ok]),resid(y[ok],dd$S[ok],dd$G[ok])))} }
 for(gene in GENES$RESPONSE){resp<-paste0('RESPONSE_NO_',gene);for(edge in c('TYMP_TO_RESPONSE','RESPONSE_TO_SPP1','RESPONSE_TO_IMMUNE')){x<-if(edge=='TYMP_TO_RESPONSE')d$TYMP else d[[resp]];y<-if(edge=='TYMP_TO_RESPONSE')d[[resp]] else if(edge=='RESPONSE_TO_SPP1')d$SPP1 else d$IMMUNE;ok<-is.finite(x)&is.finite(y);if(sum(ok)<MIN_CELLS)next;k<-k+1L;loo_rows[[k]]<-data.table(dataset=dataset,left_out_gene=gene,edge=edge,patient_key=key$patient_key,timepoint=key$timepoint,n_cells=sum(ok),state_adjusted_rho=rho(resid(x[ok],d$S[ok],d$G[ok],d$fine_state[ok]),resid(y[ok],d$S[ok],d$G[ok],d$fine_state[ok])))} }
}
summ<-function(dt,by){p<-dt[,.(patient_rho=median(get(if('state_adjusted_rho'%chin%names(dt))'state_adjusted_rho' else 'rho'),na.rm=TRUE)),by=c(by,'patient_key')];s<-p[,.(n_patients=.N,median_rho=median(patient_rho),positive_fraction=mean(patient_rho>0),p_value=sr(patient_rho)),by=by];s[,q_value:=p.adjust(p_value,'BH'),by=dataset];s}
unit<-rbindlist(unit_rows);strata<-rbindlist(strata_rows);loo<-rbindlist(loo_rows)
us<-summ(unit,c('dataset','relation_id','expected'));us[,decision:=fifelse(expected=='two_sided',fifelse(q_value<0.05,'ASSOCIATED','NOT_SUPPORTED'),fifelse(q_value<0.05 & median_rho>0 & positive_fraction>=2/3,'SUPPORTED','NOT_SUPPORTED'))]
ssum<-summ(strata,c('dataset','state_class','relation_id','expected'));ssum[,decision:=fifelse(q_value<0.05 & median_rho>0 & positive_fraction>=2/3,'SUPPORTED','NOT_SUPPORTED')]
lsum<-summ(loo,c('dataset','left_out_gene','edge'));lsum[,decision:=fifelse(q_value<0.05 & median_rho>0 & positive_fraction>=2/3,'SUPPORTED','NOT_SUPPORTED')]
fwrite(unit,file.path(part,paste0(dataset,'_state_adjusted_units_v1.tsv')),sep='\t');fwrite(us,file.path(part,paste0(dataset,'_state_adjusted_summary_v1.tsv')),sep='\t');fwrite(strata,file.path(part,paste0(dataset,'_state_stratified_units_v1.tsv.gz')),sep='\t');fwrite(ssum,file.path(part,paste0(dataset,'_state_stratified_summary_v1.tsv')),sep='\t');fwrite(loo,file.path(part,paste0(dataset,'_response_loo_units_v1.tsv.gz')),sep='\t');fwrite(lsum,file.path(part,paste0(dataset,'_response_loo_summary_v1.tsv')),sep='\t');fwrite(data.table(dataset=dataset,selected_cells=nrow(meta),patients=uniqueN(meta$patient_key),fine_states=uniqueN(meta$fine_state),status='PASS_STATE_CONFOUNDING_AUDIT_V1'),file.path(part,paste0(dataset,'_state_audit_status_v1.tsv')),sep='\t');file.rename(part,out);cat(sprintf('{"status":"PASS_STATE_CONFOUNDING_AUDIT_V1","dataset":"%s"}\n',dataset))
