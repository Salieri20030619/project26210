suppressPackageStartupMessages({library(monocle3);library(Matrix);library(data.table);library(jsonlite);library(SingleCellExperiment)})
R <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1'
I <- file.path(R,'07_trajectory/matrices_v1'); O <- file.path(R,'07_trajectory/monocle3_v3')
stopifnot(file.exists(file.path(I,'STATUS.json')),!dir.exists(O));dir.create(O)
set.seed(262148);options(mc.cores=1)
scale_partition <- function(pt,part) {
 out <- rep(NA_real_,length(pt))
 for (pa in unique(part)) { k<-which(part==pa & is.finite(pt));if(length(k)>1) { rr<-range(pt[k]);out[k]<-if(diff(rr)>0)(pt[k]-rr[1])/diff(rr) else 0 } }
 out
}
qc<-list(); roots_out<-list(); outall<-list()
for(mode in c('primary','without_stress_cycle')) for(study in sort(list.dirs(file.path(I,mode),recursive=FALSE,full.names=FALSE))) {
 cat('START',mode,study,format(Sys.time()),'\n');flush.console()
 p<-file.path(I,mode,study);od<-file.path(O,mode,study);dir.create(od,recursive=TRUE)
 md<-fread(file.path(p,'metadata.tsv.gz'));latent<-fread(file.path(p,'pca.tsv.gz'),header=TRUE);genes<-fread(file.path(p,'genes.tsv'),header=FALSE)[[1]]
 x<-as(readMM(gzfile(file.path(p,'counts.mtx.gz'))),'CsparseMatrix');cells<-md$cell_id
 stopifnot(ncol(x)==length(cells),nrow(x)==length(genes),identical(latent[[1]],cells));colnames(x)<-cells;rownames(x)<-genes
 cm<-as.data.frame(md);rownames(cm)<-cells;gm<-data.frame(gene_short_name=genes,row.names=genes)
 cds<-new_cell_data_set(x,cell_metadata=cm,gene_metadata=gm);pcs<-as.matrix(latent[,-1]);rownames(pcs)<-cells;reducedDims(cds)$PCA<-pcs
 cds<-reduce_dimension(cds,reduction_method='UMAP',preprocess_method='PCA',umap.metric='cosine',umap.min_dist=.1,umap.n_neighbors=30L,umap.fast_sgd=FALSE,cores=1)
 cds<-cluster_cells(cds,reduction_method='UMAP',k=30L,cluster_method='leiden',random_seed=262148L)
 cds<-learn_graph(cds,use_partition=TRUE,close_loop=FALSE,learn_graph_control=list(minimal_branch_len=10L),verbose=FALSE)
 part<-as.character(partitions(cds)[cells]);rootsets<-list(primary=character(),root_q25=character(),root_q50=character())
 for(pa in unique(part)) {
  kk<-which(part==pa & md$root_eligible==TRUE)
  if(length(kk)<20)next
  center<-apply(pcs[kk,,drop=FALSE],2,median);distance<-rowSums(sweep(pcs[kk,,drop=FALSE],2,center,'-')^2);ordered<-kk[order(distance)]
  for(v in names(rootsets)) {
   frac<-switch(v,primary=0,root_q25=.25,root_q50=.5);ri<-ordered[1+floor(frac*(length(ordered)-1))];rootsets[[v]]<-c(rootsets[[v]],cells[ri]);roots_out[[length(roots_out)+1]]<-data.table(mode=mode,study=study,partition=pa,root_version=v,root_cell=cells[ri],root_candidates=length(kk))
  }
 }
 out<-copy(md);out[,`:=`(mode=mode,partition=part)]
 for(v in names(rootsets)) {
  pt<-rep(NA_real_,length(cells))
  if(length(rootsets[[v]])>0) { temp<-order_cells(cds,root_cells=rootsets[[v]]);pt<-pseudotime(temp)[cells] }
  out[,(paste0('pt_',v)):=scale_partition(pt,part)]
 }
 if(length(rootsets$primary)>0)cds<-order_cells(cds,root_cells=rootsets$primary)
 fwrite(out,file.path(od,'CELL_PSEUDOTIME.tsv.gz'),sep='\t');saveRDS(cds,file.path(od,'MONOCLE3_CDS.rds'),compress='gzip')
 qt<-out[,.(cells=.N,patients=uniqueN(curated_patient_id),root_eligible=sum(root_eligible),C4_cells=sum(candidate=='4'),finite_primary=sum(is.finite(pt_primary)),C4_finite=sum(candidate=='4' & is.finite(pt_primary))),by=.(mode,study,partition)]
 qc[[length(qc)+1]]<-qt;outall[[length(outall)+1]]<-out;cat('DONE',mode,study,'partitions',uniqueN(part),'finite',sum(is.finite(out$pt_primary)),'\n');flush.console();rm(cds,temp,x);gc()
}
fwrite(rbindlist(qc),file.path(O,'PARTITION_QC.tsv'),sep='\t');fwrite(rbindlist(roots_out),file.path(O,'ROOTS.tsv'),sep='\t');fwrite(rbindlist(outall),file.path(O,'ALL_CELL_PSEUDOTIME.tsv.gz'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(O,'R_SESSION_INFO.txt'))
write_json(list(status='COMPLETE_PARTITIONED_TRAJECTORIES',figures_created=FALSE,models=length(outall),use_partition=TRUE,rootless_partitions_retained=TRUE),file.path(O,'STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
cat('ALL_MONOCLE3_DONE\n')
