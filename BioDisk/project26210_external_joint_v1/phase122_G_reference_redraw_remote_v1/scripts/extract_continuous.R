args<-commandArgs(TRUE)
base<-'/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1'
pts<-c('GSE207493__RCC104','GSE207493__RCC86','GSE156632__RCC4','GSE156632__RCC5')
for(pt in pts){
 f<-file.path(base,'03_malignancy_cnv/infercnv_v2/patients',pt,'infercnv_run/preliminary.infercnv_obj')
 o<-readRDS(f);m<-o@expr.data
 stopifnot(all(is.finite(m)),all(m>0))
 v<-colMeans(abs(log2(m)))
 d<-data.frame(exact_cell_id=colnames(m),patient=pt,cnv_score=as.numeric(v))
 write.table(d,file.path(args[1],paste0(pt,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 cat(pt,nrow(m),ncol(m),head(colnames(m),1),'\n')
 rm(o,m,d,v);gc()
}
