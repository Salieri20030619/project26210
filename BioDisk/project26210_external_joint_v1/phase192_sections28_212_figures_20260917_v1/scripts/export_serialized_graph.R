root <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1'
out <- file.path(root,'phase192_sections28_212_figures_20260917_v1/inputs/trajectory')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(mode in c('primary','without_stress_cycle'))for(study in c('GSE156632','GSE159115','GSE178481','GSE207493','GSE242299')){
 x<-readRDS(file.path(root,'part2_rebuild_v1/07_trajectory/monocle3_v3',mode,study,'MONOCLE3_CDS.rds'))
 rd<-attr(attr(x,'int_colData'),'listData')$reducedDims
 xy<-attr(rd,'listData')$UMAP
 aux<-attr(x,'principal_graph_aux')$UMAP
 nodes<-t(aux$dp_mst);g<-attr(x,'principal_graph')$UMAP;vn<-g[[9]][[3]]$name
 stopifnot(length(vn)==nrow(nodes),setequal(vn,rownames(nodes)),length(g[[3]])==length(g[[4]]),all(g[[3]]>=0),all(g[[4]]>=0),all(g[[3]]<length(vn)),all(g[[4]]<length(vn)))
 edges<-data.frame(from=vn[g[[3]]+1],to=vn[g[[4]]+1])
 closest<-aux$pr_graph_cell_proj_closest_vertex
 od<-file.path(out,mode,study);dir.create(od,recursive=TRUE,showWarnings=FALSE)
 write.table(data.frame(cell_id=rownames(xy),x=xy[,1],y=xy[,2]),file.path(od,'cells.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(node=rownames(nodes),x=nodes[,1],y=nodes[,2]),file.path(od,'nodes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(edges,file.path(od,'edges.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(cell_id=rownames(closest),node=paste0('Y_',closest[,1])),file.path(od,'closest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 cat(mode,study,nrow(xy),nrow(nodes),nrow(edges),'\n')
}
