suppressPackageStartupMessages({library(monocle3);library(igraph);library(SingleCellExperiment)})
root <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1'
out <- file.path(root,'phase192_sections28_212_figures_20260917_v1/inputs/trajectory')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for (mode in c('primary','without_stress_cycle')) for(study in c('GSE156632','GSE159115','GSE178481','GSE207493','GSE242299')) {
 path <- file.path(root,'part2_rebuild_v1/07_trajectory/monocle3_v3',mode,study,'MONOCLE3_CDS.rds')
 cds <- readRDS(path)
 xy <- reducedDims(cds)$UMAP
 nodes <- t(cds@principal_graph_aux[['UMAP']]$dp_mst)
 ed <- igraph::as_data_frame(principal_graph(cds)[['UMAP']],what='edges')
 od <- file.path(out,mode,study);dir.create(od,recursive=TRUE,showWarnings=FALSE)
 write.table(data.frame(cell_id=rownames(xy),x=xy[,1],y=xy[,2]),file.path(od,'cells.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(node=rownames(nodes),x=nodes[,1],y=nodes[,2]),file.path(od,'nodes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(ed,file.path(od,'edges.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 cat(mode,study,nrow(xy),nrow(nodes),nrow(ed),'\n')
}
writeLines(capture.output(sessionInfo()),file.path(out,'R_SESSION_INFO.txt'))
