p <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1/07_trajectory/monocle3_v3/primary/GSE156632/MONOCLE3_CDS.rds'
x <- readRDS(p)
cat('SLOTS',names(attributes(x)),'\n')
ic <- attr(x,'int_colData')
cat('INTCOL',names(attributes(ic)),'\n')
ld <- attr(ic,'listData')
cat('LISTDATA',names(ld),'\n')
rd <- ld$reducedDims
cat('REDUCED',names(attributes(rd)),names(attr(rd,'listData')),'\n')
pg <- attr(x,'principal_graph')$UMAP
cat('GRAPH_TYPE',class(pg),'LEN',length(pg),'\n')
for(i in seq_along(pg))cat(i,typeof(pg[[i]]),length(pg[[i]]),'\n')
cat('NODE_MATRIX',dim(attr(x,'principal_graph_aux')$UMAP$dp_mst),'\n')
