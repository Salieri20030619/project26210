x <- readRDS('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1/07_trajectory/monocle3_v3/primary/GSE156632/MONOCLE3_CDS.rds')
g<-attr(x,'principal_graph')$UMAP
for(i in seq_along(g[[9]]))cat('ATTR',i,names(g[[9]][[i]]),'\n')
print(head(g[[9]][[3]]$name))
print(head(attr(x,'principal_graph_aux')$UMAP$pr_graph_cell_proj_closest_vertex))
