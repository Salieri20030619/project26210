library(data.table)
p<-'part2_rebuild_v1/07_trajectory/matrices_v1/primary/GSE156632'
a<-fread(file.path(p,'metadata.tsv.gz'));b<-fread(file.path(p,'pca.tsv.gz'))
print(names(a));print(names(b));print(head(a$cell_id));print(head(b[[1]]));print(class(a$cell_id));print(class(b[[1]]));print(dim(a));print(dim(b))
