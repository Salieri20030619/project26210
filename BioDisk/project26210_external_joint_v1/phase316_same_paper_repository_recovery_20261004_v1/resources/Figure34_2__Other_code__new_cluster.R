# new cluster
sc_new<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/sc_CRC_merge_cluster_major.rds")
subset(sc_new,cells = )
data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
sc_B_T_Plasma<-sc_new[,sc_new$new_Major %in% c("TNKILC","B","Plasma")]
rm(sc_new)
table(sc_B_T_Plasma$new_Major)

new_sc_T<-sc_B_T_Plasma[,sc_B_T_Plasma$new_Major == "TNKILC"]
new_sc_T@assays

library(dplyr)
new_sc_T <- RunUMAP(new_sc_T, reduction = "harmony", dims = 1:30)
new_sc_T <- FindNeighbors(new_sc_T, reduction = "harmony", dims = 1:30)
new_sc_T <- FindClusters(new_sc_T,resolution = 0.1)
##作图
#group_by_cluster
library(RColorBrewer)
p1<-DimPlot(new_sc_T, reduction = "umap", label=T,cols = brewer.pal(10,"Paired"),group.by = "seurat_clusters") 
p1
p2<-DimPlot(new_sc_T[,!is.na(new_sc_T$new_ident)], reduction = "umap", group.by = "new_ident",na.value = alpha("white",0),shuffle = T,cols=brewer.pal(12,"Paired")) 
p2
p1+p2

new_sc_T@meta.data$class_2<-"Remain"
new_sc_T@meta.data$class_2[new_sc_T@meta.data$seurat_clusters == 0]<-"CD4 T"
new_sc_T@meta.data$class_2[new_sc_T@meta.data$seurat_clusters == 1]<-"CD8 T"
new_sc_T@meta.data$class_2[new_sc_T@meta.data$seurat_clusters == 2]<-"Mast"
table(new_sc_T@meta.data$class_2)
FeaturePlot(new_sc_T,features = "NKG7")

new_sc_T


unique(ggplot_build(  p2  )$data[[1]][,1])

DefaultAssay(sc_DC) <- "integrated"
#sc_DC<-sc_ref_rename[,sc_ref_rename@meta.data$cluster_last_version %in% c("DC","pDC","Macrophage","Monocyte")]


sc_DC <- ScaleData(sc_DC, verbose = FALSE)
sc_DC <- RunPCA(sc_DC, npcs = 50, verbose = FALSE)
ElbowPlot(sc_DC)
#VizDimLoadings(sc_DC, dims = 25:30, reduction = "pca")
DefaultAssay(sc_DC) <- "integrated"
sc_DC <- RunUMAP(sc_DC, reduction = "pca", dims = 1:20)
sc_DC <- FindNeighbors(sc_DC, reduction = "pca", dims = 1:20)
sc_DC <- FindClusters(sc_DC, resolution = 0.5)
#sc_DC$seurat_clusters %>% table()
#FeaturePlot(sc_DC,features = c("MKI67","MCM5","PCNA"))
DimPlot(sc_DC,group.by = "seurat_clusters",cols = c(brewer.pal(10,"Paired")))


