library(harmony)
library(Seurat)
library(tidyverse)
library(patchwork)
library(pacman)
library(RColorBrewer)
p_unload(Seurat)
p_load(Seurat)

source('/Path/scripts/Dotplot.R')

########## integrated using Seurat ------
CRC_obj<-readRDS("/Path/data/CRC_sc_list.RDS")
all.genes<-readRDS("/Path/data/all.genes.RDS")
features <- SelectIntegrationFeatures(object.list = CRC_obj,nfeatures = 2000)
CRC_ICB.anchors <- FindIntegrationAnchors(object.list = CRC_obj, anchor.features = features)
CRC_ICB.combined <- IntegrateData(anchorset = CRC_ICB.anchors，normalization.method = "SCT", verbose = TRUE, features.to.integrate = all_genes)
saveRDS(CRC_ICB.combined,file='/Path/data/integrate_v1.RDS')

########## integrated using Seurat ------
integrated_sc <- readRDS('/Path/data/integrate_v1.RDS')

integrated_sc <- ScaleData(integrated_sc, verbose = FALSE)
integrated_sc <- RunPCA(integrated_sc, npcs = 50, verbose = FALSE)
ElbowPlot(integrated_sc)
DefaultAssay(integrated_sc) <- "integrated"
integrated_sc <- RunUMAP(integrated_sc, reduction = "pca", dims = 1:20)
integrated_sc <- FindNeighbors(integrated_sc, reduction = "pca", dims = 1:20)
integrated_sc <- FindClusters(integrated_sc, resolution = 1)
library(RColorBrewer)

integrated_sc@meta.data$first_around_anno<-"Remain"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 0] <- "CD8_gd_NK"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 1] <- "CD4"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 2] <- "B"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 3] <- "Treg"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 4] <- "Plasma"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 5] <- "Plasma"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 6] <- "Mast"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 7] <- "Epithelial c1"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 8] <- "APOE-Fibroblast"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 9] <- "Mac-like"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 10] <- "Th1-like"

integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 11] <- "Epithelial c2"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 12] <- "Epithelial c3"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 13] <- "Myofib"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 14] <- "Endo"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 15] <- "Fibroblast c3"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 16] <- "OGN Fib"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 17] <- "Mac+ cDC2 Neu"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 18] <- "Goblet-like"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 19] <- "Prolif T"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 20] <- "Fib (pericyte-like)"

integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 21] <- "pericyte"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 22] <- "Endo c2"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 23] <- "CD45-CD11b-CD14+ cell + prolif Myeloid"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 24] <- "cDC2"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 25] <- "ILC"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 26] <- "Glial"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 27] <- "Endo_3"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 28] <- "Fib c4"

DimPlot(integrated_sc,group.by = 'first_around_anno',label = T,repel = T)



### myeloid ------------
integrated_sc_myeloid <- subset(integrated_sc,subset = first_around_anno %in%
                                  c("CD45-CD11b-CD14+ cell + prolif Myeloid","Mac-like","cDC2","Mac+ cDC2 Neu"))
DimPlot(integrated_sc_myeloid,group.by = 'first_around_anno',label = T,repel = T)
VlnPlot(integrated_sc_myeloid,features = c('LYZ','C1QA','S100A8','CD1C'),
        pt.size = 0,group.by = 'first_around_anno')

DefaultAssay(integrated_sc_myeloid) <- 'RNA'
FeaturePlot(integrated_sc_myeloid,features = c('PTPRC','MKI67','C1QA','HLA-A','ITGAM','CD14','CLEC9A','ITGAX'))
DimPlot(integrated_sc_myeloid,group.by = 'first_around_anno',label = T,repel = T)+NoLegend()

#### recluster

DefaultAssay(integrated_sc_myeloid) <- 'integrated'
integrated_sc_myeloid <- FindNeighbors(integrated_sc_myeloid, dims = 1:30)
integrated_sc_myeloid <- FindClusters(integrated_sc_myeloid, resolution = 0.3)
integrated_sc_myeloid <- RunUMAP(integrated_sc_myeloid, dims = 1:30)

table(integrated_sc_myeloid$integrated_snn_res.0.3)

DimPlot(integrated_sc_myeloid,group.by = 'integrated_snn_res.0.3')

DefaultAssay(integrated_sc_myeloid) <- 'RNA'
VlnPlot(integrated_sc_myeloid,group.by = 'integrated_snn_res.0.3',pt.size = 0,
        features = c('PTPRC','MKI67','LAMP3','S100A8','ITGAM','CCR7','CLEC9A',
                     'ITGAX','CD1C','SPP1','CSF3R','LILRB2'))



VlnPlot(integrated_sc_myeloid,group.by = 'integrated_snn_res.0.3',pt.size = 0,
        features = c('CD163','MRC1','SPP1','ITGAM','HLA-A','CD14','PDPN',
                     'ARG1','VEGFA','IL1B','C1QA','CCR2','IGF1','ITGAX','LYVE1'))

FeaturePlot(integrated_sc_myeloid,features = c('CSF3R','CXCR2','FCGR3A'))
FeaturePlot(integrated_sc_myeloid,features = c('CD163','MRC1','SPP1','PDPN','IGF1','C1QA','IL1B'))
FeaturePlot(integrated_sc_myeloid,features = c('MKI67'))
FeaturePlot(integrated_sc_myeloid,features = c('nCount_RNA','nFeature_RNA'),max.cutoff = 1000)
FeaturePlot(integrated_sc_myeloid,features = c('PTPRC','C1QA'))


Idents(integrated_sc_myeloid) <- integrated_sc_myeloid$integrated_snn_res.0.3

markers_myeloid <- FindAllMarkers(integrated_sc_myeloid,assay = 'RNA',
                                  only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

integrated_sc_myeloid$myeloid_Minor <- integrated_sc_myeloid$integrated_snn_res.0.3

levels(integrated_sc_myeloid$myeloid_Minor)[1] <- 'FOLR2+ Mac' # (M2-like) FOLR2
levels(integrated_sc_myeloid$myeloid_Minor)[2] <- 'Monocyte'
levels(integrated_sc_myeloid$myeloid_Minor)[3] <- 'cDC2'       # CD1C FCER1A
levels(integrated_sc_myeloid$myeloid_Minor)[4] <- 'SPP1+ Mac'  # SPP1
levels(integrated_sc_myeloid$myeloid_Minor)[6] <- 'cDC1&cDC3'  # CLEC9A, LAMP3
levels(integrated_sc_myeloid$myeloid_Minor)[7] <- 'Neutrophil' # CSF3R

levels(integrated_sc_myeloid$myeloid_Minor)[5] <- 'CD45- C1QA+'               #--> filter
levels(integrated_sc_myeloid$myeloid_Minor)[8] <- 'proliferation myeloid'     #--> filter

levels(integrated_sc_myeloid$myeloid_Minor)

DimPlot(integrated_sc_myeloid,group.by = 'myeloid_Minor')
integrated_sc_myeloid_fil <- subset(integrated_sc_myeloid,subset=myeloid_Minor%in%
                                      levels(integrated_sc_myeloid$myeloid_Minor)[c(1:4,6,7)])
DimPlot(integrated_sc_myeloid_fil,group.by = 'myeloid_Minor')

DefaultAssay(integrated_sc_myeloid_fil) <- 'RNA'
Dotplot(DotPlot(integrated_sc_myeloid_fil,group.by = 'myeloid_Minor',
                features = c('FOLR2','CD163L1','FCN1','THBS1',
                             'CD1C','FCER1A', 'CLEC10A','SPP1','FN1',
                             'CLEC9A','XCR1','LAMP3','CCR7',
                             'CSF3R','S100A8'))$data)

FeaturePlot(integrated_sc_myeloid_fil,
            features = c('FOLR2','FCN1','THBS1',
                         'CD1C','FCER1A','SPP1','FN1',
                         'CSF3R','S100A8'))

Dotplot(DotPlot(integrated_sc_myeloid_fil,group.by = 'myeloid_Minor',
                features = c('ITGAM','C1QA','C1QC','CD163','FCGR3A','APOE'))$data)


FeaturePlot(integrated_sc_myeloid,features = c('FCGR3A','CD163','APOE','FCGR3B'))
FeaturePlot(integrated_sc_myeloid,features = c('FCN1', 'CD14', 'CD36', 'SELL', 'S100A8', 'S100A12',
                                               'CLEC12A', 'MS4A6A', 'CXCL14'))
saveRDS(integrated_sc_myeloid_fil,'integrated_sc_myeloid_fil_temp.rds')


DimPlot(integrated_sc_myeloid_fil,group.by = 'myeloid_Minor')+
  FeaturePlot(integrated_sc_myeloid_fil,
              features = c('FOLR2','FCN1','THBS1',
                           'CD1C', 'CLEC10A','SPP1','FN1',
                           'CLEC9A','LAMP3','CCR7',
                           'CSF3R','S100A8'))

ggsave('integrated_sc_myeloid_fil_umap_feature.pdf',height = 8,width = 20)



## subset cDC1&cDC3
sc_cDC13 <- subset(integrated_sc_myeloid_fil,subset=myeloid_Minor=='cDC1&cDC3')
FeaturePlot(sc_cDC13,features = c('CCR7','CLEC9A','LAMP3'))


DefaultAssay(sc_cDC13) <- 'integrated'
sc_cDC13 <- FindNeighbors(sc_cDC13, dims = 1:30)
sc_cDC13 <- FindClusters(sc_cDC13, resolution = 0.3)
sc_cDC13 <- RunUMAP(sc_cDC13, dims = 1:30)

DimPlot(sc_cDC13,group.by = 'integrated_snn_res.0.3')

DefaultAssay(sc_cDC13) <- 'RNA'
FeaturePlot(sc_cDC13,features = c('CLEC9A','CADM1','CCR7','FSCN1','LAMP3','GZMB','IL3RA','CLEC4C','LILRB4'))
FeaturePlot(sc_cDC13,features = c('CADM1','BLTA','ITGAX'))

markers_cDC13 <- FindAllMarkers(sc_cDC13,assay = 'RNA',only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

sc_cDC13$DC_Minor <- sc_cDC13$integrated_snn_res.0.3
levels(sc_cDC13$DC_Minor)[1] <- 'cDC1'  # CLEC9A CADM1 XCR1
levels(sc_cDC13$DC_Minor)[2] <- 'cDC3'  # CCR7 LAMP3  FSCN1
levels(sc_cDC13$DC_Minor)[3] <- 'pDC'  # IL3RA LILRB4 CLEC4C GZMB

DimPlot(sc_cDC13,group.by = 'DC_Minor',label = T,repel = T)

Dotplot(DotPlot(sc_cDC13,group.by = 'DC_Minor',
                features =
                  c('CLEC9A','CADM1','XCR1','CCR7','FSCN1','LAMP3','GZMB','IL3RA','CLEC4C','LILRB4'))$data)
saveRDS(sc_cDC13,'sc_cDC.rds')

### combine all myeloid idents:
common_cells <- intersect(colnames(sc_cDC13), colnames(integrated_sc_myeloid_fil))
integrated_sc_myeloid_fil$myeloid_Minor <- as.character(integrated_sc_myeloid_fil$myeloid_Minor)
integrated_sc_myeloid_fil$myeloid_Minor[common_cells] <- as.character(sc_cDC13$DC_Minor)
table(integrated_sc_myeloid_fil$myeloid_Minor)
integrated_sc_myeloid_fil$myeloid_Minor <- factor(integrated_sc_myeloid_fil$myeloid_Minor,
                                                  levels = c('cDC1','cDC2','cDC3','pDC',
                                                             'FOLR2+ Mac','SPP1+ Mac','Monocyte','Neutrophil'))
table(integrated_sc_myeloid_fil$myeloid_Minor)
levels(integrated_sc_myeloid_fil$myeloid_Minor)
DimPlot(integrated_sc_myeloid_fil,group.by = 'myeloid_Minor',cols = brewer.pal(n=8,name = 'Paired'))+
  labs(title = NULL)
ggsave('UMAP_sc_myeloid_fil.pdf',height = 6,width = 7)
saveRDS(integrated_sc_myeloid_fil,'integrated_sc_myeloid_fil.rds')




common_cells <- intersect(colnames(sc_cDC13), colnames(integrated_sc_myeloid))
integrated_sc_myeloid$myeloid_Minor <- as.character(integrated_sc_myeloid$myeloid_Minor)
integrated_sc_myeloid$myeloid_Minor[common_cells] <- as.character(sc_cDC13$DC_Minor)
table(integrated_sc_myeloid$myeloid_Minor)
integrated_sc_myeloid$myeloid_Minor <- factor(integrated_sc_myeloid$myeloid_Minor,
                                              levels = c('cDC1','cDC2','cDC3','pDC',
                                                         'FOLR2+ Mac','SPP1+ Mac','Monocyte','Neutrophil',
                                                         'proliferation myeloid','CD45- C1QA+'))
table(integrated_sc_myeloid$myeloid_Minor)
levels(integrated_sc_myeloid$myeloid_Minor)
DimPlot(integrated_sc_myeloid,group.by = 'myeloid_Minor',cols = brewer.pal(n=10,name = 'Paired'))+
  labs(title = NULL)
ggsave('UMAP_sc_myeloid.pdf',height = 6,width = 7)
saveRDS(integrated_sc_myeloid,'integrated_sc_myeloid.rds')



Dotplot(DotPlot(integrated_sc_myeloid_fil,group.by = 'myeloid_Minor',
                features =
                  c('CLEC9A','CADM1','XCR1','CD1C','FCER1A',
                    'CCR7','FSCN1','LAMP3','GZMB','IL3RA','CLEC4C',
                    'FOLR2','CD163L1','SPP1','FN1','FCN1','THBS1','CSF3R','S100A8'))$data)
ggsave('Dotplot_sc_marker_myeloid_fil.pdf',height = 3.5,width = 7.5)






### Endothelial ------------
integrated_sc_endo <- subset(integrated_sc,subset = first_around_anno %in%
                               c("Endo",'Endo c2','Endo_3'))
DimPlot(integrated_sc_endo,group.by = 'first_around_anno',label = T,repel = T)
VlnPlot(integrated_sc_endo,features = c('PECAM1','ENG','PROX1','CCL21','HEY1','IGFBP3','CD36','CA4','ACKR1'),
        pt.size = 0,group.by = 'first_around_anno')

FeaturePlot(integrated_sc_endo,features = c('PECAM1','ENG','PROX1','CCL21','HEY1','IGFBP3','CD36','CA4','ACKR1','Cd31'))

FeaturePlot(integrated_sc_endo,features = c('CDH5','NR2F2'))

integrated_sc_endo$Endo_major <- integrated_sc_endo$first_around_anno

integrated_sc_endo$Endo_major <- as.factor(integrated_sc_endo$Endo_major)
levels(integrated_sc_endo$Endo_major)
levels(integrated_sc_endo$Endo_major)[1] <- 'VECs'
levels(integrated_sc_endo$Endo_major)[2] <- 'VECs'
levels(integrated_sc_endo$Endo_major)[2] <- 'LECs'
levels(integrated_sc_endo$Endo_major)

DimPlot(integrated_sc_endo,group.by = 'Endo_major',label = T,repel = T)
saveRDS(integrated_sc_endo,'integrated_sc_endo.rds')





### Fibroblast ---------

integrated_sc@meta.data$second_around_anno[integrated_sc@meta.data$first_around_anno %in% c("APOE-Fibroblast","OGN Fib")]<-"Fibroblast"
integrated_sc@meta.data$second_around_anno[integrated_sc@meta.data$first_around_anno %in% c("Fib (pericyte-like)","Fib c4")]<-"Smooth Muscle"
integrated_sc@meta.data$second_around_anno[integrated_sc@meta.data$first_around_anno == "Fibroblast c3"]<-"POSTN+ Myofib"
integrated_sc@meta.data$second_around_anno[integrated_sc@meta.data$first_around_anno == "Myofib"]<-"FAP+ Myofib"
integrated_sc@meta.data$second_around_anno[integrated_sc@meta.data$first_around_anno == "pericyte"]<-"Pericyte"

### Epithelial -------
integrated_sc@meta.data$second_around_anno[integrated_sc@meta.data$first_around_anno %in%
                                             c("Epithelial c1","Epithelial c2",
                                               'Epithelial c3','Goblet-like')]<-"Epithelium"


### load T recluster results from Zurui (second_around_anno_new) -------

integrated_sc_T <- readRDS('/Path/data/Tcell.RDS')



#### Combine UCASpatial_clus_v1 --------

integrated_sc$UCASpatial_clus_v1 <- integrated_sc$first_around_anno

### T
common_cells <- intersect(colnames(integrated_sc), colnames(integrated_sc_T))
integrated_sc$UCASpatial_clus_v1 <- as.character(integrated_sc$UCASpatial_clus_v1)
integrated_sc$UCASpatial_clus_v1[common_cells] <- as.character(integrated_sc_T$second_around_anno_new)
table(integrated_sc$UCASpatial_clus_v1)

### myeloid
common_cells <- intersect(colnames(integrated_sc), colnames(integrated_sc_myeloid))
integrated_sc$UCASpatial_clus_v1 <- as.character(integrated_sc$UCASpatial_clus_v1)
integrated_sc$UCASpatial_clus_v1[common_cells] <- as.character(integrated_sc_myeloid$myeloid_Minor)
table(integrated_sc$UCASpatial_clus_v1)

### endo
common_cells <- intersect(colnames(integrated_sc), colnames(integrated_sc_endo))
integrated_sc$UCASpatial_clus_v1 <- as.character(integrated_sc$UCASpatial_clus_v1)
integrated_sc$UCASpatial_clus_v1[common_cells] <- as.character(integrated_sc_endo$Endo_major)
table(integrated_sc$UCASpatial_clus_v1)

### Epi
integrated_sc@meta.data$UCASpatial_clus_v1[integrated_sc@meta.data$first_around_anno %in%
                                             c("Epithelial c1","Epithelial c2",
                                               'Epithelial c3','Goblet-like')]<-"Epithelium"
table(integrated_sc$UCASpatial_clus_v1)

### Fibro/Myofibro
integrated_sc@meta.data$UCASpatial_clus_v1[integrated_sc@meta.data$first_around_anno %in% c("APOE-Fibroblast","OGN Fib")]<-"Fibroblast"
integrated_sc@meta.data$UCASpatial_clus_v1[integrated_sc@meta.data$first_around_anno %in% c("Fib (pericyte-like)","Fib c4")]<-"Smooth Muscle"
integrated_sc@meta.data$UCASpatial_clus_v1[integrated_sc@meta.data$first_around_anno == "Fibroblast c3"]<-"POSTN+ Myofib"
integrated_sc@meta.data$UCASpatial_clus_v1[integrated_sc@meta.data$first_around_anno == "Myofib"]<-"FAP+ Myofib"
integrated_sc@meta.data$UCASpatial_clus_v1[integrated_sc@meta.data$first_around_anno == "pericyte"]<-"Pericyte"
table(integrated_sc$UCASpatial_clus_v1)

DimPlot(integrated_sc,group.by = 'UCASpatial_clus_v1')

### filter proliferation T , proliferation myeloid and CD45- C1QA+
integrated_sc_fil <- subset(integrated_sc,subset=UCASpatial_clus_v1%in%
                              names(table(integrated_sc$UCASpatial_clus_v1))[c(1:4,6:26,29:32)])
table(integrated_sc_fil$UCASpatial_clus_v1)
DimPlot(integrated_sc_fil,group.by = 'UCASpatial_clus_v1',label = T,repel = T)

names(table(integrated_sc_fil$UCASpatial_clus_v1))
integrated_sc_fil$UCASpatial_clus_v1 <- factor(integrated_sc_fil$UCASpatial_clus_v1,
                                               levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                          "CD4 Naive","CD4 Tem","CD4 Treg","Th1/Th17",
                                                          "NK","ILC",
                                                          "Mast",
                                                          "Monocyte","Neutrophil","SPP1+ Mac","FOLR2+ Mac",
                                                          "cDC1","cDC2","cDC3","pDC",
                                                          "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                          "Smooth Muscle","Pericyte",
                                                          "VECs","LECs","Glial",
                                                          "B","Plasma","Epithelium"))
DimPlot(integrated_sc_fil,group.by = 'UCASpatial_clus_v1',label = T,repel = T)

saveRDS(integrated_sc_fil,'integrated_sc_fil_v1.rds')


#### UCASpatial_clus_v2 ----------
### divide epithelium
integrated_sc_epi <- subset(integrated_sc_fil,subset = UCASpatial_clus_v1 %in%
                              c("Epithelium"))
DimPlot(integrated_sc_epi,group.by = 'UCASpatial_clus_v1',label = T,repel = T)

DefaultAssay(integrated_sc_epi) <- 'RNA'
FeaturePlot(integrated_sc_epi,features = c('PTPRC','MKI67','C1QA','HLA-A','ITGAM','CD14','CLEC9A','ITGAX'))
DimPlot(integrated_sc_epi,group.by = 'first_around_anno',label = T,repel = T)+NoLegend()

#### recluster

DefaultAssay(integrated_sc_epi) <- 'integrated'
integrated_sc_epi <- FindNeighbors(integrated_sc_epi, dims = 1:30)
integrated_sc_epi <- FindClusters(integrated_sc_epi, resolution = 0.04)
integrated_sc_epi <- RunUMAP(integrated_sc_epi, dims = 1:30)

table(integrated_sc_epi$integrated_snn_res.0.03)
DimPlot(integrated_sc_epi,group.by = 'integrated_snn_res.0.03')

p1 <- DimPlot(integrated_sc_epi,group.by = 'integrated_snn_res.0.03')

p1+DimPlot(integrated_sc_epi,group.by = 'orig.ident')

DefaultAssay(integrated_sc_epi) <- 'RNA'
Idents(integrated_sc_epi) <- integrated_sc_epi$integrated_snn_res.0.03
mk_epi <- FindAllMarkers(integrated_sc_epi,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)

FeaturePlot(integrated_sc_epi,features = c('nCount_RNA','nFeature_RNA'),max.cutoff = 1000)
FeaturePlot(integrated_sc_epi,features = c('NNMT','DPEP1','LEFTY1'))

integrated_sc_epi$Epi_Minor <- integrated_sc_epi$integrated_snn_res.0.03
table(integrated_sc_epi$Epi_Minor)

levels(integrated_sc_epi$Epi_Minor)[1] <- 'DPEP1+ Epithelium c1' ## DPEP1
levels(integrated_sc_epi$Epi_Minor)[2] <- 'LEFTY1+ Epithelium c2' ## LEFTY1
levels(integrated_sc_epi$Epi_Minor)[3] <- 'NNMT+ Epithelium c3' ## NNMT

DimPlot(integrated_sc_epi,group.by = 'Epi_Minor')



### Combine

integrated_sc_fil$UCASpatial_clus_v2 <- integrated_sc_fil$UCASpatial_clus_v1

### epi
common_cells <- intersect(colnames(integrated_sc_fil), colnames(integrated_sc_epi))
integrated_sc_fil$UCASpatial_clus_v2 <- as.character(integrated_sc_fil$UCASpatial_clus_v2)
integrated_sc_fil$UCASpatial_clus_v2[common_cells] <- as.character(integrated_sc_epi$Epi_Minor)
table(integrated_sc_fil$UCASpatial_clus_v2)

DimPlot(integrated_sc_fil,group.by = 'UCASpatial_clus_v2',label = T,repel = T)

names(table(integrated_sc_fil$UCASpatial_clus_v2))
integrated_sc_fil$UCASpatial_clus_v2 <- factor(integrated_sc_fil$UCASpatial_clus_v2,
                                               levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                          "CD4 Naive","CD4 Tem","CD4 Treg","Th1/Th17",
                                                          "NK","ILC",
                                                          "Mast",
                                                          "Monocyte","Neutrophil","SPP1+ Mac","FOLR2+ Mac",
                                                          "cDC1","cDC2","cDC3","pDC",
                                                          "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                          "Smooth Muscle","Pericyte",
                                                          "VECs","LECs","Glial",
                                                          "B","Plasma",
                                                          "DPEP1+ Epithelium c1","LEFTY1+ Epithelium c2","NNMT+ Epithelium c3"))
DimPlot(integrated_sc_fil,group.by = 'UCASpatial_clus_v2',label = T,repel = T)

saveRDS(integrated_sc_fil,'integrated_sc_fil_v2.rds')






#### UCASpatial_clus_v3 -----
integrated_sc_T_v2 <- readRDS('/Path/data/T_new.RDS')

DimPlot(integrated_sc_T_v2,label = T,group.by = 'second_around_anno_new_2')

FeaturePlot(integrated_sc_T_v2,features = 'FOXP3')

levels(integrated_sc_T_v2$second_around_anno_new_2)
integrated_sc_T_v2@meta.data$second_around_anno_new_2[integrated_sc_T_v2@meta.data$second_around_anno_new_2 %in% c("Th1")]<-"Treg"
DimPlot(integrated_sc_T_v2,label = T,group.by = 'second_around_anno_new_2')

### Combine
integrated_sc_fil_v3 <- integrated_sc_fil
integrated_sc_fil_v3$UCASpatial_clus_v3 <- integrated_sc_fil_v3$UCASpatial_clus_v2
DimPlot(integrated_sc_fil_v3,group.by = 'UCASpatial_clus_v3')

### T
common_cells <- intersect(colnames(integrated_sc_fil_v3), colnames(integrated_sc_T_v2))
integrated_sc_fil_v3$UCASpatial_clus_v3 <- as.character(integrated_sc_fil_v3$UCASpatial_clus_v3)
all(names(integrated_sc_fil_v3$UCASpatial_clus_v3[common_cells])==names(integrated_sc_T_v2$second_around_anno_new_2[common_cells]))
integrated_sc_fil_v3$UCASpatial_clus_v3[common_cells] <- as.character(integrated_sc_T_v2$second_around_anno_new_2[common_cells])

table(integrated_sc_fil_v3$UCASpatial_clus_v3)

DimPlot(integrated_sc_fil_v3,group.by = 'UCASpatial_clus_v3',label = T,repel = T)

names(table(integrated_sc_fil_v3$UCASpatial_clus_v3))
integrated_sc_fil_v3$UCASpatial_clus_v3 <- factor(integrated_sc_fil_v3$UCASpatial_clus_v3,
                                                  levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                             "CD4 Tn","CD4 Tm","CD4 X","Tfh","Th17",'Treg',
                                                             'gdT',"NK","ILC",
                                                             "Mast",
                                                             "Monocyte","Neutrophil","SPP1+ Mac","FOLR2+ Mac",
                                                             "cDC1","cDC2","cDC3","pDC",
                                                             "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                             "Smooth Muscle","Pericyte",
                                                             "VECs","LECs","Glial",
                                                             "B","Plasma",
                                                             "DPEP1+ Epithelium c1","LEFTY1+ Epithelium c2","NNMT+ Epithelium c3"))
DimPlot(integrated_sc_fil_v3,group.by = 'UCASpatial_clus_v3',label = T,repel = T)

saveRDS(integrated_sc_fil_v3,'integrated_sc_fil_v3.rds')







#### UCASpatial_clus_v4 -----

# merge T, merge endo
integrated_sc_fil_v4 <- integrated_sc_fil_v3
integrated_sc_fil_v4$UCASpatial_clus_v4 <- integrated_sc_fil_v4$UCASpatial_clus_v3
DimPlot(integrated_sc_fil_v4,group.by = 'UCASpatial_clus_v4')

levels(integrated_sc_fil_v4$UCASpatial_clus_v4)
levels(integrated_sc_fil_v4$UCASpatial_clus_v4)[6] <- "CD4 Tm"
levels(integrated_sc_fil_v4$UCASpatial_clus_v4)[9] <- "CD8 Trm"
levels(integrated_sc_fil_v4$UCASpatial_clus_v4)[7] <- "Th1/Th17"
levels(integrated_sc_fil_v4$UCASpatial_clus_v4)[7] <- "Th1/Th17"
levels(integrated_sc_fil_v4$UCASpatial_clus_v4)[25:26] <- "Endothelium"

table(integrated_sc_fil_v4$UCASpatial_clus_v4)

DimPlot(integrated_sc_fil_v4,group.by = 'UCASpatial_clus_v4',label = T,repel = T)

names(table(integrated_sc_fil_v4$UCASpatial_clus_v4))
integrated_sc_fil_v4$UCASpatial_clus_v4 <- factor(integrated_sc_fil_v4$UCASpatial_clus_v4,
                                                  levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                             "CD4 Tn","CD4 Tm","Tfh","Th1/Th17",'Treg',
                                                             "NK","ILC",
                                                             "Mast",
                                                             "Monocyte","Neutrophil","SPP1+ Mac","FOLR2+ Mac",
                                                             "cDC1","cDC2","cDC3","pDC",
                                                             "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                             "Smooth Muscle","Pericyte",
                                                             "Endothelium","Glial",
                                                             "B","Plasma",
                                                             "DPEP1+ Epithelium c1","LEFTY1+ Epithelium c2","NNMT+ Epithelium c3"))
DimPlot(integrated_sc_fil_v4,group.by = 'UCASpatial_clus_v4',label = T,repel = T)

saveRDS(integrated_sc_fil_v4,'integrated_sc_fil_v4.rds')




#### UCASpatial_clus_v5 -----

# recluster epi

DefaultAssay(integrated_sc_epi) <- 'integrated'
integrated_sc_epi <- FindNeighbors(integrated_sc_epi, dims = 1:30)
integrated_sc_epi <- FindClusters(integrated_sc_epi, resolution = 0.05)
integrated_sc_epi <- RunUMAP(integrated_sc_epi, dims = 1:30)

table(integrated_sc_epi$integrated_snn_res.0.05)
DimPlot(integrated_sc_epi,group.by = 'integrated_snn_res.0.05')

p1 <- DimPlot(integrated_sc_epi,group.by = 'integrated_snn_res.0.05')

p1+DimPlot(integrated_sc_epi,group.by = 'orig.ident')

DefaultAssay(integrated_sc_epi) <- 'RNA'
Idents(integrated_sc_epi) <- integrated_sc_epi$integrated_snn_res.0.05
mk_epi <- FindAllMarkers(integrated_sc_epi,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)

FeaturePlot(integrated_sc_epi,features = c('nCount_RNA','nFeature_RNA'),max.cutoff = 1000)
p1+FeaturePlot(integrated_sc_epi,features = c('LEFTY1','TM4SF1','SPINK4','MUC2','REG1B','REG1A'))

integrated_sc_epi$Epi_Minor_2 <- integrated_sc_epi$integrated_snn_res.0.05
table(integrated_sc_epi$Epi_Minor_2)

levels(integrated_sc_epi$Epi_Minor_2)[1] <- 'Epithelium c1' ## DPEP1
levels(integrated_sc_epi$Epi_Minor_2)[2] <- 'LEFTY1+ Epithelium c2' ## LEFTY1
levels(integrated_sc_epi$Epi_Minor_2)[3] <- 'TM4SF1+ Epithelium c3' ## NNMT
levels(integrated_sc_epi$Epi_Minor_2)[4] <- 'MUC2+ Epithelium c4' ## LEFTY1
levels(integrated_sc_epi$Epi_Minor_2)[5] <- 'REG1A+ Epithelium c5' ## NNMT

DimPlot(integrated_sc_epi,group.by = 'Epi_Minor_2',label = T,repel = T)+
  FeaturePlot(integrated_sc_epi,features = c('LEFTY1','TM4SF1','SPINK4','MUC2','REG1B','REG1A'))
ggsave('integrated_sc_epi_5clusters_umap_feature.pdf',height = 6,width = 12)

saveRDS(integrated_sc_epi,'integrated_sc_epi.rds')


#
integrated_sc_fil_v5 <- integrated_sc_fil_v4
integrated_sc_fil_v5$UCASpatial_clus_v5 <- integrated_sc_fil_v5$UCASpatial_clus_v4
DimPlot(integrated_sc_fil_v5,group.by = 'UCASpatial_clus_v5')

### epi
common_cells <- intersect(colnames(integrated_sc_fil_v5), colnames(integrated_sc_epi))
integrated_sc_fil_v5$UCASpatial_clus_v5 <- as.character(integrated_sc_fil_v5$UCASpatial_clus_v5)
integrated_sc_fil_v5$UCASpatial_clus_v5[common_cells] <- as.character(integrated_sc_epi$Epi_Minor_2[common_cells])
table(integrated_sc_fil_v5$UCASpatial_clus_v5)

DimPlot(integrated_sc_fil_v5,group.by = 'UCASpatial_clus_v5',label = T,repel = T)

names(table(integrated_sc_fil_v5$UCASpatial_clus_v5))
integrated_sc_fil_v5$UCASpatial_clus_v5 <- factor(integrated_sc_fil_v5$UCASpatial_clus_v5,
                                                  levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                             "CD4 Tn","CD4 Tm","Tfh","Th1/Th17",'Treg',
                                                             "NK","ILC",
                                                             "Mast",
                                                             "Monocyte","Neutrophil","SPP1+ Mac","FOLR2+ Mac",
                                                             "cDC1","cDC2","cDC3","pDC",
                                                             "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                             "Smooth Muscle","Pericyte",
                                                             "Endothelium","Glial",
                                                             "B","Plasma",
                                                             'Epithelium c1','LEFTY1+ Epithelium c2','TM4SF1+ Epithelium c3',
                                                             'MUC2+ Epithelium c4','REG1A+ Epithelium c5'))
DimPlot(integrated_sc_fil_v5,group.by = 'UCASpatial_clus_v5',label = T,repel = T)






#### re-integrate myeloid ------


# .libPaths(c("/data/software/R/R_package_4","/usr/local/lib/R/site-library","/usr/local/lib/R/library"))
# # library(Seurat)
CRC_obj<-readRDS("/Path/data/CRC_myeloid_list.RDS")
features <- SelectIntegrationFeatures(object.list = CRC_obj,nfeatures = 2000)
CRC_ICB.anchors <- FindIntegrationAnchors(object.list = CRC_obj, anchor.features = features)
CRC_ICB.combined <- IntegrateData(anchorset = CRC_ICB.anchors)
saveRDS(CRC_ICB.anchors,file = "/Path/data/myeloid_integrated_anchors.RDS")
saveRDS(CRC_ICB.combined,file = "/Path/data/myeloid_integrated_data.RDS")

reintegrated_sc_myeloid <- readRDS("/Path/data/myeloid_integrated_data.RDS")
DefaultAssay(reintegrated_sc_myeloid) <- 'integrated'

#### recluster
sc_CRC_mereintegrated_sc_myeloidrge <- FindVariableFeatures(reintegrated_sc_myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(reintegrated_sc_myeloid)
reintegrated_sc_myeloid <- ScaleData(reintegrated_sc_myeloid, features = all.genes)
reintegrated_sc_myeloid <- RunPCA(reintegrated_sc_myeloid, features = VariableFeatures(object = sc_CRC_merge))

reintegrated_sc_myeloid <- FindNeighbors(reintegrated_sc_myeloid, dims = 1:30)
reintegrated_sc_myeloid <- FindClusters(reintegrated_sc_myeloid, resolution = seq(0.1,1,0.1))
reintegrated_sc_myeloid <- RunUMAP(reintegrated_sc_myeloid, dims = 1:30)

table(reintegrated_sc_myeloid$integrated_snn_res.0.3)

DimPlot(reintegrated_sc_myeloid,group.by = 'integrated_snn_res.0.4')

FeaturePlot(reintegrated_sc_myeloid,
            features = c('FOLR2','CD163L1','FCN1','THBS1',
                         'CD1C','FCER1A', 'CLEC10A','SPP1','FN1',
                         'CLEC9A','XCR1','LAMP3','CCR7',
                         'CSF3R','S100A8','PDPN','CXCL9','CXCL10','IFNAR1'))

Idents(reintegrated_sc_myeloid) <- reintegrated_sc_myeloid$integrated_snn_res.0.4
mk_mye_v2 <- FindAllMarkers(reintegrated_sc_myeloid,min.pct = 0.2,logfc.threshold = 0.25,min.diff.pct = 0.1,
                            only.pos = T)

table(reintegrated_sc_myeloid$myeloid_Minor,reintegrated_sc_myeloid$integrated_snn_res.0.4)

DefaultAssay(reintegrated_sc_myeloid) <- 'RNA'
VlnPlot(reintegrated_sc_myeloid,group.by = 'integrated_snn_res.0.4',pt.size = 0,
        features = c('PTPRC','MKI67','LAMP3','S100A8','ITGAM','CCR7','CLEC9A',
                     'ITGAX','CD1C','SPP1','CSF3R','LILRB2'))

DimPlot(reintegrated_sc_myeloid)
DotPlot(reintegrated_sc_myeloid,features = c('MKI67','PTPRC'))


reintegrated_sc_myeloid$myeloid_Minor <- droplevels(reintegrated_sc_myeloid$integrated_snn_res.0.4)

levels(reintegrated_sc_myeloid$myeloid_Minor)[1] <- 'C1QC+ Mac'
levels(reintegrated_sc_myeloid$myeloid_Minor)[2] <- 'Monocyte'
levels(reintegrated_sc_myeloid$myeloid_Minor)[3] <- 'cDC2'
levels(reintegrated_sc_myeloid$myeloid_Minor)[4] <- 'CD45- myeloid'      #--> filter
levels(reintegrated_sc_myeloid$myeloid_Minor)[5] <- 'FN1+ Mac'
levels(reintegrated_sc_myeloid$myeloid_Minor)[6] <- 'PDPN+ Mac'
levels(reintegrated_sc_myeloid$myeloid_Minor)[7] <- 'proliferation myeloid' #--> filter
levels(reintegrated_sc_myeloid$myeloid_Minor)[8] <- 'cDC1'
levels(reintegrated_sc_myeloid$myeloid_Minor)[9] <- 'Neutrophil'
levels(reintegrated_sc_myeloid$myeloid_Minor)[10] <- 'cDC3'
levels(reintegrated_sc_myeloid$myeloid_Minor)[11] <- 'pDC'

levels(reintegrated_sc_myeloid$myeloid_Minor)
DimPlot(reintegrated_sc_myeloid,group.by = 'myeloid_Minor',label = T)
DotPlot(reintegrated_sc_myeloid,group.by = 'myeloid_Minor',features = c('MKI67','PTPRC'))


# filter prolif myeloid and CD45- myeloid
reintegrated_sc_myeloid_fil <- subset(reintegrated_sc_myeloid,subset=
                                        integrated_snn_res.0.4 %in% c(0,1,2,4,5,7:10))
reintegrated_sc_myeloid_fil <- RunUMAP(reintegrated_sc_myeloid_fil, dims = 1:30)
DimPlot(reintegrated_sc_myeloid_fil,group.by = 'integrated_snn_res.0.4')

mk_mye_v2 <- FindAllMarkers(reintegrated_sc_myeloid_fil,min.pct = 0.2,logfc.threshold = 0.25,min.diff.pct = 0.1,
                            only.pos = T)

VlnPlot(reintegrated_sc_myeloid_fil,group.by = 'integrated_snn_res.0.4',pt.size = 0,
        features = c('PTPRC','MKI67','LAMP3','S100A8','ITGAM','CCR7','CLEC9A',
                     'ITGAX','CD1C','SPP1','CSF3R','LILRB2'))

Dotplot(DotPlot(reintegrated_sc_myeloid_fil,group.by = 'integrated_snn_res.0.4',
                features = c('SELENOP','IGF1','FOLR2','CD163L1','FCN1','THBS1',
                             'CD1C','FCER1A', 'CLEC10A','SPP1','FN1',
                             'CLEC9A','XCR1','LAMP3','CCR7','GZMB','IL3RA' ,'LILRB4' ,'CLEC4C',
                             'CSF3R','S100A8','PDPN','CXCL9','CXCL10','IFNAR1'))$data)

FeaturePlot(reintegrated_sc_myeloid,features = c('CSF3R','CXCR2','FCGR3A'))
FeaturePlot(reintegrated_sc_myeloid,features = c('CD163','MRC1','SPP1','PDPN','IGF1','C1QA','IL1B'))
FeaturePlot(reintegrated_sc_myeloid_fil,features = c('LYVE1','APOE','C1QC','SPP1','FN1','SEPP1',
                                                     'CXCL9','CXCL10','FOLR2'))
FeaturePlot(reintegrated_sc_myeloid,features = c('nCount_RNA','nFeature_RNA'),max.cutoff = 1000)
FeaturePlot(reintegrated_sc_myeloid,features = c('PTPRC','C1QA'))

reintegrated_sc_myeloid_fil$myeloid_Minor <- droplevels(reintegrated_sc_myeloid_fil$integrated_snn_res.0.4)

levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[1] <- 'C1QC+ Mac'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[2] <- 'Monocyte'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[3] <- 'cDC2'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[4] <- 'FN1+ Mac'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[5] <- 'PDPN+ Mac'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[6] <- 'cDC1'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[7] <- 'Neutrophil'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[8] <- 'cDC3'
levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[9] <- 'pDC'

levels(reintegrated_sc_myeloid_fil$myeloid_Minor)


FeaturePlot(reintegrated_sc_myeloid_fil,
            features = c('IGF1','FN1','PDPN','FOLR2','SPP1',
                         'THBS1','CD1C','FCER1A',
                         'CLEC9A','CSF3R','S100A8','LAMP3','CCR7','GZMB',
                         'CXCL9','CXCL10'))
ggsave('myeloid_reintegrated_feature.pdf',height = 12,width = 16)

p1 <- DimPlot(reintegrated_sc_myeloid_fil,group.by = 'myeloid_Minor',label = T)
p2 <- FeaturePlot(reintegrated_sc_myeloid_fil,
                  features = c('IGF1','FN1','PDPN','FOLR2','SPP1',
                               'THBS1','CD1C','FCER1A',
                               'CLEC9A','CSF3R','S100A8','LAMP3','CCR7','GZMB',
                               'CXCL9','CXCL10'))

ggsave('myeloid_reintegrated_umap_feature.pdf',height = 8,width = 20)

DefaultAssay(reintegrated_sc_myeloid_fil) <- 'RNA'
Dotplot(DotPlot(reintegrated_sc_myeloid_fil,group.by = 'myeloid_Minor',
                features = c('IGF1','FOLR2','FCN1','THBS1',
                             'CD1C','FCER1A', 'CLEC10A','SPP1','FN1','PDPN',
                             'CLEC9A','XCR1','CSF3R','S100A8','LAMP3','CCR7','GZMB',
                             'IL3RA' ,'LILRB4' ,'CLEC4C',
                             'CXCL9','CXCL10','IFNAR1'))$data)
ggsave('myeloid_reintegrated_umap_dotplot.pdf',height = 5,width = 8)

saveRDS(reintegrated_sc_myeloid_fil,'reintegrated_sc_myeloid_fil.RDS')



#### UCASpatial_clus_v7 -------

integrated_sc_v7 <- integrated_sc
integrated_sc_v7$UCASpatial_clus_v7 <- integrated_sc_v7$UCASpatial_clus_v1
DimPlot(integrated_sc_v7,group.by = 'UCASpatial_clus_v7')

### total tag
common_cells <- intersect(colnames(integrated_sc_fil_v5), colnames(integrated_sc_v7))
integrated_sc_v7$UCASpatial_clus_v7 <- as.character(integrated_sc_v7$UCASpatial_clus_v7)
integrated_sc_v7$UCASpatial_clus_v7[common_cells] <-
  as.character(integrated_sc_fil_v5$UCASpatial_clus_v5[common_cells])
table(integrated_sc_v7$UCASpatial_clus_v7)

### myeloid
common_cells <- intersect(colnames(integrated_sc_v7), colnames(reintegrated_sc_myeloid))
integrated_sc_v7$UCASpatial_clus_v7 <- as.character(integrated_sc_v7$UCASpatial_clus_v7)
integrated_sc_v7$UCASpatial_clus_v7[common_cells] <-
  as.character(reintegrated_sc_myeloid$myeloid_Minor[common_cells])
table(integrated_sc_v7$UCASpatial_clus_v7)

DimPlot(integrated_sc_v7,group.by = 'UCASpatial_clus_v7',label = T,repel = T)

### filter proliferation myeloid and CD45- myeloid
integrated_sc_fil_v7 <- subset(integrated_sc_v7,subset=UCASpatial_clus_v7%in%
                                 names(table(integrated_sc_v7$UCASpatial_clus_v7))[c(1:3,5:29,32:37)])
table(integrated_sc_fil_v7$UCASpatial_clus_v7)

DimPlot(integrated_sc_fil_v7,group.by = 'UCASpatial_clus_v7',label = T,repel = T)


names(table(integrated_sc_fil_v7$UCASpatial_clus_v7))

integrated_sc_fil_v7$UCASpatial_clus_v7 <- factor(integrated_sc_fil_v7$UCASpatial_clus_v7,
                                                  levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                             "CD4 Tn","CD4 Tm","Tfh","Th1/Th17",'Treg',
                                                             "NK","ILC",
                                                             "Mast",
                                                             "Monocyte","Neutrophil",
                                                             "C1QC+ Mac",'FN1+ Mac','PDPN+ Mac',
                                                             "cDC1","cDC2","cDC3","pDC",
                                                             "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                             "Smooth Muscle","Pericyte",
                                                             "Endothelium","Glial",
                                                             "B","Plasma",
                                                             'Epithelium c1','Epithelium c2','Epithelium c3',
                                                             'Epithelium c4','Epithelium c5'))
DimPlot(integrated_sc_fil_v7,group.by = 'UCASpatial_clus_v7',label = T,repel = T)

saveRDS(integrated_sc_fil_v7,'integrated_sc_fil_v7.rds')

integrated_sc_fil_v8<-readRDS('integrated_sc_fil_v7.rds')
levels(integrated_sc_fil_v8$UCASpatial_clus_v7) <- c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                     "CD4 Tn","CD4 Tm","CD4 Tfh","Th1/Th17",'CD4 Treg',
                                                     "NK","ILC",
                                                     "Mast",
                                                     "cMo","Neutrophil",
                                                     "HES1hi C1QChi Mac",'TREM2+ APOC1hi Mac','TREM2+ PDPN+ Mac',
                                                     "cDC1","cDC2","CCR7+ DC","pDC",
                                                     "CD34+ iFib","FAP+ Myofib","POSTN+ Myofib",
                                                     "Smooth Muscle","Pericyte",
                                                     "Endothelial","Glial",
                                                     "B","Plasma",
                                                     'Epithelial c1','Epithelial c2','Epithelial c3',
                                                     'Epithelal c4','Epithelal c5')