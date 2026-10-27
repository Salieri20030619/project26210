#Fig3b
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


#Fig3c
library(pacman)
p_unload(ggplot2)
p_load(ggplot2)
C2L_Seurat<-readRDS("/Path/data/C2L_results_ICB.rds")
C2L_Seurat<-C2L_Seurat[[1]]
RCTD_Seurat<-readRDS("/Path/data/RCTD_results_ICB.rds")
RCTD_Seurat<-RCTD_Seurat[[1]]
SPOTlight_Seurat<-readRDS("/Path/data/SPOTlight_results_ICB.rds")
SPOTlight_Seurat<-SPOTlight_Seurat[[1]]
CARD_Seurat<-readRDS("/Path/data/CARD_results_ICB.rds")
CARD_Seurat<-CARD_Seurat[[1]]
stereoscope_Seurat<-readRDS("/Path/data/stereoscope_results_ICB.rds")
stereoscope_Seurat<-stereoscope_Seurat[[1]]
ICB_sliceA<-readRDS("/Path/data/ICB_section_A_region.RDS")
Spotiphy_Seurat<-readRDS("/Path/data/Spotiphy_results_ICB.rds")
Spotiphy_Seurat<-Spotiphy_Seurat[[1]]
plot_data<-readRDS("/Path/data/ICB_region.RDS")
total_panel_statistics<- readRDS(file = "/Path/data/total_panel_statistics.RDS")


colnames(C2L_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L_Seurat@meta.data))))))))))))
colnames(RCTD_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD_Seurat@meta.data))))))))))))
colnames(SPOTlight_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight_Seurat@meta.data))))))))))))
colnames(CARD_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD_Seurat@meta.data))))))))))))
colnames(stereoscope_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope_Seurat@meta.data))))))))))))
colnames(ICB_sliceA@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(ICB_sliceA@meta.data))))))))))))
colnames(Spotiphy_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy_Seurat@meta.data))))))))))))

Tumor1_spot<-rownames(plot_data)[plot_data$tumor_region == "Tumor1"]
Tumor2_spot<-rownames(plot_data)[plot_data$tumor_region == "Tumor2"]
TLS_spot<-rownames(plot_data)[plot_data$tumor_region == "TLS"]
Stroma_spot<-rownames(plot_data)[plot_data$tumor_region == "Stroma"]

C2L_Seurat@meta.data$tumor_region<-"Remaining"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% TLS_spot]<-"TLS"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

RCTD_Seurat@meta.data$tumor_region<-"Remaining"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% TLS_spot]<-"TLS"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

SPOTlight_Seurat@meta.data$tumor_region<-"Remaining"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% TLS_spot]<-"TLS"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

CARD_Seurat@meta.data$tumor_region<-"Remaining"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% TLS_spot]<-"TLS"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

stereoscope_Seurat@meta.data$tumor_region<-"Remaining"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% TLS_spot]<-"TLS"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

Spotiphy_Seurat@meta.data$tumor_region<-"Remaining"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% TLS_spot]<-"TLS"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

Ablation_Seurat@meta.data$tumor_region

ICB_sliceA@meta.data$tumor_region<-"Remaining"
ICB_sliceA@meta.data$tumor_region[rownames(ICB_sliceA@meta.data) %in% Tumor1_spot]<-"Tumor1"
ICB_sliceA@meta.data$tumor_region[rownames(ICB_sliceA@meta.data) %in% Tumor2_spot]<-"Tumor2"

ICB_sliceA@meta.data$TLS_score<-rowMeans(ICB_sliceA@meta.data[,c("B","Tfh","cDC3","CD4+ Tn")])-rowMeans(ICB_sliceA@meta.data[,c("CD8+ Trm","FAP+ Myofib","FN1+ Mac","NK","Epithelial c3","Epithelial c4","Epithelial c5")])
#ICB_sliceA<-UpdateSeuratObject(ICB_sliceA)

ICB_sliceA@meta.data$TLS_defination<-"Non_TLS"
ICB_sliceA@meta.data$TLS_defination[ICB_sliceA@meta.data$TLS_score > 0.07]<-"TLS"
ICB_sliceA@meta.data$tumor_region[ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] > 12000 & ICB_sliceA@meta.data$TLS_defination == "TLS"]<- "TLS"

ICB_sliceA@meta.data$stroma_defination<-"Remain"
ICB_sliceA@meta.data$stroma_defination[ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] > 11600 &
                                         ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] < 12300 &
                                         ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,1] > 19700 &
                                         ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,1] < 20700]<-"Stroma_selected"
ICB_sliceA@meta.data$tumor_region[ICB_sliceA@meta.data$stroma_defination == "Stroma_selected"]<-"Stroma_selected2"

DefaultAssay(ICB_sliceA) <- "Annotation"
SpatialPlot(ICB_sliceA,group.by = "stroma_defination",pt.size.factor = 2)
#SpatialPlot(ICB_sliceA[,ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] > 12000],group.by = "TLS_defination",pt.size.factor = 3)
SpatialPlot(ICB_sliceA,group.by = "tumor_region",pt.size.factor = 3,cols = c("white",""))

C2L_Seurat@meta.data$TREM2_mac<-rowSums(C2L_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
RCTD_Seurat@meta.data$TREM2_mac<-rowSums(RCTD_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
SPOTlight_Seurat@meta.data$TREM2_mac<-rowSums(SPOTlight_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
CARD_Seurat@meta.data$TREM2_mac<-rowSums(CARD_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
stereoscope_Seurat@meta.data$TREM2_mac<-rowSums(stereoscope_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
ICB_sliceA@meta.data$TREM2_mac<-rowSums(ICB_sliceA@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
Spotiphy_Seurat@meta.data$TREM2_mac<-rowSums(Spotiphy_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])

ST_tumor_combine_C2L<-aggregate(C2L_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                        "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                list(Region = C2L_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_C2L<-ST_tumor_combine_C2L[2:4,]
ST_tumor_combine_C2L<-reshape2::melt(ST_tumor_combine_C2L)

ST_tumor_combine_RCTD<-aggregate(RCTD_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                          "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                 list(Region = RCTD_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_RCTD<-ST_tumor_combine_RCTD[2:4,]
ST_tumor_combine_RCTD<-reshape2::melt(ST_tumor_combine_RCTD)

ST_tumor_combine_SPOTlight<-aggregate(SPOTlight_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                                    "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                      list(Region = SPOTlight_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_SPOTlight<-ST_tumor_combine_SPOTlight[2:4,]
ST_tumor_combine_SPOTlight<-reshape2::melt(ST_tumor_combine_SPOTlight)

ST_tumor_combine_UCAS<-aggregate(ICB_sliceA@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                         "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                 list(Region = ICB_sliceA@meta.data$tumor_region),mean)

ST_tumor_combine_UCAS<-ST_tumor_combine_UCAS[3:5,]
ST_tumor_combine_UCAS<-reshape2::melt(ST_tumor_combine_UCAS)

ST_tumor_combine_CARD<-aggregate(CARD_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                          "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                 list(Region = CARD_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_CARD<-ST_tumor_combine_CARD[2:4,]
ST_tumor_combine_CARD<-reshape2::melt(ST_tumor_combine_CARD)

ST_tumor_combine_stereoscope<-aggregate(stereoscope_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                                        "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                        list(Region = stereoscope_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_stereoscope<-ST_tumor_combine_stereoscope[2:4,]
ST_tumor_combine_stereoscope<-reshape2::melt(ST_tumor_combine_stereoscope)

ST_tumor_combine_Spotiphy<-aggregate(Spotiphy_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                                  "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                     list(Region = Spotiphy_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_Spotiphy<-ST_tumor_combine_Spotiphy[2:4,]
ST_tumor_combine_Spotiphy<-reshape2::melt(ST_tumor_combine_Spotiphy)

ST_tumor_combine_Ablation


ST_tumor_combine_UCAS$Region == ST_tumor_combine_C2L$Region
ST_tumor_combine_3<-cbind(ST_tumor_combine_UCAS,ST_tumor_combine_C2L[,3],ST_tumor_combine_RCTD[,3],
                          ST_tumor_combine_SPOTlight[,3],ST_tumor_combine_CARD[,3],ST_tumor_combine_stereoscope[,3],ST_tumor_combine_Spotiphy[,3])

levels(ST_tumor_combine_3$variable)
gsub("C1QC mac","C1QC+ mac",gsub("CD4+ Treg","Treg",gsub("CD4+ Tfh","Tfh",gsub("CD4","CD4+",gsub("CD8","CD8+",gsub(" cells","",total_panel_statistics$celltype))))))

levels(ST_tumor_combine_3$variable)[c(8,13)]<-c("TREM2+ mac","CCR7 DC")
levels(ST_tumor_combine_3$variable)[c(9)]<-"C1QC+ mac"
ST_tumor_combine_3<-ST_tumor_combine_3[ST_tumor_combine_3$variable != "CD4+ Tm",]
ST_tumor_combine_3$variable == total_panel_statistics$celltype
ST_tumor_combine_3$Region == total_panel_statistics$Region

#ST_tumor_combine_3$variable <- mIHC_tumor_combine$celltype
ST_tumor_combine_3$proportion_mIHC <- total_panel_statistics$Proportion
colnames(ST_tumor_combine_3)[3:9]<-c("proportion_UCAS","proportion_C2L","proportion_RCTD",
                                     "proportion_SPOTlight","proportion_CARD","proportion_stereoscope","proportion_Spotiphy")
colnames(ST_tumor_combine_3)[2]<-"celltype"
ST_tumor_combine_3$proportion_Ablation = ST_tumor_combine_Ablation$value[ST_tumor_combine_Ablation$Region != "Stroma"]
ST_tumor_combine_Ablation<-ST_tumor_combine_Ablation[ST_tumor_combine_Ablation$variable != "CD4+ Tm",]


#ST_tumor_combine_2$proportion_mIHC = mIHC_tumor_combine$proportion
#ST_tumor_combine_2$proportion_mIHC<-mIHC_tumor_combine$proportion
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)

ST_tumor_combine_3 <- readRDS("/Path/data/ST_tumor_combine_3.RDS")

ST_tumor_combine_3_tumor<-ST_tumor_combine_3[ST_tumor_combine_3$Region %in% c("Tumor1","Tumor2"),]
p1<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_UCAS*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (UCASpatial)",y="Cell Proportion%DAPI (mIHC)",title = "UCASpatial")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))
p1
p2<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_C2L*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (cell2location)",y="Cell Proportion%DAPI (mIHC)",title = "cell2location")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))

p3<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_RCTD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (RCTD)",y="Cell Proportion%DAPI (mIHC)",title = "RCTD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))

p4<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_SPOTlight*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (SPOTlight)",y="Cell Proportion%DAPI (mIHC)",title = "SPOTlight")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,15),ylim = c(0,15))

p5<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_CARD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (CARD)",y="Cell Proportion%DAPI (mIHC)",title = "CARD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,40),ylim = c(0,40))

p6<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_stereoscope*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (stereoscope)",y="Cell Proportion%DAPI (mIHC)",title = "stereoscope")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))

p7<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_Spotiphy*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Spotiphy)",y="Cell Proportion%DAPI (mIHC)",title = "Spotiphy")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,7),ylim = c(0,7))

p8<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_Ablation*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Ablation)",y="Cell Proportion%DAPI (mIHC)",title = "Ablation")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,10),ylim = c(0,10))

library(patchwork)
p1+p8+p2+p3+p4+p5+p6+p7+plot_layout(ncol=4)
library(Metrics)
PCC_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                   c(cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_UCAS)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_C2L)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_RCTD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_SPOTlight)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_CARD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_stereoscope)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_Spotiphy)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_Ablation)[4]$estimate))
colnames(PCC_df)<-c("method","PCC")
PCC_df$method<-factor(PCC_df$method,levels = c("UCASpatial","RCTD","Ablation","Spotiphy","CARD","cell2location","stereoscope","SPOTlight"))
levels(PCC_df$method)[3] <- "UCASpatial_Ablation"
p1<-ggplot(PCC_df,aes(x=method,y=PCC,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#E6D77A","#CBEBC5","#d9d9d9","#add26a","#b37cb3","#7eabca"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1

RMSE_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                    c(rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_UCAS*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_C2L*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_RCTD*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_SPOTlight*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_CARD*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_stereoscope*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_Spotiphy*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_Ablation*100)))
colnames(RMSE_df)<-c("method","RMSE")
RMSE_df$method<-factor(RMSE_df$method,levels = c("UCASpatial","RCTD","stereoscope","cell2location","Spotiphy","Ablation","SPOTlight","CARD"))
p2<-ggplot(RMSE_df,aes(x=method,y=RMSE,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#b37cb3","#add26a","#CBEBC5","#E6D77A","#7eabca","#d9d9d9"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
p1+p2


ST_tumor_combine_3_TLS<-ST_tumor_combine_3[ST_tumor_combine_3$Region %in% c("TLS"),]
p1<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_UCAS*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (UCASpatial)",y="Cell Proportion%DAPI (mIHC)",title = "UCASpatial")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))
p1
p2<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_C2L*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (cell2location)",y="Cell Proportion%DAPI (mIHC)",title = "cell2location")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,60),ylim = c(0,60))

p3<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_RCTD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (RCTD)",y="Cell Proportion%DAPI (mIHC)",title = "RCTD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p4<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_SPOTlight*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (SPOTlight)",y="Cell Proportion%DAPI (mIHC)",title = "SPOTlight")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,30),ylim = c(0,30))

p5<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_CARD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (CARD)",y="Cell Proportion%DAPI (mIHC)",title = "CARD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p6<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_stereoscope*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (stereoscope)",y="Cell Proportion%DAPI (mIHC)",title = "stereoscope")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p7<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_Spotiphy*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Spotiphy)",y="Cell Proportion%DAPI (mIHC)",title = "Spotiphy")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p8<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_Ablation*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Ablation)",y="Cell Proportion%DAPI (mIHC)",title = "Ablation")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

library(patchwork)
p1+p8+p2+p3+p4+p5+p6+p7+plot_layout(ncol=4)
library(Metrics)
PCC_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                   c(cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_UCAS)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_C2L)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_RCTD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_SPOTlight)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_CARD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_stereoscope)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_Spotiphy)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_Ablation)[4]$estimate))
colnames(PCC_df)<-c("method","PCC")
PCC_df$method<-factor(PCC_df$method,levels = c("UCASpatial","Spotiphy","cell2location","RCTD","stereoscope","Ablation","SPOTlight","CARD"))
levels(PCC_df$method)[6] <- "UCASpatial_Ablation"
p1<-ggplot(PCC_df,aes(x=method,y=PCC,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#CBEBC5","#add26a","#f4b062","#b37cb3","#E6D77A","#7eabca","#d9d9d9"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1

RMSE_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                    c(rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_UCAS*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_C2L*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_RCTD*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_SPOTlight*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_CARD*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_stereoscope*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_Spotiphy*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_Ablation*100)))
colnames(RMSE_df)<-c("method","RMSE")
RMSE_df$method<-factor(RMSE_df$method,levels = c("UCASpatial","Spotiphy","RCTD","stereoscope","Ablation","CARD","SPOTlight","cell2location"))
p2<-ggplot(RMSE_df,aes(x=method,y=RMSE,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#CBEBC5","#f4b062","#b37cb3","#E6D77A","#d9d9d9","#7eabca","#add26a"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
p1+p2


# T cell in tumor
ST_tumor_combine_Tcell<-ST_tumor_combine_3_tumor[ST_tumor_combine_3_tumor$celltype %in% c("CD8+ Tem","CD8+ Tex","CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg"),]
PCC_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                   c(cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_UCAS)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_C2L)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_RCTD)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_SPOTlight)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_CARD)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_stereoscope)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_Spotiphy)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_Ablation)[4]$estimate))
colnames(PCC_df)<-c("method","PCC")
PCC_df$method<-factor(PCC_df$method,levels = c("UCASpatial","RCTD","cell2location","stereoscope","Spotiphy","CARD","SPOTlight","Ablation"))
levels(PCC_df$method)[8] <- "UCASpatial_Ablation"
p1<-ggplot(PCC_df,aes(x=method,y=PCC,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#add26a","#b37cb3","#CBEBC5","#d9d9d9","#7eabca","#E6D77A"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1

RMSE_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                    c(rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_UCAS*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_C2L*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_RCTD*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_SPOTlight*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_CARD*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_stereoscope*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_Spotiphy*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_Ablation*100)))
colnames(RMSE_df)<-c("method","RMSE")
RMSE_df$method<-factor(RMSE_df$method,levels = c("UCASpatial","RCTD","stereoscope","cell2location","Spotiphy","Ablation","SPOTlight","CARD"))
levels(RMSE_df$method)[6] <- "UCASpatial_Ablation"
p2<-ggplot(RMSE_df,aes(x=method,y=RMSE,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#b37cb3","#add26a","#CBEBC5","#E6D77A","#7eabca","#d9d9d9"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
p1+p2

#Fig3d
library(ggthemes)
theme_base()
theme_set(theme_cowplot())
Spotiphy<-readRDS("/Path/data/Spotiphy_results.rds")


SPOTlight<-readRDS("/Path/data/SPOTlight_results.rds")
colnames(SPOTlight[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[1]]@meta.data))))))))))))
colnames(SPOTlight[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[2]]@meta.data))))))))))))
colnames(SPOTlight[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[3]]@meta.data))))))))))))
colnames(SPOTlight[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[4]]@meta.data))))))))))))
colnames(SPOTlight[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[5]]@meta.data))))))))))))
colnames(SPOTlight[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[6]]@meta.data))))))))))))
colnames(SPOTlight[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[7]]@meta.data))))))))))))
colnames(SPOTlight[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[8]]@meta.data))))))))))))

UCASpatial<-readRDS("/Path/data/UCASpatial_results.rds")
colnames(UCASpatial[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[1]]@meta.data))))))))))))
colnames(UCASpatial[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[2]]@meta.data))))))))))))
colnames(UCASpatial[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[3]]@meta.data))))))))))))
colnames(UCASpatial[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[4]]@meta.data))))))))))))
colnames(UCASpatial[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[5]]@meta.data))))))))))))
colnames(UCASpatial[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[6]]@meta.data))))))))))))
colnames(UCASpatial[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[7]]@meta.data))))))))))))
colnames(UCASpatial[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[8]]@meta.data))))))))))))

C2L<-readRDS("/Path/data/C2L_results.rds")
colnames(C2L[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[1]]@meta.data))))))))))))
colnames(C2L[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[2]]@meta.data))))))))))))
colnames(C2L[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[3]]@meta.data))))))))))))
colnames(C2L[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[4]]@meta.data))))))))))))
colnames(C2L[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[5]]@meta.data))))))))))))
colnames(C2L[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[6]]@meta.data))))))))))))
colnames(C2L[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[7]]@meta.data))))))))))))
colnames(C2L[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[8]]@meta.data))))))))))))

CARD<-readRDS("/Path/data/CARD_results.rds")
colnames(CARD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[1]]@meta.data))))))))))))
colnames(CARD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[2]]@meta.data))))))))))))
colnames(CARD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[3]]@meta.data))))))))))))
colnames(CARD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[4]]@meta.data))))))))))))
colnames(CARD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[5]]@meta.data))))))))))))
colnames(CARD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[6]]@meta.data))))))))))))
colnames(CARD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[7]]@meta.data))))))))))))
colnames(CARD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[8]]@meta.data))))))))))))

RCTD<-readRDS("/Path/data/RCTD_results.rds")
colnames(RCTD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[1]]@meta.data))))))))))))
colnames(RCTD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[2]]@meta.data))))))))))))
colnames(RCTD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[3]]@meta.data))))))))))))
colnames(RCTD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[4]]@meta.data))))))))))))
colnames(RCTD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[5]]@meta.data))))))))))))
colnames(RCTD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[6]]@meta.data))))))))))))
colnames(RCTD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[7]]@meta.data))))))))))))
colnames(RCTD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[8]]@meta.data))))))))))))


stereoscope<-readRDS("/Path/data/stereoscope_results.rds")
colnames(stereoscope[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[1]]@meta.data))))))))))))
colnames(stereoscope[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[2]]@meta.data))))))))))))
colnames(stereoscope[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[3]]@meta.data))))))))))))
colnames(stereoscope[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[4]]@meta.data))))))))))))
colnames(stereoscope[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[5]]@meta.data))))))))))))
colnames(stereoscope[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[6]]@meta.data))))))))))))
colnames(stereoscope[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[7]]@meta.data))))))))))))
colnames(stereoscope[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[8]]@meta.data))))))))))))

Spotiphy<-readRDS("/Path/data/Spotiphy_results.rds")
colnames(Spotiphy[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[1]]@meta.data))))))))))))
colnames(Spotiphy[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[2]]@meta.data))))))))))))
colnames(Spotiphy[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[3]]@meta.data))))))))))))
colnames(Spotiphy[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[4]]@meta.data))))))))))))
colnames(Spotiphy[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[5]]@meta.data))))))))))))
colnames(Spotiphy[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[6]]@meta.data))))))))))))
colnames(Spotiphy[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[7]]@meta.data))))))))))))
colnames(Spotiphy[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy[[8]]@meta.data))))))))))))

p1<-SpatialPlot(UCASpatial[[3]],features = "GZMK",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FEFAF8',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#FEFAF8")+
  labs(fill ='GZMK')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p2<-SpatialPlot(UCASpatial[[3]],features = "CXCL13",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#F9FCFE',brewer.pal(n = 9, name = "Blues")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#F9FCFE")+
  labs(fill ='CXCL13')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(UCASpatial[[3]],features = "XCL1",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#F9F8FB',brewer.pal(n = 9, name = "Purples")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#F9F8FB")+
  labs(fill ='XCL1')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(UCASpatial[[3]],features = "FAP",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FBF2DE',brewer.pal(n = 9, name = "#E2BC3C")))(n=100),
                       na.value = "#FBF2DE")+
  labs(fill ='FAP')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4
p1+p2+p3+p4


p_list <- list(p1,p2,p3,p4)
marker_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


DotPlot(sc_new,features = c("GZMK","CXCL13","KLRD1"))



p1<-SpatialPlot(UCASpatial[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(UCASpatial[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(UCASpatial[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(UCASpatial[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FBF2DE',brewer.pal(n = 9, name = "#E2BC3C")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p_list <- list(p1,p2,p3,p4)
UCASpatial_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')




p1<-SpatialPlot(CARD[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p1
p2<-SpatialPlot(CARD[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(CARD[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(CARD[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))


p_list <- list(p1,p2,p3,p4)
CARD_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')



p1<-SpatialPlot(SPOTlight[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(SPOTlight[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(SPOTlight[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(SPOTlight[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))




p_list <- list(p1,p2,p3,p4)
SPOTlight_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')



p1<-SpatialPlot(C2L[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(C2L[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(C2L[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(C2L[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p1+p2+p3


p_list <- list(p1,p2,p3,p4)
C2L_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialPlot(stereoscope[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(stereoscope[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(stereoscope[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(stereoscope[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p_list <- list(p1,p2,p3,p4)
stereoscope_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialPlot(RCTD[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(RCTD[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(RCTD[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(RCTD[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p1<-SpatialPlot(Spotiphy[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 8,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+theme(aspect.ratio = 0.7)+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(Spotiphy[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 8,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+theme(aspect.ratio = 0.7)+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(Spotiphy[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 8,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+theme(aspect.ratio = 0.7)+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(Spotiphy[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 8,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+theme(aspect.ratio = 0.7)+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))


p_list <- list(p1,p2,p3,p4)
Spotiphy_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
Spotiphy_p


plot_grid(marker_p,UCASpatial_p,CARD_p,C2L_p,RCTD_p,SPOTlight_p,stereoscope_p,
          align = 'hv',nrow = 1)


p1<-SpatialPlot(UCASpatial[[3]],features = "CD3G",image.alpha = 0,stroke = 0.01,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "white")+
  labs(fill ='CD3G')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p2<-SpatialPlot(UCASpatial[[3]],features = "CD8A",image.alpha = 0,stroke = 0.01,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "white")+
  labs(fill ='CD8A')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p_list <- list(p1,p2)
plot_grid(plotlist = p_list,ncol = 1,align = 'hv')

#Fig3e
source("/Path/summ_dist_V4_statistic.R")
CRC3a_st<-readRDS("/Path/data/CRC3a_distance_Seurat.RDS")
coord <- CRC3a_st@images$sliceD1@coordinates
stromal_name_list <- c("CD8.Tem","CD8.Tex","CD8.Trm") %>% as.list()
y="CD8.Tem"

dist_stromal <- lapply(stromal_name_list, function(y){
  weight.matr <- dplyr::select(CRC3a_st@meta.data,y,Epithelial_sum)
  dist <- summ_dist_V2(object = CRC3a_st,
                       barcode_set.2 = colnames(CRC3a_st)[which(CRC3a_st$Region == "Clone1")],
                       barcode_set.1 = colnames(CRC3a_st)[which(CRC3a_st@meta.data[[y]] > 0.03 )],
                       coord,weight = 'double',weight.matr = weight.matr)
  dist
})
dist_stromal



Tem<-data.frame(celltype="CD8+ Tem",dis=dist_stromal[[1]][[1]] * dist_stromal[[1]][[2]])
Tex<-data.frame(celltype="CD8+ Tex",dis=dist_stromal[[2]][[1]] * dist_stromal[[2]][[2]])
Trm<-data.frame(celltype="CD8+ Trm",dis=dist_stromal[[3]][[1]] * dist_stromal[[3]][[2]])
T_distribution<-rbind(Tem,Tex,Trm)
my_comparisons <- list( c("CD8+ Tem", "CD8+ Tex"), c("CD8+ Tex", "CD8+ Trm"), c("CD8+ Tem", "CD8+ Trm") )
p3<-ggplot(T_distribution[T_distribution$dis >= 0,],aes(celltype,dis,color=celltype))+geom_boxplot(outliers = FALSE)+
  scale_color_manual(values = c("#e6513a","#4a91c1","#7066a7"))+labs(x="Subpopulations",y="Weighted Distance to tumor")+
  guides(color=F)
p3
plot_data<-CRC3a_st@meta.data[CRC3a_st@meta.data$Region== "Clone1",c("CD8.Tem","CD8.Tex","CD8.Trm")]

plot_data<-reshape2::melt(plot_data)
p4<-ggplot(plot_data,aes(x=variable,y=value,fill=variable))+stat_summary(fun.y="mean",geom="bar",size=3) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.05)+scale_fill_manual(values = c("#e6513a","#4a91c1","#7066a7"))+guides(fill=F)+
  labs(x="Celltype",y="Proportion (%)")+coord_cartesian(ylim = c(0,0.02))+stat_compare_means(ref.group = "CD8.Tem",method = "t.test",label = "p.signif",label.y = 0.02)
p1+p2+p3+p4

#Fig3g
#distance calculation
theme_set(theme_cowplot())
dis1<-read.csv("/Path/data/distance_data1.txt",sep="\t")
dis2<-read.csv("/Path/data/distance_data2.txt",sep="\t")
dis3<-read.csv("/Path/data/distance_data3.txt",sep="\t")
dis4<-read.csv("/Path/data/distance_data4.txt",sep="\t")
dis5<-read.csv("/Path/data/distance_data5.txt",sep="\t")
dis6<-read.csv("/Path/data/distance_data6.txt",sep="\t")
dis7<-read.csv("/Path/data/distance_data7.txt",sep="\t")
dis8<-read.csv("/Path/data/distance_data8.txt",sep="\t")
dis9<-read.csv("/Path/data/distance_data9.txt",sep="\t")
dis10<-read.csv("/Path/data/distance_data10.txt",sep="\t")
distance_combine<-rbind(dis1,dis2,dis3,dis4,dis5,dis6,dis7,dis8,dis9,dis10)
distance_combine$Signed.distance.to.annotation.Tumor.µm[distance_combine$Signed.distance.to.annotation.Tumor.µm<0] = 0
CD8_name <- names(table(distance_combine$Classification))
CD8_name<-CD8_name[grep("CD8",CD8_name)]
distance_combine$Celltype<-"Others"
distance_combine$Celltype[distance_combine$Classification %in% c("CD8: CD103: PD-1","CD8: PD-1","GZMK: CD8: CD103: PD-1","GZMK: CD8: PD-1")]<-"CD8 Tex"
distance_combine$Celltype[distance_combine$Classification %in% c("CD8: CD103","GZMK: CD8: CD103")]<-"CD8 Trm"
distance_combine$Celltype[distance_combine$Classification %in% c("GZMK: CD8")]<-"CD8 Tem"
my_comparisons <- list( c("CD8 Tem", "CD8 Tex"), c("CD8 Tex", "CD8 Trm"), c("CD8 Tem", "CD8 Trm") )
p1<-ggplot(distance_combine[distance_combine$Celltype != "Others",],
       aes(x=Celltype,y=Signed.distance.to.annotation.Tumor.µm,color=Celltype))+geom_boxplot(outliers = F)+
  labs(x="Subpopulation",y="Distance to tumor (um)",title="mIHC(N=10)")+scale_color_manual(values = c("#cf4e39","#4888b5","#675e97"))+
  guides(color=F)+theme(plot.title = element_text(hjust=0.5))
p1
table(distance_combine$Celltype,distance_combine$Signed.distance.to.annotation.Tumor.µm == 0)
proportion<-data.frame(cell=c("CD8 Tem", "CD8 Tex","CD8 Trm"),value=c(122/37216 * 100,294/37216 * 100,396/37216 * 100))
p2<-ggplot(proportion,aes(x=cell,y=value,fill=cell))+geom_bar(stat = "identity",position = "dodge")+
  labs(x="Subpopulation",y="Proportion (%)",title="mIHC(N=10)")+scale_fill_manual(values = c("#cf4e39","#4888b5","#675e97"))+
  guides(fill=F)+theme(plot.title = element_text(hjust=0.5))+coord_cartesian(ylim = c(0,1.5))
p2
p1+p2

#Fig3i
# 4. Celltype infiltration and clonalTME phenotyping
library(dplyr)

CRC1_st<-readRDS("/Path/data/CRC1_st_clone.RDS")
CRC2_st<-readRDS("/Path/data/CRC2_st_clone.RDS")
CRC3a_st<-readRDS("/Path/data/CRC3a_st_clone.RDS")
CRC3b_st<-readRDS("/Path/data/CRC3b_st_clone.RDS")
CRC4_st<-readRDS("/Path/data/CRC4_st_clone.RDS")
CRC5_st<-readRDS("/Path/data/CRC5_st_clone.RDS")
CRC6_st<-readRDS("/Path/data/CRC6_st_clone.RDS")
CRC7_st<-readRDS("/Path/data/CRC7_st_clone.RDS")

Celltype_num=34

inhouse_cohort<-rbind(CRC1_st@meta.data[,c(4:37,40)],CRC2_st@meta.data[,c(4:37,40)],CRC3a_st@meta.data[,c(4:37,40)],CRC3b_st@meta.data[,c(4:37,40)])
public_cohort<-rbind(CRC4_st@meta.data[,c(5:38,41)],CRC5_st@meta.data[,c(4:37,40)],CRC6_st@meta.data[,c(4:37,40)],CRC7_st@meta.data[,c(4:37,40)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(CRC1_st@meta.data)),rep("CRC2",nrow(CRC2_st@meta.data)),rep("CRC3a",nrow(CRC3a_st@meta.data)),rep("CRC3b",nrow(CRC3b_st@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(CRC4_st@meta.data)),rep("CRC5",nrow(CRC5_st@meta.data)),rep("CRC6",nrow(CRC6_st@meta.data)),rep("CRC7",nrow(CRC7_st@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)

total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]


#total_sample_tumor_prop[,1:Celltype_num]


# cutoff: proportion > 0.02
#table(total_sample_tumor_prop$slice %in% c("CRC1","CRC2","CRC3a","CRC3b")
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)
#saveRDS(tumor_prop_avg,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/clonalTME_prop_avg.RDS")

#clonalTME phenotyping
colnames(tumor_prop_avg)
#tumor_prop_avg<-as.data.frame(tumor_prop_avg[,2:ncol(tumor_prop_avg)])
set.seed(1024)

km.res<-kmeans(scale(tumor_prop_avg[,c(2:6,28:30)]),centers = 2,nstart = 1024)
library(factoextra)
library(cowplot)
theme_set(theme_cowplot())
fviz_cluster(km.res, data = scale(tumor_prop_avg[,c(2:6,28:30)]),
             frame.alpha = 0, frame.level = 0.7, geom = "text",main = "clonalTME Clustering")+theme_cowplot()


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))

colnames(tumor_prop_avg)
my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)


tumor_prop_avg[,c(36,37)]
tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"
levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)]),cluster_cols = F,cluster_rows = F,border_color = NA,
                   gaps_col = c(3,8),gaps_row = c(3,8,10,20,22,25,29),
                   scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                   annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8)
ggsave("/data/huangzr/Spatial/Kras/version5/Fig_new/subTME_totalcell_heatmap.pdf",width = 6,height = 15)

#celltype difference between each TME phenotype
library(tidyr)
library(dplyr)

tumor_prop_avg_df <- as.data.frame(tumor_prop_avg)
tumor_prop_avg_df$Sample <- rownames(tumor_prop_avg_df)
colnames(tumor_prop_avg)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(tumor_prop_avg))))))))))))
colnames(total_sample_tumor_prop)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(total_sample_tumor_prop))))))))))))


celltype_list<-colnames(tumor_prop_avg)[c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)]

#calculate Roe
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"
library(dplyr)
total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()

Calculate_Roe <- function(avg_df)
{
  Ro_e_ex <- avg_df %>% dplyr::mutate(RowSum = rowSums(avg_df))
  NRow <- Ro_e_ex$RowSum
  Ro_e_ex <- Ro_e_ex%>% t() %>% as.data.frame() %>% dplyr::mutate(ColSum = colSums(Ro_e_ex))
  NCol <- Ro_e_ex$ColSum
  Roe_df <- avg_df
  N <- Ro_e_ex['RowSum','ColSum']
  for (i in rownames(avg_df)) {
    for (j in colnames(avg_df)) {
      Roe_df[i,j] = N*Roe_df[i,j] / (NRow[i] * NCol[j])
    }
  }
  return(Roe_df)
}


Roe_v_inhouse <- Calculate_Roe(Roe_total)
celltype_list[1:29]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>%
  select(celltype_list[1:29]) %>% t() %>% as.data.frame()
# Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>%
#   select(names(sort(kmeans_Roe_v_inhouse[["cluster"]]))) %>% t() %>% as.data.frame()

#colnames(Roe_v_inhouse) <- c('Dessert','Inflamed','Suppression')
max(Roe_v_inhouse)
min(Roe_v_inhouse)
library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-0.8,0.8))+labs(x="Ro/e",y="Subset")

ggsave("/Path/data/Celltype_enrich_subTME.pdf",height = 10,width = 5)

pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                   gaps_col = c(3,8),gaps_row = c(13,26),
                   scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                   annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8)
#Fig3j

# neighbour-based cellular community
st_neigbour<-readRDS("/Path/data/st.list.neighbor.rds")
length(st_neigbour)
CRC1_prop<-as.data.frame(as.matrix(t(st_neigbour[[1]])))
CRC1_prop$barcode<-paste("CRC1",rownames(CRC1_prop),sep = "_")
CRC2_prop<-as.data.frame(as.matrix(t(st_neigbour[[2]])))
CRC2_prop$barcode<-paste("CRC2",rownames(CRC2_prop),sep = "_")
CRC3a_prop<-as.data.frame(as.matrix(t(st_neigbour[[3]])))
CRC3a_prop$barcode<-paste("CRC3a",rownames(CRC3a_prop),sep = "_")
CRC3b_prop<-as.data.frame(as.matrix(t(st_neigbour[[4]])))
CRC3b_prop$barcode<-paste("CRC3b",rownames(CRC3b_prop),sep = "_")
CRC4_prop<-as.data.frame(as.matrix(t(st_neigbour[[5]])))
CRC4_prop$barcode<-paste("CRC4",rownames(CRC4_prop),sep = "_")
CRC5_prop<-as.data.frame(as.matrix(t(st_neigbour[[6]])))
CRC5_prop$barcode<-paste("CRC5",rownames(CRC5_prop),sep = "_")
CRC6_prop<-as.data.frame(as.matrix(t(st_neigbour[[7]])))
CRC6_prop$barcode<-paste("CRC6",rownames(CRC6_prop),sep = "_")
CRC7_prop<-as.data.frame(as.matrix(t(st_neigbour[[8]])))
CRC7_prop$barcode<-paste("CRC7",rownames(CRC7_prop),sep = "_")

neighbour_prop<-rbind(CRC1_prop,CRC2_prop,CRC3a_prop,CRC3b_prop,CRC4_prop,CRC5_prop,CRC6_prop,CRC7_prop)
colnames(neighbour_prop)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(neighbour_prop))))))))))))
colnames(neighbour_prop)<-gsub("Fibroblast","CD34+ iFib",colnames(neighbour_prop))
colnames(neighbour_prop)[35]<-"section_barcode"

CRC1_st<-readRDS("/Path/data/CRC1_st_clone.RDS")
CRC2_st<-readRDS("/Path/data/CRC2_st_clone.RDS")
CRC3a_st<-readRDS("/Path/data/CRC3a_st_clone.RDS")
CRC3b_st<-readRDS("/Path/data/CRC3b_st_clone.RDS")
CRC4_st<-readRDS("/Path/data/CRC4_st_clone.RDS")
CRC5_st<-readRDS("/Path/data/CRC5_st_clone.RDS")
CRC6_st<-readRDS("/Path/data/CRC6_st_clone.RDS")
CRC7_st<-readRDS("/Path/data/CRC7_st_clone.RDS")

Celltype_num=34

inhouse_cohort<-rbind(CRC1_st@meta.data[,c(4:37,40)],CRC2_st@meta.data[,c(4:37,40)],CRC3a_st@meta.data[,c(4:37,40)],CRC3b_st@meta.data[,c(4:37,40)])
public_cohort<-rbind(CRC4_st@meta.data[,c(5:38,41)],CRC5_st@meta.data[,c(4:37,40)],CRC6_st@meta.data[,c(4:37,40)],CRC7_st@meta.data[,c(4:37,40)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(CRC1_st@meta.data)),rep("CRC2",nrow(CRC2_st@meta.data)),rep("CRC3a",nrow(CRC3a_st@meta.data)),rep("CRC3b",nrow(CRC3b_st@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(CRC4_st@meta.data)),rep("CRC5",nrow(CRC5_st@meta.data)),rep("CRC6",nrow(CRC6_st@meta.data)),rep("CRC7",nrow(CRC7_st@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)

total_sample_porp<-rbind(public_cohort,inhouse_cohort)


total_sample_region_info<-total_sample_porp[,c("sample","slice","Region")]
total_sample_region_info$section_barcode <- paste(total_sample_region_info$slice,rownames(total_sample_region_info),sep = "_")

neighbour_prop$section_barcode<-gsub("-","_",substr(neighbour_prop$section_barcode, 1, 23))
total_sample_region_info$section_barcode<-gsub("-","_",substr(total_sample_region_info$section_barcode, 1, 23))

neighbour_prop$slice<-sub('^([^_]+)_.*', '\\1', neighbour_prop$section_barcode)
table(total_sample_region_info$section_barcode %in% neighbour_prop$section_barcode)
table(neighbour_prop$section_barcode %in% total_sample_region_info$section_barcode )

neighbour_prop<-merge(neighbour_prop,total_sample_region_info,by="section_barcode")

table(neighbour_prop$Region)
neighbour_prop_tumor<-neighbour_prop[!(neighbour_prop$Region %in% c("Remaining","Normal_Epithelial")),]
colnames(neighbour_prop_tumor)
tumor_sample_prop_immune_neighbour<-neighbour_prop_tumor[,2:35]
tumor_sample_prop_immune_neighbour<-tumor_sample_prop_immune_neighbour[,!(colnames(tumor_sample_prop_immune_neighbour) %in% c("Epithelial c1","Epithelial c2","Epithelial c3","Epithelial c4","Epithelial c5"))]

PCC_matrix <- matrix(NA, 29, 29)
for(i in 1:29) {
  for(j in 1:29) {
    temp <- cor.test(tumor_sample_prop_immune_neighbour[,i], tumor_sample_prop_immune_neighbour[,j],method = "pearson")
    PCC_matrix[i,j] <- temp$estimate
  }
}

rownames(PCC_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])
colnames(PCC_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])

#PCC_matrix[PCC_matrix<0.1]<-0

sig_matrix <- matrix(NA, 29, 29)
for(i in 1:29) {
  for(j in 1:29) {
    temp <- cor.test(tumor_sample_prop_immune_neighbour[,i], tumor_sample_prop_immune_neighbour[,j],method = "pearson")
    sig_matrix[i,j] <- temp$p.value
  }
}

rownames(sig_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])
colnames(sig_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])


sig_matrix

convert_p_to_stars <- function(matrix_of_p_values) {
  stars <- matrix("", nrow = nrow(matrix_of_p_values), ncol = ncol(matrix_of_p_values))

  for (i in 1:nrow(matrix_of_p_values)) {
    for (j in 1:ncol(matrix_of_p_values)) {
      p <- matrix_of_p_values[i, j]
      if (!is.na(p)) {
        if (p < 0.0001) {
          stars[i, j] <- "****"
        } else if (p < 0.001) {
          stars[i, j] <- "***"
        } else if (p < 0.01) {
          stars[i, j] <- "**"
        } else if (p < 0.05) {
          stars[i, j] <- "*"
        } else {
          stars[i, j] <- ""
        }
      } else {
        stars[i, j] <- ""
      }
    }
  }

  return(stars)
}

# 示例用法
# 假设 matrix_of_p_values 是你的 p 值矩阵
stars_matrix <- convert_p_to_stars(sig_matrix)

stars_matrix[PCC_matrix<0.2]<-""

#cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value==1),]
#cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value <= 0.2),]

cell_cor_long_matrix<-reshape2::melt(PCC_matrix)
cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value==1),]
cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value <= 0.3),]

write.csv(cell_cor_long_matrix,file = "/Path/data/cor_matrix_0.3.csv")

pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#7FB1D4","white","white","white","#FEBA78","#FE7F02"))(100),
                   cellwidth = 10,cellheight = 10,fontsize = 8,breaks = seq(-0.45, 0.45, length.out = 100),
                   border_color = "black",cutree_cols = 5,cutree_rows = 5,
                   clustering_method = "ward.D2",display_numbers = stars_matrix)

#Fig3k
#community ClonalTIME association
library(dplyr)

CRC1_st<-readRDS("/Path/data/CRC1_st_clone.RDS")
CRC2_st<-readRDS("/Path/data/CRC2_st_clone.RDS")
CRC3a_st<-readRDS("/Path/data/CRC3a_st_clone.RDS")
CRC3b_st<-readRDS("/Path/data/CRC3b_st_clone.RDS")
CRC4_st<-readRDS("/Path/data/CRC4_st_clone.RDS")
CRC5_st<-readRDS("/Path/data/CRC5_st_clone.RDS")
CRC6_st<-readRDS("/Path/data/CRC6_st_clone.RDS")
CRC7_st<-readRDS("/Path/data/CRC7_st_clone.RDS")
Celltype_num=34

inhouse_cohort<-rbind(CRC1_st@meta.data[,c(4:37,40)],CRC2_st@meta.data[,c(4:37,40)],CRC3a_st@meta.data[,c(4:37,40)],CRC3b_st@meta.data[,c(4:37,40)])
public_cohort<-rbind(CRC4_st@meta.data[,c(5:38,41)],CRC5_st@meta.data[,c(4:37,40)],CRC6_st@meta.data[,c(4:37,40)],CRC7_st@meta.data[,c(4:37,40)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(CRC1_st@meta.data)),rep("CRC2",nrow(CRC2_st@meta.data)),rep("CRC3a",nrow(CRC3a_st@meta.data)),rep("CRC3b",nrow(CRC3b_st@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(CRC4_st@meta.data)),rep("CRC5",nrow(CRC5_st@meta.data)),rep("CRC6",nrow(CRC6_st@meta.data)),rep("CRC7",nrow(CRC7_st@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)

total_sample_porp<-rbind(public_cohort,inhouse_cohort)
#total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

total_sample_porp[,1:Celltype_num] <- apply(total_sample_porp[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))
colnames(total_sample_porp)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(total_sample_porp))))))))))))
colnames(total_sample_porp)<-gsub("Fibroblast","CD34+ iFib",colnames(total_sample_porp))

#total_sample_porp_immune<-total_sample_porp[,c(1:34,35:37,81)]
total_sample_porp_immune<-total_sample_porp[,!(colnames(total_sample_porp) %in% c("Epithelial c1","Epithelial c2","Epithelial c3","Epithelial c4","Epithelial c5"))]
total_sample_porp_immune$total<-rowSums(total_sample_porp_immune[,1:29])

mat_normalized <- t(apply(total_sample_porp_immune[,1:29], MARGIN = 1, FUN = function(row) row / sum(row)))

print(mat_normalized)
colnames(mat_normalized)<-paste(colnames(mat_normalized),"normalized",sep = "_")
upper_quantile <- apply(mat_normalized, 2, quantile, probs = 0.995,na.rm=T)
lower_quantile <- apply(mat_normalized, 2, quantile, probs = 0.005,na.rm=T)
scaled_data <- scale(mat_normalized, center = lower_quantile, scale = upper_quantile - lower_quantile)
scaled_data<-as.data.frame(scaled_data)
colnames(scaled_data)<-gsub("normalized","scaled",colnames(scaled_data))
total_sample_porp_immune<-cbind(total_sample_porp_immune,mat_normalized,scaled_data)

total_sample_porp_immune$community1<-rowMeans(total_sample_porp_immune[,c("C1QC+ Mac_scaled","CD8+ Tem_scaled","CD34+ iFib_scaled","CD8+ Tex_scaled")])
total_sample_porp_immune$community2<-rowMeans(total_sample_porp_immune[,c("FN1+ Mac_scaled","Glial_scaled","Endothelial_scaled","PDPN+ Mac_scaled","FAP+ Myofib_scaled","Pericyte_scaled","Smooth Muscle_scaled")])
total_sample_porp_immune$community3<-rowMeans(total_sample_porp_immune[,c("Monocyte_scaled","Neutrophil_scaled")])


total_sample_porp_immune2<-total_sample_porp_immune[,c(32,92:94)]
total_sample_porp_immune2<-na.omit(total_sample_porp_immune2)
tumor_prop_avg<-aggregate(total_sample_porp_immune2[,2:4],by = list(total_sample_porp_immune2[,1]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)
tumor_prop_avg<-tumor_prop_avg[grep("Clone",rownames(tumor_prop_avg)),]
tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
library(ggpubr)
theme_set(theme_cowplot())
p1<-ggplot(tumor_prop_avg,aes(x=subTME_category,y=community1,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5,method = "t.test")+
  labs(x="ClonalTIME",y="Community Score",title = "Community1")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 0.13))
p2<-ggplot(tumor_prop_avg,aes(x=subTME_category,y=community2,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5,method = "t.test")+
  labs(x="ClonalTIME",y="Community Score",title = "Community2")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 0.18))
p3<-ggplot(tumor_prop_avg,aes(x=subTME_category,y=community3,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5,method = "t.test")+
  labs(x="ClonalTIME",y="Community Score",title = "Community3")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 0.16))

library("patchwork")
p1+p2+p3+plot_layout(ncol=2)

saveRDS(total_sample_porp_immune2,file="/Path/data/total_sample_porp_immune2.RDS")

#Fig3l
# PCC between each community

total_sample_porp_immune2<-readRDS("/Path/data/total_sample_porp_immune2.RDS")
PCC_matrix <- matrix(NA, 3, 3)
for(i in 1:3) {
  for(j in 1:3) {
    temp <- cor.test(total_sample_porp_immune2[,2:4][,i], total_sample_porp_immune2[,2:4][,j],method = "pearson")
    PCC_matrix[i,j] <- temp$estimate
  }
}

rownames(PCC_matrix) = colnames(total_sample_porp_immune2[,2:4])
colnames(PCC_matrix) = colnames(total_sample_porp_immune2[,2:4])

sig_matrix <- matrix(NA, 3, 3)
for(i in 1:3) {
  for(j in 1:3) {
    temp <- cor.test(total_sample_porp_immune2[,2:4][,i], total_sample_porp_immune2[,2:4][,j],method = "pearson")
    sig_matrix[i,j] <- temp$p.value
  }
}

rownames(sig_matrix) = colnames(total_sample_porp_immune2[,2:4])
colnames(sig_matrix) = colnames(total_sample_porp_immune2[,2:4])

sig_matrix
stars_matrix <- convert_p_to_stars(sig_matrix)
stars_matrix[abs(PCC_matrix)<0.2]<-""
pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","white","#FED9B4","#FE7F02"))(100),
                   cellwidth = 15,cellheight = 15,fontsize = 8,breaks = seq(-0.6, 0.6, length.out = 100),
                   border_color = "black",display_numbers = stars_matrix,cluster_rows = F,cluster_cols = F)

#Fig3m
total_sample_porp_immune<-readRDS("/Path/data/total_sample_porp_immune.RDS")
total_sample_porp_immune$Barcode<-substr(rownames(total_sample_porp_immune), 1, 18)

CRC4_st@meta.data<-merge(CRC4_st@meta.data,total_sample_porp_immune[total_sample_porp_immune$slice == "CRC4",c("Barcode","community1","community2","community3")],by="Barcode")
rownames(CRC4_st@meta.data)<-CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
CRC4_st<-UpdateSeuratObject(CRC4_st)

plot_data<-CRC4_st@meta.data[,c("Region","community1","community2","community3")]
theme_set(theme_cowplot())
p1<-ggplot(plot_data[plot_data$Region != "Remaining",],aes(x=Region,y=log2(community1),fill=Region,color=Region))+geom_boxplot()+
  scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+scale_color_manual(values = c("#378B39","#1F6C9F","#B62D28"))+
  labs(x="",y="log2(Community Score)",title = "Community1")+theme(plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)+
  stat_compare_means(ref.group = "Clone1",method = "wilcox.test",label = "p.signif")
p2<-ggplot(plot_data[plot_data$Region != "Remaining",],aes(x=Region,y=log2(community2),fill=Region,color=Region))+geom_boxplot()+
  scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+scale_color_manual(values = c("#378B39","#1F6C9F","#B62D28"))+
  labs(x="",y="log2(Community Score)",title = "Community2")+theme(plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)+
  stat_compare_means(ref.group = "Clone1",method = "wilcox.test",label = "p.signif")
p3<-ggplot(plot_data[plot_data$Region != "Remaining",],aes(x=Region,y=log2(community3),fill=Region,color=Region))+geom_boxplot()+
  scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+scale_color_manual(values = c("#378B39","#1F6C9F","#B62D28"))+
  labs(x="",y="log2(Community Score)",title = "Community3")+theme(plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)+
  stat_compare_means(ref.group = "Clone1",method = "wilcox.test",label = "p.signif")
library("patchwork")
p1+p2+p3+plot_layout(ncol=2)

CRC4_st_tumor<-CRC4_st[,CRC4_st$Node != "Remaining"]
SpatialPlot(CRC4_st_tumor,group.by ="Node",stroke = F,pt.size=1.2,cols = c("#A7CC83","#9BBFD0","#ED7B6C"))

#Fig3n
total_sample_porp_immune<-readRDS("/Path/data/total_sample_porp_immune.RDS")
total_sample_porp_immune$Barcode<-substr(rownames(total_sample_porp_immune), 1, 18)

CRC4_st@meta.data<-merge(CRC4_st@meta.data,total_sample_porp_immune[total_sample_porp_immune$slice == "CRC4",c("Barcode","community1","community2","community3")],by="Barcode")
rownames(CRC4_st@meta.data)<-CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
CRC4_st<-UpdateSeuratObject(CRC4_st)

colnames(CRC4_st@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CRC4_st@meta.data))))))))))))
colnames(CRC4_st@meta.data)<-gsub("Fibroblast","CD34+ iFib",colnames(CRC4_st@meta.data))

CRC4_st@meta.data[,c("community1","community2","community3")][is.na(CRC4_st@meta.data[,c("community1","community2","community3")])]<-0


CRC4_st@meta.data[CRC4_st@meta.data$Region == "Remaining",]<-NA

p1<-SpatialPlot(CRC4_st,features = c("community1"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "community1",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p2<-SpatialPlot(CRC4_st,features = c("CD34+ iFib"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD34+ iFib",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p3<-SpatialPlot(CRC4_st,features = c("CD8+ Tem"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8+ Tem",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p4<-SpatialPlot(CRC4_st,features = c("CD8+ Tex"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8+ Tex",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p5<-SpatialPlot(CRC4_st,features = c("C1QC+ Mac"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "C1QC+ Mac",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p1+p5+p2+p4+p3+plot_layout(ncol=5)
