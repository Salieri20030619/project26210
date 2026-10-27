#########Figure 4&FigureS3
processed_data <- map(integrated_data, basic_qc)

# Optionally save each processed object or continue with downstream analysis
for (i in seq_along(processed_data)) {
  saveRDS(processed_data[[i]], file = paste0("processed_data_", i, ".rds"))
}
rm(integrated_data)
# Check each Seurat object's project attribute
library(purrr)
library(Seurat)

# Ensure each Seurat object has a project name
processed_data <- map(processed_data, function(x) {
  if (is.null(Project(x)) || length(Project(x)) != 1) {
    Project(x) <- "DefaultProject"  # Set a default project name if missing
  }
  return(x)
})
def_assay <- "RNA"
hvg_list <- processed_data %>%
  map(~{
    DefaultAssay(.x) <- def_assay
    .x <- FindVariableFeatures(.x, selection.method = "vst", nfeatures = 3000)
    .x@assays[[def_assay]]@var.features
  }) %>% unlist() %>% table() %>% sort(decreasing = TRUE)


gene_selection_plt <- enframe(hvg_list) %>%
  group_by(value) %>%
  summarize(ngenes = length(name)) %>%
  ggplot(aes(x = value, y = ngenes)) + geom_bar(stat = "identity")

gene_selection <- names(hvg_list)[1:3000]
processed_data <- purrr::reduce(processed_data, merge, merge.data = TRUE)
cat("Merged data. Object size:", object.size(processed_data), "\n")
# Set default assay and preprocess data
DefaultAssay(processed_data) <- def_assay
processed_data <- processed_data %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = gene_selection, npcs = 30, verbose = FALSE)
# Run Harmony for batch correction
processed_data <- RunHarmony(processed_data, 
                             group.by.vars = "orig.ident", 
                             plot_convergence = TRUE, 
                             assay.use = def_assay, 
                             max.iter.harmony = 20)
subscRNA <- FindNeighbors(processed_data, reduction = "harmony", dims = 1:15)
subscRNA <- FindClusters(subscRNA, 
                         resolution = seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.2))
subscRNA2 <- RunUMAP(subscRNA, reduction = "harmony", dims = 1:15)
DimPlot(subscRNA2,raster = T)
DimPlot(subscRNA2,group.by = 'orig.ident',raster = T)
scRNA3<-RunUMAP(subscRNA,reduction="pca",dims=1:20,verbose=FALSE)
DimPlot(scRNA3,raster = T,label = T,group.by = 'RNA_snn_res.0.9')
DimPlot(scRNA3,group.by = 'orig.ident',raster = T)

library(SCP)
FeatureDimPlot(scRNA3,features = c('C1QB','C1QA'))
FeatureDimPlot(scRNA3,features = c("EPCAM", "KRT8", "KRT17", "KRT18", "KRT19"))
FeatureDimPlot(scRNA3,features = c('CD3D','CD3E'))
FeatureDimPlot(scRNA3,features = c('CD79A','CD79B'))
FeatureDimPlot(scRNA3,features = c("TPSAB1", "CPA3"))
FeatureDimPlot(scRNA3,features = c('LUM' , 'FGF7', 'MME','ACTA2'))
FeatureDimPlot(scRNA3,features = c("CD1C", "CD1E", "FCER1A"))
FeatureDimPlot(scRNA3,features = c("VWF", "PECAM1"))
FeatureDimPlot(scRNA3,features = c("NKG7", "KLRB1"))
FeatureDimPlot(scRNA3,features = c("CD4", "CD8A"))
FeatureDimPlot(scRNA3,features = c("S100A9", "S100A8"))
FeatureDimPlot(scRNA3,features = c("CD14", "CD16",'CD68','CD168'))
FeatureDimPlot(scRNA3,features = c("APOE", "SPP1",'CXCL9','FAP'))
CellDimPlot(scRNA3,group.by = 'RNA_snn_res.0.1')
DimPlot(scRNA3,label = T)
library(clustree)
clustree(scRNA3)
gene_list <- list(
  Mast = c("TPSAB1", "CPA3"),
  B = c("MS4A1", "CD79A", "CD79B"),
  NK = c("NKG7", "KLRB1"),
  Tcells = c("CD3D", "CD3E", "CD3G"),
  Monocytes = c("S100A9", "S100A8"),
  Dendritic = c("CD1C", "CD1E", "FCER1A"),
  Macrophage = c("CD14",'CD68', "C1QA", "C1QB"),
  Monocle = c('S100A9','S100A8'),
  Epithelial = c("EPCAM", "KRT8", "KRT17", "KRT18", "KRT19"),
  Endothelial = c("VWF", "PECAM1"),
  SMC = c("ACTA2", "MYH11"),
  Pericyte = c("ABCC9", "PDGFRB", "RGS5")
)
library(Seurat)
library(Matrix)
seurat_to_adata(scRNA3,Dimension='UMAP',path = './seurat_')
DotPlot(sce.all,features = c( "MS4A1", "BLNK", "CD19", "CD79A", "CD79B"),group.by = 'seurat_clusters')

table(sce.all$seurat_clusters)


sce <- subset(scRNA3, RNA_snn_res.0.9 %in% c('2','17','10','18','57','61'))


genes_to_check = c('APOE','CXCL9','SPP1')

p = DotPlot(sce, features = unique(genes_to_check),
            assay='RNA' ,group.by = 'RNA_snn_res.0.9' )  + coord_flip()+ RotatedAxis()
p

VlnPlot(sce, features = genes_to_check, pt.size = 0, layer = "scale.data")



table(scRNA3$celltype)

CellDimPlot3D(scRNA3,group.by = 'celltype',reduction = 'umap')
CellDensityPlot(scRNA3,features = 'APOE')
CellDimPlot(scRNA3,group.by = c('RNA_snn_res.0.9','celltype'))


ht1 <- GroupHeatmap(scRNA3,
                    features = genes_to_check,feature_split = 'celltype',
                    group.by = c( "celltype")
)
print(ht1$plot)



mac <- subset(scRNA3, celltype %in% c('Mph'))
# Convert to factor firs
mac@meta.data$celltype <- as.factor(mac@meta.data$celltype)

# Now you can use droplevels
mac@meta.data$celltype <- droplevels(mac@meta.data$celltype, 
                                     exclude = setdiff(levels(mac@meta.data$celltype),
                                                       unique(mac@meta.data$celltype)))

# Check the table
table(mac$celltype)


mac <- ScaleData(mac, vars.to.regress = c("nFeature_RNA", "percent_mito",'percent_ribo','percent_hb','cohort','orig.ident'), verbose = FALSE)
mac <- FindVariableFeatures(mac, nfeatures = 4000)
mac <- RunPCA(mac, npcs = 50, verbose = FALSE)
mac <- FindNeighbors(mac, reduction = "pca", dims = 1:50)
mac <- FindClusters(mac, 
                    resolution = seq(from = 0.1, 
                                     to = 1.0, 
                                     by = 0.2))
mac <- RunUMAP(mac, reduction = "pca", dims = 1:50)
mac <- RunTSNE(mac, reduction = "pca", dims = 1:50)
library(clustree)
clustree(mac)
Idents(mac) <- "RNA_snn_res.0.1"
mac$seurat_clusters <- mac@active.ident
DimPlot(mac, label = T,pt.size = 1)
DimPlot(mac, group.by = 'cohort',pt.size =0.0000000000000001)

## 重新harmony
library(harmony)
set.seed(1000)
mac2 <- RunHarmony(mac, group.by.vars = c("orig.ident",'cohort', 'Histology', 'Treatment', 'Respponse'))
DimPlot(mac2, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(mac2,reduction = 'harmony')
mac2 <- FindNeighbors(mac2, reduction = 'harmony',dims = 1:20)
mac2 <- FindClusters(mac2, 
                     resolution = seq(from = 0.1, 
                                      to = 1, 
                                      by = 0.1))
mac2 <- RunUMAP(mac2, reduction = 'harmony',dims = 1:20)
clustree(mac2)
Idents(mac2) <- "RNA_snn_res.0.1"
DimPlot(mac2, label = T,pt.size = 0.0001)

Epi <- FindVariableFeatures(Epi, nfeatures = 4000)
Epi <- RunPCA(Epi, npcs = 50, verbose = FALSE)
Epi <- FindNeighbors(Epi, reduction = "pca", dims = 1:50)
Epi <- FindClusters(Epi, 
                    resolution = seq(from = 0.1, 
                                     to = 1.0, 
                                     by = 0.2))
Epi <- RunUMAP(Epi, reduction = "pca", dims = 1:50)
Epi <- RunTSNE(Epi, reduction = "pca", dims = 1:50)
library(clustree)
clustree(Epi)
Idents(Epi) <- "RNA_snn_res.0.1"
Epi$seurat_clusters <- Epi@active.ident
DimPlot(Epi, label = T,pt.size = 1)
DimPlot(Epi, group.by = 'cohort',pt.size =0.0000000000000001)
library(harmony)
set.seed(1000)
Epi2 <- RunHarmony(Epi, group.by.vars = c("orig.ident"))
DimPlot(Epi2, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(Epi2,reduction = 'harmony')
Epi2 <- FindNeighbors(Epi2, reduction = 'harmony',dims = 1:20)
Epi2 <- FindClusters(Epi2, 
                     resolution = seq(from = 0.1, 
                                      to = 1, 
                                      by = 0.1))
Epi2 <- RunUMAP(Epi2, reduction = 'harmony',dims = 1:20)
clustree(Epi2)
Idents(Epi2) <- "RNA_snn_res.0.3"
DimPlot(Epi2, label = T)

DimPlot(Epi2, label = T,group.by = 'Histology')



library(infercnv)
library(Seurat)
library(SeuratObject)
sce_CNV <- Epi2
sce_CNV.tumor <- sce_CNV[,sce_CNV@meta.data$Histology %in% c("Tumor")]
sce_CNV.normal <- sce_CNV[,sce_CNV@meta.data$Histology %in% c("Normal")]
sce_CNV.tumor <- sce_CNV[,sce_CNV@meta.data$orig.ident %in% c( "GSM112169", "GSM6112170", "GSM7028034", "GSM7028035", "GSM7028036",
                                                               "GSM7028037", "GSM7028039")]
sce_CNV.normal <- sce_CNV[, sce_CNV@meta.data$orig.ident %in% c('RCC1.Normal',
                                                                'RCC2.Normal',
                                                                'RCC3.Normal',
                                                                'RCC4.Normal',
                                                                'RCC5.Normal'
)]

sce_CNV<-merge(sce_CNV.normal,sce_CNV.tumor)
table(sce_CNV$orig.ident)
saveRDS(as.matrix(sce_CNV[["RNA"]]@counts), "sce_CNV_counts.matrix")
matrix_counts <- readRDS("sce_CNV_counts.matrix")

sce_CNV$newcelltype <- paste0(sce_CNV$Histology,"_",sce_CNV$celltype)
unique(sce_CNV$newcelltype)


write.table(sce_CNV$newcelltype, "celltype.label.txt", sep = "\t", quote = F, col.names = F)
saveRDS(infercnv_obj,file = "./infercnv_obj_inputfor_run.rds")

options(scipen = 100)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="EE3",  
                             no_prelim_plot = T,
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             min_cells_per_gene = 10,
                             num_threads=5,
                             write_expr_matrix = T
)


