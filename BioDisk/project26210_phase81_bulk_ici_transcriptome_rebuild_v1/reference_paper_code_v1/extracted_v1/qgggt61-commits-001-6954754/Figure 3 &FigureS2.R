#########Figure 3 &FigureS2
library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(batchelor)
library(Matrix.utils)
library(SPATA2)
library(SPATA2)
library(SPATAData)
library(tidyverse)
# devtools::install_github(repo = "kueckelj/confuns")
# devtools::install_github(repo = "theMILOlab/SPATAData")
# devtools::install_github(repo = "theMILOlab/SPATA2")

setwd('H:/2025年文章/stMTMC.RCC/Figure 2/rowdata/GSM5924033')
object <- 
  initiateSpataObjectVisium(
    sample_name = "UKF275T", 
    directory_visium = "H:/2025年文章/stMTMC.RCC/Figure 2/rowdata/GSM5924033/" # adjust to your liking 
  )
object <- identifySpatialOutliers(object, method = "obs")

plot_with_outliers <- plotSurface(object, color_by = "sp_outlier", clrp_adjust = c("TRUE" = "blue"))

# remove where sp_outlier == TRUE
object <- removeSpatialOutliers(object)

plot_without_outliers <- plotSurface(object, color_by = "sp_outlier")

# left plot
plot_with_outliers

# right plot
plot_without_outliers

# 4.1. 数据清理
# before
nGenes(object)
## [1] 33538

# removes stress genes
object <- removeGenesStress(object)

# removes genes that were not detected in any of the observations
object <- removeGenesZeroCounts(object)

# afterwards
nGenes(object)
## [1] 21556




# before
nObs(object)
## [1] 3213

# check for and remove observations with zero counts
object <- removeObsZeroCounts(object)

# afterwards 
nObs(object)
## [1] 3213

object <- computeMetaFeatures(object)

# plot left
plotSurface(object, color_by = "n_counts_gene")

# plot right
plotSurface(object, color_by = "n_distinct_gene")



# obtain matrix names prior to normalization
getMatrixNames(object)
## [1] "counts"

plot_before <- 
  plotSurface(object, color_by = "MAG") + labs(color = "MAG\n(Counts)")

# create log normalized matrix
object <- normalizeCounts(object, method = "LogNormalize")
## Normalizing layer: counts
## 02:03:11 Active matrix in assay 'gene': 'LogNormalize'

# obtain matrix names after normalization
getMatrixNames(object)
## [1] "counts"       "LogNormalize"

# check active matrix 
activeMatrix(object)
## [1] "LogNormalize"

plot_afterwards <- 
  plotSurface(object, color_by = "MAG") + labs(color = "MAG\n(logNorm)")

# left plot
plot_before

# right plot
plot_afterwards

# 5. 空间可变基因
# results are stored inside the SPATA2 object
object <- runSPARKX(object, verbose = FALSE)

# get genes with a p-value < 0.01
sparkx_genes <- getSparkxGenes(object, threshold_pval = 0.01)

str(sparkx_genes)
# visualize in space
plotSurfaceComparison(object, color_by = head(sparkx_genes, 6), nrow = 2)

plotSpatialAnnotations(object, ids = ids, nrow = 1)


plotImage(object_t313)

# uses the results of computeMetaFeatures(), count distribution is heavily right skewed -> transform
plotSurface(object, color_by = "n_counts_gene", transform_with = log10) + 
  labs(color = "counts\n(log10)")



# run the pipeline
object_t313 <- 
  runBayesSpaceClustering(
    object = object, 
    name = "bayes_space", # the name of the output grouping variable
    qs = 5
  )

# run PCA based on which clustering is conducted
object_t313 <- runPCA(object_t313, n_pcs = 20)

object_t313 <- 
  runKmeansClustering(
    object = object_t313, 
    ks = c(7, 8), 
    methods_kmeans = "Lloyd"
  )

# results are immediately stored in the objects feature data
getGroupingOptions(object_t313)


# left plot
plotSurface(
  object = object_t313, 
  color_by = "bayes_space", 
  pt_clrp = "egr"
)
ggsave('bayes_space.total.pdf',height = 5,width = 5)

# right plot
plotSurface(
  object = object_t313, 
  color_by = "Lloyd_k7",
  pt_clrp = "jco"
)

plotSurface(object = object_t269, pt_alpha = 0)


object_t269 <-
  runCNV(
    object = object_t313,
    # the directory must exist in your file system
    directory_cnv_folder = "H:/2025年文章/stMTMC.RCC/Figure 2/rowdata/GSM5924033//infercnv", 
    cnv_prefix = "Chr"
  )
plotCnvHeatmap(object = object_t269, across = "bayes_space", clrp = "npg")
ggsave('cnvheatmap.total.pdf',height = 5,width = 10)



# a more complex set up
plotCnvHeatmap(
  object = object_t269, 
  across = "bayes_space", 
  # across_subset = c("B1", "B5"), # dont show the transition part
  meta_vars = c("MET", "AURKA", "EGFR", "TPX2", "CDK6", "DNMT3B", "BRAF", "ZBP1"),  # visualize SNAP25 expression on the left
  # meta_vars_clrs='BuPu',
  # meta_vars_clrs = c("SNAP25" = "inferno",
  #                    "AURKA" = "inferno"
  #                    ), # with the inferno color spectrum 
  chrom_subset = c("1", '3',"6", "7", "12", "16", "20"), # only show these chromosomes
  ggpLayers = list(arm = list(legendNone())) # remove the chrom arm legend
)
ggsave('cnvheatmap.select.pdf',height = 5,width = 10)

# right plot
plotSurface(object_t269, color_by = "SETD2", pt_clrsp = "inferno")



plotSurface(
  object = object_t269, 
  color_by = "Chr20", 
  pt_clrsp = "BuPu"
)
ggsave('Chr20.pdf',height = 4,width = 4.5)

plotSurface(
  object = object_t269, 
  color_by = "Chr7", 
  pt_clrsp = "BuPu"
)
ggsave('Chr7.pdf',height = 4,width = 4.5)
plotSurface(
  object = object_t269, 
  color_by = "Chr3", 
  pt_clrsp = "BuPu"
)
ggsave('Chr3.pdf',height = 4,width = 4.5)


plotCnvLineplot(
  object = object_t269,
  across = "bayes_space", 
  n_bins_bcsp = 1000,
  n_bins_genes = 1000,
  nrow = 3
)
ggsave('bayes_space.line.pdf',height = 9,width = 12)



plotSurface( 
  object = object_t269, 
  color_by = "Chr5", 
  pt_clrsp = "Oslo" 
)
ggsave('Chr5.pdf',height = 4,width = 4.5)

plotSurface(object_t269, color_by = "MET", pt_clrsp = "inferno")
ggsave('MET.pdf',height = 4,width = 4.5)
plotSurface(object_t269, color_by = "AURKA", pt_clrsp = "inferno")
ggsave('AURKA.pdf',height = 4,width = 4.5)
plotSurface(object_t269, color_by = "EGFR", pt_clrsp = "inferno")
ggsave('EGFR.pdf',height = 4,width = 4.5)

plotSurface(object_t269, color_by = "TPX2", pt_clrsp = "inferno")
ggsave('TPX2.pdf',height = 4,width = 4.5)

plotSurface(object_t269, color_by = "CDK6", pt_clrsp = "inferno")
ggsave('CDK6.pdf',height = 4,width = 4.5)

plotSurface(object_t269, color_by = "DNMT3B", pt_clrsp = "inferno")
ggsave('DNMT3B.pdf',height = 4,width = 4.5)

plotSurface(object_t269, color_by = "BRAF", pt_clrsp = "inferno")
ggsave('BRAF.pdf',height = 4,width = 4.5)

plotSurface(object_t269, color_by = "ZBP1", pt_clrsp = "inferno")
ggsave('BRAF.pdf',height = 4,width = 4.5)

saveRDS(object_t269,'object_t269.rds')






object_t313 <- createImageAnnotations(object)

# results are stored in the processed object that can be downloaded 
plotSpatialAnnotations(
  object = object_t313,
  # ids = c("necrotic_area", "necrotic_edge", "necrotic_edge2"), 
  nrow = 1, 
  unit = "px", 
  fill = ggplot2::alpha("orange", 0.25)
)
ggsave('region.pdf',height = 4,width = 4.5)



# object_t313 <- loadExampleObject("UKF313T", process = TRUE, meta = TRUE)


necrotic_area <- 
  getSpatialAnnotation(
    object = object_t313,
    id = c("img_ann_1"), 
    add_image = T
  )

# print summary 
necrotic_area
slotNames(necrotic_area)
# visualize from the SPATA2 object, with expand = 0
plotSpatialAnnotations(object_t313, ids = "img_ann_1", expand = 0, fill = NA)

# expand = 0 is how the image is extracted when extracting the spatial annotation
# (in SPATA2, images are plotted upside down to fit into the cartesian coordinate system)
plot(EBImage::flip(necrotic_area@image))


# visualize from the SPATA2 object, with expand = 0
plotSpatialAnnotations(
  object = object_t313,
  ids = "img_ann_1", 
  expand = "500um", 
  fill = NA, 
  sb_dist = "500um", 
  sb_pos = c("4.5mm", "3mm")
)

necrotic_area_expanded <- 
  getSpatialAnnotation(
    object = object_t313,
    id = "img_ann_1", 
    add_image = T, 
    expand = "500um"
  )


necrotic_ids <- c("img_ann_1", "img_ann_2",'img_ann_3')

# obtain distance data relative to spatial annotations (SA) with getCoordsDfSA()
coords_df <- 
  getCoordsDfSA(object_t313, ids = necrotic_ids, unit = "mm", binwidth = "200um", core0 = T) 

# add distance to necrosis as meta feature to the SPATA2 object
# to make them accessible for other functions
object_t313 <- 
  addFeatures(object_t313, feature_df = coords_df[,c("barcodes", "dist", "bins_dist")], overwrite = T)

# create ggproto layer for further surface plots
necrotic_outline <- 
  ggpLayerSpatAnnOutline(object_t313, ids = necrotic_ids, incl_edge = T, fill = "grey")

sas_screning <- 
  ggpLayerScreeningDirectionSAS(object_t313, ids = necrotic_ids, line_size = 0.75)

# show plots 
plotImage(object_t313) + 
  necrotic_outline + 
  sas_screning

plotSurface(object_t313, color_by = "dist", pt_clrsp = "plasma") + 
  necrotic_outline + 
  sas_screning + 
  labs(color = "Dist [mm]")
ggsave('Dist.pdf',height = 4,width = 4.5)


# 定义所有 gene sets
gene_sets <- list(
  POSTNCAF = c("POSTN", "FAP", "PDGFRA", "PDGFRB", "ACTA2", "COL1A1", "COL3A1", "TAGLN", "VIM", "SPARC", "LOXL1", "MMP2", "MMP14"),
  Collagen = c("COL4A1", "COL4A2", "COL11A1", "COL11A2", "COL12A1"),
  LOXfamily = c("LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"),
  MMP = c("MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP11", "MMP13", "MMP14"),
  APOEMph = c("APOE", "CD68", "CD163", "MRC1", "C1QA", "C1QB", "C1QC", "MARCO", "TREM2", "LGALS3", "MSR1", "IL10", "S100A9", "CSF1R"),
  ImmunosuppressiveCytokines = c("TGFB1", "IL10", "IL4"),
  GrowthFactors = c("TGFA", "TGFB1", "TGFB2", "TGFB3", "PDGFA", "PDGFB", "EGFR", "FGF1", "FGF2"),
  ImmuneCheckpoints = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "TIGIT", "HAVCR2", "LAG3"),
  # MyeloidSignals = c("CD47", "SIRPA", "CD276", "IDO1", "ARG1", "CSF1", "CSF1R"),
  Costimulatory = c("TNFRSF14", "BTLA", "CD160", "VSIR"),
  CCL = c("CCL22", "CCL2", "CCL5", "CCL17", "CXCL2", "CXCL5", "CXCL1"),
  Chemokines = c("CCL22", "CCL2", "CXCL8", "CXCL12", "CCR4", "CXCR2")
  
)

# 假设你已经定义了 gene_sets 列表
for (name in names(gene_sets)) {
  object_t313 <- addGeneSet(
    object = object_t313,
    genes = gene_sets[[name]],
    name = name,       # e.g., "Collagen"
    class = "Custom",  # 避免下划线！不能包含 _
    overwrite = TRUE
  )
}

plotSurface(object_t313, color_by = "dist", pt_clrsp = "plasma") + 
  # necrotic_outline + 
  sas_screning + 
  labs(color = "Dist [mm]")
ggsave('Distance Dist.pdf',height = 4,width = 4.5)

plotSasHeatmap(
  object = object_t313,
  variables = paste0("Custom_", names(gene_sets)),  # 拼接完整变量名
  ids = necrotic_ids,
  clrsp = "BuPu"
) + labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance to Tumor.total5.pdf',height = 4,width = 4.5)
saveRDS(object_t313,'object_t313.rds')

plotSurface(object_t313, color_by = "Custom_POSTNCAF", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_POSTNCAF Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_POSTNCAF", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_POSTNCAF Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_POSTNCAF", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_POSTNCAF.pdf',height = 4,width = 6)


plotSurface(object_t313, color_by = "Custom_Collagen", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_Collagen Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_Collagen", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_Collagen Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_Collagen", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_Collagen.pdf',height = 4,width = 6)


plotSurface(object_t313, color_by = "Custom_LOXfamily", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_LOXfamily Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_LOXfamily", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_LOXfamily Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_LOXfamily", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_LOXfamily.pdf',height = 4,width = 6)



plotSurface(object_t313, color_by = "Custom_MMP", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_MMP object_t313 Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_MMP", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_MMP Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_MMP", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_MMP.pdf',height = 4,width = 6)



plotSurface(object_t313, color_by = "Custom_APOEMph", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_APOEMph Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_APOEMph", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_APOEMph Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_APOEMph", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_APOEMph.pdf',height = 4,width = 6)



plotSurface(object_t313, color_by = "Custom_ImmunosuppressiveCytokines", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_ImmunosuppressiveCytokines Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_ImmunosuppressiveCytokines", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_ImmunosuppressiveCytokines Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_ImmunosuppressiveCytokines", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_ImmunosuppressiveCytokines.pdf',height = 4,width = 6)


plotSurface(object_t313, color_by = "Custom_GrowthFactors", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_GrowthFactors Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_GrowthFactors", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_GrowthFactors Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_GrowthFactors", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_GrowthFactors.pdf',height = 4,width = 6)



plotSurface(object_t313, color_by = "Custom_ImmuneCheckpoints", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_ImmuneCheckpoints Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_ImmuneCheckpoints", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_ImmuneCheckpoints Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_ImmuneCheckpoints", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_ImmuneCheckpoints.pdf',height = 4,width = 6)


plotSurface(object_t313, color_by = "Custom_Costimulatory", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_Costimulatory Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_Costimulatory", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_Costimulatory Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_Costimulatory", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_Costimulatory.pdf',height = 4,width = 6)


plotSurface(object_t313, color_by = "Custom_CCL", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_CCL Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_CCL", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_CCL Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_CCL", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_CCL.pdf',height = 4,width = 6)


plotSurface(object_t313, color_by = "Custom_Chemokines", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_Chemokines Distance Dist.pdf',height = 4,width = 4.5)
plotSurface(object_t313, color_by = "Custom_Chemokines", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  # necrotic_outline + 
  sas_screning
ggsave('Custom_Chemokines Distance Dist.pdf',height = 4,width = 4.5)

plotSasLineplot(object_t313, ids = necrotic_ids, variables = "Custom_Chemokines", line_color = "blue") + 
  labs(x = "Distance to Tumor Borders [mm]")
ggsave('Distance Dist line Custom_Chemokines.pdf',height = 4,width = 6)

