# CNV evaluation in CRC4
.libPaths(c("/data/huangzr/R_package/R-3.6.1/R_package","/data/software/R/R_package","/data/xy/R/R_package/R-3.6.1/R_package","/data/software/R/R_news/R-3.6.1/lib/R/library"))

library("devtools")
library("BiocManager")
library("infercnv")
library("tidyverse")
library("Seurat")
library("phylogram")
library("ape")
library("hdf5r")
library("SpatialInferCNV")
library(ape)
setwd("/data/huangzr/Spatial/Kras/version5/Region_devide_2")
# Create the inferCNV object from the 3 files, without the reference group and mitochondria
NC_CRC4_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./CRC4_count_matrix.tsv",
                                                 gene_order_file="./siCNV_GeneOrderFile.tsv",
                                                 annotations_file="./Region_CRC4.tsv",
                                                 delim="\t",
                                                 ref_group_names="Remaining",
                                                 chr_exclude = c("chrM","chrX","chrY"))
saveRDS(NC_CRC4_infCNV,file="./CRC4_infCNV_clone_group_2.RDS")

#Run inferCNV
NC_CRC4_infCNV_res_spot = infercnv::run(NC_CRC4_infCNV,
                                        cutoff=0.1, #(see infercnv::run documentation)
                                        out_dir="./SPCNVoutputs_CRC4_clone_group_2",
                                        cluster_by_groups=TRUE, #unsupervised analysis
                                        HMM = TRUE,k_obs_groups = 3,
                                        denoise=TRUE) #denoising applies noise reduction for the plot
saveRDS(NC_CRC4_infCNV_res_spot,file="./CRC4_infCNV_clone_group_2_res_cluster_level_3.RDS")