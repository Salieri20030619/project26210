suppressPackageStartupMessages(library(Seurat))
p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase318_original_method_reconstruction_20261005_v1"
x<-readRDS("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase316_same_paper_repository_recovery_20261004_v1/downloads/Example_Data_Human_CRC_ST_data.rds")
cat("ASSAYS",names(x@assays),"\nMETA\n");print(colnames(x@meta.data));cat("COMMANDS\n");print(x@commands);cat("MISC\n");str(x@misc,max.level=2)
write.table(x@meta.data,file.path(p,"tables/P7_STORED_METADATA.tsv"),sep="\t",quote=FALSE,col.names=NA)
