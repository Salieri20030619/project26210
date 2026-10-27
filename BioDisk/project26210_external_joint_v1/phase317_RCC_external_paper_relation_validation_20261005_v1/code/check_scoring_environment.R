p<-"/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase317_RCC_external_paper_relation_validation_20261005_v1"
suppressPackageStartupMessages(library(Seurat))
writeLines(deparse(Seurat::AddModuleScore),file.path(p,"resources/INSTALLED_AddModuleScore.R"))
writeLines(capture.output(sessionInfo()),file.path(p,"resources/SCORING_R_SESSION.txt"))
cat("Seurat",as.character(packageVersion("Seurat")),"SeuratObject",as.character(packageVersion("SeuratObject")),"\n")
print(formals(Seurat::AddModuleScore))
