suppressPackageStartupMessages(library(clusterProfiler))
d <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase200_union_C4_C10_C15_reassessment_20260918_v1"
gmt <- read.gmt("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1/06_functional_programs/run_v1/resources/KEGG_2021_Human.gmt")
rr <- read.delim(file.path(d,"gsea/RANK_REGISTRY.tsv"))
all <- list()
for(model in rr$model) for(version in c("complete","without_SPP1_TYMP")) {
 a <- read.delim(file.path(d,paste0("gsea/",model,"_rank.tsv")))
 if(version=="without_SPP1_TYMP") a <- a[!a$gene %in% c("SPP1","TYMP"),]
 stats <- setNames(a$score,a$gene)
 set.seed(262200)
 fit <- GSEA(stats,TERM2GENE=gmt,minGSSize=10,maxGSSize=500,exponent=1,pvalueCutoff=1,pAdjustMethod="BH",by="fgsea",eps=0,nPermSimple=10000,seed=TRUE,verbose=FALSE)
 z <- as.data.frame(fit@result)
 sizes <- tapply(gmt$gene %in% names(stats),gmt$term,sum)
 eligible <- names(sizes)[sizes>=10 & sizes<=500]
 z$p.adjust <- p.adjust(z$pvalue,method="BH",n=length(eligible))
 fit@result$p.adjust <- z$p.adjust
 full <- merge(data.frame(ID=eligible,intersection_size=sizes[eligible]),z,by="ID",all.x=TRUE,sort=FALSE)
 full$model <- model;full$version <- version;full$patients <- rr$n[rr$model==model]
 write.table(full,file.path(d,paste0("gsea/",model,"__",version,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 saveRDS(fit,file.path(d,paste0("gsea/",model,"__",version,".rds")))
 all[[paste(model,version)]] <- full
 cat("DONE",model,version,"family",length(eligible),"missing",sum(is.na(full$pvalue)),"\n")
}
out <- do.call(rbind,all)
write.table(out,file.path(d,"gsea/ALL_GSEA_RESULTS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
focus <- out[out$ID %in% c("Cholesterol metabolism","Lysosome","Phagosome","Antigen processing and presentation","Cell cycle"),]
write.table(focus,file.path(d,"gsea/FOCUS_GSEA_RESULTS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(d,"qc/GSEA_SESSION_INFO.txt"))
writeLines("COMPLETE",file.path(d,"qc/GSEA_COMPLETE.txt"))
