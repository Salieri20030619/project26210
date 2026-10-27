options(stringsAsFactors=FALSE, bitmapType="cairo")
suppressPackageStartupMessages({
 library(clusterProfiler)
 library(enrichplot)
 library(ggplot2)
})
root <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase199_fig3_marker_GSEA_20260918_v5"
setwd(root)
ranktab <- read.delim("tables/RANKED_PATIENT_GENES.tsv", check.names=FALSE)
gmt_lines <- strsplit(readLines("inputs/KEGG_2021_Human.gmt"), "\t", fixed=FALSE)
sets <- lapply(gmt_lines, function(x) unique(x[-c(1,2)]))
names(sets) <- vapply(gmt_lines, function(x) x[1], character(1))
sets <- lapply(sets,function(x) x[nzchar(x)])
term2gene <- data.frame(term=rep(names(sets),lengths(sets)),gene=unlist(sets,use.names=FALSE))
term2gene <- unique(term2gene)
focus <- c("Cholesterol metabolism","Lysosome","Phagosome","Antigen processing and presentation")
stopifnot(all(focus %in% names(sets)),!anyDuplicated(ranktab$gene))
objects <- list()
for (mode in c("primary","exclude_TYMP_SPP1")) {
 a <- ranktab
 if (mode=="exclude_TYMP_SPP1") a <- a[!a$gene %in% c("TYMP","SPP1"),]
 stats <- setNames(a$ranking_score,a$gene)
 stopifnot(all(is.finite(stats)),all(diff(stats)<=0))
 set.seed(20260918)
 obj <- GSEA(geneList=stats, TERM2GENE=term2gene, minGSSize=10, maxGSSize=500,
             exponent=1, pvalueCutoff=1, pAdjustMethod="BH", verbose=FALSE,
             seed=TRUE, by="fgsea", eps=0, nPermSimple=10000)
 
 tab <- as.data.frame(obj@result)
 overlaps <- vapply(sets,function(x)sum(x %in% names(stats)),integer(1))
 eligible <- names(overlaps)[overlaps>=10 & overlaps<=500]
 tab$p.adjust <- p.adjust(tab$pvalue,method="BH",n=length(eligible))
 obj@result$p.adjust <- tab$p.adjust
 coverage <- data.frame(ID=eligible,setSize=overlaps[eligible])
 fulltab <- merge(coverage,tab[,setdiff(names(tab),"setSize"),drop=FALSE],by="ID",all.x=TRUE,sort=FALSE)
 write.table(fulltab,file=paste0("tables/GSEA_KEGG_FULL_FAMILY_",mode,".tsv"),sep="\t",row.names=FALSE,quote=FALSE)
 cat("ELIGIBLE_FAMILY",length(eligible),"NON_ESTIMABLE",sum(is.na(fulltab$pvalue)),"\n")

 write.table(tab,file=paste0("tables/GSEA_KEGG_",mode,".tsv"),sep="\t",row.names=FALSE,quote=FALSE)
 chosen <- tab[match(focus,tab$ID),]
 stopifnot(!anyNA(chosen$ID),all(is.finite(chosen$NES)),all(is.finite(chosen$p.adjust)))
 write.table(chosen,file=paste0("tables/GSEA_FOUR_PROGRAMS_",mode,".tsv"),sep="\t",row.names=FALSE,quote=FALSE)
 saveRDS(obj,paste0("tables/GSEA_",mode,".rds"))
 objects[[mode]] <- obj
 cat("GSEA_COMPLETE",mode,"genes",length(stats),"pathways",nrow(tab),"\n")
 print(as.data.frame(chosen[,c("ID","setSize","NES","pvalue","p.adjust")]))
}
obj <- objects$primary
obj@result$Description[obj@result$ID=="Antigen processing and presentation"] <- "Antigen presentation"
colors <- c("#B9852D","#A13B80","#368899","#4662A6")
p <- gseaplot2(obj,geneSetID=focus,title="C4 versus other mature macrophages",
              base_size=11,rel_heights=c(1.5,.4,1),subplots=1:3,
              pvalue_table=FALSE,ES_geom="line",color=colors)
pdf("panels/GSEA_curves_raw.pdf",width=13,height=4.8,useDingbats=FALSE)
print(p)
dev.off()
saveRDS(p,"panels/GSEA_plot.rds")
capture.output(sessionInfo(),file="qc/R_SESSION_INFO.txt")
writeLines(c("seed=20260918","rank=mean paired delta / (SD paired delta / sqrt(26))",
 "patients=26","minGSSize=10","maxGSSize=500","exponent=1","eps=0","nPermSimple=10000",
 "BH family=all eligible KEGG pathways per version","GSEA gene-set null; not patient-label permutation"),
 "qc/GSEA_PARAMETERS.txt")
cat("GSEA_FIGURE_COMPLETE\n")
