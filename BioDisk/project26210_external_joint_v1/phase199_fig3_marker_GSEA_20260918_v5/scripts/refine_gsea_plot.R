suppressPackageStartupMessages({library(enrichplot);library(ggplot2)})
d <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase199_fig3_marker_GSEA_20260918_v5"
obj <- readRDS(file.path(d,"tables/GSEA_primary.rds"))
obj@result$Description[obj@result$ID=="Antigen processing and presentation"] <- "Antigen presentation"
focus <- c("Cholesterol metabolism","Lysosome","Phagosome","Antigen processing and presentation")
colors <- c("Cholesterol metabolism"="#B9852D","Lysosome"="#A13B80","Phagosome"="#368899","Antigen presentation"="#4662A6")
p <- gseaplot2(obj,geneSetID=focus,title="C4 versus other mature macrophages",base_size=11,
 rel_heights=c(1.5,.4,1),subplots=1:3,pvalue_table=FALSE,ES_geom="line",color=colors)
pdf(file.path(d,"panels/GSEA_curves_raw.pdf"),width=13,height=4.8,useDingbats=FALSE)
print(p)
dev.off()
saveRDS(p,file.path(d,"panels/GSEA_plot.rds"))
write.table(data.frame(pathway=names(colors),color=unname(colors)),file.path(d,"tables/GSEA_COLOR_MAP.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("GSEA_COLORS_REFINED\n")
