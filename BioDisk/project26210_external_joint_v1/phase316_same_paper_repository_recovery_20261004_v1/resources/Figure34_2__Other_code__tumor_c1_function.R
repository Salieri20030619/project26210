.libPaths("/data/xy/R/R_package/R4/R_package_4")

sc_marker<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/PCC_in_real/Seurat_cluster_markers.rds")
names(table(sc_marker$cluster))

CD8_Tem_marker<-sc_marker$gene[sc_marker$cluster == "CD8 Tem"]
CD8_Tex_marker<-sc_marker$gene[sc_marker$cluster == "CD8 Tex"]
C1QC_mac_marker<-sc_marker$gene[sc_marker$cluster == "IGF1+ Mac"]
Fibroblast_marker<-sc_marker$gene[sc_marker$cluster == "Fibroblast"]

CD8_Tem_GO<-enrichGO(CD8_Tem_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
CD8_Tex_GO<-enrichGO(CD8_Tex_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
C1QC_mac_GO<-enrichGO(C1QC_mac_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Fibroblast_GO<-enrichGO(Fibroblast_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")


CD8_Tem_pathway<-data.frame(cell="CD8 Tem",pathway=CD8_Tem_GO@result$Description[CD8_Tem_GO@result$p.adjust < 0.001],pathwayID=CD8_Tem_GO@result$ID[CD8_Tem_GO@result$p.adjust < 0.001])
CD8_Tex_pathway<-data.frame(cell="CD8 Tex",pathway=CD8_Tex_GO@result$Description[CD8_Tex_GO@result$p.adjust < 0.001],pathwayID=CD8_Tex_GO@result$ID[CD8_Tex_GO@result$p.adjust < 0.001])
C1QC_mac_pathway<-data.frame(cell="C1QC mac",pathway=C1QC_mac_GO@result$Description[C1QC_mac_GO@result$p.adjust < 0.001],pathwayID=C1QC_mac_GO@result$ID[C1QC_mac_GO@result$p.adjust < 0.001])
Fibroblast_pathway<-data.frame(cell="Fibroblast",pathway=Fibroblast_GO@result$Description[Fibroblast_GO@result$p.adjust < 0.001],pathwayID=Fibroblast_GO@result$ID[Fibroblast_GO@result$p.adjust < 0.001])

community1_pathway<-rbind(CD8_Tem_pathway,CD8_Tex_pathway,C1QC_mac_pathway,Fibroblast_pathway)
community1_pathway<-community1_pathway[,1:2]
result <- community1_pathway %>%  
  # 创建一个值为1的变量，表示激活  
  mutate(value = 1) %>%  
  # 使用pivot_wider将数据转换为宽格式  
  # names_from是我们要转换为列名的变量，values_fill用于填充缺失值（即0）  
  pivot_wider(names_from = pathway, values_from = value, values_fill = 0) 
community1_tmp<-as.matrix(result[2:ncol(result)])
community1_tmp[1,1]
rownames(community1_tmp)<-c("CD8 Tem","CD8 Tex","C1QC mac","Fibroblast")
community1_tmp<-t(community1_tmp)

community1_filtered<-community1_tmp[rowSums(community1_tmp)>=2,]
nrow(community1_filtered)

community1_pathway<-rbind(CD8_Tem_pathway,CD8_Tex_pathway,C1QC_mac_pathway,Fibroblast_pathway)
c1_pathway_id<- community1_pathway[,2:3]
c1_pathway_id<-unique(c1_pathway_id$pathwayID[c1_pathway_id$pathway %in% rownames(community1_filtered)])

pathway_gene_pair_c1 <- GO_data$PATHID2EXTID[c1_pathway_id]
pathway_gene_pair_c1<-reshape2::melt(pathway_gene_pair_c1)
colnames(pathway_gene_pair_c1)<-c("GENE","ID")
result <- pathway_gene_pair_c1 %>%  
  # 创建一个值为1的变量，表示激活  
  mutate(value = 1) %>%  
  # 使用pivot_wider将数据转换为宽格式  
  # names_from是我们要转换为列名的变量，values_fill用于填充缺失值（即0）  
  pivot_wider(names_from = GENE, values_from = value, values_fill = 0) 

rownames(community1_pathway)

community1_pathway<-as.matrix(result[2:ncol(result)])
rownames(community1_pathway)<-result$ID
library(stringr)
community1_pathway<-cor(t(community1_pathway))

d = dist(community1_pathway, method = 'euclidean')
tree = hclust(d, method = 'complete')
v = cutree(tree, 6)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)
gene.cluster <- as.data.frame(v)
gene.cluster$gene <- str_split(rownames(gene.cluster),'[.]',simplify = T)[,1]
table(gene.cluster$v)

#community1_pathway_cluster<-as.data.frame(p2$kmeans$cluster)

community1_pathway<-rbind(CD8_Tem_pathway,CD8_Tex_pathway,C1QC_mac_pathway,Fibroblast_pathway)
c1_pathway_id <- community1_pathway[,2:3]
c1_pathway_id<-unique(c1_pathway_id)
c1_pathway_id<-c1_pathway_id[c1_pathway_id$pathwayID %in% rownames(gene.cluster),]
gene.cluster$pathwayID<-rownames(gene.cluster)

community1_pathway_cluster<-merge(gene.cluster,c1_pathway_id,by="pathwayID")
View(community1_pathway_cluster)

community1_pathway_cluster_cand<-community1_pathway_cluster
community1_pathway_cluster_cand<-community1_pathway_cluster_cand[order(community1_pathway_cluster_cand$v),]
library(RColorBrewer)
brewer.pal(12,"Paired")
y_colour = list(cluster = c("cytokine production" = brewer.pal(12,"Paired")[1], "regulation of leukocyte cell-cell adhesion and proliferation"= brewer.pal(12,"Paired")[2], 
                            "cell killing"= brewer.pal(12,"Paired")[3], "leukocyte chemotaxis and  migration"= brewer.pal(12,"Paired")[4],
                            "homeostasis" = brewer.pal(12,"Paired")[5], "cell-matrix adhesion"= brewer.pal(12,"Paired")[6]))
annotation_row <- data.frame(  
  cluster = community1_pathway_cluster_cand$v
)  
rownames(annotation_row)<-community1_pathway_cluster_cand$pathway
rownames(community1_pathway_cluster_cand)<-community1_pathway_cluster_cand$pathway
community1_filtered
p1<-pheatmap::pheatmap(community1_filtered[community1_pathway_cluster_cand$pathway,],cluster_rows = F,cluster_cols = T,show_rownames = F,color = c("lightgrey","#5A80BE"),
                       annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA)
p1
annotation_row$cluster[annotation_row$cluster == 1]<-"cytokine production"
annotation_row$cluster[annotation_row$cluster == 2]<-"regulation of leukocyte cell-cell adhesion and proliferation"
annotation_row$cluster[annotation_row$cluster == 3]<-"cell killing"
annotation_row$cluster[annotation_row$cluster == 4]<-"leukocyte chemotaxis and  migration"
annotation_row$cluster[annotation_row$cluster == 5]<-"homeostasis"
annotation_row$cluster[annotation_row$cluster == 6]<-"cell-matrix adhesion"

community1_01_selected<-community1_filtered[community1_pathway_cluster_cand$pathway,]
p1<-pheatmap::pheatmap(community1_01_selected,cluster_rows = F,cluster_cols = T,show_rownames = F,color = c("lightgrey","#5A80BE"),
                       annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA)
p1

CD8_Tem_GO<-enrichGO(CD8_Tem_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
CD8_Tex_GO<-enrichGO(CD8_Tex_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
C1QC_mac_GO<-enrichGO(C1QC_mac_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Fibroblast_GO<-enrichGO(Fibroblast_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
CD8_Tem_GO<-CD8_Tem_GO[,c(2,5)]
CD8_Tem_GO$pvalue <- (-log10(CD8_Tem_GO$pvalue))
colnames(CD8_Tem_GO)[2]<-"CD8_Tem_p"
CD8_Tex_GO<-CD8_Tex_GO[,c(2,5)]
CD8_Tex_GO$pvalue <- (-log10(CD8_Tex_GO$pvalue))
colnames(CD8_Tex_GO)[2]<-"CD8_Tex_p"
C1QC_mac_GO<-C1QC_mac_GO[,c(2,5)]
C1QC_mac_GO$pvalue <- (-log10(C1QC_mac_GO$pvalue))
colnames(C1QC_mac_GO)[2]<-"C1QC_mac_p"
Fibroblast_GO<-Fibroblast_GO[,c(2,5)]
Fibroblast_GO$pvalue <- (-log10(Fibroblast_GO$pvalue))
colnames(Fibroblast_GO)[2]<-"Fibroblast_p"
merged_p<-merge(merge(merge(CD8_Tem_GO,CD8_Tex_GO,by="Description",all=T),C1QC_mac_GO,by="Description",all=T),Fibroblast_GO,by="Description",all=T)
merged_p[is.na(merged_p)]<-0
rownames(merged_p) <- merged_p$Description
merged_p_seleced<-merged_p[rownames(community1_01_selected),2:ncol(merged_p)]
as.numeric(merged_p_seleced)
community1_01_selected<-cbind(merged_p_seleced,annotation_row)
community1_01_selected$count<-rowSums(community1_01_selected[1:(ncol(community1_01_selected)-1)])
community1_01_selected<-community1_01_selected[order(community1_01_selected$cluster,community1_01_selected$count,decreasing = T),]
pheatmap::pheatmap(merged_p_seleced[rownames(community1_01_selected),],cluster_rows = F,cluster_cols = T,show_rownames = F,colorRampPalette(c("white","#5A80BE"))(50),
                       annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA)

