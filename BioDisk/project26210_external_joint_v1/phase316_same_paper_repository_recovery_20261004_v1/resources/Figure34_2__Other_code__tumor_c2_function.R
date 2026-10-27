.libPaths("/data/xy/R/R_package/R4/R_package_4")

sc_marker<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/PCC_in_real/Seurat_cluster_markers.rds")

FAP_Myofib_marker<-sc_marker$gene[sc_marker$cluster == "FAP+ Myofib"]
PDPN_mac_marker<-sc_marker$gene[sc_marker$cluster == "PDPN+ Mac"]
Smooth_Muscle_marker<-sc_marker$gene[sc_marker$cluster == "Smooth Muscle"]
Pericyte_marker<-sc_marker$gene[sc_marker$cluster == "Pericyte"]
Endothelium_marker<-sc_marker$gene[sc_marker$cluster == "Endothelium"]
Glial_marker<-sc_marker$gene[sc_marker$cluster == "Glial"]

FAP_Myofib_GO<-enrichGO(FAP_Myofib_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
PDPN_mac_GO<-enrichGO(PDPN_mac_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Smooth_Muscle_GO<-enrichGO(Smooth_Muscle_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Pericyte_GO<-enrichGO(Pericyte_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Endothelium_GO<-enrichGO(Endothelium_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Glial_GO<-enrichGO(Glial_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")

FAP_Myofib_pathway<-data.frame(cell="FAP Myofib",pathway=FAP_Myofib_GO@result$Description[FAP_Myofib_GO@result$p.adjust < 0.001],pathwayID=FAP_Myofib_GO@result$ID[FAP_Myofib_GO@result$p.adjust < 0.001])
PDPN_mac_pathway<-data.frame(cell="PDPN mac",pathway=PDPN_mac_GO@result$Description[PDPN_mac_GO@result$p.adjust < 0.001],pathwayID=PDPN_mac_GO@result$ID[PDPN_mac_GO@result$p.adjust < 0.001])
Smooth_Muscle_pathway<-data.frame(cell="Smooth Muscle",pathway=Smooth_Muscle_GO@result$Description[Smooth_Muscle_GO@result$p.adjust < 0.001],pathwayID=Smooth_Muscle_GO@result$ID[Smooth_Muscle_GO@result$p.adjust < 0.001])
Pericyte_pathway<-data.frame(cell="Pericyte",pathway=Pericyte_GO@result$Description[Pericyte_GO@result$p.adjust < 0.001],pathwayID=Pericyte_GO@result$ID[Pericyte_GO@result$p.adjust < 0.001])
Endothelium_pathway<-data.frame(cell="Endothelium",pathway=Endothelium_GO@result$Description[Endothelium_GO@result$p.adjust < 0.001],pathwayID=Endothelium_GO@result$ID[Endothelium_GO@result$p.adjust < 0.001])
Glial_pathway<-data.frame(cell="Glial",pathway=Glial_GO@result$Description[Glial_GO@result$p.adjust < 0.001],pathwayID=Glial_GO@result$ID[Glial_GO@result$p.adjust < 0.001])
community3_pathway<-rbind(FAP_Myofib_pathway,PDPN_mac_pathway,Smooth_Muscle_pathway,Pericyte_pathway,Endothelium_pathway,Glial_pathway,FN1_Mac_pathway)
community3_pathway<-community3_pathway[,1:2]
result <- community3_pathway %>%  
  # 创建一个值为1的变量，表示激活  
  mutate(value = 1) %>%  
  # 使用pivot_wider将数据转换为宽格式  
  # names_from是我们要转换为列名的变量，values_fill用于填充缺失值（即0）  
  pivot_wider(names_from = pathway, values_from = value, values_fill = 0) 
community3_tmp<-as.matrix(result[2:ncol(result)])
community3_tmp[1,1]
rownames(community3_tmp)<-c("FAP Myofib","PDPN mac","Smooth Muscle","Pericyte","Endothelium","Glial","FN1 mac")
community3_tmp<-t(community3_tmp)

community3_filtered<-community3_tmp[rowSums(community3_tmp)>=3.5,]
nrow(community3_filtered)

community3_pathway<-rbind(FAP_Myofib_pathway,PDPN_mac_pathway,Smooth_Muscle_pathway,Pericyte_pathway,Endothelium_pathway,Glial_pathway,FN1_Mac_pathway)
c3_pathway_id<- community3_pathway[,2:3]
c3_pathway_id<-unique(c3_pathway_id$pathwayID[c3_pathway_id$pathway %in% rownames(community3_filtered)])

pathway_gene_pair_c3 <- GO_data$PATHID2EXTID[c3_pathway_id]
pathway_gene_pair_c3<-reshape2::melt(pathway_gene_pair_c3)
colnames(pathway_gene_pair_c3)<-c("GENE","ID")
result <- pathway_gene_pair_c3 %>%  
  # 创建一个值为1的变量，表示激活  
  mutate(value = 1) %>%  
  # 使用pivot_wider将数据转换为宽格式  
  # names_from是我们要转换为列名的变量，values_fill用于填充缺失值（即0）  
  pivot_wider(names_from = GENE, values_from = value, values_fill = 0) 

rownames(community3_pathway)

community3_pathway<-as.matrix(result[2:ncol(result)])
rownames(community3_pathway)<-result$ID
library(stringr)
community3_pathway<-cor(t(community3_pathway))

d = dist(community3_pathway, method = 'euclidean')
tree = hclust(d, method = 'complete')
v = cutree(tree, 10)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)
gene.cluster <- as.data.frame(v)
gene.cluster$gene <- str_split(rownames(gene.cluster),'[.]',simplify = T)[,1]
table(gene.cluster$v)

#community3_pathway_cluster<-as.data.frame(p2$kmeans$cluster)

community3_pathway<-rbind(FAP_Myofib_pathway,PDPN_mac_pathway,Smooth_Muscle_pathway,Pericyte_pathway,Endothelium_pathway,Glial_pathway,FN1_Mac_pathway)
c3_pathway_id <- community3_pathway[,2:3]
c3_pathway_id<-unique(c3_pathway_id)
c3_pathway_id<-c3_pathway_id[c3_pathway_id$pathwayID %in% rownames(gene.cluster),]
gene.cluster$pathwayID<-rownames(gene.cluster)

community3_pathway_cluster<-merge(gene.cluster,c3_pathway_id,by="pathwayID")
View(community3_pathway_cluster)

community3_pathway_cluster_cand<-community3_pathway_cluster
community3_pathway_cluster_cand<-community3_pathway_cluster_cand[order(community3_pathway_cluster_cand$v),]
library(RColorBrewer)
brewer.pal(12,"Paired")
y_colour = list(cluster = c("extracellular matrix organization" = brewer.pal(12,"Paired")[1], "cell-matrix adhesion"= brewer.pal(12,"Paired")[2], 
                            "TGFbeta signaling"= brewer.pal(12,"Paired")[3], "wound healing"= brewer.pal(12,"Paired")[4],
                            "angiogenesis" = brewer.pal(12,"Paired")[5], "tissue remodeling"= brewer.pal(12,"Paired")[6], 
                            "WNT signaling"= brewer.pal(12,"Paired")[7], "actin filament organization"= brewer.pal(12,"Paired")[8],
                            "symbiotic interaction" = brewer.pal(12,"Paired")[9], "neurogenesis"= brewer.pal(12,"Paired")[10]))
annotation_row <- data.frame(  
  cluster = community3_pathway_cluster_cand$v
)  
rownames(annotation_row)<-community3_pathway_cluster_cand$pathway
rownames(community3_pathway_cluster_cand)<-community3_pathway_cluster_cand$pathway
annotation_row$cluster[annotation_row$cluster == 1]<-"extracellular matrix organization"
annotation_row$cluster[annotation_row$cluster == 2]<-"cell-matrix adhesion"
annotation_row$cluster[annotation_row$cluster == 3]<-"TGFbeta signaling"
annotation_row$cluster[annotation_row$cluster == 4]<-"wound healing"
annotation_row$cluster[annotation_row$cluster == 5]<-"angiogenesis"
annotation_row$cluster[annotation_row$cluster == 6]<-"tissue remodeling"
annotation_row$cluster[annotation_row$cluster == 7]<-"WNT signaling"
annotation_row$cluster[annotation_row$cluster == 8]<-"actin filament organization"
annotation_row$cluster[annotation_row$cluster == 9]<-"symbiotic interaction"
annotation_row$cluster[annotation_row$cluster == 10]<-"neurogenesis"
p3<-pheatmap::pheatmap(community3_filtered[community3_pathway_cluster_cand$pathway,],cluster_rows = F,cluster_cols = T,show_rownames = F,color = c("lightgrey","#5A80BE"),
                       annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA)

(FAP_Myofib_pathway,PDPN_mac_pathway,Smooth_Muscle_pathway,Pericyte_pathway,Endothelium_pathway,Glial_pathway,FN1_Mac_pathway)
community3_01_selected<-community3_filtered[community3_pathway_cluster_cand$pathway,]
FAP_Myofib_GO<-enrichGO(FAP_Myofib_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
PDPN_mac_GO<-enrichGO(PDPN_mac_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Smooth_Muscle_GO<-enrichGO(Smooth_Muscle_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Pericyte_GO<-enrichGO(Pericyte_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Endothelium_GO<-enrichGO(Endothelium_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Glial_GO<-enrichGO(Glial_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
FN1_Mac_GO<-enrichGO(FN1_Mac_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)

FAP_Myofib_GO<-FAP_Myofib_GO[,c(2,5)]
FAP_Myofib_GO$pvalue <- (-log10(FAP_Myofib_GO$pvalue))
colnames(FAP_Myofib_GO)[2]<-"FAP_Myofib_p"
PDPN_mac_GO<-PDPN_mac_GO[,c(2,5)]
PDPN_mac_GO$pvalue <- (-log10(PDPN_mac_GO$pvalue))
colnames(PDPN_mac_GO)[2]<-"PDPN_mac_p"
Smooth_Muscle_GO<-Smooth_Muscle_GO[,c(2,5)]
Smooth_Muscle_GO$pvalue <- (-log10(Smooth_Muscle_GO$pvalue))
colnames(Smooth_Muscle_GO)[2]<-"Smooth_Muscle_p"
Pericyte_GO<-Pericyte_GO[,c(2,5)]
Pericyte_GO$pvalue <- (-log10(Pericyte_GO$pvalue))
colnames(Pericyte_GO)[2]<-"Pericyte_p"
Endothelium_GO<-Endothelium_GO[,c(2,5)]
Endothelium_GO$pvalue <- (-log10(Endothelium_GO$pvalue))
colnames(Endothelium_GO)[2]<-"Endothelium_p"
Glial_GO<-Glial_GO[,c(2,5)]
Glial_GO$pvalue <- (-log10(Glial_GO$pvalue))
colnames(Glial_GO)[2]<-"Glial_p"
FN1_Mac_GO<-FN1_Mac_GO[,c(2,5)]
FN1_Mac_GO$pvalue <- (-log10(FN1_Mac_GO$pvalue))
colnames(FN1_Mac_GO)[2]<-"FN1_Mac_p"

merged_p<-merge(merge(merge(merge(merge(merge(FAP_Myofib_GO,PDPN_mac_GO,by="Description",all=T),Smooth_Muscle_GO,by="Description",all=T),Pericyte_GO,by="Description",all=T),
                                  Endothelium_GO,by="Description",all=T),Glial_GO,by="Description",all=T),FN1_Mac_GO,by="Description",all=T)
merged_p[is.na(merged_p)]<-0
rownames(merged_p) <- merged_p$Description
merged_p_seleced<-merged_p[rownames(community3_01_selected),2:ncol(merged_p)]
community2_01_selected<-cbind(merged_p_seleced,annotation_row)
community2_01_selected$count<-rowSums(community2_01_selected[1:(ncol(community2_01_selected)-1)])
community2_01_selected<-community2_01_selected[order(community2_01_selected$cluster,community2_01_selected$count,decreasing = T),]
pheatmap::pheatmap(merged_p_seleced[rownames(community2_01_selected),],cluster_rows = F,cluster_cols = T,show_rownames = F,colorRampPalette(c("white","#5A80BE"))(50),
                   annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA)
