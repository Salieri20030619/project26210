.libPaths("/data/xy/R/R_package/R4/R_package_4")

sc_marker<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/PCC_in_real/Seurat_cluster_markers.rds")
names(table(sc_marker$cluster))

plot_data$community[plot_data$celltype %in% c("cDC3","Endothelial","PDPN+ Mac","FAP+ Myofib","Monocyte","Smooth Muscle")]<-"community1"
plot_data$community[plot_data$celltype %in% c("FN1+ Mac","CD4+ Tn","Treg","Monocyte","Neutrophil","cDC3","CD8+ Trm","Plasma","pDC","cDC1","B","ILC")]<-"community4"
plot_data$community[plot_data$celltype %in% c("C1QC+ Mac","CD8+ Tem","Fibroblast","CD8+ Tex","cDC2")]<-"community3"
plot_data$community[plot_data$celltype %in% c("NK","Th1/Th17","Mast","CD4+ Tm","Tfh")]<-"community4"



FN1_Mac_marker<-sc_marker$gene[sc_marker$cluster == "FN1+ Mac"]
CD4_Tn_marker<-sc_marker$gene[sc_marker$cluster == "CD4 Tn"]
Treg_marker<-sc_marker$gene[sc_marker$cluster == "Treg"]
Monocyte_marker<-sc_marker$gene[sc_marker$cluster == "Monocyte"]
Neutrophil_marker<-sc_marker$gene[sc_marker$cluster == "Neutrophil"]
cDC3_marker<-sc_marker$gene[sc_marker$cluster == "cDC3"]
CD8_Trm_marker<-sc_marker$gene[sc_marker$cluster == "CD8 Trm"]
Plasma_marker<-sc_marker$gene[sc_marker$cluster == "Plasma"]
pDC_marker<-sc_marker$gene[sc_marker$cluster == "pDC"]
cDC1_marker<-sc_marker$gene[sc_marker$cluster == "cDC1"]
B_marker<-sc_marker$gene[sc_marker$cluster == "B"]
ILC_marker<-sc_marker$gene[sc_marker$cluster == "ILC"]

FN1_Mac_GO<-enrichGO(FN1_Mac_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
CD4_Tn_GO<-enrichGO(CD4_Tn_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Treg_GO<-enrichGO(Treg_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Monocyte_GO<-enrichGO(Monocyte_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Neutrophil_GO<-enrichGO(Neutrophil_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
cDC3_GO<-enrichGO(cDC3_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
CD8_Trm_GO<-enrichGO(CD8_Trm_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
Plasma_GO<-enrichGO(Plasma_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
pDC_GO<-enrichGO(pDC_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
cDC1_GO<-enrichGO(cDC1_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
B_GO<-enrichGO(B_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
ILC_GO<-enrichGO(ILC_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")

FN1_Mac_pathway<-data.frame(cell="FN1 Mac",pathway=FN1_Mac_GO@result$Description[FN1_Mac_GO@result$p.adjust < 0.001],pathwayID=FN1_Mac_GO@result$ID[FN1_Mac_GO@result$p.adjust < 0.001])
CD4_Tn_pathway<-data.frame(cell="CD4 Tn",pathway=CD4_Tn_GO@result$Description[CD4_Tn_GO@result$p.adjust < 0.001],pathwayID=CD4_Tn_GO@result$ID[CD4_Tn_GO@result$p.adjust < 0.001])
Treg_pathway<-data.frame(cell="Treg",pathway=Treg_GO@result$Description[Treg_GO@result$p.adjust < 0.001],pathwayID=Treg_GO@result$ID[Treg_GO@result$p.adjust < 0.001])
Monocyte_pathway<-data.frame(cell="Monocyte",pathway=Monocyte_GO@result$Description[Monocyte_GO@result$p.adjust < 0.001],pathwayID=Monocyte_GO@result$ID[Monocyte_GO@result$p.adjust < 0.001])
Neutrophil_pathway<-data.frame(cell="Neutrophil",pathway=Neutrophil_GO@result$Description[Neutrophil_GO@result$p.adjust < 0.001],pathwayID=Neutrophil_GO@result$ID[Neutrophil_GO@result$p.adjust < 0.001])
cDC3_pathway<-data.frame(cell="cDC3",pathway=cDC3_GO@result$Description[cDC3_GO@result$p.adjust < 0.001],pathwayID=cDC3_GO@result$ID[cDC3_GO@result$p.adjust < 0.001])
CD8_Trm_pathway<-data.frame(cell="CD8 Trm",pathway=CD8_Trm_GO@result$Description[CD8_Trm_GO@result$p.adjust < 0.001],pathwayID=CD8_Trm_GO@result$ID[CD8_Trm_GO@result$p.adjust < 0.001])
Plasma_pathway<-data.frame(cell="Plasma",pathway=Plasma_GO@result$Description[Plasma_GO@result$p.adjust < 0.001],pathwayID=Plasma_GO@result$ID[Plasma_GO@result$p.adjust < 0.001])
pDC_pathway<-data.frame(cell="pDC",pathway=pDC_GO@result$Description[pDC_GO@result$p.adjust < 0.001],pathwayID=pDC_GO@result$ID[pDC_GO@result$p.adjust < 0.001])
cDC1_pathway<-data.frame(cell="cDC1",pathway=cDC1_GO@result$Description[cDC1_GO@result$p.adjust < 0.001],pathwayID=cDC1_GO@result$ID[cDC1_GO@result$p.adjust < 0.001])
B_pathway<-data.frame(cell="B",pathway=B_GO@result$Description[B_GO@result$p.adjust < 0.001],pathwayID=B_GO@result$ID[B_GO@result$p.adjust < 0.001])
ILC_pathway<-data.frame(cell="ILC",pathway=ILC_GO@result$Description[ILC_GO@result$p.adjust < 0.001],pathwayID=ILC_GO@result$ID[ILC_GO@result$p.adjust < 0.001])

community4_pathway<-rbind(Monocyte_pathway,Neutrophil_pathway)
community4_pathway<-community4_pathway[,1:2]
result <- community4_pathway %>%  
  # 创建一个值为1的变量，表示激活  
  mutate(value = 1) %>%  
  # 使用pivot_wider将数据转换为宽格式  
  # names_from是我们要转换为列名的变量，values_fill用于填充缺失值（即0）  
  pivot_wider(names_from = pathway, values_from = value, values_fill = 0) 
community4_tmp<-as.matrix(result[2:ncol(result)])
community4_tmp[1,1]
rownames(community4_tmp)<-c("Monocyte","Neutrophil")
community4_tmp<-t(community4_tmp)

community4_filtered<-community4_tmp[rowSums(community4_tmp)>=2,]
nrow(community4_filtered)

community4_pathway<-rbind(Monocyte_pathway,Neutrophil_pathway)
c4_pathway_id<- community4_pathway[,2:3]
c4_pathway_id<-unique(c4_pathway_id$pathwayID[c4_pathway_id$pathway %in% rownames(community4_filtered)])

pathway_gene_pair_c4 <- GO_data$PATHID2EXTID[c4_pathway_id]
pathway_gene_pair_c4<-reshape2::melt(pathway_gene_pair_c4)
colnames(pathway_gene_pair_c4)<-c("GENE","ID")
result <- pathway_gene_pair_c4 %>%  
  # 创建一个值为1的变量，表示激活  
  mutate(value = 1) %>%  
  # 使用pivot_wider将数据转换为宽格式  
  # names_from是我们要转换为列名的变量，values_fill用于填充缺失值（即0）  
  pivot_wider(names_from = GENE, values_from = value, values_fill = 0) 

rownames(community4_pathway)

community4_pathway<-as.matrix(result[2:ncol(result)])
rownames(community4_pathway)<-result$ID
library(stringr)
community4_pathway<-cor(t(community4_pathway))

d = dist(community4_pathway, method = 'euclidean')
tree = hclust(d, method = 'complete')
v = cutree(tree, 6)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)
gene.cluster <- as.data.frame(v)
gene.cluster$gene <- str_split(rownames(gene.cluster),'[.]',simplify = T)[,1]
table(gene.cluster$v)

#community4_pathway_cluster<-as.data.frame(p2$kmeans$cluster)

community4_pathway<-rbind(Monocyte_pathway,Neutrophil_pathway)
c4_pathway_id <- community4_pathway[,2:3]
c4_pathway_id<-unique(c4_pathway_id)
c4_pathway_id<-c4_pathway_id[c4_pathway_id$pathwayID %in% rownames(gene.cluster),]
gene.cluster$pathwayID<-rownames(gene.cluster)

community4_pathway_cluster<-merge(gene.cluster,c4_pathway_id,by="pathwayID")
View(community4_pathway_cluster)

community4_pathway_cluster_cand<-community4_pathway_cluster
community4_pathway_cluster_cand<-community4_pathway_cluster_cand[order(community4_pathway_cluster_cand$v),]
library(RColorBrewer)
brewer.pal(12,"Paired")
y_colour = list(cluster = c("cytokine production" = brewer.pal(12,"Paired")[1], "cytokine-mediated signaling pathway"= brewer.pal(12,"Paired")[2], 
                            "regulation of leukocyte cell-cell adhesion and proliferation"= brewer.pal(12,"Paired")[3], 
                            "response to external stimulus"= brewer.pal(12,"Paired")[4],
                            "regulation of myeloid cell differentiation" = brewer.pal(12,"Paired")[5], "acute inflammatory response"= brewer.pal(12,"Paired")[6]))
annotation_row <- data.frame(  
  cluster = community4_pathway_cluster_cand$v
)  
rownames(annotation_row)<-community4_pathway_cluster_cand$pathway
rownames(community4_pathway_cluster_cand)<-community4_pathway_cluster_cand$pathway

annotation_row$cluster[annotation_row$cluster == 1]<-"cytokine production"
annotation_row$cluster[annotation_row$cluster == 2]<-"cytokine-mediated signaling pathway"
annotation_row$cluster[annotation_row$cluster == 3]<-"regulation of leukocyte cell-cell adhesion and proliferation"
annotation_row$cluster[annotation_row$cluster == 4]<-"response to external stimulus"
annotation_row$cluster[annotation_row$cluster == 5]<-"regulation of myeloid cell differentiation"
annotation_row$cluster[annotation_row$cluster == 6]<-"acute inflammatory response"

p2<-pheatmap::pheatmap(community4_filtered[community4_pathway_cluster_cand$pathway,],show_rownames = F,color = c("#5A80BE"),
                       annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA,breaks = c(0,1))



community4_01_selected<-community4_filtered[community4_pathway_cluster_cand$pathway,]
c("Monocyte","Neutrophil")

Monocyte_GO<-enrichGO(Monocyte_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Neutrophil_GO<-enrichGO(Neutrophil_marker,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 0)
Monocyte_GO<-Monocyte_GO[,c(2,5)]
Monocyte_GO$pvalue <- (-log10(Monocyte_GO$pvalue))
colnames(Monocyte_GO)[2]<-"Monocyte_p"
Neutrophil_GO<-Neutrophil_GO[,c(2,5)]
Neutrophil_GO$pvalue <- (-log10(Neutrophil_GO$pvalue))
colnames(Neutrophil_GO)[2]<-"Neutrophil_p"

merged_p<-merge(Monocyte_GO,Neutrophil_GO,by="Description",all=T)
merged_p[is.na(merged_p)]<-0
rownames(merged_p) <- merged_p$Description
merged_p_seleced<-merged_p[rownames(community4_01_selected),2:ncol(merged_p)]
as.numeric(merged_p_seleced)
community4_01_selected<-cbind(merged_p_seleced,annotation_row)
community4_01_selected$count<-rowSums(community4_01_selected[1:(ncol(community4_01_selected)-1)])
community4_01_selected<-community4_01_selected[order(community4_01_selected$cluster,community4_01_selected$count,decreasing = T),]
pheatmap::pheatmap(merged_p_seleced[rownames(community4_01_selected),],cluster_rows = F,cluster_cols = T,show_rownames = F,colorRampPalette(c("white","#5A80BE"))(50),
                   annotation_row = annotation_row,annotation_colors = y_colour,border_color = NA)

