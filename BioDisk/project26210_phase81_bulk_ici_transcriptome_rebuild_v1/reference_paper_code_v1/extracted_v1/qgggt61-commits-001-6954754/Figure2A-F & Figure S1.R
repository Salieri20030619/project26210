##########Figure2A-F & Figure S1
sample_objects <- purrr::map(1:length(sample_list), function(x) {  
  one_dir <- samples_dir[x]  
  sample_id <- sample_list[x]  
  slice_id <- sample_names[x]  
  
  sample_object <- Load10X_Spatial(
    data.dir = one_dir,  
    filename = "filtered_feature_bc_matrix.h5",  
    assay = "Spatial",  
    slice = slice_id,  
    filter.matrix = TRUE  
  )  
  
  sample_object@project.name <- sample_id  
  sample_object@meta.data$orig.ident <- slice_id  
  sample_object <- RenameCells(object = sample_object, add.cell.id = slice_id)
  
  return(sample_object)
})
sample_objects <- lapply(sample_objects, 
                         SCTransform, 
                         assay = "Spatial", 
                         method = "poisson")
ST <- merge(sample_objects[[1]], y = sample_objects[[2]])
DefaultAssay(ST) <- "SCT"

VariableFeatures(ST) <- c(VariableFeatures(sample_objects[[1]]),  
                          VariableFeatures(sample_objects[[2]]))

ST <- RunPCA(ST, assay = "SCT", verbose = FALSE)
ST <- RunHarmony(ST, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
ST <- FindNeighbors(ST, reduction = "harmony", dims = 1:30) %>%  
  FindClusters(resolution = 0.5)

ST <- RunUMAP(ST, dims = 1:30)
display.brewer.all()
cols <- colorRampPalette(brewer.pal(9, 'Set1'))(26)

DimPlot(ST, reduction = "umap", cols = cols)
SpatialPlot(ST, group.by = 'seurat_clusters', ncol = 2,label = T) & 
  scale_fill_manual(values = cols)
SpatialFeaturePlot(ST, features = c("SPP1","APOE","CD163"), alpha = c(0.1, 1))

ST$seurat_clusters <- factor(ST$seurat_clusters) 
Idents(ST) <- ST$seurat_clusters
ST <- PrepSCTFindMarkers(ST)

markers <- FindAllMarkers(
  ST, 
  only.pos = TRUE,  
  min.pct = 0.25,   
  logfc.threshold = 0.25  
)
top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) 
write.csv(markers, file = "markers.csv")
write.csv(top10, file = "top10.csv")
top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) 
write.csv(top10, file = "top20.csv")

library(scRNAtoolVis)
cellRatioPlot(object = ST, 
              sample.name = "orig.ident", 
              celltype.name = "seurat_clusters", 
              flow.curve = 0.5, 
              fill.col = cols) + 
  RotatedAxis() 

library(clusterProfiler)
library(dplyr)
gid <- bitr(unique(markers$gene), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
markers <- full_join(markers, gid, by = c('gene' = 'SYMBOL'))

dotplot(KEGG, label_format = 40) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient(high = "#4b5cc4", low = "#FE8D3C")

ggsave('KEGG.pdf',height = 10,width = 10)
dev.off()

av <- AverageExpression(ST, group.by = "seurat_clusters")
av <- av[[1]] 
head(av)
cg <- names(tail(sort(apply(av, 1, sd)), 1000))
cor_matrix <- cor(av[cg,], method = 'spearman')
pheatmap::pheatmap(cor_matrix)
ggsave('Correlation.pdf',height = 10,width = 10)
dev.off()



genes <- c("FAP", "POSTN", "COL1A1", "ACTA2", "CD3D", "CD3E", "CA9", 
           'CD68','CD163','CD14','CD16','APOE','SPP1','CXCL9',
           'KRT19','EPCAM'
           
           
)
p_all_markers <- DotPlot(
  ST,
  features = genes,
  scale = TRUE,
  
  group.by = "seurat_clusters"
) +
  theme_bw() + 
  scale_color_continuous(low = "grey", high = "#1F77B4FF") +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 0, unit = 'cm'),
    axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 12),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black")
  )
p_all_markers
ggsave('check_paper_markers_cluster2.pdf', height = 2, width = 14)
SpatialDimPlot(ST, 
               ncol = 1,
               cells.highlight = CellsByIdentities(object = ST, idents = c('9')))
SpatialDimPlot(ST, 
               ncol = 1,
               cells.highlight = CellsByIdentities(object = ST, idents = c('15')))

new.cluster.ids <- c("0"="Malignant_immune", 
                     "1"="CAF_Endo", 
                     "2"="B_CAF", 
                     "3"="TAM_CAF_Tumor", 
                     "4"="APOE+TAM_POSTN+CAF", 
                     "5"="Malignant_Endo", 
                     "6"="Malignant_Endo", 
                     "7"="Ma_CAF_TAM", 
                     "8"="Malignant_immune",
                     "9"="immune", 
                     "10"="immune",
                     "11"="SMC_CAF",
                     "12"="Malignant_Myeloid", 
                     "13"="SMC_CAF_Endo",
                     "14"="immune", 
                     "15"="Malignant_immune"
)
ST <- RenameIdents(ST, new.cluster.ids)                        
ST$celltype <- ST@active.ident
DimPlot(ST, label = T,pt.size = 1,group.by = "celltype")
SpatialPlot(ST, group.by = 'celltype', ncol = 2,label = F) & 
  scale_fill_manual(values = cols)

SpatialDimPlot(ST, 
               ncol = 1,
               cells.highlight = CellsByIdentities(object = ST, idents = c('APOE+TAM_POSTN+CAF')))
colors <- c('#FFEBCD','#AA8282','#FFD700','#FFC0CB','#DC143C','#5000B8','#008080','#9CED5F',
            '#87CEEB','#007FFF','#FFBF00')
SpatialPlot(ST, group.by = 'celltype', ncol = 2,label = F) & 
  scale_fill_manual(values = colors)

SpatialPlot(ST, group.by = 'celltype', ncol = 2, label = F) & 
  scale_fill_manual(values = colors) &
  theme(legend.position = "right") & 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) &
  theme(legend.title = element_blank(), legend.text = element_text(size = 12)) &
  guides(fill = guide_legend(override.aes = list(stroke = 0)))
ggsave('ST.plot.pdf',height = 3,width = 8)

ggsave('ST.plot2.pdf',height = 8,width = 15)
save.image('20250218.Figure2.Rdata')
ST$celltype <- factor(ST$celltype)  
Idents(ST) <- ST$celltype
ST$region<-ifelse(ST$celltype=='APOE+TAM_POSTN+CAF','APOE+TAM_POSTN+CAF','other')
Idents(ST)<-'region'
ST <- PrepSCTFindMarkers(ST)
deg1up <- rownames(markers[which(markers$avg_log2FC > 0.1 & markers$p_val_adj  < 0.5),])
deg1dw <- rownames(markers[which(markers$avg_log2FC < -0.1& markers$p_val_adj  < 0.5),])
genes<-c(deg1up,deg1dw)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
library(dplyr)

alldiff<-markers[deg1up,]

alldiff$id <- rownames(alldiff)
all<- alldiff[,c(2,8)]%>%
  arrange(desc(alldiff$avg_log2FC))

all_id <- bitr(all$id,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")%>%
  inner_join(all,by= c('SYMBOL'='id'))%>%
  arrange(desc(avg_log2FC))  
geneList1 <- all_id$avg_log2FC
names(geneList1) <- all_id$ENTREZID
geneList1 <- sort(geneList1, decreasing = TRUE)
print(geneList1)
ego3 <- gseGO(geneList     = geneList1,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              verbose      = FALSE)
volcanoGsea(data = ego3)
library(aPEAR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(cols4all)
dt<- geneList1
head(dt) 
set.seed(123)
enrich<- gseGO(dt, OrgDb = org.Hs.eg.db, ont = 'BP')

geneList <- data.frame(geneList1)
genes <-deg1up
genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
goresult <- data.frame(go_res@result)
go_select <- goresult %>% 
  filter(ONTOLOGY == "MF" & pvalue < 0.05) %>%  
  arrange(pvalue) %>%  
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))  
p <- ggplot(go_select, aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", fill = "#7f4ea850", width = 0.8) +  
  geom_text(aes(x = 0.2, label = Description), hjust = 0, vjust = 1, size = 4.5) + 
  theme_minimal() +  
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.75),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "-log10(pvalue)", y = "") +  
  coord_cartesian(clip = 'off')  
p
ggsave("gomf.pdf", p, width = 4, height = 4.5)
go_select$pval_log <- -log10(go_select$pvalue)
go_select$RichFactor <- as.numeric(gsub("/.*", "", go_select$GeneRatio)) / 
  as.numeric(gsub("/.*", "", go_select$BgRatio))
ggsave("gomf-1.pdf", p, width = 5.2, height = 4.5)

go_select <- goresult %>% 
  filter(ONTOLOGY == "CC" & pvalue < 0.05) %>%  
  arrange(pvalue) %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description)))) 
p <- ggplot(go_select, aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", fill = "#7f4ea850", width = 0.8) +  
  geom_text(aes(x = 0.2, label = Description), hjust = 0, vjust = 1, size = 4.5) +  
  theme_minimal() + 
  axis.text.y = element_blank(),
axis.title.x = element_text(margin = margin(t = 10)),
axis.title.y = element_text(margin = margin(r = 10)),
plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
axis.line = element_line(linewidth = 0.75, color = "black"),
axis.ticks.x = element_line(linewidth = 0.75),
axis.ticks.y = element_blank(),
axis.text = element_text(size = 12, color = "black"),
axis.title = element_text(size = 12, color = "black"),
plot.title = element_text(size = 12, color = "black"),
panel.grid = element_blank()
) +
  labs(x = "-log10(pvalue)", y = "") +  
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = 'off')  
ggsave("goBP.pdf", p, width = 4, height = 4.5)
go_select$pval_log <- -log10(go_select$pvalue)
go_select$RichFactor <- as.numeric(gsub("/.*", "", go_select$GeneRatio)) / 
  as.numeric(gsub("/.*", "", go_select$BgRatio))
head(go_select)
scale_factor <- max(go_select$pval_log) / max(go_select$RichFactor)
p3 <- ggplot(go_select, aes(y = reorder(Description, pval_log))) +
  geom_bar(aes(x = pval_log), stat = "identity", fill = "#FDBF6F") + 
  geom_path(aes(x = RichFactor * scale_factor, y = reorder(Description, pval_log), group = 1),
            color = "black", size = 0.75) +  
  geom_point(aes(x = RichFactor * scale_factor),
             shape = 21, fill = "#FDBF6F", color = "black", size = 3) +  
  scale_x_continuous(
    name = "-log10(pvalue)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Rich Factor")
  ) +
  labs(y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
print(p3)
keggresult <- data.frame(kegg_res@result)
write.csv(keggresult,'keggresult.csv')
keggresult <- keggresult %>% 
  arrange(pvalue) %>% 
  head(10) %>%  
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))  
go_select <- keggresult %>% 
  filter( pvalue < 0.05) %>%  
  arrange(pvalue) %>%  
  head(10) %>%  
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))  
p <- ggplot(go_select, aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", fill = "#7f4ea850", width = 0.8) +  
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.75),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "-log10(pvalue)", y = "") +  
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = 'off') 
p
ggsave("goBP.pdf", p, width = 4, height = 4.5)
go_select$pval_log <- -log10(go_select$pvalue)
go_select$RichFactor <- as.numeric(gsub("/.*", "", go_select$GeneRatio)) / 
  as.numeric(gsub("/.*", "", go_select$BgRatio))
head(go_select)
scale_factor <- max(go_select$pval_log) / max(go_select$RichFactor)
p4 <- ggplot(go_select, aes(y = reorder(Description, pval_log))) +
  geom_bar(aes(x = pval_log), stat = "identity", fill = "#8DD3C7") + 
  geom_path(aes(x = RichFactor * scale_factor, y = reorder(Description, pval_log), group = 1),
            color = "black", size = 0.75) +  
  geom_point(aes(x = RichFactor * scale_factor),
             shape = 21, fill = "#8DD3C7", color = "black", size = 3) + 
  scale_x_continuous(
    name = "-log10(pvalue)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Rich Factor")
  ) +
  labs(y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
print(p4)
p <- p2 | p3| p1 | p4
ggsave("go.pdf", p, width = 25, height = 4.5)

library(tidyverse)
library(igraph)
library(ggraph)
gsea_df <- as.data.frame(enrich)  
library(tidyr)
pathway_gene <- gsea_df %>%
  dplyr::select(ID, Description, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/") 
library(igraph)
pathway_gene_list <- split(pathway_gene$core_enrichment, pathway_gene$Description)
jaccard_matrix <- similarity(graph_from_data_frame(pathway_gene), method = "jaccard", mode = "all")
jaccard_df <- as_data_frame(graph.adjacency(jaccard_matrix, weighted = TRUE), what = "edges") %>%
  filter(weight > 0.1 & weight < 1)  
nodes <- gsea_df %>%
  dplyr::select(Description, NES, p.adjust) %>%
  dplyr::rename(name = Description)
edges <- jaccard_df %>%
  dplyr::rename(from = from, to = to)
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)


library(ggraph)
library(scales)
ggraph(net, layout = "fr") +  
  geom_edge_link(aes(width = weight), alpha = 0.3) +
  geom_node_point(aes(size = abs(NES), color = -log10(p.adjust))) +
  geom_node_text(aes(label = name), repel = TRUE)
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggforce)
V(net)$cluster <- cluster_louvain(net)$membership
layout_df <- create_layout(net, layout = "fr")
ggraph(layout_df) +
  geom_edge_link(aes(width = weight), edge_colour = "grey80", alpha = 0.3) +
  geom_mark_hull(
    aes(x = x, y = y, group = cluster, fill = as.factor(cluster), label = NULL),
    concavity = 5, expand = unit(2, "mm"), alpha = 0.2, color = "black"
  ) +
  geom_node_point(aes(size = abs(NES), color = -log10(p.adjust)), show.legend = TRUE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 0.1) +
  scale_size_continuous(range = c(0.01, 0.1)) +
  scale_color_gradient(low = "skyblue", high = "firebrick") +
  scale_fill_brewer(palette = "Set2") +
  theme_graph() +
  theme(legend.position = "right")
View(enrich@result)
set.seed(1)

p <- enrichmentNetwork(
  enrich@result,
  colorBy = 'NES',
  colorType = c("nes"),
  nodeSize = 'setSize',
  fontSize = 5,
  plotOnly = FALSE,
  verbose = TRUE
)
p
set.seed(1)
p1<- enrichmentNetwork(enrich@result,
                       colorBy= 'NES',
                       colorType= c("nes"),
                       nodeSize= 'setSize',
                       fontSize= 5,
                       drawEllipses= TRUE, 
                       verbose= TRUE)
p1
set.seed(1)
p2<- enrichmentNetwork(enrich@result,
                       colorBy= 'pvalue',
                       colorType= c("pval"), 
                       nodeSize= 'setSize',
                       fontSize= 4,
                       drawEllipses= TRUE,
                       pCutoff= -10,
                       verbose= TRUE)
p2
library(ggplot2)
p1 + scale_color_gradientn(
  colors = c("#C6D7ED", "white", "#A78FB3"),
  limits = c(-2, 2),  
  name = "NES"
)
top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) 
write.csv(markers, file = "markers.csv")
write.csv(top10, file = "top10.csv")
library(ggrepel)
library(dplyr)
library(ggrepel)
log2FC_cutoff <- 0.585
padj_cutoff <- 0.05
markers <- markers %>%
  mutate(
    change = case_when(
      p_val_adj < padj_cutoff & avg_log2FC > log2FC_cutoff ~ "Up",
      p_val_adj < padj_cutoff & avg_log2FC < -log2FC_cutoff ~ "Down",
      TRUE ~ "NotSig"
    )
  )
markers1<-markers
markers <- markers %>%
  mutate(sign = ifelse(p_val_adj < 0.005 & abs(avg_log2FC) > 2, rownames(markers), NA))
ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = change), alpha = 0.4, size = 3) +
  scale_color_manual(
    name = "Expression",
    values = c("Up" = "#5000B8", "Down" = "#008080", "NotSig" = "grey"),
    labels = c(
      paste0("Down (", n_down, ")"),
      "Not Sig",
      paste0("Up (", n_up, ")")
    )
  ) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  geom_label_repel(
    aes(label = Sign),
    fontface = "bold", color = "black",
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "grey50"
  ) +
  labs(x = "log2(Fold Change)", y = "-log10(FDR adjusted P-value)") +
  theme_bw()


top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) 
write.csv(top10, file = "top20.csv")
DoHeatmap(
  ST, 
  features = top10$gene,  
) + NoLegend()  
top10
DotPlot(ST,features = top10$gene_name)

library(scRNAtoolVis)
library(Seurat)
pdf('Degs.pdf',width = 9,height = 7)
jjVolcano(diffData = all.markers2,
          tile.col = corrplot::COL2('RdBu', 15)[1:12],
          size  = 3.5,
          fontface = 'italic',
          polar = T)
jjVolcano(diffData = all.markers2,
          tile.col = corrplot::COL2('PuOr', 15)[1:12],
          size  = 3.5,
          fontface = 'italic',
          legend.position = c(0.8,0.2),
          flip = T)
dev.off()

genes <- c(
  "SPARC", "APOE", "COL1A2", "COL1A1", "HLA-DPA1", "SULF1", "TIMP1",
  "COL4A1", "BGN", "COL3A1", "COL5A1", "LYZ", "C1QB", "POSTN", "C1QC", "COL6A3",
  "CXCL9", "CCL19", "HSPA1B"
)
DotPlot(ST,features = genes,group.by = 'celltype')+RotatedAxis()
DotPlot(ST,features = genes)+RotatedAxis()
SCpubr::do_BeeSwarmPlot(sample = ST, 
                        feature_to_rank = "POSTN", 
                        continuous_feature = TRUE)
Idents(ST)<-'region'
p2 <- SCpubr::do_DotPlot(sample = ST,
                         features = genes, 
                         plot.title = "DEGs")

p2
library(patchwork)
gene_list <- genes
plot_list <- map(gene_list, function(gene) {
  SCpubr::do_RidgePlot(
    sample = ST,
    feature = gene,
    continuous_scale = TRUE,
    compute_quantiles = TRUE,
    color_by_probabilities = TRUE
  ) + ggplot2::ggtitle(gene)  
})
combined_plot <- wrap_plots(plot_list, ncol = 2)
print(combined_plot)
ggsave('combined.plot2.pdf',height = 40,width = 10)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
markers <- markers %>%
  mutate(Difference = pct.1 - pct.2) %>%
  rownames_to_column(var = "gene_name")  
library(ggplot2)
library(ggrepel)
markers <- markers %>%
  mutate(
    Difference = pct.1 - pct.2,
    sig_star = ifelse((avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05) |
                        (avg_log2FC <= -1 & Difference <= -0.1 & p_val_adj <= 0.05), TRUE, FALSE)
  )
ggplot(markers, aes(x = Difference, y = avg_log2FC)) +
  geom_point(data = subset(markers, !sig_star), aes(color = change), size = 0.1, shape = 16) +
  geom_point(data = subset(markers, sig_star), aes(color = change), size = 2.5, shape = 15) +
  scale_color_manual(
    name = 'change',
    labels = c(
      paste0("down (", table(markers$change)[[1]], ")"),
      "ns",
      paste0("up (", table(markers$change)[[3]], ")")
    ),
    values = c("Up" = "#5000B8", "Down" = "#008080", "NotSig" = "grey")
  ) +
  geom_label_repel(
    data = subset(markers, avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05),
    aes(label = gene),
    color = "black",
    segment.colour = "black",
    label.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    max.overlaps = 20
  ) +
  geom_label_repel(
    data = subset(markers, avg_log2FC <= -1 & Difference <= -0.1 & p_val_adj <= 0.05),
    aes(label = gene),
    color = "black",
    segment.colour = "black",
    label.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    max.overlaps = 20
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "△ Percentage difference", y = "Log-Fold change") +
  theme_bw()
ggsave("volcanoplot_new.pdf", width = 5, height = 3)

