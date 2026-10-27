#################Figure6-7
meta =sce_Resistant@meta.data # a dataframe with rownames containing cell mata data
data_input <- as.matrix(sce_Resistant@assays$RNA@data)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")

#-------Add cell information into meta slot of the object (Optional)
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

cellchat@data.signaling


CellChatDB.use <- CellChatDB# use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,file = "resistant/cellchat_resistant.rds")
meta =sce_Response@meta.data # a dataframe with rownames containing cell mata data
data_input <- as.matrix(sce_Response@assays$RNA@data)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")
#-------Add cell information into meta slot of the object (Optional)
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
cellchat@data.signaling
CellChatDB.use <- CellChatDB# use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,file = "expose/cellchat_expose.rds")

library(CellChat)
cellchat_resistant <- readRDS("resistant/cellchat_resistant.rds")
cellchat_expose <- readRDS("expose/cellchat_expose.rds")
object.list <- list(Resistant=cellchat_resistant,Expose=cellchat_expose)
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c("Resistant","Expose"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c("Resistant","Expose"), measure = "weight")
gg1 + gg2
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c("Resistant","Expose"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c("Resistant","Expose"))
mat <- cellchat@net$Resistant$weight
groupSize <- as.numeric(table(cellchat_resistant@idents)) 
pdf("1.netvisual.cricle.single.Resistant.pdf",height = 20,width =20)
par(mfrow = c(6,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

mat <- cellchat@net$Expose$weight
groupSize <- as.numeric(table(cellchat_expose@idents)) 
pdf("1.netvisual.cricle.single.Expose.pdf",height = 20,width =20)
par(mfrow = c(6,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
gg1 <- netVisual_heatmap(cellchat,comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2))
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,2))
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 25)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 25)
pdf("netanalysis.heatmap.pdf",height = 15,width =8)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
cellchat_resistant <- computeNetSimilarity(cellchat_resistant, type = "functional")
cellchat_resistant <- netEmbedding(cellchat_resistant, type = "functional")
cellchat_resistant <- netClustering(cellchat_resistant, type = "functional",do.parallel = FALSE)
# figure19
netVisual_embedding(cellchat_resistant, type = "functional", label.size = 3.5)
cellchat_resistant <- computeNetSimilarity(cellchat_resistant, type = "structural")
cellchat_resistant <- netEmbedding(cellchat_resistant, type = "structural")
cellchat_resistant <- netClustering(cellchat_resistant, type = "structural",do.parallel = FALSE)
#figure20
netVisual_embedding(cellchat_resistant, type = "structural", label.size = 3.5)
a<-as.data.frame(ht2@matrix)
pathways.show <- as.character(rownames(a))
vertex.receiver <- seq(1, 4)
group.cellType <- c("Macro", "Tcell", "Tcell", "Endo", "Tcell", "DC", "CAF", "Tumor")
names(group.cellType) <- c("APOE_TAM", "CD8_Tem", "CD8_Tex", "Endo", "FOXP3_Treg", "moDC", "POSTN_CAF", "Tumor")
for (pathway in pathways.show) {
  pdf(file = paste0(save_path, pathway, "_hierarchy.pdf"))
  netVisual_aggregate(cellchat_resistant, signaling = pathway, vertex.receiver = vertex.receiver)
  netVisual_aggregate(cellchat_resistant, signaling = pathway, layout = "chord")
  dev.off()
  pdf(file = paste0(save_path, pathway, "_chord.pdf"))
  netVisual_chord_cell(cellchat_resistant, signaling = pathway, group = group.cellType, title.name = paste0(pathway, " signaling network"))
  dev.off()
  pdf(file = paste0(save_path, pathway, "_circle.pdf"), height = 10, width = 10)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat_resistant, signaling = pathway, layout = "circle")
  dev.off()
  pdf(file = paste0(save_path, pathway, "_network_centrality.pdf"))
  netAnalysis_signalingRole_network(cellchat_resistant, signaling = pathway, width = 15, height = 2.5, font.size = 10)
  gg1 <- netAnalysis_signalingRole_scatter(cellchat_resistant, signaling = pathway)
  gg2 <- netAnalysis_signalingRole_scatter(cellchat_resistant, signaling = pathway)
  print(gg1 + gg2)
  dev.off()
  pdf(file = paste0(save_path, pathway, "_heatmap.pdf"))
  p1 = netVisual_heatmap(cellchat_resistant, signaling = pathway, color.heatmap = "Reds")
  p2 = netVisual_heatmap(cellchat_resistant, color.heatmap = "Reds")
  print(p1 + p2)
  dev.off()
  pdf(file = paste0(save_path, pathway, "_contribution.pdf"))
  netAnalysis_contribution(cellchat_resistant, signaling = pathway)
  dev.off()
  pairLR <- extractEnrichedLR(cellchat_resistant, signaling = pathway, geneLR.return = FALSE)
  pdf(file = paste0(save_path, pathway, "_individual.pdf"))
  netVisual_individual(cellchat_resistant, signaling = pathway, pairLR.use = pairLR[1,], layout = "circle")
  dev.off()
}
pathways.show <- c(
  "MHC-I", "MHC-II", "MIF", "CCL", "SPP1", "CD99", "COLLAGEN", "IL2", "CXCL",
  "GALECTIN", "LCK", "ANGPTL", "PECAM1", "MK", "FN1", "ICAM", "ITGB2",
  "LAMININ", "CD70", "IL16", "LT", "COMPLEMENT", "PTN", "VCAM", "IL4", "PARs",
  "VEGF", "THY1", "TGFb", "ESAM", "ANNEXIN", "EPHA", "IGF", "THBS", "CD40",
  "CD96", "LIGHT", "JAM", "BTLA", "ANGPT", "IFN-II", "VISFATIN", "ADGRE5", "NOTCH",
  "TIGIT", "SELL", "CSF", "CD34", "GRN", "APRIL", "CD45", "SELPLG",
  "EGF", "TWEAK", "PVR", "IL6", "IL10", "BAFF", "CHEMERIN", "EPHB",
  "GAS", "PTPRM", "CDH5", "CD46", "HGF", "CALCR", "PDGF", "CDH1", "CD22",
  "TRAIL", "SEMA7", "ALCAM", "CD6", "SEMA4", "OSM",
  "CD86", "L1CAM", "IL1", "FLT3", "MPZ", "FGF", "CD23",
  "NRG", "BMP", "OX40", "CDH", "SN", "EDN", "CD80",
  "LIFR", "FASLG", "NT", "CD48", "CX3C",
  "NPR2", "NGF", "VTN", "ACTIVIN",
  "OCLN", "CD39", "GP1BA", "MSTN"
)

IL
vertex.receiver <- seq(1, 4)
for (pathway in pathways.show) {
  pdf(file = paste0(save_path, pathway, "_hierarchy.pdf"))
  netVisual_aggregate(cellchat_expose, signaling = pathway, vertex.receiver = vertex.receiver)
  netVisual_aggregate(cellchat_expose, signaling = pathway, layout = "chord")
  dev.off()
  pdf(file = paste0(save_path, pathway, "_chord.pdf"))
  netVisual_chord_cell(cellchat_expose, signaling = pathway, group = group.cellType, title.name = paste0(pathway, " signaling network"))
  dev.off()
  pdf(file = paste0(save_path, pathway, "_circle.pdf"), height = 10, width = 10)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat_expose, signaling = pathway, layout = "circle")
  dev.off()
  pdf(file = paste0(save_path, pathway, "_network_centrality.pdf"))
  netAnalysis_signalingRole_network(cellchat_expose, signaling = pathway, width = 15, height = 2.5, font.size = 10)
  gg1 <- netAnalysis_signalingRole_scatter(cellchat_expose, signaling = pathway)
  gg2 <- netAnalysis_signalingRole_scatter(cellchat_expose, signaling = pathway)
  print(gg1 + gg2)
  dev.off()
  pdf(file = paste0(save_path, pathway, "_heatmap.pdf"))
  p1 = netVisual_heatmap(cellchat_expose, signaling = pathway, color.heatmap = "Reds")
  p2 = netVisual_heatmap(cellchat_expose, color.heatmap = "Reds")
  print(p1 + p2)
  dev.off()
  pdf(file = paste0(save_path, pathway, "_contribution.pdf"))
  netAnalysis_contribution(cellchat_expose, signaling = pathway)
  dev.off()
  pairLR <- extractEnrichedLR(cellchat_expose, signaling = pathway, geneLR.return = FALSE)
  pdf(file = paste0(save_path, pathway, "_individual.pdf"))
  netVisual_individual(cellchat_expose, signaling = pathway, pairLR.use = pairLR[1,], layout = "circle")
  dev.off()
}


library(nichenetr)


seurat_obj_receiver <- subset(sce, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = "Non-Response", 
                                  ident.2 = "Response",
                                  group.by = "Respponse",
                                  min.pct = 0.1) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>% drop_na()
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 



#=========================================================================================
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>% drop_na()

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
DotPlot(subset(sce, celltype %in% sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")



celltype_order <- levels(Idents(sce)) 
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = sce,
  condition_colname = "Respponse",
  condition_oi = "Non-Response",
  condition_reference = "Response",
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% column_to_rownames("gene") 
vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc


# compare rankings between the sender-agnostic and sender-focused approach
(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p1 = make_heatmap_ggplot(vis_ligand_aupr,"Prioritized ligands","Ligand activity",legend_title = "AUPR", color ="#EFBD37")+
  theme_bw()+
  theme(axis.text.x.top = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10,face ="bold",hjust = 1))+
  xlab("Ligand\nactivity")

sender=c("Tumor")

#p2-expression sender

p2 <- DotPlot(
  subset(sce, celltype %in% sender),
  features = rev(best_upstream_ligands),    # 
  cols = "RdYlBu",
  group.by = "Respponse"                     # 
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  theme_bw() +
  geom_point(shape = 1, aes(size = pct.exp), stroke = 0.5) +  # 
  theme(
    axis.text.x = element_text(colour = 'black', size = 9, angle = 90, hjust = 0),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank()
  ) +
  ylab("Ligand Expression\nin Sender")


p2
#p3-ligands-target heatmap
# Target gene plot
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.92) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p3 <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                          color = '#810f7c', legend_title = "Regulatory potential") +
  theme(axis.text.x.top = element_text(colour = 'black', face = "italic", angle = 90, hjust = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10,face ="bold"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


pp=p1+p2+p3+plot_layout(ncol = 3, widths  = c(2, length(sender)+7, ncol(vis_ligand_target)),guides = 'collect')
pp & theme(plot.margin = margin(1,2,2,0),
           legend.position = "bottom",
           legend.direction = 'horizontal',
           # legend.justification=c(0,1),
           legend.key.width=unit(0.5,"cm"),
           legend.key.height = unit(0.3, "cm"),
           legend.title.position = 'top')




pp <- p1 + p2 + p3 + 
  plot_layout(
    ncol = 3,
    widths = c(2, length(sender) + 7, ncol(vis_ligand_target) * 1.5),  
    guides = 'collect'
  )

pp & theme(
  plot.margin = margin(1, 2, 2.5, 0),
  legend.position = "bottom",
  legend.direction = 'horizontal',
  legend.key.width = unit(0.5, "cm"),
  legend.key.height = unit(0.3, "cm"),
  legend.title.position = 'top'
)






