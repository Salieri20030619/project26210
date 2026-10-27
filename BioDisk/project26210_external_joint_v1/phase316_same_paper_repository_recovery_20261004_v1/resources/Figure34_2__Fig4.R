#Fig4a
#CNV profile
CNV_clonal_region<-readRDS("/Path/data/CNV_clonal_region.RDS")
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
clone_order<-gsub(" ","_",clone_order)
pheatmap::pheatmap(CNV_clonal_region[c("chr20q","chr20p","chr7p","chr7q","chr13q","chr5p","chr9q","chr8q","chr12p","chr22q","chr8p_2","chr12q_1","chr12q_2",
                                       "chr8p_1","chr1p_2"),clone_order],cluster_rows = F,cluster_cols = F,border_color = "black",
                   color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),gaps_col = c(3),cellwidth = 10,cellheight = 10,fontsize = 8)
#Odds ratio calculating
CNV_clonal_region<-CNV_clonal_region[rowSums(CNV_clonal_region == "0") != 13,]
CNV_clonal_region<-CNV_clonal_region[c("chr20q","chr20p","chr7p","chr7q","chr13q","chr5p","chr9q","chr8q","chr12p","chr22q","chr8p_2","chr12q_1","chr12q_2",
                                       "chr8p_1","chr1p_2"),clone_order]

group1 <- abs(CNV_clonal_region[, 1:3])
group2 <- abs(CNV_clonal_region[, 4:ncol(CNV_clonal_region)])
ncol(group2)

OR_df<-matrix(nrow=15,ncol=4)
rownames(OR_df)<-rownames(CNV_clonal_region)
colnames(OR_df) <- c("OR","lower","upper","p")
program <- c('CNV', 'WT')
outcome <- c('Low', 'High')
library(epitools)
for (i in 0:14){
  i=i+1
  b<-sum(group1[i,])
  c<-sum(group2[i,])
  data <- matrix(c(c, b, 10-c, 3-b), nrow=2, ncol=2, byrow=TRUE)
  dimnames(data) <- list('Program'=program, 'Outcome'=outcome)
  OR_df[i,1:4]<-c(oddsratio(data,method = "small")$measure[2],oddsratio(data,method = "small")$measure[4],
                  oddsratio(data,method = "small")$measure[6],oddsratio(data,method = "small")$p.value[6])
}

OR_df<-as.data.frame(OR_df)
#OR_df=OR_df[order(OR_df[,1],decreasing = F),]
OR_df$CNV<-rownames(OR_df)
#OR_df<-OR_df[c(1:4,6,7,9,10,11,5,8,12,13,15,14),]
OR_df$CNV <- factor(OR_df$CNV,levels = rev(OR_df$CNV))
library(ggplot2)
#OR_df$OR = OR_df$OR + 0.01
theme_set(theme_cowplot())

colnames(OR_df)[5]<-"CNV"
ggplot(OR_df, aes(x = OR, y = CNV)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = `upper`, xmin = `lower`), size = 0.6,
                 height = 0.4,color="black") +
  geom_point(size = 2,color="black") +
  ylab("") + xlab("Odds ratio")+coord_trans(x = 'log10')+
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(0.1,1,10,100))


#Fig4b
library("devtools")
library("BiocManager")
library("infercnv")
library("tidyverse")
library("Seurat")
library("phylogram")
library("ape")
library("hdf5r")
library("SpatialInferCNV")
library(ape)
setwd("/data/huangzr/Spatial/Kras/version5/Region_devide_2")
# Create the inferCNV object from the 3 files, without the reference group and mitochondria
NC_CRC4_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./CRC4_count_matrix.tsv",
                                                 gene_order_file="./siCNV_GeneOrderFile.tsv",
                                                 annotations_file="./Region_CRC4.tsv",
                                                 delim="\t",
                                                 ref_group_names="Remaining",
                                                 chr_exclude = c("chrM","chrX","chrY"))
saveRDS(NC_CRC4_infCNV,file="./CRC4_infCNV_clone_group_2.RDS")

#Run inferCNV
NC_CRC4_infCNV_res_spot = infercnv::run(NC_CRC4_infCNV,
                                        cutoff=0.1, #(see infercnv::run documentation)
                                        out_dir="./SPCNVoutputs_CRC4_clone_group_2",
                                        cluster_by_groups=TRUE, #unsupervised analysis
                                        HMM = TRUE,k_obs_groups = 3,
                                        denoise=TRUE) #denoising applies noise reduction for the plot
saveRDS(NC_CRC4_infCNV_res_spot,file="./CRC4_infCNV_clone_group_2_res_cluster_level_3.RDS")

library("devtools")
library("BiocManager")
library("infercnv")
library("tidyverse")
library("Seurat")
library("phylogram")
library("ape")
library("hdf5r")
library("SpatialInferCNV")
library(ape)
library(pacman)
p_unload(Seurat)
p_load(Seurat)
CRC4_st<-readRDS("/Path/data/CRC4_region.RDS")
#Use read.dendrogram() to import the dendogram file
CRC4_st_for_clustering <- read.dendrogram(file = "/Path/data/infercnv.21_denoised.observations_dendrogram.txt")

#Convert to a phylo object using as.phylo()
CRC4_st_for_clustering_phylo <- as.phylo(CRC4_st_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(CRC4_st_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes
png("CRC4_st_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(CRC4_st_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:CRC4_st_for_clustering_phylo$Nnode,node=1:CRC4_st_for_clustering_phylo$Nnode+Ntip(CRC4_st_for_clustering_phylo))
dev.off()



#Selecting Node 1
Node1 <- SelectingSubTreeData(my.subtrees, 1988)
Node5 <- SelectingSubTreeData(my.subtrees, 993)

#Selecting Node 2
Node2 <- SelectingSubTreeData(my.subtrees, 816)

#Selecting Node 3
Node3 <- SelectingSubTreeData(my.subtrees, 2)

#Selecting Node 3
Node4 <- SelectingSubTreeData(my.subtrees, 2086)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2,Node3,Node4,Node5)
#MergedNodes <- Node1

#Note, since the "CRC4itional spot" for Clone 3 did not have a node, identify it by joining the original annotations to the MergedNodes dataframe
#MergedNodes <- full_join(MergedNodes, Celltype_Anno_CRC4_st)


#Then drop the Histology column
#MergedNodes <- MergedNodes[,1:2]

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_993", "Clone1",
                           ifelse(MergedNodes$Node == "Node_816", "Clone2",
                                  ifelse(MergedNodes$Node == "Node_2", "Clone3",
                                         ifelse(MergedNodes$Node == "Node_2086", "Remaining",
                                                ifelse(MergedNodes$Node == "Node_1988", "Clone1",MergedNodes$Node)))))

#Show the first part of the finalized dataframe for this step
head(MergedNodes)


#Copy the MergedNodes dataframe to a new dataframe called "ForLoupeBrowser"
ForLoupeBrowser <- MergedNodes

#Replace the "." in the current barcode dataframe to "-" (required for LoupeBrowser)
#ForLoupeBrowser$Barcode <- gsub("\\.", "\\-", ForLoupeBrowser$Barcode)
#ForLoupeBrowser$Barcode <- trimws(substr(ForLoupeBrowser$Barcode, 8, 100))
#save Loupe object
write.csv(ForLoupeBrowser, "CRC4_st_Clones_ForLoupeBrowser.csv", row.names = FALSE)
#CRC4_st@meta.data<-CRC4_st@meta.data
CRC4_st@meta.data$Barcode<-rownames(CRC4_st@meta.data)
CRC4_st@meta.data<-merge(CRC4_st@meta.data,ForLoupeBrowser,by="Barcode",all = T)
rownames(CRC4_st@meta.data) <- CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
rownames(CRC4_st@meta.data) <- CRC4_st@meta.data$Barcode
SpatialPlot(CRC4_st,group.by = c("Node"))
library(Seurat)
#CRC4_st$Node.y[CRC4_st$Node.y %in% c("Clone4","Clone5")] <- "Clone4"
#CRC4_st$Region[CRC4_st$Node %in% c("Clone1","Clone2")] <- "Clone1"
CRC4_st$Region[CRC4_st$Node %in% c("Clone1")] <- "Clone1"
CRC4_st$Region[CRC4_st$Node %in% c("Clone2")] <- "Clone2"
CRC4_st$Region[CRC4_st$Node %in% c("Clone3")] <- "Clone3"
CRC4_st$Region[CRC4_st$Node %in% c("Remaining")] <- "Remaining"
#CRC4_st$Region[CRC4_st$Node %in% c("Clone4")] <- "Clone3"
pdf("/Path/data/clone_distribution_CRC4_st.pdf",width = 6,height = 5)
SpatialPlot(CRC4_st,group.by = c("Region"),cols = c("#F6B272","#96C2D6","#AA2E2E","white"),image.alpha = 0)
dev.off()
saveRDS(CRC4_st,file = "/Path/data/CRC4_st_clone.RDS")

#Fig4d
#signature calculation

TCGA_COAD_origin_tumor<-readRDS("/Path/data/TCGA_COAD_origin_tumor.rds")
sc_marker<-readRDS("/Path/data/Seurat_cluster_markers.rds")
top_genes <- sc_marker %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

for (cluster in clusters) {
  TCGA_COAD_origin_tumor[[paste0(cluster, "_signature")]] <- rowMeans(
    as.data.frame(scale(TCGA_COAD_origin_tumor[, colnames(TCGA_COAD_origin_tumor) %in% top_genes$gene[top_genes$cluster == cluster]])) %>%dplyr::select(where(~ !any(is.na(.))))
  )
}


# split gain and loss CNV
adj_cols <- grep("_adj", colnames(TCGA_COAD_origin_tumor))
adj_CNV<-TCGA_COAD_origin_tumor[,adj_cols]
mat_gain <- apply(adj_CNV, 2, function(x) ifelse(x > 2, x, NA))
mat_loss <- apply(adj_CNV, 2, function(x) ifelse(x < 2, x, NA))

mat_new <- cbind(mat_gain, mat_loss)
colnames(mat_new) <- c(paste0(colnames(adj_CNV), "_gain"), paste0(colnames(adj_CNV), "_loss"))

TCGA_COAD_origin_tumor<-cbind(TCGA_COAD_origin_tumor,mat_new)


#TCGA matrix
colnames(TCGA_COAD_origin_tumor)[adj_cols]
adj_cols <- grep("_adj_", colnames(TCGA_COAD_origin_tumor))
signature_cols <- grep("signature", colnames(TCGA_COAD_origin_tumor))

cor_matrix <- matrix(0, nrow = length(adj_cols), ncol = length(signature_cols))
rownames(cor_matrix) <- colnames(TCGA_COAD_origin_tumor)[adj_cols]
colnames(cor_matrix) <- colnames(TCGA_COAD_origin_tumor)[signature_cols]

for (i in 1:length(adj_cols)) {
  for (j in 1:length(signature_cols)) {
    cor_matrix[i, j] <- cor(TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), adj_cols[i]], TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), signature_cols[j]], method = "pearson")
  }
}

cor_matrix<-cor_matrix[,4:ncol(cor_matrix)]
colnames(cor_matrix)<-gsub("_signature","",colnames(cor_matrix))
rownames(cor_matrix)<-gsub("_adj","",rownames(cor_matrix))
bk = unique(c(seq(-0.4,0.4, length=100)))
cor_matrix[abs(cor_matrix) < 0.2] = 0
cor_matrix[grep("loss",rownames(cor_matrix)),]<- (-cor_matrix[grep("loss",rownames(cor_matrix)),])
p1<-pheatmap::pheatmap(cor_matrix[c("chr20q_gain","chr7p_gain","chr20p_gain","chr7q_gain","chr17p_loss","chr8q_gain","chr8p_gain","chr1p_2_loss"),1:29],
                       color = colorRampPalette(c("#3384BB","#A6CEE3","#E4F0F6","white","white","white","white","white","white","#FEECD4","#FDBF6F","#FE7F02"))(100),cluster_rows = F,
                       cluster_cols = F,breaks = bk,border_color = "#262626",angle_col = 90,cellwidth = 15,cellheight = 15,
                       gaps_col = c(10,20,25,27,29),gaps_row = c(4,5,7))


library(reshape2)
CNV_Frequency<-readRDS("/Path/data/CNV_Frequency.RDS")
CNV_Frequency<-data.frame(loss_frequency=colSums(CNV_module_filtered_1 < (-0.35))/nrow(CNV_module_filtered_1),gain_frequency=colSums(CNV_module_filtered_1 > (0.35))/nrow(CNV_module_filtered_1))
CRC4_CNV_freq<-CNV_Frequency[c("chr20q","chr7p","chr7q","chr20p","chr17p","chr17q","chr8p","chr8q","chr1p_2","chr14q"),]
CRC4_CNV_freq$region<-rownames(CRC4_CNV_freq)
CRC4_CNV_freq<-melt(CRC4_CNV_freq)
CRC4_CNV_freq$region<-factor(CRC4_CNV_freq$region,levels=rev(c("chr20q","chr20p","chr7p","chr7q","chr17p","chr17q","chr8q","chr8p","chr14q","chr1p_2")))
theme_set(theme_cowplot())

CRC4_CNV_freq_2<-CRC4_CNV_freq[c(11,12,14,13,5,18,17,9),]
CRC4_CNV_freq$region<-factor(CRC4_CNV_freq$region,levels=rev(c("chr20q","chr7p","chr20p","chr7q","chr17p","chr17q","chr8q","chr8p","chr14q","chr1p_2")))
ggplot(CRC4_CNV_freq_2,aes(x=value,y=region,fill=variable))+geom_bar(stat = "identity",width = 0.7)+scale_fill_manual(values = c("#97C3D8","#F6B373"))
p2<-TCGA_COAD_origin_tumor$chr20p_status

#Fig4e
# CNV co-occurance
copy_number_CNV_module_sep<-readRDS(copy_number_CNV_module_sep,"/Path/data/copy_number_CNV_module_sep.RDS")
CNV_sep_2<-copy_number_CNV_module_sep[,1:27]
CNV_sep_2_01_matrix<-(CNV_sep_2 > 2.549121 | CNV_sep_2 <1.569168)
CNV_sep_2_01_matrix[is.na(CNV_sep_2_01_matrix)]<-FALSE
CNV_sep_2_filtered<-CNV_sep_2[,colSums(CNV_sep_2_01_matrix)/466 > 0.05]
CNV_cor_matrix_2<-cor(CNV_sep_2)
cor_mat_2<-cor(CNV_cor_matrix_2,use = "pairwise.complete.obs")
colnames(CNV_cor_matrix_2)<-gsub("_adj","",colnames(CNV_cor_matrix_2))
rownames(CNV_cor_matrix_2)<-gsub("_adj","",rownames(CNV_cor_matrix_2))
colnames(cor_mat_2)<-gsub("_adj","",colnames(cor_mat_2))
rownames(cor_mat_2)<-gsub("_adj","",rownames(cor_mat_2))

pval_mat <- matrix(NA, ncol(CNV_sep_2), ncol(CNV_sep_2))
for(i in 1:ncol(CNV_sep_2)) {
  for(j in 1:ncol(CNV_sep_2)) {
    temp <- cor.test(CNV_sep_2[,i], CNV_sep_2[,j])
    pval_mat[i,j] <- temp$p.value
  }
}
sig_mat <- ifelse(pval_mat < 0.001, "***", ifelse(pval_mat < 0.01, "**", ifelse(pval_mat < 0.05, "*", "")))

rownames(sig_mat) = rownames(CNV_cor_matrix_2)
colnames(sig_mat) = colnames(CNV_cor_matrix_2)
sig_mat<-sig_mat[c("chr20p","chr20q","chr7p","chr7q","chr17p"),c("chr20p","chr20q","chr7p","chr7q","chr17p")]

library(pheatmap)
pheatmap(cor_mat_2[c("chr20p","chr20q","chr7p","chr7q","chr17p"),
                          c("chr20p","chr20q","chr7p","chr7q","chr17p")],
         color = colorRampPalette(c("#3384BB","#C3DAEB","white","white","white","white","#FED9B4","#FE7F02"))(100),
         border_color = "black",display_numbers = sig_mat,cluster_rows = F,cluster_cols = F)

#Fig4f
p1<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=log2(CD8A),fill=chr20q_20p_status,color=chr20q_20p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p2<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=log2(CD8A),fill=chr20q_7p_status,color=chr20q_7p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p3<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=log2(CD8A),fill=chr20q_7q_status,color=chr20q_7q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p4<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=log2(CD8A),fill=chr20q_17p_status,color=chr20q_17p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p1+p2+p3+p4


#Fig4g
CRC4_st <- readRDS("/Path/data/CRC4_st.rds")
library(pacman)
library(Seurat)
library(patchwork)
library(ggrepel)
library(clusterProfiler)

p_unload(Seurat)
p_load(Seurat)
DefaultAssay(CRC4_st) <- "Spatial"
CRC4_st_tumor<-CRC4_st[,CRC4_st$Region %in% c("Clone1","Clone2","Clone3")]

CRC4_st_tumor <- readRDS("/Path/data/CRC4_st_tumor.rds")
Idents(CRC4_st_tumor)<-CRC4_st_tumor@meta.data$Region
diff_ST_subclone<-FindAllMarkers(CRC4_st_tumor,assay = "Spatial",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,
                                 min.cells.feature = 0)
diff_ST_subclone$p_val_adj <- -log10(diff_ST_subclone$p_val_adj)
diff_ST_subclone$change<-"Remain"
diff_ST_subclone$change[diff_ST_subclone$avg_logFC > 0.1 & diff_ST_subclone$p_val_adj >1.30103]<-"Up"
diff_ST_subclone$change[diff_ST_subclone$avg_logFC < (-0.1) & diff_ST_subclone$p_val_adj >1.30103]<-"Down"
diff_ST_subclone$change <- factor(diff_ST_subclone$change,levels = c("Up","Down","Remain"))
label<-c("IRF7","IRF1","ISG15","IFI27","IFI16","CXCL9","CXCL10","CXCL11","B2M","PSMB9","PSMB8","TAP1","TAP2",
         "HLA-A","HLA-B","STAT1","OLFM4","FCGBP","MUC2","RPL37","RPL31","RPS5")
label_data_down<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3" & diff_ST_subclone$gene %in% label & diff_ST_subclone$change == "Down",]
label_data_up<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3" & diff_ST_subclone$gene %in% label & diff_ST_subclone$change == "Up",]
theme_set(theme_cowplot())

ggplot(diff_ST_subclone[diff_ST_subclone$cluster == "Clone3",],aes(x=avg_logFC,y=p_val_adj,color=change))+geom_point()+
  scale_color_manual(values = c("#FE7F00","#1F78B4","lightgrey"))+coord_cartesian(xlim=c(-2,2))+
  labs(x="log2(Clone3 vs Other Clones)",y="-log10(adj.P)",title = "Differential Expression Genes in Clone3")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "top",legend.text = element_text(hjust = 0.5))+
  geom_vline(xintercept = -0.1,lty=2)+geom_vline(xintercept = 0.1,lty=2)+
  geom_hline(yintercept = 2,lty=2)+geom_text_repel(data=label_data_down,aes(color=change, label=gene),
                                                   point.padding = unit(0.2, "lines"), box.padding = unit(0.4, "lines"),max.overlaps = 12,
                                                   nudge_y = 0.1,ylim = c(5,300),xlim = c(-2,-1))+
  geom_text_repel(data=label_data_up,aes(color=change, label=gene),
                  point.padding = unit(0.2, "lines"), box.padding = unit(0.4, "lines"),max.overlaps = 10,
                  nudge_y = 0.1,ylim = c(5,300),xlim=c(1,2))+guides(color=F)

Sub3_up<-diff_ST_subclone[(diff_ST_subclone$change == "Up") & (diff_ST_subclone$cluster == "Clone3"),]
Sub3_down<-diff_ST_subclone[(diff_ST_subclone$change == "Down") & (diff_ST_subclone$cluster == "Clone3"),]
#Sub3_up<-readRDS("/data/huangzr/Spatial/Kras/version5/Example_data/Sub3_up_gene.RDS")
#Sub3_down<-readRDS("/data/huangzr/Spatial/Kras/version5/Example_data/Sub3_down_gene.RDS")
up_GO<-enrichGO(Sub3_up,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
down_GO<-enrichGO(Sub3_down,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
up_GO@result$pvalue=-log10(up_GO@result$pvalue)
down_GO@result$pvalue=-log10(down_GO@result$pvalue)
selected_GO_up<-up_GO@result[c(9,8,6,1),]
selected_GO_down<-down_GO@result[c(1133,169,105,2),]
selected_GO_down$Description<-factor(selected_GO_down$Description,levels=selected_GO_down$Description)
selected_GO_up$Description<-factor(selected_GO_up$Description,levels=selected_GO_up$Description)
theme_set(theme_cowplot())
p1<-ggplot(selected_GO_down,aes(x=Description,y=pvalue,fill=Description))+
  geom_bar(stat = "identity",width=0.8)+coord_flip()+guides(fill=F)+
  labs(x="",y=expression(-log[10]~italic(P)),title = "Clone3 Down")+theme(axis.text.x =element_text(size=16),aspect.ratio = 1)+
  ylim(0,20)+scale_fill_manual(values = c(rep("#A5CDE2",6)))+
  theme(plot.title = element_text(hjust = 0.5))
p2<-ggplot(selected_GO_up,aes(x=Description,y=pvalue,fill=Description))+
  geom_bar(stat = "identity",width=0.8)+coord_flip()+guides(fill=F)+
  labs(x="",y=expression(-log[10]~italic(P)),title = "Clone3 Up")+theme(axis.text.x =element_text(size=16),aspect.ratio = 1)+
  ylim(0,70)+scale_fill_manual(values = c(rep("#FDBE6D",6)))+
  theme(plot.title = element_text(hjust = 0.5))
p1+p2+plot_layout(ncol=1)

#Fig4i
#diff TE analysis
library(Seurat)
library(cowplot)
p_unload(Seurat)
p_load(Seurat)

CRC4_TE <- readRDS("/Path/data/CRC4_st_with_TE.RDS")
DefaultAssay(CRC4_TE)<-"TE"
CRC4_TE@meta.data<-CRC4_TE@meta.data[colnames(CRC4_TE),]
table(colnames(CRC4_TE) == rownames(CRC4_TE@meta.data))

CRC4_TE_Tumor<-CRC4_TE[,CRC4_TE$Region %in% c("Clone1","Clone2","Clone3")]
Idents(CRC4_TE_Tumor)<-CRC4_TE_Tumor$Region
diff_TE<-FindAllMarkers(CRC4_TE_Tumor,assay = "TE",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,only.pos = F)
repeats_info<-read.csv("/Path/data/repeats_family_info.csv",header = F)
CRC_specific_TE<-read.csv("/Path/data/GCP_COAD_tumor_vs_normal_Diff_TE.csv")
CRC_specific_TE$change_tumor<-"Remain"
CRC_specific_TE$change_tumor[CRC_specific_TE$logFC>0] <- "Up"
CRC_specific_TE$change_tumor[CRC_specific_TE$logFC<0] <- "Down"
CRC_specific_TE<-CRC_specific_TE$X[CRC_specific_TE$change_tumor == "Up"]
diff_CRC_TE<-diff_TE[diff_TE$gene %in% CRC_specific_TE,]
Clone3_TE<-diff_TE[diff_TE$cluster == "Clone3",]
Clone3_TE$p_val_adj <- -log10(Clone3_TE$p_val_adj)
Clone3_TE$change<-"Remain"
Clone3_TE$change[Clone3_TE$p_val_adj > 2 & Clone3_TE$avg_logFC > 0.25] = "Up"
Clone3_TE$change[Clone3_TE$p_val_adj > 2 & Clone3_TE$avg_logFC < (-0.25)] = "Down"
Clone3_TE$group<-"Remain"
Clone3_TE$group[Clone3_TE$gene %in% c("LTR7","LTR7Y","LTR7B","HERVH-int")]<-"HERVH"

theme_set(theme_cowplot())
library(ggrepel)
for_label <-  Clone3_TE[(Clone3_TE$change != "Remain") & (Clone3_TE$gene %in% CRC_specific_TE),]
for_label_2 <-  Clone3_TE[(Clone3_TE$group == "HERVH") & (Clone3_TE$gene %in% CRC_specific_TE),]
theme_set(theme_cowplot())

ggplot(Clone3_TE[Clone3_TE$gene %in% CRC_specific_TE,],aes(x=avg_logFC,y=p_val_adj,fill=change,color=group))+geom_point(shape=21,stroke=0.5,size=1.5)+
  scale_color_manual(values = c("black",alpha("white",0),rep(alpha("white",0),9)))+
  scale_fill_manual(values = c("#1E72A9","lightgrey","#D1332D"))+
  labs(x="log2(Clone3 vs Other Subclones)",y="-log2(adj.P)",title = "Tumor Specific TE in CRC")+
  theme(plot.title = element_text(size=12,hjust=0.5))+geom_vline(xintercept = 0.25,lty=2)+geom_vline(xintercept = -0.25,lty=2)+
  geom_hline(yintercept = 2,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+
  geom_text_repel(data = for_label_2[for_label_2$avg_logFC<0,],aes(label = gene),alpha = 1,color="#1F78B4",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(-0.3,-1.2),size=4)+
  geom_text_repel(data = for_label[for_label$avg_logFC>0,],aes(label = gene),alpha = 1,color="#D1332D",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(0.3,1.2),size=4)

#Fig4j
library(Seurat)
library(cowplot)
p_unload(Seurat)
p_load(Seurat)
DefaultAssay(CRC4_TE)<-"TE"
CRC4_TE@meta.data<-CRC4_TE@meta.data[colnames(CRC4_TE),]
table(colnames(CRC4_TE) == rownames(CRC4_TE@meta.data))
CRC4_TE@meta.data$HERVH<-colMeans(CRC4_TE@assays$TE@data[c("LTR7Y","LTR7","LTR7B","HERVH-int"),])
DefaultAssay(CRC4_TE) <- "Spatial"
p1<-SpatialPlot(CRC4_TE,features = "HERVH",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HERV-H",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))

Idents(CRC4_TE)<-"Region"
HERVH_expression<-CRC4_TE@meta.data[,c("Region","HERVH")]
p2<-ggplot(HERVH_expression[HERVH_expression$Region %in% c("Clone1","Clone2","Clone3"),],aes(x=Region,y=HERVH,fill=Region,color=Region))+
  geom_violin()+scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+
  scale_color_manual(values = c("#378B39","#1C6797","#B62D28"))+
  labs(x="Clones",y="Expression",title = "HERV-H")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p1+p2

#Fig4k
TCGA_COAD_total_with_20q<-readRDS("/data/huangzr/Spatial/Kras/version5/Example_data/TCGA_COAD_total_with_20q.RDS")
TCGA_COAD_total_with_20q$chr20q_status.y<-factor(TCGA_COAD_total_with_20q$chr20q_status.y,levels=c("Normal","WT","Gain"))
ggplot(TCGA_COAD_total_with_20q[!is.na(TCGA_COAD_total_with_20q$chr20q_status.y),],aes(x=chr20q_status.y,y=log2(HERVH_expression),color=chr20q_status.y))+
  geom_boxplot(size=1)+scale_color_manual(values = c("#BFBDBD","#1C6797","#DD7E2D"))+labs(x="",y="log2(RPM)",title="HERV-H")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust=0.5))+geom_signif(comparisons = list(c("Normal", "WT"),c("Normal", "Gain"),c("WT", "Gain")),y_position=c(8, 9, 8))

#Fig4l
# ICB response
library(survival)
library(ggplot2)
library(survminer)
MSK_cohort <- readRDS("/Path/data/MSK_cohort.rds")
# ICB survival CNV
MSK_cohort$ID<-gsub("-",".",MSK_cohort$Sample.ID)
ICB_CNV_all<-read.csv("/Path/data/all_data_by_genes.txt",sep = "\t")
info<-ICB_CNV_all[,1:3]
ICB_CNV_all<-ICB_CNV_all[,colnames(ICB_CNV_all) %in% MSK_cohort$ID]
ICB_CNV_all<-cbind(info,ICB_CNV_all)
ICB_CNV_all$arm <- sub("(p|q).*", "\\1", ICB_CNV_all$Cytoband)

ICB_CNV_all_df<-data.frame(chr20q_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "20q",4:74]),
                           chr20p_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "20p",4:74]),
                           chr17p_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "17p",4:74]),
                           chr7p_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "7p",4:74]),
                           chr7q_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "7q",4:74]))

ICB_CNV_all_df_new <- apply(ICB_CNV_all_df, 2, function(x) {
  # 使用ifelse函数根据条件选择值
  ifelse(x > 0.35, "gain", ifelse(x < -0.35, "loss", "WT"))
})
ICB_CNV_all_df_new<-as.data.frame(ICB_CNV_all_df_new)
ICB_CNV_all_df_new$ID <- rownames(ICB_CNV_all_df_new)

MSK_cohort_with_CNV_status<-merge(MSK_cohort,ICB_CNV_all_df_new,by="ID")
MSK_cohort_with_CNV_status$CNV_status<-paste(MSK_cohort_with_CNV_status[,52],MSK_cohort_with_CNV_status[,53],MSK_cohort_with_CNV_status[,54],
                                             MSK_cohort_with_CNV_status[,55],MSK_cohort_with_CNV_status[,56])
#survival
MSK_cohort_with_CNV_status$status<-1
MSK_cohort_with_CNV_status$status[MSK_cohort_with_CNV_status$Overall.Survival.Status.x == "1:DECEASED"]<-2
colnames(MSK_cohort_with_CNV_status)[15]<-"time"
#msk_ICB_cohort
fit <- survfit(Surv(time,status) ~ chr20q_CNV,  # 创建生存对象
               data = MSK_cohort_with_CNV_status)
p1<-ggsurvplot(fit, data = MSK_cohort_with_CNV_status,pval = TRUE,palette = c("#FF7F00","#226C9E"))+
  labs(x="Time post first dose of ICB (Months)",title = "CRC ICB cohort (MSK-IMPACT)")
p1

#TCGA survival
clinical <- read.table('/Path/data/COAD_clinical_filter.txt',sep="\t",header=T)
clinical$vital_status <- ifelse(clinical$vital_status=="Alive", 0, 1)
clinical$days_to_death <- gsub('--',-1,clinical$days_to_death)
clinical$days_to_last_follow_up <- gsub('--',-1,clinical$days_to_last_follow_up)
clinical$time <- ifelse(clinical$vital_status==0, clinical$days_to_last_follow_up, clinical$days_to_death)
clinical_filter <- clinical[-which(clinical$time < 10),]
clinical$time <- as.numeric(clinical$time)/365
clinical_filter$time <- as.numeric(clinical_filter$time)/365
colnames(clinical_filter)[1]<-"case"
clinical_COAD_CNV<-merge(clinical_filter,TCGA_COAD_origin_tumor,by="case")

TCGA_COAD_origin_tumor_surv <- survfit(Surv(clinical_COAD_CNV$time,event=clinical_COAD_CNV$vital_status)~chr20q_status,
                                       data=clinical_COAD_CNV)
p2<-ggsurvplot(TCGA_COAD_origin_tumor_surv, data = clinical_COAD_CNV,pval = TRUE,palette = c("#FF7F00","#226C9E"))+
  labs(x="Time (Months)",title = "CRC non-ICB cohort (TCGA-COAD)")
p2

