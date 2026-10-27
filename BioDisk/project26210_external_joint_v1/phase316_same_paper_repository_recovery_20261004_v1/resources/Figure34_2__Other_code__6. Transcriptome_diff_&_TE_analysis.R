# Transcriptome analysis -> diff gene expression analysis
p_unload(Seurat)
p_load(Seurat)
DefaultAssay(CRC4_st) <- "Spatial"
CRC4_st_tumor<-CRC4_st[,CRC4_st$Region %in% c("Clone1","Clone2","Clone3")]


Idents(CRC4_st_tumor)<-CRC4_st_tumor@meta.data$Region
diff_ST_subclone<-FindAllMarkers(CRC4_st_tumor,assay = "Spatial",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,
                           min.cells.feature = 0)
table(diff_ST_subclone$change[diff_ST_subclone$cluster == "Clone3"])
diff_ST_subclone$change<-"Remain"
diff_ST_subclone$change[diff_ST_subclone$avg_logFC > 0.1 & diff_ST_subclone$p_val_adj >1.30103]<-"Up"
diff_ST_subclone$change[diff_ST_subclone$avg_logFC < (-0.1) & diff_ST_subclone$p_val_adj >1.30103]<-"Down"
diff_ST_subclone$change <- factor(diff_ST_subclone$change,levels = c("Up","Down","Remain"))

library(ggrepel)
label<-c("IRF7","IRF1","ISG15","IFI27","IFI16","CXCL9","CXCL10","CXCL11","B2M","PSMB9","PSMB8","TAP1","TAP2",
         "HLA-A","HLA-B","STAT1","OLFM4","RPL36A","RPS18","FCGBP","MUC2","RPL37","RPL31","RPS27","RPS5","RPS9","RPL38")
label_data_down<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3" & diff_ST_subclone$gene %in% label & diff_ST_subclone$change == "Down",]
label_data_up<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3" & diff_ST_subclone$gene %in% label & diff_ST_subclone$change == "Up",]

#diff_ST_subclone$p_val_adj <- -log10(diff_ST_subclone$p_val_adj)
library(ggrepel)
theme_set(theme_cowplot())

ggplot(diff_ST_subclone[diff_ST_subclone$cluster == "Clone3",],aes(x=avg_logFC,y=p_val_adj,color=change))+geom_point()+
  scale_color_manual(values = c("#FE7F00","#1F78B4","lightgrey"))+coord_cartesian(xlim=c(-2,2))+
  labs(x="log2(Clone3 vs Other Clone)",y="-log10(adj.P)",title = "")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "top",legend.text = element_text(hjust = 0.5))+
  geom_vline(xintercept = -0.1,lty=2)+geom_vline(xintercept = 0.1,lty=2)+
  geom_hline(yintercept = 2,lty=2)+geom_text_repel(data=label_data_down,aes(color=change, label=gene), 
                                                   point.padding = unit(0.2, "lines"), box.padding = unit(0.4, "lines"),max.overlaps = 12, 
                                                   nudge_y = 0.1,ylim = c(5,300),xlim = c(-2,-1))+
  geom_text_repel(data=label_data_up,aes(color=change, label=gene), 
                   point.padding = unit(0.2, "lines"), box.padding = unit(0.4, "lines"),max.overlaps = 10, 
                   nudge_y = 0.1,ylim = c(5,300),xlim=c(1,2))+guides(color=F)


#enrichment
View(diff_ST_subclone)
Clone3_up<-diff_ST_subclone$gene[diff_ST_subclone$cluster== "Clone3" & diff_ST_subclone$change == "Up"]
Clone3_down<-diff_ST_subclone$gene[diff_ST_subclone$cluster== "Clone3" & diff_ST_subclone$change == "Down"]
saveRDS(Clone3_up,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/Clone3_up_gene.RDS")
saveRDS(Clone3_down,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/Clone3_down_gene.RDS")
library(clusterProfiler)
up_GO<-enrichGO(Clone3_up,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
down_GO<-enrichGO(Clone3_down,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
View(up_GO@result)
View(down_GO@result)
up_GO@result$pvalue=-log10(up_GO@result$pvalue)
down_GO@result$pvalue=-log10(down_GO@result$pvalue)
selected_GO_up<-up_GO@result[c(8,7,4,1),]

which(down_GO@result$Description == "type I interferon production")
down_GO@result

selected_GO_down<-down_GO@result[c(617,222,57,18),]

selected_GO_down$Description<-factor(selected_GO_down$Description,levels=selected_GO_down$Description)
selected_GO_up$Description<-factor(selected_GO_up$Description,levels=selected_GO_up$Description)
#both_down_GO@result$Description == "regulation of T cell activation"
#library(ggpubr)
#library(ggrepel)
#library(scales)
#theme_set(theme_cowplot())
theme_set(theme_cowplot())
p1<-ggplot(selected_GO_down,aes(x=Description,y=pvalue,fill=Description))+
  geom_bar(stat = "identity",width=0.8)+coord_flip()+guides(fill=F)+
  labs(x="",y=expression(-log[10]~italic(P)),title = "Clone3 Down")+theme(axis.text.x =element_text(size=16),aspect.ratio = 1)+
  ylim(0,16)+scale_fill_manual(values = c(rep("#A5CDE2",6)))+
  theme(plot.title = element_text(hjust = 0.5))
p2<-ggplot(selected_GO_up,aes(x=Description,y=pvalue,fill=Description))+
  geom_bar(stat = "identity",width=0.8)+coord_flip()+guides(fill=F)+
  labs(x="",y=expression(-log[10]~italic(P)),title = "Clone3 Up")+theme(axis.text.x =element_text(size=16),aspect.ratio = 1)+
  ylim(0,40)+scale_fill_manual(values = c(rep("#FDBE6D",6)))+
  theme(plot.title = element_text(hjust = 0.5))
p2
p1+p2+plot_layout(ncol=1)

# statistic
Expression_df<-as.data.frame(t(as.data.frame(CRC4_st_tumor@assays$Spatial@data[c("STAT1","DDX58","IRF7","IRF1","IRF3","IRF9","ISG15","PSMB8","PSMB9","TAP1","TAP2","CXCL9","CXCL10","CXCL11",
                                                                            "MX1","OASL","MX2"),])))
Expression_df$Barcode <- rownames(Expression_df)
df2<-CRC4_st_tumor@meta.data[,c("Region","Barcode")]
Expression_df<-merge(Expression_df,df2,by="Barcode")

theme_set(theme_cowplot())
p1<-ggplot(Expression_df,aes(x=Region,y=PSMB9,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "PSMB9")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")

p2<-ggplot(Expression_df,aes(x=Region,y=PSMB8,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "PSMB8")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p3<-ggplot(Expression_df,aes(x=Region,y=TAP1,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "TAP1")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p4<-ggplot(Expression_df,aes(x=Region,y=TAP2,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "TAP2")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p5<-ggplot(Expression_df,aes(x=Region,y=STAT1,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "STAT1")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p6<-ggplot(Expression_df,aes(x=Region,y=IRF1,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "IRF1")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p9<-ggplot(Expression_df,aes(x=Region,y=IRF9,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "IRF9")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p10<-ggplot(Expression_df,aes(x=Region,y=MX1,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "MX1")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p11<-ggplot(Expression_df,aes(x=Region,y=CXCL9,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "CXCL9")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p12<-ggplot(Expression_df,aes(x=Region,y=CXCL10,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "CXCL10")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p13<-ggplot(Expression_df,aes(x=Region,y=CXCL11,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "CXCL11")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p14<-ggplot(Expression_df,aes(x=Region,y=ISG15,fill=Region,color=Region))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "ISG15")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p1+p2+p3+p4+p5+p6+p9+p10+p11+p12+p13+p14


ggplot(Expression_df,aes(x=Region,y=ERVH,fill=Region,color=Region))+
  geom_violin(width=1)+
  geom_boxplot(size=1,width=0.5,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "ERVH")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")+
  guides(fill=F,color=F)


p1<-SpatialPlot(CRC4_st,features = "HLA-A",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HLA-A",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))
p2<-SpatialPlot(CRC4_st,features = "HLA-B",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HLA-B",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))
p3<-SpatialPlot(CRC4_st,features = "TAP1",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "TAP1",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))
p4<-SpatialPlot(CRC4_st,features = "TAP2",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "TAP2",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))
p5<-SpatialPlot(CRC4_st,features = "B2M",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "B2M",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))
p1+p2+p5

#STAT1 binding genes
DEG_Clone3<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3",]

DefaultAssay(CRC4_TE_Tumor)<-"Spatial"
DEG<-FindAllMarkers(CRC4_TE_Tumor,assay = "Spatial",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,only.pos = F)
STAT1_target_genes<-read.csv("/data/huangzr/Spatial/Kras/version4/STAT1_target/STAT1_target_list.csv",header = T)
DEG_Clone3<-DEG[DEG$cluster == "Clone3",]
STAT1_DEG<-DEG_Clone3[DEG_Clone3$gene %in% STAT1_target_genes$X.6,]
STAT1_DEG$p_val_adj <- -log10(STAT1_DEG$p_val_adj)

STAT1_DEG$change<-"Remain"
STAT1_DEG$change[STAT1_DEG$p_val_adj > 2 & STAT1_DEG$avg_logFC > 0.25] = "Up"
STAT1_DEG$change[STAT1_DEG$p_val_adj > 2 & STAT1_DEG$avg_logFC < (-0.25)] = "Down"
for_label <- STAT1_DEG[STAT1_DEG$gene %in% c("TAP1","CXCL9","HLA-E","IFI27","PSMB9","CXCL10","IRF1","PSMB8","TAP2","IFI16"),]
ggplot(STAT1_DEG,aes(x=avg_logFC,y=p_val_adj,color=change))+geom_point(size=1.5)+
  scale_color_manual(values = c("#2C77A9","lightgrey","#D1332D"))+
  labs(x="log2(Clone3 vs Other clones)",y="-log2(P.adj)",title = "STAT1 binding genes")+
  theme(plot.title = element_text(size=12,hjust=0.5))+geom_vline(xintercept = 0.25,lty=2)+geom_vline(xintercept = -0.25,lty=2)+
  geom_hline(yintercept = 2,lty=2)+coord_cartesian(xlim = c(-1,1))+
  geom_text_repel(data = for_label[for_label$avg_logFC<0,],aes(label = gene),alpha = 1,color="#1F78B4",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(-0.5,-1.1),ylim = c(2,100),size=3)


#diff TE analysis
library(Seurat)
p_unload(Seurat)
p_load(Seurat)
CRC4_TE <- readRDS("/data/huangzr/Spatial/Kras/version5/Fig/NC_CRC2_TE.RDS")
as.matrix(CRC4_TE@assays$RNA@counts)
CRC4_st@assays$TE<-CreateAssayObject(data =  as.matrix(CRC4_TE@assays$RNA@data[rownames(CRC4_TE@assays$RNA@data) %in% repeats_info$V1,]))
saveRDS(CRC4_st,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/CRC4_st_with_TE.RDS")
table(colnames(CRC4_st@assays$TE) == rownames(CRC4_st@meta.data))
CRC4_TE<-CRC4_st
DefaultAssay(CRC4_TE)<-"TE"
CRC4_TE@meta.data<-CRC4_TE@meta.data[colnames(CRC4_TE),]

table(colnames(CRC4_TE) == rownames(CRC4_TE@meta.data))

CRC2_st2@assays$TE

table(colnames(CRC4_TE@assays$TE@data) == rownames(CRC4_TE@meta.data))
CRC4_TE@meta.data$HERVH<-colSums(CRC4_TE@assays$TE@data[c("LTR7Y","LTR7","HERVH-int"),])

CRC4_TE_Tumor@assays$TE@data[1:100,1:100]
CRC4_TE_Tumor<-CRC4_TE[,CRC4_TE$Region %in% c("Clone1","Clone2","Clone3")]
Idents(CRC4_TE_Tumor)<-CRC4_TE_Tumor$Region
diff_TE<-FindAllMarkers(CRC4_TE_Tumor,assay = "TE",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,only.pos = F)

View(diff_TE)
repeats_info<-read.csv("/data/huangzr/Spatial/Kras/spaceranger/repeats_family_info.csv",header = F)
Idents(CRC4_TE)
pre_TE_diff<-readRDS("/data/huangzr/Spatial/Kras/version5/Fig/TE_matrix_withGCP_TE.RDS")
pre_TE_diff
CRC_specific_TE<-read.csv("/data/huangzr/Spatial/Kras/version5/Fig/GCP_COAD_tumor_vs_normal_Diff_TE.csv")
CRC_specific_TE$change_tumor<-"Remain"
CRC_specific_TE$change_tumor[CRC_specific_TE$logFC>0] <- "Up"
CRC_specific_TE$change_tumor[CRC_specific_TE$logFC<0] <- "Down"
CRC_specific_TE<-CRC_specific_TE$X[CRC_specific_TE$change_tumor == "Up"]
saveRDS(diff_TE,"/data/huangzr/Spatial/Kras/version5/Fig_new/clone_TE_diff_new.RDS")

diff_TE


diff_CRC_TE<-diff_TE[diff_TE$gene %in% CRC_specific_TE,]

Clone3_TE<-diff_TE[diff_TE$cluster == "Clone3",]
Clone3_TE$p_val_adj <- -log10(Clone3_TE$p_val_adj)
Clone3_TE$change<-"Remain"
Clone3_TE$change[Clone3_TE$p_val_adj > 2 & Clone3_TE$avg_logFC > 0.25] = "Up"
Clone3_TE$change[Clone3_TE$p_val_adj > 2 & Clone3_TE$avg_logFC < (-0.25)] = "Down"
Clone3_TE$group<-"Remain"
Clone3_TE$group[Clone3_TE$gene %in% c("LTR7","LTR7Y","HERVH-int")]<-"HERVH"
#Clone3_TE$p_val_adj <- (-log10(Clone3_TE$p_val_adj))

theme_set(theme_cowplot())
library(ggrepel)
for_label <-  Clone3_TE[Clone3_TE$change != "Remain",]
for_label_2 <-  Clone3_TE[Clone3_TE$group == "HERVH",]
theme_set(theme_cowplot())

ggplot(Clone3_TE,aes(x=avg_logFC,y=p_val_adj,fill=change,color=group))+geom_point(shape=21,stroke=0.5,size=1.5)+
  scale_color_manual(values = c("black",alpha("white",0),rep(alpha("white",0),9)))+ 
  scale_fill_manual(values = c("#1E72A9","lightgrey","#D1332D"))+
  labs(x="log2(Clone3 vs Other Subclones)",y="-log2(P value)",title = "Differential Expressed Repeats(Clone3/Others)")+
  theme(plot.title = element_text(size=12,hjust=0.5))+geom_vline(xintercept = 0.25,lty=2)+geom_vline(xintercept = -0.25,lty=2)+
  geom_hline(yintercept = 2,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+
  geom_text_repel(data = for_label_2[for_label_2$avg_logFC<0,],aes(label = gene),alpha = 1,color="#1F78B4",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(-0.3,-1.2),size=4)+
  geom_text_repel(data = for_label[for_label$avg_logFC>0,],aes(label = gene),alpha = 1,color="#D1332D",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(0.3,1.2),size=4)


DefaultAssay(CRC4_TE) <- "Spatial"
SpatialPlot(CRC4_TE,features = "HERVH",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HERV-H",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))
Idents(CRC4_TE)<-"Region"
VlnPlot(CRC4_TE[,CRC4_TE$Region %in% c("Clone1","Clone2","Clone3")],features = "HERVH",pt.size = 0)

HERVH_expression<-CRC4_TE@meta.data[,c("Region","HERVH")]
ggplot(HERVH_expression[HERVH_expression$Region %in% c("Clone1","Clone2","Clone3"),],aes(x=Region,y=HERVH,fill=Region,color=Region))+geom_violin()+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F","grey"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00","lightgrey"))+
  labs(x="",y="Expression",title = "HERVH")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
CRC4_st
DefaultAssay(CRC4_st)<-"Spatial"
ERVH_element_exp<-t(CRC4_st@assays$TE[c("LTR7","LTR7Y","LTR7B","HERVH-int"),])
table(rownames(ERVH_element_exp) == rownames(CRC4_st@meta.data))
ERVH_element_exp<-as.data.frame(ERVH_element_exp)
ERVH_element_exp<-ERVH_element_exp[rownames(CRC4_st@meta.data),]
ERVH_element_exp$Barcode<-rownames(ERVH_element_exp)
CRC4_st@meta.data<-merge(ERVH_element_exp[,c(3,5)],CRC4_st@meta.data,by="Barcode")
rownames(CRC4_st@meta.data) <- CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
p1<-SpatialPlot(CRC4_st,features = c("LTR7Y"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "LTR7Y",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p2<-SpatialPlot(CRC4_st,features = c("LTR7B"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "LTR7B",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p3<-SpatialPlot(CRC4_st,features = c("LTR7"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "LTR7",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p4<-SpatialPlot(CRC4_st,features = c("HERVH-int"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HERVH-int",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p5<-SpatialPlot(CRC4_st,features = c("HLA-A"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HLA-A",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p6<-SpatialPlot(CRC4_st,features = c("HLA-B"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HLA-B",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p7<-SpatialPlot(CRC4_st,features = c("B2M"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "B2M",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))



CRC4_st@meta.data$HERVH<-rowMeans(CRC4_st@meta.data[,c("LTR7","LTR7Y","LTR7B","HERVH-int")])


SpatialPlot(CRC4_st,features = c("HERVH"),image.alpha = 0,stroke = NA,ncol = 1,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors =c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HERV-H",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p_unload(Seurat)
p_load(Seurat)
Idents(CRC4_st)<-factor(Idents(CRC4_st),levels=c("Clone1","Clone2","Clone3","Remaining"))
Idents(CRC4_st)<-"Region"
VlnPlot(CRC4_st[,CRC4_st$Region %in% c("Clone1","Clone2","Clone3")],features = "HERVH",pt.size = 0,
        cols = c("#A7CC84","#9FC4D6","#F5BA6F"))


p5+p6+p7+p1+p2+p3+p4+plot_layout(ncol=7)
