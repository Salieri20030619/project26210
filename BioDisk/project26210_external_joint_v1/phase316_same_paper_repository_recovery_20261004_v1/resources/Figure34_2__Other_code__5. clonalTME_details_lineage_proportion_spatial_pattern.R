#cell composition for each lineage
#group by clonalTME phenotype and clone(CRC4 slide)

colnames(tumor_prop_avg)
T_cell_cols<-c(4,5,6,2,3,28:30)
CD8_T_cell_cols<-c(4,5,6)
CD4_T_cell_cols<-c(2,3,28:30)
myeloid_cell_cols<-c(10,13,17,18,20,21,31:34)
mac_cell_cols<-c(18,10,13,21)
CD8_proportion<-tumor_prop_avg[,c(CD8_T_cell_cols)]/rowSums(tumor_prop_avg[,c(CD8_T_cell_cols)])
CD4_proportion<-tumor_prop_avg[,c(CD4_T_cell_cols)]/rowSums(tumor_prop_avg[,c(CD4_T_cell_cols)])
myeloid_proportion<-tumor_prop_avg[,c(myeloid_cell_cols)]/rowSums(tumor_prop_avg[,c(myeloid_cell_cols)])
mac_proportion<-tumor_prop_avg[,c(mac_cell_cols)]/rowSums(tumor_prop_avg[,c(mac_cell_cols)])
T_proportion<-tumor_prop_avg[,c(T_cell_cols)]/ropheatmap如何改colorkey的标题wSums(tumor_prop_avg[,c(T_cell_cols)])
CD8_proportion$clone<-rownames(CD8_proportion)
CD4_proportion$clone<-rownames(CD4_proportion)
mac_proportion$clone<-rownames(mac_proportion)
myeloid_proportion$clone<-rownames(myeloid_proportion)

CD8_proportion<-reshape2::melt(CD8_proportion)
CD4_proportion<-reshape2::melt(CD4_proportion)
mac_proportion<-reshape2::melt(mac_proportion)
myeloid_proportion<-reshape2::melt(myeloid_proportion)
mac_proportion<-mac_proportion[,c(1,2,4)]
ggplot(mac_proportion[mac_proportion$clone %in% c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),],
       aes(x=clone,y=value,fill=variable))+geom_bar(stat = "identity")


T_proportion$clone<-rownames(T_proportion)
T_proportion$clone<-rownames(T_proportion)
T_proportion$clone<-rownames(T_proportion)
CD8_proportion<-reshape2::melt(CD8_proportion)
CD4_proportion<-reshape2::melt(CD4_proportion)
myeloid_proportion<-reshape2::melt(myeloid_proportion)
mac_proportion<-reshape2::melt(mac_proportion)
T_proportion<-reshape2::melt(T_proportion)
CD8_proportion$Var1 <- factor(CD8_proportion$Var1,levels = c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2","CRC1 Clone2",
                                                                     "CRC2 Clone1","CRC3a Clone1","CRC3b Clone1",
                                                                     "CRC4 Clone3","CRC5 Clone1","CRC5 Clone2","CRC5 Clone3",
                                                                     "CRC6 Clone1","CRC7 Clone1"))
CD4_proportion$Var1 <- factor(CD4_proportion$Var1,levels = c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2","CRC1 Clone2",
                                                             "CRC2 Clone1","CRC3a Clone1","CRC3b Clone1",
                                                             "CRC4 Clone3","CRC5 Clone1","CRC5 Clone2","CRC5 Clone3",
                                                             "CRC6 Clone1","CRC7 Clone1"))
T_proportion$clone <- factor(T_proportion$clone,levels = c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2","CRC1 Clone2",
                                                             "CRC2 Clone1","CRC3a Clone1","CRC3b Clone1",
                                                             "CRC4 Clone3","CRC5 Clone1","CRC5 Clone2","CRC5 Clone3",
                                                             "CRC6 Clone1","CRC7 Clone1"))


p1<-ggplot(CD8_proportion,aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p1
p2<-ggplot(CD4_proportion,aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p2
p1+p2
myeloid_proportion$Var1 <- factor(myeloid_proportion$Var1,levels = c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2","CRC1 Clone2",
                                                                     "CRC2 Clone1","CRC3a Clone1","CRC3b Clone1",
                                                                     "CRC4 Clone3","CRC5 Clone1","CRC5 Clone2","CRC5 Clone3",
                                                                     "CRC6 Clone1","CRC7 Clone1"))
mac_proportion$Var1 <- factor(mac_proportion$Var1,levels = c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2","CRC1 Clone2",
                                                                     "CRC2 Clone1","CRC3a Clone1","CRC3b Clone1",
                                                                     "CRC4 Clone3","CRC5 Clone1","CRC5 Clone2","CRC5 Clone3",
                                                                     "CRC6 Clone1","CRC7 Clone1"))
library(RColorBrewer)
p3<-ggplot(myeloid_proportion,aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p3
p4<-ggplot(mac_proportion,aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p4

p1+p2+p3

# T cell abundance phenotype
clone_T_total <- rowMeans(tumor_prop_avg[,c(4,5,6,2,3,28:30)])
data.frame(T_cell=clone_T_total,group=tumor_prop_avg$subTME_category)

#myeloid abundance
CRC4_mye_total <- rowMeans(tumor_prop_avg[,c(10,13,17,18,20,21,31:34)])[6:8]
CRC4_mye_total<-as.data.frame(CRC4_mye_total)
CRC4_mye_total$clone<-rownames(CRC4_mye_total)
CRC4_mye_total<-reshape2::melt(CRC4_mye_total)
p1<-ggplot(CRC4_mye_total,aes(x=clone,y=value,fill=clone))+geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = c("#96C2D6","#90CA68","#F6B272"))+labs(x="clone",y="Proportion",title = "Myeloid")+
  theme(plot.title = element_text(hjust = 0.5))+guides(fill=F)

CRC4_T_total <- rowMeans(tumor_prop_avg[,c(4,5,6,2,3,28:30)])[6:8]
CRC4_T_total<-as.data.frame(CRC4_T_total)
CRC4_T_total$clone<-rownames(CRC4_T_total)
CRC4_T_total<-reshape2::melt(CRC4_T_total)
p2<-ggplot(CRC4_T_total,aes(x=clone,y=value,fill=clone))+geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = c("#96C2D6","#90CA68","#F6B272"))+labs(x="clone",y="Proportion",title = "T cell")+
  theme(plot.title = element_text(hjust = 0.5))+guides(fill=F)
library(patchwork)
p1+p2+plot_layout(ncol=1)
#cell proportion of CRC4 slide
p1<-SpatialPlot(CRC4_st,features = c("CD8.Tex"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.2,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8 Tex",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p1
p2<-SpatialPlot(CRC4_st,features = c("CD8.Tem"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.2,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8 Tem",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p2
p3<-SpatialPlot(CRC4_st,features = c("CD8.Trm"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.25,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8 Trm",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p3

p4<-SpatialPlot(CRC4_st,features = c("NK"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.15,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "NK",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p4
p5<-SpatialPlot(CRC4_st,features = c("Treg"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.2,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD4 Treg",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p5
p6<-SpatialPlot(CRC4_st,features = c("FN1..Mac"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.25,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "FN1..Mac",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p6
p7<-SpatialPlot(CRC4_st,features = c("IGF1..Mac"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.25,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "IGF1..Mac",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p7
p8<-SpatialPlot(CRC4_st,features = c("PDPN..Mac"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 1,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors = c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "PDPN..Mac",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p8
library(patchwork)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol=4)

#Roe calculate for each clone

#total_sample_tumor_prop$subTME<-"T cell-excluded"
#total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"
CRC4_subset<-total_sample_tumor_prop[total_sample_tumor_prop$slice == "CRC4",]
library(dplyr)
CRC4_subset %>% group_by(sample)
Roe_total<-aggregate(CRC4_subset[,1:34],by = list(CRC4_subset[,37]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)
celltype_list[1:29]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  select(celltype_list[1:29]) %>% t() %>% as.data.frame()
# Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
#   select(names(sort(kmeans_Roe_v_inhouse[["cluster"]]))) %>% t() %>% as.data.frame()

#colnames(Roe_v_inhouse) <- c('Dessert','Inflamed','Suppression')
max(Roe_v_inhouse)
min(Roe_v_inhouse)
library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "CRC4 Clone1",],`Ro/e`)$features.plot)
#arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


#dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-0.8,0.8))+labs(x="Ro/e",y="Subset")

ggsave("/data/huangzr/Spatial/Kras/version5/Fig_new/Celltype_enrich_subTME.pdf",height = 8,width = 5)


p1<-ggplot(dotplot_Roe_v_inhouse,
           aes(y=features.plot, x=id, fill= `Ro/e`, size = Group))+
  geom_point(shape=21) + theme_bw()+ labs(x='', y='')+
  theme(axis.text = element_text(color = 'black', size=14),
        legend.key.size = unit(3, 'mm'),
        panel.border = element_rect(color = 'black', size= 1),
        axis.ticks.length = unit(2, 'mm'),
        axis.ticks = element_line(colour = 'black'))+
  scale_fill_gradientn(limits = c(0,2),colours = colorRampPalette(c("#4A91C1","white","#FF7F00"))(100), 
                       guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'))
p1

#Proportion plot for each clone
library(reshape2)
CRC4_clone_prop_avg<-tumor_prop_avg[c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),]
CRC4_clone_prop_avg_df<-reshape2::melt(CRC4_clone_prop_avg)



library(gridExtra)
# 获取Var2列的所有唯一因子
factors <- unique(CRC4_clone_prop_avg_df$Var2)

# 初始化一个空列表来存储所有的图形
plot_list <- list()
i=1
CRC4_clone_prop_avg_df$value = CRC4_clone_prop_avg_df$value * 100
CRC4_clone_prop_avg_df$Var1 = gsub("CRC4 ","",CRC4_clone_prop_avg_df$Var1)
theme_set(theme_cowplot())
# 对每一个因子，创建一个图形
for (i in seq_along(factors)) {
  factor <- factors[i]
  
  # 对每一个因子，筛选出对应的数据子集
  data_subset <- CRC4_clone_prop_avg_df[CRC4_clone_prop_avg_df$Var2 == factor, ]
  
  # 创建图形
  p <- ggplot(data_subset, aes(x = Var1, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Clone", y = "Proportion(%)", title = factor) +
    scale_fill_manual(values = c("#F6B272","#96C2D6","#AA2E2E")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5))+guides(fill=F)
  
  # 将图形添加到列表中
  plot_list[[i]] <- p
}
p
# 使用grid.arrange函数将所有的图形拼在一起
do.call(grid.arrange, c(plot_list, ncol = 5))

#lineage transition for CRC4
p1<-ggplot(CD8_proportion[CD8_proportion$Var1 %in% c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),],aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p1
p2<-ggplot(CD4_proportion[CD4_proportion$Var1 %in% c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),],aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p2
library(RColorBrewer)
myeloid_proportion$Var2 <- factor(myeloid_proportion$Var2,levels=c("IGF1..Mac","FN1..Mac","PDPN..Mac","Monocyte","cDC1","cDC2","cDC3","pDC","Mast","Neutrophil"))
CRC4_clone_myeloid <-myeloid_proportion[myeloid_proportion$Var1 %in% c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),]
CRC4_clone_myeloid$clone <- as.character(CRC4_clone_myeloid$Var1)
CRC4_clone_myeloid$clone<-gsub("CRC4 ","",CRC4_clone_myeloid$Var1)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())
p3<-ggplot(CRC4_clone_myeloid,aes(x=clone,y=value,fill=Var2))+
  geom_bar(stat="identity",width=0.8)+scale_fill_manual(values =brewer.pal(10, "Paired"))+labs(x="Clone",y="Proportion(% total Myeloid)",fill="Subset")
p3
p4<-ggplot(mac_proportion[mac_proportion$Var1 %in% c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),],aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_manual(values =brewer.pal(12, "Paired"))
p4
p3+p4
CRC4_clone_T <-T_proportion[T_proportion$clone %in% c("CRC4 Clone1","CRC4 Clone2","CRC4 Clone3"),]
CRC4_clone_T$clone <- as.character(CRC4_clone_T$clone)
CRC4_clone_T$clone<-gsub("CRC4 ","",CRC4_clone_T$clone)
p5<-ggplot(CRC4_clone_T,aes(x=clone,y=value,fill=variable))+
  geom_bar(stat="identity",width=0.8)+scale_fill_manual(values =c(brewer.pal(8, "Set2"),brewer.pal(12, "Paired")[12]))+labs(x="Clone",y="Proportion(% total T cell)",fill="Subset")
p5
p3+p5


#Clonal distribution CRC4
library(Seurat)
library(pacman)
p_unload(Seurat)
p_load(Seurat)
SpatialPlot(CRC4_st,group.by = "Region",cols = c("#050404","#706F70","#D71219",alpha("white",0)),stroke = NA,pt.size.factor = 1.2)
SpatialPlot(CRC4_st,group.by = "Region",cols = c("#A6CEE3","#B2DF8A","#FB9A99",alpha("white",0)),stroke = NA,pt.size.factor = 1.4)


tumor_prop_avg_UCASpatial[]
