# Spatial Figure Process
library(pacman)
p_unload(Seurat)
p_load(Seurat)

load("/data/huangzr/Spatial/Kras/version4/public_CRC_bulk/TCGA/workspace_230710.RData")
NC_CRC2_st@meta.data

CRC2_TE_exp<-as.data.frame(NC_CRC2_TE@assays$RNA@data)

ERV1<-repeats_info$repeats[repeats_info$family == "ERV1"]

ERV1_level<-colMeans(CRC2_TE_exp[rownames(CRC2_TE_exp) %in% ERV1,])

ERV1_level<-as.data.frame(ERV1_level)
ERV1_level$barcode<-rownames(ERV1_level)
colnames(ERV1_level)[1]<-"ERV1"
NC_CRC2_st@meta.data<-merge(NC_CRC2_st@meta.data,ERV1_level,by="barcode",all.x = T)
rownames(NC_CRC2_st@meta.data)<-NC_CRC2_st@meta.data$barcode
NC_CRC2_st@meta.data<-NC_CRC2_st@meta.data[colnames(NC_CRC2_st),]

SpatialPlot(NC_CRC2_st,features = c("LTR7Y"))

NC_CRC2_st
library(Seurat)
VlnPlot(NC_CRC2_st[,NC_CRC2_st$Region=="Tumor"],features = "ERV1",group.by = "clone_region_update",
        cols = c("#A6CEE3","#B2DF8A","#FDBF6F"),pt.size = 0)
table(rownames(NC_CRC2_st@meta.data) %in% names(ERV1_level))


ERV1_clone<-NC_CRC2_st@meta.data[NC_CRC2_st@meta.data$Region=="Tumor",c(50,92,93,94,95)]
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
ERV1_clone$clone_region_update <- gsub("_tumor","",ERV1_clone$clone_region_update)

NC_CRC2_st@assays$TE<-CreateAssayObject(data=as.matrix(CRC2_TE_exp)[,rownames(NC_CRC2_st@meta.data)])
NC_CRC2_st@assays$TE@data[1:5,1:5]
DefaultAssay(NC_CRC2_st)<-"Spatial"
NC_CRC2_st@meta.data<-cbind(NC_CRC2_st@meta.data,t(NC_CRC2_st@assays$TE[c("LTR7Y","LTR7","HERVH-int"),]))
NC_CRC2_st@meta.data<-NC_CRC2_st@meta.data[colnames(NC_CRC2_st@assays$TE),]

NC_CRC2_st@meta.data$ERVH<-rowMeans(NC_CRC2_st@meta.data[,c("LTR7Y","LTR7","HERVH-int")])
VlnPlot(NC_CRC2_st[,NC_CRC2_st@meta.data$Region=="Tumor"],features = c("LTR7Y","LTR7","HERVH-int"),group.by = "clone_region_update",
        cols = c("#A6CEE3","#B2DF8A","#FDBF6F"),pt.size = 0)
table(colnames(NC_CRC2_st@assays$TE@data) == rownames(NC_CRC2_st@meta.data))


STAT1_act<-rowMeans(as.data.frame(t(as.matrix(NC_CRC2_st@assays$Spatial@data[rownames(NC_CRC2_st@assays$Spatial@data) %in% STAT1_list,]))))
NC_CRC2_tumor@meta.data$STAT1_activity<-STAT1_act[rownames(NC_CRC2_tumor@meta.data)]
VlnPlot(NC_CRC2_tumor,features = c("STAT1_activity"),group.by = "clone_region_update",
        cols = c("#A6CEE3","#B2DF8A","#FDBF6F"),pt.size = 0)
SpatialPlot(NC_CRC2_st,features = c("STAT1_activity"),image.alpha = 0,
            alpha = c(0,1))

SpatialPlot(NC_CRC2_st,features = c("CD3G"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "ERVH",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
CRC2_STAT1_act<-NC_CRC2_tumor@meta.data[,c("clone_region_update","STAT1_activity")]
CRC2_STAT1_act$clone_region_update <- gsub("_tumor","",CRC2_STAT1_act$clone_region_update)
ggplot(CRC2_STAT1_act,aes(x=clone_region_update,y=STAT1_activity,fill=clone_region_update,color=clone_region_update))+geom_violin()+
  geom_boxplot(size=1,width=0.4,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F)+guides(fill=F)+
  labs(x="",y="Activity",title = "STAT1 Activity")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,0.24))
compare_means(STAT1_activity~clone_region_update,CRC2_STAT1_act,method = "t.test")
ERV1_clone<-ERV1_clone[,c(1,3,4,5)]
ERV1_clone<-reshape2::melt(ERV1_clone)



ggplot(ERV1_clone,aes(x=variable,y=value,fill=clone_region_update,color=clone_region_update))+
  geom_boxplot(size=1,width=0.6,outlier.alpha = 0,alpha=0.2,position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F)+
  labs(x="",y="Expression",title = "ERVH in Tumor Region")+theme(plot.title = element_text(hjust = 0.5))
library(ggpubr)
compare_means(value~clone_region_update,ERV1_clone,method = "t.test",group.by = "variable")

c("ANXA1","KLF4","FAP","FPR2"))

rownames(NC_CRC2_st@assays$Spatial)[grep("IRF",rownames(NC_CRC2_st@assays$Spatial))]
Idents(NC_CRC2_st)<-NC_CRC2_st$clone_region_update
VlnPlot(NC_CRC2_st,features = c("IRF1","IRF9","STAT1","ISG15"),group.by = "clone_region_update")
VlnPlot(NC_CRC2_tumor,features = c("IRF1","IRF9","STAT1","ISG15"),
        group.by = "clone_region_update",pt.size = 0,cols = c("#A6CEE3","#B2DF8A","#FDBF6F"))

VlnPlot(NC_CRC2_tumor,features = rownames(NC_CRC2_st@assays$AUCell)[grep("RIG",rownames(NC_CRC2_st@assays$AUCell))][1:4],
        group.by = "clone_region_update",pt.size = 0,cols = c("#A6CEE3","#B2DF8A","#FDBF6F"))
SpatialPlot(NC_CRC2_st,features = c("IRF1","IRF3","IRF7","IRF9"),image.alpha = 0,
            alpha = c(0,1))
Pathway_Activity<-as.data.frame(t(NC_CRC2_st@assays$AUCell[c("RESPONSE-TO-VIRUS","CELLULAR-RESPONSE-TO-VIRUS","INTERFERON-MEDIATED-SIGNALING-PATHWAY",
                           "RESPONSE-TO-TYPE-I-INTERFERON","RESPONSE-TO-TYPE-II-INTERFERON","MDA-5-SIGNALING-PATHWAY","REGULATION-OF-MDA-5-SIGNALING-PATHWAY"),]))

Pathway_Activity<-Pathway_Activity[rownames(NC_CRC2_st@meta.data),]
table(rownames(Pathway_Activity) == rownames(NC_CRC2_st@meta.data))
Pathway_Activity<-cbind(NC_CRC2_st@meta.data,Pathway_Activity)[,c(50,97:103)]
Pathway_Activity<-reshape2::melt(Pathway_Activity)

Pathway_Activity$clone_region_update<-gsub("_tumor","",Pathway_Activity$clone_region_update)
ggplot(Pathway_Activity[Pathway_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Pathway_Activity$variable == "INTERFERON-MEDIATED-SIGNALING-PATHWAY",],
       aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin()+
  geom_boxplot(size=1,width=0.3,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Pathway Activity",title = "Interferon mediated signal")+theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 0.070805420,lty=2,size=1,alpha=0.6)
compare_means(value~clone_region_update,Pathway_Activity[Pathway_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Pathway_Activity$variable == "MDA-5-SIGNALING-PATHWAY",],method = "t.test",group.by = "variable")

quantile(Pathway_Activity$value[Pathway_Activity$variable == "INTERFERON-MEDIATED-SIGNALING-PATHWAY"])
ggplot(Pathway_Activity[Pathway_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Pathway_Activity$variable == "RESPONSE-TO-VIRUS",],
       aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+
  geom_density_ridges_gradient(aes(x=clone_region_update,y=value))
library(ggridges)

Pathway_Activity$clone_region_update <- factor(Pathway_Activity$clone_region_update)
p1<-ggplot(Pathway_Activity[Pathway_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Pathway_Activity$variable == "RESPONSE-TO-VIRUS",],
       aes(x=value,y=clone_region_update,group=clone_region_update)) +  
  geom_density_ridges_gradient(aes(fill = clone_region_update,color = clone_region_update), scale = 2, size = 0.5,alpha=0.5)+
  scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  geom_vline(xintercept = 0.0865,lty=2,size=0.5)+labs(x="Pathway Activity",y="",title = "Response to virus")+
  theme(plot.title = element_text(hjust = 0.5))


p2<-ggplot(Pathway_Activity[Pathway_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Pathway_Activity$variable == "INTERFERON-MEDIATED-SIGNALING-PATHWAY",],
       aes(x=value,y=clone_region_update,group=clone_region_update)) +  
  geom_density_ridges_gradient(aes(fill = clone_region_update,color = clone_region_update), scale = 2, size = 0.5,alpha=0.5)+
  scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  geom_vline(xintercept = 0.073,lty=2,size=0.5)+labs(x="Pathway Activity",y="",title = "Interferon mediated signaling")+
  theme(plot.title = element_text(hjust = 0.5))
theme_set(theme_cowplot())
p3<-ggplot(Pathway_Activity[Pathway_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Pathway_Activity$variable == "MDA-5-SIGNALING-PATHWAY",],
           aes(x=value,y=clone_region_update,group=clone_region_update)) +  
  geom_density_ridges_gradient(aes(fill = clone_region_update,color = clone_region_update), scale = 2, size = 0.5,alpha=0.5)+
  scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  geom_vline(xintercept = 0.073,lty=2,size=0.5)+labs(x="Pathway Activity",y="",title = "MDA-5 signaling")+
  theme(plot.title = element_text(hjust = 0.5))
p3
p1/p3


Exp_Activity<-as.data.frame(t(as.matrix(NC_CRC2_st@assays$Spatial[c("IRF1","IRF9","STAT1","ISG15","IFI16","IFI30"),])))

Exp_Activity<-Exp_Activity[rownames(NC_CRC2_st@meta.data),]
table(rownames(Exp_Activity) == rownames(NC_CRC2_st@meta.data))
Exp_Activity<-cbind(NC_CRC2_st@meta.data,Exp_Activity)[,c(50,97:102)]
Exp_Activity<-reshape2::melt(Exp_Activity)

Exp_Activity$clone_region_update<-gsub("_tumor","",Exp_Activity$clone_region_update)
p1<-ggplot(Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Exp_Activity$variable == "IRF1",],
       aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin(width=1.2)+
  geom_boxplot(size=1,width=0.12,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Expression",title = "IRF1")+theme(plot.title = element_text(hjust = 0.5))
p2<-ggplot(Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Exp_Activity$variable == "IRF9",],
           aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin(width=1.2)+
  geom_boxplot(size=1,width=0.12,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Expression",title = "IRF9")+theme(plot.title = element_text(hjust = 0.5))
p3<-ggplot(Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Exp_Activity$variable == "STAT1",],
           aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin(width=1.2)+
  geom_boxplot(size=1,width=0.12,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Expression",title = "STAT1")+theme(plot.title = element_text(hjust = 0.5))
p4<-ggplot(Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Exp_Activity$variable == "ISG15",],
           aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin(width=1.2)+
  geom_boxplot(size=1,width=0.12,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Expression",title = "ISG15")+theme(plot.title = element_text(hjust = 0.5))
p5<-ggplot(Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Exp_Activity$variable == "IFI16",],
           aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin(width=1.2)+
  geom_boxplot(size=1,width=0.12,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Expression",title = "IFI16")+theme(plot.title = element_text(hjust = 0.5))
p6<-ggplot(Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & Exp_Activity$variable == "IFI30",],
           aes(x=clone_region_update,y=value,fill=clone_region_update,color=clone_region_update))+geom_violin(width=1.2)+
  geom_boxplot(size=1,width=0.12,outlier.alpha = 0,position = position_dodge(width=0.9),fill="white")+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F,fill=F)+
  labs(x="",y="Expression",title = "IFI30")+theme(plot.title = element_text(hjust = 0.5))

library(patchwork)
p1+theme(axis.text.x=element_blank())+p3+theme(axis.text.x=element_blank())+p5+theme(axis.text.x=element_blank())+
  p2+theme(axis.text.x=element_blank())+p4+theme(axis.text.x=element_blank())+p6+theme(axis.text.x=element_blank())+plot_layout(ncol = 3)


compare_means(value~clone_region_update,
              Exp_Activity[Exp_Activity$clone_region_update %in% c("Clone1","Clone2","Clone3") & 
                             Exp_Activity$variable == "IFI30",],method = "t.test",group.by = "variable")


Clone3_STAT_lfc$change <- factor(Clone3_STAT_lfc$change,levels = c("up_regulated","down_regulated","remain"))
library(ggrepel)
label<-label[!(label$gene %in% c("DTX3E","DTX3L","ICAM1","CTSS","BST2")),]


ggplot(Clone3_STAT_lfc,aes(x=avg_logFC,y=p_val_adj,color=change))+geom_point()+
  scale_color_manual(values = c("#FE7F00","#1F78B4","lightgrey"))+coord_cartesian(xlim=c(-1,1))+
  labs(x="log2(Clone3 vs Other Clone)",y="-log10(adj.P)",title = "STAT1 Binding Genes")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "top",legend.text = element_text(hjust = 0.5))+geom_vline(xintercept = -0.25,lty=2)+geom_vline(xintercept = 0.25,lty=2)+
  geom_hline(yintercept = 2,lty=2)+geom_text_repel(data=label,aes(color=change, label=gene), 
                                                   point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), 
                                                   nudge_y = 0.1,xlim = c(-1.1,-0.6))+guides(color=F)

p1<-SpatialPlot(NC_CRC2_st,features = c("HLA-A"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HLA-A",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p2<-SpatialPlot(NC_CRC2_st,features = c("HLA-B"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HLA-B",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p3<-SpatialPlot(NC_CRC2_st,features = c("B2M"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "B2M",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p4<-SpatialPlot(NC_CRC2_st,features = c("PSMB9"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "PSMB9",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p5<-SpatialPlot(NC_CRC2_st,features = c("TAP1"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "TAP1",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p1+p2+p3+p4+p5+plot_layout(ncol = 5)
COAD_TPM$HERVH_expression
tail(colnames(COAD_TPM))
COAD_TPM$stage2 <- gsub("Stage ","",COAD_TPM$stage2)
COAD_TPM$stage2<-factor(COAD_TPM$stage2,levels = c("Normal","I","II","III","IV"))
COAD_TPM$metastasis <- "remain"
COAD_TPM$metastasis[COAD_TPM$stage2 == "Normal"]="Normal"
COAD_TPM$metastasis[COAD_TPM$stage2 %in% c("I","II","III")]="Non-metastasis"
COAD_TPM$metastasis[COAD_TPM$stage2 %in% c("IV")]="Metastasis"
COAD_TPM$metastasis<-factor(COAD_TPM$metastasis,levels = c("Normal","Non-metastasis","Metastasis"))
ggplot(COAD_TPM[!(is.na(COAD_TPM$stage2)),],aes(x=metastasis,y=HERVH_expression,color=metastasis,fill=metastasis))+
  geom_boxplot(size=1,outlier.alpha = 0)+scale_color_manual(values = c("#D3D3D3",rep("#1F78B4",1),rep("#FE7F00",1)))+
  scale_fill_manual(values = c("lightgrey",rep("#EDF5F9",1),rep("#FEF2E2",1)))+coord_cartesian(ylim = c(0,80))+
  labs(x=" ",y="Expression(CPM)",title="ERVH(TCGA-COAD)")+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust=1))+guides(fill=F,color=F)

COAD_TPM$HERVH_level<-COAD_TPM$HERVH_expression>5
table(COAD_TPM$HERVH_level[!(is.na(COAD_TPM$stage2))],COAD_TPM$stage2[!(is.na(COAD_TPM$stage2))])
df<-data.frame(stage=c("Normal","I","II","III","IV"),value=c(5.555556,68.75,64.7619,66.66667,57.14286))
ggplot(df,aes(x=stage,y=value))+geom_bar(stat = "identity")
compare_means(HERVH_expression~metastasis,COAD_TPM[!(is.na(COAD_TPM$stage2)),],method = "wilcox.test")


table(COAD_TPM_tumor_2$HERVH_group)
COAD_TPM_tumor_2<-COAD_TPM_tumor_2[order(COAD_TPM_tumor_2$HERVH_expression),]
COAD_TPM_tumor_2$HERVH_bin <- c(rep("Low",71),rep("Mid",140),rep("High",71))
COAD_TPM_tumor_2$HERVH_bin <- factor(COAD_TPM_tumor_2$HERVH_bin,levels = c("Low","Mid","High"))
theme_set(theme_cowplot())
ggplot(COAD_TPM_tumor_2,aes(x=HERVH_bin,y=CD))+geom_boxplot()
plot_list<-COAD_TPM_tumor_2[,c("HERVH_bin","CD3G","CD8A","CD4","FOXP3","LYZ","ITGAM","ITGAX","C1QA","C1QB","GZMA","GZMB",
                    "NKG7","PRF1","IL2RA","TNFRSF9","IFNG","STAT1","IRF1","IFI16","B2M","HLA_A","HLA_B","B2M","TAP1",
                    "TAP2","PSMB8","PSMB9")]
plot_list<-reshape2::melt(plot_list)
colnames(plot_list) <- c("ERVH_group","Gene","Exp")
ggplot(plot_list[plot_list$Gene=="CD8A",],aes(x=ERVH_group,y=Exp,fill= ERVH_group,color=ERVH_group))+
  stat_summary(fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7,position = position_dodge( .9))+
  stat_summary(fun.data='mean_sdl',fun.args = list(mult=1),geom="errorbar",width=0.2)

COAD_TPM_tumor_2$HERVH_bin <- factor(COAD_TPM_tumor_2$HERVH_bin,levels = c("High","Mid","Low"))
p1<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_bin,y=CD3G,fill=HERVH_bin,color=HERVH_bin))+
  geom_boxplot(size=1.5,width=0.8,outlier.alpha = 0)+scale_fill_manual(values = c("#EDF5F9","#EFF8E7","#FEF2E2"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F)+guides(fill=F)+
  labs(x="",y="Expression(TPM)",title = "CD3G")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,10))
p2<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_bin,y=STAT1,fill=HERVH_bin,color=HERVH_bin))+
  geom_boxplot(size=1.5,width=0.8,outlier.alpha = 0)+scale_fill_manual(values = c("#EDF5F9","#EFF8E7","#FEF2E2"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F)+guides(fill=F)+
  labs(x="",y="Expression(TPM)",title = "STAT1")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,300))
p3<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_bin,y=HLA_A,fill=HERVH_bin,color=HERVH_bin))+
  geom_boxplot(size=1.5,width=0.8,outlier.alpha = 0)+scale_fill_manual(values = c("#EDF5F9","#EFF8E7","#FEF2E2"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F)+guides(fill=F)+
  labs(x="",y="Expression(TPM)",title = "HLA-A")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,4000))


load("/data/huangzr/Spatial/Kras/version4/Figure/workspace_230722.RData")
NC_CRC2_tumor@meta.data

CRC2_ERVH<-NC_CRC2_st@meta.data[,c(1,92:96)]
NC_CRC2_tumor@meta.data<-merge(NC_CRC2_tumor@meta.data,CRC2_ERVH,by="barcode",all.x=T)
rownames(NC_CRC2_tumor@meta.data)<-NC_CRC2_tumor@meta.data$barcode
VlnPlot(NC_CRC2_tumor,features = "ERVH",group.by = "clone_region_update")
clone_ERVH<-NC_CRC2_tumor@meta.data[,c(50,101)]


antigen_gene<-as.data.frame(t(as.matrix(NC_CRC2_tumor@assays$Spatial[c("TAP1","TAP2","PSMB8","PSMB9"),])))
antigen_gene$barcode<-rownames(antigen_gene)
clone_ERVH$barcode<-rownames(clone_ERVH)
clone_ERVH<-merge(clone_ERVH,antigen_gene,by="barcode")

theme_set(theme_cowplot())
clone_ERVH$clone_region_update <- gsub("_tumor","",clone_ERVH$clone_region_update)
clone_ERVH<-clone_ERVH[,c(2,4:7)]

p1<-ggplot(clone_ERVH,aes(x=clone_region_update,y=TAP1,fill=clone_region_update,color=clone_region_update))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "TAP1")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,3))+
  theme(axis.text.x = element_blank())+guides(color=F)
p2<-ggplot(clone_ERVH,aes(x=clone_region_update,y=TAP2,fill=clone_region_update,color=clone_region_update))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "TAP2")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,2))+
  theme(axis.text.x = element_blank())+guides(color=F)
p3<-ggplot(clone_ERVH,aes(x=clone_region_update,y=PSMB8,fill=clone_region_update,color=clone_region_update))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "PSMB8")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,3))+
  theme(axis.text.x = element_blank())+guides(color=F)
p4<-ggplot(clone_ERVH,aes(x=clone_region_update,y=PSMB9,fill=clone_region_update,color=clone_region_update))+
  geom_boxplot(size=1,width=0.7,outlier.alpha = 0,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+
  labs(x="",y="Expression",title = "PSMB9")+theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(ylim = c(0,3))+
  theme(axis.text.x = element_blank())+guides(color=F)

p1+p2+p3+p4
library(ggpubr)
compare_means(PSMB8~clone_region_update,clone_ERVH,method = "t.test")

DefaultAssay(NC_CRC2_st) <- "TE"
NC_CRC2_st_tumor<-NC_CRC2_st[,NC_CRC2_st$Tumor == "Tumor"]
Idents(NC_CRC2_st_tumor)<-NC_CRC2_st_tumor$clone_region_update
DefaultAssay(NC_CRC2_st_tumor)<-"TE2"
TE_diff<-FindAllMarkers(NC_CRC2_st_tumor,logfc.threshold = 0)
NC_CRC2_st_tumor@assays$TE@data[1:5,1:5]

NC_CRC2_st_tumor@assays$TE2<- CreateAssayObject(data = NC_CRC2_st_tumor@assays$TE@data[rownames(NC_CRC2_st_tumor@assays$TE@data) %in% repeats_info$repeats[!(repeats_info$class %in% c("Simple_repeat","Satellite"))],])
View(TE_diff)
TE_diff$p_val_adj <- -log10(TE_diff$p_val_adj)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
TE_diff$change<-"Remain"
TE_diff$change[TE_diff$p_val_adj > 2 & TE_diff$avg_logFC < (-0.25)]<-"Down"
TE_diff$change[TE_diff$p_val_adj > 2 & TE_diff$avg_logFC > (0.25)]<-"Up"
TE_diff$Group<-"Remain"
TE_diff$Group[TE_diff$gene %in% c("HERVH-int","LTR7",'LTR7Y')]<-"ERVH"
ggplot(TE_diff[TE_diff$cluster == "Clone3_tumor",],aes(x=avg_logFC,y=p_val_adj,fill=change,color=Group))+
  geom_point(shape=21,size=3,stroke=1)+
  coord_cartesian(xlim = c(-1.5,1.5))+labs(x="log2(Clone3 vs Other Clone)",y="-log10(adj.P)",title = "Differential TE in Clone3")+
  scale_fill_manual(values = c("#1F78B4","#D3D3D3","#FE7F01"))+geom_hline(yintercept = 2,lty=2)+
  geom_vline(xintercept = 0.25,lty=2)+geom_vline(xintercept = -0.25,lty=2)+
  scale_color_manual(values = c("black",alpha("white",0)))+
  geom_text_repel(data=label,aes(label=gene),point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"),
                  nudge_y = 0.1,xlim = c(-1.5,-0.3))+theme(plot.title = element_text(hjust=0.5))+guides(fill=F,color=F)




library(ggrepel)
label<-TE_diff[(TE_diff$gene %in% c("HERVH-int","LTR7",'LTR7Y')) & TE_diff$cluster == "Clone3_tumor",]


ggplot(Clone3_STAT_lfc,aes(x=avg_logFC,y=p_val_adj,color=change))+geom_point()+
  scale_color_manual(values = c("#FE7F00","#1F78B4","lightgrey"))+coord_cartesian(xlim=c(-1,1))+
  labs(x="log2(Clone3 vs Other Clone)",y="-log10(adj.P)",title = "STAT1 Binding Genes")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "top",legend.text = element_text(hjust = 0.5))+geom_vline(xintercept = -0.25,lty=2)+geom_vline(xintercept = 0.25,lty=2)+
  geom_hline(yintercept = 2,lty=2)+geom_text_repel(data=label,aes(color=change, label=gene), 
                                                   point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), 
                                                   nudge_y = 0.1,xlim = c(-1.1,-0.6))+guides(color=F)


SpatialPlot(NC_CRC2_st,features = c("STAT1"),image.alpha = 0,stroke = NA,alpha = c(0,1),ncol = 1)+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "STAT1",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

table(names(NC_CRC2_st@assays$Spatial@data["STAT1",]) == rownames(NC_CRC2_st@meta.data))
NC_CRC2_st@meta.data$STAT1<-NC_CRC2_st@assays$Spatial@data["STAT1",]
STAT1_clone<-NC_CRC2_st@meta.data[,50:51]
library(cowplot)
theme_set(theme_cowplot())
ggplot(STAT1_clone[STAT1_clone$clone_region_update %in% c("Clone1_tumor","Clone2_tumor","Clone3_tumor"),],aes(x=clone_region_update,y=STAT1,fill=clone_region_update,color=clone_region_update))+
  geom_violin(width=1.1)+geom_boxplot(size=1,width=0.25,outlier.alpha = 0,alpha=1,fill="white",position = position_dodge(width=0.8))+scale_fill_manual(values = c("#A6CEE3","#B2DF8A","#FDBF6F"))+
  scale_color_manual(values = c("#1F78B4","#33A02C","#FE7F00"))+guides(color=F)+
  labs(x="",y="Expression",title = "STAT1")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.x = element_blank())
library(ggpubr)
compare_means(value~clone_region_update,ERV1_clone,method = "t.test",group.by = "variable")

#CRC2 CNV
CRC2_CNV<-read.csv("/data/huangzr/Spatial/Kras/version4/Tumor_region_identified_by_CNV/CRC2_CNV.csv")
View(CRC2_CNV)

#TCGA
TCGA_ERV_data<-readRDS("/data/huangzr/Spatial/Kras/version4/public_CRC_bulk/TCGA/COAD_data_v2.RDS")
load("/data/huangzr/Spatial/Kras/version4/public_CRC_bulk/TCGA/workspace_230710.RData")
colnames(TCGA_ERV_data)[59000:ncol(TCGA_ERV_data)]

COAD_TPM_tumor_2$HERVH_group<-"Mid"
COAD_TPM_tumor_2<-COAD_TPM_tumor_2[order(COAD_TPM_tumor_2$HERVH_expression),]

COAD_TPM_tumor_2$HERVH_group[(282-69):282]="High"
COAD_TPM_tumor_2$HERVH_group[1:70]="Low"
table(COAD_TPM_tumor_2$HERVH_group)

COAD_TPM_tumor_2$HERVH_expression[HERVH_rank=113]
COAD_TPM_tumor_2$HERVH_group<-factor(COAD_TPM_tumor_2$HERVH_group,levels = c("Low","Mid","High"))
p1<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=STAT1,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "STAT1")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,300))
p1
p2<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=IRF1,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "IRF1")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,120))
p2
p3<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=IRF7,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "IRF7")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,90))
p3
p4<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=ISG15,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "ISG15")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,600))
p4
p5<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=HLA_B,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "HLA-B")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,6000))
p5
p6<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=TAP1,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "TAP1")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,500))
p6
p7<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=IFNG,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "IFNG")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,4))
p7
p8<-ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=PSMB9,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "PSMB9")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,300))
p8
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol =4)
compare_means(STAT1~HERVH_group,COAD_TPM_tumor_2)


ggplot(COAD_TPM_tumor_2,aes(x=HERVH_group,y=IFNG,color=HERVH_group,fill=HERVH_group))+geom_boxplot(width=0.7,lwd=1,outlier.alpha = 0)+
  scale_color_manual(values = c("#FE7F00","#33A02C","#1F78B4"))+scale_fill_manual(values = c("#FEF2E2","#EFF8E7","#EDF5F9"))+
  labs(x="ERVH",y="Expression(TPM)",title = "IFNG")+guides(fill=F,color=F)+theme(plot.title = element_text(hjust=0.5))+
  coord_cartesian(ylim = c(0,4))


compare_means(PSMB9~HERVH_group,COAD_TPM_tumor_2)

VlnPlot(NC_CRC2_tumor,features = c("CASP8","CASP3"),
        group.by = "clone_region_update",pt.size = 0,ncol=4,cols = c("#A6CEE3","#B2DF8A","#FDBF6F"))

#CNV profie in spatial
CRC2_CNV

gc()
gc()
