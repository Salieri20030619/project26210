library(pacman)
p_unload(ggplot2)
p_load(ggplot2)
C2L_Seurat<-readRDS("/Path/data/C2L_results_ICB.rds")
C2L_Seurat<-C2L_Seurat[[1]]
RCTD_Seurat<-readRDS("/Path/data/RCTD_results_ICB.rds")
RCTD_Seurat<-RCTD_Seurat[[1]]
SPOTlight_Seurat<-readRDS("/Path/data/SPOTlight_results_ICB.rds")
SPOTlight_Seurat<-SPOTlight_Seurat[[1]]
CARD_Seurat<-readRDS("/Path/data/CARD_results_ICB.rds")
CARD_Seurat<-CARD_Seurat[[1]]
stereoscope_Seurat<-readRDS("/Path/data/stereoscope_results_ICB.rds")
stereoscope_Seurat<-stereoscope_Seurat[[1]]
ICB_sliceA<-readRDS("/Path/data/ICB_section_A_region.RDS")
Spotiphy_Seurat<-readRDS("/Path/data/Spotiphy_results_ICB.rds")
Spotiphy_Seurat<-Spotiphy_Seurat[[1]]
plot_data<-readRDS("/Path/data/ICB_region.RDS")
total_panel_statistics<- readRDS(file = "/Path/data/total_panel_statistics.RDS")


colnames(C2L_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L_Seurat@meta.data))))))))))))
colnames(RCTD_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD_Seurat@meta.data))))))))))))
colnames(SPOTlight_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight_Seurat@meta.data))))))))))))
colnames(CARD_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD_Seurat@meta.data))))))))))))
colnames(stereoscope_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope_Seurat@meta.data))))))))))))
colnames(ICB_sliceA@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(ICB_sliceA@meta.data))))))))))))
colnames(Spotiphy_Seurat@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(Spotiphy_Seurat@meta.data))))))))))))

Tumor1_spot<-rownames(plot_data)[plot_data$tumor_region == "Tumor1"]
Tumor2_spot<-rownames(plot_data)[plot_data$tumor_region == "Tumor2"]
TLS_spot<-rownames(plot_data)[plot_data$tumor_region == "TLS"]
Stroma_spot<-rownames(plot_data)[plot_data$tumor_region == "Stroma"]

C2L_Seurat@meta.data$tumor_region<-"Remaining"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% TLS_spot]<-"TLS"
C2L_Seurat@meta.data$tumor_region[rownames(C2L_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

RCTD_Seurat@meta.data$tumor_region<-"Remaining"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% TLS_spot]<-"TLS"
RCTD_Seurat@meta.data$tumor_region[rownames(RCTD_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

SPOTlight_Seurat@meta.data$tumor_region<-"Remaining"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% TLS_spot]<-"TLS"
SPOTlight_Seurat@meta.data$tumor_region[rownames(SPOTlight_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

CARD_Seurat@meta.data$tumor_region<-"Remaining"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% TLS_spot]<-"TLS"
CARD_Seurat@meta.data$tumor_region[rownames(CARD_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

stereoscope_Seurat@meta.data$tumor_region<-"Remaining"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% TLS_spot]<-"TLS"
stereoscope_Seurat@meta.data$tumor_region[rownames(stereoscope_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

Spotiphy_Seurat@meta.data$tumor_region<-"Remaining"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% Tumor1_spot]<-"Tumor1"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% Tumor2_spot]<-"Tumor2"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% TLS_spot]<-"TLS"
Spotiphy_Seurat@meta.data$tumor_region[rownames(Spotiphy_Seurat@meta.data) %in% Stroma_spot]<-"Stroma"

Ablation_Seurat@meta.data$tumor_region

ICB_sliceA@meta.data$tumor_region<-"Remaining"
ICB_sliceA@meta.data$tumor_region[rownames(ICB_sliceA@meta.data) %in% Tumor1_spot]<-"Tumor1"
ICB_sliceA@meta.data$tumor_region[rownames(ICB_sliceA@meta.data) %in% Tumor2_spot]<-"Tumor2"

ICB_sliceA@meta.data$TLS_score<-rowMeans(ICB_sliceA@meta.data[,c("B","Tfh","cDC3","CD4+ Tn")])-rowMeans(ICB_sliceA@meta.data[,c("CD8+ Trm","FAP+ Myofib","FN1+ Mac","NK","Epithelial c3","Epithelial c4","Epithelial c5")])
#ICB_sliceA<-UpdateSeuratObject(ICB_sliceA)

ICB_sliceA@meta.data$TLS_defination<-"Non_TLS"
ICB_sliceA@meta.data$TLS_defination[ICB_sliceA@meta.data$TLS_score > 0.07]<-"TLS"
ICB_sliceA@meta.data$tumor_region[ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] > 12000 & ICB_sliceA@meta.data$TLS_defination == "TLS"]<- "TLS"

ICB_sliceA@meta.data$stroma_defination<-"Remain"
ICB_sliceA@meta.data$stroma_defination[ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] > 11600 &
                                         ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] < 12300 &
                                         ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,1] > 19700 &
                                         ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,1] < 20700]<-"Stroma_selected"
ICB_sliceA@meta.data$tumor_region[ICB_sliceA@meta.data$stroma_defination == "Stroma_selected"]<-"Stroma_selected2"

DefaultAssay(ICB_sliceA) <- "Annotation"
SpatialPlot(ICB_sliceA,group.by = "stroma_defination",pt.size.factor = 2)
#SpatialPlot(ICB_sliceA[,ICB_sliceA@images$sliceA1@boundaries$centroids@coords[,2] > 12000],group.by = "TLS_defination",pt.size.factor = 3)
SpatialPlot(ICB_sliceA,group.by = "tumor_region",pt.size.factor = 3,cols = c("white",""))

C2L_Seurat@meta.data$TREM2_mac<-rowSums(C2L_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
RCTD_Seurat@meta.data$TREM2_mac<-rowSums(RCTD_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
SPOTlight_Seurat@meta.data$TREM2_mac<-rowSums(SPOTlight_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
CARD_Seurat@meta.data$TREM2_mac<-rowSums(CARD_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
stereoscope_Seurat@meta.data$TREM2_mac<-rowSums(stereoscope_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
ICB_sliceA@meta.data$TREM2_mac<-rowSums(ICB_sliceA@meta.data[,c("FN1+ Mac","PDPN+ Mac")])
Spotiphy_Seurat@meta.data$TREM2_mac<-rowSums(Spotiphy_Seurat@meta.data[,c("FN1+ Mac","PDPN+ Mac")])

ST_tumor_combine_C2L<-aggregate(C2L_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                        "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                list(Region = C2L_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_C2L<-ST_tumor_combine_C2L[2:4,]
ST_tumor_combine_C2L<-reshape2::melt(ST_tumor_combine_C2L)

ST_tumor_combine_RCTD<-aggregate(RCTD_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                          "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                 list(Region = RCTD_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_RCTD<-ST_tumor_combine_RCTD[2:4,]
ST_tumor_combine_RCTD<-reshape2::melt(ST_tumor_combine_RCTD)

ST_tumor_combine_SPOTlight<-aggregate(SPOTlight_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                                    "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                      list(Region = SPOTlight_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_SPOTlight<-ST_tumor_combine_SPOTlight[2:4,]
ST_tumor_combine_SPOTlight<-reshape2::melt(ST_tumor_combine_SPOTlight)

ST_tumor_combine_UCAS<-aggregate(ICB_sliceA@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                         "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                 list(Region = ICB_sliceA@meta.data$tumor_region),mean)

ST_tumor_combine_UCAS<-ST_tumor_combine_UCAS[3:5,]
ST_tumor_combine_UCAS<-reshape2::melt(ST_tumor_combine_UCAS)

ST_tumor_combine_CARD<-aggregate(CARD_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                          "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                 list(Region = CARD_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_CARD<-ST_tumor_combine_CARD[2:4,]
ST_tumor_combine_CARD<-reshape2::melt(ST_tumor_combine_CARD)

ST_tumor_combine_stereoscope<-aggregate(stereoscope_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                                        "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                        list(Region = stereoscope_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_stereoscope<-ST_tumor_combine_stereoscope[2:4,]
ST_tumor_combine_stereoscope<-reshape2::melt(ST_tumor_combine_stereoscope)

ST_tumor_combine_Spotiphy<-aggregate(Spotiphy_Seurat@meta.data[,c("CD8+ Tem", "CD8+ Tex", "CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg",
                                                                  "TREM2_mac","C1QC+ Mac","B","cDC1","cDC2","cDC3")],
                                     list(Region = Spotiphy_Seurat@meta.data$tumor_region),mean)

ST_tumor_combine_Spotiphy<-ST_tumor_combine_Spotiphy[2:4,]
ST_tumor_combine_Spotiphy<-reshape2::melt(ST_tumor_combine_Spotiphy)

ST_tumor_combine_Ablation


ST_tumor_combine_UCAS$Region == ST_tumor_combine_C2L$Region
ST_tumor_combine_3<-cbind(ST_tumor_combine_UCAS,ST_tumor_combine_C2L[,3],ST_tumor_combine_RCTD[,3],
                          ST_tumor_combine_SPOTlight[,3],ST_tumor_combine_CARD[,3],ST_tumor_combine_stereoscope[,3],ST_tumor_combine_Spotiphy[,3])

levels(ST_tumor_combine_3$variable)
gsub("C1QC mac","C1QC+ mac",gsub("CD4+ Treg","Treg",gsub("CD4+ Tfh","Tfh",gsub("CD4","CD4+",gsub("CD8","CD8+",gsub(" cells","",total_panel_statistics$celltype))))))

levels(ST_tumor_combine_3$variable)[c(8,13)]<-c("TREM2+ mac","CCR7 DC")
levels(ST_tumor_combine_3$variable)[c(9)]<-"C1QC+ mac"
ST_tumor_combine_3<-ST_tumor_combine_3[ST_tumor_combine_3$variable != "CD4+ Tm",]
ST_tumor_combine_3$variable == total_panel_statistics$celltype
ST_tumor_combine_3$Region == total_panel_statistics$Region

#ST_tumor_combine_3$variable <- mIHC_tumor_combine$celltype
ST_tumor_combine_3$proportion_mIHC <- total_panel_statistics$Proportion
colnames(ST_tumor_combine_3)[3:9]<-c("proportion_UCAS","proportion_C2L","proportion_RCTD",
                                     "proportion_SPOTlight","proportion_CARD","proportion_stereoscope","proportion_Spotiphy")
colnames(ST_tumor_combine_3)[2]<-"celltype"
ST_tumor_combine_3$proportion_Ablation = ST_tumor_combine_Ablation$value[ST_tumor_combine_Ablation$Region != "Stroma"]
ST_tumor_combine_Ablation<-ST_tumor_combine_Ablation[ST_tumor_combine_Ablation$variable != "CD4+ Tm",]


#ST_tumor_combine_2$proportion_mIHC = mIHC_tumor_combine$proportion
#ST_tumor_combine_2$proportion_mIHC<-mIHC_tumor_combine$proportion
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)

ST_tumor_combine_3 <- readRDS("/Path/data/ST_tumor_combine_3.RDS")

ST_tumor_combine_3_tumor<-ST_tumor_combine_3[ST_tumor_combine_3$Region %in% c("Tumor1","Tumor2"),]
p1<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_UCAS*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (UCASpatial)",y="Cell Proportion%DAPI (mIHC)",title = "UCASpatial")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))
p1
p2<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_C2L*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (cell2location)",y="Cell Proportion%DAPI (mIHC)",title = "cell2location")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))

p3<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_RCTD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (RCTD)",y="Cell Proportion%DAPI (mIHC)",title = "RCTD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))

p4<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_SPOTlight*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (SPOTlight)",y="Cell Proportion%DAPI (mIHC)",title = "SPOTlight")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,15),ylim = c(0,15))

p5<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_CARD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (CARD)",y="Cell Proportion%DAPI (mIHC)",title = "CARD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,40),ylim = c(0,40))

p6<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_stereoscope*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (stereoscope)",y="Cell Proportion%DAPI (mIHC)",title = "stereoscope")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,5),ylim = c(0,5))

p7<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_Spotiphy*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Spotiphy)",y="Cell Proportion%DAPI (mIHC)",title = "Spotiphy")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,7),ylim = c(0,7))

p8<-ggplot(ST_tumor_combine_3_tumor,aes(x=proportion_Ablation*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Ablation)",y="Cell Proportion%DAPI (mIHC)",title = "Ablation")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,10),ylim = c(0,10))

library(patchwork)
p1+p8+p2+p3+p4+p5+p6+p7+plot_layout(ncol=4)
library(Metrics)
PCC_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                   c(cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_UCAS)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_C2L)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_RCTD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_SPOTlight)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_CARD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_stereoscope)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_Spotiphy)[4]$estimate,
                     cor.test(ST_tumor_combine_3_tumor$proportion_mIHC,ST_tumor_combine_3_tumor$proportion_Ablation)[4]$estimate))
colnames(PCC_df)<-c("method","PCC")
PCC_df$method<-factor(PCC_df$method,levels = c("UCASpatial","RCTD","Ablation","Spotiphy","CARD","cell2location","stereoscope","SPOTlight"))
levels(PCC_df$method)[3] <- "UCASpatial_Ablation"
p1<-ggplot(PCC_df,aes(x=method,y=PCC,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#E6D77A","#CBEBC5","#d9d9d9","#add26a","#b37cb3","#7eabca"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1

RMSE_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                    c(rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_UCAS*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_C2L*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_RCTD*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_SPOTlight*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_CARD*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_stereoscope*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_Spotiphy*100),
                      rmse(ST_tumor_combine_3_tumor$proportion_mIHC, ST_tumor_combine_3_tumor$proportion_Ablation*100)))
colnames(RMSE_df)<-c("method","RMSE")
RMSE_df$method<-factor(RMSE_df$method,levels = c("UCASpatial","RCTD","stereoscope","cell2location","Spotiphy","Ablation","SPOTlight","CARD"))
p2<-ggplot(RMSE_df,aes(x=method,y=RMSE,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#b37cb3","#add26a","#CBEBC5","#E6D77A","#7eabca","#d9d9d9"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
p1+p2


ST_tumor_combine_3_TLS<-ST_tumor_combine_3[ST_tumor_combine_3$Region %in% c("TLS"),]
p1<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_UCAS*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (UCASpatial)",y="Cell Proportion%DAPI (mIHC)",title = "UCASpatial")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))
p1
p2<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_C2L*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (cell2location)",y="Cell Proportion%DAPI (mIHC)",title = "cell2location")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,60),ylim = c(0,60))

p3<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_RCTD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (RCTD)",y="Cell Proportion%DAPI (mIHC)",title = "RCTD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p4<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_SPOTlight*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (SPOTlight)",y="Cell Proportion%DAPI (mIHC)",title = "SPOTlight")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,30),ylim = c(0,30))

p5<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_CARD*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (CARD)",y="Cell Proportion%DAPI (mIHC)",title = "CARD")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p6<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_stereoscope*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (stereoscope)",y="Cell Proportion%DAPI (mIHC)",title = "stereoscope")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p7<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_Spotiphy*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Spotiphy)",y="Cell Proportion%DAPI (mIHC)",title = "Spotiphy")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

p8<-ggplot(ST_tumor_combine_3_TLS,aes(x=proportion_Ablation*100,y=proportion_mIHC,color=celltype,shape=Region))+
  geom_point(size=3)+scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2")))+
  labs(x="Cell Proportion% (Ablation)",y="Cell Proportion%DAPI (mIHC)",title = "Ablation")+
  geom_abline(slope = 1,intercept = 0,lty="dashed")+coord_cartesian(xlim = c(0,25),ylim = c(0,25))

library(patchwork)
p1+p8+p2+p3+p4+p5+p6+p7+plot_layout(ncol=4)
library(Metrics)
PCC_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                   c(cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_UCAS)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_C2L)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_RCTD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_SPOTlight)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_CARD)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_stereoscope)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_Spotiphy)[4]$estimate,
                     cor.test(ST_tumor_combine_3_TLS$proportion_mIHC,ST_tumor_combine_3_TLS$proportion_Ablation)[4]$estimate))
colnames(PCC_df)<-c("method","PCC")
PCC_df$method<-factor(PCC_df$method,levels = c("UCASpatial","Spotiphy","cell2location","RCTD","stereoscope","Ablation","SPOTlight","CARD"))
levels(PCC_df$method)[6] <- "UCASpatial_Ablation"
p1<-ggplot(PCC_df,aes(x=method,y=PCC,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#CBEBC5","#add26a","#f4b062","#b37cb3","#E6D77A","#7eabca","#d9d9d9"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1

RMSE_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                    c(rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_UCAS*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_C2L*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_RCTD*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_SPOTlight*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_CARD*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_stereoscope*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_Spotiphy*100),
                      rmse(ST_tumor_combine_3_TLS$proportion_mIHC, ST_tumor_combine_3_TLS$proportion_Ablation*100)))
colnames(RMSE_df)<-c("method","RMSE")
RMSE_df$method<-factor(RMSE_df$method,levels = c("UCASpatial","Spotiphy","RCTD","stereoscope","Ablation","CARD","SPOTlight","cell2location"))
p2<-ggplot(RMSE_df,aes(x=method,y=RMSE,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#CBEBC5","#f4b062","#b37cb3","#E6D77A","#d9d9d9","#7eabca","#add26a"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
p1+p2


# T cell in tumor
ST_tumor_combine_Tcell<-ST_tumor_combine_3_tumor[ST_tumor_combine_3_tumor$celltype %in% c("CD8+ Tem","CD8+ Tex","CD8+ Trm","CD4+ Tm","Tfh","Th1/Th17","Treg"),]
PCC_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                   c(cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_UCAS)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_C2L)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_RCTD)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_SPOTlight)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_CARD)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_stereoscope)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_Spotiphy)[4]$estimate,
                     cor.test(ST_tumor_combine_Tcell$proportion_mIHC,ST_tumor_combine_Tcell$proportion_Ablation)[4]$estimate))
colnames(PCC_df)<-c("method","PCC")
PCC_df$method<-factor(PCC_df$method,levels = c("UCASpatial","RCTD","cell2location","stereoscope","Spotiphy","CARD","SPOTlight","Ablation"))
levels(PCC_df$method)[8] <- "UCASpatial_Ablation"
p1<-ggplot(PCC_df,aes(x=method,y=PCC,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#add26a","#b37cb3","#CBEBC5","#d9d9d9","#7eabca","#E6D77A"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1

RMSE_df<-data.frame(method=c("UCASpatial","cell2location","RCTD","SPOTlight","CARD","stereoscope","Spotiphy","Ablation"),
                    c(rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_UCAS*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_C2L*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_RCTD*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_SPOTlight*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_CARD*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_stereoscope*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_Spotiphy*100),
                      rmse(ST_tumor_combine_Tcell$proportion_mIHC, ST_tumor_combine_Tcell$proportion_Ablation*100)))
colnames(RMSE_df)<-c("method","RMSE")
RMSE_df$method<-factor(RMSE_df$method,levels = c("UCASpatial","RCTD","stereoscope","cell2location","Spotiphy","Ablation","SPOTlight","CARD"))
levels(RMSE_df$method)[6] <- "UCASpatial_Ablation"
p2<-ggplot(RMSE_df,aes(x=method,y=RMSE,fill=method))+geom_bar(stat="identity")+guides(fill=F)+
  scale_fill_manual(values = c("#ee7c6e","#f4b062","#b37cb3","#add26a","#CBEBC5","#E6D77A","#7eabca","#d9d9d9"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
p1+p2
