



total_sample_porp_immune<-readRDS("/Path/data/total_sample_porp_immune.RDS")
total_sample_porp_immune$Barcode<-substr(rownames(total_sample_porp_immune), 1, 18)

CRC4_st@meta.data<-merge(CRC4_st@meta.data,total_sample_porp_immune[total_sample_porp_immune$slice == "CRC4",c("Barcode","community1","community2","community3")],by="Barcode")
rownames(CRC4_st@meta.data)<-CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
CRC4_st<-UpdateSeuratObject(CRC4_st)

plot_data<-CRC4_st@meta.data[,c("Region","community1","community2","community3")]
theme_set(theme_cowplot())
p1<-ggplot(plot_data[plot_data$Region != "Remaining",],aes(x=Region,y=log2(community1),fill=Region,color=Region))+geom_boxplot()+
  scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+scale_color_manual(values = c("#378B39","#1F6C9F","#B62D28"))+
  labs(x="",y="log2(Community Score)",title = "Community1")+theme(plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)+
  stat_compare_means(ref.group = "Clone1",method = "wilcox.test",label = "p.signif")
p2<-ggplot(plot_data[plot_data$Region != "Remaining",],aes(x=Region,y=log2(community2),fill=Region,color=Region))+geom_boxplot()+
  scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+scale_color_manual(values = c("#378B39","#1F6C9F","#B62D28"))+
  labs(x="",y="log2(Community Score)",title = "Community2")+theme(plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)+
  stat_compare_means(ref.group = "Clone1",method = "wilcox.test",label = "p.signif")
p3<-ggplot(plot_data[plot_data$Region != "Remaining",],aes(x=Region,y=log2(community3),fill=Region,color=Region))+geom_boxplot()+
  scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+scale_color_manual(values = c("#378B39","#1F6C9F","#B62D28"))+
  labs(x="",y="log2(Community Score)",title = "Community3")+theme(plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)+
  stat_compare_means(ref.group = "Clone1",method = "wilcox.test",label = "p.signif")
library("patchwork")
p1+p2+p3+plot_layout(ncol=2)

CRC4_st_tumor<-CRC4_st[,CRC4_st$Node != "Remaining"]
SpatialPlot(CRC4_st_tumor,group.by ="Node",stroke = F,pt.size=1.2,cols = c("#A7CC83","#9BBFD0","#ED7B6C"))