#distance calculation
theme_set(theme_cowplot())
dis1<-read.csv("/Path/data/distance_data1.txt",sep="\t")
dis2<-read.csv("/Path/data/distance_data2.txt",sep="\t")
dis3<-read.csv("/Path/data/distance_data3.txt",sep="\t")
dis4<-read.csv("/Path/data/distance_data4.txt",sep="\t")
dis5<-read.csv("/Path/data/distance_data5.txt",sep="\t")
dis6<-read.csv("/Path/data/distance_data6.txt",sep="\t")
dis7<-read.csv("/Path/data/distance_data7.txt",sep="\t")
dis8<-read.csv("/Path/data/distance_data8.txt",sep="\t")
dis9<-read.csv("/Path/data/distance_data9.txt",sep="\t")
dis10<-read.csv("/Path/data/distance_data10.txt",sep="\t")
distance_combine<-rbind(dis1,dis2,dis3,dis4,dis5,dis6,dis7,dis8,dis9,dis10)
distance_combine$Signed.distance.to.annotation.Tumor.µm[distance_combine$Signed.distance.to.annotation.Tumor.µm<0] = 0
CD8_name <- names(table(distance_combine$Classification))
CD8_name<-CD8_name[grep("CD8",CD8_name)]
distance_combine$Celltype<-"Others"
distance_combine$Celltype[distance_combine$Classification %in% c("CD8: CD103: PD-1","CD8: PD-1","GZMK: CD8: CD103: PD-1","GZMK: CD8: PD-1")]<-"CD8 Tex"
distance_combine$Celltype[distance_combine$Classification %in% c("CD8: CD103","GZMK: CD8: CD103")]<-"CD8 Trm"
distance_combine$Celltype[distance_combine$Classification %in% c("GZMK: CD8")]<-"CD8 Tem"
my_comparisons <- list( c("CD8 Tem", "CD8 Tex"), c("CD8 Tex", "CD8 Trm"), c("CD8 Tem", "CD8 Trm") )
p1<-ggplot(distance_combine[distance_combine$Celltype != "Others",],
       aes(x=Celltype,y=Signed.distance.to.annotation.Tumor.µm,color=Celltype))+geom_boxplot(outliers = F)+
  labs(x="Subpopulation",y="Distance to tumor (um)",title="mIHC(N=10)")+scale_color_manual(values = c("#cf4e39","#4888b5","#675e97"))+
  guides(color=F)+theme(plot.title = element_text(hjust=0.5))
p1
table(distance_combine$Celltype,distance_combine$Signed.distance.to.annotation.Tumor.µm == 0)
proportion<-data.frame(cell=c("CD8 Tem", "CD8 Tex","CD8 Trm"),value=c(122/37216 * 100,294/37216 * 100,396/37216 * 100))
p2<-ggplot(proportion,aes(x=cell,y=value,fill=cell))+geom_bar(stat = "identity",position = "dodge")+
  labs(x="Subpopulation",y="Proportion (%)",title="mIHC(N=10)")+scale_fill_manual(values = c("#cf4e39","#4888b5","#675e97"))+
  guides(fill=F)+theme(plot.title = element_text(hjust=0.5))+coord_cartesian(ylim = c(0,1.5))
p2
p1+p2
