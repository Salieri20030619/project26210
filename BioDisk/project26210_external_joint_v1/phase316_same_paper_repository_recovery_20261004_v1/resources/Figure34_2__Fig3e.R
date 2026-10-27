source("/Path/summ_dist_V4_statistic.R")
CRC3a_st<-readRDS("/Path/data/CRC3a_distance_Seurat.RDS")
coord <- CRC3a_st@images$sliceD1@coordinates
stromal_name_list <- c("CD8.Tem","CD8.Tex","CD8.Trm") %>% as.list()
y="CD8.Tem"

dist_stromal <- lapply(stromal_name_list, function(y){
  weight.matr <- dplyr::select(CRC3a_st@meta.data,y,Epithelial_sum)
  dist <- summ_dist_V2(object = CRC3a_st,
                       barcode_set.2 = colnames(CRC3a_st)[which(CRC3a_st$Region == "Clone1")],
                       barcode_set.1 = colnames(CRC3a_st)[which(CRC3a_st@meta.data[[y]] > 0.03 )],
                       coord,weight = 'double',weight.matr = weight.matr)
  dist
})
dist_stromal



Tem<-data.frame(celltype="CD8+ Tem",dis=dist_stromal[[1]][[1]] * dist_stromal[[1]][[2]])
Tex<-data.frame(celltype="CD8+ Tex",dis=dist_stromal[[2]][[1]] * dist_stromal[[2]][[2]])
Trm<-data.frame(celltype="CD8+ Trm",dis=dist_stromal[[3]][[1]] * dist_stromal[[3]][[2]])
T_distribution<-rbind(Tem,Tex,Trm)
my_comparisons <- list( c("CD8+ Tem", "CD8+ Tex"), c("CD8+ Tex", "CD8+ Trm"), c("CD8+ Tem", "CD8+ Trm") )
p3<-ggplot(T_distribution[T_distribution$dis >= 0,],aes(celltype,dis,color=celltype))+geom_boxplot(outliers = FALSE)+
  scale_color_manual(values = c("#e6513a","#4a91c1","#7066a7"))+labs(x="Subpopulations",y="Weighted Distance to tumor")+
  guides(color=F)
p3
plot_data<-CRC3a_st@meta.data[CRC3a_st@meta.data$Region== "Clone1",c("CD8.Tem","CD8.Tex","CD8.Trm")]

plot_data<-reshape2::melt(plot_data)
p4<-ggplot(plot_data,aes(x=variable,y=value,fill=variable))+stat_summary(fun.y="mean",geom="bar",size=3) + 
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.05)+scale_fill_manual(values = c("#e6513a","#4a91c1","#7066a7"))+guides(fill=F)+
  labs(x="Celltype",y="Proportion (%)")+coord_cartesian(ylim = c(0,0.02))+stat_compare_means(ref.group = "CD8.Tem",method = "t.test",label = "p.signif",label.y = 0.02)
p1+p2+p3+p4
