total_sample_porp_immune<-readRDS("/Path/data/total_sample_porp_immune.RDS")
total_sample_porp_immune$Barcode<-substr(rownames(total_sample_porp_immune), 1, 18)

CRC4_st@meta.data<-merge(CRC4_st@meta.data,total_sample_porp_immune[total_sample_porp_immune$slice == "CRC4",c("Barcode","community1","community2","community3")],by="Barcode")
rownames(CRC4_st@meta.data)<-CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
CRC4_st<-UpdateSeuratObject(CRC4_st)

colnames(CRC4_st@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CRC4_st@meta.data))))))))))))
colnames(CRC4_st@meta.data)<-gsub("Fibroblast","CD34+ iFib",colnames(CRC4_st@meta.data))
                                  
CRC4_st@meta.data[,c("community1","community2","community3")][is.na(CRC4_st@meta.data[,c("community1","community2","community3")])]<-0


CRC4_st@meta.data[CRC4_st@meta.data$Region == "Remaining",]<-NA

p1<-SpatialPlot(CRC4_st,features = c("community1"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "community1",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p2<-SpatialPlot(CRC4_st,features = c("CD34+ iFib"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD34+ iFib",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p3<-SpatialPlot(CRC4_st,features = c("CD8+ Tem"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8+ Tem",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p4<-SpatialPlot(CRC4_st,features = c("CD8+ Tex"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8+ Tex",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p5<-SpatialPlot(CRC4_st,features = c("C1QC+ Mac"),image.alpha = 0,stroke = NA,ncol = 1,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "C1QC+ Mac",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))

p1+p5+p2+p4+p3+plot_layout(ncol=5)
