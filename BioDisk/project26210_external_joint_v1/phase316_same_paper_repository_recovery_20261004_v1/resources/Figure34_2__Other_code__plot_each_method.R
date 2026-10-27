# plot each method TME (total celltype)
setwd("/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/")
celltype_all<-colnames(C2L[[1]]@meta.data)[4:37]

CRC1_st<-readRDS("/data/huangzr/Spatial/Kras/version5/Region_devide_2/clone_region/CRC1_st_clone.RDS")
CRC1_st@meta.data[,c(4:37,40)]

#prepare clonal region
CRC1_Region<-CRC1_st@meta.data[,c("Barcode","Region")]
CRC2_Region<-CRC2_st@meta.data[,c("Barcode","Region")]
CRC3a_Region<-CRC3a_st@meta.data[,c("Barcode","Region")]
CRC3b_Region<-CRC3b_st@meta.data[,c("Barcode","Region")]
CRC4_Region<-CRC4_st@meta.data[,c("Barcode","Region")]
CRC5_Region<-CRC5_st@meta.data[,c("Barcode","Region")]
CRC6_Region<-CRC6_st@meta.data[,c("Barcode","Region")]
CRC7_Region<-CRC7_st@meta.data[,c("Barcode","Region")]
CRC4_st@meta.data[,c("Barcode","Region")]

#clone region annotation
RCTD[[1]]@meta.data<-merge(RCTD[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(RCTD[[1]]@meta.data)<-RCTD[[1]]@meta.data$Barcode
RCTD[[1]]@meta.data<-RCTD[[1]]@meta.data[colnames(RCTD[[1]]),]

RCTD[[2]]@meta.data<-merge(RCTD[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(RCTD[[2]]@meta.data)<-RCTD[[2]]@meta.data$Barcode
RCTD[[2]]@meta.data<-RCTD[[2]]@meta.data[colnames(RCTD[[2]]),]

RCTD[[3]]@meta.data<-merge(RCTD[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(RCTD[[3]]@meta.data)<-RCTD[[3]]@meta.data$Barcode
RCTD[[3]]@meta.data<-RCTD[[3]]@meta.data[colnames(RCTD[[3]]),]

RCTD[[4]]@meta.data<-merge(RCTD[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(RCTD[[4]]@meta.data)<-RCTD[[4]]@meta.data$Barcode
RCTD[[4]]@meta.data<-RCTD[[4]]@meta.data[colnames(RCTD[[4]]),]

RCTD[[5]]@meta.data<-merge(RCTD[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(RCTD[[5]]@meta.data)<-RCTD[[5]]@meta.data$Barcode
RCTD[[5]]@meta.data<-RCTD[[5]]@meta.data[colnames(RCTD[[5]]),]

colnames(RCTD[[6]]@meta.data)[1]<-"Barcode"
RCTD[[6]]@meta.data<-merge(RCTD[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(RCTD[[6]]@meta.data)<-RCTD[[6]]@meta.data$Barcode
RCTD[[6]]@meta.data<-RCTD[[6]]@meta.data[colnames(RCTD[[6]]),]

RCTD[[7]]@meta.data<-merge(RCTD[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(RCTD[[7]]@meta.data)<-RCTD[[7]]@meta.data$Barcode
RCTD[[7]]@meta.data<-RCTD[[7]]@meta.data[colnames(RCTD[[7]]),]

RCTD[[8]]@meta.data<-merge(RCTD[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(RCTD[[8]]@meta.data)<-RCTD[[8]]@meta.data$Barcode
RCTD[[8]]@meta.data<-RCTD[[8]]@meta.data[colnames(RCTD[[8]]),]


inhouse_cohort<-rbind(RCTD[[1]]@meta.data[,c(4:38)],RCTD[[4]]@meta.data[,c(4:38)],RCTD[[2]]@meta.data[,c(4:38)],RCTD[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(RCTD[[6]]@meta.data[,c(4:38)],RCTD[[5]]@meta.data[,c(4:38)],RCTD[[7]]@meta.data[,c(4:38)],RCTD[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(RCTD[[1]]@meta.data)),rep("CRC2",nrow(RCTD[[4]]@meta.data)),rep("CRC3a",nrow(RCTD[[2]]@meta.data)),rep("CRC3b",nrow(RCTD[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(RCTD[[6]]@meta.data)),rep("CRC5",nrow(RCTD[[5]]@meta.data)),rep("CRC6",nrow(RCTD[[7]]@meta.data)),rep("CRC7",nrow(RCTD[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.02的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)
tumor_prop_avg_RCTD<-tumor_prop_avg
saveRDS(tumor_prop_avg_RCTD,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_RCTD.RDS")
my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"

#Roe calculating
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"

library(dplyr)
sele
total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)

Roe_celltype_list<-celltype_all[c(-11,-18,-21,-29,-33)]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  dplyr::select(Roe_celltype_list) %>% t() %>% as.data.frame()

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
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

p1<-ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+labs(x="Ro/e",y="Subset",title = "RCTD")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/RCTD_ROE.pdf",width = 5,height=8)

levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
p1<-pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                   gaps_col = c(3,8),scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                   annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8,main = "RCTD")
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/RCTD_heatmap.pdf",width = 5,height=8)

long_df <- pivot_longer(tumor_prop_avg, cols = 1:34, names_to = "Variable", values_to = "Value")
var=celltype_all[2]
plot_list <- lapply(celltype_all, function(var) {
  ggplot(subset(long_df, Variable == var),aes(x=subTME_category,y=Value,fill=subTME_category,color=subTME_category))+geom_boxplot(size=0.6)+
    scale_fill_manual(values = c("#A6CEE3","#FDBF6F"))+scale_color_manual(values = c("#418DBF","#FD8B18"))+
    labs(x="",y="Proportion",title = var)+
    stat_compare_means(label="p.signif",ref.group = "T cell-excluded")+
    theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)
})
ggsave(plot = do.call(gridExtra::grid.arrange, c(plot_list, ncol = 6)),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/RCTD_box.pdf",width = 20,height=40)



#C2L
C2L<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/C2L_results.rds")
colnames(C2L[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[1]]@meta.data))))))))))))
colnames(C2L[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[2]]@meta.data))))))))))))
colnames(C2L[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[3]]@meta.data))))))))))))
colnames(C2L[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[4]]@meta.data))))))))))))
colnames(C2L[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[5]]@meta.data))))))))))))
colnames(C2L[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[6]]@meta.data))))))))))))
colnames(C2L[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[7]]@meta.data))))))))))))
colnames(C2L[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[8]]@meta.data))))))))))))

SpatialFeaturePlot(C2L[[6]],features = celltype_all[1])

#clone region annotation
C2L[[1]]@meta.data<-merge(C2L[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(C2L[[1]]@meta.data)<-C2L[[1]]@meta.data$Barcode
C2L[[1]]@meta.data<-C2L[[1]]@meta.data[colnames(C2L[[1]]),]

C2L[[2]]@meta.data<-merge(C2L[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(C2L[[2]]@meta.data)<-C2L[[2]]@meta.data$Barcode
C2L[[2]]@meta.data<-C2L[[2]]@meta.data[colnames(C2L[[2]]),]

C2L[[3]]@meta.data<-merge(C2L[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(C2L[[3]]@meta.data)<-C2L[[3]]@meta.data$Barcode
C2L[[3]]@meta.data<-C2L[[3]]@meta.data[colnames(C2L[[3]]),]

C2L[[4]]@meta.data<-merge(C2L[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(C2L[[4]]@meta.data)<-C2L[[4]]@meta.data$Barcode
C2L[[4]]@meta.data<-C2L[[4]]@meta.data[colnames(C2L[[4]]),]

C2L[[5]]@meta.data<-merge(C2L[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(C2L[[5]]@meta.data)<-C2L[[5]]@meta.data$Barcode
C2L[[5]]@meta.data<-C2L[[5]]@meta.data[colnames(C2L[[5]]),]

colnames(C2L[[6]]@meta.data)[1]<-"Barcode"
C2L[[6]]@meta.data<-merge(C2L[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(C2L[[6]]@meta.data)<-C2L[[6]]@meta.data$Barcode
C2L[[6]]@meta.data<-C2L[[6]]@meta.data[colnames(C2L[[6]]),]

C2L[[7]]@meta.data<-merge(C2L[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(C2L[[7]]@meta.data)<-C2L[[7]]@meta.data$Barcode
C2L[[7]]@meta.data<-C2L[[7]]@meta.data[colnames(C2L[[7]]),]

C2L[[8]]@meta.data<-merge(C2L[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(C2L[[8]]@meta.data)<-C2L[[8]]@meta.data$Barcode
C2L[[8]]@meta.data<-C2L[[8]]@meta.data[colnames(C2L[[8]]),]


inhouse_cohort<-rbind(C2L[[1]]@meta.data[,c(4:38)],C2L[[4]]@meta.data[,c(4:38)],C2L[[2]]@meta.data[,c(4:38)],C2L[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(C2L[[6]]@meta.data[,c(4:38)],C2L[[5]]@meta.data[,c(4:38)],C2L[[7]]@meta.data[,c(4:38)],C2L[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(C2L[[1]]@meta.data)),rep("CRC2",nrow(C2L[[4]]@meta.data)),rep("CRC3a",nrow(C2L[[2]]@meta.data)),rep("CRC3b",nrow(C2L[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(C2L[[6]]@meta.data)),rep("CRC5",nrow(C2L[[5]]@meta.data)),rep("CRC6",nrow(C2L[[7]]@meta.data)),rep("CRC7",nrow(C2L[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.02的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)

tumor_prop_avg_C2L<-tumor_prop_avg
saveRDS(tumor_prop_avg_C2L,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_C2L.RDS")

my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"

#Roe calculating
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"

library(dplyr)
sele
total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)

Roe_celltype_list<-celltype_all[c(-11,-18,-21,-29,-33)]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  dplyr::select(Roe_celltype_list) %>% t() %>% as.data.frame()

library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
Roe_v_inhouse[is.na(Roe_v_inhouse)]<-0
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

p1<-ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+labs(x="Ro/e",y="Subset",title = "C2L")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/C2L_ROE.pdf",width = 5,height=8)

levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
p1<-pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                       gaps_col = c(3,8),scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                       annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8,main = "C2L")
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/C2L_heatmap.pdf",width = 5,height=8)

long_df <- pivot_longer(tumor_prop_avg, cols = 1:34, names_to = "Variable", values_to = "Value")

var=celltype_all[2]
plot_list <- lapply(celltype_all, function(var) {
  ggplot(subset(long_df, Variable == var),aes(x=subTME_category,y=Value,fill=subTME_category,color=subTME_category))+geom_boxplot(size=0.6)+
    scale_fill_manual(values = c("#A6CEE3","#FDBF6F"))+scale_color_manual(values = c("#418DBF","#FD8B18"))+
    labs(x="",y="Proportion",title = var)+
    stat_compare_means(label="p.signif",ref.group = "T cell-excluded")+
    theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)
})
ggsave(plot = do.call(gridExtra::grid.arrange, c(plot_list, ncol = 6)),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/C2L_box.pdf",width = 20,height=40)

#stereoscope
stereoscope<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/stereoscope_results.rds")
colnames(stereoscope[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[1]]@meta.data))))))))))))
colnames(stereoscope[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[2]]@meta.data))))))))))))
colnames(stereoscope[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[3]]@meta.data))))))))))))
colnames(stereoscope[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[4]]@meta.data))))))))))))
colnames(stereoscope[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[5]]@meta.data))))))))))))
colnames(stereoscope[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[6]]@meta.data))))))))))))
colnames(stereoscope[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[7]]@meta.data))))))))))))
colnames(stereoscope[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[8]]@meta.data))))))))))))

#SpatialFeaturePlot(stereoscope[[6]],features = celltype_all[1])

#clone region annotation
stereoscope[[1]]@meta.data<-merge(stereoscope[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(stereoscope[[1]]@meta.data)<-stereoscope[[1]]@meta.data$Barcode
stereoscope[[1]]@meta.data<-stereoscope[[1]]@meta.data[colnames(stereoscope[[1]]),]

stereoscope[[2]]@meta.data<-merge(stereoscope[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(stereoscope[[2]]@meta.data)<-stereoscope[[2]]@meta.data$Barcode
stereoscope[[2]]@meta.data<-stereoscope[[2]]@meta.data[colnames(stereoscope[[2]]),]

stereoscope[[3]]@meta.data<-merge(stereoscope[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(stereoscope[[3]]@meta.data)<-stereoscope[[3]]@meta.data$Barcode
stereoscope[[3]]@meta.data<-stereoscope[[3]]@meta.data[colnames(stereoscope[[3]]),]

stereoscope[[4]]@meta.data<-merge(stereoscope[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(stereoscope[[4]]@meta.data)<-stereoscope[[4]]@meta.data$Barcode
stereoscope[[4]]@meta.data<-stereoscope[[4]]@meta.data[colnames(stereoscope[[4]]),]

stereoscope[[5]]@meta.data<-merge(stereoscope[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(stereoscope[[5]]@meta.data)<-stereoscope[[5]]@meta.data$Barcode
stereoscope[[5]]@meta.data<-stereoscope[[5]]@meta.data[colnames(stereoscope[[5]]),]

colnames(stereoscope[[6]]@meta.data)[1]<-"Barcode"
stereoscope[[6]]@meta.data<-merge(stereoscope[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(stereoscope[[6]]@meta.data)<-stereoscope[[6]]@meta.data$Barcode
stereoscope[[6]]@meta.data<-stereoscope[[6]]@meta.data[colnames(stereoscope[[6]]),]

stereoscope[[7]]@meta.data<-merge(stereoscope[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(stereoscope[[7]]@meta.data)<-stereoscope[[7]]@meta.data$Barcode
stereoscope[[7]]@meta.data<-stereoscope[[7]]@meta.data[colnames(stereoscope[[7]]),]

stereoscope[[8]]@meta.data<-merge(stereoscope[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(stereoscope[[8]]@meta.data)<-stereoscope[[8]]@meta.data$Barcode
stereoscope[[8]]@meta.data<-stereoscope[[8]]@meta.data[colnames(stereoscope[[8]]),]


inhouse_cohort<-rbind(stereoscope[[1]]@meta.data[,c(4:38)],stereoscope[[4]]@meta.data[,c(4:38)],stereoscope[[2]]@meta.data[,c(4:38)],stereoscope[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(stereoscope[[6]]@meta.data[,c(4:38)],stereoscope[[5]]@meta.data[,c(4:38)],stereoscope[[7]]@meta.data[,c(4:38)],stereoscope[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(stereoscope[[1]]@meta.data)),rep("CRC2",nrow(stereoscope[[4]]@meta.data)),rep("CRC3a",nrow(stereoscope[[2]]@meta.data)),rep("CRC3b",nrow(stereoscope[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(stereoscope[[6]]@meta.data)),rep("CRC5",nrow(stereoscope[[5]]@meta.data)),rep("CRC6",nrow(stereoscope[[7]]@meta.data)),rep("CRC7",nrow(stereoscope[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.02的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)

tumor_prop_avg_steroscope<-tumor_prop_avg
saveRDS(tumor_prop_avg_steroscope,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_steroscope.RDS")

my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"

#Roe calculating
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"

library(dplyr)
sele
total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)

Roe_celltype_list<-celltype_all[c(-11,-18,-21,-29,-33)]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  dplyr::select(Roe_celltype_list) %>% t() %>% as.data.frame()

library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
Roe_v_inhouse[is.na(Roe_v_inhouse)]<-0
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

p1<-ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+labs(x="Ro/e",y="Subset",title = "stereoscope")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/stereoscope_ROE.pdf",width = 5,height=8)

levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
p1<-pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                       gaps_col = c(3,8),scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                       annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8,main = "stereoscope")
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/stereoscope_heatmap.pdf",width = 5,height=8)

long_df <- pivot_longer(tumor_prop_avg, cols = 1:34, names_to = "Variable", values_to = "Value")

var=celltype_all[2]
plot_list <- lapply(celltype_all, function(var) {
  ggplot(subset(long_df, Variable == var),aes(x=subTME_category,y=Value,fill=subTME_category,color=subTME_category))+geom_boxplot(size=0.6)+
    scale_fill_manual(values = c("#A6CEE3","#FDBF6F"))+scale_color_manual(values = c("#418DBF","#FD8B18"))+
    labs(x="",y="Proportion",title = var)+
    stat_compare_means(label="p.signif",ref.group = "T cell-excluded")+
    theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)
})
ggsave(plot = do.call(gridExtra::grid.arrange, c(plot_list, ncol = 6)),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/stereoscope_box.pdf",width = 20,height=40)


#CARD
CARD<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/CARD_results.rds")
colnames(CARD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[1]]@meta.data))))))))))))
colnames(CARD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[2]]@meta.data))))))))))))
colnames(CARD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[3]]@meta.data))))))))))))
colnames(CARD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[4]]@meta.data))))))))))))
colnames(CARD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[5]]@meta.data))))))))))))
colnames(CARD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[6]]@meta.data))))))))))))
colnames(CARD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[7]]@meta.data))))))))))))
colnames(CARD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[8]]@meta.data))))))))))))

#SpatialFeaturePlot(CARD[[6]],features = celltype_all[1])

#clone region annotation
CARD[[1]]@meta.data<-merge(CARD[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(CARD[[1]]@meta.data)<-CARD[[1]]@meta.data$Barcode
CARD[[1]]@meta.data<-CARD[[1]]@meta.data[colnames(CARD[[1]]),]

CARD[[2]]@meta.data<-merge(CARD[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(CARD[[2]]@meta.data)<-CARD[[2]]@meta.data$Barcode
CARD[[2]]@meta.data<-CARD[[2]]@meta.data[colnames(CARD[[2]]),]

CARD[[3]]@meta.data<-merge(CARD[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(CARD[[3]]@meta.data)<-CARD[[3]]@meta.data$Barcode
CARD[[3]]@meta.data<-CARD[[3]]@meta.data[colnames(CARD[[3]]),]

CARD[[4]]@meta.data<-merge(CARD[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(CARD[[4]]@meta.data)<-CARD[[4]]@meta.data$Barcode
CARD[[4]]@meta.data<-CARD[[4]]@meta.data[colnames(CARD[[4]]),]

CARD[[5]]@meta.data<-merge(CARD[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(CARD[[5]]@meta.data)<-CARD[[5]]@meta.data$Barcode
CARD[[5]]@meta.data<-CARD[[5]]@meta.data[colnames(CARD[[5]]),]

colnames(CARD[[6]]@meta.data)[1]<-"Barcode"
CARD[[6]]@meta.data<-merge(CARD[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(CARD[[6]]@meta.data)<-CARD[[6]]@meta.data$Barcode
CARD[[6]]@meta.data<-CARD[[6]]@meta.data[colnames(CARD[[6]]),]

CARD[[7]]@meta.data<-merge(CARD[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(CARD[[7]]@meta.data)<-CARD[[7]]@meta.data$Barcode
CARD[[7]]@meta.data<-CARD[[7]]@meta.data[colnames(CARD[[7]]),]

CARD[[8]]@meta.data<-merge(CARD[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(CARD[[8]]@meta.data)<-CARD[[8]]@meta.data$Barcode
CARD[[8]]@meta.data<-CARD[[8]]@meta.data[colnames(CARD[[8]]),]


inhouse_cohort<-rbind(CARD[[1]]@meta.data[,c(4:38)],CARD[[4]]@meta.data[,c(4:38)],CARD[[2]]@meta.data[,c(4:38)],CARD[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(CARD[[6]]@meta.data[,c(4:38)],CARD[[5]]@meta.data[,c(4:38)],CARD[[7]]@meta.data[,c(4:38)],CARD[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(CARD[[1]]@meta.data)),rep("CRC2",nrow(CARD[[4]]@meta.data)),rep("CRC3a",nrow(CARD[[2]]@meta.data)),rep("CRC3b",nrow(CARD[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(CARD[[6]]@meta.data)),rep("CRC5",nrow(CARD[[5]]@meta.data)),rep("CRC6",nrow(CARD[[7]]@meta.data)),rep("CRC7",nrow(CARD[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.02的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)

tumor_prop_avg_CARD<-tumor_prop_avg
saveRDS(tumor_prop_avg_CARD,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_CARD.RDS")

my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"

#Roe calculating
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"

library(dplyr)

total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)

Roe_celltype_list<-celltype_all[c(-11,-18,-21,-29,-33)]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  dplyr::select(Roe_celltype_list) %>% t() %>% as.data.frame()

library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
Roe_v_inhouse[is.na(Roe_v_inhouse)]<-0
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

p1<-ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+labs(x="Ro/e",y="Subset",title = "CARD")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/CARD_ROE.pdf",width = 5,height=8)

levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
p1<-pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                       gaps_col = c(3,8),scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                       annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8,main = "CARD")
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/CARD_heatmap.pdf",width = 5,height=8)

long_df <- pivot_longer(tumor_prop_avg, cols = 1:34, names_to = "Variable", values_to = "Value")

var=celltype_all[2]
plot_list <- lapply(celltype_all, function(var) {
  ggplot(subset(long_df, Variable == var),aes(x=subTME_category,y=Value,fill=subTME_category,color=subTME_category))+geom_boxplot(size=0.6)+
    scale_fill_manual(values = c("#A6CEE3","#FDBF6F"))+scale_color_manual(values = c("#418DBF","#FD8B18"))+
    labs(x="",y="Proportion",title = var)+
    stat_compare_means(label="p.signif",ref.group = "T cell-excluded")+
    theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)
})
ggsave(plot = do.call(gridExtra::grid.arrange, c(plot_list, ncol = 6)),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/CARD_box.pdf",width = 20,height=40)

#SPOTlight
SPOTlight<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/SPOTlight_results.rds")
colnames(SPOTlight[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[1]]@meta.data))))))))))))
colnames(SPOTlight[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[2]]@meta.data))))))))))))
colnames(SPOTlight[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[3]]@meta.data))))))))))))
colnames(SPOTlight[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[4]]@meta.data))))))))))))
colnames(SPOTlight[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[5]]@meta.data))))))))))))
colnames(SPOTlight[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[6]]@meta.data))))))))))))
colnames(SPOTlight[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[7]]@meta.data))))))))))))
colnames(SPOTlight[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[8]]@meta.data))))))))))))

#SpatialFeaturePlot(SPOTlight[[6]],features = celltype_all[1])

#clone region annotation
SPOTlight[[1]]@meta.data<-merge(SPOTlight[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(SPOTlight[[1]]@meta.data)<-SPOTlight[[1]]@meta.data$Barcode
SPOTlight[[1]]@meta.data<-SPOTlight[[1]]@meta.data[colnames(SPOTlight[[1]]),]

SPOTlight[[2]]@meta.data<-merge(SPOTlight[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(SPOTlight[[2]]@meta.data)<-SPOTlight[[2]]@meta.data$Barcode
SPOTlight[[2]]@meta.data<-SPOTlight[[2]]@meta.data[colnames(SPOTlight[[2]]),]

SPOTlight[[3]]@meta.data<-merge(SPOTlight[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(SPOTlight[[3]]@meta.data)<-SPOTlight[[3]]@meta.data$Barcode
SPOTlight[[3]]@meta.data<-SPOTlight[[3]]@meta.data[colnames(SPOTlight[[3]]),]

SPOTlight[[4]]@meta.data<-merge(SPOTlight[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(SPOTlight[[4]]@meta.data)<-SPOTlight[[4]]@meta.data$Barcode
SPOTlight[[4]]@meta.data<-SPOTlight[[4]]@meta.data[colnames(SPOTlight[[4]]),]

SPOTlight[[5]]@meta.data<-merge(SPOTlight[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(SPOTlight[[5]]@meta.data)<-SPOTlight[[5]]@meta.data$Barcode
SPOTlight[[5]]@meta.data<-SPOTlight[[5]]@meta.data[colnames(SPOTlight[[5]]),]

colnames(SPOTlight[[6]]@meta.data)[1]<-"Barcode"
SPOTlight[[6]]@meta.data<-merge(SPOTlight[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(SPOTlight[[6]]@meta.data)<-SPOTlight[[6]]@meta.data$Barcode
SPOTlight[[6]]@meta.data<-SPOTlight[[6]]@meta.data[colnames(SPOTlight[[6]]),]

SPOTlight[[7]]@meta.data<-merge(SPOTlight[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(SPOTlight[[7]]@meta.data)<-SPOTlight[[7]]@meta.data$Barcode
SPOTlight[[7]]@meta.data<-SPOTlight[[7]]@meta.data[colnames(SPOTlight[[7]]),]

SPOTlight[[8]]@meta.data<-merge(SPOTlight[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(SPOTlight[[8]]@meta.data)<-SPOTlight[[8]]@meta.data$Barcode
SPOTlight[[8]]@meta.data<-SPOTlight[[8]]@meta.data[colnames(SPOTlight[[8]]),]


inhouse_cohort<-rbind(SPOTlight[[1]]@meta.data[,c(4:38)],SPOTlight[[4]]@meta.data[,c(4:38)],SPOTlight[[2]]@meta.data[,c(4:38)],SPOTlight[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(SPOTlight[[6]]@meta.data[,c(4:38)],SPOTlight[[5]]@meta.data[,c(4:38)],SPOTlight[[7]]@meta.data[,c(4:38)],SPOTlight[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(SPOTlight[[1]]@meta.data)),rep("CRC2",nrow(SPOTlight[[4]]@meta.data)),rep("CRC3a",nrow(SPOTlight[[2]]@meta.data)),rep("CRC3b",nrow(SPOTlight[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(SPOTlight[[6]]@meta.data)),rep("CRC5",nrow(SPOTlight[[5]]@meta.data)),rep("CRC6",nrow(SPOTlight[[7]]@meta.data)),rep("CRC7",nrow(SPOTlight[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.03的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)

tumor_prop_avg_SPOTlight<-tumor_prop_avg
saveRDS(tumor_prop_avg_SPOTlight,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_SPOTlight.RDS")

my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"

#Roe calculating
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"

library(dplyr)

total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)

Roe_celltype_list<-celltype_all[c(-11,-18,-21,-29,-33)]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  dplyr::select(Roe_celltype_list) %>% t() %>% as.data.frame()

library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
Roe_v_inhouse[is.na(Roe_v_inhouse)]<-0
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

p1<-ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+labs(x="Ro/e",y="Subset",title = "SPOTlight")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/SPOTlight_ROE.pdf",width = 5,height=8)

levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
p1<-pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                       gaps_col = c(3,8),scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                       annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8,main = "SPOTlight")
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/SPOTlight_heatmap.pdf",width = 5,height=8)

long_df <- pivot_longer(tumor_prop_avg, cols = 1:34, names_to = "Variable", values_to = "Value")

var=celltype_all[2]
plot_list <- lapply(celltype_all, function(var) {
  ggplot(subset(long_df, Variable == var),aes(x=subTME_category,y=Value,fill=subTME_category,color=subTME_category))+geom_boxplot(size=0.6)+
    scale_fill_manual(values = c("#A6CEE3","#FDBF6F"))+scale_color_manual(values = c("#418DBF","#FD8B18"))+
    labs(x="",y="Proportion",title = var)+
    stat_compare_means(label="p.signif",ref.group = "T cell-excluded")+
    theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)
})
ggsave(plot = do.call(gridExtra::grid.arrange, c(plot_list, ncol = 6)),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/SPOTlight_box.pdf",width = 20,height=40)

#UCASpatial
UCASpatial<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/UCASpatial_results.rds")
colnames(UCASpatial[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[1]]@meta.data))))))))))))
colnames(UCASpatial[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[2]]@meta.data))))))))))))
colnames(UCASpatial[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[3]]@meta.data))))))))))))
colnames(UCASpatial[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[4]]@meta.data))))))))))))
colnames(UCASpatial[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[5]]@meta.data))))))))))))
colnames(UCASpatial[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[6]]@meta.data))))))))))))
colnames(UCASpatial[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[7]]@meta.data))))))))))))
colnames(UCASpatial[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[8]]@meta.data))))))))))))

#SpatialFeaturePlot(UCASpatial[[6]],features = celltype_all[1])

#clone region annotation
UCASpatial[[1]]@meta.data<-merge(UCASpatial[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(UCASpatial[[1]]@meta.data)<-UCASpatial[[1]]@meta.data$Barcode
UCASpatial[[1]]@meta.data<-UCASpatial[[1]]@meta.data[colnames(UCASpatial[[1]]),]

UCASpatial[[2]]@meta.data<-merge(UCASpatial[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(UCASpatial[[2]]@meta.data)<-UCASpatial[[2]]@meta.data$Barcode
UCASpatial[[2]]@meta.data<-UCASpatial[[2]]@meta.data[colnames(UCASpatial[[2]]),]

UCASpatial[[3]]@meta.data<-merge(UCASpatial[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(UCASpatial[[3]]@meta.data)<-UCASpatial[[3]]@meta.data$Barcode
UCASpatial[[3]]@meta.data<-UCASpatial[[3]]@meta.data[colnames(UCASpatial[[3]]),]

UCASpatial[[4]]@meta.data<-merge(UCASpatial[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(UCASpatial[[4]]@meta.data)<-UCASpatial[[4]]@meta.data$Barcode
UCASpatial[[4]]@meta.data<-UCASpatial[[4]]@meta.data[colnames(UCASpatial[[4]]),]

UCASpatial[[5]]@meta.data<-merge(UCASpatial[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(UCASpatial[[5]]@meta.data)<-UCASpatial[[5]]@meta.data$Barcode
UCASpatial[[5]]@meta.data<-UCASpatial[[5]]@meta.data[colnames(UCASpatial[[5]]),]

colnames(UCASpatial[[6]]@meta.data)[1]<-"Barcode"
UCASpatial[[6]]@meta.data<-merge(UCASpatial[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(UCASpatial[[6]]@meta.data)<-UCASpatial[[6]]@meta.data$Barcode
UCASpatial[[6]]@meta.data<-UCASpatial[[6]]@meta.data[colnames(UCASpatial[[6]]),]

UCASpatial[[7]]@meta.data<-merge(UCASpatial[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(UCASpatial[[7]]@meta.data)<-UCASpatial[[7]]@meta.data$Barcode
UCASpatial[[7]]@meta.data<-UCASpatial[[7]]@meta.data[colnames(UCASpatial[[7]]),]

UCASpatial[[8]]@meta.data<-merge(UCASpatial[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(UCASpatial[[8]]@meta.data)<-UCASpatial[[8]]@meta.data$Barcode
UCASpatial[[8]]@meta.data<-UCASpatial[[8]]@meta.data[colnames(UCASpatial[[8]]),]


inhouse_cohort<-rbind(UCASpatial[[1]]@meta.data[,c(4:38)],UCASpatial[[4]]@meta.data[,c(4:38)],UCASpatial[[2]]@meta.data[,c(4:38)],UCASpatial[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(UCASpatial[[6]]@meta.data[,c(4:38)],UCASpatial[[5]]@meta.data[,c(4:38)],UCASpatial[[7]]@meta.data[,c(4:38)],UCASpatial[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(UCASpatial[[1]]@meta.data)),rep("CRC2",nrow(UCASpatial[[4]]@meta.data)),rep("CRC3a",nrow(UCASpatial[[2]]@meta.data)),rep("CRC3b",nrow(UCASpatial[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(UCASpatial[[6]]@meta.data)),rep("CRC5",nrow(UCASpatial[[5]]@meta.data)),rep("CRC6",nrow(UCASpatial[[7]]@meta.data)),rep("CRC7",nrow(UCASpatial[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.02的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)

tumor_prop_avg_UCASpatial<-tumor_prop_avg
saveRDS(tumor_prop_avg_UCASpatial,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_UCASpatial.RDS")


my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"

#Roe calculating
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"

library(dplyr)

total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()
Roe_v_inhouse <- Calculate_Roe(Roe_total)

Roe_celltype_list<-celltype_all[c(-8,-15,-16,-25,-27)]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  dplyr::select(Roe_celltype_list) %>% t() %>% as.data.frame()

library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
Roe_v_inhouse[is.na(Roe_v_inhouse)]<-0
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

p1<-ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+labs(x="Ro/e",y="Subset",title = "UCASpatial")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/UCASpatial_ROE.pdf",width = 5,height=8)

levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
p1<-pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                       gaps_col = c(3,8),scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                       annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8,main = "UCASpatial")
ggsave(plot = p1,filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/UCASpatial_heatmap.pdf",width = 5,height=8)

long_df <- pivot_longer(tumor_prop_avg, cols = 1:34, names_to = "Variable", values_to = "Value")

var=celltype_all[2]
plot_list <- lapply(celltype_all, function(var) {
  ggplot(subset(long_df, Variable == var),aes(x=subTME_category,y=Value,fill=subTME_category,color=subTME_category))+geom_boxplot(size=0.6)+
    scale_fill_manual(values = c("#A6CEE3","#FDBF6F"))+scale_color_manual(values = c("#418DBF","#FD8B18"))+
    labs(x="",y="Proportion",title = var)+
    stat_compare_means(label="p.signif",ref.group = "T cell-excluded")+
    theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+guides(fill=F,color=F)
})
ggsave(plot = do.call(gridExtra::grid.arrange, c(plot_list, ncol = 6)),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/UCASpatial_box.pdf",width = 20,height=40)




#SpatialFeaturePlot(UCASpatial_without_ent[[6]],features = celltype_all[1])

#clone region annotation
UCASpatial_without_ent[[1]]@meta.data<-merge(UCASpatial_without_ent[[1]]@meta.data,CRC1_Region,by="Barcode")
rownames(UCASpatial_without_ent[[1]]@meta.data)<-UCASpatial_without_ent[[1]]@meta.data$Barcode
UCASpatial_without_ent[[1]]@meta.data<-UCASpatial_without_ent[[1]]@meta.data[colnames(UCASpatial_without_ent[[1]]),]

UCASpatial_without_ent[[2]]@meta.data<-merge(UCASpatial_without_ent[[2]]@meta.data,CRC3a_Region,by="Barcode")
rownames(UCASpatial_without_ent[[2]]@meta.data)<-UCASpatial_without_ent[[2]]@meta.data$Barcode
UCASpatial_without_ent[[2]]@meta.data<-UCASpatial_without_ent[[2]]@meta.data[colnames(UCASpatial_without_ent[[2]]),]

UCASpatial_without_ent[[3]]@meta.data<-merge(UCASpatial_without_ent[[3]]@meta.data,CRC3b_Region,by="Barcode")
rownames(UCASpatial_without_ent[[3]]@meta.data)<-UCASpatial_without_ent[[3]]@meta.data$Barcode
UCASpatial_without_ent[[3]]@meta.data<-UCASpatial_without_ent[[3]]@meta.data[colnames(UCASpatial_without_ent[[3]]),]

UCASpatial_without_ent[[4]]@meta.data<-merge(UCASpatial_without_ent[[4]]@meta.data,CRC2_Region,by="Barcode")
rownames(UCASpatial_without_ent[[4]]@meta.data)<-UCASpatial_without_ent[[4]]@meta.data$Barcode
UCASpatial_without_ent[[4]]@meta.data<-UCASpatial_without_ent[[4]]@meta.data[colnames(UCASpatial_without_ent[[4]]),]

UCASpatial_without_ent[[5]]@meta.data<-merge(UCASpatial_without_ent[[5]]@meta.data,CRC5_Region,by="Barcode")
rownames(UCASpatial_without_ent[[5]]@meta.data)<-UCASpatial_without_ent[[5]]@meta.data$Barcode
UCASpatial_without_ent[[5]]@meta.data<-UCASpatial_without_ent[[5]]@meta.data[colnames(UCASpatial_without_ent[[5]]),]

colnames(UCASpatial_without_ent[[6]]@meta.data)[1]<-"Barcode"
UCASpatial_without_ent[[6]]@meta.data<-merge(UCASpatial_without_ent[[6]]@meta.data,CRC4_Region,by="Barcode")
rownames(UCASpatial_without_ent[[6]]@meta.data)<-UCASpatial_without_ent[[6]]@meta.data$Barcode
UCASpatial_without_ent[[6]]@meta.data<-UCASpatial_without_ent[[6]]@meta.data[colnames(UCASpatial_without_ent[[6]]),]

UCASpatial_without_ent[[7]]@meta.data<-merge(UCASpatial_without_ent[[7]]@meta.data,CRC6_Region,by="Barcode")
rownames(UCASpatial_without_ent[[7]]@meta.data)<-UCASpatial_without_ent[[7]]@meta.data$Barcode
UCASpatial_without_ent[[7]]@meta.data<-UCASpatial_without_ent[[7]]@meta.data[colnames(UCASpatial_without_ent[[7]]),]

UCASpatial_without_ent[[8]]@meta.data<-merge(UCASpatial_without_ent[[8]]@meta.data,CRC7_Region,by="Barcode")
rownames(UCASpatial_without_ent[[8]]@meta.data)<-UCASpatial_without_ent[[8]]@meta.data$Barcode
UCASpatial_without_ent[[8]]@meta.data<-UCASpatial_without_ent[[8]]@meta.data[colnames(UCASpatial_without_ent[[8]]),]


inhouse_cohort<-rbind(UCASpatial_without_ent[[1]]@meta.data[,c(4:38)],UCASpatial_without_ent[[4]]@meta.data[,c(4:38)],UCASpatial_without_ent[[2]]@meta.data[,c(4:38)],UCASpatial_without_ent[[3]]@meta.data[,c(4:38)])
public_cohort<-rbind(UCASpatial_without_ent[[6]]@meta.data[,c(4:38)],UCASpatial_without_ent[[5]]@meta.data[,c(4:38)],UCASpatial_without_ent[[7]]@meta.data[,c(4:38)],UCASpatial_without_ent[[8]]@meta.data[,c(4:38)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(UCASpatial_without_ent[[1]]@meta.data)),rep("CRC2",nrow(UCASpatial_without_ent[[4]]@meta.data)),rep("CRC3a",nrow(UCASpatial_without_ent[[2]]@meta.data)),rep("CRC3b",nrow(UCASpatial_without_ent[[3]]@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(UCASpatial_without_ent[[6]]@meta.data)),rep("CRC5",nrow(UCASpatial_without_ent[[5]]@meta.data)),rep("CRC6",nrow(UCASpatial_without_ent[[7]]@meta.data)),rep("CRC7",nrow(UCASpatial_without_ent[[8]]@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)
total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

# 将小于0.02的值赋值为0
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


total_sample_tumor_prop<-na.omit(total_sample_tumor_prop)
tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)
tumor_prop_avg_UCASpatial_without_ent<-tumor_prop_avg
saveRDS(tumor_prop_avg_UCASpatial_without_ent,"/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/tumor_prop_avg_UCASpatial_without_ent.RDS")
my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"






# PCC mac and T (each method)

tumor_prop_avg_RCTD





cor.test.result <- cor.test(tumor_prop_avg$`C1QC+ Mac`, tumor_prop_avg$`CD8+ Tem`)
p1<-ggplot(tumor_prop_avg,aes(x=`C1QC+ Mac`,y=`CD8+ Tem`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)

cor.test.result <- cor.test(tumor_prop_avg$`C1QC+ Mac`, tumor_prop_avg$`CD8+ Tex`)
p2<-ggplot(tumor_prop_avg,aes(x=`C1QC+ Mac`,y=`CD8+ Tex`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)
cor.test.result <- cor.test(tumor_prop_avg$`C1QC+ Mac`, tumor_prop_avg$`CD8+ Trm`)
p3<-ggplot(tumor_prop_avg,aes(x=`C1QC+ Mac`,y=`CD8+ Trm`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


cor.test.result <- cor.test(tumor_prop_avg$`FN1+ Mac`, tumor_prop_avg$`CD8+ Tem`)
p4<-ggplot(tumor_prop_avg,aes(x=`FN1+ Mac`,y=`CD8+ Tem`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)

cor.test.result <- cor.test(tumor_prop_avg$`FN1+ Mac`, tumor_prop_avg$`CD8+ Tex`)
p5<-ggplot(tumor_prop_avg,aes(x=`FN1+ Mac`,y=`CD8+ Tex`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)
cor.test.result <- cor.test(tumor_prop_avg$`FN1+ Mac`, tumor_prop_avg$`CD8+ Trm`)
p6<-ggplot(tumor_prop_avg,aes(x=`FN1+ Mac`,y=`CD8+ Trm`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


cor.test.result <- cor.test(tumor_prop_avg$`FN1+ Mac`, tumor_prop_avg$`CD8+ Tem`)
p7<-ggplot(tumor_prop_avg,aes(x=`PDPN+ Mac`,y=`CD8+ Tem`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)

cor.test.result <- cor.test(tumor_prop_avg$`PDPN+ Mac`, tumor_prop_avg$`CD8+ Tex`)
p8<-ggplot(tumor_prop_avg,aes(x=`PDPN+ Mac`,y=`CD8+ Tex`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)
cor.test.result <- cor.test(tumor_prop_avg$`PDPN+ Mac`, tumor_prop_avg$`CD8+ Trm`)
p9<-ggplot(tumor_prop_avg,aes(x=`PDPN+ Mac`,y=`CD8+ Trm`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


p1+p2+p3+p4+p5+p6+p7+p8+p9+plot_layout(ncol=3)




ggsave(plot = Seurat::SpatialFeaturePlot(
  object = test,
  features = cell_types_all,stroke = NA,alpha = c(0.3,1),
  min.cutoff = 0.03,ncol=9),filename = paste('UCASpatial_',z,'cut0.03.pdf',sep = ''),
  height = 12,width = 20)

SpatialFeaturePlot(object = RCTD[[6]],features = celltype_all[1],stroke = NA,alpha = c(0.1,1),
  min.cutoff = 0.03)+theme(legend.position = "right")+scale_fill_viridis(direction = 1)
celltype_all<-celltype_all[c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)]
library(gridExtra)
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(RCTD[[6]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = -1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})




# 使用do.call和grid.arrange来组合这些图形
library(gridExtra)



ggsave(plot = do.call(grid.arrange, c(plots, ncol = 6)),filename = 'RCTD_CRC4_cut0.03.pdf',
  height = 40,width = 25)


SpatialFeaturePlot(RCTD[[6]],features = celltype_all[1])