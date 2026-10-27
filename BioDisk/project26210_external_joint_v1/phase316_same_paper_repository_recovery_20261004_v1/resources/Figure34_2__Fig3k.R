#community ClonalTIME association
library(dplyr)

CRC1_st<-readRDS("/Path/data/CRC1_st_clone.RDS")
CRC2_st<-readRDS("/Path/data/CRC2_st_clone.RDS")
CRC3a_st<-readRDS("/Path/data/CRC3a_st_clone.RDS")
CRC3b_st<-readRDS("/Path/data/CRC3b_st_clone.RDS")
CRC4_st<-readRDS("/Path/data/CRC4_st_clone.RDS")
CRC5_st<-readRDS("/Path/data/CRC5_st_clone.RDS")
CRC6_st<-readRDS("/Path/data/CRC6_st_clone.RDS")
CRC7_st<-readRDS("/Path/data/CRC7_st_clone.RDS")
Celltype_num=34

inhouse_cohort<-rbind(CRC1_st@meta.data[,c(4:37,40)],CRC2_st@meta.data[,c(4:37,40)],CRC3a_st@meta.data[,c(4:37,40)],CRC3b_st@meta.data[,c(4:37,40)])
public_cohort<-rbind(CRC4_st@meta.data[,c(5:38,41)],CRC5_st@meta.data[,c(4:37,40)],CRC6_st@meta.data[,c(4:37,40)],CRC7_st@meta.data[,c(4:37,40)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(CRC1_st@meta.data)),rep("CRC2",nrow(CRC2_st@meta.data)),rep("CRC3a",nrow(CRC3a_st@meta.data)),rep("CRC3b",nrow(CRC3b_st@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(CRC4_st@meta.data)),rep("CRC5",nrow(CRC5_st@meta.data)),rep("CRC6",nrow(CRC6_st@meta.data)),rep("CRC7",nrow(CRC7_st@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)

total_sample_porp<-rbind(public_cohort,inhouse_cohort)
#total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]

total_sample_porp[,1:Celltype_num] <- apply(total_sample_porp[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))
colnames(total_sample_porp)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(total_sample_porp))))))))))))
colnames(total_sample_porp)<-gsub("Fibroblast","CD34+ iFib",colnames(total_sample_porp))

#total_sample_porp_immune<-total_sample_porp[,c(1:34,35:37,81)]
total_sample_porp_immune<-total_sample_porp[,!(colnames(total_sample_porp) %in% c("Epithelial c1","Epithelial c2","Epithelial c3","Epithelial c4","Epithelial c5"))]
total_sample_porp_immune$total<-rowSums(total_sample_porp_immune[,1:29])

mat_normalized <- t(apply(total_sample_porp_immune[,1:29], MARGIN = 1, FUN = function(row) row / sum(row)))

print(mat_normalized)
colnames(mat_normalized)<-paste(colnames(mat_normalized),"normalized",sep = "_")
upper_quantile <- apply(mat_normalized, 2, quantile, probs = 0.995,na.rm=T)
lower_quantile <- apply(mat_normalized, 2, quantile, probs = 0.005,na.rm=T)
scaled_data <- scale(mat_normalized, center = lower_quantile, scale = upper_quantile - lower_quantile)
scaled_data<-as.data.frame(scaled_data)
colnames(scaled_data)<-gsub("normalized","scaled",colnames(scaled_data))
total_sample_porp_immune<-cbind(total_sample_porp_immune,mat_normalized,scaled_data)

total_sample_porp_immune$community1<-rowMeans(total_sample_porp_immune[,c("C1QC+ Mac_scaled","CD8+ Tem_scaled","CD34+ iFib_scaled","CD8+ Tex_scaled")])
total_sample_porp_immune$community2<-rowMeans(total_sample_porp_immune[,c("FN1+ Mac_scaled","Glial_scaled","Endothelial_scaled","PDPN+ Mac_scaled","FAP+ Myofib_scaled","Pericyte_scaled","Smooth Muscle_scaled")])
total_sample_porp_immune$community3<-rowMeans(total_sample_porp_immune[,c("Monocyte_scaled","Neutrophil_scaled")])


total_sample_porp_immune2<-total_sample_porp_immune[,c(32,92:94)]
total_sample_porp_immune2<-na.omit(total_sample_porp_immune2)
tumor_prop_avg<-aggregate(total_sample_porp_immune2[,2:4],by = list(total_sample_porp_immune2[,1]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)
tumor_prop_avg<-tumor_prop_avg[grep("Clone",rownames(tumor_prop_avg)),]
tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))
library(ggpubr)
theme_set(theme_cowplot())
p1<-ggplot(tumor_prop_avg,aes(x=subTME_category,y=community1,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5,method = "t.test")+
  labs(x="ClonalTIME",y="Community Score",title = "Community1")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 0.13))
p2<-ggplot(tumor_prop_avg,aes(x=subTME_category,y=community2,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5,method = "t.test")+
  labs(x="ClonalTIME",y="Community Score",title = "Community2")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 0.18))
p3<-ggplot(tumor_prop_avg,aes(x=subTME_category,y=community3,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5,method = "t.test")+
  labs(x="ClonalTIME",y="Community Score",title = "Community3")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 0.16))

library("patchwork")
p1+p2+p3+plot_layout(ncol=2)

saveRDS(total_sample_porp_immune2,file="/Path/data/total_sample_porp_immune2.RDS")

