# TLS region
library(pacman)
p_unload(Seurat)
p_load(Seurat)

Celltype_num=34

CRC3a_st_TLS<-readRDS("/data/huangzr/Spatial/Kras/version4/boudary/Pt3_region.RDS")
SpatialPlot(CRC3a_st_TLS,group.by = "Region")
CRC3a_TLS_spot<-rownames(CRC3a_st_TLS@meta.data)[CRC3a_st_TLS@meta.data$Region == "TLS"]

CRC4_st_TLS<-readRDS("/data/huangzr/Spatial/Kras/version4/boudary/CRC2_region.RDS")
SpatialPlot(CRC4_st_TLS,group.by = "Region")
CRC4_TLS_spot<-rownames(CRC4_st_TLS@meta.data)[CRC4_st_TLS@meta.data$Region == "TLS"]

CRC5_st_TLS<-readRDS("/data/huangzr/Spatial/Kras/version4/boudary/CRC1_region.RDS")
SpatialPlot(CRC5_st_TLS,group.by = "Region")
CRC5_TLS_spot<-rownames(CRC5_st_TLS@meta.data)[CRC5_st_TLS@meta.data$Region == "TLS"]

CRC6_st_TLS<-readRDS("/data/huangzr/Spatial/Kras/version4/boudary/CRC3_region.RDS")
SpatialPlot(CRC6_st_TLS,group.by = "Region")
CRC6_TLS_spot<-rownames(CRC6_st_TLS@meta.data)[CRC6_st_TLS@meta.data$Region == "TLS"]

CRC7_st_TLS<-readRDS("/data/huangzr/Spatial/Kras/version4/boudary/CRC4_region.RDS")
SpatialPlot(CRC7_st_TLS,group.by = "Region")
CRC7_TLS_spot<-rownames(CRC7_st_TLS@meta.data)[CRC7_st_TLS@meta.data$Region == "TLS"]

#UCASpatial TLS Composition
UCASpatial<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/UCASpatial_results.rds")
colnames(UCASpatial[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[1]]@meta.data))))))))))))
colnames(UCASpatial[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[2]]@meta.data))))))))))))
colnames(UCASpatial[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[3]]@meta.data))))))))))))
colnames(UCASpatial[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[4]]@meta.data))))))))))))
colnames(UCASpatial[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[5]]@meta.data))))))))))))
colnames(UCASpatial[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[6]]@meta.data))))))))))))
colnames(UCASpatial[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[7]]@meta.data))))))))))))
colnames(UCASpatial[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[8]]@meta.data))))))))))))

CRC2_TLS_UCASpatial <- UCASpatial[[4]]@meta.data[rownames(UCASpatial[[4]]@meta.data) %in% CRC2_TLS_spot,4:37]
CRC2_TLS_UCASpatial$section="CRC2"
CRC3a_TLS_UCASpatial <- UCASpatial[[2]]@meta.data[rownames(UCASpatial[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_UCASpatial$section="CRC3a"
CRC4_TLS_UCASpatial <- UCASpatial[[6]]@meta.data[rownames(UCASpatial[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_UCASpatial$section="CRC4"
CRC5_TLS_UCASpatial <- UCASpatial[[5]]@meta.data[rownames(UCASpatial[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_UCASpatial$section="CRC5"
CRC6_TLS_UCASpatial <- UCASpatial[[7]]@meta.data[rownames(UCASpatial[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_UCASpatial$section="CRC6"
CRC7_TLS_UCASpatial <- UCASpatial[[8]]@meta.data[rownames(UCASpatial[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_UCASpatial$section="CRC7"

UCASpatial_TLS<-rbind(CRC2_TLS_UCASpatial,CRC3a_TLS_UCASpatial,CRC4_TLS_UCASpatial,CRC5_TLS_UCASpatial,CRC6_TLS_UCASpatial,CRC7_TLS_UCASpatial)

# 将小于0.02的值赋值为0
UCASpatial_TLS[,1:Celltype_num] <- apply(UCASpatial_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
UCASpatial_TLS[,1:Celltype_num] <- t(apply(UCASpatial_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


UCASpatial_TLS<-na.omit(UCASpatial_TLS)
UCASpatial_TLS_avg<-aggregate(UCASpatial_TLS[,1:Celltype_num],by = list(UCASpatial_TLS[,Celltype_num+1]),mean)
rownames(UCASpatial_TLS_avg)<- UCASpatial_TLS_avg$Group.1
UCASpatial_TLS_avg<-UCASpatial_TLS_avg[,2:ncol(UCASpatial_TLS_avg)] %>% t() %>% as.data.frame()
UCASpatial_TLS_avg<-t(UCASpatial_TLS_avg)
UCASpatial_TLS_avg_total<-colMeans(UCASpatial_TLS[,1:Celltype_num])
UCASpatial_TLS_avg<-reshape2::melt(UCASpatial_TLS_avg)
colnames(UCASpatial_TLS_avg)<-c("section","subpopulation","proportion")
UCASpatial_TLS_avg$method="UCASpatial"

#CARD TLS Composition
CARD<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/CARD_results.rds")
colnames(CARD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[1]]@meta.data))))))))))))
colnames(CARD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[2]]@meta.data))))))))))))
colnames(CARD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[3]]@meta.data))))))))))))
colnames(CARD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[4]]@meta.data))))))))))))
colnames(CARD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[5]]@meta.data))))))))))))
colnames(CARD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[6]]@meta.data))))))))))))
colnames(CARD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[7]]@meta.data))))))))))))
colnames(CARD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[8]]@meta.data))))))))))))


CRC3a_TLS_CARD <- CARD[[2]]@meta.data[rownames(CARD[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_CARD$section="CRC3a"
CRC4_TLS_CARD <- CARD[[6]]@meta.data[rownames(CARD[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_CARD$section="CRC4"
CRC5_TLS_CARD <- CARD[[5]]@meta.data[rownames(CARD[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_CARD$section="CRC5"
CRC6_TLS_CARD <- CARD[[7]]@meta.data[rownames(CARD[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_CARD$section="CRC6"
CRC7_TLS_CARD <- CARD[[8]]@meta.data[rownames(CARD[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_CARD$section="CRC7"

colnames(CRC3a_TLS_CARD) == colnames(CRC5_TLS_CARD)

CARD_TLS<-rbind(CRC3a_TLS_CARD,CRC4_TLS_CARD,CRC5_TLS_CARD,CRC6_TLS_CARD,CRC7_TLS_CARD)

# 将小于0.02的值赋值为0
CARD_TLS[,1:Celltype_num] <- apply(CARD_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
CARD_TLS[,1:Celltype_num] <- t(apply(CARD_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


CARD_TLS<-na.omit(CARD_TLS)
CARD_TLS_avg<-aggregate(CARD_TLS[,1:Celltype_num],by = list(CARD_TLS[,Celltype_num+1]),mean)
rownames(CARD_TLS_avg)<- CARD_TLS_avg$Group.1
CARD_TLS_avg<-CARD_TLS_avg[,2:ncol(CARD_TLS_avg)] %>% t() %>% as.data.frame()
CARD_TLS_avg<-t(CARD_TLS_avg)
CARD_TLS_avg_total<-colMeans(CARD_TLS[,1:Celltype_num])
CARD_TLS_avg<-reshape2::melt(CARD_TLS_avg)
colnames(CARD_TLS_avg)<-c("section","subpopulation","proportion")
CARD_TLS_avg$method="CARD"

#C2L TLS Composition
C2L<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/C2L_results.rds")
colnames(C2L[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[1]]@meta.data))))))))))))
colnames(C2L[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[2]]@meta.data))))))))))))
colnames(C2L[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[3]]@meta.data))))))))))))
colnames(C2L[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[4]]@meta.data))))))))))))
colnames(C2L[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[5]]@meta.data))))))))))))
colnames(C2L[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[6]]@meta.data))))))))))))
colnames(C2L[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[7]]@meta.data))))))))))))
colnames(C2L[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[8]]@meta.data))))))))))))

CRC3a_TLS_C2L <- C2L[[2]]@meta.data[rownames(C2L[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_C2L$section="CRC3a"
CRC4_TLS_C2L <- C2L[[6]]@meta.data[rownames(C2L[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_C2L$section="CRC4"
CRC5_TLS_C2L <- C2L[[5]]@meta.data[rownames(C2L[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_C2L$section="CRC5"
CRC6_TLS_C2L <- C2L[[7]]@meta.data[rownames(C2L[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_C2L$section="CRC6"
CRC7_TLS_C2L <- C2L[[8]]@meta.data[rownames(C2L[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_C2L$section="CRC7"

C2L_TLS<-rbind(CRC3a_TLS_C2L,CRC4_TLS_C2L,CRC5_TLS_C2L,CRC6_TLS_C2L,CRC7_TLS_C2L)

# 将小于0.02的值赋值为0
C2L_TLS[,1:Celltype_num] <- apply(C2L_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
C2L_TLS[,1:Celltype_num] <- t(apply(C2L_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


C2L_TLS<-na.omit(C2L_TLS)
C2L_TLS_avg<-aggregate(C2L_TLS[,1:Celltype_num],by = list(C2L_TLS[,Celltype_num+1]),mean)
rownames(C2L_TLS_avg)<- C2L_TLS_avg$Group.1
C2L_TLS_avg<-C2L_TLS_avg[,2:ncol(C2L_TLS_avg)] %>% t() %>% as.data.frame()
C2L_TLS_avg<-t(C2L_TLS_avg)
C2L_TLS_avg_total<-colMeans(C2L_TLS[,1:Celltype_num])
C2L_TLS_avg<-reshape2::melt(C2L_TLS_avg)
colnames(C2L_TLS_avg)<-c("section","subpopulation","proportion")
C2L_TLS_avg$method="C2L"

#RCTD TLS Composition
RCTD<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/RCTD_results.rds")
colnames(RCTD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[1]]@meta.data))))))))))))
colnames(RCTD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[2]]@meta.data))))))))))))
colnames(RCTD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[3]]@meta.data))))))))))))
colnames(RCTD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[4]]@meta.data))))))))))))
colnames(RCTD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[5]]@meta.data))))))))))))
colnames(RCTD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[6]]@meta.data))))))))))))
colnames(RCTD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[7]]@meta.data))))))))))))
colnames(RCTD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[8]]@meta.data))))))))))))

CRC3a_TLS_RCTD <- RCTD[[2]]@meta.data[rownames(RCTD[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_RCTD$section="CRC3a"
CRC4_TLS_RCTD <- RCTD[[6]]@meta.data[rownames(RCTD[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_RCTD$section="CRC4"
CRC5_TLS_RCTD <- RCTD[[5]]@meta.data[rownames(RCTD[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_RCTD$section="CRC5"
CRC6_TLS_RCTD <- RCTD[[7]]@meta.data[rownames(RCTD[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_RCTD$section="CRC6"
CRC7_TLS_RCTD <- RCTD[[8]]@meta.data[rownames(RCTD[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_RCTD$section="CRC7"

RCTD_TLS<-rbind(CRC3a_TLS_RCTD,CRC4_TLS_RCTD,CRC5_TLS_RCTD,CRC6_TLS_RCTD,CRC7_TLS_RCTD)

# 将小于0.02的值赋值为0
RCTD_TLS[,1:Celltype_num] <- apply(RCTD_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
RCTD_TLS[,1:Celltype_num] <- t(apply(RCTD_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


RCTD_TLS<-na.omit(RCTD_TLS)
RCTD_TLS_avg<-aggregate(RCTD_TLS[,1:Celltype_num],by = list(RCTD_TLS[,Celltype_num+1]),mean)
rownames(RCTD_TLS_avg)<- RCTD_TLS_avg$Group.1
RCTD_TLS_avg<-RCTD_TLS_avg[,2:ncol(RCTD_TLS_avg)] %>% t() %>% as.data.frame()
RCTD_TLS_avg<-t(RCTD_TLS_avg)
RCTD_TLS_avg_total<-colMeans(RCTD_TLS[,1:Celltype_num])
RCTD_TLS_avg<-reshape2::melt(RCTD_TLS_avg)
colnames(RCTD_TLS_avg)<-c("section","subpopulation","proportion")
RCTD_TLS_avg$method="RCTD"

#stereoscope TLS Composition
stereoscope<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/stereoscope_results.rds")
colnames(stereoscope[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[1]]@meta.data))))))))))))
colnames(stereoscope[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[2]]@meta.data))))))))))))
colnames(stereoscope[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[3]]@meta.data))))))))))))
colnames(stereoscope[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[4]]@meta.data))))))))))))
colnames(stereoscope[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[5]]@meta.data))))))))))))
colnames(stereoscope[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[6]]@meta.data))))))))))))
colnames(stereoscope[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[7]]@meta.data))))))))))))
colnames(stereoscope[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[8]]@meta.data))))))))))))

CRC3a_TLS_stereoscope <- stereoscope[[2]]@meta.data[rownames(stereoscope[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_stereoscope$section="CRC3a"
CRC4_TLS_stereoscope <- stereoscope[[6]]@meta.data[rownames(stereoscope[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_stereoscope$section="CRC4"
CRC5_TLS_stereoscope <- stereoscope[[5]]@meta.data[rownames(stereoscope[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_stereoscope$section="CRC5"
CRC6_TLS_stereoscope <- stereoscope[[7]]@meta.data[rownames(stereoscope[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_stereoscope$section="CRC6"
CRC7_TLS_stereoscope <- stereoscope[[8]]@meta.data[rownames(stereoscope[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_stereoscope$section="CRC7"

stereoscope_TLS<-rbind(CRC3a_TLS_stereoscope,CRC4_TLS_stereoscope,CRC5_TLS_stereoscope,CRC6_TLS_stereoscope,CRC7_TLS_stereoscope)

# 将小于0.02的值赋值为0
stereoscope_TLS[,1:Celltype_num] <- apply(stereoscope_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
stereoscope_TLS[,1:Celltype_num] <- t(apply(stereoscope_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


stereoscope_TLS<-na.omit(stereoscope_TLS)
stereoscope_TLS_avg<-aggregate(stereoscope_TLS[,1:Celltype_num],by = list(stereoscope_TLS[,Celltype_num+1]),mean)
rownames(stereoscope_TLS_avg)<- stereoscope_TLS_avg$Group.1
stereoscope_TLS_avg<-stereoscope_TLS_avg[,2:ncol(stereoscope_TLS_avg)] %>% t() %>% as.data.frame()
stereoscope_TLS_avg<-t(stereoscope_TLS_avg)
stereoscope_TLS_avg_total<-colMeans(stereoscope_TLS[,1:Celltype_num])
stereoscope_TLS_avg<-reshape2::melt(stereoscope_TLS_avg)
colnames(stereoscope_TLS_avg)<-c("section","subpopulation","proportion")
stereoscope_TLS_avg$method="stereoscope"

#SPOTlight TLS Composition
SPOTlight<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/SPOTlight_results.rds")
colnames(SPOTlight[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[1]]@meta.data))))))))))))
colnames(SPOTlight[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[2]]@meta.data))))))))))))
colnames(SPOTlight[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[3]]@meta.data))))))))))))
colnames(SPOTlight[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[4]]@meta.data))))))))))))
colnames(SPOTlight[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[5]]@meta.data))))))))))))
colnames(SPOTlight[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[6]]@meta.data))))))))))))
colnames(SPOTlight[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[7]]@meta.data))))))))))))
colnames(SPOTlight[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[8]]@meta.data))))))))))))

CRC3a_TLS_SPOTlight <- SPOTlight[[2]]@meta.data[rownames(SPOTlight[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_SPOTlight$section="CRC3a"
CRC4_TLS_SPOTlight <- SPOTlight[[6]]@meta.data[rownames(SPOTlight[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_SPOTlight$section="CRC4"
CRC5_TLS_SPOTlight <- SPOTlight[[5]]@meta.data[rownames(SPOTlight[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_SPOTlight$section="CRC5"
CRC6_TLS_SPOTlight <- SPOTlight[[7]]@meta.data[rownames(SPOTlight[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_SPOTlight$section="CRC6"
CRC7_TLS_SPOTlight <- SPOTlight[[8]]@meta.data[rownames(SPOTlight[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_SPOTlight$section="CRC7"

SPOTlight_TLS<-rbind(CRC3a_TLS_SPOTlight,CRC4_TLS_SPOTlight,CRC5_TLS_SPOTlight,CRC6_TLS_SPOTlight,CRC7_TLS_SPOTlight)

# 将小于0.02的值赋值为0
SPOTlight_TLS[,1:Celltype_num] <- apply(SPOTlight_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
SPOTlight_TLS[,1:Celltype_num] <- t(apply(SPOTlight_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


SPOTlight_TLS<-na.omit(SPOTlight_TLS)
SPOTlight_TLS_avg<-aggregate(SPOTlight_TLS[,1:Celltype_num],by = list(SPOTlight_TLS[,Celltype_num+1]),mean)
rownames(SPOTlight_TLS_avg)<- SPOTlight_TLS_avg$Group.1
SPOTlight_TLS_avg<-SPOTlight_TLS_avg[,2:ncol(SPOTlight_TLS_avg)] %>% t() %>% as.data.frame()
SPOTlight_TLS_avg<-t(SPOTlight_TLS_avg)
SPOTlight_TLS_avg_total<-colMeans(SPOTlight_TLS[,1:Celltype_num])
SPOTlight_TLS_avg<-reshape2::melt(SPOTlight_TLS_avg)
colnames(SPOTlight_TLS_avg)<-c("section","subpopulation","proportion")
SPOTlight_TLS_avg$method="SPOTlight"


#UCASpatial_without_ent TLS Composition
UCASpatial_without_ent<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240517_UCASpatial_in_real_without_weight/without_weight_and_meta/UCASpatial_nd_CRC.rds")
colnames(UCASpatial_without_ent[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[1]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[2]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[3]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[4]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[5]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[6]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[7]]@meta.data))))))))))))
colnames(UCASpatial_without_ent[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial_without_ent[[8]]@meta.data))))))))))))

CRC3a_TLS_UCASpatial_without_ent <- UCASpatial_without_ent[[2]]@meta.data[rownames(UCASpatial_without_ent[[2]]@meta.data) %in% CRC3a_TLS_spot,4:37]
CRC3a_TLS_UCASpatial_without_ent$section="CRC3a"
CRC4_TLS_UCASpatial_without_ent <- UCASpatial_without_ent[[6]]@meta.data[rownames(UCASpatial_without_ent[[6]]@meta.data) %in% CRC4_TLS_spot,4:37]
CRC4_TLS_UCASpatial_without_ent$section="CRC4"
CRC5_TLS_UCASpatial_without_ent <- UCASpatial_without_ent[[5]]@meta.data[rownames(UCASpatial_without_ent[[5]]@meta.data) %in% CRC5_TLS_spot,4:37]
CRC5_TLS_UCASpatial_without_ent$section="CRC5"
CRC6_TLS_UCASpatial_without_ent <- UCASpatial_without_ent[[7]]@meta.data[rownames(UCASpatial_without_ent[[7]]@meta.data) %in% CRC6_TLS_spot,4:37]
CRC6_TLS_UCASpatial_without_ent$section="CRC6"
CRC7_TLS_UCASpatial_without_ent <- UCASpatial_without_ent[[8]]@meta.data[rownames(UCASpatial_without_ent[[8]]@meta.data) %in% CRC7_TLS_spot,4:37]
CRC7_TLS_UCASpatial_without_ent$section="CRC7"

UCASpatial_without_ent_TLS<-rbind(CRC3a_TLS_UCASpatial_without_ent,CRC4_TLS_UCASpatial_without_ent,CRC5_TLS_UCASpatial_without_ent,CRC6_TLS_UCASpatial_without_ent,CRC7_TLS_UCASpatial_without_ent)

# 将小于0.02的值赋值为0
UCASpatial_without_ent_TLS[,1:Celltype_num] <- apply(UCASpatial_without_ent_TLS[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
UCASpatial_without_ent_TLS[,1:Celltype_num] <- t(apply(UCASpatial_without_ent_TLS[,1:Celltype_num], 1, function(x) x / sum(x)))


UCASpatial_without_ent_TLS<-na.omit(UCASpatial_without_ent_TLS)
UCASpatial_without_ent_TLS_avg<-aggregate(UCASpatial_without_ent_TLS[,1:Celltype_num],by = list(UCASpatial_without_ent_TLS[,Celltype_num+1]),mean)
rownames(UCASpatial_without_ent_TLS_avg)<- UCASpatial_without_ent_TLS_avg$Group.1
UCASpatial_without_ent_TLS_avg<-UCASpatial_without_ent_TLS_avg[,2:ncol(UCASpatial_without_ent_TLS_avg)] %>% t() %>% as.data.frame()
UCASpatial_without_ent_TLS_avg<-t(UCASpatial_without_ent_TLS_avg)
UCASpatial_without_ent_TLS_avg_total<-colMeans(UCASpatial_without_ent_TLS[,1:Celltype_num])
UCASpatial_without_ent_TLS_avg<-reshape2::melt(UCASpatial_without_ent_TLS_avg)
colnames(UCASpatial_without_ent_TLS_avg)<-c("section","subpopulation","proportion")
UCASpatial_without_ent_TLS_avg$method="UCASpatial_without_ent"






colnames(UCASpatial_without_ent_TLS_avg)<-colnames(SPOTlight_TLS_avg)
section_TLS_prop<-rbind(UCASpatial_TLS_avg,UCASpatial_without_ent_TLS_avg,SPOTlight_TLS_avg,RCTD_TLS_avg,C2L_TLS_avg,CARD_TLS_avg,stereoscope_TLS_avg)
theme_set(theme_cowplot())

section_TLS_prop$method<-factor(section_TLS_prop$method,levels=c("UCASpatial","UCASpatial_without_ent","SPOTlight","RCTD","C2L","CARD","stereoscope"))
section_TLS_prop$subpopulation<-factor(section_TLS_prop$subpopulation,levels=c("CD8+ Tem","CD8+ Tex","CD8+ Trm", "CD4+ Tn" , "CD4+ Tm", "Tfh", "Th1/Th17", "Treg","NK", "ILC","B","Plasma","Mast","Monocyte", "Neutrophil",  "C1QC+ Mac", "FN1+ Mac",  "PDPN+ Mac",  "cDC1", "cDC2","cDC3" ,"pDC",  "Endothelial",  "Glial", "Fibroblast", "FAP+ Myofib", "POSTN+ Myofib","Smooth Muscle", "Pericyte", "Epithelial c1","Epithelial c2","Epithelial c3","Epithelial c4",
                                                                               "Epithelial c5"))






levels(section_TLS_prop$method)[2]<-"UCASpatial_ablation"

levels(section_TLS_prop$subpopulation)

colors=c("#F9D4C7","#F5B09B","#F58773","#F16A54","#E83F3B","#C7242D","#922629","#5A1312","#F7A06D","#E2662C","#D3B67E","#EBDAB6",
"#983D2B","#E1DEE9","#D0E1EE","#B6CEE4","#8CBED2","#62A2C6","#448AB6","#2A70AB","#1B5191","#2A3864","#82522C","#B47C3C","#CCCFE2",
"#A6ABC9","#7779B0","#64569C","#403066","#B7E0DA","#71C0B6","#408F8A","#1F655B","#143930")


theme_set(theme_cowplot())
ggplot(section_TLS_prop[section_TLS_prop$section != "CRC2",],aes(x=section,y=proportion,fill=subpopulation))+geom_bar(stat = "identity")+
  scale_fill_manual(values = colors)+facet_wrap(~method,ncol = 4)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


ggplot(section_TLS_prop,aes(x=section,y=proportion,fill=subpopulation))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c(brewer.pal(12,"Paired"),brewer.pal(12,"Set2"),brewer.pal(12,"Set1"),brewer.pal(12,"Set3")))+
  facet_wrap(~method,ncol = 3)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
CRC2_st_TLS@meta.data$cor<-F
CRC2_st_TLS@meta.data$cor[CRC2_st_TLS@images$sliceC1@coordinates$row>48]<-T
p1<-SpatialPlot(CRC2_st_TLS,group.by = "cor")

p2<-SpatialPlot(CRC2_st_TLS,group.by = "Region")
p1+p2+plot_layout(ncol=1)


SpatialPlot(SPOTlight[[5]],features = "CD8+ Tem")

min(CRC5_st_TLS[,CRC5_TLS_spot]@images$image@coordinates$col)

row:28:40
col:85:115

CRC_5_TLS_round<-rownames(CRC5_st_TLS@images$image@coordinates)[CRC5_st_TLS@images$image@coordinates$row > 27 & CRC5_st_TLS@images$image@coordinates$row < 43 & CRC5_st_TLS@images$image@coordinates$col > 90 & CRC5_st_TLS@images$image@coordinates$col < 110]
CRC_5_TLS_round %in% colnames(RCTD[[5]])

SubsetData(RCTD[[5]],CRC_5_TLS_round)

RCTD[[5]]@meta.data$TLS_around <- rownames(RCTD[[5]]@meta.data) %in% CRC_5_TLS_round
SubsetData(RCTD[[5]],subset.name = TLS_around)
Seurat::SubsetData(RCTD[[5]], CRC_5_TLS_round)
SpatialFeaturePlot(RCTD[[5]][,RCTD[[5]]@meta.data$TLS_around == T], features = "CD8+ Tem", stroke = NA,
                   min.cutoff = 0.03,alpha = c(0.1,1)) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))

SpatialFeaturePlot(RCTD[[5]][,RCTD[[5]]@images$sliceD1@coordinates$row>30], features = "CD8+ Tem", stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))

subset(RCTD[[5]], imagerow > 4000 & imagerow < 8000 & imagecol < 5000 & imagecol > 4000)
saveRDS(CRC_5_TLS_round,file = "/data/huangzr/Spatial/Kras/version5/SI_Fig/CRC5_TLS_around.rds")
library(SeuratData)
UCAS_CRC5_TLS<-SubsetData(UCASpatial[[5]],CRC_5_TLS_round)
SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                   min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))
CRC_5_TLS_round<-readRDS("/data/huangzr/Spatial/Kras/version5/SI_Fig/CRC5_TLS_around.rds")

p_unload(Seurat)
p_load(Seurat)

UCAS_CRC5_TLS<-UCASpatial[[5]][,CRC_5_TLS_round]
UCASpatial_without_ent_CRC5_TLS<-UCASpatial_without_ent[[5]][,CRC_5_TLS_round]
SPOTlight_CRC5_TLS<-SPOTlight[[5]][,CRC_5_TLS_round]
CARD_CRC5_TLS<-CARD[[5]][,CRC_5_TLS_round]
C2L_CRC5_TLS<-C2L[[5]][,CRC_5_TLS_round]
stereoscope_CRC5_TLS<-stereoscope[[5]][,CRC_5_TLS_round]


RCTD[[5]]@meta.data<-RCTD[[5]]@meta.data[colnames(RCTD[[5]]),]
RCTD_CRC5_TLS<-RCTD[[5]][,colnames(RCTD[[5]]) %in% CRC_5_TLS_round]
SPOTlight_CRC5_TLS<-SPOTlight[[5]][,colnames(SPOTlight[[5]]) %in% CRC_5_TLS_round]
CARD_CRC5_TLS<-CARD[[5]][,colnames(CARD[[5]]) %in% CRC_5_TLS_round]
C2L_CRC5_TLS<-C2L[[5]][,colnames(C2L[[5]]) %in% CRC_5_TLS_round]
stereoscope_CRC5_TLS<-stereoscope[[5]][,colnames(stereoscope[[5]]) %in% CRC_5_TLS_round]

DefaultAssay(UCAS_CRC5_TLS)<-"Spatial"
RCTD[[5]]@assays$Annotation<-NA

colnames(RCTD_CRC5_TLS@meta.data)[7]<-"CD8+ Tem"

p1<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                   min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
 theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")
p2<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")
p3<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p4<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("Cell2location")
p5<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")
p6<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")

p7<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD8A", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))
p8<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "GZMK", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))
p9<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 6) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_ablation")



p7+p1+p9+p2+p3+p8+p4+p5+p6+plot_layout(ncol=5)
p4


SpatialPlot(UCASpatial[[3]],features = "CD8+ Trm",cols = c("white","#7A040E"))


SpatialPlot(UCASpatial[[3]],features = "GZMK",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FEFAF8',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#FEFAF8")+
  labs(fill ='GZMK')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p1<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")
p2<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")
p3<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")

p4<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")
p5<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")

p6<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")
p7<-SpatialFeaturePlot(UCAS_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
UCASpatial_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")
p2<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")
p3<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")

p4<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")
p5<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")

p6<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")
p7<-SpatialFeaturePlot(CARD_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("CARD")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
CARD_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')

p1<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p2<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p3<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")

p4<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p5<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")

p6<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p7<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
RCTD_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p2<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p3<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")

p4<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p5<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")

p6<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p7<-SpatialFeaturePlot(RCTD_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("RCTD")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
RCTD_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')

p1<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")
  
  
p2<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")
p3<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")

p4<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")
p5<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")

p6<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")
p7<-SpatialFeaturePlot(SPOTlight_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("SPOTlight")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
SPOTlight_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")

p2<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")
p3<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")

p4<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")
p5<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")

p6<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")
p7<-SpatialFeaturePlot(stereoscope_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("stereoscope")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
stereoscope_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')



p1<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")

p2<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")
p3<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")

p4<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")
p5<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")

p6<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")
p7<-SpatialFeaturePlot(C2L_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("C2L")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
C2L_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "CD8+ Tem", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")

p2<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "CD8+ Tex", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")
p3<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "CD8+ Trm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")

p4<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "CD4+ Tn", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")
p5<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "CD4+ Tm", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")

p6<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "B", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")
p7<-SpatialFeaturePlot(UCASpatial_without_ent_CRC5_TLS, features = "cDC2", stroke = NA,
                       min.cutoff = 0.03,alpha = c(0.1,1),pt.size.factor = 7) +
  theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust=0.5))+ggtitle("UCASpatial_without_ent")
p_list <- list(p1,p2,p3,p4,p5,p6,p7)
UCASpatial_without_ent_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')

plot_grid(UCAS_p,UCASpatial_without_ent_p,CARD_p,C2L_p,RCTD_p,SPOTlight_p,stereoscope_p,
          align = 'hv',nrow = 1)
