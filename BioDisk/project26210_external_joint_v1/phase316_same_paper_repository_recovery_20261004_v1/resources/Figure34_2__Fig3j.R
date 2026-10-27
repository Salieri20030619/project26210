
# neighbour-based cellular community
st_neigbour<-readRDS("/Path/data/st.list.neighbor.rds")
length(st_neigbour)
CRC1_prop<-as.data.frame(as.matrix(t(st_neigbour[[1]])))
CRC1_prop$barcode<-paste("CRC1",rownames(CRC1_prop),sep = "_")
CRC2_prop<-as.data.frame(as.matrix(t(st_neigbour[[2]])))
CRC2_prop$barcode<-paste("CRC2",rownames(CRC2_prop),sep = "_")
CRC3a_prop<-as.data.frame(as.matrix(t(st_neigbour[[3]])))
CRC3a_prop$barcode<-paste("CRC3a",rownames(CRC3a_prop),sep = "_")
CRC3b_prop<-as.data.frame(as.matrix(t(st_neigbour[[4]])))
CRC3b_prop$barcode<-paste("CRC3b",rownames(CRC3b_prop),sep = "_")
CRC4_prop<-as.data.frame(as.matrix(t(st_neigbour[[5]])))
CRC4_prop$barcode<-paste("CRC4",rownames(CRC4_prop),sep = "_")
CRC5_prop<-as.data.frame(as.matrix(t(st_neigbour[[6]])))
CRC5_prop$barcode<-paste("CRC5",rownames(CRC5_prop),sep = "_")
CRC6_prop<-as.data.frame(as.matrix(t(st_neigbour[[7]])))
CRC6_prop$barcode<-paste("CRC6",rownames(CRC6_prop),sep = "_")
CRC7_prop<-as.data.frame(as.matrix(t(st_neigbour[[8]])))
CRC7_prop$barcode<-paste("CRC7",rownames(CRC7_prop),sep = "_")

neighbour_prop<-rbind(CRC1_prop,CRC2_prop,CRC3a_prop,CRC3b_prop,CRC4_prop,CRC5_prop,CRC6_prop,CRC7_prop)
colnames(neighbour_prop)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(neighbour_prop))))))))))))
colnames(neighbour_prop)<-gsub("Fibroblast","CD34+ iFib",colnames(neighbour_prop))
colnames(neighbour_prop)[35]<-"section_barcode"

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


total_sample_region_info<-total_sample_porp[,c("sample","slice","Region")]
total_sample_region_info$section_barcode <- paste(total_sample_region_info$slice,rownames(total_sample_region_info),sep = "_")

neighbour_prop$section_barcode<-gsub("-","_",substr(neighbour_prop$section_barcode, 1, 23))
total_sample_region_info$section_barcode<-gsub("-","_",substr(total_sample_region_info$section_barcode, 1, 23))

neighbour_prop$slice<-sub('^([^_]+)_.*', '\\1', neighbour_prop$section_barcode)
table(total_sample_region_info$section_barcode %in% neighbour_prop$section_barcode)
table(neighbour_prop$section_barcode %in% total_sample_region_info$section_barcode )

neighbour_prop<-merge(neighbour_prop,total_sample_region_info,by="section_barcode")

table(neighbour_prop$Region)
neighbour_prop_tumor<-neighbour_prop[!(neighbour_prop$Region %in% c("Remaining","Normal_Epithelial")),]
colnames(neighbour_prop_tumor)
tumor_sample_prop_immune_neighbour<-neighbour_prop_tumor[,2:35]
tumor_sample_prop_immune_neighbour<-tumor_sample_prop_immune_neighbour[,!(colnames(tumor_sample_prop_immune_neighbour) %in% c("Epithelial c1","Epithelial c2","Epithelial c3","Epithelial c4","Epithelial c5"))]

PCC_matrix <- matrix(NA, 29, 29)
for(i in 1:29) {
  for(j in 1:29) {
    temp <- cor.test(tumor_sample_prop_immune_neighbour[,i], tumor_sample_prop_immune_neighbour[,j],method = "pearson")
    PCC_matrix[i,j] <- temp$estimate
  }
}

rownames(PCC_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])
colnames(PCC_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])

#PCC_matrix[PCC_matrix<0.1]<-0

sig_matrix <- matrix(NA, 29, 29)
for(i in 1:29) {
  for(j in 1:29) {
    temp <- cor.test(tumor_sample_prop_immune_neighbour[,i], tumor_sample_prop_immune_neighbour[,j],method = "pearson")
    sig_matrix[i,j] <- temp$p.value
  }
}

rownames(sig_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])
colnames(sig_matrix) = colnames(tumor_sample_prop_immune_neighbour[,1:29])


sig_matrix

convert_p_to_stars <- function(matrix_of_p_values) {
  stars <- matrix("", nrow = nrow(matrix_of_p_values), ncol = ncol(matrix_of_p_values))
  
  for (i in 1:nrow(matrix_of_p_values)) {
    for (j in 1:ncol(matrix_of_p_values)) {
      p <- matrix_of_p_values[i, j]
      if (!is.na(p)) {
        if (p < 0.0001) {
          stars[i, j] <- "****"
        } else if (p < 0.001) {
          stars[i, j] <- "***"
        } else if (p < 0.01) {
          stars[i, j] <- "**"
        } else if (p < 0.05) {
          stars[i, j] <- "*"
        } else {
          stars[i, j] <- ""
        }
      } else {
        stars[i, j] <- ""
      }
    }
  }
  
  return(stars)
}

# 示例用法
# 假设 matrix_of_p_values 是你的 p 值矩阵
stars_matrix <- convert_p_to_stars(sig_matrix)

stars_matrix[PCC_matrix<0.2]<-""

#cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value==1),]
#cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value <= 0.2),]

cell_cor_long_matrix<-reshape2::melt(PCC_matrix)
cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value==1),]
cell_cor_long_matrix<-cell_cor_long_matrix[!(cell_cor_long_matrix$value <= 0.3),]

write.csv(cell_cor_long_matrix,file = "/Path/data/cor_matrix_0.3.csv")

pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#7FB1D4","white","white","white","#FEBA78","#FE7F02"))(100),
                   cellwidth = 10,cellheight = 10,fontsize = 8,breaks = seq(-0.45, 0.45, length.out = 100),
                   border_color = "black",cutree_cols = 5,cutree_rows = 5,
                   clustering_method = "ward.D2",display_numbers = stars_matrix)